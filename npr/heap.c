#include "npr/heap.h"
#include "npr/bits.h"
#include <assert.h>

#define ALIGN_UP(v,align) ((v)+((align)-1))&~((align)-1)
#define ALIGN_UP(v,align) ((v)+((align)-1))&~((align)-1)

#define ALIGN_DOWN(v, align) (v)&~((align)-1)

#ifdef _WIN32
#include <windows.h>
static int get_page_size() {
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    return si.dwPageSize;
}

static const int RWE = PAGE_EXECUTE_READWRITE;
static const int RW = PAGE_READWRITE;

static uintptr_t
alloc_page(size_t sz, int prot) {
    return (uintptr_t)VirtualAlloc(NULL, sz, MEM_COMMIT, prot);
}
static void
free_page(void *p, size_t sz) {
    VirtualFree(p, sz, MEM_FREE);
}

#else
#include <sys/mman.h>
#include <unistd.h>


static int get_page_size() {
    return sysconf(_SC_PAGE_SIZE);
}
static const int RWE = PROT_READ | PROT_WRITE | PROT_EXEC;
static const int RW = PROT_READ | PROT_WRITE;

static uintptr_t
alloc_page(size_t sz, int prot) {
    return (uintptr_t)mmap(NULL, sz, prot, MAP_ANON|MAP_PRIVATE, 0, 0);
}
static void
free_page(void *p, size_t sz) {
    munmap(p, sz);
}
#endif

static const int SIZEOF_PAGE_HEADER = 32;

typedef uint32_t chunksize_t;
static const chunksize_t NEXT_ALLOCATED_BIT = 0x80000000;
static const chunksize_t THIS_ALLOCATED_BIT = 0x40000000;
static const chunksize_t ALLOC_BIT_MASK = 0xc0000000;
static const chunksize_t SIZE_MASK = 0x3fffffff;
static const int CHUNK_ALIGN = 8;

void
npr_heap_init(struct npr_heap *h, int exec)
{
    int ps, shift;
    if (exec) {
        h->flags = RWE;
    } else {
        h->flags = RW;
    }

    ps = get_page_size();
    h->chunk_size = h->page_size = ps;
    assert(npr_popcnt32(ps) == 1);
    shift = npr_bsf32(ps);
    h->pfn_mask = ~((((uintptr_t)1)<<shift)-1);
    h->num_ptr_per_page = ps / (sizeof (void*));

    h->page_head.next = &h->page_tail;
    h->page_tail.prev = &h->page_head;
    h->free_head.next = &h->free_tail;
    h->free_tail.prev = &h->free_head;
}

static void
add_page(struct npr_heap *h, struct npr_heap_page_header *ph)
{
    ph->next = h->page_head.next;
    h->page_head.next->prev = ph;

    ph->prev = &h->page_head;
    h->page_head.next = ph;
}

static void
del_page(struct npr_heap *h, struct npr_heap_page_header *ph)
{
    ph->prev->next = ph->next;
    ph->next->prev = ph->prev;
}

static void
add_free(struct npr_heap *h, struct npr_free_chunk_header *f)
{
    h->free_tail.prev->next = f;
    f->prev = h->free_tail.prev;
    f->next = &h->free_tail;
    h->free_tail.prev = f;
}
static void
del_free(struct npr_free_chunk_header *f)
{
    f->prev->next = f->next;
    f->next->prev = f->prev;
}


static __inline size_t
calc_alloc_size(size_t n)
{
    return ALIGN_UP(n+4, CHUNK_ALIGN);
}

void *
npr_heap_alloc(struct npr_heap *h, size_t n)
{
    size_t alloc_size = calc_alloc_size(n);
    if (alloc_size >= (h->page_size-128)) {
        /* large mem */
        uintptr_t page;
        struct npr_heap_page_header *ph;
        alloc_size = ALIGN_UP(alloc_size+SIZEOF_PAGE_HEADER, h->page_size);
        page = alloc_page(alloc_size, h->flags);
        ph = (struct npr_heap_page_header*)page;
        ph->flags = NPR_LARGE_CHUNK;
        ph->page_size = alloc_size;
        add_page(h, ph);
        return (void*)(page+SIZEOF_PAGE_HEADER);
    } else {
        struct npr_free_chunk_header *f = h->free_head.next;
        struct npr_heap_page_header *ph;

        while (f != &h->free_tail) {
            if (f->sz >= alloc_size) {
                del_free(f);
                goto found;
            }
            f = f->next;
        }

        /* not found */
        {
            size_t cs = h->chunk_size;
            uintptr_t page = alloc_page(cs, h->flags);
            h->chunk_size *= 2;
            ph = (struct npr_heap_page_header*)page;
            ph->flags = 0;
            add_page(h, ph);
            f = (struct npr_free_chunk_header*)(page+SIZEOF_PAGE_HEADER);
            f->sz = cs-SIZEOF_PAGE_HEADER;
        }

    found:
        /* |<-    alloc  size    ->|
         * |           |<- 4byte ->|
         *
         * |  data     |   size    |
         *  ^           ^           ^
         *  |           |           |
         *  f_addr      size_addr   next_addr
         */

        {
            uintptr_t f_addr = (uintptr_t)f;
            uintptr_t next_addr = f_addr + alloc_size;
            struct npr_free_chunk_header *fn;
            uintptr_t size_addr = next_addr - sizeof(chunksize_t);
            *(chunksize_t*)size_addr = alloc_size | THIS_ALLOCATED_BIT;

            if (f->sz > alloc_size) {
                fn = (struct npr_free_chunk_header*)next_addr;
                fn->sz = f->sz - alloc_size;
                add_free(h, fn);
            } else {
                fn = (struct npr_free_chunk_header*)f->next;
            }

            if ((f_addr & ~(h->pfn_mask)) != SIZEOF_PAGE_HEADER) {
                uintptr_t prev_size = f_addr - sizeof(chunksize_t);
                *(chunksize_t*)prev_size |= NEXT_ALLOCATED_BIT;
            }

            return (void*)f_addr;
        }
    }
}

void
npr_heap_free(struct npr_heap *h, void *p, size_t n)
{
    size_t alloc_size = calc_alloc_size(n);
    uintptr_t free_addr;
    struct npr_free_chunk_header *cur;
    uintptr_t next_addr;
    uintptr_t size_addr;
    chunksize_t *csz;


    if (alloc_size >= (h->page_size-128)) {
        /* large mem */
        uintptr_t free_addr = (uintptr_t)p;
        struct npr_heap_page_header *ph = (struct npr_heap_page_header*)(free_addr & h->pfn_mask);
        del_page(h, ph);
        alloc_size = ALIGN_UP(alloc_size+SIZEOF_PAGE_HEADER, h->page_size);
        free_page((void*)(free_addr-SIZEOF_PAGE_HEADER), alloc_size);
        return;
    }
    free_addr = (uintptr_t)p;
    cur = (struct npr_free_chunk_header*)p;
    next_addr = (free_addr + alloc_size);
    size_addr = next_addr-sizeof(chunksize_t);
    csz = (chunksize_t*)size_addr;

    if ((next_addr & ~(h->pfn_mask)) == 0) {
        /* page last chunk */
        cur->sz = alloc_size;
    } else {
        assert(alloc_size == (*csz & SIZE_MASK));
        if (! (*csz & NEXT_ALLOCATED_BIT)) {
            /* merge next chunk */
            struct npr_free_chunk_header *next = (struct npr_free_chunk_header*)next_addr;
            size_t merge_sz;
            uintptr_t next_size_addr;
            chunksize_t *nsz;
            del_free(next);
            merge_sz = alloc_size + next->sz;
            next_size_addr = (next_addr + next->sz - sizeof(chunksize_t));
            cur->sz = merge_sz;
            nsz = (chunksize_t*)next_size_addr;
            *nsz = merge_sz | NEXT_ALLOCATED_BIT;
        } else {
            cur->sz = alloc_size;
            *csz = alloc_size | NEXT_ALLOCATED_BIT;
        }
    }

    /* printf("free : %p %x\n", free_addr, h->pfn_mask); */
    if ((free_addr & ~(h->pfn_mask)) == SIZEOF_PAGE_HEADER) {
        /* page first chunk */
        /* todo : free page */
        add_free(h, cur);
    } else {
        /* merge prev */
        uintptr_t prev_size_addr = free_addr - sizeof(chunksize_t);
        chunksize_t prev_size = *(chunksize_t*)prev_size_addr;
        if (prev_size & THIS_ALLOCATED_BIT) {
            add_free(h, cur);
            *(chunksize_t*)prev_size_addr &= ~NEXT_ALLOCATED_BIT;
        } else {
            /* merge prev chunk */
            uintptr_t prev_addr = (free_addr - (prev_size&SIZE_MASK));
            struct npr_free_chunk_header *prev = (struct npr_free_chunk_header*)prev_addr;
            size_t merge_sz = cur->sz + prev->sz;
            uintptr_t size_addr;
            prev->sz = merge_sz;
            size_addr = prev_addr + merge_sz - sizeof(chunksize_t);
            *(chunksize_t*)size_addr = merge_sz | NEXT_ALLOCATED_BIT;
        }
    }
}

void
npr_heap_dump(FILE *fp, struct npr_heap *h)
{
    struct npr_heap_page_header *ph;
    struct npr_free_chunk_header *fch;

    fprintf(fp, "=== dump heap %p ==\n", h);
    fprintf(fp, "  === pages ===\n");
    ph = h->page_head.next;
    while (ph != &h->page_tail) {
        fprintf(fp, "  page %p : flags=%08x, next=%p, prev=%p\n",
                ph, ph->flags,
                ph->next,
                ph->prev);
        assert(ph->next->prev == ph);
        ph = ph->next;
    }

    fprintf(fp, "  === free chain ===\n");
    fch = h->free_head.next;

    while (fch != &h->free_tail) {
        fprintf(fp, "  free %p : size=%5d, next=%p, prev=%p\n",
                fch,
                (int)fch->sz,
                fch->next, fch->prev);
        assert(fch->next->prev == fch);
        fch = fch->next;
    }
    
}

void
npr_heap_fini(struct npr_heap *h)
{
    struct npr_heap_page_header *ph = h->page_head.next;
    while (ph != &h->page_tail) {
        struct npr_heap_page_header *next = ph->next;
        if (ph->flags & NPR_LARGE_CHUNK) {
            free_page(ph, ph->page_size);
        } else {
            free_page(ph, h->page_size);
        }
        ph = next;
    }
}

