#include <assert.h>
#include "pcgpmem.h"

void arenaInit(Arena *a, void *mem, size_t size) {
	a->start = mem;
	a->end = a->start + size;
	a->ptr = a->start;
}

void *arenaAlloc(Arena *a, size_t size) {
	uintptr_t ptr = (uintptr_t)a->ptr;
	uintptr_t aptr = (uintptr_t)ptr;
	uintptr_t mod = ptr & (PCGPMEM_ALIGN - 1);
	if (!mod) {
		aptr += PCGPMEM_ALIGN - mod;
	}
	ptr = aptr + size;
	assert(ptr <= (uintptr_t)a->end);
	a->ptr = (void *)ptr;
	return (void *)aptr;
}

void arenaFree(Arena *a) {
	a->ptr = a->start;
}
