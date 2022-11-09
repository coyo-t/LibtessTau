from typing import Generic, TypeVar
from .heap import Heap

DEFAULT_SIZE = 32

T = TypeVar('T')
class Priority(Generic[T]):
	def __init__ (self, cmp = None, size = DEFAULT_SIZE):
		if cmp is None:
			cmp = lambda a, b: a <= b
		
		self.heap = Heap[T](size=size, cmp=cmp)
		self.items = self._realloc_items(size, setatr=False)
		self.order = list[int]()
		self.size = 0
		self.size_max = size
		self.initialized = False
		self.cmp = cmp

	def _realloc_items (self, size: int = -1, setatr = True):
		size = self.size_max if size == -1 else size
		newi: list[T] = [None] * size
		if setatr:
			newi[:len(self.items)] = self.items
		self.items = newi
		return newi

	def __getitem__ (self, order_index: int):
		return self.items[self.order[order_index]]

	def init (ds, *, debug = False, insr_sort_limit = 10):
		cmp = ds.cmp
		order = ds.order = [i for i in range(ds.size)]
		items = ds.items
		# perform an iterative impl of quicksort over recursive
		# technically, theres no reason it cant be recursive, but
		# libtess went with iterative over recursive and i agree
		# with that decision so it stays
		seed: int = 0x7830F0C3
		stack: list[tuple[int, int]] = [(0, ds.size - 1)]
		while len(stack) > 0:
			p, r = stack.pop()
			while r > (p + insr_sort_limit):
				# libtess used a u32 for seed, so im adding the bitmask
				# here to retain behavior even if it doesn't really matter.
				seed = (seed * 0x5BC19F0D + 1) & 0xFFFFFFFF
				i = p + seed % (r - p + 1)
				piv      = order[i]
				order[i] = order[p]
				order[p] = piv
				i = p - 1
				j = r + 1
				while True:
					while True:
						i += 1
						if cmp(ds[i], items[piv]):
							break
					while True:
						j -= 1
						if cmp(items[piv], ds[j]):
							break
					if not i < j:
						break
					else:
						order[i], order[j] = order[j], order[i]
				if (i - p) < (r - j):
					top = (j + 1, r)
					r = i - 1
				else:
					top = (p, i - 1)
					p = j + 1
				stack.append(top)
			# depending on how large you consider a "small" list to be based
			# on INSR_SORT_LIMIT, perform an insertion sort on
			# the partition over quicksort
			for i in range(p+1, r+1):
				piv = order[i]
				j = i
				while j > p and not cmp(items[piv], ds[j-1]):
					order[j] = order[j-1]
					j -= 1
				order[j] = piv
		ds.size_max = ds.size
		ds.initialized = True
		ds.heap.init()
		if debug:
			for i in range(ds.size - 1):
				assert cmp(ds[i+1], ds[i])

	def insert (ds, new_item: T):
		if ds.initialized:
			return ds.heap.insert(new_item)
		cur = ds.size
		ds.size += 1
		if ds.size >= ds.size_max:
			ds.size_max <<= 1
			ds._realloc_items()
		ds.items[cur] = new_item
		return -(cur + 1)

	def extract_min (ds) -> T:
		if ds.size == 0:
			return ds.heap.extract_min()
		sort_min = ds[ds.size-1]
		if not ds.heap.is_empty():
			heap_min = ds.heap.min()
			if ds.cmp(heap_min, sort_min):
				return ds.heap.extract_min()
		while True:
			ds.size -= 1
			if ds.size <= 0 or (ds[ds.size-1] is not None):
				break
		return sort_min

	@property
	def is_empty (ds):
		return ds.size == 0 and ds.heap.is_empty()

	def delete (ds, cur: int):
		if cur >= 0:
			del ds.heap[cur]
			return
		cur = -(cur + 1)
		assert cur < ds.size_max and (ds.items[cur] is not None)
		ds.items[cur] = None
		while ds.size > 0 and (ds[ds.size-1] is None):
			ds.size -= 1

	def min (ds) -> T:
		if ds.size == 0:
			return ds.heap.min()
		sort_min = ds[ds.size-1]
		if not ds.heap.is_empty():
			heap_min = ds.heap.min()
			if ds.cmp(heap_min, sort_min):
				return heap_min
		return sort_min

	def __len__ (ds):
		return ds.size

	def __delitem__ (ds, cur: int):
		ds.delete(ds, cur)

	def __iter__ (ds):
		for i, _ in zip(ds.order, range(len(ds.size))):
			yield ds.items[i]
