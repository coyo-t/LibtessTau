from typing import TypeVar, Generic
DEFAULT_SIZE = 32

T = TypeVar('T')
S = TypeVar('S')
class Heap(Generic[T]):
	class HandleElem(Generic[S]):
		__slots__ = 'node', 'item'
		def __init__ (self, item: S = None, node: int = 0):
			self.item = item
			self.node = node
	
	def __init__ (self, size = DEFAULT_SIZE, cmp = None):
		if cmp is None:
			cmp = lambda a, b: a <= b
		
		self.size = 0
		self.size_max = size
		self.free_list = 0
		self.initialized = False
		self.cmp = cmp

		self.handles = self._realloc_handles(size, setatr=False)
		self.nodes   = self._realloc_nodes(size, setatr=False)
		self.nodes[1] = 1

	def _thethree (self):
		return self.nodes, self.handles, self.cmp

	def _realloc_nodes (self, size: int = -1, setatr = True):
		size = (self.size_max if size == -1 else size) + 1
		newn = [-1] * (size + 1)
		if setatr:
			newn[:len(self.nodes)] = self.nodes
		self.nodes = newn
		return newn

	def _realloc_handles (self, size: int = -1, setatr = True) -> list[HandleElem[T]]:
		size = (self.size_max if size == -1 else size) + 1
		newn = [self.HandleElem[T]() for _ in range(size + 1)]
		if setatr:
			newn[:len(self.handles)] = self.handles
		self.handles = newn
		return newn

	def float_down (ds, cur: int):
		n, h, cmp = ds._thethree()
		h_cur = n[cur]
		while True:
			child = cur << 1
			if child < ds.size and cmp(h[n[child+1]].item, h[n[child]].item):
				child += 1
			assert child <= ds.size_max
			h_child = n[child]
			if child > ds.size or cmp(h[h_cur].item, h[h_child].item):
				n[cur] = h_cur
				h[h_cur].node = cur
				break
			n[cur] = h_child
			h[h_child].node = cur
			cur = child

	def float_up (ds, cur: int):
		n, h, cmp = ds._thethree()
		h_cur = n[cur]
		while True:
			parent = cur >> 1
			h_parent = n[parent]
			if parent == 0 or cmp(h[h_parent].item, h[h_cur].item):
				n[cur] = h_cur
				h[h_cur].node = cur
				break
			n[cur] = h_parent
			h[h_parent].node = cur
			cur = parent

	def init (ds):
		for i in range(ds.size, 0, -1):
			ds.float_down(i)
		ds.initialized = True

	def insert (ds, new_item: T):
		ds.size += 1
		cur = ds.size
		if (cur << 1) > ds.size_max:
			ds.size_max <<= 1
			ds._realloc_handles()
			ds._realloc_nodes()
		if ds.free_list == 0:
			free = cur
		else:
			free = ds.free_list
			ds.free_list = ds.handles[free].node
		ds.nodes[cur] = free
		hh = ds.handles[free]
		hh.node = cur
		hh.item = new_item
		if ds.initialized:
			ds.float_up(cur)
		return free

	def extract_min (ds) -> T:
		n, h, _ = ds._thethree()
		h_min = n[1]
		min_v = h[h_min].item
		if ds.size > 0:
			n[1] = n[ds.size]
			h[n[1]].node = 1
			h[h_min].item = None
			h[h_min].node = ds.free_list
			ds.free_list = h_min
			ds.size -= 1
			if ds.size > 0:
				ds.float_down(1)
		return min_v

	def delete (ds, h_cur: int):
		n, h, cmp = ds._thethree()
		assert (1 <= h_cur <= ds.size_max) and (h[h_cur].item is not None)
		cur = h[h_cur].node
		n[cur] = n[ds.size]
		h[n[cur]].node = cur
		ds.size -= 1
		if cur <= ds.size:
			if cur <= 1 or cmp(h[n[cur>>1]].item, h[n[cur]].item):
				ds.float_down(cur)
			else:
				ds.float_up(cur)
		h[h_cur].item = None
		h[h_cur].node = ds.free_list
		ds.free_list = h_cur

	def is_empty (ds):
		return ds.size == 0

	def min (ds) -> T:
		return ds.handles[ds.nodes[1]].item

	def __len__ (ds):
		return ds.size
	
	def __delitem__ (ds, h_cur: int):
		ds.delete(h_cur)
