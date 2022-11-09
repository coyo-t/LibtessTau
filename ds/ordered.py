from typing import TypeVar, Generic

T = TypeVar('T')
S = TypeVar('S')

class Ordered(Generic[T]):
	__slots__ = 'head', 'cmp'
	class Node(Generic[S]):
		__slots__ = 'prev', 'next', 'item'
		def __init__ (self, item: S = None):
			self.prev = self.next = self
			self.item: S = item
		
		def delete (self):
			self.prev.next = self.next
			self.next.prev = self.prev
			self.next = self.prev = self

	def _new_node (self, itm = None):
		return self.Node[T](itm)

	def __init__ (self, cmp = None):
		if cmp is None:
			cmp = lambda a, b: a <= b
		self.cmp = cmp
		self.head = self._new_node()

	def insert_before (ds, node: Node[T], item: T) -> Node[T]:
		while True:
			if ((node:=node.prev).item is None) or ds.cmp(node.item, item):
				break
		new = ds._new_node(item)
		new.next = node.next
		new.prev = node
		node.next.prev = new
		node.next = new
		return new

	def delete (ds, node: Node[T]):
		node.next.prev = node.prev
		node.prev.next = node.next
		del node

	def search (ds, item: T) -> Node[T]:
		node = ds.head
		while True:
			if ((node:=node.next).item is None) or ds.cmp(item, node.item):
				break
		return node

	def insert (ds, item: T):
		return ds.insert_before(ds.head, item)

	def to_list (ds):
		outl = list[T]()
		start = node = ds.head
		while (node:=node.next) is not start:
			outl.append(node.item)
		return outl

	def size (ds):
		outc = 0
		start = node = ds.head
		while (node:=node.next) is not start:
			outc += 1
		return outc

	def min (ds) -> Node[T]:
		return ds.head.next

	def max (ds) -> Node[T]:
		return ds.head.prev


# def __main ():
# 	f = Ordered[str]()
# 	f.insert('haha')
# 	f.insert('wow!!!')

# 	l = f.to_list()
# 	print(l)


# 	s = f.search('haha').item
# 	print(s)


# if __name__ == '__main__':
# 	__main()
