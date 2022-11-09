from collections.abc import Iterator
from typing import Self, TypeVar, Generic


def _node_link_prop (*, attr=None, getv=None, setv=None):
	if attr is not None:
		if getv is None:
			getv = lambda self: getattr(self, attr)
		if setv is None:
			setv = lambda self, v: setattr(self, attr, v)
	return property(getv, setv)


def _attr_or_call (thing, default, attr):
	if thing is None:
		return default
	if isinstance(thing, str):
		return attr
	elif callable(thing):
		return thing


class Node:
	__slots__ = '_prev_node', '_next_node'
	def __init__ (self):
		self._prev_node = self._next_node = self

	prev: Self = _node_link_prop(attr=__slots__[0])
	next: Self = _node_link_prop(attr=__slots__[1])

	left = prev
	'Alias for `.prev`, when the idea of a "Previous" and "Next" node don\'t make much sense'
	right = next
	'Same as `.left`, but for `.next`'

	start = next
	'An alias for `.next` to make DxDLists made of only a "head" node more readable'
	end   = prev
	'same as `.start`, but for `.prev` instead of `.next`'

	def insert_before (self, node: Self) -> Self:
		'''
		...prev - +node+ - self - next...
		'''
		node.next = self
		node.prev = self.prev
		self.prev.next = node
		self.prev = node
		return node

	def insert_after (self, node: Self) -> Self:
		'''
		...prev - self - +node+ - next...
		'''
		node.prev = self
		node.next = self.next
		self.next.prev = node
		self.next = node
		return node

	def remove_from_chain (self, clear_self = False):
		'''
		(...prev - self - next...) becomes (...prev - next...)
		'''
		self.prev.next = self.next
		self.next.prev = self.prev
		if clear_self:
			self.next = self.prev = self

	def __iter__ (self) -> Iterator[Self]:
		n = self
		while True:
			yield n
			if (n := n.next) is None or n is self:
				break

	def __reversed__ (self) -> Iterator[Self]:
		n = self
		while True:
			yield n
			if (n := n.prev) is None or n is self:
				break

	def walk (head, key=None, start=None, sentinel=None):
		'''
		if start is passed the same obj as head, or True,
		chain will skip head as the first node
		'''
		nextv = _attr_or_call(
			key,
			lambda self: getattr(self, 'next'),
			lambda self: getattr(self, key),
		)
		
		if (isinstance(start, bool) and start) or (start is head):
			start = nextv(head)
		elif start is None:
			start = head
		
		getend = _attr_or_call(
			sentinel,
			lambda self: (self is None) or (self is head),
			lambda self: (self is None) or (self is getattr(head, sentinel))
		)

		n = start
		while True:
			yield n
			if getend(n:=nextv(n)):
				break

	def chain_length (self) -> int:
		i = 0
		n = self
		while True:
			i += 1
			if (n:=n.next) is None or (n is self):
				return i

	def __lshift__ (self, other: Self) -> Self:
		"""
		`(self >> other) == self`

		The new chain is `... -> prev -> self -> other -> ...`
		"""
		self.prev = other.prev
		self.next = other
		other.prev.next = self
		other.prev      = self
		return self

	def __rshift__ (self, other: Self) -> Self:
		"""
		`(self << other) == other`

		The new chain is `... <- self <- other <- next <- ...`
		"""
		other.prev = self
		other.next = self.next
		self.next.prev = other
		self.next      = other
		return other

	def __pos__ (self) -> Iterator[Self]:
		"""
		A forwards-marching head iterator, IE in the chain:
		`a - b - c`, `for n in +a:` will yield `b, c`, but not `a`
		"""
		n = self
		while not (((n:=n.next) is None) or (n is self)):
			yield n

	def __neg__ (self) -> Iterator[Self]:
		"""
		The same as `+a`, except in reverse; in the chain `a - b - c` you
		will get `c, b` but not `a`
		"""
		n = self
		while not (((n:=n.prev) is None) or (n is self)):
			yield n

	def walk_head_safe (self) -> Iterator[tuple[Self, Self]]:
		'Yields the nodes `node` and `node.next`'
		n = self.next
		nx: Self
		while n is not self:
			nx = n.next
			yield n, nx
			n = nx



NT = TypeVar('NT', bound=Node)
class DxDList(Generic[NT]):
	__slots__ = '_head'
	def __init__ (self, chain: NT = None):
		self._head = Node()
		head = self._head
		if chain is not None:
			head.next = chain
			head.prev = chain.prev
			chain.prev.next = head
			chain.prev = head

	@property
	def head (self) -> NT:
		return self._head
	
	@head.setter
	def head (self, new: NT):
		self._head = new
	
	@property
	def start (self) -> NT:
		return self._head.next

	@start.setter
	def start (self, node: NT):
		h = self._head
		node.next = h.next.next
		node.prev = h
		node.next.prev = node
		h.next = node

	@property
	def end (self) -> NT:
		return self._head.prev

	@end.setter
	def end (self, node: NT):
		h = self._head
		node.prev = h.prev.prev
		node.next = h
		node.prev.next = node
		h.next = node

	def __iter__ (self) -> Iterator[NT]:
		n = start = self._head
		while not ((n:=n.next) is None or n is start):
			yield n

	def __reversed__ (self) -> Iterator[NT]:
		n = end = self._head
		while not ((n:=n.prev) is None or n is end):
			yield n

	def __delitem__ (self, node: NT):
		node.remove_from_chain()

