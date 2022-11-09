from math import sqrt
from typing import Self


def _safe_div (a, b):
	if b == 0.0:
		return 0.0
	return a / b

def _unzip_input (*args):
	if len(args) == 1 and isinstance(args[0], tuple|list):
		args = args[0]
	for value in args:
		if hasattr(value, '__iter__'):
			for v in value:
				yield v
			continue
		yield value

def _co2_leq (a0, a1, b0, b1) -> bool:
	return (a0 < b0) or (a0 == b0 and a1 <= b1)

class CoordBase:
	__slots__ = ()
	__alias__ = ()
	__name__  = ''

	@classmethod
	@property
	def ZERO (cls):
		return cls(0.0)
	
	@classmethod
	@property
	def ONE (cls):
		return cls(1.0)

	def __init__ (self, *args):
		self.set_co(0.0)
		if len(args) == 1:
			self.set_co(args[0])
			args = args[0]
		else:
			self.set_co(*args)

	def _get_i (self, index: int) -> float:
		return getattr(self, self.__slots__[index])
	def _set_i (self, index: int, value: float):
		setattr(self, self.__slots__[index], value)
	def _get_s (self, index: str) -> float:
		return getattr(self, index)
	def _set_s (self, index: str, value: float):
		setattr(self, index, value)

	def _swizzle (self, key: str):
		slots = self.__slots__
		alias = self.__alias__
		_alias = None
		for ch in key:
			if ch == '_':
				yield None
				continue
			if _alias is None:
				if ch in slots:
					_alias = slots
				else:
					for a in alias:
						if ch in a:
							_alias = a
							break
				if _alias is None:
					raise KeyError(f"Couldn't determine swizzle alias for key {key}")
			if not ch in _alias:
				KeyError(f"Can't mix swizzle sets for {key} in {_alias}")
			yield slots[_alias.index(ch)]

	def _index_list (self, indicies: tuple):
		for k in indicies:
			if (k is ...) or (k is None):
				yield None
				continue
			elif isinstance(k, str):
				yield k
			elif isinstance(k, int):
				yield self.__slots__[k]

	def __iter__ (self):
		return (self._get_s(k) for k in self.__slots__)

	def __getitem__ (self, index):
		if isinstance(index, int):
			return self._get_i(index)
		elif index is ...:
			return self.__class__(*self)
		
		if isinstance(index, str):
			if len(index) == 0:
				return self.__class__(*self)
			genr, count = self._swizzle(index), len(index)
			out = tuple[float](0.0 if k is None else self._get_s(k) for k in genr)
			if count == 1:
				return out[0]
			elif count == 2:
				return Co2(*out)
			elif count == 3:
				return Co3(*out)
			return out
		elif hasattr(index, '__iter__'):
			genr = self._index_list(index)
			out = tuple[float](0.0 if k is None else self._get_s(k) for k in genr)
			if len(out) == 1:
				return out[0]
			return out

	def __setitem__ (self, index, value):
		slr = lambda: range(len(self.__slots__))
		is_dead = lambda x: (x is None) or (x is ...)
		if hasattr(value, '__iter__'):
			v_genr = _unzip_input(value)
		else:
			v = 0.0 if is_dead(value) else value
			v_genr = (v for _ in slr())
		
		if isinstance(index, int):
			self._set_i(index, next(v_genr))
			return
		elif index is ...:
			self.set_co(v_genr)
			return
		
		if isinstance(index, str):
			if len(index) == 0:
				self.set_co(v_genr)
				return
			genr = self._swizzle(index)
		elif hasattr(index, '__iter__'):
			genr = self._index_list(index)
		for i, v, _, in zip(genr, v_genr, range(len(self.__slots__))):
			if i is None:
				continue
			self._set_s(i, v)

	def __repr__ (self):
		return f"{self.__name__}({', '.join(str(v) for v in self)})"

	def set_co (self, /, *args):
		if len(args) == 1:
			args = args[0]
			if not hasattr(args, '__iter__'):
				if args is None:
					args = 0.0
				for k in self.__slots__:
					self._set_s(k, args)
				return
		for k, v in zip(self.__slots__, _unzip_input(args)):
			if k is None:
				continue
			self._set_s(k, 0.0 if v is None else v)

	def length_squared (self) -> float:
		return sum(k ** 2 for k in self)

	def length (self) -> float:
		return sqrt(self.length_squared())

	def _inv_mag (self) -> float:
		mag = self.length_squared()
		if mag == 0.0:
			return 0.0
		return 1.0 / sqrt(mag)

	def normalize (self):
		'Normalizes the coord in-place'
		mag = self._inv_mag()
		self.set_co(k * mag for k in self)

	def normalized (self):
		mag = self._inv_mag()
		return self.__class__(k * mag for k in self)

	def dist_l1 (self, other) -> float:
		return sum(abs(i-j) for i, j in zip(self, other))

	def dot (self, other) -> float:
		return sum(a * b for a, b in zip(self, other))

	def __neg__ (self):
		return self.__class__(-v for v in self)

	def _op_uhh (self, thing):
		if hasattr(thing, '__iter__'):
			return zip(self, _unzip_input(thing))
		else:
			return zip(self, (thing for _ in range(len(self.__slots__))))

	def __eq__ (self, other) -> bool:
		other = self._op_uhh(other)
		return all(a == b for a, b in other)

	def pow (self, other) -> Self:
		return self.__class__(a ** b for a, b in self._op_uhh(other))

	def __add__ (self, other) -> Self:
		return self.__class__(a + b for a, b in self._op_uhh(other))
	def __radd__ (self, other) -> Self:
		return self.__class__(b + a for a, b in self._op_uhh(other))

	def __sub__ (self, other) -> Self:
		return self.__class__(a - b for a, b in self._op_uhh(other))
	def __rsub__ (self, other) -> Self:
		return self.__class__(b - a for a, b in self._op_uhh(other))
	
	def __mul__ (self, other) -> Self:
		return self.__class__(a * b for a, b in self._op_uhh(other))
	def __rmul__ (self, other) -> Self:
		return self.__class__(b * a for a, b in self._op_uhh(other))
	
	def __truediv__ (self, other) -> Self:
		return self.__class__(_safe_div(a, b) for a, b in self._op_uhh(other))
	def __rtruediv__ (self, other) -> Self:
		return self.__class__(_safe_div(b, a) for a, b in self._op_uhh(other))


class Co2(CoordBase):
	__slots__ = 'x', 'y'
	__alias__ = 'uv', 'st', 'ij'
	__name__ = 'Co2'
	def __init__ (self, /, *args, x = None, y = None):
		if len(args) == 0:
			args = (x, y)
		super().__init__(*args)

	def __le__ (a, b: Self) -> bool:
		return _co2_leq(a.x, a.y, b.x, b.y)

	def __gt__ (a, b: Self) -> bool:
		return not _co2_leq(a.x, a.y, b.x, b.y)

	@property
	def xy (self) -> Self:
		'Alias for `co[...]`'
		return Co2(self.x, self.y)

	@xy.setter
	def xy (self, v):
		self[...] = v

	@property
	def yx (self) -> Self:
		'alias for `self` transposed'
		return Co2(self.y, self.x)


class Co3(CoordBase):
	__slots__ = 'x', 'y', 'z'
	__alias__ = 'stp', 'uvw', 'ijk', 'rgb'
	__name__ = 'Co3'
	def __init__ (self, /, *args, x = None, y = None, z = None):
		if len(args) == 0:
			args = (x, y, z)
		super().__init__(*args)

	def cross (a, b: Self):
		return Co3(a.y*b.z - a.z*b.y,
		           a.z*b.x - a.x*b.z,
		           a.x*b.y - a.y*b.x)

	@property
	def xyz (self) -> Self:
		'An alias for `co[...]`. Just for readability I suppose.'
		return Co3(self.x, self.y, self.z)

	@xyz.setter
	def xyz (self, v):
		self[...] = v

class Co4(CoordBase):
	__slots__ = 'x', 'y', 'z', 'w'
	__alias__ = 'rgba', 'stpq', 'ijkl'
	__name__ = 'Co4'
	def __init__ (self, /, *args, r = None, g = None, b = None, a = None):
		if len(args) == 0:
			args = (r, g, b, a)
		super().__init__(*args)


# def __main ():
# 	a = Co2(10, 20)
# 	x, y = a
# 	print(x, y)
	

# if __name__ == '__main__':
# 	__main()