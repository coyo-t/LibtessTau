from ds.coord import Co2, Co3

def _eval_sign (
	u: Co2, v: Co2, w: Co2,
	transposed = False,
	*,
	eval = False,
	sign = False
) -> float:
	if transposed:
		u = u['yx']
		v = v['yx']
		w = w['yx']

	assert u <= v <= w

	gap_l = v.x - u.x
	gap_r = w.x - v.x

	if gap_l + gap_r > 0:
		if eval:
			if gap_l < gap_r:
				return (v.y - u.y) + (u.y - w.y) * (gap_l / (gap_l + gap_r))
			else:
				return (v.y - w.y) + (w.y - u.y) * (gap_r / (gap_l + gap_r))
		if sign:
			return (v.y - w.y) * gap_l + (v.y - u.y) * gap_r
	return 0.0 # Vertical line

def edge_eval (u: Co2, v: Co2, w: Co2):
	return _eval_sign(u,v,w, eval=True)
def edge_sign (u: Co2, v: Co2, w: Co2):
	return _eval_sign(u,v,w, sign=True)

def trans_eval (u: Co2, v: Co2, w: Co2):
	return _eval_sign(u,v,w,True, eval=True)
def trans_sign (u: Co2, v: Co2, w: Co2):
	return _eval_sign(u,v,w,True, sign=True)


def edge_intersect (o1: Co2, d1: Co2,
                    o2: Co2, d2: Co2,
                    v: Co2 = None) -> Co2:
	def interpolate (a: float, x: float, b: float, y: float):
		a = 0.0 if a < 0 else a
		b = 0.0 if b < 0 else b
		
		if a <= b:
			if b == 0:
				return (x + y) / 2
			else:
				return x + (y-x) * (a/(a+b))
		else:
			return y + (x-y) * (b/(a+b))

	if v is None:
		v = Co2.ZERO

	# This is certainly not the most efficient way to find the intersection
	# of two line segments, but it is very numerically stable.
	#
	# Strategy: find the two middle vertices in the VertLeq ordering,
	# and interpolate the intersection s-value from these.  Then repeat
	# using the TransLeq ordering to find the intersection t-value.
	if o1 > d1:
		o1, d1 = d1, o1

	if o2 > d2:
		o2, d2 = d2, o2

	if o1 > o2:
		o1, o2 = o2, o1
		d1, d2 = d2, d1

	if o2 > d1:
		# Technically, no intersection -- do our best
		v.x = (o2.x + d1.x) / 2
	elif d1 <= d2:
		# Interpolate between o2 and d1 */
		z1 = edge_eval(o1, o2, d1)
		z2 = edge_eval(o2, d1, d2)
		if (z1 + z2) < 0:
			z1 = -z1
			z2 = -z2
		v.x = interpolate(z1, o2.x, z2, d1.x)
	else:
		z1 = edge_sign(o1, o2, d1)
		z2 = -edge_sign(o1, d2, d1)
		if (z1 + z2) < 0:
			z1 = -z1
			z2 = -z2
		v.x = interpolate(z1, o2.x, z2, d2.x)

	# Repeat for .y
	if o1.yx > d1.yx:
		o1, d1 = d1, o1

	if o2.yx > d2.yx:
		o2, d2 = d2, o2

	if o1.yx > o2.yx:
		o1, o2 = o2, o1
		d1, d2 = d2, d1

	if o2.yx > d1.yx:
		v.y = (o2.y + d1.y) / 2
	elif d1.yx <= d2.yx:
		z1 = trans_eval(o1, o2, d1)
		z2 = trans_eval(o2, d1, d2)
		if (z1 + z2) < 0:
			z1 = -z1
			z2 = -z2
		v.y = interpolate(z1, o2.y, z2, d1.y)
	else:
		z1 = trans_sign(o1, o2, d1)
		z2 = -trans_sign(o1, d2, d1)
		if (z1 + z2) < 0:
			z1 = -z1
			z2 = -z2
		v.y = interpolate(z1, o2.y, z2, d2.y)

	return v

def vert_ccw (u: Co2, v: Co2, w: Co2) -> bool:
	return ((u.x * (v.y - w.y) +
	         v.x * (w.y - u.y) +
	         w.x * (u.y - v.y)) >= 0)

def short_axis (co: Co3):
	i = 0
	if abs(co.y) < abs(co.x):
		i = 1
	if abs(co.z) < abs(co[i]):
		i = 2
	return i, co[i]

def long_axis (co: Co3):
	i = 0
	if abs(co.y) > abs(co.x):
		i = 1
	if abs(co.z) > abs(co[i]):
		i = 2
	return i, co[i]

def in_circle (v: Co2, v0: Co2, v1: Co2, v2: Co2) -> float:
	adx = v0.x - v.x
	ady = v0.y - v.y
	bdx = v1.x - v.x
	bdy = v1.y - v.y
	cdx = v2.x - v.x
	cdy = v2.y - v.y

	abdet = adx * bdy - bdx * ady
	bcdet = bdx * cdy - cdx * bdy
	cadet = cdx * ady - adx * cdy

	alift = adx * adx + ady * ady
	blift = bdx * bdx + bdy * bdy
	clift = cdx * cdx + cdy * cdy

	return alift * bcdet + blift * cadet + clift * abdet


