from pathlib import Path

from tess import *

def read_curve_file (path):
	splines = []
	with open(path, 'r') as f:
		linez = f.readlines()
		for l in linez:
			acc = ''
			spline = []
			pt = []
			for ch in l:
				if ch == ',' and len(acc) != 0:
					pt.append(float(acc))
					if len(pt) == 2:
						spline.append(tuple(pt))
						pt.clear()
					acc = ''
				else:
					acc += ch
			splines.append(tuple(spline))
	return tuple(splines)


OUTP = Path('D:/_projects/_desktop/libtess2 FINAL DRAFT AAAAAAAAAAAAAAAAAAAA/blen tes')
curve = read_curve_file(Path(OUTP, 'bltest.txt'))

t = Tesselator()
t.normal = Co3(0.0, 0.0, 1.0)
t.process_cdt = True

polygon_to_test = curve
for c in polygon_to_test:
	t.add_contour(c)

t.tesselate(
	WindingRule.ODD,
	ElementType.POLYGONS,
	poly_size=3,
	vert_size=2
)

verts, elems = t.verticies, t.elements

out = list[str]()

LINE_V = 'v {:f} {:f} {:f}'
LINE_F = 'f {:d} {:d} {:d}'

out.append('o Haha Wow!!!')

vgen = iter(verts)
while (x:=next(vgen, None)) is not None:
	y = next(vgen)
	out.append(LINE_V.format(x, y, 0.0))

out.append('s 0')

fgen = iter(elems)
while (v1:=next(fgen, None)) is not None:
	v2 = next(fgen)
	v3 = next(fgen)

	out.append(LINE_F.format(v1+1, v2+1, v3+1))

with open(Path(OUTP, 'thingy.obj'), 'w') as f:
	f.write('\n'.join(out))
