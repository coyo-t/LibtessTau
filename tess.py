from typing import Self, LiteralString
from enum import Enum, auto as iota

from ds.coord import Co2, Co3
from ds.node import Node
from ds.ordered import Ordered
from ds.priority import Priority
import geom


class ElementType(Enum):
	POLYGONS           = iota()
	CONNECTED_POLYGONS = iota()
	BOUNDARY_CONTOURS  = iota()


class WindingRule(Enum):
	ODD         = iota(), lambda n: n & 1
	NONZERO     = iota(), lambda n: n != 0
	POS         = iota(), lambda n: n > 0
	NEG         = iota(), lambda n: n < 0
	ABS_GEQ_TWO = iota(), lambda n: (n >= 2) or (n <= -2)

	def is_inside (self, n: int) -> bool:
		return self.value[1](n)


class Vertex(Node):
	def __init__ (self):
		super().__init__()
		self.prev: Self
		self.next: Self
		self.edge: Edge = None
		
		# internal
		self.src_co = Co3()
		self.co = Co2()

		self.pq_handle: int = -1
		self.n:         int = -1
		self.index:     int = -1

	@classmethod
	def Make (cls, e_origin: 'Edge', v_next: Self, v_new: Self = None):
		if v_new is None:
			v_new = cls()

		v_prev = v_next.prev
		v_new.prev = v_prev
		v_prev.next = v_new
		v_new.next = v_next
		v_next.prev = v_new

		v_new.edge = e_origin

		e = e_origin
		while True:
			e.org = v_new
			if (e:=e.onext) is e_origin:
				break

		return v_new

	def kill (self, new_origin: Self):
		e = eStart = self.edge
		while True:
			e.org = new_origin
			if (e:=e.onext) is eStart:
				break

		self.remove_from_chain()

	def weights (self: Self, org: Self, dst: Self):
		t1 = org.co.dist_l1(self.co)
		t2 = dst.co.dist_l1(self.co)
		w0 = 0.5 * t2 / (t1 + t2)
		w1 = 0.5 * t1 / (t1 + t2)

		self.src_co.xyz += w0*org.src_co + w1*dst.src_co

		# isect['xyz'] += w0 * org + w1 * dst

		# isect.x += w0 * org.x + w1 * dst.x
		# isect.y += w0 * org.y + w1 * dst.y
		# isect.z += w0 * org.z + w1 * dst.z

	def get_intersect_data (self: Self, orgUp: Self, dstUp: Self, orgLo: Self, dstLo: Self):
		self.src_co.xyz = 0.0
		self.index = -1
		self.weights(orgUp, dstUp)
		self.weights(orgLo, dstLo)

	def __eq__ (self, other: Self):
		return self.co == other.co

	def __le__ (self, other: Self):
		return self.co <= other.co


class Face(Node):
	def __init__ (self):
		super().__init__()
		self.prev: Self
		self.next: Self
		self.edge: Edge = None

		# internal
		self.trail: Self = None
		self.n: int = -1
		
		self.marked = False
		self.inside = False

	@classmethod
	def Make (cls, e_orig: 'Edge', f_next: Self, f_new: Self = None):
		if f_new is None:
			f_new = cls()

		fPrev: Face = f_next.prev
		f_new.prev = fPrev
		fPrev.next = f_new
		f_new.next = f_next
		f_next.prev = f_new

		f_new.edge = e_orig

		f_new.inside = f_next.inside

		e = e_orig
		while True:
			e.lface = f_new
			if (e:=e.lnext) is e_orig:
				break
		return f_new

	def kill (self, new_lface: Self):
		e = eStart = self.edge
		while True:
			e.lface = new_lface
			if (e:=e.lnext) is eStart:
				break
		
		self.remove_from_chain()

	def count_verts (self) -> int:
		eCur = self.edge
		n = 0
		while True:
			n += 1
			if (eCur:=eCur.lnext) is self.edge:
				break
		return n


class Edge(Node):
	_ptr_counter  = 0
	_pair_counter = 0

	def __init__ (self):
		super().__init__()
		self.next: Self = self
		self.twin:  Self = None
		self.onext: Self = None
		self.lnext: Self = None
		self.org: Vertex = None
		self.lface: Face = None

		# internal
		self.active_region: ActiveRegion = None
		self.winding = 0
		self.mark = False

		self._ptr_index  = 0
		self._pair_index = 0

	prev: Self = property(lambda self: self.twin.next,
	                      lambda self, v: setattr(self.twin, 'next', v))
	
	dst: Vertex = property(lambda self: self.twin.org,
	                       lambda self, v: setattr(self.twin, 'org', v))
	
	rface: Face = property(lambda self: self.twin.lface,
	                       lambda self, v: setattr(self.twin, 'lface', v))

	oprev: Self = property(lambda self: self.twin.lnext,
	                       lambda self, v: setattr(self.twin, 'lnext', v))

	lprev: Self = property(lambda self: self.onext.twin,
	                       lambda self, v: setattr(self.onext, 'twin', v))
	
	dprev: Self = property(lambda self: self.lnext.twin,
	                       lambda self, v: setattr(self.lnext, 'twin', v))
	
	rprev: Self = property(lambda self: self.twin.onext,
	                       lambda self, v: setattr(self.twin, 'onext', v))
	
	dnext: Self = property(lambda self: self.rprev.twin,
	                       lambda self, v: setattr(self.rprev, 'twin', v))
	
	rnext: Self = property(lambda self: self.oprev.twin,
	                       lambda self, v: setattr(self.oprev, 'twin', v))

	@property
	def neighbour_face (self):
		if self.rface is None:
			return -1
		if not self.rface.inside:
			return -1
		return self.rface.n

	def is_locally_delaunay (e):
		return geom.in_circle(
			e.twin.lnext.lnext.org.co,
			e.lnext.org.co,
			e.lnext.lnext.org.co,
			e.org.co
		) < 0

	@property
	def goes_left (self):
		return self.dst <= self.org
	@property
	def goes_right (self):
		return self.org <= self.dst
	@property
	def is_internal (self):
		return (self.rface is not None) and self.rface.inside

	@classmethod
	def _make_twins (cls):
		e = cls()
		t = cls()
		e.twin = t
		t.twin = e

		e._pair_index = t._pair_index = cls._pair_counter
		cls._pair_counter += 1
		e._ptr_index = cls._ptr_counter
		t._ptr_index = cls._ptr_counter + 1
		cls._ptr_counter += 2

		return e, t

	@classmethod
	def Make (cls, enext: Self):
		e, t = cls._make_twins()

		# Make sure eNext points to the first edge of the edge pair
		# First as in, first allocated. I don't get it either. =_="
		if enext.twin._ptr_index < enext._ptr_index:
			enext = enext.twin


		ePrev: Edge = enext.twin.next
		t.next = ePrev
		ePrev.twin.next = e
		e.next = enext
		enext.twin.next = t


		e.onext = e
		e.lnext = t

		t.onext = t
		t.lnext = e

		return e

	def kill (self):
		if self.twin._ptr_index < self._ptr_index:
			self = self.twin

		eNext = self.next
		ePrev = self.twin.next
		eNext.twin.next = ePrev
		ePrev.twin.next = eNext

	def splice_insert (a: Self, b: Self):
		aOnext = a.onext
		bOnext = b.onext

		aOnext.twin.lnext = b
		bOnext.twin.lnext = a
		a.onext = bOnext
		b.onext = aOnext

	def splice (eOrg: Self, eDst: Self):
		if eOrg is eDst:
			return
		
		joining_loops = False
		joining_vertices = False

		if eDst.org is not eOrg.org:
			joining_vertices = True
			Vertex.kill(eDst.org, eOrg.org)

		if eDst.lface is not eOrg.lface:
			joining_loops = True
			Face.kill(eDst.lface, eOrg.lface)

		eDst.splice_insert(eOrg)

		if not joining_vertices:

			Vertex.Make(eDst, eOrg.org)
			eOrg.org.edge = eOrg

		if not joining_loops:
			Face.Make(eDst, eOrg.lface)
			eOrg.lface.edge = eOrg

	def add_edge_vertex (self):
		eNew = Edge.Make(self)
		eNewSym = eNew.twin

		eNew.splice_insert(self.lnext)
		eNew.org = self.dst

		Vertex.Make(eNewSym, eNew.org)
		eNew.lface = eNewSym.lface = self.lface
		return eNew

	def split (self):
		temp_edge = self.add_edge_vertex()

		eNew = temp_edge.twin

		self.twin.splice_insert(self.twin.oprev)
		self.twin.splice_insert(eNew)

		self.dst = eNew.org
		eNew.dst.edge = eNew.twin
		eNew.rface = self.rface
		eNew.winding = self.winding
		eNew.twin.winding = self.twin.winding

		return eNew

	def leq (e1: Self, e2: Self, event: Vertex):
		'''
		Original signature of `Tesselator tess, ActiveRegion reg1, ActiveRegion reg2`
		
		Remap by passing in the arguments `reg1.edge_up, reg2.edge_up, tess.event`
		'''
		if e1.dst is event:
			if e2.dst is event:
				if e1.org <= e2.org:
					return geom.edge_sign(e2.dst.co, e1.org.co, e2.org.co) <= 0
				return geom.edge_sign(e1.dst.co, e2.org.co, e1.org.co) >= 0
			return geom.edge_sign(e2.dst.co, event.co, e2.org.co) <= 0
		if e2.dst is event:
			return geom.edge_sign(e1.dst.co, event.co, e1.org.co) >= 0

		t1 = geom.edge_eval(e1.dst.co, event.co, e1.org.co)
		t2 = geom.edge_eval(e2.dst.co, event.co, e2.org.co)
		return t1 >= t2

	def connect (e_org: Self, e_dst: Self) -> Self:
		en = e_org.Make(e_org)
		tn = en.twin

		joining_loops = False
		if e_dst.lface is not e_org.lface:
			joining_loops = True
			Face.kill(e_dst.lface, e_org.lface)

		en.splice_insert(e_org.lnext)
		tn.splice_insert(e_dst)

		en.org = e_org.dst
		tn.org = e_dst.org
		en.lface = tn.lface = e_org.lface

		e_org.lface.edge = tn

		if not joining_loops:
			Face.Make(en, e_org.lface)
		return en

	def add_winding (dst: Self, src: Self):
		dst.winding += src.winding
		dst.twin.winding += src.twin.winding


class Mesh:
	__slots__ = 'vhead', 'fhead', 'ehead', 'thead'
	def __init__ (self):
		self.vhead = Vertex()
		self.fhead = Face()
		e = self.ehead = Edge()
		t = self.thead = Edge()
		e.twin = t
		t.twin = e

	def check_self (self):
		fhead = self.fhead
		vhead = self.vhead
		ehead = self.ehead

		f = fprev = fhead
		while (f := fprev.next) is not fhead:
			assert f.prev is fprev
			e = f.edge
			while True:
				assert e.twin is not e
				assert e.twin.twin is e
				assert e.lnext.onext.twin is e
				assert e.onext.twin.lnext is e
				assert e.lface is f
				if (e:=e.lnext) is f.edge:
					break
			fprev = f
		assert f.prev is fprev and f.edge is None

		v = vprev = vhead
		while (v := vprev.next) is not vhead:
			assert v.prev is vprev
			e = v.edge
			while True:
				assert e.twin is not e
				assert e.twin.twin is e
				assert e.lnext.onext.twin is e
				assert e.onext.twin.lnext is e
				assert e.org is v
				if (e:=e.onext) is v.edge:
					break
			vprev = v
		assert v.prev is vprev and v.edge is None

		e = eprev = ehead
		while (e:=eprev.next) is not ehead:
			assert e.twin.next is eprev.twin
			assert e.twin is not e
			assert e.twin.twin is e
			assert e.org is not None
			assert e.dst is not None
			assert e.lnext.onext.twin is e
			assert e.onext.twin.lnext is e
			eprev = e
		assert (
			(e.twin.next is eprev.twin) and
			(e.twin is self.thead) and
			(e.twin.twin is e) and
			(e.org is None and e.dst is None) and # linter breaks here, dunno why.
			(e.lface is None and e.rface is None)
		)

	def make_edge (self):
		e = Edge.Make(self.ehead)

		Vertex.Make(e, self.vhead)
		Vertex.Make(e.twin, self.vhead)

		Face.Make(e, self.fhead)
		return e

	def set_winding_number (self, value: int, keep_only_bounds: bool):
		# for e in +self.ehead:
		for (e, _) in self.ehead.walk_head_safe():
			if e.rface.inside != e.lface.inside:
				e.winding = value if e.lface.inside else -value
			else:
				if not keep_only_bounds:
					e.winding = 0
				else:
					del self[e]

	def __delitem__ (self, e_del: Edge):
		joining_loops = False

		if e_del.lface is not e_del.rface:
			joining_loops = True
			Face.kill(e_del.lface, e_del.rface)

		if e_del.onext is e_del:
			Vertex.kill(e_del.org, None)
		else:
			e_del.rface.edge = e_del.oprev
			e_del.org.edge = e_del.onext

			e_del.splice_insert(e_del.oprev)

			if not joining_loops:
				Face.Make(e_del, e_del.lface)
		
		t_del = e_del.twin
		if t_del.onext is t_del:
			Vertex.kill(t_del.org, None)
			Face.kill(t_del.lface, None)
		else:
			e_del.lface.edge = t_del.oprev
			t_del.org.edge = t_del.onext
			t_del.splice_insert(t_del.oprev)

		Edge.kill(e_del)

	def compute_normal (self, norm: Co3 = None):
		vHead = self.vhead

		v = vHead.start
		minVal = Co3()
		maxVal = Co3()
		minVert = [v] * 3
		maxVert = [v] * 3

		for v in +vHead:
			for i, c in enumerate(v.src_co):
				if c < minVal[i]:
					minVal[i] = c
					minVert[i] = v
				if c > maxVal[i]:
					maxVal[i] = c
					maxVert[i] = v

		i = 0
		if (maxVal.y - minVal.y) > (maxVal.x - minVal.x):
			i = 1
		if (maxVal.z - minVal.z) > (maxVal[i] - minVal[i]):
			i = 2

		if norm is None:
			norm = Co3()

		if minVal[i] >= maxVal[i]:
			norm.xyz = (0.0, 0.0, 1.0)
			return

		maxLen2 = 0
		v1 = minVert[i]
		v2 = maxVert[i]
		d1 = v1.src_co - v2.src_co
		
		for v in +vHead:
			d2 = v.src_co - v2.src_co
			tNorm = d1.cross(d2)
			tLen2 = tNorm.length_squared()
			if tLen2 > maxLen2:
				maxLen2 = tLen2
				norm.xyz = tNorm

		if maxLen2 <= 0:
			norm.xyz = Co3.ZERO
			norm[geom.short_axis(d1)[0]] = 1.0
		
		return norm

	def check_orientation (self):
		area = 0
		for f in +self.fhead:
			e = f.edge
			if e.winding <= 0:
				continue
			
			while True:
				area += (e.org.co.x - e.dst.co.x) * (e.org.co.y + e.dst.co.y)
				if (e := e.lnext) is f.edge:
					break

		if area < 0:
			for v in +self.vhead:
				v.co.y = - v.co.y

	def remove_degen_edges (self):
		eHead = self.ehead

		e = eHead.next
		while e is not eHead:
			eNext = e.next
			eLnext = e.lnext

			if (e.org == e.dst) and e.lnext.lnext is not e:
				Edge.splice(eLnext, e)
				del self[e]
				e = eLnext
				eLnext = e.lnext

			if eLnext.lnext is e:
				if eLnext is not e:
					if eLnext is eNext or eLnext is eNext.twin:
						eNext = eNext.next
					del self[eLnext]
				if e is eNext or e is eNext.twin:
					eNext = eNext.next
				del self[e]
			e = eNext

	def remove_degen_faces (self):
		for (f, _) in self.fhead.walk_head_safe():
			e = f.edge
			assert e.lnext is not e

			if e.lnext.lnext is e:
				Edge.add_winding(e.onext, e)
				del self[e]

	def tessellate_interior (self):
		for (f, _) in self.fhead.walk_head_safe():
			if f.inside:
				self.tessellate_mono_region(f)

	def tessellate_mono_region (self, face: Face):
		up = face.edge
		assert (up.lnext is not up) and (up.lnext.lnext is not up)

		while up.dst <= up.org:
			up = up.lprev

		while up.org <= up.dst:
			up = up.lnext

		lo = up.lprev
		while up.lnext is not lo:
			if up.dst <= lo.org:
				while (
					(lo.lnext is not up) and
					(
						lo.lnext.goes_left or
						geom.edge_sign(lo.org.co, lo.dst.co, lo.lnext.dst.co) <= 0
					)
				):
					lo = Edge.connect(lo.lnext, lo).twin
				lo = lo.lprev
			else:
				while (
					(lo.lnext is not up) and
					(
						(up.lprev.goes_right) or
						geom.edge_sign(up.dst.co, up.org.co, up.lprev.org.co) >= 0
					)
				):
					up = Edge.connect(up, up.lprev).twin
				up = up.lnext

		assert lo.lnext is not up
		while lo.lnext.lnext is not up:
			lo = Edge.connect(lo.lnext, lo).twin

	def merge_convex_faces (self, max_verts_per_face: int):
		for e in +self.ehead:
			eNext = e.next
			eSym = e.twin
			if e.twin is None:
				continue
			
			if (e.lface is None) or (not e.lface.inside):
				continue
			if (eSym.lface is None) or (not eSym.lface.inside):
				continue

			leftNv = Face.count_verts(e.lface)
			rightNv = Face.count_verts(eSym.lface)
			if (leftNv + rightNv - 2) > max_verts_per_face:
				continue

			va = e.lprev.org.co
			vb = e.org.co
			vc = e.twin.lnext.dst.co

			vd = e.twin.lprev.org.co
			ve = e.twin.org.co
			vf = e.lnext.dst.co

			if geom.vert_ccw(va, vb, vc) and geom.vert_ccw(vd, ve, vf):
				if (e is eNext) or (e is eNext.twin):
					eNext = eNext.next
				del self[e]

	def flip_edge (mesh, edge: Edge):
		a0 = edge
		a1 = a0.lnext
		a2 = a1.lnext
		b0 = edge.twin
		b1 = b0.lnext
		b2 = b1.lnext

		aOrg = a0.org
		aOpp = a2.org
		bOrg = b0.org
		bOpp = b2.org

		fa = a0.lface
		fb = b0.lface

		assert edge.is_internal
		assert a2.lnext is a0
		assert b2.lnext is b0

		a0.org = bOpp
		a0.onext = b1.twin
		b0.org = aOpp
		b0.onext = a1.twin
		a2.onext = b0
		b2.onext = a0
		b1.onext = a2.twin
		a1.onext = b2.twin

		a0.lnext = a2
		a2.lnext = b1
		b1.lnext = a0

		b0.lnext = b2
		b2.lnext = a1
		a1.lnext = b0

		a1.lface = fb
		b1.lface = fa

		fa.edge = a0
		fb.edge = b0

		if (aOrg.edge is a0):
			aOrg.edge = b1
		if (bOrg.edge is b0):
			bOrg.edge = a1

		for k in (a0, a1, a2, b0, b1, b2):
			assert k.lnext.onext.twin is k
			assert k.onext.twin.lnext is k
			assert k.org.edge.org is k.org

		for k in (aOrg, bOrg):
			assert k.edge.org is k

		assert a0.oprev.onext.org is a0.org

	def refine_delaunay (mesh):
		stack = list[Edge]()
		max_faces = 0
		for f in +mesh.fhead:
			if f.inside:
				e = f.edge
				while True:
					e.mark = e.is_internal
					if e.mark and not e.twin.mark:
						stack.append(e)
					if (e:=e.lnext) is f.edge:
						break
			max_faces += 1

		max_iter = max_faces ** 2
		cur_iter = 0

		while len(stack) > 0 and cur_iter < max_iter:
			e = stack.pop()
			e.mark = e.twin.mark = False
			if not Edge.is_locally_delaunay(e):
				mesh.flip_edge(e)
				edges = (
					e.lnext,
					e.lprev,
					e.twin.lnext,
					e.twin.lprev
				)
				for ek in edges:
					if not ek.mark and ek.is_internal:
						ek.mark = ek.twin.mark = True
						stack.append(ek)
			cur_iter += 1


class ActiveRegion:
	def __init__ (self):
		self.edge_up: Edge = None
		self.node_up: Ordered.Node[ActiveRegion] = None

		self.winding_number: int = 0

		self.is_inside   = False
		self.is_sentinel = False
		self.is_dirty    = False
		self.fix_upper   = False

	@property
	def below (self):
		return self.node_up.prev.item

	@property
	def above (self):
		return self.node_up.next.item

	def finish (self):
		e = self.edge_up
		f = e.lface

		f.inside = self.is_inside
		f.edge = e
		self.delete()

	def delete (self):
		if self.fix_upper:
			# It was created with zero winding number, so it better be
			# deleted with zero winding number (ie. it better not get merged
			# with a real edge).
			assert self.edge_up.winding == 0
		self.edge_up.active_region = None
		self.node_up.delete()


class Tesselator:
	def __init__ (self):
		self.mesh: Mesh = None

		self.normal = Co3()

		self.bounds_min = Co2()
		self.bounds_max = Co2()

		self.process_cdt = False
		self.reverse_contours = False

		self.winding_rule = WindingRule.ODD

		self.dict = Ordered[ActiveRegion]()
		self.pq = Priority[Vertex]()
		self.event: Vertex = None

		self.vertex_index_counter = 0
		# self.vertex_count
		# self.element_count
		self.verticies      = list[float]()
		self.elements       = list[int]()
		self.vertex_indices = list[int]()

		self.S_UNIT = Co2(1.0, 0.0)
		self.TOLERANCE_NONZERO = False

	def auto_set_s_unit (self, mode: str):
		if mode == 'flat':
			self.S_UNIT = Co2(1.0, 0.0)
		elif mode == 'slanted':
			# Pre-normalized
			self.S_UNIT = Co2(0.50941539564955385, 0.86052074622010633)
		elif mode == 'random':
			import random
			self.S_UNIT = Co2(
				random.random() * 2.0 - 1.0,
				random.random() * 2.0 - 1.0
			).normalized()

	def project_polygon (self, true_project = True):
		norm = Co3(self.normal)
		computed_normal = False
		if norm == Co3.ZERO:
			self.mesh.compute_normal(norm)
			computed_normal = True
		
		SUX, SUY = self.S_UNIT
		
		i = geom.long_axis(norm)[0]
		j = (i + 1) % 3
		k = (i + 2) % 3

		sUnit = Co3()
		sUnit[i, j, k] = 0.0, SUX, SUY
		tUnit = Co3()

		if true_project:
			norm.normalize()

			w = sUnit.dot(norm)
			sUnit -= w * norm
			sUnit.normalize()

			tUnit = norm.cross(sUnit).normalized()
		else:
			gt = norm[i] > 0
			tUnit[i, j, k] = (
				0.0,
				-SUY if gt else +SUY,
				+SUX if gt else -SUX
			)

		vhead = self.mesh.vhead
		for v in +vhead:
			v.co.xy = (
				v.src_co.dot(sUnit),
				v.src_co.dot(tUnit)
			)
		
		if computed_normal:
			self.mesh.check_orientation()

		first = True
		for v in +vhead:
			if first:
				self.bounds_min.x = self.bounds_max.x = v.co.x
				self.bounds_min.y = self.bounds_max.y = v.co.y
				first = False
			else:
				x, y = v.co
				self.bounds_min.x = min(self.bounds_min.x, x)
				self.bounds_min.y = min(self.bounds_min.y, y)
				self.bounds_max.x = max(self.bounds_max.x, x)
				self.bounds_max.y = max(self.bounds_max.y, y)

	def add_contour (self, vertices):
		"""
		Libtess2's signature here was:
		`(..., int size, void* verts, int stride, int num_verts)`,
		however here I've opted to omit `size`, `stride`, and `num_verts`.
		
		These are all determined by the first element of `verts`, and `verts` itself.

		`size = len(verts[0])`
		
		`num_verts = len(verts)`
		"""
		
		if self.mesh is None:
			self.mesh = Mesh()

		size = max(2, min(3, len(vertices[0])))

		e = None
		for co in vertices:
			if e is None:
				e = self.mesh.make_edge()
				Edge.splice(e, e.twin)
			else:
				Edge.split(e)
				e = e.lnext
			
			src_co = e.org.src_co
			src_co.xyz = co[:2], 0.0
			if size > 2:
				src_co.z = co[2]

			e.org.index = self.vertex_index_counter
			self.vertex_index_counter += 1

			e.winding      = -1 if self.reverse_contours else +1
			e.twin.winding = -e.winding

	def tesselate (
		self,
		winding_rule: WindingRule,
		element_type: ElementType,
		poly_size,
		vert_size,
		normal: Co3 = None
	):
		self.vertex_index_counter = 0

		if normal is not None:
			self.normal.xyz = normal

		self.winding_rule = winding_rule

		vert_size = max(2, min(3, vert_size))

		if self.mesh is None:
			return

		self.project_polygon()
		self.compute_interior()

		mesh = self.mesh

		rc = 1
		if element_type is ElementType.BOUNDARY_CONTOURS:
			rc = mesh.set_winding_number(1, True)
		else:
			rc = mesh.tessellate_interior()
			if rc != 0 and self.process_cdt:
				mesh.refine_delaunay()

		if rc == 0:
			raise NotImplementedError("I dont actually know what this means =/")

		mesh.check_self()

		if element_type is ElementType.BOUNDARY_CONTOURS:
			self.output_contours(vert_size) #FIXME
		else:
			self.output_polymesh(self.mesh, element_type, poly_size, vert_size)

		# tessMeshDeleteMesh( &tess->alloc, mesh );
		self.mesh = None

	def compute_winding (self, reg: ActiveRegion):
		reg.winding_number = reg.above.winding_number + reg.edge_up.winding
		reg.is_inside = self.winding_rule.is_inside(reg.winding_number)

	def add_sentinel (self, x_min: float, x_max: float, y: float):
		"""
		Add two sentinel edges above and below all other edges,
		to avoid special cases at the top and bottom
		"""
		e = self.mesh.make_edge()

		e.org.co.xy = x_max, y
		e.dst.co.xy = x_min, y
		self.event = e.dst

		reg = ActiveRegion()
		reg.edge_up = e

		reg.is_sentinel = True
		reg.node_up = self.dict.insert(reg)

	def init_pq (self, extra_verts = 8):
		vHead = self.mesh.vhead
		vertexCount = sum(1 for _ in +vHead) + extra_verts
		
		pq = self.pq = Priority[Vertex](
			size=vertexCount,
			cmp=lambda a, b: a <= b
		)

		for v in +vHead:
			v.pq_handle = pq.insert(v)

		pq.init()

	def init_edge_dict (self):
		def leq (a: ActiveRegion, b: ActiveRegion):
			return Edge.leq(a.edge_up, b.edge_up, self.event)

		self.dict = Ordered[ActiveRegion](cmp=leq)
		
		w = (self.bounds_max.x - self.bounds_min.x) + 0.01
		h = (self.bounds_max.y - self.bounds_min.y) + 0.01

		xmin = self.bounds_min.x - w
		ymin = self.bounds_min.y - h
		xmax = self.bounds_max.x + w
		ymax = self.bounds_max.y + h

		self.add_sentinel(xmin, xmax, ymin)
		self.add_sentinel(xmin, xmax, ymax)

	def done_edge_dict (self):
		fixed_edges = 0

		reg: ActiveRegion = None
		while (reg:=self.dict.min().item) is not None:
			if not reg.is_sentinel:
				assert reg.fix_upper
				fixed_edges += 1
				assert fixed_edges == 1
			assert reg.winding_number == 0
			reg.delete()
		self.dict = None

	def done_pq (tess):
		tess.pq = None

	def fix_upper_edge (self, reg: ActiveRegion, new: Edge):
		assert reg.fix_upper
		del self.mesh[reg.edge_up]
		reg.fix_upper = False
		reg.edge_up = new
		new.active_region = reg

	def check_for_right_splice (self, reg_up: ActiveRegion):
		regLo = reg_up.below
		e_up = reg_up.edge_up
		e_lo = regLo.edge_up

		if e_up.org <= e_lo.org:
			if geom.edge_sign(e_lo.dst.co, e_up.org.co, e_lo.org.co) > 0:
				return False

			if e_up.org != e_lo.org:
				Edge.split(e_lo.twin)
				Edge.splice(e_up, e_lo.oprev)
				reg_up.is_dirty = regLo.is_dirty = True

			elif e_up.org is not e_lo.org:
				del self.pq[e_up.org.pq_handle]
				Edge.splice(e_lo.oprev, e_up)
		else:
			if geom.edge_sign(e_up.dst.co, e_lo.org.co, e_up.org.co) <= 0:
				return False

			reg_up.above.is_dirty = reg_up.is_dirty = True
			Edge.split(e_up.twin)
			Edge.splice(e_lo.oprev, e_up)
		return True

	def check_for_left_splice (self, reg_up: ActiveRegion):
		regLo = reg_up.below
		eUp = reg_up.edge_up
		eLo = regLo.edge_up

		assert eUp.dst != eLo.dst

		if eUp.dst <= eLo.dst:
			if geom.edge_sign(eUp.dst.co, eLo.dst.co, eUp.org.co) < 0:
				return False
			reg_up.above.is_dirty = reg_up.is_dirty = True
			e = Edge.split(eUp)
			Edge.splice(eLo.twin, e)
			e.lface.inside = reg_up.is_inside
		else:
			if geom.edge_sign(eLo.dst.co, eUp.dst.co, eLo.org.co) > 0:
				return False
			reg_up.is_dirty = regLo.is_dirty = True
			e = Edge.split(eLo)
			Edge.splice(eUp.lnext, eLo.twin)
			e.rface.inside = reg_up.is_inside
		return True

	def topleft_region (self, reg: ActiveRegion):
		org = reg.edge_up.org
		while (reg:=reg.above).edge_up.org is org:
			pass

		if reg.fix_upper:
			e = Edge.connect(reg.below.edge_up.twin, reg.edge_up.lnext)
			self.fix_upper_edge(reg, e)
			reg = reg.above
		return reg

	def topright_region (self, reg: ActiveRegion):
		dst = reg.edge_up.dst
		while (reg:=reg.above).edge_up.dst is dst:
			pass
		return reg

	def add_region_below (self, above: ActiveRegion, news_up: Edge):
		new = ActiveRegion()
		new.edge_up = news_up
		new.node_up = self.dict.insert_before(above.node_up, new)

		news_up.active_region = new
		return new

	def add_right_edges (self,
	                     reg_up: ActiveRegion,
	                     first: Edge,
	                     last: Edge,
	                     topleft_edge: Edge,
	                     do_cleanup: bool):
		# ActiveRegion *reg, *regPrev;
		# TESShalfEdge *e, *ePrev;
		e = first
		while True:
			assert e.org <= e.dst
			self.add_region_below(reg_up, e.twin)
			if (e:=e.onext) is last:
				break

		if topleft_edge is None:
			topleft_edge = reg_up.below.edge_up.rprev

		regPrev = reg_up
		ePrev = topleft_edge
		firstTime = True
		while True:
			reg = regPrev.below
			e = reg.edge_up.twin
			if e.org is not ePrev.org:
				break

			if e.onext is not ePrev:
				Edge.splice(e.oprev, e)
				Edge.splice(ePrev.oprev, e)

			reg.winding_number = regPrev.winding_number - e.winding
			reg.is_inside = self.winding_rule.is_inside(reg.winding_number)

			regPrev.is_dirty = True
			if not firstTime and self.check_for_right_splice(regPrev):
				Edge.add_winding(e, ePrev)
				regPrev.delete()
				del self.mesh[ePrev]
			firstTime = False
			regPrev = reg
			ePrev = e
		regPrev.is_dirty = True
		assert (regPrev.winding_number - e.winding) == reg.winding_number

		if do_cleanup:
			self.walk_dirty_regions(regPrev)

	def connect_left_degen (self, reg_up: ActiveRegion, event: Vertex):
		e = reg_up.edge_up
		if e.org == event:
			assert self.TOLERANCE_NONZERO
			Edge.splice(e, event.edge)
			return

		if e.dst != event:
			Edge.split(e.twin)
			if reg_up.fix_upper:
				del self.mesh[e.onext]
				reg_up.fix_upper = False
			Edge.splice(event.edge, e)
			self.sweep_event(event)
			return

		assert self.TOLERANCE_NONZERO
		reg_up = self.topright_region(reg_up)
		reg = reg_up.below
		eTopRight = reg.edge_up.twin
		eTopLeft = eLast = eTopRight.onext
		if reg.fix_upper:
			assert eTopLeft is not eTopRight
			reg.delete()
			del self.mesh[eTopRight]
			eTopRight = eTopLeft.oprev
		Edge.splice(event.edge, eTopRight)
		if not eTopLeft.goes_left:
			eTopLeft = None
		self.add_right_edges(reg_up, eTopRight.onext, eLast, eTopLeft, True)

	def connect_left_vertex (self, event: Vertex):
		tmp = ActiveRegion()
		tmp.edge_up = event.edge.twin
		regUp = self.dict.search(tmp).item
		regLo = regUp.below
		if regLo is None:
			return # This may happen if the input polygon is coplanar.

		eUp = regUp.edge_up
		eLo = regLo.edge_up

		if geom.edge_sign(eUp.dst.co, event.co, eUp.org.co) == 0:
			self.connect_left_degen(regUp, event)
			return

		reg = regUp if eLo.dst <= eUp.dst else regLo

		eNew: Edge = None
		if regUp.is_inside or reg.fix_upper:
			if reg is regUp:
				eNew = Edge.connect(event.edge.twin, eUp.lnext)
			else:
				eNew = Edge.connect(eLo.dnext, event.edge).twin
			if reg.fix_upper:
				self.fix_upper_edge(reg, eNew)
			else:
				self.compute_winding(self.add_region_below(regUp, eNew))
			self.sweep_event(event)
		else:
			self.add_right_edges(regUp, event.edge, event.edge, None, True)

	def connect_right_vertex (self, reg_up: ActiveRegion, bottomleft_edge: Edge):
		eTopLeft = bottomleft_edge.onext
		regLo = reg_up.below
		eUp = reg_up.edge_up
		eLo = regLo.edge_up

		if eUp.dst is not eLo.dst:
			self.check_for_intersect(reg_up)

		degenerate = False
		if eUp.org == self.event:
			Edge.splice(eTopLeft.oprev, eUp)
			reg_up = self.topleft_region(reg_up)
			eTopLeft = reg_up.below.edge_up
			self.finish_left_regions(reg_up.below, regLo)
			degenerate = True

		if eLo.org == self.event:
			Edge.splice(bottomleft_edge, eLo.oprev)
			bottomleft_edge = self.finish_left_regions(regLo, None)
			degenerate = True
		if degenerate:
			self.add_right_edges(reg_up, bottomleft_edge.onext, eTopLeft, eTopLeft, True)
			return

		if eLo.org <= eUp.org:
			eNew = eLo.oprev
		else:
			eNew = eUp
		eNew = Edge.connect(bottomleft_edge.lprev, eNew)

		# Prevent cleanup, otherwise eNew might disappear before we've even
		# had a chance to mark it as a temporary edge.
		self.add_right_edges(reg_up, eNew, eNew.onext, eNew.onext, False)
		eNew.twin.active_region.fix_upper = True
		self.walk_dirty_regions(reg_up)

	def finish_left_regions (self, first: ActiveRegion, last: ActiveRegion):
		regPrev = first
		ePrev = first.edge_up
		while regPrev is not last:
			regPrev.fix_upper = False
			reg = regPrev.below
			e = reg.edge_up
			if e.org is not ePrev.org:
				if not reg.fix_upper:
					regPrev.finish()
					break
				e = Edge.connect(ePrev.lprev, e.twin)
				self.fix_upper_edge(reg, e)

			if ePrev.onext is not e:
				Edge.splice(e.oprev, e)
				Edge.splice(ePrev, e)
			regPrev.finish() # may change reg->eUp
			ePrev = reg.edge_up
			regPrev = reg
		return ePrev

	def check_for_intersect (self, reg_up: ActiveRegion):
		regLo = reg_up.below
		eUp = reg_up.edge_up
		eLo = regLo.edge_up
		orgUp = eUp.org
		orgLo = eLo.org
		dstUp = eUp.dst
		dstLo = eLo.dst
		# TESSreal tMinUp, tMaxLo;
		# TESSvertex isect, *orgMin;
		# TESShalfEdge *e;

		assert dstLo != dstUp
		assert geom.edge_sign(dstUp.co, self.event.co, orgUp.co) <= 0
		assert geom.edge_sign(dstLo.co, self.event.co, orgLo.co) >= 0
		assert orgUp is not self.event and orgLo is not self.event
		assert not (reg_up.fix_upper or regLo.fix_upper)

		if orgUp is orgLo:
			return False

		y_min_up = min(orgUp.co.y, dstUp.co.y)
		y_max_lo = max(orgLo.co.y, dstLo.co.y)
		if y_min_up > y_max_lo:
			return False

		if orgUp <= orgLo:
			if geom.edge_sign(dstLo.co, orgUp.co, orgLo.co) > 0:
				return False
		else:
			if geom.edge_sign(dstUp.co, orgLo.co, orgUp.co) < 0:
				return False

		# DebugEvent( tess );
		isect = geom.edge_intersect(dstUp.co, orgUp.co, dstLo.co, orgLo.co)

		assert min(orgUp.co.y, dstUp.co.y) <= isect.y
		assert isect.y <= max(orgLo.co.y, dstLo.co.y)
		assert min(dstLo.co.x, dstUp.co.x) <= isect.x
		assert isect.x <= max(orgLo.co.x, orgUp.co.x)

		if isect <= self.event.co:
			isect.xy = self.event.co
		orgMin = orgUp if orgUp <= orgLo else orgLo
		
		if orgMin.co <= isect:
			isect.xy = orgMin.co

		if isect == orgUp.co or isect == orgLo.co:
			self.check_for_right_splice(reg_up)
			return False

		if ((
				dstUp != self.event and
				geom.edge_sign(dstUp.co, self.event.co, isect) >= 0
			) or (
				dstLo != self.event and
				geom.edge_sign(dstLo.co, self.event.co, isect) <= 0
		)):
			if dstLo is self.event:
				Edge.split(eUp.twin)
				Edge.splice(eLo.twin, eUp)
				reg_up = self.topleft_region(reg_up)
				eUp = reg_up.below.edge_up
				self.finish_left_regions(reg_up.below, regLo)
				self.add_right_edges(reg_up, eUp.oprev, eUp, eUp, True)
				return True

			if dstUp is self.event:
				Edge.split(eLo.twin)
				Edge.splice(eUp.lnext, eLo.oprev)
				regLo = reg_up
				reg_up = self.topright_region(reg_up)
				e = reg_up.below.edge_up.rprev
				regLo.edge_up = eLo.oprev
				eLo = self.finish_left_regions(regLo, None)
				self.add_right_edges(reg_up, eLo.onext, eUp.rprev, e, True)
				return True
			if geom.edge_sign(dstUp.co, self.event.co, isect) >= 0:
				reg_up.above.is_dirty = reg_up.is_dirty = True
				Edge.split(eUp.twin)
				eUp.org.co.xy = self.event.co

			if geom.edge_sign(dstLo.co, self.event.co, isect) <= 0:
				reg_up.is_dirty = regLo.is_dirty = True
				Edge.split(eLo.twin)
				eLo.org.co.xy = self.event.co
			return False

		Edge.split(eUp.twin)
		Edge.split(eLo.twin)
		Edge.splice(eLo.oprev, eUp)
		eUp.org.co.xy = isect
		eUp.org.pq_handle = self.pq.insert(eUp.org)
		
		Vertex.get_intersect_data(eUp.org, orgUp, dstUp, orgLo, dstLo)
		reg_up.above.is_dirty = reg_up.is_dirty = regLo.is_dirty = True
		return False

	def walk_dirty_regions (self, reg_up: ActiveRegion):
		reg_lo = reg_up.below

		while True:
			while reg_lo.is_dirty:
				reg_up = reg_lo
				reg_lo = reg_lo.below

			if not reg_up.is_dirty:
				reg_lo = reg_up
				reg_up = reg_up.above
				if (reg_up is None) or (not reg_up.is_dirty):
					return # walked all dirty regions
			reg_up.is_dirty = False
			e_up = reg_up.edge_up
			e_lo = reg_lo.edge_up

			if e_up.dst is not e_lo.dst:
				if self.check_for_left_splice(reg_up):
					if reg_lo.fix_upper:
						reg_lo.delete()
						del self.mesh[e_lo]
						reg_lo = reg_up.below
						e_lo = reg_lo.edge_up
					elif reg_up.fix_upper:
						reg_up.delete()
						del self.mesh[e_up]
						reg_up = reg_lo.above
						e_up = reg_up.edge_up

			if e_up.org is not e_lo.org:
				if (
					(e_up.dst is not e_lo.dst) and
					not (reg_up.fix_upper or reg_lo.fix_upper) and
					((e_up.dst is self.event) or (e_lo.dst is self.event))
				):
					if self.check_for_intersect(reg_up):
						return # Called recursively; finished now.
				else:
					self.check_for_right_splice(reg_up)

			if (e_up.org is e_lo.org) and (e_up.dst is e_lo.dst):
				Edge.add_winding(e_lo, e_up)
				reg_up.delete()
				del self.mesh[e_up]
				reg_up = reg_lo.above

	def sweep_event (self, event: Vertex):
		self.event = event
		# DebugEvent( tess );

		e = event.edge
		while e.active_region is None:
			e: Edge = e.onext

			if e is event.edge:
				self.connect_left_vertex(event)
				return

		regUp: ActiveRegion = self.topleft_region(e.active_region)
		reg = regUp.below
		eTopLeft = reg.edge_up
		eBottomLeft: Edge = self.finish_left_regions(reg, None)

		if eBottomLeft.onext is eTopLeft:
			self.connect_right_vertex(regUp, eBottomLeft)
		else:
			self.add_right_edges(regUp, eBottomLeft.onext, eTopLeft, eTopLeft, True)

	def compute_interior (self):
		self.mesh.remove_degen_edges()

		self.init_pq()
		self.init_edge_dict()

		pq = self.pq
		while not pq.is_empty:
			v = pq.extract_min()
			while True:
				vNext = pq.min()
				if (vNext is None) or (vNext != v):
					break
				vNext = pq.extract_min()
				Edge.splice(v.edge, vNext.edge)

			self.sweep_event(v)

		self.event = self.dict.min().item.edge_up.org
		# DebugEvent( tess );
		self.done_edge_dict()
		self.done_pq()

		self.mesh.remove_degen_faces()
		self.mesh.check_self()

	def output_polymesh (
		self,
		mesh: Mesh,
		element_type: ElementType,
		poly_size: int,
		vert_size: int
	):
		#TODO: this should be outputting tuples, not a homogenous list?
		# tldr just make it more pythonic i guess =_="
		maxFaceCount = 0
		maxVertexCount = 0
		elements = 0

		if poly_size > 3:
			mesh.merge_convex_faces(poly_size)

		for v in +mesh.vhead:
			v.n = -1

		for f in +mesh.fhead:
			f.n = -1
			if not f.inside:
				continue

			edge = f.edge
			faceVerts = 0
			while True:
				v = edge.org
				if v.n == -1:
					v.n = maxVertexCount
					maxVertexCount += 1
				faceVerts += 1
				if (edge:=edge.lnext) is f.edge:
					break

			assert faceVerts <= poly_size

			f.n = maxFaceCount
			maxFaceCount += 1

		self.element_count = maxFaceCount
		if element_type is ElementType.CONNECTED_POLYGONS:
			maxFaceCount *= 2
		
		elements = self.elements = [0] * (maxFaceCount * poly_size)
		
		self.vertex_count = maxVertexCount
		verts = self.verticies = [0.0] * (maxVertexCount * vert_size)
		
		self.vertex_indices = [0] * (maxVertexCount)

		for v in +mesh.vhead:
			if v.n != -1:
				vert = v.n * vert_size
				verts[vert] = v.src_co.x
				verts[vert+1] = v.src_co.y
				if vert_size > 2:
					vert[vert+2] = v.src_co.z
				self.vertex_indices[v.n] = v.index

		ec = 0
		for f in +mesh.fhead:
			if not f.inside:
				continue

			edge = f.edge
			faceVerts = 0
			while True:
				v = edge.org
				elements[ec] = v.n
				ec += 1
				faceVerts += 1
				if (edge:=edge.lnext) is f.edge:
					break
			
			for _ in range(faceVerts, poly_size):
				elements[ec] = -1
				ec += 1

			if element_type == 'connected_polygons':
				edge = f.edge
				while True:
					elements[ec] = edge.neighbour_face
					ec += 1
					if (edge:=edge.lnext) is f.edge:
						break
				for _ in range(faceVerts, poly_size):
					elements[ec] = -1
					ec += 1


