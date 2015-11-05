var DefaultBezierContext = require('./default-context').DefaultBezierContext;

var N = 4;
var MAX_DEPTH = 5;

function arraycopy(src, j, dst, jd, n) {
	for (var k = 0; k < n; k++) {
		dst[jd + k] = src[j + k]
	}
}

function SpiroSeg() {
	this.x = 0;
	this.y = 0;
	this.type = 'corner';
	this.bend_th = 0;
	this.ks = [0, 0, 0, 0]
	this.seg_ch = 0;
	this.seg_th = 0;
	this.l = 0;
}

function BandMat() {
	this.a = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
	this.al = [0, 0, 0, 0, 0]
}
BandMat.prototype.copyfrom = function(from) {
	this.a = from.a.slice(0);
	this.al = from.al.slice(0);
}
BandMat.arraycopy = function(from, fromi, to, toi, nelem) {
	for (var i = 0; i < nelem; ++i)
		to[i + toi].copyfrom(from[i + fromi]);
}

function integrate_spiro(ks, xy) {
	var th1 = ks[0];
	var th2 = .5 * ks[1];
	var th3 = (1. / 6) * ks[2];
	var th4 = (1. / 24) * ks[3];
	var x, y;
	var ds = 1. / N;
	var ds2 = ds * ds;
	var ds3 = ds2 * ds;
	var k0 = ks[0] * ds;
	var k1 = ks[1] * ds;
	var k2 = ks[2] * ds;
	var k3 = ks[3] * ds;
	var i;
	var s = .5 * ds - .5;

	x = 0;
	y = 0;

	for (i = 0; i < N; i++) {
		var u, v;
		var km0, km1, km2, km3;

		km0 = (((1. / 6) * k3 * s + .5 * k2) * s + k1) * s + k0;
		km1 = ((.5 * k3 * s + k2) * s + k1) * ds;
		km2 = (k3 * s + k2) * ds2;
		km3 = k3 * ds3; {
			var t1_1 = km0;
			var t1_2 = .5 * km1;
			var t1_3 = (1. / 6) * km2;
			var t1_4 = (1. / 24) * km3;
			var t2_2 = t1_1 * t1_1;
			var t2_3 = 2 * (t1_1 * t1_2);
			var t2_4 = 2 * (t1_1 * t1_3) + t1_2 * t1_2;
			var t2_5 = 2 * (t1_1 * t1_4 + t1_2 * t1_3);
			var t2_6 = 2 * (t1_2 * t1_4) + t1_3 * t1_3;
			var t2_7 = 2 * (t1_3 * t1_4);
			var t2_8 = t1_4 * t1_4;
			var t3_4 = t2_2 * t1_2 + t2_3 * t1_1;
			var t3_6 = t2_2 * t1_4 + t2_3 * t1_3 + t2_4 * t1_2 + t2_5 * t1_1;
			var t3_8 = t2_4 * t1_4 + t2_5 * t1_3 + t2_6 * t1_2 + t2_7 * t1_1;
			var t3_10 = t2_6 * t1_4 + t2_7 * t1_3 + t2_8 * t1_2;
			var t4_4 = t2_2 * t2_2;
			var t4_5 = 2 * (t2_2 * t2_3);
			var t4_6 = 2 * (t2_2 * t2_4) + t2_3 * t2_3;
			var t4_7 = 2 * (t2_2 * t2_5 + t2_3 * t2_4);
			var t4_8 = 2 * (t2_2 * t2_6 + t2_3 * t2_5) + t2_4 * t2_4;
			var t4_9 = 2 * (t2_2 * t2_7 + t2_3 * t2_6 + t2_4 * t2_5);
			var t4_10 = 2 * (t2_2 * t2_8 + t2_3 * t2_7 + t2_4 * t2_6) + t2_5 * t2_5;
			var t5_6 = t4_4 * t1_2 + t4_5 * t1_1;
			var t5_8 = t4_4 * t1_4 + t4_5 * t1_3 + t4_6 * t1_2 + t4_7 * t1_1;
			var t5_10 = t4_6 * t1_4 + t4_7 * t1_3 + t4_8 * t1_2 + t4_9 * t1_1;
			var t6_6 = t4_4 * t2_2;
			var t6_7 = t4_4 * t2_3 + t4_5 * t2_2;
			var t6_8 = t4_4 * t2_4 + t4_5 * t2_3 + t4_6 * t2_2;
			var t6_9 = t4_4 * t2_5 + t4_5 * t2_4 + t4_6 * t2_3 + t4_7 * t2_2;
			var t6_10 = t4_4 * t2_6 + t4_5 * t2_5 + t4_6 * t2_4 + t4_7 * t2_3 + t4_8 * t2_2;
			var t7_8 = t6_6 * t1_2 + t6_7 * t1_1;
			var t7_10 = t6_6 * t1_4 + t6_7 * t1_3 + t6_8 * t1_2 + t6_9 * t1_1;
			var t8_8 = t6_6 * t2_2;
			var t8_9 = t6_6 * t2_3 + t6_7 * t2_2;
			var t8_10 = t6_6 * t2_4 + t6_7 * t2_3 + t6_8 * t2_2;
			var t9_10 = t8_8 * t1_2 + t8_9 * t1_1;
			var t10_10 = t8_8 * t2_2;
			u = 1;
			v = 0;
			v += (1. / 12) * t1_2 + (1. / 80) * t1_4;
			u -= (1. / 24) * t2_2 + (1. / 160) * t2_4 + (1. / 896) * t2_6 + (1. / 4608) * t2_8;
			v -= (1. / 480) * t3_4 + (1. / 2688) * t3_6 + (1. / 13824) * t3_8 + (1. / 67584) * t3_10;
			u += (1. / 1920) * t4_4 + (1. / 10752) * t4_6 + (1. / 55296) * t4_8 + (1. / 270336) * t4_10;
			v += (1. / 53760) * t5_6 + (1. / 276480) * t5_8 + (1. / 1.35168e+06) * t5_10;
			u -= (1. / 322560) * t6_6 + (1. / 1.65888e+06) * t6_8 + (1. / 8.11008e+06) * t6_10;
			v -= (1. / 1.16122e+07) * t7_8 + (1. / 5.67706e+07) * t7_10;
			u += (1. / 9.28973e+07) * t8_8 + (1. / 4.54164e+08) * t8_10;
			v += (1. / 4.08748e+09) * t9_10;
			u -= (1. / 4.08748e+10) * t10_10;
		}

		{
			var th = (((th4 * s + th3) * s + th2) * s + th1) * s;
			var cth = Math.cos(th);
			var sth = Math.sin(th);

			x += cth * u - sth * v;
			y += cth * v + sth * u;
			s += ds;
		}
	}

	xy[0] = x * ds;
	xy[1] = y * ds;
}

function compute_ends(ks, ends, seg_ch) {
	var xy = [0, 0];
	var ch, th;
	var l, l2, l3;
	var th_even, th_odd;
	var k0_even, k0_odd;
	var k1_even, k1_odd;
	var k2_even, k2_odd;

	integrate_spiro(ks, xy);
	ch = Math.hypot(xy[0], xy[1]);
	th = Math.atan2(xy[1], xy[0]);
	l = ch / seg_ch;

	th_even = .5 * ks[0] + (1. / 48) * ks[2];
	th_odd = .125 * ks[1] + (1. / 384) * ks[3] - th;
	ends[0][0] = th_even - th_odd;
	ends[1][0] = th_even + th_odd;
	k0_even = l * (ks[0] + .125 * ks[2]);
	k0_odd = l * (.5 * ks[1] + (1. / 48) * ks[3]);
	ends[0][1] = k0_even - k0_odd;
	ends[1][1] = k0_even + k0_odd;
	l2 = l * l;
	k1_even = l2 * (ks[1] + .125 * ks[3]);
	k1_odd = l2 * .5 * ks[2];
	ends[0][2] = k1_even - k1_odd;
	ends[1][2] = k1_even + k1_odd;
	l3 = l2 * l;
	k2_even = l3 * ks[2];
	k2_odd = l3 * .5 * ks[3];
	ends[0][3] = k2_even - k2_odd;
	ends[1][3] = k2_even + k2_odd;

	return l;
}

function compute_pderivs(s, ends, derivs, jinc) {
	var recip_d = 2e6;
	var delta = 1. / recip_d;
	var try_ks = [0, 0, 0, 0];
	var try_ends = [
		[0, 0, 0, 0],
		[0, 0, 0, 0]
	];
	var i, j, k;

	compute_ends(s.ks, ends, s.seg_ch);
	for (i = 0; i < jinc; i++) {
		for (j = 0; j < 4; j++)
			try_ks[j] = s.ks[j];
		try_ks[i] += delta;
		compute_ends(try_ks, try_ends, s.seg_ch);
		for (k = 0; k < 2; k++)
			for (j = 0; j < 4; j++)
				derivs[j][k][i] = recip_d * (try_ends[k][j] - ends[k][j]);
	}
}

function mod_2pi(th) {
	var u = th / (2 * Math.PI);
	return 2 * Math.PI * (u - Math.floor(u + 0.5));
}

function setup_path(src, n) {
	var n_seg = src[0].type === 'open' ? n - 1 : n;
	var r = [];
	for (var j = 0; j < n_seg + 1; j++) {
		r[j] = new SpiroSeg
	}
	var i;
	var ilast;

	for (i = 0; i < n_seg; i++) {
		r[i] = new SpiroSeg;
		r[i].x = src[i].x;
		r[i].y = src[i].y;
		r[i].type = src[i].type;
		r[i].af = src[i].af;
		r[i].ks = [0, 0, 0, 0]
	}
	r[n_seg] = new SpiroSeg;
	r[n_seg].x = src[n_seg % n].x;
	r[n_seg].y = src[n_seg % n].y;
	r[n_seg].type = src[n_seg % n].type;
	r[n_seg].af = src[n_seg % n].af;

	for (i = 0; i < n_seg; i++) {
		var dx = r[i + 1].x - r[i].x;
		var dy = r[i + 1].y - r[i].y;
		r[i].seg_ch = Math.hypot(dx, dy);
		r[i].seg_th = Math.atan2(dy, dx);
	}

	ilast = n_seg - 1;
	for (i = 0; i < n_seg; i++) {
		if (r[i].type === 'open' || r[i].type === 'open_end' || r[i].type === 'corner')
			r[i].bend_th = 0.;
		else
			r[i].bend_th = mod_2pi(r[i].seg_th - r[ilast].seg_th);
		ilast = i;
	}
	return r;
}

function bandec11(m, perm, n) {
	var i, j, k;
	var l;

	/* pack top triangle to the LEFT. */
	for (i = 0; i < 5; i++) {
		for (j = 0; j < i + 6; j++)
			m[i].a[j] = m[i].a[j + 5 - i];
		for (; j < 11; j++) m[i].a[j] = 0.;
	}
	l = 5;
	for (k = 0; k < n; k++) {
		var pivot = k;
		var pivot_val = m[k].a[0];
		var pivot_scale;

		l = l < n ? l + 1 : n;

		for (j = k + 1; j < l; j++)
			if (Math.abs(m[j].a[0]) > Math.abs(pivot_val)) {
				pivot_val = m[j].a[0];
				pivot = j;
			}

		perm[k] = pivot;
		if (pivot !== k) {
			for (j = 0; j < 11; j++) {
				var tmp = m[k].a[j];
				m[k].a[j] = m[pivot].a[j];
				m[pivot].a[j] = tmp;
			}
		}

		if (Math.abs(pivot_val) < 1e-12) pivot_val = 1e-12;
		pivot_scale = 1. / pivot_val;
		for (i = k + 1; i < l; i++) {
			var x = m[i].a[0] * pivot_scale;
			m[k].al[i - k - 1] = x;
			for (j = 1; j < 11; j++)
				m[i].a[j - 1] = m[i].a[j] - x * m[k].a[j];
			m[i].a[10] = 0.;
		}
	}
}

function banbks11(m, perm, v, n) {
	var i, k, l;

	/* forward substitution */
	l = 5;
	for (k = 0; k < n; k++) {
		i = perm[k];
		if (i !== k) {
			var tmp = v[k];
			v[k] = v[i];
			v[i] = tmp;
		}
		if (l < n) l++;
		for (i = k + 1; i < l; i++)
			v[i] -= m[k].al[i - k - 1] * v[k];
	}

	/* back substitution */
	l = 1;
	for (i = n - 1; i >= 0; i--) {
		var x = v[i];
		for (k = 1; k < l; k++)
			x -= m[i].a[k] * v[k + i];
		v[i] = x / m[i].a[0];
		if (l < 11) l++;
	}
}

function compute_jinc(ty0, ty1) {
	if (ty0 === 'g4' || ty1 === 'g4' || ty0 === 'right' || ty1 === 'left')
		return 4;
	else if (ty0 === 'g2' && ty1 === 'g2')
		return 2;
	else if (((ty0 === 'open' || ty0 === 'corner' || ty0 === 'left') && ty1 === 'g2') ||
		(ty0 === 'g2' && (ty1 === 'open_end' || ty1 === 'corner' || ty1 === 'right')))
		return 1;
	else
		return 0;
}

function count_vec(s, nseg) {
	var n = 0;
	for (var i = 0; i < nseg; i++)
		n += compute_jinc(s[i].type, s[i + 1].type);
	return n;
}

function add_mat_line(m, v, derivs, x, y, j, jj, jinc, nmat) {
	if (jj >= 0) {
		var joff = (j + 5 - jj + nmat) % nmat;
		if (nmat < 6) {
			joff = j + 5 - jj;
		} else if (nmat === 6) {
			joff = 2 + (j + 3 - jj + nmat) % nmat;
		}
		v[jj] += x;
		for (var k = 0; k < jinc; k++)
			m[jj].a[joff + k] += y * derivs[k];
	}
}

function spiro_iter(s, m, perm, v, n) {
	var cyclic = s[0].type !== 'open' && s[0].type !== 'corner';
	var i, j, jj;
	var nmat = count_vec(s, n);
	var norm;
	var n_invert;

	for (i = 0; i < nmat; i++) {
		v[i] = 0.;
		m[i] = new BandMat();
		m[i].a = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
		m[i].al = [0, 0, 0, 0, 0];
	}

	j = 0;
	if (s[0].type === 'g4')
		jj = nmat - 2;
	else if (s[0].type === 'g2')
		jj = nmat - 1;
	else
		jj = 0;
	for (i = 0; i < n; i++) {
		var ty0 = s[i].type;
		var ty1 = s[i + 1].type;
		var jinc = compute_jinc(ty0, ty1);
		var th = s[i].bend_th;
		var ends = [
			[0, 0, 0, 0],
			[0, 0, 0, 0]
		];
		var derivs = [
			[
				[0, 0, 0, 0],
				[0, 0, 0, 0]
			],
			[
				[0, 0, 0, 0],
				[0, 0, 0, 0]
			],
			[
				[0, 0, 0, 0],
				[0, 0, 0, 0]
			],
			[
				[0, 0, 0, 0],
				[0, 0, 0, 0]
			]
		];
		var jthl = -1,
			jk0l = -1,
			jk1l = -1,
			jk2l = -1;
		var jthr = -1,
			jk0r = -1,
			jk1r = -1,
			jk2r = -1;

		compute_pderivs(s[i], ends, derivs, jinc);

		/* constraints crossing LEFT */
		if (ty0 === 'g4' || ty0 === 'g2' || ty0 === 'left' || ty0 === 'right') {
			jthl = jj++;
			jj %= nmat;
			jk0l = jj++;
		}
		if (ty0 === 'g4') {
			jj %= nmat;
			jk1l = jj++;
			jk2l = jj++;
		}

		/* constraints on LEFT */
		if ((ty0 === 'left' || ty0 === 'corner' || ty0 === 'open' || ty0 === 'g2') &&
			jinc === 4) {
			if (ty0 !== 'g2') jk1l = jj++;
			jk2l = jj++;
		}

		/* constraints on RIGHT */
		if ((ty1 === 'right' || ty1 === 'corner' || ty1 === 'open_end' || ty1 === 'g2') &&
			jinc === 4) {
			if (ty1 !== 'g2') jk1r = jj++;
			jk2r = jj++;
		}

		/* constraints crossing RIGHT */
		if (ty1 === 'g4' || ty1 === 'g2' || ty1 === 'left' || ty1 === 'right') {
			jthr = jj;
			jk0r = (jj + 1) % nmat;
		}
		if (ty1 === 'g4') {
			jk1r = (jj + 2) % nmat;
			jk2r = (jj + 3) % nmat;
		}

		add_mat_line(m, v, derivs[0][0], th - ends[0][0], 1, j, jthl, jinc, nmat);
		add_mat_line(m, v, derivs[1][0], ends[0][1], -1, j, jk0l, jinc, nmat);
		add_mat_line(m, v, derivs[2][0], ends[0][2], -1, j, jk1l, jinc, nmat);
		add_mat_line(m, v, derivs[3][0], ends[0][3], -1, j, jk2l, jinc, nmat);
		add_mat_line(m, v, derivs[0][1], -ends[1][0], 1, j, jthr, jinc, nmat);
		add_mat_line(m, v, derivs[1][1], -ends[1][1], 1, j, jk0r, jinc, nmat);
		add_mat_line(m, v, derivs[2][1], -ends[1][2], 1, j, jk1r, jinc, nmat);
		add_mat_line(m, v, derivs[3][1], -ends[1][3], 1, j, jk2r, jinc, nmat);
		j += jinc;
	}
	if (cyclic) {
		BandMat.arraycopy(m, 0, m, nmat, nmat);
		BandMat.arraycopy(m, 0, m, 2 * nmat, nmat);
		arraycopy(v, 0, v, nmat, nmat);
		arraycopy(v, 0, v, 2 * nmat, nmat);
		n_invert = 3 * nmat;
		j = nmat;
	} else {
		n_invert = nmat;
		j = 0;
	}
	bandec11(m, perm, n_invert);
	banbks11(m, perm, v, n_invert);
	norm = 0.;
	for (i = 0; i < n; i++) {
		var ty0 = s[i].type;
		var ty1 = s[i + 1].type;
		var jinc = compute_jinc(ty0, ty1);
		for (var k = 0; k < jinc; k++) {
			var dk = v[j++];
			s[i].ks[k] += dk;
			norm += dk * dk;
		}
	}
	return norm;
}

function solve_spiro(s, nseg) {
	var m = [];
	var v = [];
	var perm = [];
	var nmat = count_vec(s, nseg);
	var n_alloc = nmat;
	var norm;

	if (nmat === 0) return 0;
	if (s[0].type !== 'open' && s[0].type !== 'corner') n_alloc *= 3;
	if (n_alloc < 5) n_alloc = 5;
	for (var j = 0; j < n_alloc; j++) {
		m[j] = new BandMat;
		v[j] = 0;
		perm[j] = 0;
	}
	for (var i = 0; i < 10; i++) {
		norm = spiro_iter(s, m, perm, v, nseg);
		if (norm < 1e-12) break;
	}
	return 0;
}
function findIntersection(p1, c1, c2, p2){
	var d1 = {x : c1.x - p1.x, y: c1.y - p1.y}
	var d2 = {x : c2.x - p2.x, y: c2.y - p2.y}
	
	var det = d2.x * d1.y - d2.y * d1.x;
	if(Math.abs(det) < 1e-6) return null;
	var u = ((p2.y - p1.y) * d2.x - (p2.x - p1.x) * d2.y) / det
	var v = ((p2.y - p1.y) * d1.x - (p2.x - p1.x) * d1.y) / det
	if(u <= 0 || v <= 0) return null;
	return {
		x: p1.x + d1.x * u,
		y: p1.y + d1.y * u
	}
}

function findquad(p1, c1, c2, p2) {
	var pt = findIntersection(p1, c1, c2, p2);
	if(!pt) pt = {
		x : (p1.x + p2.x) / 2,
		y : (p1.y + p2.y) / 2
	}
	return pt
}
function spiro_seg_to_bpath(ks, x0, y0, x1, y1, bc, depth, subdivided, isquad, af) {
	var bend = Math.abs(ks[0]) + Math.abs(.5 * ks[1]) + Math.abs(.125 * ks[2]) + Math.abs((1. / 48) * ks[3]);

	if (bend <= 1e-8) {
		bc.lineTo(x1, y1, subdivided);
	} else {
		var seg_ch = Math.hypot(x1 - x0, y1 - y0);
		var seg_th = Math.atan2(y1 - y0, x1 - x0);
		var xy = [];
		var ch, th;
		var scale, rot;
		var th_even, th_odd;
		var ul, vl;
		var ur, vr;

		integrate_spiro(ks, xy);
		ch = Math.hypot(xy[0], xy[1]);
		th = Math.atan2(xy[1], xy[0]);
		scale = seg_ch / ch;
		rot = seg_th - th;
		if (depth > MAX_DEPTH || bend < (isquad ? 0.75 : 1.0)) {
			th_even = (1. / 384) * ks[3] + (1. / 8) * ks[1] + rot;
			th_odd = (1. / 48) * ks[2] + .5 * ks[0];
			ul = (scale * (1. / 3)) * Math.cos(th_even - th_odd);
			vl = (scale * (1. / 3)) * Math.sin(th_even - th_odd);
			ur = (scale * (1. / 3)) * Math.cos(th_even + th_odd);
			vr = (scale * (1. / 3)) * Math.sin(th_even + th_odd);
			if(isquad){
				var pt = findquad({x : x0, y : y0}, {x : x0 + ul, y : y0 + vl}, {x : x1 - ur, y : y1 - vr}, {x: x1, y: y1});
				bc.curveTo(pt.x, pt.y, x1, y1, subdivided);
			} else {
				bc.cubicTo(x0 + ul, y0 + vl, x1 - ur, y1 - vr, x1, y1, subdivided);
			}
		} else {
			/* subdivide */
			var ksub = [];
			var thsub;
			var xysub = [];
			var xmid, ymid;
			var cth, sth;

			ksub[0] = .5 * ks[0] - .125 * ks[1] + (1. / 64) * ks[2] - (1. / 768) * ks[3];
			ksub[1] = .25 * ks[1] - (1. / 16) * ks[2] + (1. / 128) * ks[3];
			ksub[2] = .125 * ks[2] - (1. / 32) * ks[3];
			ksub[3] = (1. / 16) * ks[3];
			thsub = rot - .25 * ks[0] + (1. / 32) * ks[1] - (1. / 384) * ks[2] + (1. / 6144) * ks[3];
			cth = .5 * scale * Math.cos(thsub);
			sth = .5 * scale * Math.sin(thsub);
			integrate_spiro(ksub, xysub);
			xmid = x0 + cth * xysub[0] - sth * xysub[1];
			ymid = y0 + cth * xysub[1] + sth * xysub[0];
			spiro_seg_to_bpath(ksub, x0, y0, xmid, ymid, bc, depth + 1, true, isquad);
			ksub[0] += .25 * ks[1] + (1. / 384) * ks[3];
			ksub[1] += .125 * ks[2];
			ksub[2] += (1. / 16) * ks[3];
			spiro_seg_to_bpath(ksub, xmid, ymid, x1, y1, bc, depth + 1, subdivided, isquad);
		}
	};
	if (af) {
		af.call(bc, x0, y0, x1, y1)
	}
}

function run_spiro(src, n) {
	var nseg = src[0].type === 'open' ? n - 1 : n;
	var s = setup_path(src, n);
	if (nseg > 1) solve_spiro(s, nseg);
	return s;
}

function spiro_to_bpath(s, n, bc, isquad) {
	var nsegs = s[n - 1].type === 'open_end' ? n - 1 : n;
	for (var i = 0; i < nsegs; i++) {
		var x0 = s[i].x;
		var y0 = s[i].y;
		var x1 = s[i + 1].x;
		var y1 = s[i + 1].y;
		if (i === 0) {
			bc.moveTo(x0, y0);
			if(s[0].af) s[0].af.call(bc, x0, y0);
		}
		spiro_seg_to_bpath(s[i].ks, x0, y0, x1, y1, bc, 0, false, isquad, s[i + 1].af);
	}
}

function spiroToBezierOnContext(spiros, isClosed, bc, isquad) {
	var s;
	var n = spiros.length;
	if (n < 1) return;
	if (!isClosed) {
		var oldty_start = spiros[0].type;
		var oldty_end = spiros[n - 1].type;
		spiros[0].type = 'open';
		spiros[n - 1].type = 'open_end';
		s = run_spiro(spiros, n);
		spiros[n - 1].type = oldty_end;
		spiros[0].type = oldty_start;
	} else {
		s = run_spiro(spiros, n);
	}
	spiro_to_bpath(s, n, bc, isquad);
}
function spiroToBezier(spiros, isClosed){
	var c = new DefaultBezierContext();
	spiroToBezierOnContext(spiros, isClosed, c);
	return c.strands;
}

exports.spiroToBezierOnContext = spiroToBezierOnContext;
exports.spiroToBezier = spiroToBezier;