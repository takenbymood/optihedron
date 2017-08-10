def subdivide(l,n):
	return [l[x:x+n] for x in xrange(0, len(l), n)]