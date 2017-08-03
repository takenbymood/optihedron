#Credit: Rosetta code https://rosettacode.org/wiki/Gray_code#Python

def genBin(n):
	if n:
		bits = []
		while n:
			n,remainder = divmod(n, 2)
			bits.insert(0, remainder)
		return bits
	else: return [0]

def readBin(bits):
	i = 0
	for bit in bits:
		i = i * 2 + bit
	return i

def binToGray(bits):
	return bits[:1] + [i ^ ishift for i, ishift in zip(bits[:-1], bits[1:])]

def grayToBin(bits):
	b = [bits[0]]
	for nextb in bits[1:]: b.append(b[-1] ^ nextb)
	return b

def genCode(n):
	return binToGray(genBin(n))

def readCode(bits):
	return readBin(grayToBin(bits))