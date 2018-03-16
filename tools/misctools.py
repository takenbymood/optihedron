import string
import random
import os

def randomStr(N): 
	return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(N))


def removeByExtension(dirName,extension):
	test = os.listdir(dirName)
	for item in test:
		if item.endswith(extension):
			os.remove(os.path.join(dirName, item))
			print("removed {}".format(item))


def removeByPattern(dirName,pattern):
	test = os.listdir(dirName)
	for item in test:
		if pattern in item:
			os.remove(os.path.join(dirName, item))
			print("removed {}".format(item))

def toInt(s):
	try:
		return int(s)
	except ValueError:
		print('error converting '+str(s)+' to integer value')
		return s