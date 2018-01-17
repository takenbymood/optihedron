import string
import random

def randomStr(N): 
	return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(N))