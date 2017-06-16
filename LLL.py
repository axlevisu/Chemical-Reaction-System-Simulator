import numpy as np

def gram_schmidt(B):
	oB = [B[0]]
	for b in B[1:]:
		x = 0
		for ob in oB:
			x += (b.dot(ob)/(1.0*ob.dot(ob)))*ob
		oB.append(b - x)
	return np.array(oB)

def basis_reduce(C,delta = 3./4):
	B = C
	oB = gram_schmidt(B)
	mu = B.dot(oB.T)/(np.sum(oB*oB,axis=-1))
	n = B.shape[0]
	k=1
	while k<n:
		j = k-1
		while j>=0:
			if np.abs(mu[k][j]) > 0.5:
				B[k] = B[k] - round(mu[k][j])*B[j]
				oB = gram_schmidt(B)
				mu = B.dot(oB.T)/(np.sum(oB*oB,axis=-1))
			j-=1
		if oB[k].dot(oB[k]) >= (delta - mu[k][k-1])*oB[k].dot(oB[k]):
			k+=1
		else:
			B[k] = B[k-1] + B[k]
			B[k-1] = B[k] - B[k-1]
			B[k] = B[k] - B[k-1]
			oB = gram_schmidt(B)
			mu = B.dot(oB.T)/(np.sum(oB*oB,axis=-1))
			k = max(k-1,1)
	return B

def main():
	A = np.array([[1,4,2,5],[2,2,6,3]])*1.
	print A
	print LLL(A)

if __name__ == '__main__':
	main()