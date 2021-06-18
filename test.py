from qr_factorisation import HouseHolder_decomp,givens_rotation,Gram_Schmidt_process
import numpy as np

if __name__ == "__main__":
	A=np.random.rand(5,3)
	val = input("Enter your value: 1-HouseHolder Process 2-Givens Process 3-Gram Schmidt Process	")
	val=int(val)
	if val==1:
		Q,R=HouseHolder_decomp(A)
		print('HouseHolder')
		print('Q=')
		print(Q)
		print('R=')
		print(R)
	else:
		if val==2:
			Q,R=givens_rotation(A)
			print('Givens Rotation')
			print('Q=')
			print(Q)
			print('R=')
			print(R)
		else:
			if val==3:
				Q,R=Gram_Schmidt_process(A)
				print('Gram Schmidt')
				print('Q=')
				print(Q)
				print('R=')
				print(R)
			else:
				val=input("Enter your value: 1-HouseHolder Process 2-Givens Process 3-Gram Schmidt Process")
