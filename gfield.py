# GF arithemtics class
# GF(2^m) operations class in Python 
# Basic Information :
# Author : Colin
# Start Date : Nov 15 2013
# Function: a GF class defining finite field arithemtics and functions for general polynomial operations 
# Update note: Nov 18: added substitute method , forshadowing the future syndrome calculation
# TODO: -- some methods are better turned to standalone functions;
		# -- better documentation and commenting



from itertools import izip 
import copy

#primitive polynomial table , represented in integer form
primpoly = {
	1:0b11,
	2:0b111,
	3:11,
	4:19,
	5:37,
	6:67,
	7:137,
	8:285,
	9:529,
	10:1033
	}
# Genator of logarithm table of GF(2^m) with m as an argument 
def logtable(m=4):
	"create a discrete logarithm table as per m in GF(2^m)"
	order = 2**m 
	gf_exp = [0]*(2*order)
	gf_log = [0]*order 
	
	gf_exp[0] = 1 
	
	x = 1 
	for i in xrange(1,order-1):
		x <<= 1
		if x & order:
			x ^= primpoly[m]
		gf_exp[i] = x 
		gf_log[x] = i 
	for i in xrange(order-1,2*order):
		gf_exp[i] = gf_exp[i-order+1] 
	table = {'exp':gf_exp, 'log':gf_log}

	return table 
			
							
# ---------------------------- GENERAL POLYNOMIAL OPERATION FUNCTIONS added Nov 22 ALL INPUTS ARE INTEGERS-----------------------

def vec2num(v):
	"input: v as a vector and convert to a number output (binary)"
	assert v 
	s = ''.join( map(str,v) ) 
	return int(s,2)   

def poly_deg(v):
	"find the degree of input polynomial represented as a integer or binary"
	if v:
		r = -1
		while v:
			v = v >> 1
			r += 1
		return r 
	else:
		return 0

def poly_mul(x,y,vector_out=False): 
		"general polynomial multiplication without reducing to primitive polynomial, return a integer form of poly" 
		z = 0
		i = 0
		while (y>>i) > 0:
			if y & (1<<i):
				z ^= x<<i 
			i += 1
		if vector_out:
			return map(int,bin(z)[2:]) 
		else:
			return z 

def poly_add(x,y,vector_out=False):
	"polynomial addition. input are as numbers "

	z = x ^ y 

	if vector_out:
		assert z 
		return map(int,bin(z)[2:] )
	else:
		return z 


def poly_div(N, D, vector_out=False):
	"General polynomial division.  input N,D as numbers ; output: r=N mod D "
	if D==0:
			raise ZeroDivisionError() 
		#nested helper function to see the degree of a ascending poly vector
	def degree(ascending_poly):
		while ascending_poly and ascending_poly[-1]==0:
			ascending_poly.pop()
		return len(ascending_poly)-1 
	# ---------------------------------------------------------
	dD = poly_deg(D)
	dN = poly_deg(N)  # degree of D and N 
	Dvec = map(int, bin(D)[2:] )[::-1]
	Nvec = map(int, bin(N)[2:] )[::-1] # ascending repr of poly
	
 	q = [0]*dN 
 	while dN>=dD:
 		d = [0]*(dN - dD) + Dvec
 		q[dN-dD] = 1 
 		
 		Nvec = [(coeffN ^ coeffd) for coeffN,coeffd in izip(Nvec,d)] 
 		dN = degree(Nvec)
 	r = Nvec[::-1]
 	if r == []:   # D|N  
 		return 0  
 	q = q[::-1]  # switch back to descending poly vector
	if not vector_out: # if vector out is off. turn out integer form 
		r = int( ''.join(map(str, r)), 2) # convert to numerical output 
	return r






# --------------------------------------------------END OF POLYNOMIAL OPERATION FUNCTIONS ,follows the class FiniteField ------
class FiniteField(object):
	"""The finite field class does the follows:
	first, use order of the generator polynomial as a argument m, and construct a GF(2^m) using logarithm table 
	as a internal data attributes. 
	multiplication and division as well as xor operations are defined. It suits various value of m using a dictionary
	to swtich to the right order
	the methods are :
	----------------------------------------------------------
	add   i.e. xor 
	subtract i.e. xor equivalent to add in GF 
	mul :  multiplication modulo generator polynomial
	div : division 
	showpoly  : show the polynomial representation of input number
	showvector: show code vector representation of input number

	--------some private member functions for internal use and debugging only ------
	_poly_add: general polynomial addition without reduction 
	_poly_mul:  x,y as integer , output a z integer according to general poly mult 
	_poly_div : N,D as integer , output the N mod D in integer form 
	_reduce: take a integer/binary number x and gives x mod primitive_polynomial 
	_vec2num: take a descending order poly in list form L, and output converted binary/integer 
	    	 example: _vec2num([1,1,0,0]) returns 12
	--------------------------------------------------------------
	"""
	def __init__(self, m):
		super(FiniteField, self).__init__()
		self.m = m 
		self.primpoly = primpoly[m] #primitive poly in number repr
		self.order=2**m
		self.gf_exp = logtable(m)['exp'] 
		self.gf_log = logtable(m)['log']


	def _poly_mul(self,x,y): 
		"general polynomial multiplication without reducing to primitive polynomial, return a integer form of poly" 
		z = 0
		i = 0
		while (y>>i) > 0:
			if y & (1<<i):
				z ^= x<<i 
			i += 1
		return z 

	def _poly_div(self,N,D):  # remains to test , Nov 15 2013. use some test cases. general poly division over GF(2) 
		"""long division: input descending poly (number) and output q,r in a descending poly vector form,
		 r is again converted back to number 
		 ----------------------------------------------------------------
		 input : N(int), D(int) 
		 return r= N mod D ; int form still 
		 ----------------------------------------------------------------
		 """
		if D==0:
			raise ZeroDivisionError() 
		#nested helper function to see the degree of a ascending poly vector
		def degree(ascending_poly):
			while ascending_poly and ascending_poly[-1]==0:
				ascending_poly.pop()
			return len(ascending_poly)-1 
		# ---------------------------------------------------------
		dD = self.finddegree(D)
		dN = self.finddegree(N)
		Dvec = self.showvector(D)[::-1]
  		Nvec = self.showvector(N)[::-1] # ascending repr of poly
		
	 	q = [0]*dN 
	 	while dN>=dD:
	 		d = [0]*(dN - dD) + Dvec
	 		q[dN-dD] = 1 
	 		
	 		Nvec = [(coeffN ^ coeffd) for coeffN,coeffd in izip(Nvec,d)] 
	 		dN = degree(Nvec)
	 	r = Nvec[::-1]
	 	if r == []:   
	 		return 0  
	 	q = q[::-1]  # switch back to descending poly vector
		r = self._vec2num(r) # convert to numerical output 
		return r 

	def _vec2num(self,x):
		"convert a list repr of poly to a binary number, added on Nov 15"

		# convert the x , e.g. [0,0,1,0,1] to str  ing '00101' and then strip leading 0s 
		# get '101' and convert it to a binary integer int(x,2) base 2 
		assert x != []  # x is non-empty	
		string_rep = ''.join(map(str, x))
		# string_rep = string_rep.lstrip('0') 
		out_num = int(string_rep, 2) 
		return out_num 

	
					
	def _reduce(self,x):
		"take a input x(number) and modulo primitive poly, output a number repr of poly"
		remainder = self._poly_div(x,self.primpoly) 
		return remainder



	def _poly_add(self,x,y):
		"general polynomial addition over GF(2), without modulo primpoly"
		summ = x^y
		return summ 
		

		

# ------------------------all the operations in the GF(2^m) are modulo primpoly -------------------------

	def add(self,x,y):
		assert x>=0 and y>=0 
		summ = x ^ y
		reduced_sum = self._reduce(summ)
		return reduced_sum 

	def subtract(self,x,y):
		return self.add(x, y)

	def mul(self,x,y):
		assert y>=0 and x >= 0 
		x = self._reduce(x) 
		y = self._reduce(y) #reduce inputs to modulo primpoly
		if x==0 or y==0:
			return 0
		return self.gf_exp[self.gf_log[x] + self.gf_log[y]]

	def div(self,x,y):
		assert x>=0 and y>=0 
		x = self._reduce(x) 
		y = self._reduce(y) 
		if y==0:
			raise ZeroDivisionError()
		if x==0:
			return 0
		return self.gf_exp[self.gf_log[x] + self.order-1 - self.gf_log[y]]

 	def finddegree(self,v):
 		"find the degree of input polynomial represented as a integer or binary"
 		if v:
 			r = -1
 			while v:
 				v = v >> 1
 				r += 1
 			return r 
 		else:
 			return 0

 	def showpoly(self,f):
 		"string representation of a input polynomial"
 		deg = self.finddegree(f)
 		s = ''

 		if f==0:
 			return '0'
 		for i in range(deg,0,-1):
 			if (1<<i) & f:
 				s = s + (' a^'+ repr(i))
 		if 1 & f:
 			s = s + ' 1'
 		return s.strip().replace(' ','+') 

 	def showvector(self,f):
 		"f is a integer number repr of polynomial, output a vector "
 		# assert f>=0 
 		bstring = bin(f) 
 		bstring = bstring[2:] #bstring[0] = 1 
 		vec = map(int,bstring)
 		return vec 


 	def substitute(self,fin,powr,exp_out=False):
 		"""given a received vector descending polynomial, substitute alpha^i to X ,used for syndrome calculation later
 		------------------------------------
 		input : powr is the power of alpha i in alpha^pow fin(alpha^power) 
 		input: fin is a descending polynomial vector or a integer number 
 		output: a integer or a exponential of alpha if exp_out is True 
 		CAUTION: COPY the fin first. what if fin is 0 ?  
 		-------------------------------------"""
 		# assert fin    # fin cannot be zero or [] 
 		if fin==0 and not exp_out:  # fin(X) = 0 
 			return 0 

 		fin2=copy.copy(fin)

 		if not isinstance(fin2,list):  # if fin2 is not a vector, convert it to a vector first
 			fin2 = self.showvector(fin2)
 			N = len(fin2)
 		else:
 			N = len(fin2)
 		result = 0
 		for i in xrange(0,N):
 			if fin2[i] :
 				result = self.add(self.gf_exp[(N-1-i)*powr], result) #modulo addition

 		if exp_out:
 			assert result != 0 # log(0) is invalid 
 			return self.gf_log[result] # return the logarithm of the result polynomial 
 		return result  







if __name__ == '__main__':
	a = FiniteField(4)
	a.showpoly(15)
	a.showvector(10)
	print bin(a.div(0b1001,0b1110))

	print 'numerize a vector ', a._vec2num([0,0,1,0,1])






