
# Author:Colin 
# Date Nov 18 2013 
# Testing simple Hamming code for C(39,32),C(72,64),C(137,128),C(266,256) 







import numpy
import random
import sys

 
def createMessage(length):
    msg = []
    for i in range(length):
        letter = random.choice([0,1])
        msg.append(letter)
    return numpy.array(msg, dtype = int)
 
def encode(m, g):
    "encoder by Generator matrix g.  m * G = c "
    enc = numpy.dot(m, g)%2
    return enc
 
def decode(received, h):
    "syndrome calculation syn = r x H' or H x r' = syn' "
    syndrome = numpy.dot(h, received)%2
    return syndrome
 
def noise(m, error, bit):
    "error is error rate, bit is total #corrupted bits,m is the input vector; insert noise in the message"
    noisy_msg = numpy.copy(m)
    cont = 0
    for i in range(len(noisy_msg)):
        e = random.random() # return a random number in [0,1)
        if e <= error:
            if noisy_msg[i] == 0:
                noisy_msg[i] = 1
                cont +=1
            else:
                noisy_msg[i] = 0
                cont +=1
            if cont == bit:  # bit is  #corrupted bits
                break
    return noisy_msg #return noise mask 
 
def findError(synd, H_matrix):
    # need to altered by Colin , return the location index of error 
    "error locator function. return the position of error"
    if all(synd==0):
        return -1 # no error 
    block_length = int(H_matrix.shape[1]) # H is a rxn matrix, so n is the code length 
    e = numpy.zeros(block_length,int) 
    for i in xrange(0,block_length):
        e1 = numpy.copy(e)  #copy original zeros vector
        e1[i] = 1
        compare = numpy.dot(e1,H_matrix.transpose()) 
        if numpy.array_equal(compare,synd):
            return i # locator index 
        else:
            continue 
    else:
        return 1000 # beyond correction capability        
 
def correct(noisy, n):
    "flip bits at n-th position"
    if noisy[n] == 0:
        noisy[n] = 1
    else:# noisy[n] == 1:
        noisy[n] = 0
    #print noisy
    return noisy
 
def hamming(length, n, error, bit):
    g =  numpy.array([[1, 0, 0, 0, 0, 1, 1],[0, 1, 0, 0, 1, 0, 1],[0, 0, 1, 0, 1, 1, 0],[0, 0, 0, 1, 1, 1, 1]])
    h = numpy.array([[0, 0, 0, 1, 1, 1, 1],[0, 1, 1, 0, 0, 1, 1],[1, 0, 1, 0, 1, 0, 1]]) # colin. remove a extra , comma
    # corrected = 0
    # uncorrected = 0
    for i in range(n):
        msg = createMessage(length)
        enc = encode(msg, g)
        print "original message: ", msg
        print "encoded msg: ", enc
        noisy = noise(enc, error, bit)
        print 'noisy/received mesage: ', noisy
        dec = decode(noisy, h)
        print 'syndrome vector: ', dec
        error_index = findError(dec, h)
        if error_index >= 0 and error_index < 1000:
            corrected_vector = correct(noisy,error_index)
            print 'corrected vector is ',corrected_vector 
        else:
            print 'no error detected!' 
    #         else:
    #             uncorrected+=1
    # print 'corrected: ', corrected
    # print 'uncorrected: ', uncorrected
        
        
if __name__ == '__main__':
    
    length = int(raw_input('bits length: '))
    n = int(raw_input('Repetitions: '))
    error = float(raw_input('Error percentage: '))
    bit = int(raw_input('Bits with error: '))
    hamming(length,n,error, bit)