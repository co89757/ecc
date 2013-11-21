
# Author:Colin 
# Date Nov 18 2013 
# Testing simple Hamming code for C(39,32),C(72,64),C(137,128),C(266,256) 
# UPDATE: Nov 20 2013: Include DED functionality by overall parity bit ;
# use error_status_code to switch to error scenarios detected 
#TODO: larger G/H for wider information word length 







import numpy
import random
# import sys

 
def createMessage(length):
    msg = []
    for i in range(length):
        letter = random.choice([0,1])
        msg.append(letter)
    return numpy.array(msg, dtype = int)
 
def encode(m, g):
    "encoder by Generator matrix g.  m * G = c , works on SEC part ; Return a numpy array "
    enc = numpy.dot(m, g)%2
    return enc
 
def syndrome(received, h):
    "syndrome calculation syn = r x H' or H x r' = syn' , return a numpy array"
    synd = numpy.dot(h, received)%2
    return synd
 
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
    # need to altered by Colin , return the location index of error, works on the SEC part  
    "error locator function. return the position of error"
    if all(synd==0):
        return -1 # no error , clean syndrome returns -1 
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
    for i in range(n): # n incoming received message repetitions 
        msg = createMessage(length) # generate a random information vector u(x)
        enc = encode(msg, g) 

        # ------- check parity right after encoder : parity before corruption ----------
        op_before =  list(enc).count(1) % 2 
        print "original message: ", msg, " parity: ",op_before 
        print "encoded msg: ", enc

        # ------introduce corruption --------------
        noisy = noise(enc, error, bit)

        #------check parity after corruption --------
        op_after = list(noisy).count(1) % 2 
        # --------DED_flag denotes the status of DED part, 1 means detected parity mismatch -- --
        DED_flag = op_before ^ op_after 
        print 'noisy/received mesage: ', noisy," parity(after pollution): ", op_after 
        # -----compute syndrome , set the SEC flag ------------------------------
        syndromes = syndrome(noisy, h)
        print 'syndrome vector: ', syndromes
         # -------- error pattern / error location returned if only 1 error bit, if clean, error_index = -1 ----
        error_index = findError(syndromes, h) 

        SEC_flag = int(error_index >= 0)  # SEC_flag 1 denotes error alert in SEC part 

        # create a error status code , a tuple of SEC and DED flags 
        ERROR_STATUS_CODE = (SEC_flag,DED_flag)  

        #---------------CORRECTION ACTION AS PER THE STATUS OF (SEC_FLAG,DED_FLAG) ----------
        # ------------------------------------------------------------------------------------
        if ERROR_STATUS_CODE == (0,0):
            print 'Clean! No errors !'
        elif ERROR_STATUS_CODE == (1,1):
            corrected_vector = correct(noisy, error_index)
            print 'corrected vector : ', corrected_vector
        elif ERROR_STATUS_CODE == (1,0):
            print "2 errors detected! Unable to correct! " 
        else:
            print "The impossible occurs, something wrong! "


        


        # if error_index >= 0 and error_index < 1000:
        #     corrected_vector = correct(noisy,error_index)
        #     print 'corrected vector is ',corrected_vector 
        # else: # if synd = 0 
        #     print 'no error detected!' 





        
        
if __name__ == '__main__':
    
    length = int(raw_input('bits length: '))
    n = int(raw_input('Repetitions: '))
    error = float(raw_input('Error percentage: '))
    bit = int(raw_input('Bits with error: '))
    hamming(length,n,error, bit)