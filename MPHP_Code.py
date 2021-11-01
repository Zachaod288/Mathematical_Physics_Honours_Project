import numpy as np
import scipy
from scipy.special import gamma
import time
start = time.process_time()

def As(A, B_i, i):
    Ak = []

    for l in range(len(A)):
        for m in range(len(B_i[i-1])):
            if l == 0:
                Ak.append(A[l]*B_i[i-1-l][-(m+1)])
            if l> 0:
                if 2*l+m>=len(Ak):
                    break
                temp = Ak[2*l+m]
                Ak[2*l+m] += A[l]*B_i[i-1-l][-(m+1)]


   
    return Ak 

def approx(n): #Definition returning the approximation given
    return (-1)**(n+1) * (6/(scipy.pi**3))**(0.5) * gamma(n+0.5)*(3**n)


def main():
    num = int(input("How many A values to calculate?: "))
    start = time.process_time()
    #Will do A1 and A2 outside of the loop as they don't depend on previous A's
    B_i = [1]    # Polynomial coefficients will go here
    A = []
    b_i = np.zeros(2)
    b_div = np.zeros(2)
    b_2div = np.zeros(2)

    for j in range(1, len(b_i)+1):

        b_i[j-1] = (-1)*(1/(2**j)) 
        b_div[j-1] = b_i[j-1]*(2*j) 
        b_2div[j-1] = b_div[j-1]*(2*j-1)


    b_temp = []
    b_vals = (1/b_div[0])*b_i[1]
    b_temp.append(b_vals)
    b_vals = (1/b_div[1])*(b_2div[-1]*b_vals)
    b_temp.append(b_vals)


    A.append(b_vals)
    B_i.append(np.multiply(b_temp[::-1], b_i))


    #Now for i>2

    for i in range(2,num+1): # upper number gives A_(n+1) values
        A_vals = As(A, B_i, i)

        b_i = np.zeros(2*i)
        b_div = np.zeros(2*i) # First derivative
        b_2div = np.zeros(2*i) # Second derivative

        for j in range(1, len(b_i)+1):
            b_i[j-1] = (-1)**i*(1/(2**j)) 
            b_div[j-1] = b_i[j-1]*(2*j)    # Generally coefficients of the B_i,j term you're solving for
            b_2div[j-1] = b_div[j-1]*(2*j-1) # Generally coefficient of the B_i,(j+1) term you've just found
 

        b_vals = (1/b_div[-1])*(1/4)*(-B_i[i-1][-1]) # Initialsie the B_i_n value
        b_temp = []   # Will contain the B_i_j values
        b_temp.append(b_vals)

        b_vals = (1/b_div[-2])*(b_2div[-1]*b_vals - (1/4)*B_i[i-1][-2])
        b_temp.append(b_vals)

        
        for k in range(3, 2*i+1):   # leaves the last two values. All these terms contain A's
            if 2*i-k>1:
                #K is the kth highest ordered term.   e.g if k=1 at x^16, then k=3 at x^12

                b_vals = (1/b_div[-k])*(b_2div[-(k-1)]*b_vals - (1/4)*B_i[i-1][-(k)] + A_vals[k-3])
                #        from xB'_n(x)    -B''_n(x)              (1/4)x^4*B_(n-1)
                b_temp.append(b_vals)
            elif 2*i-k<=1:
                b_vals = (1/b_div[-k])*(b_2div[-(k-1)]*b_vals + A_vals[k-3])
                b_temp.append(b_vals)


        if i%2== 0:
            A.append(-b_vals)
        elif i%2 == 1:
            A.append(b_vals)
        B_i.append(np.multiply(b_temp[::-1], b_i))  # reverse the order of the list and append it to B_i

   
    #print(Approx)
    print("A values are: ")
    print("n, A_n")
    for i in range(len(A)):
        print(i+1,  A[i],)
    print("\n")
    Approx = []
    print("Approximations are: ")
    for i in range(len(A)):
        Approx.append(approx(i))
        print(i+1, Approx[i])
        
    print("Time taken :", time.process_time() - start)


    return
main()
