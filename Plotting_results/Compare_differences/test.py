ximport numpy as np

def test_fun(x):
    #x=(x-np.min(x))/(np.max(x)-np.min(x))
    x=x+1
    

def test_fun2(x):
    #x[:]=(x-np.min(x))/(np.max(x)-np.min(x))
    x[:]=x+1
    
a=np.array([[1.0,2.0,3.0],[4.0,5.0,6.0]])
print("id(a)=",id(a))
print("a=",a)

test_fun(a)
print("id(a)=",id(a))
print("a=",a)

test_fun2(a)
print("id(a)=",id(a))
print(a)
