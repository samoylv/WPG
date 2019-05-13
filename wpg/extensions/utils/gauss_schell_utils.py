def modes2D(m,N):
    """
    return N pairs of Transverse Gauss-Hermite Mode Order pairs with order up to m. 

    """
    A=list()
    for i in range(1,m+1):
       A.append ( list(itertools.product([0,i],repeat=2)))
       
    # combine   
    A=list(itertools.chain.from_iterable(A))
    
    # remove duplicates
    temp = []
    for a,b in A:
        if (a,b) not in temp: #to check for the duplicate tuples
            temp.append((a,b))
    
    return temp[0:N]# -*- coding: utf-8 -*-

