def create_field(nm,q,x,y):
    import pykat as pk
    TEM = pk.optics.gaussian_beams.HG_mode(q,n=nm[0],m=nm[1])
    field = TEM.Unm(y,x)
    return field


def overlap(Field_in,Field_nm,x,y):
    import numpy as np
    #numerical overlap integral.
    norm = (x[len(x)-1]-x[0])*(y[len(y)-1]-y[0])/(len(x)*len(y)) #normilization factor
    ccField= np.conj(Field_nm)
    C = np.sum(Field_in*ccField) * norm
    return C


def mode_indices(max_order):
    mode_list = []
    for i in range(max_order+1):
        for j in range(i+1):
            mode_list.append((j,i-j))
    return mode_list


def HG_mode_content(Field_in,q,x,y,max_order=None):
    if max_order==None:
        max_order=10
        
    Var_nm = mode_indices(max_order)
    number_of_coeffs=int((max_order+1)*(max_order+2)/2)
    mode_content=[]
    
    for i in range(number_of_coeffs):
        Field_nm=create_field(Var_nm[i],q,x,y) #create HGnm mode field
        Cnm = overlap(Field_in,Field_nm,x,y) #calculate HGnm content #mode indicies and complex weight
        mode_content.append((Var_nm[i],Cnm)) #mode indicies and complex weight 
        
    return mode_content


def recompose(Cnm,mode,q,x,y): #how to structure variables
    # Cnm complex list of weights
    # mode list of indixes (n,m)
    # q beam parameter 
    # x,y position vectors defining grid size 
    
    import numpy as np
    
    recomposed_field, trash = np.meshgrid(x,y)
    recomposed_field.astype(complex)
    
    for i in range(len(mode)):
        Unm = create_field(mode[i],q,x,y)
        recomposed_field = recomposed_field + Unm * Cnm[i]
        
    return recomposed_field