__all__ = ['assemble_stiffnessmatrix1D',
           'assemble_massmatrix1D',
           'assemble_matrix_ex01',
           'assemble_matrix_ex02',
           'assemble_matrix_ex03',
           'assemble_vector_ex01',
           'assemble_norm_ex01'
]

from pyccel.decorators import types

# assembles stiffness matrix 1D
#==============================================================================
@types('int', 'int', 'int[:]', 'double[:,:,:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]')
def assemble_stiffnessmatrix1D(ne, degree, spans, basis, weights, points,  matrix):

    # ... sizes
    k1 = weights.shape[1]
    # ... build matrices
    for ie1 in range(0, ne):
            i_span_1 = spans[ie1]        
            # evaluation dependant uniquement de l'element

            for il_1 in range(0, degree+1):
                i1 = i_span_1 - degree + il_1
                for il_2 in range(0, degree+1):
                            i2 = i_span_1 - degree + il_2
                            v = 0.0
                            for g1 in range(0, k1):
                                
                                    bi_x = basis[ie1, il_1, 1, g1]
                                    bj_x = basis[ie1, il_2, 1, g1]
                                    
                                    wvol = weights[ie1, g1]
                                    
                                    v += bi_x * bj_x * wvol

                            matrix[ degree+ i1, degree+ i2-i1]  += v
    # ...

# assembles mass matrix 1D
#==============================================================================
@types('int', 'int', 'int[:]', 'double[:,:,:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]')
def assemble_massmatrix1D(ne, degree, spans, basis, weights, points, matrix):

    # ... sizes
    k1 = weights.shape[1]
    # ... build matrices
    for ie1 in range(0, ne):
            i_span_1 = spans[ie1]        
            # evaluation dependant uniquement de l'element

            for il_1 in range(0, degree+1):
                i1 = i_span_1 - degree + il_1
                for il_2 in range(0, degree+1):
                            i2 = i_span_1 - degree + il_2
                            v = 0.0
                            for g1 in range(0, k1):
                                
                                    bi_0 = basis[ie1, il_1, 0, g1]
                                    bj_0 = basis[ie1, il_2, 0, g1]
                                    
                                    wvol = weights[ie1, g1]
                                    
                                    v += bi_0 * bj_0 * wvol

                            matrix[degree+i1, degree+ i2-i1]  += v
    # ...

#==============================================================================
#---1 : Assembles stiffness matrix 3D
@types('int', 'int', 'int', 'int', 'int', 'int', 'int[:]', 'int[:]', 'int[:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'real', 'double[:,:,:,:,:,:]')
def assemble_matrix_ex01(ne1, ne2, ne3 ,
                        p1, p2, p3,
                        spans_1, spans_2, spans_3, 
                        basis_1, basis_2, basis_3,
                        weights_1, weights_2, weights_3,
                        points_1, points_2, points_3,
                        Re_Pe, matrix):

    # ... sizes
    k1 = weights_1.shape[1]
    k2 = weights_2.shape[1]
    k3 = weights_3.shape[1]
    # ...
    v_b1 = 1.
    v_b2 = 1.
    v_b3 = 1.

    # ... build matrices
    for ie1 in range(0, ne1):
      i_span_1 = spans_1[ie1]
      for ie2 in range(0, ne2):
        i_span_2 = spans_2[ie2]
        for ie3 in range(0, ne3):
           i_span_3 = spans_3[ie3]

           for il_1 in range(0, p1+1):
             for il_2 in range(0, p2+1):
                for il_3 in range(0, p3+1):

                   for jl_1 in range(0, p1+1):  
                      for jl_2 in range(0, p2+1):
                         for jl_3 in range(0, p3+1):
                            i1 = i_span_1 - p1 + il_1
                            j1 = i_span_1 - p1 + jl_1
        
                            i2 = i_span_2 - p2 + il_2
                            j2 = i_span_2 - p2 + jl_2

                            i3 = i_span_3 - p3 + il_3
                            j3 = i_span_3 - p3 + jl_3

                            v = 0.0
                            for g1 in range(0, k1):
                               for g2 in range(0, k2):
                                  for g3 in range(0, k3):

                                     bi_0 = basis_1[ie1, il_1, 0, g1] * basis_2[ie2, il_2, 0, g2] * basis_3[ie3, il_3, 0, g3]
                                     bi_x = basis_1[ie1, il_1, 1, g1] * basis_2[ie2, il_2, 0, g2] * basis_3[ie3, il_3, 0, g3]
                                     bi_y = basis_1[ie1, il_1, 0, g1] * basis_2[ie2, il_2, 1, g2] * basis_3[ie3, il_3, 0, g3]
                                     bi_z = basis_1[ie1, il_1, 0, g1] * basis_2[ie2, il_2, 0, g2] * basis_3[ie3, il_3, 1, g3]

                                     bj_0 = basis_1[ie1, jl_1, 0, g1] * basis_2[ie2, jl_2, 0, g2] * basis_3[ie3, jl_3, 0, g3]
                                     bj_x = basis_1[ie1, jl_1, 1, g1] * basis_2[ie2, jl_2, 0, g2] * basis_3[ie3, jl_3, 0, g3]
                                     bj_y = basis_1[ie1, jl_1, 0, g1] * basis_2[ie2, jl_2, 1, g2] * basis_3[ie3, jl_3, 0, g3]
                                     bj_z = basis_1[ie1, jl_1, 0, g1] * basis_2[ie2, jl_2, 0, g2] * basis_3[ie3, jl_3, 1, g3]
                                     
                                     # ....
                                     wvol = weights_1[ie1, g1] * weights_2[ie2, g2] * weights_3[ie3, g3]
                                     
                                     # ... Assembles diffusion term
                                     v +=  (bi_x * bj_x + bi_y * bj_y + bi_z * bj_z) * wvol /Re_Pe
                                     # ... Assembles convection term
                                     v +=  ( v_b1 * bj_x + v_b2 * bj_y + v_b3 * bj_z) * bi_0 * wvol
                                     
                            matrix[p1+i1, p2+i2, p3+i3, p1+j1-i1, p2+j2-i2, p3+j3-i3]  += v
    # ...
    

#==============================================================================
#---1 : Assembles mass matrix 3D
@types('int', 'int', 'int', 'int', 'int', 'int', 'int[:]', 'int[:]', 'int[:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:,:,:,:,:]')
def assemble_matrix_ex02(ne1, ne2, ne3 ,
                        p1, p2, p3,
                        spans_1, spans_2, spans_3, 
                        basis_1, basis_2, basis_3,
                        weights_1, weights_2, weights_3,
                        points_1, points_2, points_3,
                        matrix):

    # ... sizes
    k1 = weights_1.shape[1]
    k2 = weights_2.shape[1]
    k3 = weights_3.shape[1]
    # ... build matrices
    for ie1 in range(0, ne1):
      i_span_1 = spans_1[ie1]
      for ie2 in range(0, ne2):
        i_span_2 = spans_2[ie2]
        for ie3 in range(0, ne3):
           i_span_3 = spans_3[ie3]

           for il_1 in range(0, p1+1):
             for il_2 in range(0, p2+1):
                for il_3 in range(0, p3+1):

                   for jl_1 in range(0, p1+1):  
                      for jl_2 in range(0, p2+1):
                         for jl_3 in range(0, p3+1):
                            i1 = i_span_1 - p1 + il_1
                            j1 = i_span_1 - p1 + jl_1
        
                            i2 = i_span_2 - p2 + il_2
                            j2 = i_span_2 - p2 + jl_2

                            i3 = i_span_3 - p3 + il_3
                            j3 = i_span_3 - p3 + jl_3

                            v = 0.0
                            for g1 in range(0, k1):
                               for g2 in range(0, k2):
                                  for g3 in range(0, k3):

                                     bi_0 = basis_1[ie1, il_1, 0, g1] * basis_2[ie2, il_2, 0, g2] * basis_3[ie3, il_3, 0, g3]
                                     # ...
                                     bj_0 = basis_1[ie1, jl_1, 0, g1] * basis_2[ie2, jl_2, 0, g2] * basis_3[ie3, jl_3, 0, g3]

                                     wvol = weights_1[ie1, g1] * weights_2[ie2, g2] * weights_3[ie3, g3]

                                     v += bi_0 * bj_0 * wvol

                            matrix[p1+i1, p2+i2, p3+i3, p1+j1-i1, p2+j2-i2, p3+j3-i3]  += v
    # ...

#==============================================================================
#---1 : Assembles Advection matrix 3D
@types('int', 'int', 'int', 'int', 'int', 'int', 'int[:]', 'int[:]', 'int[:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:,:,:,:,:]')
def assemble_matrix_ex03(ne1, ne2, ne3 ,
                        p1, p2, p3,
                        spans_1, spans_2, spans_3, 
                        basis_1, basis_2, basis_3,
                        weights_1, weights_2, weights_3,
                        points_1, points_2, points_3,
                        matrix):

    # ... sizes
    k1 = weights_1.shape[1]
    k2 = weights_2.shape[1]
    k3 = weights_3.shape[1]
    # ...
    v_b1 = 1.
    v_b2 = 1.
    v_b3 = 1.
    # ... build matrices
    for ie1 in range(0, ne1):
      i_span_1 = spans_1[ie1]
      for ie2 in range(0, ne2):
        i_span_2 = spans_2[ie2]
        for ie3 in range(0, ne3):
           i_span_3 = spans_3[ie3]

           for il_1 in range(0, p1+1):
             for il_2 in range(0, p2+1):
                for il_3 in range(0, p3+1):

                   for jl_1 in range(0, p1+1):  
                      for jl_2 in range(0, p2+1):
                         for jl_3 in range(0, p3+1):
                            i1 = i_span_1 - p1 + il_1
                            j1 = i_span_1 - p1 + jl_1
        
                            i2 = i_span_2 - p2 + il_2
                            j2 = i_span_2 - p2 + jl_2

                            i3 = i_span_3 - p3 + il_3
                            j3 = i_span_3 - p3 + jl_3

                            v = 0.0
                            for g1 in range(0, k1):
                               for g2 in range(0, k2):
                                  for g3 in range(0, k3):

                                     bi_0 = basis_1[ie1, il_1, 0, g1] * basis_2[ie2, il_2, 0, g2] * basis_3[ie3, il_3, 0, g3]
                                     bi_x = basis_1[ie1, il_1, 1, g1] * basis_2[ie2, il_2, 0, g2] * basis_3[ie3, il_3, 0, g3]
                                     bi_y = basis_1[ie1, il_1, 0, g1] * basis_2[ie2, il_2, 1, g2] * basis_3[ie3, il_3, 0, g3]
                                     bi_z = basis_1[ie1, il_1, 0, g1] * basis_2[ie2, il_2, 0, g2] * basis_3[ie3, il_3, 1, g3]

                                     bj_0 = basis_1[ie1, jl_1, 0, g1] * basis_2[ie2, jl_2, 0, g2] * basis_3[ie3, jl_3, 0, g3]
                                     bj_x = basis_1[ie1, jl_1, 1, g1] * basis_2[ie2, jl_2, 0, g2] * basis_3[ie3, jl_3, 0, g3]
                                     bj_y = basis_1[ie1, jl_1, 0, g1] * basis_2[ie2, jl_2, 1, g2] * basis_3[ie3, jl_3, 0, g3]
                                     bj_z = basis_1[ie1, jl_1, 0, g1] * basis_2[ie2, jl_2, 0, g2] * basis_3[ie3, jl_3, 1, g3]
                                     
                                     # ....
                                     wvol = weights_1[ie1, g1] * weights_2[ie2, g2] * weights_3[ie3, g3]
                                     
                                     # ... Assembles convection term
                                     v +=  ( v_b1 * bj_x + v_b2 * bj_y + v_b3 * bj_z) * bi_0 * wvol

                            matrix[p1+i1, p2+i2, p3+i3, p1+j1-i1, p2+j2-i2, p3+j3-i3]  += v
    # ...
    
#==============================================================================
#---1 : Assemble rhs Poisson
@types('int', 'int', 'int', 'int', 'int', 'int', 'int[:]', 'int[:]', 'int[:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:,:]', 'real', 'real', 'real', 'double[:,:,:]')
def assemble_vector_ex01(ne1, ne2, ne3, p1, p2, p3, spans_1, spans_2, spans_3,  basis_1, basis_2, basis_3,  weights_1, weights_2, weights_3, points_1, points_2, points_3, vector_u0, Re_Pe, time, dt, rhs):

    from numpy import exp
    from numpy import pi
    from numpy import sin
    from numpy import arctan2
    from numpy import cos, cosh
    from numpy import sqrt
    from numpy import zeros

    # ... sizes
    k1        = weights_1.shape[1]
    k2        = weights_2.shape[1]
    k3        = weights_3.shape[1]
    
    # ...
    lcoeffs_u = zeros((p1+1,p2+1,p3+1))
    # ...
    lvalues_u = zeros((k1, k2, k3))
    # ...
    v_b1 = 1.
    v_b2 = 1.
    v_b3 = 1.
    # ... build rhs
    for ie1 in range(0, ne1):
        i_span_1 = spans_1[ie1]
        for ie2 in range(0, ne2):
            i_span_2 = spans_2[ie2]
            for ie3 in range(0, ne3):
                i_span_3 = spans_3[ie3]

                lcoeffs_u[ : , : , : ]   = vector_u0[i_span_1 : i_span_1+p1+1, i_span_2 : i_span_2+p2+1, i_span_3 : i_span_3+p3+1]
                for g1 in range(0, k1):
                   for g2 in range(0, k2):
                      for g3 in range(0, k3):

                            wvol  = weights_1[ie1, g1]*weights_2[ie2, g2]*weights_3[ie3, g3]
                            
                            # ... Assembles the initial solution
                            u_0   = 0.0
                            for il_1 in range(0, p1+1):
                                for il_2 in range(0, p2+1):
                                    for il_3 in range(0, p3+1):

                                        bj_0    = basis_1[ie1,il_1,0,g1] * basis_2[ie2,il_2,0,g2] * basis_3[ie3,il_3,0,g3]
                                        
                                        coeff_u = lcoeffs_u[il_1,il_2,il_3]
                                        
                                        u_0    +=  coeff_u * bj_0
                                        
                            x     =  points_1[ie1, g1]
                            y     =  points_2[ie2, g2]
                            z     =  points_3[ie3, g3]
                            t     =  time

                            #.. Test 1
                            f     = pi*cos(pi*t)*sin(pi*x)*sin(pi*y)*sin(pi*z) + (3.*pi**2*sin(pi*t)*sin(pi*x)*sin(pi*y)*sin(pi*z))/Re_Pe
                            f    += v_b1*pi*sin(pi*t)*cos(pi*x)*sin(pi*y)*sin(pi*z) + v_b2*pi*sin(pi*t)*sin(pi*x)*cos(pi*y)*sin(pi*z) + v_b3*pi*sin(pi*t)*sin(pi*x)*sin(pi*y)*cos(pi*z)

                            lvalues_u[g1,g2,g3]  = u_0 + dt * f

                for il_1 in range(0, p1+1):
                  for il_2 in range(0, p2+1):
                    for il_3 in range(0, p3+1):

                      i1 = i_span_1 - p1 + il_1
                      i2 = i_span_2 - p2 + il_2
                      i3 = i_span_3 - p3 + il_3

                      v = 0.0
                      for g1 in range(0, k1):
                         for g2 in range(0, k2):
                            for g3 in range(0, k3):

                              bi_0 = basis_1[ie1, il_1, 0, g1] * basis_2[ie2, il_2, 0, g2] * basis_3[ie3, il_3, 0, g3]

                              wvol = weights_1[ie1, g1] * weights_2[ie2, g2] * weights_3[ie3, g3]

                              u    = lvalues_u[g1,g2,g3]
                              v   += bi_0 * u * wvol

                      rhs[i1+p1,i2+p2,i3+p3] += v   
    # ...

#==============================================================================
#---1 : Assemble rhs of L2 projection
@types('int', 'int', 'int', 'int', 'int', 'int', 'int[:]', 'int[:]', 'int[:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'real', 'double[:,:,:]')
def assemble_vector_ex02(ne1, ne2, ne3, p1, p2, p3, spans_1, spans_2, spans_3,  basis_1, basis_2, basis_3,  weights_1, weights_2, weights_3, points_1, points_2, points_3, time, rhs):

    from numpy import exp
    from numpy import pi
    from numpy import sin
    from numpy import arctan2
    from numpy import cos, cosh
    from numpy import sqrt
    from numpy import zeros

    # ... sizes
    k1        = weights_1.shape[1]
    k2        = weights_2.shape[1]
    k3        = weights_3.shape[1]

    # ...
    lvalues_u = zeros((k1, k2, k3))
    # ...
    # ... build rhs
    for ie1 in range(0, ne1):
        i_span_1 = spans_1[ie1]
        for ie2 in range(0, ne2):
            i_span_2 = spans_2[ie2]
            for ie3 in range(0, ne3):
                i_span_3 = spans_3[ie3]

                for g1 in range(0, k1):
                   for g2 in range(0, k2):
                      for g3 in range(0, k3):

                            wvol  = weights_1[ie1, g1]*weights_2[ie2, g2]*weights_3[ie3, g3]

                            x     =  points_1[ie1, g1]
                            y     =  points_2[ie2, g2]
                            z     =  points_3[ie3, g3]
                            t     =  time

                            #.. Test 1
                            f     = sin(pi*t)*sin(pi*x)*sin(pi*y)*sin(pi*z)

                            lvalues_u[g1,g2,g3]  = f

                for il_1 in range(0, p1+1):
                  for il_2 in range(0, p2+1):
                    for il_3 in range(0, p3+1):

                      i1 = i_span_1 - p1 + il_1
                      i2 = i_span_2 - p2 + il_2
                      i3 = i_span_3 - p3 + il_3

                      v = 0.0
                      for g1 in range(0, k1):
                         for g2 in range(0, k2):
                            for g3 in range(0, k3):

                              bi_0 = basis_1[ie1, il_1, 0, g1] * basis_2[ie2, il_2, 0, g2] * basis_3[ie3, il_3, 0, g3]

                              wvol = weights_1[ie1, g1] * weights_2[ie2, g2] * weights_3[ie3, g3]

                              u    = lvalues_u[g1,g2,g3]
                              v   += bi_0 * u * wvol

                      rhs[i1+p1,i2+p2,i3+p3] += v   
    # ...

#==============================================================================
#---1 : Assemble rhs Poisson
@types('int', 'int', 'int', 'int', 'int', 'int', 'int[:]', 'int[:]', 'int[:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'real', 'int', 'double[:,:,:]')
def assemble_vector_ex03(ne1, ne2, ne3, p1, p2, p3, spans_1, spans_2, spans_3,  basis_1, basis_2, basis_3,  weights_1, weights_2, weights_3, points_1, points_2, points_3, Re_Pe, term, rhs):

    from numpy import exp
    from numpy import pi
    from numpy import sin
    from numpy import arctan2
    from numpy import cos, cosh
    from numpy import sqrt
    from numpy import zeros

    # ... sizes
    k1        = weights_1.shape[1]
    k2        = weights_2.shape[1]
    k3        = weights_3.shape[1]
    
    # ...
    lvalues_u = zeros((k1, k2, k3))
    # ...
    v_b1 = 1.
    v_b2 = 1.
    v_b3 = 1.
    # ... build rhs
    for ie1 in range(0, ne1):
        i_span_1 = spans_1[ie1]
        for ie2 in range(0, ne2):
            i_span_2 = spans_2[ie2]
            for ie3 in range(0, ne3):
                i_span_3 = spans_3[ie3]

                for g1 in range(0, k1):
                   for g2 in range(0, k2):
                      for g3 in range(0, k3):

                            wvol  = weights_1[ie1, g1]*weights_2[ie2, g2]*weights_3[ie3, g3]
                                        
                            x     =  points_1[ie1, g1]
                            y     =  points_2[ie2, g2]
                            z     =  points_3[ie3, g3]
                            #.. Test 1
                            if term == 0 :
                                 f     = pi*sin(pi*x)*sin(pi*y)*sin(pi*z) 
                            else :
                                  f     = (3.*pi**2*sin(pi*x)*sin(pi*y)*sin(pi*z))/Re_Pe
                                  f    += v_b1*pi*cos(pi*x)*sin(pi*y)*sin(pi*z) + v_b2*pi*sin(pi*x)*cos(pi*y)*sin(pi*z) + v_b3*pi*sin(pi*x)*sin(pi*y)*cos(pi*z)

                            lvalues_u[g1,g2,g3]  = f

                for il_1 in range(0, p1+1):
                  for il_2 in range(0, p2+1):
                    for il_3 in range(0, p3+1):

                      i1 = i_span_1 - p1 + il_1
                      i2 = i_span_2 - p2 + il_2
                      i3 = i_span_3 - p3 + il_3

                      v = 0.0
                      for g1 in range(0, k1):
                         for g2 in range(0, k2):
                            for g3 in range(0, k3):

                              bi_0 = basis_1[ie1, il_1, 0, g1] * basis_2[ie2, il_2, 0, g2] * basis_3[ie3, il_3, 0, g3]

                              wvol = weights_1[ie1, g1] * weights_2[ie2, g2] * weights_3[ie3, g3]

                              u    = lvalues_u[g1,g2,g3]
                              v   += bi_0 * u * wvol

                      rhs[i1+p1,i2+p2,i3+p3] += v   
    # ...
#==============================================================================
#---1 : Assemble l2 and H1 error norm
#==============================================================================
@types('int', 'int', 'int', 'int', 'int', 'int', 'int[:]', 'int[:]', 'int[:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:,:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:]', 'double[:,:,:]', 'real', 'double[:,:,:]')
def assemble_norm_ex01(ne1, ne2, ne3, p1, p2, p3, spans_1, spans_2, spans_3,  basis_1, basis_2, basis_3,  weights_1, weights_2, weights_3, points_1, points_2, points_3, vector_u, time, rhs):

    from numpy import exp
    from numpy import cos
    from numpy import sin
    from numpy import pi
    from numpy import sqrt
    from numpy import zeros

    # ... sizes
    k1         = weights_1.shape[1]
    k2         = weights_2.shape[1]
    k3         = weights_3.shape[1]

    # ...
    lcoeffs_u  = zeros((p1+1,p2+1,p3+1))
    # ...
    lvalues_u  = zeros((k1, k2, k3))
    lvalues_ux = zeros((k1, k2, k3))
    lvalues_uy = zeros((k1, k2, k3))
    lvalues_uz = zeros((k1, k2, k3))

    # ...
    norm_H1    = 0.
    norm_l2    = 0.
    for ie1 in range(0, ne1):
       i_span_1 = spans_1[ie1]
       for ie2 in range(0, ne2):
          i_span_2 = spans_2[ie2]
          for ie3 in range(0, ne3):
             i_span_3 = spans_3[ie3]

             lvalues_u[ : , : , :]  = 0.0
             lvalues_ux[ : , : , :] = 0.0
             lvalues_uy[ : , : , :] = 0.0
             lvalues_uz[ : , : , :] = 0.0
             lcoeffs_u[ : , : , :]  = vector_u[i_span_1 : i_span_1+p1+1, i_span_2 : i_span_2+p2+1, i_span_3 : i_span_3+p3+1]
             for il_1 in range(0, p1+1):
                for il_2 in range(0, p2+1):
                   for il_3 in range(0, p3+1):
                      coeff_u = lcoeffs_u[il_1,il_2,il_3]

                      for g1 in range(0, k1):
                        b1  = basis_1[ie1,il_1,0,g1]
                        db1 = basis_1[ie1,il_1,1,g1]
                        for g2 in range(0, k2):
                            b2  = basis_2[ie2,il_2,0,g2]
                            db2 = basis_2[ie2,il_2,1,g2]
                            for g3 in range(0, k3):
                               b3  = basis_3[ie3,il_3,0,g3]
                               db3 = basis_3[ie3,il_3,1,g3]

                               lvalues_u[g1,g2,g3]  += coeff_u*b1*b2*b3
                               lvalues_ux[g1,g2,g3] += coeff_u*db1*b2*b3
                               lvalues_uy[g1,g2,g3] += coeff_u*b1*db2*b3
                               lvalues_uz[g1,g2,g3] += coeff_u*b1*b2*db3

             v = 0.0
             w = 0.0
             for g1 in range(0, k1):
               for g2 in range(0, k2):
                 for g3 in range(0, k3):
                    wvol = weights_1[ie1, g1] * weights_2[ie2, g2] * weights_3[ie3, g3]

                    x  = points_1[ie1, g1]
                    y  = points_2[ie2, g2]
                    z  = points_3[ie3, g3]
                    t  = time

                    # ... TEST 1
                    u    = sin(pi*t)*sin(pi*x)*sin(pi*y)*sin(pi*z)
                    ux   = pi*sin(pi*t)*cos(pi*x)*sin(pi*y)*sin(pi*z)
                    uy   = pi*sin(pi*t)*sin(pi*x)*cos(pi*y)*sin(pi*z)
                    uz   = pi*sin(pi*t)*sin(pi*x)*sin(pi*y)*cos(pi*z)

                    uh  = lvalues_u[g1,g2,g3]
                    uhx = lvalues_ux[g1,g2,g3]
                    uhy = lvalues_uy[g1,g2,g3]
                    uhz = lvalues_uz[g1,g2,g3]

                    v  += ((ux-uhx)**2+(uy-uhy)**2+(uz-uhz)**2) * wvol
                    w  += (u-uh)**2 * wvol

             norm_H1 += v
             norm_l2 += w

    norm_H1 = sqrt(norm_H1)
    norm_l2 = sqrt(norm_l2)

    rhs[p1,p2,p3]   = norm_l2
    rhs[p1,p2,p3+1] = norm_H1
    # ...
