from simplines import compile_kernel, apply_dirichlet

from simplines import SplineSpace
from simplines import TensorSpace
from simplines import StencilMatrix
from simplines import StencilVector
from simplines import pyccel_sol_field_3d
from simplines import prolongation_matrix

#---In Poisson equation
from gallery_section_04 import assemble_vector_ex01
from gallery_section_04 import assemble_vector_ex02
from gallery_section_04 import assemble_vector_ex03

from gallery_section_04 import assemble_matrix_ex01
from gallery_section_04 import assemble_matrix_ex02
from gallery_section_04 import assemble_matrix_ex03
from gallery_section_04 import assemble_matrix_ex04
from gallery_section_04 import assemble_norm_ex01

assemble_rhs              = compile_kernel(assemble_vector_ex01, arity=1)
assemble_in_rhs           = compile_kernel(assemble_vector_ex02, arity=1)
assemble_f_rhs            = compile_kernel(assemble_vector_ex03, arity=1)

assemble_stiffness2D      = compile_kernel(assemble_matrix_ex01, arity=2)
assemble_mass2D           = compile_kernel(assemble_matrix_ex02, arity=2)
assemble_advec2D          = compile_kernel(assemble_matrix_ex03, arity=2)
assemble_mass2D_in_center = compile_kernel(assemble_matrix_ex04, arity=2)

assemble_norm_l2          = compile_kernel(assemble_norm_ex01, arity=1)

#---mae : In one dimension
from gallery_section_04 import assemble_stiffnessmatrix1D
from gallery_section_04 import assemble_massmatrix1D

assemble_stiffness1D = compile_kernel( assemble_stiffnessmatrix1D, arity=2)
assemble_mass1D      = compile_kernel( assemble_massmatrix1D, arity=2)

#from matplotlib.pyplot import plot, show
import matplotlib.pyplot            as     plt
from   mpl_toolkits.axes_grid1      import make_axes_locatable
from   mpl_toolkits.mplot3d         import axes3d
from   matplotlib                   import cm
from   mpl_toolkits.mplot3d.axes3d  import get_test_data
from   matplotlib.ticker            import LinearLocator, FormatStrFormatter
#..
from   scipy.sparse                 import kron
from   scipy.sparse                 import csr_matrix
from   scipy.sparse                 import csc_matrix, linalg as sla
from   numpy                        import zeros, linalg, asarray
from   numpy                        import cos, sin, pi, exp, sqrt, arctan2
from   tabulate                     import tabulate
import numpy                        as     np
import timeit
import time

from tabulate import tabulate


#==============================================================================
#.......Picard BFO ALGORITHM
#==============================================================================
class Convection_diffusion(object):
    def __init__(self,V1, V2, V3, V, Re_Pe, dt):

       mass                = assemble_mass2D(V)
       self.m              = mass.tosparse()
       # ... after apllying direchlet
       # ... comlpute matrix E : mass matrix
       mass                = apply_dirichlet(V, mass)
       M                   = mass.tosparse()
       np.savetxt('E.m', M.toarray(), fmt='%.2e')

       # ...
       # ... comlpute matrix A : stiffness matrix
       stiffness           = assemble_stiffness2D(V, value = [Re_Pe])
       stiffness           = apply_dirichlet(V, stiffness)
       S                   = stiffness.tosparse()
       np.savetxt('A.m', S.toarray(), fmt='%.2e')
       # ...
       #adv_mat             = assemble_advec2D(V)
       #adv_mat             = apply_dirichlet(V, adv_mat)
       #A                   = adv_mat.tosparse()
       # ...Global matrix
       G_mat               = M + dt*S

       # ... 
       self.lu             = sla.splu(csc_matrix(G_mat))       
       self.V              = V
       self.Re_Pe          = Re_Pe
       self.dt             = dt

       
    def Proj(self):
             
       sol_in_st           = StencilVector(self.V.vector_space)
       #--Assembles a right hand side of Poisson equation
       rhs                 = assemble_in_rhs( self.V)
       b                   = rhs.toarray()*sin(pi*0.)
       # ...
       lum                 = sla.splu(csc_matrix(self.m))
       x                   = lum.solve(b)
       x                   = x.reshape(self.V.nbasis)
       sol_in_st.from_array(self.V, x)

       #--Computes error l2 and H1
       Norm                = assemble_norm_l2(self.V, fields=[sol_in_st], value = [0.])
       norm                = Norm.toarray()
       l2_norm             = norm[0]
       H1_norm             = norm[1]
       return sol_in_st, x, l2_norm, H1_norm

    def compute_matrix_B_C(self, Vc):
       #--Assembles a right hand side of Poisson equation
       rhs                 = assemble_in_rhs( self.V )
       b0                  = rhs.toarray()
       # ...
       rhs                 = assemble_f_rhs( self.V , value = [self.Re_Pe])
       b1                  = rhs.toarray()
       # ...
       B                   = np.zeros((self.V.nbasis[0]*self.V.nbasis[1]*self.V.nbasis[2],2)) 
       B[:,0]              = b0[:]
       B[:,1]              = b1[:]
       
       # ... computes C matrix for outputs
       mass                = assemble_mass2D_in_center(Vc)
       m_center            = mass.tosparse()
       # ...
       C                   = np.zeros((2, self.V.nbasis[0]*self.V.nbasis[1]*self.V.nbasis[2])) 
       C[0,:]              = np.ones(self.V.nbasis[0]*self.V.nbasis[1]*self.V.nbasis[2]).dot(m_center.toarray())
       C[1,:]              = np.ones(self.V.nbasis[0]*self.V.nbasis[1]*self.V.nbasis[2]).dot(self.m.toarray())
       return B, C
       
    def solve(self, u_n = None, time = None, B = None):
             
       sol_in_st               = StencilVector(self.V.vector_space)

       if B is None:

           #--Assembles a right hand side of Poisson equation
           rhs                 = assemble_rhs( self.V, fields=[u_n], value = [self.Re_Pe, time, self.dt])
           rhs                 = apply_dirichlet(self.V, rhs)
           b                   = rhs.toarray()

           # ...
           x                   = self.lu.solve(b)

       else :
           # ... The inputs of underlying model 
           u                   = np.asarray([pi*cos(pi*t), sin(pi*t)])
           # ...
           b                   = B.dot(u)
           # ...
           x                   = u_n.toarray() + dt*self.lu.solve(b)
           # ... The outputs of underlying model
           y                   = C.dot(x)
       # ...
       x                   = x.reshape(self.V.nbasis)
       sol_in_st.from_array(V, x)
       #--Computes error l2 and H1
       Norm                = assemble_norm_l2(self.V, fields=[sol_in_st], value = [time])
       norm                = Norm.toarray()
       l2_norm             = norm[0]
       H1_norm             = norm[1]
       return sol_in_st, x, l2_norm, H1_norm

       
# ...
nbr_freedm    = 202140
degree        = 2
nelements_x   = 12 - degree
nelements_y   = 15 - degree
nelements_z   = 1123 - degree

'''
if False :
 nbr_freedm = 202140    
 p=10 
 q=10 
 n=10 
 while n*p*q <=nbr_freedm: 
   while n*p*q <=nbr_freedm: 
     while n*p*q <=nbr_freedm: 
         if n*p*q == nbr_freedm : 
              print(n,p,q) 
         q +=1 
     q = 10 
     p+=1 
   p = 10 
   n +=1
   nelemenent_x = n-degree
   nelemenent_y = p-degree
   nelemenent_z = q-degree
'''

'''   
# ... test
degree        = 2
nelements_x   = 8
nelements_y   = 8
nelements_z   = 8
'''
n             = (nelements_x + degree) * (nelements_y + degree) *(nelements_z + degree)
print('total number of freedom = ', n)
# ...
Re_Pe       = 1.
dt          = 1e-4
t_f         = 0.02
t           = 0. # initial time

#----------------------
#..... Initialisation and computing optimal mapping for 16*16
#----------------------
grid_x  = np.linspace( 0., 2., nelements_x+1)
grid_y  = np.linspace(-2., 2., nelements_y+1)
grid_z  = np.linspace(-3., 3., nelements_z+1)

# create the spline space for each direction
V1   = SplineSpace(degree=degree, nelements= nelements_x, grid= grid_x, nderiv = 1)
V2   = SplineSpace(degree=degree, nelements= nelements_y, grid= grid_y, nderiv = 1)
V3   = SplineSpace(degree=degree, nelements= nelements_z, grid= grid_z, nderiv = 1)
V    = TensorSpace(V1, V2, V3)

# ... for outputes
# create the spline space for each direction
Vc1   = SplineSpace(degree=degree, nelements= nelements_x, grid= grid_x, nderiv = 1, sharing_grid = np.linspace( 0.5, 1.5, nelements_x+1))
Vc2   = SplineSpace(degree=degree, nelements= nelements_y, grid= grid_y, nderiv = 1, sharing_grid = np.linspace( -0.5, 0.5, nelements_y+1))
Vc3   = SplineSpace(degree=degree, nelements= nelements_z, grid= grid_z, nderiv = 1, sharing_grid = np.linspace( 0.5, 0.5, nelements_z+1))
Vc    = TensorSpace(Vc1, Vc2, Vc3)
# ... Initialization
Cd   = Convection_diffusion(V1, V2, V3, V, Re_Pe, dt)

# ... Computes the L2 projection of the Initial Solution 
u_0, xuh_0, l2_error, H1_error = Cd.Proj()

# ... compute matrix B
B, C = Cd.compute_matrix_B_C(Vc)
np.savetxt('B.m', B, fmt='%.2e')
np.savetxt('C.m', C, fmt='%.2e')

print('l2_error = {}  H1_error = {} at the time = {}'.format(l2_error, H1_error, t))


t   += dt
while (t < t_f):
	u_0, xuh, l2_error, H1_error = Cd.solve( u_n = u_0, time = t, B=B)
	print('l2_error = {}  H1_error = {} at the time = {}\n\n'.format(l2_error, H1_error, t))
	t   += dt

# # ........................................................
if True :
	nbpts = 50
	#---Compute a solution
	u, ux, uy, uz, X, Y, Z      = pyccel_sol_field_3d((nbpts, nbpts, nbpts),  xuh, V.knots, V.degree)
		
	np.save('mesh_X.npy',X )
	np.save('mesh_Y.npy',Y )
	np.save('mesh_Z.npy',Z )
	np.save('ux.npy',ux )
	#np.save('uy.npy',uy )
	#np.save('uz.npy',uz )
	np.save('un_sol.npy',u)
	
	figtitle  = 'Solution and domain'
	fig, axes = plt.subplots( 1, 3, figsize=[14,12], gridspec_kw={'width_ratios': [2,2, 2]} , num=figtitle )
	for ax in axes:
	   ax.set_aspect('equal')

	for i in range(nbpts):
	     phidx = X[:,i,25]
	     phidy = Y[:,i,25]
	  
	     axes[0].plot(phidx, phidy, '-k', linewidth = 1.)
	for i in range(nbpts):
	     phidx = X[i,:,25]
	     phidy = Y[i,:,25]
	  
	     axes[0].plot(phidx, phidy, '-k', linewidth = 1.)
	#axes[0].axis('off')
	axes[0].margins(0,0)

	for i in range(nbpts):
	     phidx = X[:,40,i]
	     phidy = Z[:,40,i]
	  
	     axes[1].plot(phidx, phidy, '-k', linewidth = 1.)
	for i in range(nbpts):
	     phidx = X[i,40,:]
	     phidy = Z[i,40,:]
	  
	     axes[1].plot(phidx, phidy, '-k', linewidth = 1.)
	#axes[1].axis('off')
	axes[1].margins(0,0)

	for i in range(nbpts):
	     phidx = Y[40,:,i]
	     phidy = Z[40,:,i]
	  
	     axes[2].plot(phidx, phidy, '-k', linewidth = 1.)
	for i in range(nbpts):
	     phidx = Y[40,i,:]
	     phidy = Z[40,i,:]
	  
	     axes[2].plot(phidx, phidy, '-k', linewidth = 1.)
	#axes[2].axis('off')
	axes[2].margins(0,0)
	  
	plt.show()
