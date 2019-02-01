import numpy as np
from sympy import Matrix,lcm
import matplotlib.pyplot as plt
from fractions import Fraction
from CRN.ReactionSystem import MassActionSystem, arraypow
from timeit import default_timer
from RBM.RBM import RBM
from LLL import basis_reduce

def Lyapunov(X,Y):
	return np.sum(Y - X + (X*np.log(X/Y)))

random_seed =6
np.random.seed(random_seed)

# Number of timesteps
ts = 20000
t =  25

A = [[2,1,0],[0,1,2]]
X = [0.2,0.42,0.38]

A = np.array(A)
X = np.array(X)
theta = np.ones(A.shape[0])
species =[]
for i in xrange(A.shape[0]):
	species.append("\\theta_"+str(i+1))

reactions =[]
rates =[]
St = theta.shape[0] # Number of theta species

for i in xrange(A.shape[1]):
	column = A.T[i].tolist()
	complexl = column 
	complexr = [0]*St
	reactions.append([complexl[:],complexr[:]])
	reactions.append([complexr,complexl])
	rates = rates +[1.,X[i]]

reactions = np.array(reactions)
rates = np.array(rates)
system = MassActionSystem(reactions,rates,species)

print "The Reactions are:"
print system.display_reactions()
system.set_concentrations(theta)

output = system.run(t=t,ts=ts, plot=False)
theta_final = system.current_concentrations()
thetas = system.concentrations_profile()
dthetas = system.dydt_profile()

D_infinity = Lyapunov(X,arraypow(theta_final,A))
Y = np.array(map(lambda f: arraypow(f,A), thetas))
D_profile = np.array(map(lambda y: Lyapunov(X,y), Y))
temp = dthetas**2/thetas
D_dot = temp[:,0] + temp[:,1]

t_index = np.linspace(0, t, ts)
plt.plot(t_index[:-1], (D_dot/(D_profile-D_infinity))[:-1] )
plt.legend(loc='best')
plt.xlabel('Time (t)')
plt.ylabel('-D_dot/(D - D_infinity)')
plt.grid()
plt.show() 