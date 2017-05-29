# MassActionSystem.py
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def arraypow(x,A):
	#Computes(\ theta^A)
	return np.prod(x**(A.T),axis=-1)

def odes(y,t,R,k):
	l = R[:,0] # l, left complex
	r = R[:,1] # r, right complex
	return ((r-l).T).dot(arraypow(y,l.T)*k) 


class MassActionSystem(object):
	"""docstring for MassActionSystem"""
	def __init__(self, reactions, rates, species =None):
		# super(MassActionSystem, self).__init__()
		"""Each reaction in reactions is 
		a pair of complex (l,r) in l --> r"""
		reactions = np.array(reactions)
		rates = np.array(rates)*1.0
		self.reactions = reactions
		self.rates = rates
		if species is None:
			species =['S_'+str(i+1) for i in xrange(reactions.shape[2])]
		self.species = species
#### TODO: Add a method to add more reactions

	def display_reactions(self):
		"""
		Puts together the set of reations in a string and returns it
		"""
		reaction_set =""
		for i in xrange(self.reactions.shape[0]):
			l = self.reactions[i][0]
			r = self.reactions[i][0]
			S = self.species
			left_complex =""
			right_complex=""
			for j in xrange(len(S)):
				left_complex+= " " + str(l[j]) +S[j] + " +"
				right_complex+=" " + str(r[j]) +S[j] + " +"
			reaction = left_complex[:-1] + " ------> " + right_complex[:-1] + "rate: "+ str(self.rates[i]) + "\n\n"
			reaction_set +=reaction 
		return reaction_set[:-1]

	def set_concentrations(self,concentrations):
		"""
		Sets the concentration of the species
		"""
		concentrations = 1.0*np.array(concentrations)
		if concentrations.shape[0] == self.reactions.shape[2]:
			self.concentrations = concentrations
		else:
			print "Wrong Initialization: Shapes of concentrations and reactions dont match ", concentrations.shape[0], "!=", self.reactions.shape[1]
		return self.concentrations

	def dydt(self):
		"""
		Returns rates of concentrations change at the current concentration
		"""
		return odes(self.concentrations,0,self.reactions,self.rates)

	def current_concentrations(self):
		"""
		Returns the concentrations after the latest run
		"""
		return self.concentrations

	def run(self,t=10000000,ts=100000000,plot=False):
		"""
		Run it for time t with ts number of time steps
		Outputs the concentration profile for the run
		default values are 100000 and 1000000
		"""	
		t_index = np.linspace(0, t, ts)
		y = self.concentrations
		output = odeint(odes, y, t_index, args= (self.reactions,self.rates))
		self.concentrations = output[-1,:]
		if plot:
			for i in xrange(output.shape[1]):
				label = self.species[i]
				plt.plot(t_index, output[:, i], label=label)
			plt.legend(loc='best')
			plt.xlabel('t')
			plt.grid()
			plt.show()
		return output

	def run_till(self,delta_time=0.001,eps = 10**-12,every=10, stop=100000):
		"""
		Runs till the gradients are less than eps
		returns time taken with array of concentrations
		default values of time_step and eps are 10**-3 and 10**-12 
		"""
		y = self.concentrations
		output = [y]
		t=0
	 	t_index = np.linspace(0,every,int(every/(delta_time)))
		if every>10*delta_time:
			 while True:
			 	gy = odes(y,t,self.reactions,self.rates)
				if (np.abs(gy) <eps).all():
					break
			 	o = odeint(odes, y, t_index, args= (self.reactions,self.rates))
			 	y = o[-1,:]
			 	t+=every
			 	output = output + o[1:]
			 	if t>stop:
			 		break 
		else:
			while True:
				gy = odes(y,t,self.reactions,self.rates)
				if (np.abs(gy) <eps).all():
					break
				y+= gy*delta_time
				t+=delta_time
				output.append(y)
				if t>stop:
					break

		output = np.array(output)
		self.concentrations =y
		return output,t
		
	

def main():
	reactions = [[[1,0],[0,1]], [[0,1],[1,0]]]
	rates = [1,1]
	system = MassActionSystem(reactions,rates)
	system.set_concentrations([2.,4.])
	system.run()
	print system.current_concentrations()
	print system.dydt()
	return

if __name__ == '__main__':
	main()