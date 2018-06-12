# evaluate each sample 

import sys
import numpy as np
import cantera as ct

# Public vars
# ===========================================

if len(sys.argv) > 1:
	p = float(sys.argv[1])*ct.one_atm
	T = float(sys.argv[2])
	COMP = sys.argv[3]
	mech = sys.argv[4]
	n_start = int(sys.argv[5])
	n_end = int(sys.argv[6])
	UFfile = sys.argv[7]

else:
	p = 1*ct.one_atm
	T = 1015.4
	COMP = 'H2:0.0009904917;O2:0.1641540720;N2:0.7728429367;H2O:0.0620124996'
	mech = 'mech/h2_li_19.xml'
	n_start = 1
	n_end = 1
	UFfile = 'mech/h2_li_UF_Weiqi.txt'

print( 'run MC for idt gradient' )
print(p/ct.one_atm, T, COMP, mech, n_start, n_end)

gas = ct.Solution(mech) 
gas.TPY = T, p, COMP

# Get equilibrium temperature for ignition break
gas.equilibrate('HP')
T_equi = gas.T

nsp = gas.n_species
nrxn = gas.n_reactions
uq_input = np.loadtxt('data/mc_in_raw.txt')

# Dimension of the input
m = uq_input.shape[1]

# Number of samples
N = uq_input.shape[0]

# perturbation factor
perturbation = 1.E-2
sigma = np.log(np.loadtxt(UFfile)[:,1])/3.

def ign_uq(factor):
	gas.set_multiplier(1.0) # reset all multipliers
	for i in range(m):
		gas.set_multiplier(factor[i],i)
	
	gas.TPY = T,p,COMP
	r = ct.IdealGasConstPressureReactor(gas)
	sim = ct.ReactorNet([r])
	t_end = 10;
	time = []
	temp = []
	while sim.time < t_end and r.T < T_equi - 10 :
		sim.step()
		time.append(sim.time)
		temp.append(r.T)

	time = np.array(time)
	temp = np.array(temp)

	diff_temp = np.diff(temp)/np.diff(time)
	
	if 'CH3OCH3' in COMP:
		ref_pos = np.argmin( np.abs(temp - T_equi + 300) )
		ign_pos = np.argmax( diff_temp[ref_pos:-1] ) + ref_pos
	else:
		ign_pos = np.argmax( diff_temp )

	if 0:
		# print( temp[-1], T_equi )
		import matplotlib.pyplot as plt
		# %matplotlib notebook

		# pressure and temperature
		fig, ax1 = plt.subplots()
		ax1.plot(time[0:-1], diff_temp, 'b-o')
		ax1.set_ylabel('dT/dt')

		ax2 = ax1.twinx()
		ax2.plot(time, temp, 'r--')
		ax2.set_ylabel('temp')
		plt.show()
	if np.max(temp) > T+10:
		idt = time[ign_pos]
	else:
		idt = 1.E8
	# Return the log10 normal ignition delay
	return np.log10(idt)

def gradient( factor ):
	
	f = ign_uq(factor)
	
	df = np.zeros((1, m))
	dk = np.zeros((m, 1))

	for i in range(m):
		if sigma[i] == 0:
			df[:,i] = 0
		else:
			dk[:] = factor[:,None]
			dk[i] = factor[i]*(1+perturbation)
			df[:,i] = (ign_uq(dk)-f)/np.log(1+perturbation)*sigma[i]

	# return idt in log10.
	return np.append(f, df)

if __name__ == '__main__':

	out_list = np.zeros((n_end-n_start+1, m+1 ))
	for i in range(n_start-1,n_end):
		# out_list[i-n_start+1, :] = gradient( np.divide(uq_input[i,:],uq_input[i,:]) )
		out_list[i-n_start+1, :] = gradient( uq_input[i,:] )

		print('Sample {:5d} idt = {:3f} s T_id = {:5f}\n'.format(i-n_start+1,out_list[i-n_start+1,0],out_list[i-n_start+1, 1]))
        
		if 0:
			import matplotlib.pyplot as plt
			print( 'perturbation = {}, S for R16, R9 = {}  {}'.format(perturbation, out_list[i-n_start+1, 16], out_list[i-n_start+1, 9] ) )
			# %matplotlib notebook
			plt.bar(range(1, nrxn+1), out_list[ i-n_start+1, 1:nrxn+1], align='center', alpha=0.5)
			plt.ylabel('df/ln(1+dk)')
			plt.show()

	np.savetxt( "data/mc_out_idt_gradient.txt"+"_"+str(n_start)+"_"+str(n_end), out_list )