import numpy as np
import math
import matplotlib.pyplot as plt
plt.switch_backend('Qt4Agg')
import re
import glob, os
import errno
import datetime
import copy


def getdate():
    now = datetime.datetime.now()
    return now.strftime("%d-%m-%y")

def getname(mydir,prefix,N,Nm,sigma,w):
    today=getdate()
    #return mydir+prefix+"_"+today+"_"+str(N)+"_"+str(Nm)+"_"+"{:.3f}".format(sigma)+"_"+"{:.2f}".format(w)+".txt"
    return mydir+prefix+"_"+today+"_"+str(N)+"_"+"{:.3f}".format(sigma)+"_"+"{:.2f}".format(w)+".txt"

def read_oldformat(filename):
    with open(filename) as fp:
        for i, line in enumerate(fp):
            if i==0:
                sigma0 = re.findall("\d+\.\d+",line)
                sigma = float(sigma0[0])
            if i==1:
                Ns0 = re.findall("\d+", line)
                Ns = int(Ns0[0])
            if i==2:
                chi0 = re.findall("\d+",line)
                chi = float(chi0[0])
            if i==3:
                d0 = re.findall("\d+\.\d+",line)
                d = float(d0[0])
            if i>6:
                break

    fdata = np.genfromtxt(filename,skip_header=7, skip_footer=1, usecols=(0,1))
    newdata = fdata[~np.isnan(fdata).any(axis=1)]
    cen = newdata[:Ns+1,1]                        # total monomer density
    #ge = newdata[Ns+1:2*Ns+2,1]                  # end monomer density
    #gm = newdata[2*Ns+2:,1]                      # middle monomer density
    #z = newdata[:Ns+1,0]                         # z coordinate
    return Ns, sigma, chi, d, cen

def write_newformat(mydir, myprefix, Ns, nk, sigma, wall, phi, Gtot,nflag):
    name = getname(mydir,myprefix, Ns, nk, sigma, abs(wall))
    file = open(name, "w+")
    file.write("# sigma: %f\n" % (sigma))
    file.write("# Ns: %d\n" % (Ns))
    file.write("# chi: %f\n" % (chi))
    file.write("# d: %f\n" % (d))
    file.write("# Nm: %d\n" % (nk))
    file.write("# z\tPHI(z)\tEM(z)\tMM(z)n=%d\n" % (int((nk-1)/2)))
    
    if (nflag == 0):
        for k in range(nk):
            file.write("%d\t%f\t%f\t%f\n" % (k, phi[k]/np.sum(phi), Gtot[nk-1][k]/sigma, Gtot[int((nk-1)/2)][k]/sigma))
    elif (nflag == 1):
        for k in range(nk):
            file.write("%d\t%f\t%f\t%f\n" % (k, phi[k], Gtot[nk-1][k], Gtot[int((nk-1)/2)][k]))
    file.close()
    return 0

def write_oldformat(mydir, myprefix, Ns, nk, sigma, wall, phi, Gtot, V0, norm, nflag):
    name = getname(mydir,myprefix, Ns, nk, sigma, abs(wall))
    file = open(name, "w+")
    file.write("sigmaexact: %f\n" % (sigma))
    file.write("Ns: %d\n" % (Ns))
    file.write("chi: %f\n" % (chi))
    file.write("d: %f\n" % (d))
    file.write("Nm: %d\n" % (nk))
    file.write("\n")
    file.write("z\tDichteverteilung(z)=cen(z)\n")

    if (nflag == 0):
        for k in range(nk):
            file.write("%d\t%f\n" % (k, phi[k]/np.sum(phi)))
    elif (nflag == 1):
        for k in range(nk):
            file.write("%d\t%f\n" % (k, phi[k]))

    file.write("z\tEndmonomerverteilung(z)\n")
    
    if (nflag == 0):
        for k in range(nk):
            file.write("%d\t%f\n" % (k, Gtot[nk-1][k]/sigma))
    elif (nflag == 1):
        for k in range(nk):
            file.write("%d\t%f\n" % (k, Gtot[nk-1][k]))

    file.write("Monomer Nr.: %d\n" % (int((nk-1)/2)))
    file.write("z\t Verteilung mittleres_Monomer\n")

    if (nflag == 0):
        for k in range(nk):
            file.write("%d\t%f\n" % (k,Gtot[int((nk-1)/2)][k]/sigma))
    elif (nflag == 1):
        for k in range(nk):
            file.write("%d\t%f\n" % (k,Gtot[int((nk-1)/2)][k]))

    file.write("\n")
    file.write("V0: %f norm: %f\n" % (V0, norm))
    file.close()
    return 0

def write_fe(mydir, myprefix,Ns, nk, sigma, wall,Gtot):
    name = getname(mydir,myprefix, Ns, nk, sigma, abs(wall))
    file = open(name, "w+")
    file.write("# sigma: %f\n" % (sigma))
    file.write("# Ns: %d\n" % (Ns))
    file.write("# chi: %f\n" % (chik))
    file.write("# d: %f\n" % (d))
    file.write("# Nm: %d\n" % (nk))
    file.write("# z\tFreeEnergy(log_of_endmon)\n")
    for k in range(nk):
        #file.write("%d\t%f\n" % (k, -1.0*np.log(Gtot[nk-1][k])))
        file.write("%d\t%f\n" %(k, -1.0*np.log(Gtot[nk-1][k]) if Gtot[nk-1][k] != 0 else 0.0))
    file.close()
    return 0

def write_adsorption(mydir,Ns,sigma, Nk, dk_werte, endmon_wand, vflag):
    if (vflag == 0):
        myprefix="poly_endm_vs_Nk_"
    elif (vflag == 1):
        myprefix="poly_endm_vs_d_"

    filename=mydir+myprefix+str(Ns)+"_"+"{:.3f}".format(sigma)+".txt"
    file = open(filename, "w+")
    
    if (vflag == 0):
        file.write("Nk")
        for dk in dk_werte:
            file.write(" %f" % dk)
        file.write("\n")
        file.write("\n")
        file.write("\n")
        for t, nk in enumerate(Nk):
            file.write("%d" % nk)
            for i, dk in enumerate(dk_werte):
                file.write(" %f" % endmon_wand[t][i])
            file.write("\n")
        file.close()    

    elif (vflag == 1):
        file.write("d")
        for nk in Nk:
            file.write(" %d" % nk)
        file.write("\n")
        file.write("\n")
        file.write("\n")
        for i, dk in enumerate(dk_werte):
            file.write("%f" % dk)
            for t, nk in enumerate(Nk):
                file.write(" %f" % endmon_wand[t][i])
            file.write("\n")
        file.close()
    return 0


def Vakt(r,ww):
    # r - monomer density = phi(z)
    # type = Auswahlparameter fuer Art des Potentials
    # ww - adsorption energy = chi parameter
    return -np.log(1.0-r) - 2.0*ww*r


#targetdir="./polydisp/"
#targetdir="./"
#os.chdir(targetdir)

#if not os.path.exists("mb"):
#    os.mkdir("mb")
#if not os.path.exists("nnmb"):
#        os.mkdir("nnmb")
#if not os.path.exists("fe"):
#        os.mkdir("fe")
#if not os.path.exists("MB"):
#        os.mkdir("MB")
#if not os.path.exists("NNMB"):
#        os.mkdir("NNMB")


      
#Ns - Kettenlaenge = Anzahl Monomere
#Nz - Gittergroesse
#sig - Graftingdichte
#sigma - ??   zwischen 0,025, 0,05  0,1, 0,2 0,3 0,4 0,5 (0,7)
#varV - ?
#Vakt - Interaction model
#V0 - ?
#cen[] - Polymerdichte
#chi - 
#w0, w1 - Wahrscheinlichkeit fuer Schritt in z, oder x-y richtung
#norm1 - ??
#F - ? 
#G - Green-Funktionen? 
#jacob - Jacobimatrix - Was beschreibt diese?
#epsi - Schrittlaenge fuer Iteration;
#wall - 0 bis man 99% der endmon an die wand kriegt (fuer sigma bis 0,4)
      
#Endmonomere an Ort 0 (als array fuer versch wall, sigma etc speichern und in extra datei ausgeben)
#-> Abhaengigkeit von Sigma bei festem Wall, versch wall bei einem Sigma usw
#Fluktuation des mittleren Monomers (quadr. Abweichung vom Mittelwert)
#verschiedene Sigma für feste Ns, wall
#versch. wall für feste Ns, Sigma
#versch Ns, dann z/Ns auftragen
#mittlere Höhe der Bürste (als Integral) und maximale Höhe (ab Wert wo Dichte < als bestimmtes epsilon etc)
      
#evtl: zustandssumme aus Summe von G ueber z -> - ln Z = F (sinnvoll erstmal fuer neutrale wand)
#chemisches potential je flaeche: my = dF/dsigma = F(sigma+epsilon)-F(sigma) / epsilon

Ns_werte=[256]
#Ns_werte=[192]
#sigma_werte=[0.025,0.050,0.075,0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800,0.900]
sigma_werte=[0.300]
#d_werte=[0.00]
d_werte=np.arange(0.0,145.0,step=0.5)
print(d_werte)


for u, Ns in enumerate(Ns_werte):
    for v, sigma in enumerate(sigma_werte):
        for p, d in enumerate(d_werte):

            print(Ns, sigma, d)
            wall = -float(2*d)  # wall adsorption (dont know why i doubled d here)
            Nz = Ns+3           # size of lattice
            tmax = 300          # max. number of iterations
            epsi = 0.00001      # Abbruchbedingung der Schleife 
            w0 = 4.0/6.0        # transition probability on
            w1 = 1.0/6.0        # simple cubic lattice
            norm1 = 1.0
            norm = 0.0
            V0 = 0.0
            chi = 0.0
            
            # greens functions
            #G1 - forward, G2 - backward, Gtot = full
            G1 = np.zeros((Ns,Ns+2))
            G2 = np.zeros((Ns,Ns+2))
            Gtot = np.zeros((Ns,Ns+1))
            #phi=np.zeros(Ns+1)

            varV=np.zeros((Ns+1,Ns+1))
            cen=np.zeros(Ns+1)
            jacob=np.zeros((Ns,Ns))
            func=np.zeros(Ns)
            F=np.zeros((Ns,Ns+1))
            #pivot=np.zeros(Nz-8)

            #Randbedingungen: das erste Monomer ist gegraftet und immer an der Wand, also gilt:
            #G1[0][1]=float(1)
            #daher Umstellen: arrayindex s (erster) startet bei 0 bis Ns-1
            #arrayindex z (zweiter) startet auch bei 0 bis Ns-1 aber z=0 ist die Wand, daher ist 1 das erste Monomer

            #endmon_wand=np.zeros((N_anz,d_anz))

            # outer loop of iterative method
            for t in range(tmax):
                print("Iteration-Nr. %d" % t)
                
                # in beginning norm1 = 1.0, so this is zero
                #add on top of potential if its values are too small
                V0 += float(int(np.log10(norm1)/150.0))
                
                # create potential from varied densities, cen[j] = old concentration of monomers
                # 2nd index is spatial pos., cen[j] should not be >= 1, otherwise nan
                # vary potential on diagonal
                for h in range(Ns+1):
                    varV[h,1:Ns+1] = np.exp(-1.*(Vakt(cen[1:Ns+1],float(chi)/10.0)+V0))
                    varV[h,h] = np.exp(-1.*(Vakt(cen[h]+epsi,float(chi)/10.0)+V0))

                # set matrix elements to zero 
                F.fill(0)

                # inner loop
                for h in range(Ns+1):
                    
                    # reset all greens functions to zero before new run
                    G1.fill(0)
                    G2.fill(0)
                    Gtot.fill(0)
                    
                    # first, calculate forward G, i.e. G1
                    # wall is at 0, so monomers can only be at 1,2,3...
                    # boundary condition: at wall G1 = 1 for first monomer (grafting)
                    G1[0,1]=float(1)
                    
                    # calculate full G1 from forward moving condition
                    # connects previous steps with step probabilities
                    # and exposure to potential V
                    # (pg. 111, Polymers at Interfaces)
                    for n in range(1,Ns):
                        G1[n,1:n+2] = w1*(G1[n-1,:n+1]+G1[n-1,2:n+3])+w0*G1[n-1,1:n+2]
                        G1[n,1:n+2] *= varV[h,1:n+2]                

                    # apply additiponal potential at surface (attractive)
                    # end-modification =  only if terminal monomer (Ns-1) is at wall (1)
                    G1[Ns-1,1] *= np.exp(-wall)

                    # norm1 = sum of all probabilities for the terminal monomer (sum over position)
                    # connection condition for backward G (=G2) 
                    # together, G1*G2 comes out to be sigma (grafting density)
                    # norm = normation in grafting density units
                    norm1 = np.sum(G1[Ns-1,:])
                    norm = sigma/norm1

                    # G2 (backwards propagating) begins with terminal monomer
                    # set inital values of G2
                    # if at wall, gain of adsorption energy
                    G2[Ns-1,1:Ns+1] = norm*varV[h,1:Ns+1]
                    G2[Ns-1,1] *= np.exp(-wall)
                    
                    # calculate all remaining values of G2 (from last to first monomer)
                    # monomer k can only reach up to position k (no need to calculate pos. j for k if j>k)
                    for n in range(Ns-1):
                        G2[Ns-2-n,1:Ns-n] = w1*(G2[Ns-1-n,:Ns-n-1]+G2[Ns-1-n,2:Ns-n+1])+w0*G2[Ns-1-n,1:Ns-n]
                        G2[Ns-2-n,1:Ns-n] *= varV[h,1:Ns-n]
                
                    # obtain full Gtot by multiplication of G1*G2
                    # (see Polymers at Interfaces for equation)
                    Gtot = np.multiply(G1,G2)
                    Gtot[:Ns,1:Ns+1] /= varV[h,1:Ns+1] 
                   
                    # energy gain at all appears twice (in G1 and G2)
                    # it has to be eliminated once to get the correct Gtot
                    Gtot[Ns-1,1] *= np.exp(wall)

                    # generate matrix elements of A for Ax=b
                    # A = jacobi matrix (1st derivative at old positions)
                    # b = function value at old position
                    # x = difference old to new position
                    # Solution x will be stored in func (is on input b and on output x)
                    # sum over segment index s
                    # this gives density phi(z) under premise that there is some fluctuation for the h-entries
                    # new density minus old to see the difference (cen[j] contains the old values of phi(z)) 
                    for j in range(1,Ns+1):
                        F[j-1,h] += np.sum(Gtot[:Ns,j]) - cen[j]
                        if j == h:
                            F[j-1,h] -= epsi

                # here, inner loop over h is ending
                # for iterative method we have to move towards self-consistent solution!
                # next, use multidimensional Newton-Raphson method to find root of matrix equation
                # here: jacob * x = func
                # we search the root of -func = -F[j][0]
                    
                # first column of F is target vector b: Ax = b
                # why this F with this index? => it is the unchanged one (no epsi subtracted)
                # func[j] = Difference between Gtot[j+1] and cen[j] (new and old)
                # CAREFUL: now the spatial coordinate is shifted (such that arry begins at 0)
                
                # CAREFUL: Python uses call by reference
                # make a separate copy of cen in memory
                func[:Ns] = F[:Ns,0] 
                cen2 = copy.deepcopy(cen)

                # generate jacobian matrix (matrix of first derivatives with resp. to x)
                for j in range(Ns):
                    jacob[j,:Ns] = (F[j,1:Ns+1]-F[j,0])/epsi    
    
                # solve matrix equation Ax=b 
                # store solution again in func (on input, func=b, on output func=x)
                x = np.linalg.solve(jacob, func)
                func = copy.deepcopy(x)

                # parameter myhelp is a measure for convergence
                # parameter a is control parameter to prevent overshoots during updating
                myhelp = np.sum(np.absolute(func))
                myhelp *= 1.5
                alpha = 1.0/(myhelp+1.5)

                # next, calculate new density field
                # avoid densities > 1, otherwise VarV will become nan becuase of log(1-chi)
                # ensure positive densities at the same time
                # reminder: func is now solution x of above matrix equation

                cen[cen <= 1.0 + alpha*np.roll(np.append(func,[0]),1)] += -alpha*np.roll(np.append(func,[0]),1)
                cen[cen < 0.] = np.abs(cen[cen < 0.])
                cen[0] = 0.
        
                norm = np.sum(cen)
                print("momentary myhelp: %f" % myhelp)
        
                # Stop iteration if myhelp is sufficiently small
                if myhelp < 0.00001:
                    print("Iteration break.")
                    print("momentary myhelp: %f" % myhelp)            
                    print("momentary iteration: %d" % t)
            
                    # write brush profile to file
                    write_oldformat("./", "b", Ns, Ns, sigma, wall, cen, Gtot, V0, norm, 1)

                    # write also adsorption curves?
                    break
            # here, outer loop over t is ending


    #phi=np.sum(Gtot, axis=0)        
    #endmon_wand[t,i]=Gtot[nk-1,1]/sigma
    # write out files with density profiles etc. 
    #write_oldformat("./mb/", "mb", Ns, nk, sigma, wall, phi, Gtot, V0, norm, 0)
    #write_oldformat("./nnmb/", "nnmb", Ns, nk, sigma, wall, phi, Gtot, V0, norm, 1)
    #write_newformat("./MB/", "MB",Ns, nk, sigma, wall, phi, Gtot, 0)
    #write_newformat("./NNMB/", "NNMB",Ns, nk, sigma, wall, phi, Gtot, 1)
    #write_fe("./fe/", "fe", Ns, nk, sigma, wall, Gtot)
        
    #print("Adsorptionskurven schreiben")
    #write_adsorption("./",Ns,sigma, Nk, dk_werte, endmon_wand, 0)
    #write_adsorption("./",Ns,sigma, Nk, dk_werte, endmon_wand, 1)
    #print("Program end")
