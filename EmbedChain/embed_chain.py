import numpy as np
import math
import matplotlib.pyplot as plt
plt.switch_backend('Qt4Agg')
import re
import glob, os
import errno
import datetime
from scipy.signal import argrelextrema

def getdate():
    now = datetime.datetime.now()
    return now.strftime("%d-%m-%y")

def getname(mydir,prefix,N,Nm,sigma,w):
    today=getdate()
    return mydir+prefix+"_"+today+"_"+str(N)+"_"+str(Nm)+"_"+"{:.3f}".format(sigma)+"_"+"{:.2f}".format(w)+".txt"

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
    file.write("# chi: %f\n" % (chik))
    file.write("# d: %f\n" % (d))
    file.write("# Nm: %d\n" % (nk))
    file.write("# z\tPHI(z)\tEM(z)\tMM(z)n=%d\n" % (int((nk-1)/2)))
    
    if (nflag == 0):
        for k in range(nk):
            file.write("%d\t%f\t%f\t%f\n" % (k, phi[k]/np.sum(phi), Gtot[nk-1,k]/sigma, Gtot[int((nk-1)/2),k]/sigma))
    elif (nflag == 1):
        for k in range(nk):
            file.write("%d\t%f\t%f\t%f\n" % (k, phi[k], Gtot[nk-1,k], Gtot[int((nk-1)/2),k]))
    file.close()
    return 0

def write_oldformat(mydir, myprefix, Ns, nk, sigma, wall, phi, Gtot, V0, norm, nflag):
    name = getname(mydir,myprefix, Ns, nk, sigma, abs(wall))
    file = open(name, "w+")
    file.write("sigmaexact: %f\n" % (sigma))
    file.write("Ns: %d\n" % (Ns))
    file.write("chi: %f\n" % (chik))
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
            file.write("%d\t%f\n" % (k, Gtot[nk-1,k]/sigma))
    elif (nflag == 1):
        for k in range(nk):
            file.write("%d\t%f\n" % (k, Gtot[nk-1,k]))

    file.write("Monomer Nr.: %d\n" % (int((nk-1)/2)))
    file.write("z\t Verteilung mittleres_Monomer\n")

    if (nflag == 0):
        for k in range(nk):
            file.write("%d\t%f\n" % (k,Gtot[int((nk-1)/2),k]/sigma))
    elif (nflag == 1):
        for k in range(nk):
            file.write("%d\t%f\n" % (k,Gtot[int((nk-1)/2),k]))

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
        #file.write("%d\t%f\n" % (k, -1.0*np.log(Gtot[nk-1,k])))
        file.write("%d\t%f\n" %(k, -1.0*np.log(Gtot[nk-1,k]) if Gtot[nk-1,k] != 0 else 0.0))
    file.close()
    return 0


def write_fe_norm1(mydir, myprefix,Ns, nk, sigma, wall,Gtot):
    name = getname(mydir,myprefix, Ns, nk, sigma, abs(wall))
    file = open(name, "w+")
    file.write("# sigma: %f\n" % (sigma))
    file.write("# Ns: %d\n" % (Ns))
    file.write("# chi: %f\n" % (chik))
    file.write("# d: %f\n" % (d))
    file.write("# Nm: %d\n" % (nk))
    file.write("# z\tFreeEnergy(log_of_endmon)\n")
    for k in range(nk):
        #file.write("%d\t%f\n" % (k, -1.0*np.log(Gtot[nk-1,k])))
        file.write("%d\t%f\n" %(k, -1.0*np.log(Gtot[nk-1,k]/np.sum(Gtot[nk-1,:])) if Gtot[nk-1,k] != 0 else 0.0))
    file.close()
    return 0




def write_adsorption(mydir,Ns,sigma, Nk, dk_werte, endmon_wand, vflag, phrase):
    if (vflag == 0):
        myprefix="poly_"+phrase+"_vs_Nk_"
    elif (vflag == 1):
        myprefix="poly_"+phrase+"_vs_d_"

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
                file.write(" %f" % endmon_wand[t,i])
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
                file.write(" %f" % endmon_wand[t,i])
            file.write("\n")
        file.close()
    return 0



def write_zaverage(mydir,Ns,sigma, Nk, dk_werte, z_end_avg, vflag):
    if (vflag == 0):
        myprefix="fullcurve_vs_Nk_"
    elif (vflag == 1):
        myprefix="fullcurve_vs_d_"

    filename=mydir+myprefix+str(Ns)+"_"+"{:.3f}".format(sigma)+".txt"
    file = open(filename, "w+")

    if (vflag == 0):
        file.write("# Nk")
        for dk in dk_werte:
            file.write(" %f" % dk*2)
        file.write("\n")
        #file.write("\n")
        #file.write("\n")
        for t, nk in enumerate(Nk):
            file.write("%d" % nk)
            for i, dk in enumerate(dk_werte):
                file.write(" %f" % z_end_avg[t,i])
            file.write("\n")
        file.close()

    elif (vflag == 1):
        file.write("# epsilon")
        for nk in Nk:
            file.write(" %d" % nk)
        file.write("\n")
        #file.write("\n")
        #file.write("\n")
        for i, dk in enumerate(dk_werte):
            file.write("%f" % (2*float(dk)))
            for t, nk in enumerate(Nk):
                file.write(" %f" % z_end_avg[t,i])
            file.write("\n")
        file.close()
    return 0


def write_zaverage_custom(mydir,Ns,sigma, Nk, dk_werte, z_end_avg, vflag, phrase):
    if (vflag == 0):
        myprefix="z"+phrase+"_vs_Nk_"
    elif (vflag == 1):
        myprefix="z"+phrase+"_vs_d_"

    filename=mydir+myprefix+str(Ns)+"_"+"{:.3f}".format(sigma)+".txt"
    file = open(filename, "w+")

    if (vflag == 0):
        file.write("# Nk")
        for dk in dk_werte:
            file.write(" %f" % dk*2)
        file.write("\n")
        #file.write("\n")
        #file.write("\n")
        for t, nk in enumerate(Nk):
            file.write("%d" % nk)
            for i, dk in enumerate(dk_werte):
                file.write(" %f" % z_end_avg[t,i])
            file.write("\n")
        file.close()

    elif (vflag == 1):
        file.write("# epsilon")
        for nk in Nk:
            file.write(" %d" % nk)
        file.write("\n")
        #file.write("\n")
        #file.write("\n")
        for i, dk in enumerate(dk_werte):
            file.write("%f" % (2*float(dk)))
            for t, nk in enumerate(Nk):
                file.write(" %f" % z_end_avg[t,i])
            file.write("\n")
        file.close()
    return 0
#def write_fdiff("./",Ns,sigma, Nk, dk_werte, fe_wall, fe_h_val, fe_h_pos, fe_diff, vflag):
def write_fdiff(mydir, Ns,sigma, Nk, dk_werte, fe_diff, vflag):
    myprefix="fdiff_"
    filename=mydir+myprefix+str(Ns)+"_"+"{:.3f}".format(sigma)+".txt"
    file = open(filename, "w+")

    file.write("# epsilon")
    for nk in Nk:
        file.write(" %d" % nk)
    file.write("\n")
    for i, dk in enumerate(dk_werte):
        file.write("%f" % (2*float(dk)))
        for t, nk in enumerate(Nk):
            file.write(" %f" % fe_diff[t,i])
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

if not os.path.exists("mb"):
        os.mkdir("mb")
#if not os.path.exists("nnmb"):
#        os.mkdir("nnmb")
if not os.path.exists("fe"):
        os.mkdir("fe")
if not os.path.exists("fe1"):
        os.mkdir("fe1")
#if not os.path.exists("MB"):
#        os.mkdir("MB")
#if not os.path.exists("NNMB"):
#        os.mkdir("NNMB")

#avail_list=["b_21-08-13_128_0.100_0.00.txt"]
avail_list=["b_21-08-13_128_0.100_0.00.txt"]

for avail in avail_list:

    Ns, sigma, chi, d, cen = read_oldformat(avail) # read in file of eq. brush - cen here is list, not np.array

    #print(Ns,sigma,chi,d)
    #plt.plot(range(len(cen)), cen)
    #plt.show()
    #exit()
    # jetzt mit Dichte cen[j] das Potential VarV rekonstruieren (ohne variation mit epsilon)
    # und in diesem Potential in einem Durchlauf mal die Greensfunktion G1 (evtl bis G2 und Gtot)
    # berechnen lassen für verschiedene Ns, d; sigma zu verändern widerspräche 
    # aber unter der Annahme dass es nur eine oder wenige Einzelketten sind...

    N_start = 128
    N_end = int(2.5*128) #512 #1024
    N_increment=8 #16
    #Nk=np.arange(N_start, N_end+N_increment, step=N_increment)
    #Nk=np.append(Nk,[128*,1280])
    #Nk=np.arange(8,128, step=8, dtype=int)
    #Nk=np.append(Nk, np.arange(N_start, N_end+N_increment, step=N_increment, dtype=int))
    #Nk=np.append(Nk, np.arange(int(N_start*3.0),int(N_start*6.5), step=int(N_start*0.5), dtype=int))
    #Nk=np.append(Nk, np.arange(int(N_start*7.0),int(N_start*11.0), step=int(N_start), dtype=int))
    #Nk=np.array([128, 136, 144, 152, 160, 168, 192, 256, 320, 336, 384, 448, 512, 640, 768, 1024, 1280],dtype=int)
    #["128","136","144","152","160","168","192","224","256","320","384","512","768","1024","1280"]
    Nk=np.array([167],dtype=int)
    #Nk=np.array([128, 136, 144, 152, 160, 168, 192, 224, 256, 288, 320, 352, 384, 448, 512, 640, 768, 1024, 1280],dtype=int)
    #Nk=np.array([8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,124,128,132,136,140,144,148,152,156,160,164,168,172,176,180,184,188,192,196,200,204,208,212,216,220,224,228,232,236,240,244,248,252,256,260,264,268,272,276,280,284,288,292,296,300,304,308,312,316,320,324,328,332,336,340,344,348,352,356,360,364,368,372,376,380,384],dtype=int)
    #Nk=np.array([128, 136, 144, 152, 160], dtype=int)
    N_end=np.max(Nk)
    N_anz=len(Nk)
    print(Nk)
    print("Nk(from,to,No): ",N_start,N_end,N_anz)

    
    d_start = 0.00 #  0.00
    d_end = 4.25 #100.0 #25.0 #50.0 #100.0 #100
    d_increment = 0.005
    dk_werte=np.arange(d_start, d_end+d_increment,step=d_increment)
    if max(dk_werte)>d_end:
        dk_werte=np.delete(dk_werte, -1)
    d_anz = len(dk_werte)
    #print(dk_werte)
    print("dk(from,to,No): ",d_start,d_end,d_anz)
    #print(d_end, max(dk_werte))
    print("Nm, d, sigma")

    chik=0
    w0=4.0/6.0                     # transition probability 
    w1=1.0/6.0                     # on simple cubic lattice
    norm1=0
    norm=0
    endmon_wand = np.zeros((N_anz,d_anz))
    loblock_wand = np.zeros((N_anz,d_anz))
    hiblock_wand = np.zeros((N_anz,d_anz))
    allmon_wand = np.zeros((N_anz,d_anz))

    z_end_avg = np.zeros((N_anz,d_anz))
    z_loblock_avg  = np.zeros((N_anz,d_anz))
    z_hiblock_avg  = np.zeros((N_anz,d_anz))
    z_all_avg = np.zeros((N_anz,d_anz))

    fe_wall  = np.zeros((N_anz,d_anz))
    fe_h_pos = np.zeros((N_anz,d_anz))
    fe_h_val = np.zeros((N_anz,d_anz))
    fe_diff = np.zeros((N_anz,d_anz))
    
    V0=0
    V_=[]
    for i in range(Ns-1):
        V_.append(np.exp(-1.*Vakt(cen[i], float(chik)/10.0) + V0))
    for i in range(N_end+2-Ns):
        V_.append(np.exp(-1.*0.0))                   # at this point V is a list, not an np.array... convert to np.array???
    V=np.array(V_)
    #for i in range(len(V)):
    #    V[i]=np.exp(-V[i])

    for i, dk in enumerate(dk_werte):
        wall=-float(2*dk)
        #wall=-2
        ############wall=-2
        print("MYENERGY: ", np.exp(wall))
        #energy=np.ones(max(Nk))
        Ewall = np.ones(max(Nk)+1)
        #for n in range(1,max(Nk)+1):
        #   Ewall[n] *= ((np.exp(-wall)-1.)/(1.*n*n) + 1.)
        #    Ewall[n] *= ( np.exp(-wall/(1.*n)) )
        #    Ewall[n] = ((np.exp(2.)-1.)/(1.*n) + 1.)
        Ewall[1] *= np.exp(-wall)
        Ewall_inv = 1./Ewall 

        #plt.plot(np.arange(1,(max(Nk)+2)), np.exp(2.)*1.*np.ones(max(Nk)+1), ls="--")
        ##plt.plot(np.arange(1,(max(Nk)+2)),1.*np.ones(max(Nk)+1), ls=":")
        #plt.plot(np.arange(1,(max(Nk)+2)), Ewall, c="r")
        #plt.plot(np.arange(1,(max(Nk)+2)), Ewall_inv, c="c")
        #plt.plot(np.arange(1,(max(Nk)+2)), Ewall*Ewall_inv, c="y", ls="-.")
        #print(Ewall, Ewall_inv, Ewall*Ewall_inv)    
        #plt.show()
        #exit()

        # nu - upper block, nl - lower block
        nl = 0; # nl=1 means the lowest adsorbed block is adsorbing, nl=2 means adsorbed + 1 extra
        nu = 25; 
        
        #if i>0:
        #    exit()
        for t, nk in enumerate(Nk):
            print(nk, dk, sigma)
            
            G1=np.zeros((nk,nk+2))
            G2=np.zeros((nk,nk+2))
            Gtot=np.zeros((nk,nk+1))
            phi=np.zeros(nk+1)
            # G[mon-index,z-position]
            G1[0,1]=1.0                           # calculate forward Gs ; an der Wand ist G1 fuer alle Gitterpkte = 1 (alle gegraftet)
            #for n in range(1,nk):                     # calculate forward Gs ; an der Wand ist G1 fuer alle Gitterpkte = 1 (alle gegraftet)
                #for j in range(1,n+2):                    # Monomere starten bei 2 = wand aber erstes mon wird als null gezählt in der laenge der kette
                    #G1[n,j] = w1*(G1[n-1,j-1]+G1[n-1,j+1])+w0*G1[n-1,j]      # S.111
                    #G1[n,j] *= V[j] #np.exp(-V[j])
            for n in range(1,nk):
                G1[n,1:n+2] = w1*(G1[n-1,:n+1]+G1[n-1,2:n+3])+w0*G1[n-1,1:n+2]
                G1[n,1:n+2] *= V[1:n+2]
                #print(n)
                if n < nl:
                    #print("G1,l",n)
                    G1[n,1:nk+1] *= Ewall[1:nk+1]
                if n > nk-1-nu:
                    #print(G1[n,1:nk+1])
                    #print("G1,u",n)
                    G1[n,1:nk+1] *= Ewall[1:nk+1]
                    #print(G1[n,1:nk+1])
            #exit()

                
            #for n in range(1,nk+1): # 1/r potential
            #    G1[nk-1,n] *= Ewall[n] #((np.exp(-wall)-1.)/(1.*n) + 1.)                   # additional potential at surface, belohnung energetisch //letztes Monomer ist an Wand
            #print(G1[nk-1,1:nk+1])
            #print(nk-1)
            #G1[nk-1,1:nk+1] *= Ewall[1:nk+1]
            #print(G1[nk-1,1:nk+1])
            #exit()

            norm1 = np.sum(G1[nk-1,:])                 # norm1 = Summe aller Wahrsch. fuer das letzte Monomer  (aufsummieren fuer alle Orte z)
            norm = sigma/norm1                        # anschlussbedingung fuer die G2, zusammen ergibt dann G1*G2 bei Ns-1 wieder sigma

            #for k in range(1,nk+1):                     # initial values for backward Gs Green function rueckw erzeugen
            #    G2[nk-1,k] = norm*V[k] #norm*np.exp(-V[k])       # Variation des Pot., Auswirkung auf Pot bei Ort k wenn ich an Ort h Pot variiere
            G2[nk-1,1:nk+1] = norm*V[1:nk+1]


            #for n in range(1,nk+1):# 1/r potential
            #    G2[nk-1,n] *= Ewall[n] #((np.exp(-wall)-1.)/(1.*n) + 1.)               # falls letztes Monomer an Wand
            #print("G2,u", nk-1)
            G2[nk-1,1:nk+1] *= Ewall[1:nk+1]


            for n in range(nk-1):                    # calculate backward Gs
            #    for j in range(1,nk-n):
            #        G2[nk-2-n,j] = w1*(G2[nk-1-n,j-1]+G2[nk-1-n,j+1])+w0*G2[nk-1-n,j] 
            #        G2[nk-2-n,j] *= V[j] #np.exp(-V[j])
                #print(G2[nk-2-n,1:nk-n])
                G2[nk-2-n,1:nk-n] = w1*(G2[nk-1-n,:nk-n-1]+G2[nk-1-n,2:nk-n+1])+w0*G2[nk-1-n,1:nk-n]
                G2[nk-2-n,1:nk-n] *= V[1:nk-n]
                #print(G2[nk-2-n,1:nk-n])
                if nk-2-n < nl and nk-2-n > 0:
                    #print("G2,l",nk-2-n)
                    G2[nk-2-n,1:nk-n] *= Ewall[1:nk-n]
                if nk-2-n > nk-1-nu:
                    #print("G2,u",nk-2-n)
                    G2[nk-2-n,1:nk-n] *= Ewall[1:nk-n]
            #exit()
                #exit()


            #for n in range(nk):                      # calculate total Gs
            #    for k in range(1,nk+1):                     
            #        Gtot[n,k] = G1[n,k]*G2[n,k]*(1./V[k])  #*np.exp(V[k])      # e⁻V ist das #1/G(z), ist bereits das phi(s,z)
            Gtot = np.multiply(G1,G2)
            #for n in range(nk):    
                #for k in range(1,nk+1):
                #    Gtot[n,k] /=  V[k]
            Gtot[:nk,1:nk+1] /=  V[1:nk+1]

            #for n in range(1,nk+1):# 1/r potential
            #    Gtot[nk-1,n] *= Ewall_inv[n] #((np.exp(wall)-1.)/(1.*n) + 1.)                # faktor exp(-wall) kommt durch die beiden Gs 1x zu oft vor daher, wieder rausmultipl.
            #Gtot[nk-1,1:nk+1] *= Ewall_inv[1:nk+1]
            #print("EWALL_INV")
            #print(Ewall_inv[1:nk+1])
            #print(Gtot[nk-1, 1:nk+1])
            #Gtot[nk-nu:nk, 1:nk+1] *= Ewall_inv[1:nk+1]
            #print(*range(nk-nu,nk))
            #print(Gtot[nk-1:nk, 1:nk+1])

            #print(Gtot[:nl,1])
            #print(Gtot[nk-nu:nk,1])
            for n in range(1,nl):
                #print("Gtot",n)
                Gtot[n, 1:nk+1] *=Ewall_inv[1:nk+1]
            for n in range(nk-nu, nk):
                #print("Gtot",n)
                Gtot[n, 1:nk+1] *=Ewall_inv[1:nk+1]

            #print(Gtot[:nl,1])
            #print(Gtot[nk-nu:nk,1])
            #print(np.sum(Gtot[nk-nu:nk,1])/(sigma*nu))
            #print(Gtot[0,:])
            #exit()

            phi=np.sum(Gtot, axis=0)
            
            endmon_wand[t,i] = Gtot[nk-1,1]/sigma
            loblock_wand[t,i] = np.sum(Gtot[:nl,1])/(nl*sigma)
            #print(np.sum(Gtot[:nl,1])/(nl*sigma))
            hiblock_wand[t,i] = np.sum(Gtot[nk-nu:nk,1])/(nu*sigma)
            allmon_wand[t,i] = np.sum(Gtot[:nk,1])/(nk*sigma)

            #print(allmon_wand[t,i], loblock_wand[t,i], hiblock_wand[t,i], endmon_wand[t,i])

            #print(np.sum(Gtot[:nl,1])/(nl*sigma))
            #print(np.sum(Gtot[:nk,1])/(nk*sigma))
            #exit()

            #print(phi)
            #print(Gtot[nk-1,:]/sigma)
            #print(np.sum(phi))
            #print(np.sum(Gtot[nk-1,:]))
            #print(sigma)
            #exit()
            
            #print(np.shape(Gtot[nk-1,:]))
            #print(np.shape(np.linspace(0,nk,nk+1)))
            #print(np.linspace(0,nk+1,nk+2))
            #print(Gtot[nk-1,:]/sigma)

            #print(Gtot[nk-1,:])
            #print(np.linspace(0,nk+1,nk+2))
            #print(np.dot(np.linspace(0,nk+1,nk+2),Gtot[nk-1,:])/np.sum(Gtot[nk-1,:]))
            z_end_avg[t,i] = np.dot(np.linspace(0,nk+1,nk+2),Gtot[nk-1,:])/np.sum(Gtot[nk-1,:])
            #print("z_end_avg = ",z_end_avg[t,i])


            z_loblock_avg[t,i] = 0
            for n in range(0,nl):
                z_loblock_avg[t,i] += np.dot(np.linspace(0,nk+1,nk+2),Gtot[n,:])/np.sum(Gtot[n,:])
            z_loblock_avg[t,i] /= nl
            #print("z_loblock_avg = ",z_loblock_avg[t,i])
            

            z_hiblock_avg[t,i] = 0
            for n in range(nk-nu,nk):
                z_hiblock_avg[t,i] += np.dot(np.linspace(0,nk+1,nk+2),Gtot[n,:])/np.sum(Gtot[n,:])
            z_hiblock_avg[t,i] /= nu
            #print("z_hiblock_avg = ",z_hiblock_avg[t,i])

            z_all_avg[t,i] = 0
            for n in range(nk):
                z_all_avg[t,i] += np.dot(np.linspace(0,nk+1,nk+2),Gtot[n,:])/np.sum(Gtot[n,:])
            z_all_avg[t,i] /= nk
            #print("z_all_avg = ",z_all_avg[t,i])

            #exit()
            #z_loblock_avg[t,i] = 
            #z_hiblock_avg[t,i] =
            #z_all_avg[t,i] = 

            #fe_loc = np.array([-1.0*np.log(Gtot[nk-1,k]/np.sum(Gtot[nk-1,:])) if Gtot[nk-1,k] != 0 else 0.0 for k in range(nk)])
            fe_loc = np.array([-1.0*np.log(item/np.sum(Gtot[nk-1,:])) if item != 0 else 0.0 for index, item in enumerate(Gtot[nk-1,:])])
            fe_loc_max = argrelextrema(fe_loc, np.greater) 
            fe_loc_min = argrelextrema(fe_loc, np.less)
            fe_h_pos[t,i] = int(fe_loc_min[0][0]) if (len(fe_loc_max[0])>0 and len(fe_loc_min[0])>0) else int(fe_loc_max[0][0])
            #fe_wall[t,i] = -1.0*np.log(Gtot[nk-1,1]/np.sum(Gtot[nk-1,:]))
            #fe_h_pos[t,i] = int(argrelextrema(fe_loc, np.less)[0][0])
            fe_h_val[t,i] = fe_loc[int(fe_h_pos[t,i])]
            #fe_diff[t,i] = fe_wall[t,i] - fe_h_val[t,i]
            fe_diff[t,i] = fe_loc[1]-fe_h_val[t,i]
            
            # write out files with density profiles etc. 
            write_oldformat("./mb/", "mb", Ns, nk, sigma, wall, phi, Gtot, V0, norm, 0)
            #write_oldformat("./nnmb/", "nnmb", Ns, nk, sigma, wall, phi, Gtot, V0, norm, 1)
            #write_newformat("./MB/", "MB",Ns, nk, sigma, wall, phi, Gtot, 0)
            #write_newformat("./NNMB/", "NNMB",Ns, nk, sigma, wall, phi, Gtot, 1)
            write_fe("./fe/", "fe", Ns, nk, sigma, wall, Gtot)
            write_fe_norm1("./fe1/", "fe1", Ns, nk, sigma, wall, Gtot)

        

    print("Adsorptionskurven schreiben")
    write_adsorption("./",Ns,sigma, Nk, dk_werte, endmon_wand, 0, "endm") # vs Nk
    write_adsorption("./",Ns,sigma, Nk, dk_werte, endmon_wand, 1, "endm") # vs d

    write_adsorption("./",Ns,sigma, Nk, dk_werte, loblock_wand, 1, "loblock") # vs d
    write_adsorption("./",Ns,sigma, Nk, dk_werte, hiblock_wand, 1, "hiblock") # vs d
    write_adsorption("./",Ns,sigma, Nk, dk_werte, allmon_wand, 1, "allm") # vs d

    #write_zaverage("./",Ns,sigma, Nk, dk_werte, z_end_avg, 0)
    write_zaverage("./",Ns,sigma, Nk, dk_werte, z_end_avg, 1)
    write_zaverage_custom("./",Ns,sigma, Nk, dk_werte, z_end_avg, 1, "end")
    write_zaverage_custom("./",Ns,sigma, Nk, dk_werte, z_loblock_avg, 1, "loblock")
    write_zaverage_custom("./",Ns,sigma, Nk, dk_werte, z_hiblock_avg, 1, "hiblock")
    write_zaverage_custom("./",Ns,sigma, Nk, dk_werte, z_all_avg, 1, "all")

    #write_fdiff("./",Ns,sigma, Nk, dk_werte, fe_wall, fe_h_val, fe_h_pos, fe_diff, 1)
    write_fdiff("./",Ns,sigma, Nk, dk_werte,  fe_diff, 1)
    print("Program end")
