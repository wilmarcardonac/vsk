import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as py 

l,vlsmica = np.loadtxt('./vsk_angular_power_spectrum/vl_smica.txt',unpack=True,usecols=[0,1])

slsmica = np.loadtxt('./vsk_angular_power_spectrum/sl_smica.txt',unpack=True,usecols=[1])

klsmica = np.loadtxt('./vsk_angular_power_spectrum/kl_smica.txt',unpack=True,usecols=[1])

#vlnilc = np.loadtxt('../inpainted-nilc/vsk/vsk_angular_power_spectrum/vl_nilc.txt',unpack=True,usecols=[1])

#klnilc = np.loadtxt('../inpainted-nilc/vsk/vsk_angular_power_spectrum/kl_nilc.txt',unpack=True,usecols=[1])

#slnilc = np.loadtxt('../inpainted-nilc/vsk/vsk_angular_power_spectrum/sl_nilc.txt',unpack=True,usecols=[1])

multipoles = len(l)

number_maps = 1000

vlmean = np.zeros(multipoles)

slmean = np.zeros(multipoles)

klmean = np.zeros(multipoles)

vl = np.zeros((multipoles,number_maps))

sl = np.zeros((multipoles,number_maps))

kl = np.zeros((multipoles,number_maps))

for index in range(1,number_maps+1):

    vl[:,index-1] = np.loadtxt('./vsk_angular_power_spectrum/vl_'+str(index).zfill(4)+'.txt',unpack=True,usecols=[1])

    sl[:,index-1] = np.loadtxt('./vsk_angular_power_spectrum/sl_'+str(index).zfill(4)+'.txt',unpack=True,usecols=[1])

    kl[:,index-1] = np.loadtxt('./vsk_angular_power_spectrum/kl_'+str(index).zfill(4)+'.txt',unpack=True,usecols=[1])


#68%

vl68 = np.zeros((2,multipoles))

sl68 = np.zeros((2,multipoles))

kl68 = np.zeros((2,multipoles))

for index in range(len(l)):

    vl68[0,index] = np.min(vl[index,:])

    vl68[1,index] = np.max(vl[index,:])

    step = abs((vl68[1,index]-vl68[0,index])/1000000.)

    vl68[0,index] = vl68[0,index] - step

    vl68[1,index] = vl68[1,index] + step

    counter = np.compress((vl68[0,index]<vl[index,:]) & ( vl68[1,index]>vl[index,:] ),vl[index,:]).size

    while (counter > 0.68*number_maps):

        vl68[0,index] = vl68[0,index] + step

        vl68[1,index] = vl68[1,index] - step

        counter = np.compress(( vl68[0,index]<vl[index,:]) & ( vl68[1,index] >vl[index,:]),vl[index,:]).size

    sl68[0,index] = np.min(sl[index,:])

    sl68[1,index] = np.max(sl[index,:])

    step = abs((sl68[1,index]-sl68[0,index])/1000000.)

    sl68[0,index] = sl68[0,index] - step

    sl68[1,index] = sl68[1,index] + step

    counter = np.compress((sl68[0,index]<sl[index,:]) & ( sl68[1,index]>sl[index,:] ), sl[index,:]).size

    while (counter > 0.68*number_maps):

        sl68[0,index] = sl68[0,index] + step

        sl68[1,index] = sl68[1,index] - step

        counter = np.compress(( sl68[0,index]<sl[index,:]) & ( sl68[1,index]>sl[index,:] ), sl[index,:]).size


    kl68[0,index] = np.min(kl[index,:])

    kl68[1,index] = np.max(kl[index,:])

    step = abs((kl68[1,index]-kl68[0,index])/1000000.)

    kl68[0,index] = kl68[0,index] - step

    kl68[1,index] = kl68[1,index] + step

    counter = np.compress((kl68[0,index]<kl[index,:]) & ( kl68[1,index]>kl[index,:] ), kl[index,:]).size

    while (counter > 0.68*number_maps):

        kl68[0,index] = kl68[0,index] + step

        kl68[1,index] = kl68[1,index] - step

        counter = np.compress(( kl68[0,index]<kl[index,:]) & ( kl68[1,index]>kl[index,:] ), kl[index,:]).size


#95%

vl95 = np.zeros((2,multipoles))

sl95 = np.zeros((2,multipoles))

kl95 = np.zeros((2,multipoles))

for index in range(len(l)):

    vl95[0,index] = np.min(vl[index,:])

    vl95[1,index] = np.max(vl[index,:])

    step = abs((vl95[1,index]-vl95[0,index])/1000000.)

    vl95[0,index] = vl95[0,index] - step

    vl95[1,index] = vl95[1,index] + step

    counter = np.compress((vl95[0,index]<vl[index,:]) & ( vl95[1,index]>vl[index,:] ),vl[index,:]).size

    while (counter > 0.95*number_maps):

        vl95[0,index] = vl95[0,index] + step

        vl95[1,index] = vl95[1,index] - step

        counter = np.compress(( vl95[0,index]<vl[index,:]) & ( vl95[1,index] >vl[index,:]),vl[index,:]).size

    sl95[0,index] = np.min(sl[index,:])

    sl95[1,index] = np.max(sl[index,:])

    step = abs((sl95[1,index]-sl95[0,index])/1000000.)

    sl95[0,index] = sl95[0,index] - step

    sl95[1,index] = sl95[1,index] + step

    counter = np.compress((sl95[0,index]<sl[index,:]) & ( sl95[1,index]>sl[index,:] ), sl[index,:]).size

    while (counter > 0.95*number_maps):

        sl95[0,index] = sl95[0,index] + step

        sl95[1,index] = sl95[1,index] - step

        counter = np.compress(( sl95[0,index]<sl[index,:]) & ( sl95[1,index]>sl[index,:] ), sl[index,:]).size


    kl95[0,index] = np.min(kl[index,:])

    kl95[1,index] = np.max(kl[index,:])

    step = abs((kl95[1,index]-kl95[0,index])/1000000.)

    kl95[0,index] = kl95[0,index] - step

    kl95[1,index] = kl95[1,index] + step

    counter = np.compress((kl95[0,index]<kl[index,:]) & ( kl95[1,index]>kl[index,:] ), kl[index,:]).size

    while (counter > 0.95*number_maps):

        kl95[0,index] = kl95[0,index] + step

        kl95[1,index] = kl95[1,index] - step

        counter = np.compress(( kl95[0,index]<kl[index,:]) & ( kl95[1,index]>kl[index,:] ), kl[index,:]).size

    vlmean[index] = np.mean(vl[index,:])

    slmean[index] = np.mean(sl[index,:])

    klmean[index] = np.mean(kl[index,:])

print "MEANS :"

print "V ", vlmean

print "S ", slmean

print "K ", klmean

py.figure(1)

py.plot(l[1:],vlsmica[1:],label=r'INPAINTED SMICA',linestyle=':')

#py.plot(l[1:],vlnilc[1:],label=r'INPAINTED NILC',linestyle='dashed')

py.plot(l[1:],vlmean[1:],label=r'simul. mean')

py.fill_between(l[1:],vl68[0,1:],vl68[1,1:],facecolor='gray',alpha=0.4)

py.fill_between(l[1:],vl95[0,1:],vl95[1,1:],facecolor='gray',alpha=0.3)

py.xlabel(r'$\ell$',fontsize='large')

py.xlim(0,multipoles)

py.xticks(np.arange(multipoles))

# Label for y-axis

py.ylabel(r'$V_{\ell}\quad[\mu K_{cmb}^4]$',fontsize='large')

#py.ylim(1.e-25,1.e-17)

py.yscale('log')

py.legend(loc=0,ncol=2)#'lower right')

#py.show()

py.savefig("./figures/Vl_ut78.pdf")

py.close(1)

py.figure(2)

py.plot(l[1:],slsmica[1:],label=r'INPAINTED SMICA',linestyle=':')

#py.plot(l[1:],slnilc[1:],label=r'INPAINTED NILC',linestyle='dashed')

py.plot(l[1:],slmean[1:],label=r'simul. mean')

py.fill_between(l[1:],sl68[0,1:],sl68[1,1:],facecolor='gray',alpha=0.4)

py.fill_between(l[1:],sl95[0,1:],sl95[1,1:],facecolor='gray',alpha=0.3)

py.xlabel(r'$\ell$',fontsize='large')

py.xlim(0,multipoles)

py.xticks(np.arange(multipoles))

# Label for y-axis

py.ylabel(r'$S_{\ell}$',fontsize='large')

py.ylim(4.e-5,5.e-2)

py.yscale('log')

py.legend(loc=0,ncol=2)#'lower right')

#py.show()

py.savefig("./figures/Sl_ut78.pdf")

py.close(2)

py.figure(3)

py.plot(l[1:],klsmica[1:],label=r'INPAINTED SMICA',linestyle=':')

#py.plot(l[1:],klnilc[1:],label=r'INPAINTED NILC',linestyle='dashed')

py.plot(l[1:],klmean[1:],label=r'simul. mean')

py.fill_between(l[1:],kl68[0,1:],kl68[1,1:],facecolor='gray',alpha=0.4)

py.fill_between(l[1:],kl95[0,1:],kl95[1,1:],facecolor='gray',alpha=0.3)

py.xlabel(r'$\ell$',fontsize='large')

py.xlim(0,multipoles)

py.xticks(np.arange(multipoles))

# Label for y-axis

py.ylabel(r'$K_{\ell}$',fontsize='large')

py.ylim(1.e-5,1.e-1)

py.yscale('log')

py.legend(loc=0,ncol=2)#'lower right')

#py.show()

py.savefig("./figures/Kl_ut78.pdf")

py.close(3)

for index in range(number_maps+1,2*number_maps+1):

    vl[:,index-1-number_maps] = np.loadtxt('./vsk_angular_power_spectrum/vl_'+str(index).zfill(4)+'.txt',unpack=True,usecols=[1])

    sl[:,index-1-number_maps] = np.loadtxt('./vsk_angular_power_spectrum/sl_'+str(index).zfill(4)+'.txt',unpack=True,usecols=[1])

    kl[:,index-1-number_maps] = np.loadtxt('./vsk_angular_power_spectrum/kl_'+str(index).zfill(4)+'.txt',unpack=True,usecols=[1])

for index in range(len(l)):

    vlmean[index] = np.mean(vl[index,:])

    slmean[index] = np.mean(sl[index,:])

    klmean[index] = np.mean(kl[index,:])

print "MEANS :"

print "V ", vlmean

print "S ", slmean

print "K ", klmean

# Arrays to store angular power spectra minus mean values

vvar = np.zeros((multipoles,number_maps))
svar = np.zeros((multipoles,number_maps))
kvar = np.zeros((multipoles,number_maps))

# We fill the arrays 

for index1 in range(multipoles):
    for index2 in range(number_maps):
        vl[index1,index2] = vl[index1,index2] - vlmean[index1] 
        sl[index1,index2] = sl[index1,index2] - slmean[index1] 
        kl[index1,index2] = kl[index1,index2] - klmean[index1] 

# Computing covariance matrices and their inverses 

vcov = np.zeros((multipoles,multipoles))

scov = np.zeros((multipoles,multipoles))

kcov = np.zeros((multipoles,multipoles))

vcov = np.cov(vl,rowvar=1)

scov = np.cov(sl,rowvar=1)

kcov = np.cov(kl,rowvar=1)

vinvcov = np.zeros((multipoles,multipoles))

sinvcov = np.zeros((multipoles,multipoles))

kinvcov = np.zeros((multipoles,multipoles))

vinvcov = np.linalg.inv(vcov)

sinvcov = np.linalg.inv(scov)

kinvcov = np.linalg.inv(kcov)

# We define and compute chi^2 (S-K-V estimators) for a set of Gaussian maps 

vchi2 = np.zeros(number_maps)

schi2 = np.zeros(number_maps)

kchi2 = np.zeros(number_maps)

for index3 in range(number_maps):
    for index1 in range(multipoles):
        for index2 in range(multipoles):
            vchi2[index3] += vl[index1,index3]*vinvcov[index1,index2]*vl[index2,index3]
            schi2[index3] += sl[index1,index3]*sinvcov[index1,index2]*sl[index2,index3]
            kchi2[index3] += kl[index1,index3]*kinvcov[index1,index2]*kl[index2,index3]

smicavchi2 = 0 

smicaschi2 = 0

smicakchi2 = 0

nilcvchi2 = 0 

nilcschi2 = 0

nilckchi2 = 0

for index1 in range(multipoles):
    for index2 in range(multipoles):
        smicavchi2 += (vlsmica[index1] - vlmean[index1])*vinvcov[index1,index2]*(vlsmica[index2] - vlmean[index2])
        smicaschi2 += (slsmica[index1] - slmean[index1])*sinvcov[index1,index2]*(slsmica[index2] - slmean[index2])
        smicakchi2 += (klsmica[index1] - klmean[index1])*kinvcov[index1,index2]*(klsmica[index2] - klmean[index2])        
#        nilcvchi2 += (vlnilc[index1] - vlmean[index1])*vinvcov[index1,index2]*(vlnilc[index2] - vlmean[index2])
#        nilcschi2 += (slnilc[index1] - slmean[index1])*sinvcov[index1,index2]*(slnilc[index2] - slmean[index2])
#        nilckchi2 += (klnilc[index1] - klmean[index1])*kinvcov[index1,index2]*(klnilc[index2] - klmean[index2])        


# lower tail probability 

def lower_tail_probability(distchi2,chi2):
    counter = len(distchi2[distchi2 < chi2])
    low_tail_pro = float(float(counter)/float(len(distchi2)))
    return low_tail_pro
        
# We compute histograms for the chi^2 vectors 

py.figure(4)

vhisto = py.hist(vchi2,fill=False)

py.ylabel('Frequency')

py.xlabel(r'$\chi^2_{V}$',fontsize='large')

ymax = 1000

py.ylim(0,ymax)

py.vlines(smicavchi2,0,ymax,linestyle=':',label=r'INPAINTED SMICA')

#py.vlines(nilcvchi2,0,ymax,linestyle='dashed',label=r'INPAINTED NILC')

py.legend(loc=0)

py.savefig("./figures/vchi2.pdf")

py.close(4)

#py.show()

py.figure(5)

shisto = py.hist(schi2,fill=False)

py.ylabel('Frequency')

py.xlabel(r'$\chi^2_{S}$',fontsize='large')

py.ylim(0,ymax)

py.vlines(smicaschi2,0,ymax,linestyle=':',label=r'INPAINTED SMICA')

#py.vlines(nilcschi2,0,ymax,linestyle='dashed',label=r'INPAINTED NILC')

py.legend(loc=0)

py.savefig("./figures/schi2.pdf")

py.close(5)

#py.show()

py.figure(6)

khisto = py.hist(kchi2,fill=False)

py.ylabel('Frequency')

py.xlabel(r'$\chi^2_{K}$',fontsize='large')

py.ylim(0,ymax)

py.vlines(smicakchi2,0,ymax,linestyle=':',label=r'INPAINTED SMICA')

#py.vlines(nilckchi2,0,ymax,linestyle='dashed',label=r'INPAINTED NILC')

py.legend()

py.savefig("./figures/kchi2.pdf")

py.close(6)

#py.show()

print " INPAINTED SMICA:lower tail probability variance, skewness, kurtosis "
print lower_tail_probability(vchi2,smicavchi2)
print lower_tail_probability(schi2,smicaschi2)
print lower_tail_probability(kchi2,smicakchi2)

#print " INPAINTED NILC:lower tail probability variance, skewness, kurtosis "
#print lower_tail_probability(vchi2,nilcvchi2)
#print lower_tail_probability(schi2,nilcschi2)
#print lower_tail_probability(kchi2,nilckchi2)

exit()
