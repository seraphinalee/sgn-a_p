from matplotlib import pyplot as plt
import csv



RealNumber = RealField(100)
print(parent(1.0))


## bound
bd = 10000000
s = 0.5

d = 1


noncm = RealNumber(4/pi)
hascm = RealNumber(8/(3*pi))

lst = ## put list of labels of elliptic curves here

count0 = 0
counts = ''
exp_val = 1


for lbl in lst:
    E=EllipticCurve(lbl).minimal_model()
    
    if E.has_cm():
        exp_val = hascm
    else:
        exp_val = noncm
    
    B_new1=1.00;
    B_new2=1.00;
    B_new3=1.00;
    B_new4=1.00;
    B_new5=1.00;
    
    B_old=1.00;

    ## euler products
    BSDlist=[];
    MBSDlist1=[];
    MBSDlist2=[];
    MBSDlist3=[];
    MBSDlist4=[];
    MBSDlist5=[];

    ## x-coordinates
    xpoints=[];
    
    
    p=2


    for j in range(bd):
        ap = E.ap(p)
        
        ## one more factor in the euler product, for various values of c
        B_new1 = B_new1*((p + 1 - 0.5 * sgn(ap) * (p**(0.5)))/p)
        B_new2 = B_new2*((p + 1 - 0.7 * sgn(ap) * (p**(0.5)))/p)
        B_new3 = B_new3*((p + 1 - sgn(ap) * (p**(0.5)))/p)
        B_new4 = B_new4*((p + 1 - exp_val * sgn(ap) * (p**(0.5)))/p)
        B_new5 = B_new5*((p + 1 - 1.5 * sgn(ap) * (p**(0.5)))/p)
        
        B_old = B_old*(1-(ap * 1.0/ p**(2*s)) + p**(1-2*(2*s)) ) ## = Np/p
        
        
        ## only plot once every 10000
        if j % 10000 == 0 and j > 0:
            MBSDlist1.append(log(B_new1)/log(log(1.0*p)))
            MBSDlist2.append(log(B_new2)/log(log(1.0*p)))
            MBSDlist3.append(log(B_new3)/log(log(1.0*p)))
            MBSDlist4.append(log(B_new4)/log(log(1.0*p)))
            MBSDlist5.append(log(B_new5)/log(log(1.0*p)))
            BSDlist.append(log(B_old)/log(log(1.0*p)))
            xpoints.append(p)
        if j %10000 == 0:
            print(j)
        p=next_prime(p)
    
    ## example in how to save file names
    rank = '?'
    if lbl in lst:
        rank = '0'
        counts = str(count0)
        count0 += 1


    print("rank="+rank)
    print('===')
    
    
    plt.plot(xpoints, MBSDlist1, 'bo', markersize = 1)
    plt.plot(xpoints, MBSDlist2, 'go', markersize = 1)
    plt.plot(xpoints, MBSDlist3, 'mo', markersize = 1)
    plt.plot(xpoints, MBSDlist4, 'ro', markersize = 1)
    plt.plot(xpoints, MBSDlist5, 'ko', markersize = 1)
    plt.plot(xpoints, BSDlist, 'co', markersize = 1)
    
    
    plt.xlabel('p')
    plt.ylabel('log(prod N_{p,c}/p)/log log p')
    plt.title(str(lbl)+', Rank = '+rank+', CM='+str(E.has_cm()))
    plt.savefig('/Users/user/Desktop/data/'+rank+'_'+counts+'.pdf')
    plt.clf()
    plt.cla()
    plt.close()
    
    MBSDlist1.insert(0, 'c = 0.5')
    MBSDlist2.insert(0, 'c = 0.7')
    MBSDlist3.insert(0, 'c = 1')
    MBSDlist4.insert(0, 'c = exp_val')
    MBSDlist5.insert(0, 'c = 1.5')
    BSDlist.insert(0, 'BSD')
    with open('/Users/user/Desktop/data/'+rank+'_'+counts+'.csv', 'w') as f1:
        writefile = csv.writer(f1)
        writefile.writerows([lbl, MBSDlist1, MBSDlist2, MBSDlist3, MBSDlist4, MBSDlist5, BSDlist])
    
    