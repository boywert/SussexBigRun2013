import numpy
import os
import math
import struct
import copy
import convert_config
import time as libtime

global SNAPfile
global AHFdir
global AHFprefix

global G
global m2Mpc
global m2km
global kpc2Mpc
global Msun2Gadget
global kg2Msun
global NFILES
global spin_model


NFILES = 8
G = 6.67384e-11 # m^3/(kgs^2)
m2Mpc = 1./3.08567758e22
m2km = 0.001
kpc2Mpc = 0.001
Msun2Gadget = 1.e-10
Msun2kg = 1.989e30
kg2Msun = 1./Msun2kg
#change to (Mpc/h) (km/s)^2 / (1e10Msun/h)
G = G*m2Mpc*m2km**2./(Msun2Gadget*kg2Msun)
print G

AHFdir = convert_config.AHFdir
AHFprefix = convert_config.AHFprefix
SUSSINGtree = convert_config.SUSSINGtree
SNAPfile = convert_config.SNAPfile
FileOut = convert_config.FileOut
FileOut2 = convert_config.FileOut2
spin_model = convert_config.spin_model

def readAHFascii():
    halocat = {}
    timesnap = numpy.loadtxt(SNAPfile)
    for time in timesnap:
        zstring = "%.3f" % (time[2])
        #print zstring[len(zstring)-1]
        filename = "%s/%s_%03d.z%s.AHF_halos" % (AHFdir, AHFprefix, time[0], zstring)
        #print "checking "+filename
        if os.path.isfile(filename) == False:
            zstring = "%.3f" % (time[2]+0.00001)
            #print zstring[len(zstring)-1]
            filename = "%s/%s_%03d.z%s.AHF_halos" % (AHFdir, AHFprefix, time[0], zstring)
        if os.path.isfile(filename) == False:
            zstring = "%.3f" % (time[2]-0.00001)
            #print "checking "+filename
            #print zstring[len(zstring)-1]
            filename = "%s/%s_%03d.z%s.AHF_halos" % (AHFdir, AHFprefix, time[0], zstring)
        if os.path.isfile(filename) == False:
            print "checking error"+filename
            exit()
        print "Reading "+filename
        stat = os.stat(filename)
        #print stat.st_size
        if(stat.st_size > 384):
            data = numpy.loadtxt(filename)
            shape = data.shape
            if (len(shape) == 1) & (shape[0] > 1):
                data_tmp = []
                data_tmp.append(data)
                data = data_tmp

            for halo in data:
                hid = long(halo[0])
                #print hid
                halocat[hid] = {}
                halocat[hid]["Redshift"] = time[2];
                halocat[hid]["UID"] = long(halo[0])
                halocat[hid]["ID"] = hid
                halocat[hid]["Mvir"] = halo[3]*Msun2Gadget
                halocat[hid]["Rvir"] = halo[11]*kpc2Mpc
                if(halo[44] < 1.0e34):
                    halocat[hid]["M_200Mean"] = halo[44]*Msun2Gadget
                else:
                    halocat[hid]["M_200Mean"] = 0.
                if(halo[45] < 1.0e34):
                    halocat[hid]["M_200Crit"] = halo[45]*Msun2Gadget
                else:
                    halocat[hid]["M_200Crit"] = 0.
                if(halo[46] < 1.0e34):
                    halocat[hid]["M_TopHat"] =  halo[46]*Msun2Gadget
                else:
                    halocat[hid]["M_200Crit"] = 0.
                    
                halocat[hid]["Len"] = halo[4]
                halocat[hid]["Pos"] = (halo[5]*kpc2Mpc,halo[6]*kpc2Mpc,halo[7]*kpc2Mpc)
                halocat[hid]["Vel"] = (halo[8],halo[9],halo[10])
                halocat[hid]["Vmax"] = halo[16]
                halocat[hid]["VelDisp"] = halo[18]
                lambda_bullock = halo[19]
                # 
                # 
                
                if(spin_model == 99): # Boyd's stupid model
                    lambda_bullock = 0.02

                
                # use Peebles lambdaE definition to find angular momentum
                #halocat[hid]["Ep"] = halo[38]
                
                #halocat[hid]["Ek"] = halo[39]
                # total_energy = math.fabs((halo[38] + halo[39])*Msun2Gadget)
                #halocat[hid]["LambdaE"] = halo[20]
                # J = halo[20]*G*halocat[hid]["Mvir"]**(3./2.)/total_energy**(0.5)
                J = lambda_bullock*numpy.sqrt(2.0*G*halocat[hid]["Mvir"]*halocat[hid]["Rvir"])
                
                #halocat[hid]["TotalEnergy"] = total_energy
                halocat[hid]["Spin"] = (halo[21]*J,halo[22]*J,halo[23]*J)
                halocat[hid]["FirstProgenitor"] = -1
                halocat[hid]["NextProgenitor"] = -1
                halocat[hid]["Descendant"] = -1
                halocat[hid]["HostHalo"] = long(halo[1])
                if(halocat[hid]["HostHalo"] == 0):
                    halocat[hid]["HostHalo"] = -1
                halocat[hid]["NextHalo"] = -1
                halocat[hid]["SnapNum"] = long(time[0])
                #halocat[hid]["NextinTree"] = -1
                halocat[hid]["HaloNr"] = -1
                halocat[hid]["TreeNr"] = -1
                halocat[hid]["movetonew"] = -1
                # print halocat[hid]['Spin'] 
    print "Make host-sub structures ..."
    for haloc in halocat.iterkeys():
        #print haloc
        halo = halocat[haloc]
        hosthalo= halo["HostHalo"]
        upperhost = halo["HostHalo"]
        ref = upperhost
        while upperhost > -1:
            ref = upperhost
            upperhost = halocat[upperhost]["HostHalo"]
        if(hosthalo != ref):
            hosthalo = ref
        halocat[haloc]["MainHalo"] = hosthalo
        if(halocat[haloc]["MainHalo"] > -1):
            cursub = halocat[haloc]["MainHalo"]
            while cursub > -1:
                curid = halocat[cursub]["ID"]
                cursub = halocat[cursub]["NextHalo"]
            #print curid, "change",halocat[curid]["NextHalo"],"to",haloc
            halocat[curid]["NextHalo"] = haloc

    return halocat

def readSussingtree(SUSSINGtree,halocat):
    print "Copy halo catalogue ..."
    halocopy = copy.copy(halocat)
    print "Read tree file ..."
    f = open(SUSSINGtree)
    line = f.read().splitlines()
    count = 0;
    for (i,item) in enumerate(line):
        if(i == 2):
            ngals_read = item.split()
            totalhalo = long(ngals_read[0])
            print "tree",totalhalo,"halocat",len(halocat)
        if(i >= 3):
            col = item.split()
            if(count == 0):
                if(col[0] == "END"):
                    print "finish reading ",SUSSINGtree
                    return halocopy
                if(len(col) != 2):
                    print "line",i,"has error"
                    exit()
                else:
                    haloid = long(col[0])
                    nprog = long(col[1])
                    count = nprog
            else:
                if(len(col) != 1):
                    print "line",i,"has error"
                    exit()
                else:
                    progid = long(col[0])
                    if(nprog==count):
                        halocopy[haloid]["FirstProgenitor"] = progid
                        halocopy[progid]["Descendant"] = haloid
                        # print haloid,"=>",progid
                    else:
                        halocopy[prevhalo]["NextProgenitor"] = progid
                        halocopy[progid]["Descendant"] = haloid
                        # print prevhalo,"=>",progid
                    prevhalo = progid
                    count -= 1

def treecrowler(hid,halocat,treenr,fulltree):
    halocat[hid]["TreeNr"] = treenr
    halocat[hid]["HaloNr"] = len(fulltree[treenr])
    fulltree[treenr].append(hid)
    progid = halocat[hid]["FirstProgenitor"]
    # print progid
    if progid > -1:
        (halocat,fulltree) = treecrowler(progid,halocat,treenr,fulltree)
    nextprog = halocat[hid]["NextProgenitor"]
    # print "nextprog",nextprog
    # time.sleep(1)
    if nextprog > -1:
        (halocat,fulltree) = treecrowler(nextprog,halocat,treenr,fulltree)
    return (halocat,fulltree)

def outputtrees(halocat2,fileout,fileout2):
    halocat = copy.copy(halocat2)
    ntrees = 0
    cumnhalo = []
    fulltree = {}
    print "start outputting trees"
    for haloid in halocat.iterkeys():
        halo = halocat[haloid]
        if(halo["SnapNum"] == 61) & (halo["MainHalo"] == -1) & (halo["FirstProgenitor"] > -1):
            curid = haloid
            fulltree[ntrees] = []
            while curid > -1:
                (halocat,fulltree) = treecrowler(curid,halocat,ntrees,fulltree)
                curid = halocat[curid]["NextHalo"]
            if len(fulltree[ntrees]) > 0:
                ntrees += 1

    #GROUP TREES INTO BUSHES
    newfulltree = {}
    newntrees = 0
    oldmergetonew = []
    newmergetonew = {}

    for tree in range(ntrees):
        oldmergetonew.append(-1)

    for tree in range(ntrees):
        #print "tree:",tree
        if(len(fulltree[tree]) > 0):
            checked = 0
            newfulltree[newntrees] = []
            newmergetonew[newntrees] = -1
            for hids in fulltree[tree]:
                newfulltree[newntrees].append(hids)
            fulltree[tree] = []
            oldmergetonew[tree] = newntrees
            # print "move oldtree:",tree,"=>",newntrees
            countrep = 1
            while checked == 0:
               #print "Check tree ..."
                count = 0
                maptree = {}
                maptree[-1] = -1
                insidecheck = 1 # 1 if it's OK
                for hid in newfulltree[newntrees]:
                    maptree[hid] = count
                    count += 1
                # print tree,": trial:",countrep
                if len(maptree) != len(newfulltree[newntrees])+1:
                    print "INSIDE: dupplicate"
                    exit()
                for hid in newfulltree[newntrees]:
                    halo = halocat[hid]
                    if halo["MainHalo"] not in maptree: # -1 is in maptree
                        target = halo["MainHalo"]
                        oldtree = halocat[target]["TreeNr"]
                        if(oldtree > -1):
                            if(oldmergetonew[oldtree] == -1):
                                for hids in fulltree[oldtree]:
                                    newfulltree[newntrees].append(hids)
                                fulltree[oldtree] = []
                                oldmergetonew[oldtree] = newntrees
                                # print "move oldtree:",oldtree,"=>",newntrees
                            else:
                                srctree = oldmergetonew[oldtree]
                                reftree = srctree
                                while srctree > -1:
                                    reftree = srctree
                                    srctree = newmergetonew[srctree]
                                srctree = reftree
                                for hids in newfulltree[srctree]:
                                    newfulltree[newntrees].append(hids)
                                # print "move newtree:",srctree,"=>",newntrees
                                newmergetonew[srctree] = newntrees
                                newfulltree[srctree] = []
                        else:
                            # print "remove main halo"
                            halo["NextHalo"] = -1
                            halo["MainHalo"] = -1
                        insidecheck = 0
                        break
                    if halo["NextHalo"] not in maptree: # -1 is in maptree
                        target = halo["NextHalo"]
                        oldtree = halocat[target]["TreeNr"]
                        if(oldtree > -1):
                            if(oldmergetonew[oldtree] == -1):
                                for hids in fulltree[oldtree]:
                                    newfulltree[newntrees].append(hids)
                                fulltree[oldtree] = []
                                oldmergetonew[oldtree] = newntrees
                                print "move oldtree:",oldtree,"=>",newntrees
                            else:
                                srctree = oldmergetonew[oldtree]
                                reftree = srctree
                                while srctree > -1:
                                    reftree = srctree
                                    srctree = newmergetonew[srctree]
                                srctree = reftree
                                for hids in newfulltree[srctree]:
                                    newfulltree[newntrees].append(hids)
                                print "move newtree:",srctree,"=>",newntrees
                                newmergetonew[srctree] = newntrees
                                newfulltree[srctree] = []
                        else:
                            # print "forward halo"
                            halo["NextHalo"] = halocat[target]["NextHalo"]
                        insidecheck = 0
                        break
                countrep += 1
                if(insidecheck == 1):
                    checked = 1
            newntrees += 1

    nhalopertree = []
    nhalos = 0
    ntrees = 0
    fulltree = {}
    for tree in range(newntrees):
        if len(newfulltree[tree]) > 0:
            nhalopertree.append(len(newfulltree[tree]))
            fulltree[ntrees] = newfulltree[tree]
            nhalos += len(newfulltree[tree])
            ntrees += 1

    # check uniqueness
    count = 0
    unique = {}
    for tree in range(ntrees):
        for hid in fulltree[tree]:
            unique[hid] = count
            count += 1

    if count != len(unique):
        print "trees are not isolated"
        exit()

    fp = open(fileout,"wb")
    fp2 = open(fileout2,"wb")
    print "Ntrees:",ntrees
    buffer = struct.pack("i",int(ntrees))
    fp.write(buffer)
    print "Nhalos:",nhalos
    buffer = struct.pack("i",int(nhalos))
    fp.write(buffer)
    for tree in range(ntrees):
        #print tree,":",nhalopertree[tree]
        buffer = struct.pack("i",int(nhalopertree[tree]))
        fp.write(buffer)

    maptreeall = {}
    maptreeall[-1] = -1
    countall = 0
    for tree in range(ntrees):
        count = 0
        maptree = {}
        maptree[-1] = -1
        for hid in fulltree[tree]:
            maptree[hid] = count
            maptreeall[hid] = countall
            count += 1
            countall += 1
        if len(fulltree[tree])+1 != len(maptree):
            print "There are dupplicated entries"
            print "Exit"
            exit()
        # print "tree",tree,"total halos",len(fulltree[tree])
        for hid in fulltree[tree]:
            halo = halocat[hid]
            # if(halo["NextProgenitor"] != -1):
            #     print "prog:",int(maptree[halo["FirstProgenitor"]])
            #     print "next:",int(maptree[halo["NextProgenitor"]])
            #     time.sleep(1)
            buffer = struct.pack("i",int(maptree[halo["Descendant"]]))
            fp.write(buffer)
            buffer = struct.pack("i",int(maptree[halo["FirstProgenitor"]]))
            fp.write(buffer)
            buffer = struct.pack("i",int(maptree[halo["NextProgenitor"]]))
            fp.write(buffer)
            if(halo["MainHalo"] > -1):
                buffer = struct.pack("i",int(maptree[halo["MainHalo"]]))
            else:
                buffer = struct.pack("i",int(maptree[halo["ID"]]))
            fp.write(buffer)
            if halo["NextHalo"] not in maptree:
                print tree
                print halo
                print halocat[halo["NextHalo"]]
                print fulltree[tree]
                exit()
            buffer = struct.pack("i",int(maptree[halo["NextHalo"]]))
            fp.write(buffer)
            buffer = struct.pack("i",halo["Len"])
            fp.write(buffer)
            buffer = struct.pack("f",halo["M_200Mean"])
            fp.write(buffer)
            buffer = struct.pack("f",halo["M_200Crit"])
            fp.write(buffer)
            buffer = struct.pack("f",halo["M_TopHat"])
            fp.write(buffer)
            buffer = struct.pack("fff",halo["Pos"][0],halo["Pos"][1],halo["Pos"][2])
            fp.write(buffer)
            buffer = struct.pack("fff",halo["Vel"][0],halo["Vel"][1],halo["Vel"][2])
            fp.write(buffer)
            buffer = struct.pack("f",halo["VelDisp"])
            fp.write(buffer)
            buffer = struct.pack("f",halo["Vmax"])
            fp.write(buffer)
            buffer = struct.pack("fff",halo["Spin"][0],halo["Spin"][1],halo["Spin"][2])
            fp.write(buffer)
            buffer = struct.pack("q",0) #Mostboundid
            fp.write(buffer)
            buffer = struct.pack("i",halo["SnapNum"])
            fp.write(buffer)
            buffer = struct.pack("i",0) #FileNr
            fp.write(buffer)
            buffer = struct.pack("i",0) #Subhaloindex
            fp.write(buffer)
            buffer = struct.pack("f",0.0) #Subhalfmass
            fp.write(buffer)


    fp.close()
    for tree in range(ntrees):
        lastprogid = fulltree[tree][len(fulltree[tree])-1]
        for hid in fulltree[tree]:
            halo = halocat[hid]
            buffer = struct.pack("q",halo["UID"])
            fp2.write(buffer)
            buffer = struct.pack("q",0)  #FileNr
            fp2.write(buffer)
            buffer = struct.pack("q",maptreeall[halo["FirstProgenitor"]])
            fp2.write(buffer)
            buffer = struct.pack("q",maptreeall[lastprogid])
            fp2.write(buffer)
            buffer = struct.pack("q",maptreeall[halo["NextProgenitor"]])
            fp2.write(buffer)
            buffer = struct.pack("q",maptreeall[halo["Descendant"]])
            fp2.write(buffer)
            if(halo["MainHalo"] > -1):
                buffer = struct.pack("q",maptreeall[halo["MainHalo"]])
            else:
                buffer = struct.pack("q",maptreeall[halo["ID"]])
            fp2.write(buffer)
            buffer = struct.pack("q",maptreeall[halo["NextHalo"]])
            fp2.write(buffer)
            buffer = struct.pack("d",halo["Redshift"])
            fp2.write(buffer)
            buffer = struct.pack("i",0) #Peano keys
            fp2.write(buffer)
            buffer = struct.pack("i",0) #dummy
            fp2.write(buffer)
            
    fp2.close()  


#;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


halo = readAHFascii()
ahf = readSussingtree(SUSSINGtree,halo)
outputtrees(ahf,FileOut,FileOut2)

