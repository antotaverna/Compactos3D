import illustris_python as il
import numpy  as np

basePath = '/big4/simul/Illustris-TNG300/output'

#leo SubHalos
fields = ['SubhaloFlag',
          'SubhaloGrNr',
          'SubhaloCM',
          'SubhaloMass',
          'SubhaloHalfmassRad',
          'SubhaloHalfmassRadType',
          'SubhaloParent',
          'SubhaloPos',
          'SubhaloVel',
          'SubhaloVelDisp',
          'SubhaloSFR',
          'SubhaloStellarPhotometrics',
          'SubhaloStellarPhotometricsMassInRad',
          'SubhaloStellarPhotometricsRad',
          'SubhaloMassInRadType',
          ]

subhalos = il.groupcat.loadSubhalos(basePath,99,fields=fields)
print('subhalos keys:',subhalos.keys())
print('subhalos.shape=',subhalos['SubhaloMass'].shape)


#leo Halos
#Ger: En los "fields" del load.Halos le pones el "GroupFirstSub" que te la el nombre (o posición en la lista) del primer Subhalo del grupo (que es la glx central)
#Y "GroupNsub" que te dice el número de subhalos que tiene ese halo, Que serían algo así como el número de galaxias en el halo

fields_h = ['GroupFirstSub',
            'GroupNsubs',
            'GroupCM',
            'Group_M_Crit200',
            'Group_M_Mean200',
            'GroupPos',
            'GroupVel',
            'Group_R_Crit200',
            'Group_R_Mean200']

halos = il.groupcat.loadHalos(basePath,99,fields=fields_h)
halos['Group_M_Crit200'].shape

n_subs = halos['GroupNsubs']  # Number subhalos
first = halos['GroupFirstSub'] # Glx central- first subhalos
grnr = subhalos['SubhaloGrNr'] # index into the Group table of the FOF host/parent of this Subhalo.
mag = subhalos['SubhaloStellarPhotometrics']

#link --> GroupFirstSub(halo) + SubhaloGrNr(sub halo)
f = open ("galaxycatalogue_cut_Mst.dat", "w")
#for i in range(len(first[:2])):
for i in range(len(first)):
    if(i%10000==0):
        print('****', i)

    #print(len(grnr[ii]),n_subs[i],i,grnr[ii][0],first[i])

    ii=np.where((subhalos['SubhaloGrNr']==i) & (subhalos['SubhaloStellarPhotometricsMassInRad']>0.001))
    for j in range(len(subhalos['SubhaloCM'][ii])):

        cadena=str(subhalos['SubhaloCM'][ii][j][0])+'  '+str(subhalos['SubhaloCM'][ii][j][1])+'  '+str(subhalos['SubhaloCM'][ii][j][2])+\
            '  '+str(subhalos['SubhaloHalfmassRad'][ii][j])+'  '+str(subhalos['SubhaloMass'][ii][j])+\
            '  '+str(subhalos['SubhaloStellarPhotometricsMassInRad'][ii][j])+'  '+str(subhalos['SubhaloStellarPhotometricsRad'][ii][j])+\
            '  '+str(subhalos['SubhaloParent'][ii][j])+'  '+str(subhalos['SubhaloGrNr'][ii][j])+\
            '  '+str(subhalos['SubhaloPos'][ii][j][0])+'  '+str(subhalos['SubhaloPos'][ii][j][1])+'  '+str(subhalos['SubhaloPos'][ii][j][2])+\
            '  '+str(subhalos['SubhaloSFR'][ii][j])+'  '+str(subhalos['SubhaloVelDisp'][ii][j])+'  '+str(subhalos['SubhaloVel'][ii][j][0])+\
            '  '+str(subhalos['SubhaloVel'][ii][j][1])+'  '+str(subhalos['SubhaloVel'][ii][j][2])+\
            '  '+str(mag[ii][j][0])+'  '+str(mag[ii][j][1])+'  '+str(mag[ii][j][2])+'  '+str(mag[ii][j][3])+\
            '  '+str(mag[ii][j][4])+'  '+str(mag[ii][j][5])+'  '+str(mag[ii][j][6])+'  '+str(mag[ii][j][7])+\
            '  '+str(first[i])+'  '+str(n_subs[i])+'  '+str(halos['Group_M_Crit200'][i])+'  '+str(halos['Group_M_Mean200'][i])+\
            '  '+str(halos['Group_R_Crit200'][i])+'  '+str(halos['Group_R_Mean200'][i])+\
            '  '+str(halos['GroupCM'][i][0])+'  '+str(halos['GroupCM'][i][1])+'  '+str(halos['GroupCM'][i][2])+\
            '  '+str(halos['GroupPos'][i][0])+'  '+str(halos['GroupPos'][i][1])+'  '+str(halos['GroupPos'][i][2])+\
            '  '+str(halos['GroupVel'][i][0])+'  '+str(halos['GroupVel'][i][1])+'  '+str(halos['GroupVel'][i][2])+'\n'
        f.write(cadena)
f.close()
