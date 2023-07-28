#! /bin/csh

##################################################################
# En principio corro todo para z=0
##################################################################
# Calculo b0 para 3 cosmologias y todos los z (over=fija?)
gfortran -o overdensity overdensity.f

# Identifico grupos fof en z0/
cc -lm -o idenanalit idenanalit.c

# Genero box con grupos sin repetidas
gfortran -o readgroups2 readgroups2.f

# Calculo las propiedades de los grupos
gfortran -mcmodel=large -o masavir2 masavir2.f

# DIFUSOS/
# readgroups2
# masavir2

# check subestructura
gfortran -o cross_check_sub_new_euge cross_check_sub_new_euge.f

# Genero catalogo 3-D CGs
gfortran -o aislados aislados.f

# Genero mock catalogo 3-D CGs
gfortran -mcmodel=large -o mock_ais mock_ais.f

# Propiedades catalogo 3-D CGs (oct mock y oct mock mlim +3, +1)
gfortran -mcmodel=large -o props_ais props_ais.f

# Propiedades catalogo 3-D densos y no-substr (oct mock y oct mock mlim +3)
#gfortran -mcmodel=large -o props_densos props_densos.f

############################ overfensity.f
echo --------------------------------------------
echo --------------------------------------------
echo overdensity ncosmo=1
#./overdensity 1
echo overdensity ncosmo=2
#./overdensity 2
echo overdensity ncosmo=3
#./overdensity 3
echo overdensity ncosmo=4
#./overdensity 4
echo overdensity ncosmo=7
#./overdensity 7

############################ idenanalit.c
echo --------------------------------------------
echo --------------------------------------------
echo --------------------------------------------
echo idenanalit ncosmo=1, z=0 ---------------------
#./idenanalit 1 0
echo idenanalit ncosmo=2, z=0 ---------------------
#./idenanalit 2 0
echo idenanalit ncosmo=3, z=0 ---------------------
#./idenanalit 3 0
echo idenanalit ncosmo=4, z=0 ---------------------
#./idenanalit 4 0
echo idenanalit ncosmo=7, z=0 ---------------------
#./idenanalit 7 0

############################ readgroups2.f
echo --------------------------------------------
echo --------------------------------------------
echo --------------------------------------------
echo readgroups2 ncosmo=1, z=0 ---------------------
#./readgroups2 1 0 
echo readgroups2 ncosmo=2, z=0 ---------------------
#./readgroups2 2 0 
echo readgroups2 ncosmo=3, z=0 ---------------------
#./readgroups2 3 0 
echo readgroups2 ncosmo=4, z=0 ---------------------
#./readgroups2 4 0 
echo readgroups2 ncosmo=7, z=0 ---------------------
#./readgroups2 7 0 

############################ masavir2.f
echo --------------------------------------------
echo --------------------------------------------
echo --------------------------------------------
echo masavir2 ncosmo=1, z=0 ---------------------
#./masavir2 1 0
echo masavir2 ncosmo=2, z=0 ---------------------
#./masavir2 2 0
echo masavir2 ncosmo=3, z=0 ---------------------
#./masavir2 3 0
echo masavir2 ncosmo=4, z=0 ---------------------
#./masavir2 4 0
echo masavir2 ncosmo=7, z=0 ---------------------
#./masavir2 7 0
echo masavir2 ncosmo=8, z=0 ---------------------
#./masavir2 8 0

############################ cross_check_sub_new_euge.f
echo --------------------------------------------
echo --------------------------------------------
echo --------------------------------------------
echo check substr ncosmo=1, z=0 ---------------------
#./cross_check_sub_new_euge 1
echo check substr ncosmo=2, z=0 ---------------------
#./cross_check_sub_new_euge 2
echo check substr ncosmo=3, z=0 ---------------------
#./cross_check_sub_new_euge 3
echo check substr ncosmo=4, z=0 ---------------------
#./cross_check_sub_new_euge 4
echo check substr ncosmo=7, z=0 ---------------------
#./cross_check_sub_new_euge 7
echo check substr ncosmo=8, z=0 ---------------------
#./cross_check_sub_new_euge 8

############################ aislados.f
echo --------------------------------------------
echo --------------------------------------------
echo --------------------------------------------
echo Aislados ncosmo=1, z=0 ---------------------
#./aislados 1 0
echo Aislados ncosmo=2, z=0 ---------------------
#./aislados 2 0
echo Aislados ncosmo=3, z=0 ---------------------
#./aislados 3 0
echo Aislados ncosmo=4, z=0 ---------------------
#./aislados 4 0
echo Aislados ncosmo=7, z=0 ---------------------
#./aislados 7 0
echo Aislados ncosmo=8, z=0 ---------------------
#./aislados 8 0

############################ mock_ais.f
echo --------------------------------------------
echo --------------------------------------------
echo --------------------------------------------
echo mock_aislados_y_densos ncosmo=1, z=0, wh_gr=0-1-2-3------------------
#./mock_ais 1 0 0 #ais
#./mock_ais 1 0 1 #dens
##./mock_ais 1 0 2 #dens_no_subs
#./mock_ais 1 0 3 #difusos
echo mock_aislados_y_densos ncosmo=2, z=0, wh_gr=0-1-2-3------------------
#./mock_ais 2 0 0
#./mock_ais 2 0 1 #dens
##./mock_ais 2 0 2 #dens_no_subs
#./mock_ais 2 0 3 #difusos
echo mock_aislados_y_densos ncosmo=3, z=0, wh_gr=0-1-2-3------------------
#./mock_ais 3 0 0
#./mock_ais 3 0 1 #dens
##./mock_ais 3 0 2 #dens_no_subs
#./mock_ais 3 0 3 #difusos
echo mock_aislados_y_densos ncosmo=4, z=0, wh_gr=0-1-2-3-------------------
#./mock_ais 4 0 0
#./mock_ais 4 0 1 #dens
##./mock_ais 4 0 2 #dens_no_subs
#./mock_ais 4 0 3 #difusos
echo mock_aislados_y_densos ncosmo=7, z=0, wh_gr=0-1-2-3------------------
#./mock_ais 7 0 0
#./mock_ais 7 0 1 #dens
##./mock_ais 7 0 2 #dens_no_subs
#./mock_ais 7 0 3 #difusos
echo mock_aislados_y_densos ncosmo=8, z=0, wh_gr=0-1-2-3------------------
#./mock_ais 8 0 0
#./mock_ais 8 0 1 #dens
##./mock_ais 8 0 2 #dens_no_subs
#./mock_ais 8 0 3 #difusos

############################
echo --------------------------------------------
echo --------------------------------------------
echo --------------------------------------------
echo props_aislados ncosmo=1, z=0, mem_lim=0-1-3 ---------------------
#./props_ais 1 0 0 #oct mock
./props_ais 1 0 1 #oct mock mlim +1
./props_ais 1 0 3 #oct mock mlim +3
echo props_aislados ncosmo=2, z=0, mem_lim=0-1-3 ---------------------
#./props_ais 2 0 0 #oct mock
./props_ais 2 0 1 #oct mock mlim +1
./props_ais 2 0 3 #oct mock mlim +3
echo props_aislados ncosmo=3, z=0, mem_lim=0-1-3 ---------------------
#./props_ais 3 0 0 #oct mock
./props_ais 3 0 1 #oct mock mlim +1
./props_ais 3 0 3 #oct mock mlim +3
echo props_aislados ncosmo=4, z=0, mem_lim=0-1-3 ---------------------
#./props_ais 4 0 0 #oct mock
./props_ais 4 0 1 #oct mock mlim +1
./props_ais 4 0 3 #oct mock mlim +3
echo props_aislados ncosmo=7, z=0, mem_lim=0-1-3 ---------------------
#./props_ais 7 0 0 #oct mock
./props_ais 7 0 1 #oct mock mlim +1
./props_ais 7 0 3 #oct mock mlim +3
echo props_aislados ncosmo=8, z=0, mem_lim=0-1-3 ---------------------
#./props_ais 7 0 0 #oct mock
./props_ais 8 0 1 #oct mock mlim +1
./props_ais 8 0 3 #oct mock mlim +3

############################
echo --------------------------------------------
echo --------------------------------------------
echo --------------------------------------------
echo props_densos_y_no_sbrt ncosmo=1, wh_gr=1-2, mem_lim=0-3 ---------------------
#./props_densos 1 1 0 #oct mock densos
#./props_densos 1 1 3 #oct mock densos mlim +3
#./props_densos 1 2 0 #oct mock nosubst
#./props_densos 1 2 3 #oct mock nosust mlim +3
echo props_densos_y_no_sbrt ncosmo=2, wh_gr=1-2, mem_lim=0-3 ---------------------
#./props_densos 2 1 0 #oct mock densos
#./props_densos 2 1 3 #oct mock densos mlim +3
#./props_densos 2 2 0 #oct mock nosubst
#./props_densos 2 2 3 #oct mock nosust mlim +3
echo props_densos_y_no_sbrt ncosmo=3, wh_gr=1-2, mem_lim=0-3 ---------------------
#./props_densos 3 1 0 #oct mock densos
#./props_densos 3 1 3 #oct mock densos mlim +3
#./props_densos 3 2 0 #oct mock nosubst
#./props_densos 3 2 3 #oct mock nosust mlim +3
echo props_densos_y_no_sbrt ncosmo=4, wh_gr=1-2, mem_lim=0-3 ---------------------
#./props_densos 4 1 0 #oct mock densos
#./props_densos 4 1 3 #oct mock densos mlim +3
#./props_densos 4 2 0 #oct mock nosubst
#./props_densos 4 2 3 #oct mock nosust mlim +3
echo props_densos_y_no_sbrt ncosmo=7, wh_gr=1-2, mem_lim=0-3 ---------------------
#./props_densos 7 1 0 #oct mock densos
#./props_densos 7 1 3 #oct mock densos mlim +3
#./props_densos 7 2 0 #oct mock nosubst
#./props_densos 7 2 3 #oct mock nosust mlim +3


############################
#source("pru_dij_g11.r")






##############################################################
##### PARA TODOS LOS Z:[0-27]
##############################################################

# Calculo b0 para 3 cosmologias y todos los z

#gfortran -o overdensity overdensity.f
#./overdensity 1
#./overdensity 2
#./overdensity 3
#################################

# Genero directorios

#@ k=0
#while ($k < 28)
#	rm -rf guo_11/z$k
#	rm -rf guo_13/z$k
#	rm -rf hen_15/z$k
#
#	mkdir guo_11/z$k
#	mkdir guo_13/z$k
#	mkdir hen_15/z$k
#@ k = $k + 1
#end
#################################


#@ j=1
#@ i=0
#while ($j < 4)
#
#   echo overdensity ncosmo= $j
#   ./overdensity $j
#  
#   while ($i < 28)
#
#	echo idenanalit ncosmo= $j, z= $i ---------------------
#	./idenanalit $j $i
#
#	echo readgroups2 ncosmo= $j, z= $i ---------------------
#        echo -------- limpieza repetidas
#	./readgroups2 $j $i   #utiliza archivos "****" 
#
#	echo masavir2 ncosmo= $j, z= $i ---------------------
#	./masavir2 $j $i
#
#	echo aislados ncosmo= $j, z= $i ---------------------
#	./aislados $j $i
#
#   @ i = $i + 1
#@ j = $j + 1
#end

