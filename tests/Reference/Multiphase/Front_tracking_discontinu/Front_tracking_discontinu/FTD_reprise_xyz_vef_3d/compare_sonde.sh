# Chutte d'une goutte
echo Chute d\'une goutte sous gravite g=[0. 0. -10.]
cat <<END | gnuplot
set output "compare_sonde.eps"
set terminal postscript color eps
plot \
-10.*x title "vitesse theorique",\
"GOUTTE_reference.son" using 1:4 title "calcul de reference" ps 2,\
"prepare_V_GOUTTE_reference.son" using 1:4 title "reference prepare" ps 3,\
"prepare_V_GOUTTE.son" using 1:4 title "resultat prepare.data",\
"FTD_reprise_xyz_vef_3d_V_GOUTTE_reference.son" using 1:4 title "reference seq", \
"FTD_reprise_xyz_vef_3d_V_GOUTTE.son" using 1:4 title "resultat seq"
END
echo Pour voir le fichier resultat: gv compare_sonde.eps
gv compare_sonde.eps
