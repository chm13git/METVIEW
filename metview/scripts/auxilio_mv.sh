#!/bin/bash -x

# Testando se o paramentro referente o horario da rodada foi passado

if [ $# -ne 1 ];then

        echo "Entre com a data e com o horario de referencia (00 ou 12)."
        echo "Ex.: auxilio_mv.sh 00"
        exit

fi

DATA=`date +%Y%m%d`
HH=$1
#DATA=20210628

echo Inicio: `date`
echo
echo Rodando script auxilio MV para a rodada de $HH de $DATA!
sleep 2

cd /home/operador/metview/neris/work

metview -b /home/operador/metview/neris/nebrasil_emborg_cosmo.mv $DATA $HH
metview -b /home/operador/metview/neris/nebrasil_opeanv_cosmo.mv $DATA $HH
#metview -b nebrasil_emborg_cosmo.mv $DATA $HH
#metview -b nebrasil_opeanv_cosmo.mv $DATA $HH

echo "Transferindo figuras para o Supervisor (DPNT02B e DPNT02C)"
scp auxilio_mv_nebrasil_*_${HH}_*.png supervisor@10.13.200.5:~/grads/gif/auxilio$HH/cosmomet/ww3met
scp auxilio_mv_nebrasil_*_${HH}_*.png supervisor@10.13.200.8:~/grads/gif/auxilio$HH/cosmomet/ww3met

echo Final: `date`
