#!/bin/bash

topDir=jobs
jobFile=SRValidData_HI2023_R375790_cfg.py
#jobFile2=SRValidData_HI_R362297_TEMPLATE_cfg.py
dbFile=srCondWrite_TEMPLATE_cfg.py

#FR SETTINGS
# doSRARR=(1 1 1 1)
# srLTHARR=(10.0 10.0 8.0 10.0)
# srHTHARR=(14.0 14.0 10.0 14.0)

# doZSARR=(1 1 0 0)
# doZSCenterARR=(0 0 0 0)
# zsMIEBARR=(8.0 8.0 8.0 8.0)
# zsMIEEARR=(8.0 8.0 8.0 8.0)

# zsHIEBARR=(2.5 3.5 0.0 0.0)
# zsHIEEARR=(3.0 4.0 0.0 0.0)

#Recent 2023.09.16
#doSRARR=(1 1 1 1)
#srLTHARR=(8.0 8.0 8.0 8.0)
#srHTHARR=(10.0 10.0 10.0 10.0)

#doZSARR=(1 1 1 1)
#doZSCenterARR=(0 0 0 0)
#zsMIEBARR=(6.0 7.0 8.0 10.0)
#zsMIEEARR=(6.0 7.0 8.0 10.0)

#zsHIEBARR=(6.0 7.0 8.0 10.0)
#zsHIEEARR=(6.0 7.0 8.0 10.0)

#END FR SETTINGS

#doSRARR=(1 1 1 1)
#srLTHARR=(4.0 6.0 7.0 8.0) # this is 2xGeV so 2 GeV
#srHTHARR=(8.0 12.0 14.0 16.0) # this is 2xGeV so 4 GeV

#doSRARR=(1 1 1 1)
#srLTHARR=(5.0 5.0 5.0 5.0)
#srHTHARR=(10.0 10.0 10.0 10.0)

#doZSARR=(1 1 1 1)
#zsMIEBARR=(6.0 8.0 12.0 16.0)
#zsMIEEARR=(6.0 8.0 12.0 16.0)

#doZSARR=(0 0 0 0)
#zsMIEBARR=(5.5 5.5 5.5 5.5)
#zsMIEEARR=(6.0 6.0 6.0 6.0)

# test only 1 conf
# doSRARR=(1)
# srLTHARR=(4.0)
# srHTHARR=(8.0)

# doZSARR=(0 0 0 0)
# doZSCenterARR=(0 0 0 0)
# zsMIEBARR=(8.0 8.0 8.0 8.0)
# zsMIEEARR=(8.0 8.0 8.0 8.0)

# zsHIEBARR=(2.5 3.5 0.0 0.0)
# zsHIEEARR=(3.0 4.0 0.0 0.0)

# Create TARGET and SAFE
# doSRARR=(1 1 1) # 1 to DoSR
# srLTHARR=(4.0 8.0 10.0)
# srHTHARR=(8.0 10.0 14.0)

# doZSARR=(0 0 0) # 1 to DoZS 
# doZSCenterARR=(0 0 0)
# zsMIEBARR=(8.0 8.0 8.0)
# zsMIEEARR=(8.0 8.0 8.0)

# zsHIEBARR=(2.5 3.5 0.0)
# zsHIEEARR=(3.0 4.0 0.0)

# EXTRA SAFE tests
doSRARR=(1 1 1) # 1 to DoSR
srLTHARR=(4.0 8.0 10.0)
srHTHARR=(8.0 10.0 14.0)

doZSARR=(0 0 0) # 1 to DoZS 
doZSCenterARR=(0 0 0)
zsMIEBARR=(10.0 10.0 10.0)
zsMIEEARR=(10.0 10.0 10.0)

zsHIEBARR=(2.5 3.5 0.0)
zsHIEEARR=(3.0 4.0 0.0)

pos=0
for i in "${doSRARR[@]}"
do
 
    doSR=${doSRARR[$pos]}
    srLTH=${srLTHARR[$pos]}
    srHTH=${srHTHARR[$pos]}
    
    doZS=${doZSARR[$pos]}
    doZSCenter=${doZSCenterARR[$pos]}
    zsMIEB=${zsMIEBARR[$pos]}
    zsMIEE=${zsMIEEARR[$pos]}
        
    zsHIEB=${zsHIEBARR[$pos]}
    zsHIEE=${zsHIEEARR[$pos]}
    
    jobStr="NoSR"
    
    if [[ $doSR -eq 1 ]]
    then
	jobStr="DoSR_LTH"$srLTH"_HTH"$srHTH
    fi
    
    jobStr2="NoZS_MIEB"$zsMIEB"_MIEE"$zsMIEE
    
    if [[ $doZS -eq 1 ]]
    then
	jobStr2="DoZS_MIEB"$zsMIEB"_MIEE"$zsMIEE"_HIEB"$zsHIEB"_HIEE"$zsHIEE"_CenterTow"$doZSCenter
    fi
    
    jobStr=$jobStr"_"$jobStr2
    
    while [[ $jobStr == *"."* ]]
    do
	jobStr=$(echo $jobStr | sed -e "s@\.@p@g")
    done
    
    fileOut=$(echo $jobFile | sed -e "s@TEMPLATE@$jobStr@g")
#    fileOut2=$(echo $jobFile2 | sed -e "s@TEMPLATE@$jobStr@g")
    
    dirOut=$topDir/$jobStr
    fileOut=$dirOut/$fileOut
    fileOut2=$dirOut/$fileOut2
    
    dbFileNewName=$(echo $dbFile | sed -e "s@TEMPLATE@$jobStr@")
    dbFileOut=$dirOut/$dbFileNewName
    
    mkdir -p $dirOut
    cp $jobFile $fileOut
#    cp $jobFile2 $fileOut2
    cp $dbFile $dbFileOut
    
    sed -i "s@TEMPLATE@$jobStr@g" $dbFileOut
    sed -i "s@TEMPLATE@$jobStr@g" $fileOut
#    sed -i "s@TEMPLATE@$jobStr@g" $fileOut2
    
    if [[ $doSR -eq 1 ]]
    then
	echo "" >> $dbFileOut
	echo "process.writeInDB.trigPrimBypass = cms.bool(True)" >> $dbFileOut
	echo "process.writeInDB.trigPrimBypassMode = cms.int32(1)" >> $dbFileOut
	echo "process.writeInDB.trigPrimBypassHTH = cms.double($srHTH)" >> $dbFileOut
	echo "process.writeInDB.trigPrimBypassLTH = cms.double($srLTH)" >> $dbFileOut
	echo "" >> $dbFileOut
	
	echo "" >> $fileOut
	echo "simEcalDigis.trigPrimBypass = cms.bool(True)" >> $fileOut
	echo "simEcalDigis.trigPrimBypassMode = cms.int32(1)" >> $fileOut
	echo "simEcalDigis.trigPrimBypassLTH = cms.double($srLTH)" >> $fileOut
	echo "simEcalDigis.trigPrimBypassHTH = cms.double($srHTH)" >> $fileOut
	echo "" >> $fileOut

#	echo "" >> $fileOut2
#	echo "simEcalDigis.trigPrimBypass = cms.bool(True)" >> $fileOut2
#	echo "simEcalDigis.trigPrimBypassMode = cms.int32(1)" >> $fileOut2
#	echo "simEcalDigis.trigPrimBypassLTH = cms.double($srLTH)" >> $fileOut2
#	echo "simEcalDigis.trigPrimBypassHTH = cms.double($srHTH)" >> $fileOut2
#	echo "" >> $fileOut2
    fi
    
    if [[ $doZS -eq 1 ]]
    then
	echo "" >> $dbFileOut
	echo "process.writeInDB.srpBarrelHighInterestChannelZS = cms.double($zsHIEB*0.035)" >> $dbFileOut
	echo "process.writeInDB.srpEndcapHighInterestChannelZS = cms.double($zsHIEE*0.06)" >> $dbFileOut
	echo "" >> $dbFileOut

	if [[ $doZSCenter -eq 1 ]]
	then
	    echo "process.writeInDB.actions = cms.vint32(1, 2, 2, 2, 5, 6, 6, 6)" >> $dbFileOut
	else
	    echo "process.writeInDB.actions = cms.vint32(1, 2, 2, 3, 5, 6, 6, 7)" >> $dbFileOut
	fi

	echo "" >> $dbFileOut
    fi

    doZSOrSR=0
    if [[ $doZS -eq 1 ]]
    then
	doZSOrSR=1
    elif [[ $doSR -eq 1 ]]
    then
	doZSOrSR=1
    fi

    if [[ $doZSOrSR -eq 1 ]]
    then
	echo "" >> $fileOut
	#    echo "process.load(\"CondCore.CondDB.CondDB_cfi\")" >> $fileOut
	echo "process.GlobalTag.toGet = cms.VPSet(" >> $fileOut
	echo " cms.PSet(record = cms.string(\"EcalSRSettingsRcd\")," >> $fileOut
	echo "   tag = cms.string(\"EcalSRSettings_$jobStr\")," >> $fileOut
	echo "   connect = cms.string(\"sqlite_file:EcalSRSettings_$jobStr.db\")" >> $fileOut
	echo " )" >> $fileOut
	echo ")" >> $fileOut
	echo "" >> $fileOut

#	echo "" >> $fileOut2
	#    echo "process.load(\"CondCore.CondDB.CondDB_cfi\")" >> $fileOut2
#	echo "process.GlobalTag.toGet = cms.VPSet(" >> $fileOut2
#	echo " cms.PSet(record = cms.string(\"EcalSRSettingsRcd\")," >> $fileOut2
#	echo "   tag = cms.string(\"EcalSRSettings_$jobStr\")," >> $fileOut2
#	echo "   connect = cms.string(\"sqlite_file:EcalSRSettings_$jobStr.db\")" >> $fileOut2
#	echo " )" >> $fileOut2
#	echo ")" >> $fileOut2
#	echo "" >> $fileOut2

	sed -i "s@EBLIZS@$zsMIEB@" $dbFileOut
	sed -i "s@EELIZS@$zsMIEE@" $dbFileOut

	cd $dirOut   
	cmsRun $dbFileNewName
	cd ../../
    fi
    
    pos=$((pos+1))
done
