editRundeck()
# -------------------------------------------------------------------
{
# Adapted from R. Ruedy's script add_regr_test_params
# This function enables regression testing of restarts. It modifies
# the RUNDECK by adding an NDISK line right before the &&END_PARAMETERS
# and an end time near the end of the RUNDECK.

  local deck=$1".R"
  local ndisk=$2
  local datee=$3
  local houre=$4

  echo $deck
  ndisk_line="ndisk=$ndisk"
  end_hour_line=" DATEE=$datee, HOURE=$houre,"
  eof1='&&END_PARAMETERS'
  eof2='/'

# this section may be omitted, once all old style rundecks are modified
  a=$( grep ^' '${eof2} $deck )
  if [[ $a = '' ]];then
     eof2='&END'
  fi

  a=$( grep -n ${eof1}     $deck | head -1 ) ; n1=${a%%:*}
  a=$( grep -n ^' '${eof2} $deck | head -1 ) ; n2=${a%%:*}

  cp ${deck} templ
  head -$(( n1-1 )) templ                   > ${deck}
  echo "${ndisk_line}"                      >> ${deck}
  tail +${n1} templ | head -$(( n2 - n1 ))  >> ${deck}
#  eline=$( tail +${n1} templ | head -$(( n2 - n1 )) | grep -i istart | tail -1 )
#  echo  ' '${eline#*ISTART*,}               >> ${deck}
  echo "${end_hour_line}"                   >> ${deck}
  echo " ISTART=2, ${end_hour_line}"         >> ${deck}
  tail +${n2} templ                         >> ${deck}
}

editRundeck $1 $2 $3 $4;
