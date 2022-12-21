echo "Uaage 1: get_theta90.sh
               Read energies.dat, sort by the second column."
echo

echo "Uaage 2: get_theta90.sh Any_argument
               Read energies_average.dat, sort by the second column."
echo

echo "Uaage 3: get_theta90.sh file_name column_number
               Read file_name, sort by the column column_number."
echo

if [ $# == 0 ]; then
  grep "^  90.0" energies.dat | sort -n -k2  > therta90.dat 
elif [ $# == 1 ]; then
  grep "^  90.0" energies_average.dat | sort -n -k2  > therta90.dat 
else
  grep "^  90.0" $1 | sort -n -k $2  > therta90.dat 
fi
