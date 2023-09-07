#!/bin/bash

usr_id_hyperion="adbh716"

ssh_hyperion_fast='ssh -f -N -J adbh716@hpc-jump00.city.ac.uk -N adbh716@hyperion.city.ac.uk'
$ssh_hyperion_fast

ssh_hyperion='ssh -J adbh716@hpc-jump00.city.ac.uk adbh716@hyperion.city.ac.uk'
sftp_hyperion='sftp -J adbh716@hpc-jump00.city.ac.uk adbh716@hyperion.city.ac.uk'

crnt_dir=$PWD

cd $PWD

case_name=$(basename "$PWD")

#====================================================================================
# Mode: either launch new case or restart downloading

echo
echo 'Press <'L'> to launch a new case or <'D'> to restart downloading.'

read -s -n 1 key  # Wait for the user to press a key
case $key in      # Check which key was pressed
    l|L)
        echo "You pressed 'L'. New case setup will commence."
        mode_script="L"
        ;;
    d|D)
        echo "You pressed 'D'. Downloading files of existing case will commence."
        mode_script="D"
        ;;
    *)
        echo "Invalid input. Please press 'L' or 'R' or 'D'."
        exit
        ;;
esac
echo

#====================================================================================

read -e -p "Please enter the number of nodes (max 8 nodes): " -i "" number_of_nodes
echo

echo 'To keep the temp files press <'Y'> or <'N'> to only keep the .plt.'
read -s -n 1 key  # Wait for the user to press a key
case $key in      # Check which key was pressed
  y|Y)
    echo "You pressed 'y'. temp files will be saved."
    save_files="T"
    ;;
  n|N)
    echo "You pressed 'n'. Only plt files will be saved."
    save_files="F"
    ;;
  *)
    echo "Invalid input. Please press 'y' or 'n'."
    ;;
esac
echo

#====================================================================================
# IF for mode
if [ "$mode_script" == "L" ]; then

  #--------------------------------------------------------------------
  # Input data for sim
  while true; do
    read -p "Please enter the name of the case for squeue (6 characters or fewer): " case_name_squeue
    if [ ${#case_name_squeue} -gt 6 ]; then
        echo "The name of the case for squeue is too long (more than 6 characters). Please try again."
    else
        break  # Exit the loop if the string length is valid
    fi
  done
  echo

  echo 'The simulation with local ID: <'$case_name'> with mode: <'$mode_script'> will launch on Hyperion with ID: <'$case_name_squeue'> on user: <'$usr_id_hyperion'> with #'$number_of_nodes 'node(s)'
  echo

  echo 'Press <'Y'> to continue or <'N'> to abort.'

  read -s -n 1 key   # Wait for the user to press a key
  case $key in       # Check which key was pressed
    y|Y)
      echo "You pressed 'y'. Continuing..."
      ;;
    n|N)
      echo "You pressed 'n'. Exiting..."
      exit 1
      ;;
    *)
      echo "Invalid input. Please press 'y' or 'n'."
      ;;
  esac
  echo
  #--------------------------------------------------------------------

  #--------------------------------------------------------------------
  # Creating local directories
  echo 'Creating local directories.'
  echo
  {
  mkdir temp
  mkdir plt

  cp sh.ct temp/
  cp sh.neu temp/

  if [ $number_of_nodes == '1' ]; then
      cp /home/ek/Desktop/cluster_runs/ForestFV_export/ForestFV_V5_export_96 temp/ForestFV_export
  elif [ $number_of_nodes == '2' ]; then
      cp /home/ek/Desktop/cluster_runs/ForestFV_export/ForestFV_V5_export_192 temp/ForestFV_export
  elif [ $number_of_nodes == '3' ]; then
      cp /home/ek/Desktop/cluster_runs/ForestFV_export/ForestFV_V5_export_288 temp/ForestFV_export
  elif [ $number_of_nodes == '4' ]; then
      cp /home/ek/Desktop/cluster_runs/ForestFV_export/ForestFV_V5_export_384 temp/ForestFV_export
  elif [ $number_of_nodes == '5' ]; then
      cp /home/ek/Desktop/cluster_runs/ForestFV_export/ForestFV_V5_export_480 temp/ForestFV_export
  elif [ $number_of_nodes == '6' ]; then
      cp /home/ek/Desktop/cluster_runs/ForestFV_export/ForestFV_V5_export_576 temp/ForestFV_export
  fi
  } &> /dev/null
  #--------------------------------------------------------------------


  #--------------------------------------------------------------------
  # Connecting to Hyperion
  echo 'Connecting to Hyperion.'
  echo
  echo 'Press <'U'> to upload new version of ForestFV or <'C'> to keep current version.'

  read -s -n 1 key  # Wait for the user to press a key
  case $key in      # Check which key was pressed
      u|U)
          echo "You pressed <'U'>. Uploading..."

          cd /home/ek/Desktop/Forest_V5/src/
{
$sftp_hyperion <<EOF
cd Forest_V5/src/
rm c/*
rm c/.vscode/*
rmdir c/.vscode/
rm c/.editorconfig
rmdir c/
put -r c/ ./
EOF
} &> /dev/null
          {
          make_status=$($ssh_hyperion 'cd Forest_V5/run; rm ForestFV_V5; cd ../src; flight start; flight env activate gridware; module add compilers/gcc; make clean; make -j 10;')
          } &> /dev/null

          #Check if successful
          {
          dummy=$($ssh_hyperion '[ "$(ls -A Forest_V5/run)" ] && dummy=1 || dummy=2; echo $dummy;')
          } &> /dev/null
          dummy="${dummy//[$'\t\r\n ']}"
          
          if [ "$dummy" -eq "1" ]; then
              echo "Upload and Compilation successful."    
          elif [ "$dummy" -eq "2" ]; then
              echo "Compilation failed, exiting..."
              exit
          fi

          ;;
      c|C)
          echo "You pressed <'C'>. Version already on cluster will be used."

          #Check if successful
          {
          dummy=$($ssh_hyperion  '[ "$(ls -A Forest_V5/run)" ] && dummy=1 || dummy=2; echo $dummy;')
          } &> /dev/null
          dummy="${dummy//[$'\t\r\n ']}"
          
          if [ "$dummy" -eq "1" ]; then
            echo "Current version of ForestFV found."
          elif [ "$dummy" -eq "2" ]; then
            echo "Current version of ForestFV not found, please compile and run again, exiting..."
            exit
          fi
          ;;
      *)
          echo "Invalid input. Please press <'U'> or <'C'>."
      ;;
  esac
  echo
  #--------------------------------------------------------------------


#--------------------------------------------------------------------
# Creating directories on Hyperion
echo 'Creating directories on Hyperion.'

path_temp=${crnt_dir#*ForestFV_V5}
path_local=${path_temp#*/}

{
$ssh_hyperion -x << EOF
cd Forest_V5/
mkdir -p $path_local
EOF
} &> /dev/null
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Upload sh.* files. [ok]
echo 'Uploading sh.* files on Hyperion.'
echo
{
$sftp_hyperion <<EOF
cd Forest_V5/$path_local
put sh.* ./
EOF
} &> /dev/null
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Copy: 1) ForestFV_V5                      [ok]
#       2) Script for sbatch                [ok]
#       3) script for filename (download)   [ok]
# Launch case                               [ok]
echo 'Launching simulation.'
echo
$ssh_hyperion /bin/bash << EOF
cd Forest_V5/$path_local

cp /users/adbh716/Forest_V5/run/ForestFV_V5 ./
cp /users/adbh716/misc/launch_n"$number_of_nodes" ./"$case_name_squeue"
cp /users/adbh716/misc/script_filename.sh ./
EOF

job_output=$($ssh_hyperion "cd Forest_V5/$path_local && sbatch $case_name_squeue")
job_id=$(echo "$job_output" | awk '{print $4}')
echo "Job submitted with ID: $job_id"
#--------------------------------------------------------------------

elif [ "$mode_script" == "D" ]; then

#-------------------------------------
# find job id

path_temp=${crnt_dir#*ForestFV_V5}
path_local=${path_temp#*/}

latest_output_file=$($ssh_hyperion "cd Forest_V5/$path_local && ls -t output.* | head -n 1")
job_id=$(echo "$latest_output_file" | awk -F'[.]' '{print $2}')

echo "Latest Output File: $latest_output_file"
echo "Job ID: $job_id"
echo

#-------------------------------------

fi

job_status_temp="T"

# Download files
while true; do

    # Check the status of the job using squeue
    {
    job_status=$($ssh_hyperion "squeue -j "$job_id" -h -o "%T"")
    } &> /dev/null

    if ! [ -n "$job_status" ]; then

      if [ "$job_status_temp" == "T" ]
        echo "Job $job_id has finished."
        job_status_temp="F"
      fi
        
        remaining_files=$($ssh_hyperion "cd Forest_V5/$path_local && ls | grep 'F_0*' | wc -l")

        echo "Remaining files found are: $remaining_files"
        if [ "$remaining_files" -eq "1" ]; then
          break
        fi
    fi

    file_name_slash=$($ssh_hyperion 'cd Forest_V5/'$path_local'; ./script_filename.sh; file_name=$(cat filename.txt); echo $file_name;')

    file_name="${file_name_slash%/*}"               

    file_name_trim="$(echo -e "${file_name}" | tr -d '[:space:]')"

    if ! [ "$file_name_trim" == "" ]; then
        echo "Downloading timestep: $file_name_trim"
{
$sftp_hyperion <<EOF
cd Forest_V5/$path_local
get -r $file_name_slash ./
rm $file_name/*
rmdir $file_name
EOF
} &> /dev/null
        mv $file_name/ temp/
        cd temp/$file_name/
        temp_timestep="${file_name%/*}"      
        echo "${temp_timestep:2}" > tstep.dat
        
        echo "Exporting to .dat, timestep: $file_name_trim"
        {
        /home/ek/Local/petsc_intallation/install/bin/mpirun -np 4 ../ForestFV_export sh
        } &> /dev/null

        cd ../
        
        file_dat="field_linked${temp_timestep:2}prt00.dat"
        file_plt="field_linked${temp_timestep:2}prt00.plt"
        
        echo "Translating to .plt, timestep: $file_name_trim"
        echo
        {
        /home/ek/Desktop/tecplot_new/360ex_2021r1/bin/preplot "$file_dat" "$file_plt"
        } &> /dev/null

        if [ "$save_files" == "F" ]; then
            rm $file_dat
        fi

        mv *.plt ../plt/
        
        cd ../

        if [ "$save_files" == "F" ]; then
            cd temp/
            rm -r $file_name
            cd ../
        fi

        

    fi

    sleep 60

done


