#!/bin/bash
 
#SBATCH --job-name=VSCODE_on_cnode

#SBATCH --time=08:00:00 #This is only an example time parameter but you can set this as you need

#SBATCH --output=%x_%j_%N.log

#SBATCH --mem=3gb #This is only an example memory parameter but you can set this as you need

#SBATCH --cpus-per-task=1 #This is only an example cpu parameter but you can set this as you need

#SBATCH --partition=cpuq  #This is only an example partition parameter but you can set this as you need
 
# Load the vscode 4.89.1 module

module load code-server/4.89.1
 
# This is a very important variable that create the code-server directory in which will be saved all the logs, the extensions and the cache of the sessions. If you don't set this variable, code-server will overwrite and use the cache of the code-server in Posit Workbench! You can insert any folder that you want but you can't choose the $HOME/.local/share/ folder. Here,  code-server will create this folder automatically in your $HOME.

export XDG_DATA_HOME=$HOME
 
PORT=37247  #The port on which code-server will run
 
# The length of the password

LENGTH=7
 
# Generate a numeric random password

PASSWORD=$(cat /dev/urandom | tr -dc '0-9' | head -c "$LENGTH")
 
echo "Your random password is: $PASSWORD"
 
# Is important to create a specific folder on which the VSCODE configuration file will be saved. Because the VSCODE inside Posit Workbench save the config file in the folder $HOME/.config/code-server folder. So to avoid overwrite the VSCODE Posit Workbench configurations is important to create this secondary folder named as code-server_cnode

DIR="$HOME/.config/code-server_cnode/"
 
#Create the configuration directory

if [ ! -d "$DIR" ]; then

  #If the directory doesn't exist create it

  mkdir -p "$DIR"

  echo "The folder '$DIR' is created."

else

  echo "The folder '$DIR' already exists."

fi
 
# The name of the code server configuration file

CONFIG_FILE="$HOME/.config/code-server_cnode/config.yaml"
 
#Write the configurations inside the CONFIG_FILE

echo "bind-addr: 0.0.0.0:$PORT  # Specify the PORT on which code-server will run" > $CONFIG_FILE

echo "auth: password             # Enable the password authentication" >> $CONFIG_FILE

echo "password: $PASSWORD   # Substitute with your password" >> $CONFIG_FILE

echo "cert: false                # Disable the certificate SSL/TLS" >> $CONFIG_FILE

echo "disable-telemetry: true    # Disable telemetry" >> $CONFIG_FILE
 
echo "Configuration file was created successfully: $CONFIG_FILE"
 
 
# Print some useful information

echo "********************************************************************"

echo "Starting code-server in Slurm"

echo "Environment information:"

echo "Date:" $(date)

echo "Allocated node:" $(hostname)

echo "Listening on:" $PORT

echo "Password to access VSCode:" $PASSWORD
 
echo "********************************************************************"
 
#Launch code server on the compute node http://localhost:37247  ssh -NL 37247:cnode42:37247 manuel.tardaguila@hpclogin.fht.org

PASSWORD=$PASSWORD code-server --config $CONFIG_FILE
```
 
