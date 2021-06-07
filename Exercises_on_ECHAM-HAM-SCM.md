# How to run ECHAM-HAM in the Single Column Model (SCM) Mode on MISTRAL

## Preparation 

* Learn vim
  * <https://www.openvim.com> 
	* <https://vim-adventures.com>
	* <http://www.vimgenius.com/lessons/vim-intro/levels/level-1>


* Test (or prepare) your machine for remote connections with ssh and sshfs
  * <https://www.digitalocean.com/community/tutorials/how-to-use-sshfs-to-mount-remote-file-systems-over-ssh>


* Test (or install) panoply
  * <https://www.giss.nasa.gov/tools/panoply/download/>


* Mistral is the supercomputer at the German Climate Computer Center DKRZ. Check infos on: <https://www.dkrz.de/up/systems/mistral>


## EXERCISE 1: First ECHAM Run on Mistral

### Login into Mistral

* What you need
  * Your account number ("MyUser")
  * Your password ("MyPASSWORD")

* Log into mistral using your user and password (platform dependent)
  * On Mac
    
      ```	
      ssh MyUser@mistral.dkrz.de
      ```

#### Now on MISTRAL

* Go to your home directory

    ```
    cd ~
    ```
    
### Prepare your environment
    
* Make a variable including your project number (bb1224 for T2 course in SS2021)

    ``` 
      export ACCOUNT=bb1224
    ``` 

* Load the required modules for ECHAM-HAM

      ``` 
        module use /pf/zmaw/m222045/m222045_modulefiles
        module load jobscript_toolkit
        module load autotools_for_echam/6.3.0
        module load intel/16.0
        module load mxm/3.4.3082
        module load bullxmpi_mlx/bullxmpi_mlx-1.2.9.2
      ```

* And set the correct language
      ```
        export LC_ALL=C
      ```

* Add the lines in your ~/.bash_profile for your next session:

    ``` 
      vi ~/.bash_profile
    ```
  * On vim add
 	
  
      ``` 
        export ACCOUNT=bb1224
        module use /pf/zmaw/m222045/m222045_modulefiles
        module load jobscript_toolkit
        module load autotools_for_echam/6.3.0
        module load intel/16.0
        module load mxm/3.4.3082
        module load bullxmpi_mlx/bullxmpi_mlx-1.2.9.2
        export LC_ALL=C
      ```
    * In short: inside vim: 
      * press 'i' for insertmode, then paste the lines 
      * press 'esc' for normal mode, 
      * type ':wq' and press 'enter' to save and quit

### Download ECHAM

* Copy ECHAM-HAM to your home directory: about 5 min

    ```
        cp -av /work/bb1224/2021_MS-COURSE/model/MyClimateModel ~/MyClimateModel
    ```
    
### Download the Patch for Pratical Work with the Single Column Model
* Download patch, configuration and boundary conditions

    ``` 
      git clone https://github.com/diegovillanueva/SCM_workshop.git
    ```

* What is on the files?
  * 'SCM_workshop/scm-patch.patch'
  * 'SCM_workshop/ioColumn.patch'
    * Patches for the source code (MyClimateModel/src/)
  * 'SCM_workshop/settings_patch_scm.generic.patch'
    * A patch for the simulation settings (MyClimateModel/my_experiments/MyRun/)
  * 'SCM_workshop/stdatm_form_diss.nc'
    * A file with boundary conditions for the SCM model.

* Patch the code for SCM
	* Go to model directory
	
      ```	
        cd ~/MyClimateModel
      ```

  * Apply source code patch

      ```
      git apply --whitespace=nowarn ~/SCM_workshop/scm-patch.patch
      patch src/mo_iocolumn.f90 ~/SCM_workshop/ioColumn.patch
      ```
    
### Compile the ECHAM Model for Fortran
* Compile and configure the model: about 30 min (make a break and grep a coffee)
	* Configure the model for this machine (only needed the first time)
		
      ```
      ./configure --with-fortran=intel
      ```

	* Compile (repeat after code changes)
	
      ```
      make -j 40    # compile files in src/
      make install  # prepare binary file (bin/echam6)
      ```

### Running the ECHAM Single Column Model
* Make a run
	* Prepare run
	
	    ``` 
        cd my_experiments
        prepare_run.sh MyRun
	    ```
	* Configure run for SCM
	
	    ```
        cd MyRun
        patch settings_MyRun ../../../SCM_workshop/settings_patch_scm.generic.patch
	    ```
     
    * Create a directory for the output and run
    
      ```
      mkdir /work/$ACCOUNT/userspace/$USER
      jobsubm_echam.sh settings_MyRun
      ```

* Check if is runing

      ```
      squeue -u $USER
      ```

* Check the content of your results 

    ```
    cd /work/$ACCOUNT/userspace/$USER/MyRun
    cdo infon -vertmean -timmean -selname,xl,relhum MyRun_200701.01_echam.nc
    ```

### Perform Data Analysis and Plotting - On your Local Computer (open a new terminal)

* Open Terminal
  * Mount mistral disk (platform dependent) 
    * On Mac

      ```
          MyUser=MyUser #replace right side with your user
          MyPass=MyPass #replace right side with your pass
          sshfs $MyUser@mistral.dkrz.de:/work/bb1224/$MyUser ./MISTRAL -ovolname=NAME -o password_stdin <<< "$MyPass"
      ```

* Plot your results (platform dependent)

  * On Mac (panoply)
  
      ```
        open MISTRAL/MyRun/MyRun_200701.01_echam.nc
      ```


## EXERCISE 2: Rerunning ECHAM on MISTRAL with Parameter Changes

* Change a parameter in your simulation and rerun

	* Go to your experiment folder

      ```
    	cd ~/MyClimateModel/my_experiments
      ```

	* Prepare run

      ```
      prepare_run.sh MyRunWithChanges
      ```

	* Configure run for SCM

    	```
      cd MyRunWithChanges
      patch settings_MyRunWithChanges ../../../SCM_workshop/settings_patch_scm.generic.patch
   		```

	* Change a parameter (see ../../include/physctl.inc)

      ```
    	vi settings_MyRunWithChanges
      ```
    * e.g., Change 'nauto=2 to nauto=1'  in line 167 (This changes the autoconversion scheme)
      * TIP: press "/", type 'nauto', press 'enter' and press 'n' repeatedly to jump  
      * type '167' and then 'gg' to jump

	* Run

	    ```
	    jobsubm_echam.sh settings_MyRunWithChanges
	    ```

### Perform Data Analysis and Plotting on your Local Computer


* Plot your results (platform dependent)

  * On Mac (panoply)
  
      ```     
      open MISTRAL/MyRunWithChanges/MyRunWithChanges_200701.01_echam.nc
      ```

* Plot differences (e.g., differences in xl):

	* On Mac (panoply): drag and drop variable into plot 

## EXERCISE 3: Source Code Changes in ECHAM on MISTRAL
* Change a the source code in the model and recompile
  * Go to source directory

      ```
      cd ~/MyClimateModel/src
      ```
      
    * Change a file: e.g., mo_cloud_micro_2m.f90 (line 2996)
    
      ```
      vi mo_cloud_micro_2m.f90
      ```
      [change 'ccraut' to '100.0_dp*ccraut' in line 2834 (increase autoconversion efficiency by a hundredfold)]

* Go back, back up previous binary and compile

    ```
    cd ~/MyClimateModel
    mv bin/echam6 bin/echam6.original
    make -j 40    # compile files in src/
    make install  # prepare binary file (bin/echam6)
    ```

* Rerun and plot differences

	* Go to experiments

		```
		cd ~/MyClimateModel/my_experiments
		```
	* Prepare run

		```
		prepare_run.sh MyRunWithChangesInCode
		```
	* Configure run for SCM

    ```
		cd MyRunWithChangesInCode
		patch settings_MyRunWithChangesInCode ../../../SCM_workshop/settings_patch_scm.generic.patch
    ```

  * Run ECHAM
  
    ```
    jobsubm_echam.sh settings_MyRunWithChangesInCode
    ```


### Perform Data Analysis and Plotting  On your Local Computer

* Plot your results (platform dependent)

	* On Mac (panoply)
  
    ```     
    open MISTRAL/MyRunWithChangesInCode/MyRunWithChangesInCode_200701.01_echam.nc    
    ```

* Plot differences (e.g., differences in xl)

	* On Mac (panoply): drag and drop variable into plot


