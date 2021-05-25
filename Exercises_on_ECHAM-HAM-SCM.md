# How to run ECHAM-HAM in the Single Column Model (SCM) Mode on MISTRAL

## Preparation 

* Learn vim
  * <https://www.openvim.com> 
	* <https://vim-adventures.com>
	* <http://www.vimgenius.com/lessons/vim-intro/levels/level-1>
	* In short: 'inside vim, press 'i' for insertmode, then paste, press 'esc' for normal mode, type ':wq' and press 'enter' to save and quit)'


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
    
* Make a variable including your project number

    ``` 
      export ACCOUNT=bb1224
    ``` 
* Add the following lines in your ~/.bash_profile:

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
      ```

### Download ECHAM

* Copy ECHAM-HAM to your home directory: about 5 min

    ```
        rsync -a --progress /work/bb1224/b380602/echam6.3-ham2.3-moz1.0/ ~/MyClimateModel
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
     
    * Run
    
      ```
      jobsubm_echam.sh settings_MyRun
      ```

* Check the content of your results 

    ```
    cd /work/$ACCOUNT/$USER/MyRun
    cdo infon -vertmean -timmean -selname,xl,relhum MyRun_200701.01_echam.nc
    ```

### Perform Data Analysis and Plotting
#### On your Local Computer

* Mount mistral disk (platform dependent) 
	* On Mac

      ```
          sshfs MyUser@mistral.dkrz.de:/work/bb1224/MyUser ./MISTRAL -ovolname=NAME -o password_stdin <<< "MyPASSWORD"
      ```

* Plot your results (platform dependent)
	* Open a new terminal, and go to were you mounted the MISTRAL disk
	* On Mac (panoply)
			open MISTRAL/MyRun/MyRun_200701.01_echam.nc


## EXERCISE 2: Rerunning ECHAM on MISTRAL with Parameter Changes

* Change a parameter in your simulation and rerun

	* Go to experiments

      ```
    	cd /pf/b/$USER/MyClimateModel/my_experiments
      ```

	* Prepare run

      ```
    	cd my_experiments
		prepare_run.sh MyRunWithChanges
      ```

	* Configure run for SCM

    	```
      cd MyRunWithChanges
      patch settings_MyRunWithChanges ../../../SCM_workshop/settings_patch_scm.generic.patch
   		```

	* Change a parameter, e.g.  (see ../../include/physctl.inc)

      ```
    	vi settings_MyRunWithChanges
      ```
    	[Change 'nauto=2 to nauto=1'  (This changes the autoconversion scheme)]

	* Run

	    ```
	    jobsubm_echam.sh settings_MyRunWithChanges
	    ```

## Perform Data Analysis and Plotting

#### On your Local Computer

* Plot your results (platform dependent)
	* Open a new terminal, and go to were you mounted the MISTRAL disk
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
      cd /pf/b/$USER/MyClimateModel/src
      ```
      
    * Change a file: e.g., mo_cloud_micro_2m.f90 (line 2996)
    
      ```
      vi mo_cloud_micro_2m.f90
      ```
      [change 'ccraut' to '100.0_dp*ccraut' (increase autoconversion efficiency by a hundredfold)]

* Go back, back up previous binary and compile

    ```
    cd ..
    mv bin/echam6 bin/echam6.original
    make -j 40    # compile files in src/
    make install  # prepare binary file (bin/echam6)
    ```

* Rerun and plot differences

	* Go to experiments

		```
		cd /pf/b/$USER/MyClimateModel/my_experiments
		```
	* Prepare run

		```
		cd my_experiments
		prepare_run.sh MyRunWithChangesInCode
		```
	* Configure run for SCM

    ```
		cd MyRunWithChangesInCode
		patch settings_MyRunWithChangesInCode ../../../SCM_workshop/settings_patch_scm.generic.patch
    ```
  * Change necessary parameters, e.g.  (see ../../include/physctl.inc)

    ```
      vi settings_MyRunWithChangesInCode
    ```
	[in vim change 'nauto=2 to nauto=1']

  * Run ECHAM
  
    ```
    jobsubm_echam.sh settings_MyRunWithChangesInCode
    ```

### Perform Data Analysis and Plotting
#### On your Local Computer

* Plot your results (platform dependent)
	* Open a new terminal, and go to were you mounted the MISTRAL disk
	* On Mac (panoply)
    ```     
    open MISTRAL/MyRunWithChangesInCode/MyRunWithChangesInCode_200701.01_echam.nc    
    ```

* Plot differences (e.g., differences in xl)
	* On Mac (panoply): drag and drop variable into plot



