# -*- tcl -*-

# NAMD run script template. Contact: giacomo.fiorin@gmail.com
# Version: 2022-04-22


if { [catch numPes] } {
    # Wrap NAMD's print command
    proc print { args } {
        puts [join ${args}]
    }
}

set mol_name popc_400
if { [info exists env(mol_name)] > 0 } {
    set mol_name $env(mol_name)
}

if { [file exists ${mol_name}.psf] == 0 } {
    print "Error: The file \"${mol_name}.psf\" is missing."
    exit 1
}


# Job control

# Define run_index based on myReplica
if { [catch numReplicas] == 0 } {
    set n_replicas [numReplicas]
    if { ${n_replicas} > 1 } {
        set env(run_index) [expr 100 + [myReplica]]
    }
}

# Extract the run label from the name of this script
set run [file rootname [file tail [info script]]]

# Handle different replicas for ensemble runs
set run_path_prefix ""
if { [info exists env(run_index)] > 0 } {
    # Set replica-specific run label
    set run "${run}-$env(run_index)"
    # Assume that each replica runs in a separate sub-folder
    set run_path_prefix "${run}/"
} else {
    if { (${run} == "us") } {
        print "Error: Undefined run_index variable for umbrella-sampling jobs"
        exit 1
    }
}

# Prefix run label with molecular system's label
set run ${mol_name}.${run}

set old ""
set job ${run}
set ijob 0


if { (${run} != "${mol_name}.min") && (${run} != "${mol_name}.test") } {
    if { [info exists env(old)] > 0 } {
        set old $env(old)
    } else {
        # Auto-generate this job's label
        set coor_files [list]
        catch {
            set coor_files [glob -type f "${run_path_prefix}${run}.\[0-9\]\[0-9\]\[0-9\]\[0-9\]\[0-9\]\[0-9\].coor"]
        } err
        if { [llength ${coor_files}] > 0 } {
            set coor_files [lsort ${coor_files}]
            set last_coor_file [lindex ${coor_files} end]
            # Get the job index from the last .coor file's name
            set last_ijob [string trimleft \
		[file extension [file rootname ${last_coor_file}]] ".0"]
            if { [file extension [file rootname ${last_coor_file}]] ==
                 ".000000" } {
                # Special case (when everything is trimmed out)
                set last_ijob 0
            }
            set ijob [expr ${last_ijob} + 1]
            set old [format "%s.%06d" ${run} ${last_ijob}]
        }
        set job [format "%s.%06d" ${run} ${ijob}]
    }
}


if { [catch numPes] } {
    # Not running NAMD: print the name of this job and exit
    print ${job}
    exit
}


# Detect force field files

if { [file exists charmm] > 0 } {
    set ff "CHARMM"
    set ff_folder "charmm"
}

if { [file exists martini] > 0 } {
    set ff "MARTINI"
    set ff_folder "martini"
}


# Set run parameters from environment variables

set ensemble "NPT"; # Choices: NVT, NPT, NPAT, NPgammaT

set ff "CHARMM"

set pbc "yes"
set pbc_aniso_xy "yes"

set cutoff 12.0
set vdwForceSwitching on
set fullElectFrequency 3

set temperature 300.0
set pressure 1.0
set pressure 1.0
set surface_tension 0.0
set langevinDamping -1.0
set langevinPistonPeriod 200.0
set langevinPistonDecay 100.0


proc update_defaults {} {

    global ff
    global ff_folder

    global timestep
    global langevinDamping

    global log_freq
    global dcd_freq
    global restart_freq

    if { [info exists ff] == 0 } {
        # Fail on first call of this function when ff is undefined
        print "ERROR: Force field undefined."
        exit 1
    }

    if { [info exists ff_folder] == 0 } {
        if { ${ff} == "CHARMM" } {
            set ff_folder "charmm"
        }
        if { ${ff} == "MARTINI" } {
            set ff_folder "martini"
        }
    }

    if { [info exists timestep] == 0 } {
        if { ${ff} == "CHARMM" } {
            set timestep 2.0
        }
        if { ${ff} == "MARTINI" } {
            set timestep 25.0
        }
    }

    if { ${langevinDamping} < 0.0 } {
        if { ${ff} == "CHARMM" } {
            set langevinDamping 1.0
        }
        if { ${ff} == "MARTINI" } {
            set langevinDamping 1.0
        }
    }

    if { [info exists log_freq] == 0 } {
        if { ${ff} == "CHARMM" } {
            set log_freq 600
        }
        if { ${ff} == "MARTINI" } {
            set log_freq 2000
        }
    }

    if { [info exists dcd_freq] == 0 } {
        set dcd_freq 2500
        if { ${ff} == "MARTINI" } {
            set dcd_freq 2000
        }
    }

    if { ([info exists restart_freq] == 0) && [info exists dcd_freq] == 1 } {
        set restart_freq [expr ${dcd_freq} * 100]
    }
}


foreach keyword { ff ff_folder ensemble timestep cutoff \
                      fullElectFrequency vdwForceSwitching \
                      pbc pbc_aniso_xy temperature pressure surface_tension \
                      langevinDamping langevinPistonPeriod langevinPistonDecay \
                      log_freq dcd_freq restart_freq \
                  } {
    if { [info exists env(${keyword})] > 0 } {
        puts "Setting ${keyword} = \"$env(${keyword})\" from environment variable (default = [set ${keyword}])."
        set ${keyword} $env(${keyword})
        update_defaults
    }
}

update_defaults

# Input

# Physical description
paraTypeCharmm          on
structure               ${mol_name}.psf
proc parameters_safe { param_file } {
    # Allow including optional files, doesn't fail when absent
    if { [file exists ${param_file}] > 0 } {
        print "Found parameters file ${param_file}."
        parameters      ${param_file}
    }
}

if { ${ff} == "CHARMM" } {

    parameters_safe     ${ff_folder}/par_all36_lipid.prm
    parameters_safe     ${ff_folder}/par_all36_na.prm
    parameters_safe     ${ff_folder}/par_all36_carb.prm
    parameters_safe     ${ff_folder}/par_all36_cgenff.prm
    parameters_safe     ${ff_folder}/par_all36_cgenff_${mol_name}.prm
    parameters_safe     ${ff_folder}/par_all36m_prot.prm
    parameters_safe     ${ff_folder}/par_water_ions.prm

    exclude             scaled1-4
    1-4scaling          1.0
    rigidBonds          all
    useSettle           on
    vdwForceSwitching   ${vdwForceSwitching}
}

if { ${ff} == "MARTINI" } {

    parameters_safe     ${ff_folder}/martini_v2.2.namd.prm

    cosAngles           on
    exclude             1-2
    1-4scaling          1.0
    switching           on
    martiniSwitching    on
    dielectric          15.0
}

# Cutoffs and approximations
switchdist              [expr ${cutoff} - 2.0]
cutoff                  ${cutoff}
pairlistdist            [expr ${cutoff} + 1.5]
if { ${ff} == "CHARMM" } {
    nonbondedFreq       1
    fullElectFrequency  3
    stepspercycle       24
}
if { ${ff} == "MARTINI" } {
    nonbondedFreq       1
    fullElectFrequency  1
    stepspercycle       10
}
margin                  2.0
if { ${pbc} == "yes" && ${ff} != "MARTINI" } {
    PME                 yes
    #PMEGridSizeX        72
    #PMEGridSizeY        72
    #PMEGridSizeZ        72
    PMEGridSpacing      0.9
    FFTWUseWisdom       yes
    FFTWWisdomFile      ${run_path_prefix}.fftw_wisdom.txt
}


# Initial data and boundary conditions
coordinates             ${run_path_prefix}${mol_name}.pdb
if { ${ijob} > 0 } {
    binCoordinates      ${run_path_prefix}${old}.coor
    binVelocities       ${run_path_prefix}${old}.vel
    if { ${pbc} == "yes" } {
    extendedSystem      ${run_path_prefix}${old}.xsc
    }
} elseif { (${run} != "${mol_name}.min") } {

    if { [info exists env(run_index)] > 0 } {
        binCoordinates ${run_path_prefix}${run}.initial.coor
    } else {
        if { [file exists ${mol_name}.eq.coor] } {
            binCoordinates  ${mol_name}.eq.coor
        } else {
            binCoordinates  ${mol_name}.min.coor
        }
    }
    seed                87654321
    temperature         ${temperature}
    if { ${pbc} == "yes" } {
        if { [info exists env(run_index)] > 0 } {
            extendedSystem ${run_path_prefix}${run}.initial.xsc
        } else {
            if { [file exists ${mol_name}.eq.xsc] } {
                extendedSystem      ${mol_name}.eq.xsc
            } elseif { [file exists ${mol_name}.min.xsc] } {
                extendedSystem      ${mol_name}.min.xsc
            } else {
                extendedSystem      ${mol_name}.xsc
            }
        }
    }
} else {
    temperature 0.0
    if { ${pbc} == "yes" } {
        extendedSystem ${mol_name}.xsc
    }
}

COMmotion no
zeroMomentum yes
if { ${pbc} == "yes" } {
    if { ${ff} == "CHARMM" } {
        wrapAll yes
    } else {
        wrapAll yes
    }
}

if { (${run} != "${mol_name}.min") } {

    # Thermodynamic ensemble options

    if { ${ensemble} != "NVE" } {

        # Per-atom Langevin temperature coupling
        langevin                on
        langevinTemp            ${temperature}
        langevinDamping         ${langevinDamping}

    }

    if { (${pbc} == "yes") && (${ensemble} != "NVE") && \
             (${ensemble} != "NVT") } {

        # Langevin piston pressure coupling

        if { ${pbc_aniso_xy} == "yes" } {
            if { (${ensemble} == "NPAT") } {
                useConstantArea yes
            } else {
                useFlexibleCell     yes
                useConstantRatio    yes
            }

            if { (${ensemble} == "NPgammaT") } {
                surfaceTensionTarget ${surface_tension}
            }
        }

        useGroupPressure        on
        LangevinPiston          on
        LangevinPistonTarget    ${pressure}
        LangevinPistonPeriod    ${langevinPistonPeriod}
        LangevinPistonDecay     ${langevinPistonDecay}
        LangevinPistonTemp      ${temperature}
    }
}


# Simulation type

if { (${run} != "${mol_name}.min") } {
    timestep ${timestep}
}

if { [info exists env(numsteps)] > 0 } {
    set numsteps        $env(numsteps)
} else {
    if { (${run} != "${mol_name}.min") } {
        # 1 ns
        set numsteps    [expr 25000 * [stepspercycle]]
    } else {
        set numsteps    [expr 100 * [stepspercycle]]
    }
}


if { [file exists ${run_path_prefix}${mol_name}.posre.pdb] > 0 } {
    # Position restraints
    constraints         on
    consRef             ${run_path_prefix}${mol_name}.posre.pdb
    consKFile           ${run_path_prefix}${mol_name}.posre.pdb
    consKCol            O
}


# Output

# Main output and restarts
outputName              ${run_path_prefix}${job}
binaryOutput            yes
if { ${restart_freq} > 0 } {
    if { (${run} == "${mol_name}.min") } {
        restartName     ${run_path_prefix}${job}
        restartFreq     100
    } else {
        restartName     ${run_path_prefix}${job}.rs
        restartFreq     ${restart_freq}
    }
    binaryRestart       yes
}

# Trajectory output
if { (${run} != "${mol_name}.min") } {
    DCDfile             ${run_path_prefix}${job}.coor.dcd
    DCDfreq             ${dcd_freq}
    if { (${pbc} == "yes") } {
        DCDUnitCell     yes
        XSTfile         ${run_path_prefix}${job}.xst
        XSTfreq         ${log_freq}
    }
}

# Standard output frequencies
outputEnergies          ${log_freq}
outputMomenta           ${log_freq}
outputPressure          ${log_freq}
outputTiming            ${log_freq}


# Colvars configuration
if { [file exists ${run_path_prefix}colvars.tcl] > 0 } {
    source ${run_path_prefix}colvars.tcl
}


# Run for the required steps
if { (${run} != "${mol_name}.min") } {

    if { ${n_replicas} == 1 } {
        run ${numsteps}
    } else {
        # Run multiple-replicas scripts
        if { [file exists replicas.tcl] > 0 } {
            source replicas.tcl
        } else {
            run ${numsteps}
        }
    }

} else {

    minimize ${numsteps}

    foreach ext { "coor" "vel" "xsc" } {
        # Cleanup unneeded backup files
        foreach backup_ext { "old" "BAK" } {
            file delete ${job}.${ext}.${backup_ext}
        }
    }
    file delete ${job}.vel
}
