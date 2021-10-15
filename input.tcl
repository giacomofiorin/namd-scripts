# -*- tcl -*-

# NAMD run script template. Contact: giacomo.fiorin@gmail.com


set mol_name ppta
if { [info exists env(mol_name)] > 0 } {
    set mol_name $env(mol_name)
}

if { [file exists ${mol_name}.psf] == 0 } {
    if { [catch numPes] } {
        puts stderr "Error: The file \"${mol_name}.psf\" is missing."
    } else {
        print "Error: The file \"${mol_name}.psf\" is missing."
    }
    exit 1
}


# Job control

# Extract the run label from the name of this script
set run [file rootname [file tail [info script]]]
if { [info exists env(run_index)] > 0 } {
    # For uncoupled ensembles, use this variable to identify the copy
    set run "${run}-$env(run_index)"
} else {
    if { (${run} == "us") } {
        # Umbrella sampling requires the flag to identify the window
        print "Missing run_index variable"
        exit 1
    }
}

# Prefix run label with molecular system's label
set run ${mol_name}.${run}

set old ""
set job ${run}
set ijob 0


if { (${run} != "${mol_name}.min") } {
    if { [info exists env(old)] > 0 } {
        set old $env(old)
    } else {
        # Auto-generate job label
        set coor_files [list]
        catch {
            set coor_files [glob -type f \
                "${run}.\[0-9\]\[0-9\]\[0-9\]\[0-9\]\[0-9\]\[0-9\].coor"]
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


if { [catch numReplicas] == 0 } {
    set n_replicas [numReplicas]
    if { ${n_replicas} > 1 } {
        # Append replica label if multiple replicas are in use
        set job [format "${job}.rep%03d" [myReplica]]
    }
}

if { [catch numPes] } {
    # Not running NAMD: print the name of this job and exit
    puts ${job}
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

set temperature 300.0
set pressure 1.0
set pressure 1.0
set surface_tension 0.0
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

    if { [info exists langevinDamping] == 0 } {
        if { ${ff} == "CHARMM" } {
            set langevinDamping 10.0
        }
        if { ${ff} == "MARTINI" } {
            set langevinDamping 1.0
        }
    }

    if { [info exists log_freq] == 0 } {
        set log_freq 500
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
    vdwForceSwitching   on
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
    stepspercycle       20
}
if { ${ff} == "MARTINI" } {
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
    FFTWWisdomFile      .fftw_wisdom.txt
}


# Initial data and boundary conditions
coordinates             ${mol_name}.pdb
if { ${ijob} > 0 } {
    binCoordinates      ${old}.coor
    binVelocities       ${old}.vel
    if { ${pbc} == "yes" } {
    extendedSystem      ${old}.xsc
    }
} elseif { (${run} != "${mol_name}.min") } {

    if { [info exists env(run_index)] > 0 } {
        binCoordinates ${run}.initial.coor
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
            extendedSystem ${run}.initial.xsc
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
    temperature         0.0
    if { ${pbc} == "yes" } {
    extendedSystem      ${mol_name}.xsc
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

    if { (${ensemble} == "NVT") || (${ensemble} == "NPT") } {

        # Per-atom Langevin temperature coupling
        langevin                on
        langevinTemp            ${temperature}
        langevinDamping         ${langevinDamping}

    }

    if { (${pbc} == "yes") && (${ensemble} != "NVT") } {

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
        set numsteps    500000
    } else {
        set numsteps    2000
    }
}


if { [file exists ${mol_name}.posre.pdb] > 0 } {
    # Position restraints
    constraints         on
    consRef             ${mol_name}.pdb
    consKFile           ${mol_name}.posre.pdb
    consKCol            O
}


# Output

# Main output and restarts
outputName              ${job}
binaryOutput            yes
if { ${restart_freq} > 0 } {
    if { (${run} == "${mol_name}.min") } {
        restartName     ${job}
        restartFreq     100
    } else {
        restartName     ${job}.rs
        restartFreq     ${restart_freq}
    }
    binaryRestart       yes
}

# Trajectory output
if { (${run} != "${mol_name}.min") } {
    DCDfile             ${job}.coor.dcd
    DCDfreq             ${dcd_freq}
    if { (${pbc} == "yes") } {
        DCDUnitCell     yes
        XSTfile         ${job}.xst
        XSTfreq         ${log_freq}
    }
}

# Standard output frequencies
outputEnergies          ${log_freq}
outputMomenta           ${log_freq}
outputPressure          ${log_freq}
outputTiming            ${log_freq}


# Colvars configuration
if { [file exists ${run}.colvars.tcl] > 0 } {
    source ${run}.colvars.tcl
} else {
    if { [file exists colvars.tcl] > 0 } {
        source colvars.tcl
    }
}


# Run for the required steps
if { (${run} != "${mol_name}.min") } {

    if { ${n_replicas} == 1 } {
        run ${numsteps}
    } else {
        # Run multiple-replicas scripts
        if { [file exists ${run}.replicas.tcl] > 0 } {
            print "RUN ${run}.replicas.tcl"
            source ${run}.replicas.tcl
        } else {
            if { [file exists replicas.tcl] > 0 } {
                source replicas.tcl
            }
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
