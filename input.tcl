# -*- tcl -*-

set mol_name ppta
if { [info exists env(mol_name)] > 0 } {
    set mol_name $env(mol_name)
}

if { [file exists ${mol_name}.psf] == 0 } {
    print "The file \"${mol_name}.psf\" is missing."
    exit 1
}

set ff "CHARMM"
if { [file exists martini] > 0 } {
    set ff "MARTINI"
}

if { (${ff} == "CHARMM") && ([file exists charmm] == 0) } {
    print "The \"charmm\" directory is missing."
    exit 1
}

# Choices: NVT, NPT, NPAT, NPgammaT
set ensemble "NPT"
if { [info exists env(ensemble)] > 0 } {
    set ensemble $env(ensemble)
}

set pbc "yes"
set pbc_aniso_xy "yes"

set temperature 300.0
set pressure 1.0
set pressure 1.0
set surface_tension 0.0

foreach keyword { pbc pbc_aniso_xy temperature pressure surface_tension } {
    if { [info exists env(${keyword})] > 0 } {
        set ${keyword} $env(${keyword})
        puts "Setting ${keyword} equal to $env(${keyword})"
    }
}

set timestep 2.0
if { ${ff} == "MARTINI" } {
    set timestep 25.0
}
if { [info exists env(timestep)] > 0 } {
    set timestep $env(timestep)
}

set cutoff 12.0
if { [info exists env(cutoff)] > 0 } {
    set cutoff $env(cutoff)
}

set log_freq 500
if { [info exists env(log_freq)] > 0 } {
    set log_freq $env(log_freq)
}

set dcd_freq 2500
if { [info exists env(dcd_freq)] > 0 } {
    set dcd_freq $env(dcd_freq)
}

set restart_freq [expr ${dcd_freq} * 10]
if { [info exists env(restart_freq)] > 0 } {
    set restart_freq $env(restart_freq)
}


# Job control
set run                 [file rootname [file tail [info script]]]
if { [info exists env(run_index)] > 0 } {
    set run "${run}-$env(run_index)"
} else {
    if { (${run} == "us") } {
        print "Missing run_index variable"
        exit 1
    }
}
set run                 ${mol_name}.${run}
set old                 ""
set job                 ${run}
set ijob                0

if { (${run} != "${mol_name}.min") } {
    if { [info exists env(old)] > 0 } {
        set old $env(old)
    } else {
        set coor_files [list]
        catch {
            set coor_files [glob -type f \
                "${run}.\[0-9\]\[0-9\]\[0-9\]\[0-9\]\[0-9\]\[0-9\].coor"]
        } err
        if { [llength ${coor_files}] > 0 } {
            set coor_files [lsort ${coor_files}]
            set last_coor_file [lindex ${coor_files} end]
            set last_ijob [string trimleft \
		[file extension [file rootname ${last_coor_file}]] ".0"]
            if { [file extension [file rootname ${last_coor_file}]] ==
                 ".000000" } {
                set last_ijob 0
            }
            set ijob [expr ${last_ijob} + 1]
            set old [format "%s.%06d" ${run} ${last_ijob}]
        }
        set job [format "%s.%06d" ${run} ${ijob}]
    }
    if { [info exists env(job)] > 0 } {
        set job $env(job)
    } 
}

if { [lindex [file split [info nameofexecutable]] end] == "tclsh" } {
    # Print the name of this job
    puts ${job}
    exit
}


# Input

# Physical description 
paraTypeCharmm          on
structure               ${mol_name}.psf
proc parameters_safe { param_file } {
    if { [file exists ${param_file}] > 0 } {
        print "Found parameters file ${param_file}."
        parameters      ${param_file}
    }
} 

if { ${ff} == "CHARMM" } {

    parameters_safe     charmm/par_all36_prot.prm
    parameters_safe     charmm/par_all36_lipid.prm
    parameters_safe     charmm/par_all36_na.prm
    parameters_safe     charmm/par_all36_carb.prm
    parameters_safe     charmm/par_all36_cgenff.prm
    parameters_safe     charmm/par_all36_cgenff_${mol_name}.prm
    parameters_safe     charmm/par_water_ions.prm

    exclude             scaled1-4
    1-4scaling          1.0
    rigidBonds          all
    useSettle           on
}

if { ${ff} == "MARTINI" } {

    parameters_safe     martini/martini_v2.2.namd.prm

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
        langevinDamping         10.0            

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
        LangevinPistonPeriod    200.0
        LangevinPistonDecay     100.0
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
if { [file exists colvars.tcl] > 0 } {
    source colvars.tcl
}
if { [file exists ${run}.colvars.tcl] > 0 } {
    source ${run}.colvars.tcl
}


# run
if { (${run} != "${mol_name}.min") } {
    run                ${numsteps}
} else {
    minimize           ${numsteps}
    # Cleanup unneeded files
    foreach ext { "coor" "vel" "xsc" } {
        foreach backup_ext { "old" "BAK" } {
            file delete ${job}.${ext}.${backup_ext}
        }
    }
    file delete ${job}.vel
}
