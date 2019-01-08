.. _inputs:

===========
Input
===========

Command Line
____________

The program is invoked on the command line:

python -m transire.py -i example.in -o output.out

  or

transire.py -i example.in -o output.out

-i  path to file with input keywords (required).

-o  name of output file (optional).  Can be specified in the input file.

Input File
__________

The input file is a plain text document comprised of keywords and values separated by an '='.
Lines without keywords or that start with '#' are ignored.
The input file is specified using the '-i' flag on the command line.

**Formatting:**

The ideal formatting for a keyword line is:

keyword = value

with a space on either side of the '=' and nothing after the value.  However, the spaces are not
generally necessary but may result in an error if omitted.

Keywords with the (array) mark in the description are given as several integers separated by spaces
such as:

keyword = 1 0 0

General inputs
______________

**crys_a_file**
            the complete or relative path to the source file for the first crystal.
            the file format can be .cif, extended .xyz, or any format supported by
            ASE that includes the lattice dimensions of the unit cell.

**crys_a_surface**
            (array) miller indices for the initial surface of the first crystal.

**crys_a_layers**
            number of layers in the z direction of the first crystal (default = 3).

**crys_b_file**
            the complete or relative path to the source file for the second crystal.
            the file format can be .cif, extended .xyz, or any format supported by
            ASE that includes the lattice dimensions of the unit cell.

**crys_b_surface**
            (array) miller indices for the initial surface of the second crystal.

**crys_b_layers**
            number of layers in the z direction of the first crystal (default = 3).

**separation**
            distance to be added between the two crystals when forming the interface.

**max_atoms**
            the maximum number of atoms that can present in the interface structure.  If
            the number of atoms is greater than max_atoms, an error message is printed 
            and the interface structure is rejected. (default = 75000)

**project_name**
            name used when creatimg output files. (default = "default_name")

**tolerance**
            error criteria when replicating unit cells to ensure periodicity and 
            and replicating supercells to match up interface.  A smaller value leads to
            less error and stress, but can result in much larger interface structures.
            (default = 5.0e-2)

**calculate_initial_energy**
            switch to enable calculating the energy at the initial interface structure.
            (default = True)

**ignore_initial_structure**
            disables the forced exit if the initial structure generation is unsuccessful.
            useful if the initial structure has a very large number of atoms.  If set to
            true and an error occurs, then "calculate_initial_energy" is set to False.
            (default = False)

**output_file**
            name of the file for printing results.  Can be specified on the command line
            with the '-o' option, though the value in the input file supersedes if it is present.
            (default = transire.out)

**remove_duplicates**
            checks constructed supercells to ensure that no atoms with identical coordinates
            have been generated.  Generally not needed and can be memory intensive.
            (default = False)

**energy_method**
            specifies which program to use for calculating the energy.  Current options are:
            lammps and cp2k (default = cp2k)

**print_debug**
            print additional output from the interface structure generation. (default = False)

**calculate_binding_energy**
            Each energy calculation will separately calculate 

**full_periodicity**
            logical switch to turn on periodicity along all three axes instead of just
            parallel to the interface. (default = False)

**z_axis_vacuum**
            adds a length of vacuum along the z-axis, split evenly between the top and
            bottom of the interface structure. (default = 0.0)


Search parameters
_________________

**search_list**
            (array) a list of searches to perform in sequential order.  The searches
            can be placed in any order and repeated. The lowest interface energy configuration 
            is passed to the next search.  Availiable searches:

                  0) Markov chain search
                  1) Twist angle search
                  2) Molecule insertion
                  3) Separation optimization

            example:  search_list = 1 0 1 2

**surface_search**
            switch to turn on performing sequence of searches in 'search_list'
            repeatedly for a range of surfaces given by 'range_surface_(h/k/l)_(a/b)'
            If False, then the single set of indices given with 'crys_(a/b)_surface' are used.
            (default = False)

**range_surface_(h/k/l)_(a/b)**
            (array) a list of values to be used in miller indices
            for generating surfaces used in the search.  Two or more
            values must be provided for each keyword.  To specify a
            single value, duplicate the value.  example:

|                  range_surface_h_a = 0 1
|                  range_surface_k_a = 0 1 2
|                  range_surface_l_a = 1 2
|                  range_surface_h_b = 0 0
|                  range_surface_k_b = 2 1
|                  range_surface_l_b = 0 1 2 3

Markov Chain parameters
_______________________

**number_of_steps**
            number of combined translations and rotations to be performed
            (default = 25)

**markov_type**
            type of movement allowed in constrained search (default = 2).

            0) Rotation only
            1) Translation only
            2) Rotation and Translation

**mc_translate_x**
            x coordinate for restart option (default = 0.0)

**mc_translate_y**
            y coordinate for restart option (default = 0.0)

**mc_rotate**
            angle for restart configuration in degrees (default = 0.0)

**mc_restart**
            logical for using a previous MC step configuration at the start of the MC search 
            (default = False)

Twist Angle parameters
______________________

**angles_list**
            (array) a list of angles in degrees to use when generating interface
            structures.  example:

            angles_to_gen = 15 27 85

**angles_stepsize**
            number of degrees between each interface structure when using 
            number_of_angles

**number_of_angles**
            number of angles to include in twist angle search when not using
            angles_to_gen.  example:

            angles_to_iter = 1
            number_of_angles = 60

**starting_angle**
            angle to start with when iterating over angles (default = 0.0)

**angle_write_energy_file**
            switch to enable printing a log file with the angles used in
            the search and the associated interface energies. (default = True)

**angle_write_coord_file**
            switch to enable printing the coordinates of all interface structures
            generated in the twist angle search in the xyz format.  (default = True)

**angle_calculate_energy**
            switch to enable calculating the interface energy of each interface
            structure generated. (default = True)

**angle_return_initial**
            switch to disable returning the lowest energy configuration from the
            twist angle search.  This is useful when performing multiple twist angle
            searches or multiple surfaces. (default = False)

**angle_optimize_separation**
            switch to enable separation optimizer after each rotation.
            (default = True)

**ras_depth**
            number of layers to perform in Reducing Angle Search that reduces the number of
            energy calculations by searching a range of rotations to find the resulting 
            interface structures that have the fewest atoms.  ras_depth = 1 is equivalent
            to a normal Twist Angle search.  Each subsequent layer involves an angle
            search around the angles that result in the smallest interface structures found
            in the previous layer using a stepsize one order of magnitude smaller than the
            previous layer.  (default = 1)

**ras_factor**
            number of angles to be used as the starting points for each layer of the ras search.
            (default = 5)

**ras_energy**
            switch to enable calculating the energy of the interfaces produced by the final
            layer of the ras search (default = True)

**ras_all_angles**
            switch to accept all of the interfaces successfully generated during the first layer.
            Overwrites 'ras_factor'. (default = False)

**read_in_ras_file**
            path to file with previous ras results that are read in to populate angles_list.

Molecule Insertion parameters
_____________________________

.. important::
            performing another search after inserting the molecules will remove the inserted
            molecules.  Insertion should be done either as the last step in the search or after
            each search.

**insert_file**
            relative or absolute path to coordinate file of atoms to be inserted.
            A file must be provided for the insertion or the program will quit with an error.
            The inserted atoms are placed between the surfaces such that the centers of the
            surfaces and the inserted atoms all align.  The surfaces are separated to make
            room for the inserted atoms.

**calc_insert_energy**
            calculate the energy of the interface after the atoms are inserted.
            If false, an input file is written but no calculation takes place.
            (default = False)

**insert_vacuum_below**
            value to add below the inserted molecule to allow for control of alignment.
            (default = 0.5)

**insert_vacuum_above**
            value to add above the inserted molecule to allow for control of alignment.
            (default = 0.0)

**rotate_insert_x**
            value in degrees to rotate the inserted molecule by around the x-axis.

**rotate_insert_y**
            value in degrees to rotate the inserted molecule by around the y-axis.

**rotate_insert_z**
            value in degrees to rotate the inserted molecule by around the z-axis.

.. note::
            multiple rotations occur in order of x -> y -> z

Separation Optimization parameters
__________________________________

**sep_guess**
            the initial guess for the optimization (default = 0.5).
            The conclusion of the optimization replaces this value with the optimal separation.

**sep_max_steps**
            max number of steps to be carried out in the optimization before returning the
            the last guess. (default = 25)

**sep_tolerance**
            convergence condition for the distance between two steps (default = 1e-5)

**sep_intial_step**
            initial change in separation for each step. (default = 0.1)

CP2K parameters
_______________

**cp2k_input**
            the complete or relative path to the file with pycp2k commands to
            generate the cp2k input file.  the same parameter can be given
            using '-c' in the command line.  The value in the input file overrides
            the command line if both are given.

**max_mpi_processes**
            the max number of processes to be passed to mpirun or similar program
            as set during the installation of pycp2k (default = 32)

**atoms_per_process**
            changes the number of processes for cp2k calculation based on the number
            of atoms in the calculation.  max_mpi_processes is used if any of these
            are true: (default = 13)

            1) atoms_per_process = 0 (default)
            2) number of atoms/atoms_per_process > max_mpi_processes
            3) number of atoms < atoms_per_process

**working_directory**
            path to directory where the cp2k input and output will be
            generated (default = "./")

LAMMPS parameters
_________________

**lammps_input**
            input file for specifying LAMMPS commands.  The format is a single list declaration
            of the form:

            commands=[
            'entry 1',
            'entry 2',
            ...
            'entry X'
            ]

            The first and last lines in the example must be present.  LAMMPS keywords related
            to defining the system are handled by ASE.  For all available keywords see:
            http://lammps.sandia.gov/doc/Section_commands.html

ET parameters
_____________

**perform_ET**
            switch to enable the electron transport calculation. (default = 'False')

**et_method**
            selects program to use for performing NEGF ET calculation.
            options include ASE and slymer. (default = 'ASE')

**number_of_layers_a**
            number of layers of crystal a that is assigned to the left lead in the
            Green's Function method.  Each layer is one unit cell thick as defined
            by the initially read in coordinate file. If using the ASE ET calculator,
            the number of layers must be a multiple of 2. (default = 2)

**number_of_layers_b**
            number of layers of crystal b that is assigned to the right lead in the
            Green's Function method.  Each layer is one unit cell thick as defined
            by the initially read in coordinate file. If using the ASE ET calculator,
            the number of layers must be a multiple of 2. (default = 2)

**ET_restart**
            switch to disable all calculations before the ET calculation step.  If set to
            'True', the output from a previous run is expected to be provided using
            "restart_path".  (default = 'False')

**restart_path**
            direct path to folder containing the ".out" files from a previous ET calculation.
            (default = './')

**restart_file**
            file name for restart file in "restart_path".  Only include the main part of the name
            and not the suffix (eg .traj)

**exclude_coupling**
            for use with debugging.  The coupling between the layers in the leads are used
            in place of the coupling between the leads and the scatter region. (default = False)

**energy_levels_ET**
            (array) specify the range of energy values in eV relative to the Fermi level.
            The first and second numbers are the lower and upper limits respectively.  The
            third number is the step size between each energy.  The default behavior is
            to only calculate transmission at the Fermi level.

**orthonormal_overlap**
            switch to replace the overlap matrix with the identity matrix. (default = False)
