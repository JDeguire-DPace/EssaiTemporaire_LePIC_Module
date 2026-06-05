module mod_intro
        use iso_fortran_env
        use mpi
    implicit none
    private

    public :: print_introduction, print_license, draw_lepic

    contains

    subroutine draw_lepic()
        implicit none
        write(*,*) '                                      '
        write(*,'(a)') '                                                                       ⣶⣄⠀⠀⢀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀'
        write(*,'(a)') ' █████                ███████████  █████   █████████              ⠀⠀⠀⠀⢸⣿⣿⣷⣴⣿⡄⠀⠀⠀⠀⠀⢀⡀⠀⠀⠀'
        write(*,'(a)') '░░███                ░░███░░░░░███░░███   ███░░░░░███     ███    ⠀⠀⠀⠀⠰⣶⣾⣿⣿⣿⣿⣿⡇⠀⢠⣷⣤⣶⣿⡇⠀⠀⠀⠀⠀'
        write(*,'(a)') ' ░███         ██████  ░███    ░███ ░███  ███     ░░░     ░███     ⠀⠀  ⠙⣿⣿⣿⣿⣿⣿⣿⣀⣿⣿⣿⣿⣿⣧⣀⠀⠀⠀'
        write(*,'(a)') ' ░███        ███░░███ ░██████████  ░███ ░███          ███████████  ⣷⣦⣀⠘⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠃⠀'
        write(*,'(a)') ' ░███       ░███████  ░███░░░░░░   ░███ ░███         ░░░░░███░░░ ⢲⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠁⠀⠀'
        write(*,'(a)') ' ░███      █░███░░░   ░███         ░███ ░░███     ███    ░███   ⠀ ⠙⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡟⠁⠀'
        write(*,'(a)') ' ███████████░░██████  █████        █████ ░░█████████     ░░░    ⠀⠀⠚⠻⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⠿⠿⠂⠀⠀'
        write(*,'(a)') '░░░░░░░░░░░  ░░░░░░  ░░░░░        ░░░░░   ░░░░░░░░░           ⠀⠀⠀⠀⠀⠀⠉⠙⢻⣿⣿⡿⠛⠉⡇⠀⠀⠀⠀⠀⠀⠀'
        write(*,'(a)') '                                                                        ⠋⠁⠀⠀⠸⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀'
        write(*,'(a)') '                                Laplace Explicit Particle-In-Cell code   '
        write(*,'(a)') '                                        (Modular Object-Oriented Edition)'
        write(*,'(a)') '                                     '

        return
    end subroutine draw_lepic

    subroutine print_license()
            implicit none
            character(len=*), parameter :: fleur = '+'
  
            write(*,'(A)') '  ' // repeat(fleur//' ',37)
            write(*,'(A)') '  ' //fleur//'                                                                            '//fleur
            write(*,'(A)') '  ' //fleur//'     3D explicit parallel particle in cell code                             '//fleur
            write(*,'(A)') '  ' //fleur//'                                                                            '//fleur
            write(*,'(A)') '  ' //fleur//'   LePIC 3D algorithm version: 3.3.6.2                                      '//fleur
            write(*,'(A)') '  ' //fleur//'   Last modifications: Sep/2024                                             '//fleur
            write(*,'(A)') '  ' //fleur//'                                                                            '//fleur
            write(*,'(A)') '  ' //fleur//'   Copyright or © or Copr. Gwenael Fubiani (2022/11/22)                     '//fleur
            write(*,'(A)') '  ' //fleur//'                                                                            '//fleur
            write(*,'(A)') '  ' //fleur//'   gwenael.fubiani@cnrs.fr                                                  '//fleur
            write(*,'(A)') '  ' //fleur//'                                                                            '//fleur
            write(*,'(A)') '  ' //fleur//'   This software is a computer program whose purpose is to model low        '//fleur
            write(*,'(A)') '  ' //fleur//'   temperature plasmas with a Particle-In-Cell algorihm in 3D-3V            '//fleur
            write(*,'(A)') '  ' //fleur//'   dimensions.                                                              '//fleur
            write(*,'(A)') '  ' //fleur//'                                                                            '//fleur
            write(*,'(A)') '  ' //fleur//'   This software is governed by the CeCILL-C license under French law and   '//fleur
            write(*,'(A)') '  ' //fleur//'   abiding by the rules of distribution of free software.  You can  use,    '//fleur
            write(*,'(A)') '  ' //fleur//'   modify and/ or redistribute the software under the terms of the CeCILL-C '//fleur
            write(*,'(A)') '  ' //fleur//'   license as circulated by CEA, CNRS and INRIA at the following URL        '//fleur
            write(*,'(A)') '  ' //fleur//'   "http://www.cecill.info".                                                '//fleur
            write(*,'(A)') '  ' //fleur//'                                                                            '//fleur
            write(*,'(A)') '  ' //fleur//'   As a counterpart to the access to the source code and  rights to copy,   '//fleur
            write(*,'(A)') '  ' //fleur//'   modify and redistribute granted by the license, users are provided only  '//fleur
            write(*,'(A)') '  ' //fleur//'   with a limited warranty and the software''s author, the holder of the     '//fleur
            write(*,'(A)') '  ' //fleur//'   economic rights,  and the successive licensors  have only  limited       '//fleur
            write(*,'(A)') '  ' //fleur//'   liability.                                                               '//fleur
            write(*,'(A)') '  ' //fleur//'                                                                            '//fleur
            write(*,'(A)') '  ' //fleur//'   In this respect, the user''s attention is drawn to the risks associated   '//fleur
            write(*,'(A)') '  ' //fleur//'   with loading,  using,  modifying and/or developing or reproducing the    '//fleur
            write(*,'(A)') '  ' //fleur//'   software by the user in light of its specific status of free software,   '//fleur
            write(*,'(A)') '  ' //fleur//'   that may mean  that it is complicated to manipulate,  and  that  also    '//fleur
            write(*,'(A)') '  ' //fleur//'   therefore means  that it is reserved for developers  and  experienced    '//fleur
            write(*,'(A)') '  ' //fleur//'   professionals having in-depth computer knowledge. Users are therefore    '//fleur
            write(*,'(A)') '  ' //fleur//'   encouraged to load and test the software''s suitability as regards their  '//fleur
            write(*,'(A)') '  ' //fleur//'   requirements in conditions enabling the security of their systems and/or '//fleur
            write(*,'(A)') '  ' //fleur//'   data to be ensured and,  more generally, to use and operate it in the    '//fleur
            write(*,'(A)') '  ' //fleur//'   same conditions as regards security.                                     '//fleur
            write(*,'(A)') '  ' //fleur//'                                                                            '//fleur
            write(*,'(A)') '  ' //fleur//'   The fact that you are presently reading this means that you have had     '//fleur
            write(*,'(A)') '  ' //fleur//'   knowledge of the CeCILL-C license and that you accept its terms.         '//fleur
            write(*,'(A)') '  ' //fleur//'                                                                            '//fleur
            write(*,'(A)') '  ' //' ' // repeat(fleur//' ',38)
            write(*,'(A)') '                                                  '

            return
    end subroutine print_license

    subroutine print_introduction()
        implicit none
        call draw_lepic()
        call print_license()
        return
    end subroutine print_introduction


end module mod_intro