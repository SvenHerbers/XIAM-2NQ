# XIAM-2NQ
XIAM-2NQ is a spectral fitting code that handles up to two quadrupolar nuclei and 4 internal rotors. It is an extension of the original XIAM code by Hartwig, available at the PROSPE (http://info.ifpan.edu.pl/~kisiel/prospe.htm) website.
The repository contains all source files which can be compiled using the make file to create a linux applications. Alternatively, the files can be compiled in steps using the ifortran compiler by intel.
XIAMi2NQ.exe is a compiled executable using the Intel(R) Fortran Intel(R) 64 Compiler Classic for applications running on IA-32 (Version 2021.13.1).
Compilation was carried out using the following command line entries: <br/> 
ifort -c iamint.f <br/> 
ifort -c iamio.f <br/> 
ifort -c mgetx.f <br/> 
ifort -c iamv.f <br/> 
ifort -c iamlib.f <br/> 
ifort -c iamfit.f <br/> 
ifort -c iamv2.f <br/> 
ifort -c iamadj.f <br/> 
ifort -c iamm.f <br/> 
ifort -c iamsys.f <br/> 
ifort -c iam.f <br/> 

ifort  -static -o XIAMi2NQ.exe iam.f iamsys.f iamm.f iamadj.f iamv2.f iamfit.f iamlib.f iamv.f mgetx.f iamio.f iamint.f <br/> 

Citation: J. Chem. Phys. 162, 234304 (2025) https://doi.org/10.1063/5.0267651

----------------------------------------------------------------------------
## Parameter Table
In this table, `Π = (Pα - ρ Pz)` refers to the relative internal angular momentum in the ρ-Axis System (RAS).

The operator `Pz` is used with a double meaning, referring either to the z axis in the RAS or in the Principal Axis System (PAS).  
In expressions where both occur, the RAS `Pz` is identifiable by being enclosed between two Wigner small-d matrices, `dt` and `d`, indicating the rotation from the RAS to the PAS.


### Semi-Rigid Rotor Parameters Hrr
Here only A reduction (S reduction also available)  
Operators are formulated in the PAS

| `Parameter` | `Operator` |
|------------|------------|
| `BJ` | `P**2` |
| `BK` | `Pz**2` |
| `B-` | `Px**2-Py**2` |
| `DJ` | `-P**4` |
| `DJK` | `-Pz**2P**2` |
| `DK` | `-Pz**4` |
| `dj` | `-2P**2 (Px**2-Py**2)` |
| `dk` | `-[Pz**2(Px**2 - Py**2)+(Px**2 - Py**2)Pz**2]` |
| `H_J` | `P**6` |
| `HJK` | `P**4 Pz**2` |
| `HKJ` | `P**2 Pz**4` |
| `H_K` | `Pz**6` |
| `h_j` | `2P**4 (Px**2 - Py**2)` |
| `hjk` | `P**2 [Pz**2 (Px**2 - Py**2) + (Px**2 - Py**2) Pz**2]` |
| `h_k` | `[Pz**4 (Px**2 - Py**2) + (Px**2 - Py**2) Pz**4]` |
| `L_J` | `P**8` |
| `LJK` | `P**4 Pz**4` |
| `LJJK` | `P**6 Pz**2` |
| `LKKJ` | `P**2 Pz**6` |
| `L_K` | `Pz**8` |
| `l_j`  | `2*P**6 (Px**2 - Py**2)` |
| `ljk` | `P**4  [Pz**2 (Px**2 - Py**2) + (Px**2 - Py**2) Pz**2]` |
| `lkj` | `P**2  [Pz**4 (Px**2 - Py**2) + (Px**2 - Py**2) Pz**4]` |
| `l_k` | `[Pz**6 (Px**2 - Py**2) + (Px**2 - Py**2) Pz**6]` |
|------------|------------|
| `Dzx` | `(PzPx + PxPz)`|
| `DzxJ` | `P**2 (PzPx + PxPz)`|
| `DzxK` | `(Pz**3Px + PxPz**3)`|



### Internal Rotor Parameters Hir
Operators are formulated in the RAS

| `Parameter` | `Operator` |
|------------|------------|
| `F` | `Π**2` |
| `rho` | `(ρ Pz**2 - 2 Pα Pz)` |
| `V1n` | `0.5(1-cos(nα))` |
| `V2n` | `0.5(1-cos(2nα))` |
| `Fmk` | `Pα ρ Pz Π**2` |
| `Fm2k2` | `Pα**2 ρ**2 Pz**2 Π**2` |
| `Fmk3` | `Pα ρ**3 Pz**3 Π**2` |
| `Fm3k` | `Pα**3 ρ Pz Π**2` |
| `Fm3k3` | `Pα**3 ρ**3 Pz**3 Π**2` |
|------------|------------|
| `mk` | `Pα ρ Pz` |
| `m2k2` | `Pα**2 ρ**2 Pz**2` |
| `mk3` | `Pα ρ**3 Pz**3` |
| `m3k` | `Pα**3 ρ Pz` |
|------------|------------|
| `DFm2` | `Pα**2 Π**2` |
| `Fk2` | `ρ Pz**2 Π**2` |

### Internal rotation - overall rotation distortion Hird
Mixed PAS/RAS formulation  
RAS operators are rotated into PAS, then multiplied with PAS operators  
Rotation is indicated by dt(O)d, with d representing the wigner small d matrix, and dt its transposed

| `Parameter` | `Operator` |
|------------|------------|
| `Dpi2J` | `2*dt(Π**2)d P**2` |
| `Dpi2K` | `[dt(Π**2)d Pz**2+ Pz**2 dt(Π**2)d]` |
| `Dpi2-` | `[dt(Π**2)d (Px**2 - Py**2)+ (Px**2 - Py**2) dt(Π**2)d]` |
| `Dp2JJ` | `2*dt(Π**2)d P**4` |
| `Dp2KK` | `[dt(Π**2)d Pz**4 + Pz**4 dt(Π**2)d]` |
| `Dp2JK` | `P**2[dt(Π**2)d Pz**2 + Pz**2 dt(Π**2)d]` |
| `Dp2-j` | `-P**2[dt(Π**2)d (Px**2 - Py**2) + (Px**2 - Py**2) dt(Π**2)d]` |
| `Dp2-k` | `-[dt(Π**2)d Pz**2(Px**2 -P_y**2)`<br>` + dt(Π**2)d (Px**2 - Py**2) Pz**2`<br>` + Pz**2 dt(Π**2)d (Px**2 - Py**2)`<br>` + (Px**2 - Py**2) dt(Π**2)d Pz**2`<br>` + Pz**2 (Px**2 - Py**2)dt(Π**2)d`<br>` + (Px**2 - Py**2) Pz**2 dt(Π**2)d]` |
| `Dp2zx` | `-[dt(Π**2)d Pz Px`<br>` + dt(Π**2)d Px Pz`<br>` + Pz dt(Π**2)d Px`<br>` + Px dt(Π**2)d Pz`<br>` + Pz Px dt(Π**2)d`<br>` + Px Pz dt(Π**2)d]` |
|------------|------------|
| `Dpi4`  | `dt(Π**4)d` |
| `Dpi4J` | `2*dt(Π**4)d P**2` |
| `Dpi4K` | `[dt(Π**4)d Pz**2+ Pz**2 dt(Π**4)d]` |
| `Dpi4-` | `[dt(Π**4)d (Px**2 - Py**2)+ (Px**2 - Py**2) dt(Π**4)d]` |
| `Dp4JJ` | `2*dt(Π**4)d P**4` |
| `Dp4KK` | `[dt(Π**4)d Pz**4 + Pz**4 dt(Π**4)d]` |
| `Dp4JK` | `P**2[dt(Π**4)d Pz**2 + Pz**2 dt(Π**4)d]` |
| `Dp4-j` | `-P**2[dt(Π**4)d (Px**2 - Py**2) + (Px**2 - Py**2) dt(Π**4)d]` |
| `Dp4-k` | `-[dt(Π**4)d Pz**2(Px**2 -P_y**2)`<br>` + dt(Π**4)d (Px**2 - Py**2) Pz**2`<br>` + Pz**2 dt(Π**4)d (Px**2 - Py**2)`<br>` + (Px**2 - Py**2) dt(Π**4)d Pz**2`<br>` + Pz**2 (Px**2 - Py**2)dt(Π**4)d`<br>` + (Px**2 - Py**2) Pz**2 dt(Π**4)d]` |
| `Dp4zx` | `-[dt(Π**4)d Pz Px`<br>` + dt(Π**4)d Px Pz`<br>` + Pz dt(Π**4)d Px`<br>` + Px dt(Π**4)d Pz`<br>` + Pz Px dt(Π**4)d`<br>` + Px Pz dt(Π**4)d]` |
|------------|------------|
| `Dpi6`  | `dt(Π**6)d` |
| `Dpi6J` | `2*dt(Π**6)d P**2` |
| `Dpi6K` | `[dt(Π**6)d Pz**2+ Pz**2 dt(Π**6)d]` |
| `Dpi6-` | `[dt(Π**6)d (Px**2 - Py**2)+ (Px**2 - Py**2) dt(Π**6)d]` |
|------------|------------|
| `rhoJ`  | `2*dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d P**2` |
| `rhoK` | `[dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d Pz**2+ Pz**2 dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d]` |
| `rho-` | `[dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d (Px**2 - Py**2)+ (Px**2 - Py**2) dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d]` |
| `rhoJJ`  | `2*dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d P**4` |
| `rhoKK` | `[dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d Pz**4+ Pz**4 dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d]` |
| `rhoJK` | `P**2 [dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d Pz**2+ Pz**2 dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d]` |
| `rho-j` | `-P**2[dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d (Px**2 - Py**2) + (Px**2 - Py**2) dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d]` |
| `rho-k` | `-[dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d Pz**2(Px**2 -P_y**2)`<br>` + dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d (Px**2 - Py**2) Pz**2`<br>` + Pz**2 dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d (Px**2 - Py**2)`<br>` + (Px**2 - Py**2) dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d Pz**2`<br>` + Pz**2 (Px**2 - Py**2)dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d`<br>` + (Px**2 - Py**2) Pz**2 dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d]` |
| `rhozx` | `-[dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d Pz Px`<br>` + dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d Px Pz`<br>` + Pz dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d Px`<br>` + Px dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d Pz`<br>` + Pz Px dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d`<br>` + Px Pz dt(F*(2.0 * ρ * Pz**2 - 2.0 * Pα * Pz))d]` |
|------------|------------|
| `mkJ` | `2*dt(Pα ρ Pz)d P**2` |
| `mkK` | `[dt(Pα ρ Pz)d Pz**2+ Pz**2 dt(Pα ρ Pz)d]` |
| `mkD` | `[dt(Pα ρ Pz)d (Px**2 - Py**2)+ (Px**2 - Py**2) dt(Pα ρ Pz)d]` |
| `mk3J` | `2*dt(Pα ρ**3 Pz**3)d P**2` |
| `mk3K` | `[dt(Pα ρ**3 Pz**3)d Pz**2+ Pz**2 dt(Pαρ**3Pz**3)d]` |
| `mk3-` | `[dt(Pα ρ**3 Pz**3)d (Px**2 - Py**2)+ (Px**2 - Py**2) dt(Pαρ**3Pz**3)d]` |
| `m3kJ` | `2*dt(Pα**3 ρ Pz)d P**2` |
| `m3kK` | `[dt(Pα**3 ρ Pz)d Pz**2+ Pz**2 dt(Pα**3 ρ Pz)d]` |
| `m3k-` | `[dt(Pα**3 ρ Pz)d (Px**2 - Py**2)+ (Px**2 - Py**2) dt(Pα**3 ρ Pz)d]` |
|------------|------------|
| `Dc3J` | `2*dt(cos(nα))d P**2` |
| `Dc3K` | `[dt(cos(nα))d Pz**2+ Pz**2 dt(cos(nα))d]` |
| `Dc3-` | `[dt(cos(nα))d (Px**2 - Py**2)+ (Px**2 - Py**2) dt(cos(nα))d]` |
| `D3JJ` | `2*dt(cos(nα))d P**4` |
| `D3KK` | `[dt(cos(nα))d Pz**4 + Pz**4 dt(cos(nα))d]` |
| `D3JK` | `P**2[dt(cos(nα))d Pz**2 + Pz**2 dt(cos(nα))d]` |
| `D3-j` | `-P**2[dt(cos(nα))d (Px**2 - Py**2) + (Px**2 - Py**2) dt(cos(nα))d]` |
| `D3-k` | `-[dt(cos(nα))d Pz**2(Px**2 -P_y**2)`<br>` + dt(cos(nα))d (Px**2 - Py**2) Pz**2`<br>` + Pz**2 dt(cos(nα))d (Px**2 - Py**2)`<br>` + (Px**2 - Py**2) dt(cos(nα))d Pz**2`<br>` + Pz**2 (Px**2 - Py**2)dt(cos(nα))d`<br>` + (Px**2 - Py**2) Pz**2 dt(cos(nα))d]` |
| `D3zx` | `-[dt(cos(nα))d Pz Px`<br>` + dt(cos(nα))d Px Pz`<br>` + Pz dt(cos(nα))d Px`<br>` + Px dt(cos(nα))d Pz`<br>` + Pz Px dt(cos(nα))d`<br>` + Px Pz dt(cos(nα))d]` |
### Top-Top coupling term Hii

| `Parameter` | `Operator` |
|------------|------------|
| `F12` | `dt1 Π1 d1 * dt2 Π2 d2 + dt2 Π2 d2 * dt1 Π1 d1` |
| `Vcc` | `dt1 cos(nα1) d1 * dt2 cos(nα2) d2` |
| `Vss` | `dt1 sin(nα1) d1 * dt2 sin(nα2) d2` |

### Other parameters

| `Parameter` | `Description` |
|------------|---------------|
| `F0` | `Structural constant to internal rotation. Can be used as input to derive ρ,β,γ`<br>`F = ℏ**2/(2rIα)= F0/r `<br>`With`<br>`r=1-Σg λg**2 Iα/Ig`<br>`and with λ being the direction cosines, and g=a,b,c in the PAS`|
| `delta` | `The angle between the internal rotation axis and the principal axis z` |
| `epsil` | `The angle between the principal axis x and the projection of the internal rotation axis onto xy-plane` |
| `beta` | `rotation about y axis` |
| `gamma` | `rotation about x axis` |
| `betJ1` `betJ2` `betJ3` `betJ4` | `Parameters to add P**2 dependence to beta.`<br/>`beta_total = beta + betJ1*(J(J+1))+ betJ2*(J(J+1))**2+ betJ3*(J(J+1))**3+ betJ4*(J(J+1))**4` |

### First order local parameters

These parameters multiply with Px,Py,Pz in the PAS, but can be defined for each Species separately (up to 11 species possible).  
They can be particularly helpful for 1st order corrections, when dealing with two low barrier coupling rotors, when the few coupling parameters available in XIAM (F12,Vcc,Vss) reach their limits.

| `Parameter` |
|-------------|
| `DxS1, DxS2,.. DxS11` |
| `DyS1, DyS2,.. DyS11` |
| `DzS1, DzS2,.. DzS11` |

### Nuclear quadrupole coupling parameters
for first or second nucleus

| `Parameter` |
|-------------|
| `chi_z, chi2_z` |
| `chi_-, chi2_-` |
| `chi_xy, chi2_xy` |
| `chi_xz, chi2_xz` |
| `chi_yz, chi2_yz` |

### Vibrational state interaction parameters
Parameter E is always available. 
For Wilson and Pickett coriolis coupling parameters, vibrational coupling mode must be activated by setting control parameter ctrl to  `ctrl 1`.
Can only be used with the “old” approximate nuclear quadrupole coupling (elements off diagonal in J neglected)  with a single quadrupolar nucleus.  


| `Parameter` | `Description` |
|------------|---------------|
| `E` | `Energy-offset assignable to the various vibrational states` |
| `Gx12,Gx34,Gy12,Gy34,Gz12,Gz34` | `Wilson type Coriolis coupling parameters for coupling between states 1&2 and 3&4. Should not be mixed with Pickett type at the moment, due to a likely phase inconsistency` |
| `Fxy12,Fxy34,Fyz12,Fyz34,Fxz12,Fxz34` | `Pickett type Coriolis coupling parameters for coupling between states 1&2 and 3&4. Should not be mixed with Wilson type at the moment, due to a likely phase inconsistency` |

## Update Notes
  Sven Herbers, 28-January-2026
  - Updated Example-Methylformate fit of the subset of v=0, Jmax=50, Kamax=20 lines. The unweighted rms of XIAM on this subset is 95 kHz; RAM36 global fits on the complete set of lines yield 71 kHz rms within this subset. 
  - Added a fit to Example-Methylformate treating the *complete* set (49 parameters, Jmax=62, Kamax=27, fmax= 668.1 GHz, fmin= 1.6 GHz, 6976 assignments) of v=0 lines from the v=0,1 dataset of lines provided in V. Ilyushin, et al. J. Mol. Spectrosc. 255, 32–38 (2009). [https://doi.org/10.1016/j.jms.2009.01.016](https://doi.org/10.1016/j.jms.2009.01.016). The unweighted rms of XIAM fits on this subset is 145 kHz; RAM36 global fits on the complete set of lines yield 75 kHz rms within this subset. 

  XIAM-2NQ v0.34 - Sven Herbers, 25-January-2026
  - Many new parameters available in Hird and Hir, parameter table will be updated in the following days.

  XIAM-2NQ v0.32 - Sven Herbers, 18-January-2026
  - Implementation of octic centrifugal distortion coefficients for Watson A and Watson S reduction (reduc 0 or reduc 1)
  - Many new parameters available in Hird and Hir
  - Updated Example-Methylformate treating the subset of v=0, Jmax=40, Kamax=15 lines. The unweighted rms of XIAM on this subset is 89 kHz; RAM36 global fits on the complete set of lines yield 66 kHz rms within this subset. 

  XIAM-2NQ v0.25 -> v0.29 - Sven Herbers, 05-January-2026
  - Many new parameters available in Hird and Hir, parameter table added to readme.md.
  - Minor code cleanup (removal/replacement of some functions)
  - Added Example-Methylformate treating the subset of v=0, Jmax=30, Kamax=10 lines. The unweighted rms of XIAM on this subset is 61 kHz; RAM36 global fits on the complete set of lines yield 54 kHz rms within this subset. 

  XIAM-2NQ v0.24b - Sven Herbers, 17-June-2025 

  - The XIAM-2NQ Publication is out! Updated go-to citation to J. Chem. Phys. 162, 234304 (2025) https://doi.org/10.1063/5.0267651
  - Removed some unnecessary comments from iam.fi

  XIAM-2NQ v0.24 - Sven Herbers, 31-Jan-2025 

  A small update to the iamint.f file:
  - An error was fixed that caused int 3 predictions to not work if no spin was present.
  - Modifications were made to allow for intensitiy predictions with J>100, if iam.fi 
    is modified accordingly and XIAM is recompiled. 

  XIAM-2NQ v0.23 -  Sven Herbers, 6-Jan-2025
  
  This is an updated version of XIAM-NQ with severale changes, the most important one
being the implementation of exact quadrupole coupling for two nuclei. 


----------------------------------------------------------------------------
## Further documentation
It follows a list of changes compared to original XIAM, the original documentation by Hartwig is printed below 
for completeness.<br/> 
New control parameters:  

     ctrl (default 3)    If spin is not 0 this switches exact quadrupole coupling on (3) or off (0) 
                           by discarding matrix elements off-diagonal in J.                       
     ctrl 0 : Switches to XIAM_mod - matrix elements off-diagonal in J are ignored. It allows for 
                analytical gradients for rigid rotor parameters and an increase in the torsional basis 
                if DIMVV > 1 is used pre-compilation.
     ctrl 1 : Sets up an experimental treatment of an uncoupled double well, without exact quadrupole 
                treatment, which is currently still being tested. Analytical gradients not implemented.          
     ctrl 2 : not used                      
     ctrl 3 : Exact quadrupole treatment, analytical gradients not implemented, increase in torsional 
                basis (e.g. "V 0 1  V 1 0") not implemented.
     qsum    (default 1000) partition function for SPCAT-type intensities.
     fsort   (default 3) This controls at which stage the F1 quantum number is assigned, relevant only
             for the two nuclei case. Case I: Separate by blocks of even and odd K as well as Wang signs
             first. Case II: use only separte blocks of even and odd K.
     fsort 1: Case I for A and E species
     fsort 2: Case II for A and E species
     fsort 3: Case I for A species and Case II for E species

     altering fsort might help if quantum number confusions are encounter (not uncommon in cases where 
     F1 is far from being a good quantum number)

New fit parameters:<br/> 
     Several new parameters, unrelated to the quadrupole treatment, have been implemented into the code. 
     These parameters are still being tested and will be described in a later version.
     
Intensity calculations <br/> 
     Previous XIAM versions did not implement intensity calculations for hyperfine components. An 
     approximation based on equation (15.142) from W. Gordy and R. L. Cook, Microwave Molecular Spectra, 
     3rd ed. (Wiley, New York, 1984), Vol. 18. has been implemented in XIAM-2NQ. 
     In this approximation, XIAM’s calculated non-hyperfine intensities are multiplied by 
     (2F1+1)x(2F1'+1)xsixj(J',F1',I1,F1,J,1)**2 x (2F+1)x(2F'+1)x sixj(F1',F',I2,F,F1,1)**2
     This approximations assumes that F1 and F are good quantum numbers. To loosen the restrictions
     on F1 a little, the F1 eigenvector components are used instead of the assigned F1 quantum number
     on ground and excited state. For each possible combination of F1 and F1' the expression is then  
     calculated and multiplied with the probability of the combination. <br/> <br/> 
     
     To evaluate the 6j symbols, a function from the wigner.f90 module, written by O. C. Gorton 
     (https://github.com/ogorton/, ogorton@sdsu.edu), was rewritten into subroutines and integrated 
     into XIAM's source code.
     This approximate treatment is expected to be accurate for single nucleus cases of 
     D, N, and Cl nuclei. For Br and I nuclei, it should still provide reliable results 
     for most transitions; however, transitions with delta J > 1 can not be predicted this way. 
     For two nuclei cases the predicted intensities will become worse if the quadrupole coupling
     of the first nucleus 1 is not about 10 times larger than the quadrupole coupling of nucleus 2
     chi_1 > 10* chi_2. In the equivalent case chi_1 = chi_2, the intensities are little reliable 
     and tend to predict more components than are present in the experimental spectrum.
     

     Some changes have been made to the intensity predictions, and new output columns have been added:
        log10_str**2: Base-10 logarithm of the reduced transition dipole moment squared.
        log10_total : Used for control parameter 'limit', this is the log10(str**2 *d_pop) with d_pop as 
                      population difference. 
        stat.w.     : M Degeneracy of lower state (2*J+1) or (2*F+1) in case of hyperfine components.
        d_pop.      : Population difference between the two involved states. Assumed are delta_M = 0 
                      transition rules, which means in the calculation of the population 
                      difference only the lower value of J or F is used to calculate the 'degeneracy'
                      of upper and lower involved levels.
        spcat       : Intensities as in Pickett’s SPCAT.

XIAM-2NQ uses 2I as input for control parameter 'spin', it also uses 2F as input for the quantum number 'F'.


XIAM_mod parameters: The additional empirical disortion parameters (Dc3K and Dc3-) available in 
                     XIAM_mod are also available in XIAM-NQ. 
Internal rotation overall - rotation distortion operator in XIAM_mod:
<pre>
  Hird =  2 Dpi2J (p_alpha - rho P_r)**2 P**2
         + Dpi2K [(p_alpha - rho P_r)**2 P_z**2 
                                  + P_z**2 (p_alpha - rho P_r)**2]
         + Dpi2- [(p_alpha - rho P_r)**2 (P_x**2 - P_y**2)
                                  + (P_x**2 - P_y**2) (p_alpha - rho P_r)**2]
         + 2 Dc3J  cos(3alpha) P**2
         + Dc3K  [cos(3alpha) P_z**2 + P_z**2 cos(3alpha)]       
         + Dc3-  [cos(3alpha) (P_x**2 - P_y**2) + (P_x**2 - P_y**2) cos(3alpha)]
</pre> 

-----------------------------------------------------------------------------
<pre>
  Program XIAM                   (Holger Hartwig, 15.November 1996)
                                 (email: hartwig@phc.uni-kiel.de  )
  Version 2.5e                   (email: phc25@rz.uni-kiel.d400.de)


  General Description
 
  XIAM can predict and fit the rotational spectrum of an asymmetric molecule
with maximal 3 symmetric internal rotors and one nucleus leading to a weak 
nuclear quadrupole coupling in the spectrum. Centrifugal distortion constants
of the overall rotation up to the 6th order and some 4th order distortion 
constants between the internal and overall rotation are included. To analyze
spectra of excited states of the internal rotation motion, some top-top 
coupling terms are used.


 Copyright 
 
 Everybody can use and copy XIAM for free. If you make changes in the source
code: (i)  mark them in the source itself 
      (ii) mention them in a README file
      (iii) give XIAM a new unique name (eg. xiam-xx)
      (iiii) send me an email (Thanks) 
 

 Method 
 
  XIAM uses the internal axes method given by Woods [1,2] and modified by 
Vacherand et. al. [3]. The centrifugal distortion constants of the 
Watson' A and S reduction [4] and the van Eijck-Typke reduction  are used. 
The nuclear quadrupole coupling is implemented with the matrix elements 
given in [5] but neglects the J off-diagonal elements. 
This leads to a similar restriction as the normal perturbation treatment 
concerning the magnitude of the nuclear quadrupole coupling. 
To analyze excited states different sets of rotational  
and kinetic internal rotation constants can be fitted simultaneously.
This program is described and used in [6] as well. Please cite [6]! 

  Hamiltonian

The hamiltonian matrix is set up in the principal axes system PAS of the whole
molecule. Therefore the centrifugal distortion constants are comparable
with rigid rotor calculations. The internal rotation operator Hi of each top
is set up in his own internal axes system (rho-system: RAS) and the resulting
eigenvalues of the diagonalisation of Hi are transformed (rotated) into
the principal axes system. This transformation is done via a rotation about 
two Euler angles beta and gamma, were beta and gamma are the angles between 
the rho system and the principal axes system. P_r is the angular momentum 
vector along the rho axis.

                        rho_x 
 gamma = arccos ---------------------------- rotation about x axis
                (rho_x**2 + rho_y**2)**0.5

                        rho_z 
 beta  = arccos ------------------------------------- rotation about y axis
                (rho_x**2 + rho_y**2 + rho_z**2)**0.5

  The top-top coupling matrix elements Hii are calculated by transforming the
operator (p_alpha - rho P_r) into the principal axes system and multiplying
the matrices then. 
  
 Rigid rotor part Hrr:

 Hrr = BJ P**2 + BK P_z**2 + B- (P_x**2 - P_y**2)


 Centrifugal distortion part Hcd (Watson A)

 Hcd = -DJ P**4 - DJK P_z**2  P**2 - DK P_z**4
       -2 dj P**2 (P_x**2 - P_y**2) 
       -dk (P_z**2 (P_x**2 - P_y**2) + (P_x**2 - P_y**2) P_z**2)
       +H_J P**6 + HJK P**4 P_z**2 + HKJ P**2 P_z**4 + H_K P_z**6
       +2 h_j P**4 (P_x**2 - P_y**2) 
       +hjk P**2 [P_z**2 (P_x**2 - P_y**2) + (P_x**2 - P_y**2) P_z**2)
       +h_k [P_z**4 (P_x**2 - P_y**2) + (P_x**2 - P_y**2) P_z**4]

 Internal rotation part of the i-th top in the rho system  
 
  Hi =   F (p_alpha - rho P_r)**2
       + 0.5 V1n (1 - cos( n alpha)) 
       + 0.5 V2n (1 - cos(2n alpha))
  where P_r is the angular momentum vector along the rho axis 

 Internal rotation distortion operator in the rho system Hid (none) 

 Internal rotation overall - rotation distortion operator in the PAS

  Hird =  2 Dpi2J (p_alpha - rho P_r)**2 P**2
         + Dpi2K [(p_alpha - rho P_r)**2 P_z**2 
                                  + P_z**2 (p_alpha - rho P_r)**2]
         + Dpi2- [(p_alpha - rho P_r)**2 (P_x**2 - P_y**2)
                                  + (P_x**2 - P_y**2) (p_alpha - rho P_r)**2]
         + Dc3J  cos(3alpha) P**2
                 

 Top-Top coupling term Hii

 Hii =   F12 ((p_alpha_1 - rho_1 P_r_1)(p_alpha_2 - rho_2 P_r_2) 
             +(p_alpha_2 - rho_2 P_r_2)(p_alpha_1 - rho_1 P_r_1))
       + Vss sin(n alpha_1) sin(n alpha_2)
       + Vcc cos(n alpha_1) cos(n alpha_2)


 The matrix elements of the nuclear quadrupole coupling are
 
   chi_z  e1 (3 K**2 - J(J+1))                             diagonal 
   chi-   e1 0.5((J(J+1)-K(K+1))(J(J+1)-(K+1)(K+2)))**0.5  2th off
   chi_xz e1 (1+2K)(J(J+1)-K(K+1))**0.5                    1th off
 i chi_yz e1 (1+2K)(J(J+1)-K(K+1))**0.5                    1th off
 i chi_xy e1 ((J(J+1)-K(K+1))(J(J+1)-(K+1)(K+2)))**0.5     2th off  
  
  e1= (0.75 G(G+1.0) - I(I+1)J(J+1))/(2I(2I-1)J(J+1)(2J-1)(2J+3))
  G = F(F+1) - I(I+1) - J(J+1)

  Spin Rotation parameters are
  C+  : C_x + C_y  
  C_z 
  C-  : C_x - C_y
 
 To fit different torsional states with different rotational constants,
 several indepentent sets of parameters can be used (maximum: DIMVB). 


  Input and Output

The input of XIAM is read from the standard input stream of the operating
system (UNIT=5), the output is written to the standard output stream (UNIT 6).
To use an input file on DOS, Unix, and OS/2 type

              xiam < inputfilename > outputfilename.

The input in a free format, it is analysed by the subroutine getx, which 
enables the user to write comments into the inputfile. They will be ignored
by XIAM. There are three types of comments:

 1) Pascal syntax:
           input_a {comment .. } input_b  
    is read  input_a input_b 

 2) Fortran 90 syntax:      
           input_a ! comment until end of line      
    is read  input_a
 
 3) Continuation character syntax:
           input_a \ comment until end of line
           input_b      
    is read  input_a input_b 

Comment 3) can be used to split one line over several lines in the 
input file. 
To get the '!', '{' or '\' characters without special meaning use 
'\!',  '\{' or '\\' instead. 

  The input of XIAM is divided into 7 blocks. Each block is separated 
from the other by an empty line. One block can not contain an empty line.
Each block has a special meaning:

  1. Block: Title and Comments. 
     This block will be written into the output file without any changes.
  
  2. Block: Control parameters.
     In this block basic numbers are given, e.g. the number of internal rotors,
     the nuclear spin I ...
  
  3. Block: Molecular parameters.
     Here the initial values of the parameters of the Hamiltonian are given.
  
  4. Block: Parameters to fit.
     The parameters (or linear combination of parameters) which are to be
     fitted are declared here.
  
  5. Block: Symmetry labeling.
     The symmetry species of the internal rotation (A, E) are defined here.
  
  6. Block: Torsional states.
     The torsional basis functions for the matrix are defined.
  
  7. Block: Transitions
     The quantum numbers and frequencies (if already known) are typed here.

Block 5 and 6 can be omitted if only a rigid rotor is calculated.

The input of each block has always the form keyword value [value [value..]] 

Some control parameters are sums of binary numbers, were each binary number
has a special meaning (e.g. woods and adj).

  List of keywords for each block:

1. Block: (no keyword)

2. Block:

  ncyc (or nzyk)
          maximum number of iteration cycles for the fit            C_NZYK   
  print   controls the output amount. (0)                          C_PRINT
  aprint  fine control of output amount (for debugging) 
  xprint  fine control of output amount (for debugging) 
  ints        intensity calculation                                     C_INTS 
          0: no intensities 
          1: intensities  of the transitions in the input-list
          2: output of all transitions whose frequencies are
             in the region 'freq_l'  to 'freq_h' and have intensities 
             more than 'limit'. 
             Output for low barrier molecules (sorted acc. symmetry label) 
          3: same as 2: but with a different sorting scheme.
             better suited for high barrier cases.
  reduc   type of reduction (0: Watson A, 1: Watson S, 2: van Eijck-Typke)
  freq_l  low boundary for int = 2 and 3 
  freq_h  upper boundary for int = 2 and 3 
  limit   intensity -limit for int = 2 and 3 
  temp    temperature for the Boltzmann distribution intensities 
  orger   use the given freqency errors to calculate the errors    C_ORGER
          of the parameters (1) or calculate weights of the 
          frequency errors and use obs.-calc. to calculate the 
          standard deviation of the parameters (0)  
  eval    create a file 'eval.out' containing the eigenvalues       C_EVAL 
          (0: no  1:yes)
  dfrq        create a file 'dfreq.out' containing the derivatives      C_DFRQ 
          of the frequencies (0: no  1:yes)
  maxm    number of basis functions for all tops. (8)
          2*maxm+1 functions will be used in the matrix 
          of operator Hi
  woods   controls the treatment of torsional integrals (abbr. TI) 
          <v_k sigma|v'_k' sigma'> for all sets.
          sum of binary values: 
          1:  calculate TI (required for the following terms)
          2:  normalize TI 
          4:  use TI in the transformation of one top (exact)
          8:  use TI in the transformation of all other tops
         32:  use TI in the rigid rotor part
         64:  use TI in the transformation of one top 
                                             (approximation Vacherand)
          if woods 0 is used TI's are completely ignored.
          normally good results can be obtained 
          with woods 33 (= 32 + 1)
  ndata       number of transitions to be read (normally not used)      C_NDATA
  nfold   potential barrier fold  (normally 3 for methyl tops)      C_NFOLD
  spin    twice nuclear spin of quadrupole coupling nucleus        C_SPIN 
          Nitrogen-14 : spin 2
  ntop    number of internal rotating tops                         C_NTOP 
  adjf    controls the value of F, F12, rho depending on the       C_ADJF 
          geometry (the values of rho and the angles)
          sum of binary values:
          1:  calculate F from rho, beta, and gamma each iteration
          2:  calculate F12 from rho, beta, and gamma each iteration
          4:  calculate F in a single top mode always (w/o F12)
          8:  calculate rho from F0, beta and gamma 
         16:  calculate rho, delta and epsilon from F0, beta and gamma 
             adfj is automatically calculated if not specified in the input.
  rofit   robust fitting control (0:off       )                 C_ROFIT
  defer       default frequency error                                 C_DEFER
  eps         tolerance limit in the SVD fit                          C_EPS  
  weigf       not used 
  convg       convergence limit (0.99 to 0.999999)                    C_CNVG 
  lambda  initial value of lambda in the Marquardt fit            C_LMBDA
  fitscl  scale the design matrix (1:no 0:yes)                    C_FITSC
  svderr  calculate the errors including all parameters  (1:no 0:yes)   
          (including parameters rejected by the SVD-fit)

4. Block

This block defines the initial values of the parameters in the Hamiltonian.
If different sets of rotational constants are used, the keywords contain
additional numbers in parenthesis, e.g. BJ(2) is the BJ value for the second
set of constants (noted as B 2 in the transition list).
If no number is given, the parameter is equal for all sets of constants. 

  BJ  *    0.5 (B_x + B_y)       B_x,y,z are the cartesian rotational constants
  BK  *    B_z - 0.5 (B_x + B_y) The representation must be chosen with the 
  B-  *    0.5 (B_x - B_y)       initial values of BJ, BK, B-.
  
  DJ  *    4th order centrifugal distortion  
  DJK *  
  DK  *  
  dj  *  
  dk or R6  *    choose R6 for van Eijck Typke's reduction. 
  
  H_J *    6th order centrifugal distortion  
  HJK *  
  HKJ *  
  H_K *  
  h_j *  
  hjk *  
  h_k *  
         
  chi_zz *   \chi zz     Quadrupole coupling constants
  chi-   *   \chi minus
  chi_xy *    
  chi_xz *  
  chi_yz *  
  
  C+  * = C_x + C_y      Spin Rotation coupling
  C_z *
  C-  * = C_x - C_y
  
  F12     Top Top interaction constants
  Vss    
  Vcc    
           
  V1n_1    V1n_2       internal rotation parameters for each top. a third top 
  V2n_1        V2n_2       can be calculated if the matrix dimensions are set  
  F_1          F_2         properly in the 'iam.fi' file before compilation.
  rho_1        rho_2  
  beta_1       beta_2 
  gamma_1      gamma_2
  Dpi2J_1      Dpi2J_2 
  Dpi2K_1      Dpi2K_2
  Dpi2-_1      Dpi2-_2
  Dc3J_1       Dc3J_2
  F0_1         F0_2       The inverse value of I_alpha in GHz (F0 = 505.379/I_alpha) 
  delta_1  delta_2    The angle between the internal rotation axis and the 
                      principal axis z (radiant: 0...2Pi).
  epsil_1  epsil_2    The angle between the principal axis x and the projection
                      of the internal rotation axis onto xy-plane.
               The delta and epsilon angles are the polar-coordinates 
               of the internal rotation axis in the principal axes system.
               lambda_z = cos(delta)
               lambda_x = sqrt((1-cos(delta)^2)) * cos(epsilon) 
               lambda_y = sqrt((1-cos(delta)^2)) * sin(epsilon) 
           cos(epsilon) = lambda_x * sign(lambda_y)/ SQRT(lambda_x^2 + lambda_y^2)   
                                                   ! SQRT correction: Mark Marshall 31.08.2021
  
  mu_x  Dipole moment components for intensity calculation.   
  mu_y  can not be fitted.
  mu_z

Note: often the value of I_alpha is predicted by the geometry (methyl top
      I_alpha = 3.18 - 3.2 uAA). To fix this value in the fit one must  
      set F0 to the appropriate value (about 159 - 157 GHz) and set 
      the 4th bit of the control parameter 'adjf' (decimal 8).
      The parameter rho can not be fitted in this case.

      To fit the angles between of internal rotation axis directly 
      (instead of indirectly with gamma and beta) the 5th bit of the 
      'adjf' parameter must be set (decimal 16). Gamma and beta can 
      not be fitted then. 
      

4. Block (fitting of parameters)

 In this block each line must begin with 'fit', 'dqu', dqx followed by the name
of the parameter to be fitted. If linear combinations of parameters are fitted,
the syntax is:
fit/dqu/dqx  parameter 1   coefficient 1   parameter 2   coefficient 2 ...
The number of lines in this block is the number of independent parameters
in the fit. 

fit   fit parameters using analytic derivatives
dqu   fit parameters using difference quotients to calculate the derivatives
dqx   fit parameters using two difference quotients to calculate the 
      derivatives 
parameter   is the keyword of the parameter (see 3th block)
coefficient is a number giving the coefficient for the linear combination.  

 The dqu/x keyword allows two optional numbers to occur before the first 
parameter. The first number gives the step width in percent for the 
differential quotient (default 0.1%), the second number is multiplied with
the variated parameter after the fit (Hartley convergence parameter).

 Not all parameters can be fitted with analytic derivatives, 
in the moment only the rigid rotor parameters have this feature.
They are denoted with a '*' in the table above. 
Maybe analytic derivatives of some other parameters will be available in 
the future, but the differential quotient method will always be necessary 
if (a) the parameters are affected by the 'adjf' control keyword and 
(b) woods 64 is used (Vacherand Method).

If only a prediction is wanted the 4 block can be replaced by an empty line.

Examples for fitting of linear combinations:
 
 dqu  V1n_1 1.0  V1n_2 1.0 
    The potential parameter V3 for two equivalent tops is fitted. 

 dqu  delta_1 1.0  delta_2  -1.0 
    Here the sum of the angles delta_1 and delta_2 is kept constant.
    This may be necessary if the molecule possesses a symmetry plane. 

5. Block (torsional symmetry)

S (or G for compatibility with old versions) followed by the sigma
  for each top.
Each line in this block defines the symmetry species of the torsional part.
Example for one three fold top 

   S  0        ! this is the A Species  (0 : sigma =0)
   S  1        ! this is the E Species  (1 : sigma =1)

Example for two different tops

  S  0 0  ! Top 1: sigma=0   Top 2: sigma=0  total symmetry AA
  S  0 1  ! Top 1: sigma=0   Top 2: sigma=1  total symmetry EiE
  S  1 0  ! Top 0: sigma=1   Top 2: sigma=1  total symmetry EEj
  S  1 1  ! Top 1: sigma=1   Top 2: sigma=1  total symmetry AE
  S -1 1  ! Top 1: sigma=-1  Top 2: sigma=1  total symmetry EA

A name of the symmetry species can be added at the beginning of a line,
the name must start with a slash:
 /A  S 0
 /E  S 1
The name is used in the output of the transitions.

6. Block 

V (only one keyword) followed by torsional state of each top.

 Each line in this block gives one set of rotational constants.
Example for one top in the ground state:

  V 0  ! matrix size is (2J+1) * (2J+1) 
 
 Example for two tops in the ground state:

  V 0 0 ! matrix size is (2J+1) * (2J+1) 
 
 Example for one top in the ground and first excited state with different 
sets of rotational constants for the states: 

  V 0   ! matrix size of ground state is (2J+1) * (2J+1) 
  V 1   ! matrix size of first excited state is (2J+1) * (2J+1) 
The constants of the second line are calculated from the constants of the
first line plus additional terms +BJ_2, +BK_2,  .., +gam_2.
In the transition list (block 7) B 1 refers to the ground state, B 2 to the 
excited state.

One line can include several V's. This will increase the size of the matrix
to include off diagonal matrix elements <v| H |v'>.
Example for one top in the ground and first excited state with the same
set of constants:  

  V 0  V 1  ! matrix size is [(2J+1)+(2J+1)] * [(2J+1)+(2J+1)]
            ! which includes ground and first excited state in one matrix 
In the transition list (block 7) V 1 refers to the ground state, V 2 to the 
excited state.

 Example for two tops in the first and second excited state, the
off diagonal elements are necessary for the top top interaction operators:

  V 0 1  V 1 0  ! matrix size is [(2J+1)+(2J+1)] * [(2J+1)+(2J+1)]


7. Block

Each line is one transition, the sequence of the keywords is arbitrary.
                   
Jup   Jlo     or  J    J quantum number

K-up  K-lo    or  K-   K minus and K plus (pseudo) quantum numbers
K+up  K+lo    or  K+   K- and K+ are only used to calculate t

=                      Frequency, if '=' not given, the predicted frequency 
                                  will be written, default GHz 

Err                    Frequency error

Sup   Slo     or  S    Torsional symmetry. the number refers to the n-th line
                       in block 5. a symmetry number of -1 indicates a rigid 
                       rotor transition

V1up  V1lo    or  V1   If more than one V keywords are given in one line of the
                       6th block, V 1 refers to the lowest, V 2 to the next
                       higher torsional level in the matrix, ....

V2up  V2lo             (not used)

Bup   Blo     or  B    the n-th set of rotational constants are used 
                       (refers to the n-th line in the 6th block). 

Fup   Flo     or  F    twice F, a F value of -1 indicates a transition without
                       nq-hfs

Tup   Tlo     or  T    tau number of the whole matrix ranging over all 
                       torsional states. tau numbers the eigenvalues in 
                       ascending order, starting with t 1
 
tup   tlo     or  t    tau number (range 1 - 2J+1) for one torsional state 

#                      splitting fit: the number (counted relatively from this 
                       transition) refers to the reference transition.
                       Example: # -1: the difference between this frequency
                                     and the previous frequency will be fitted 
                       If this transition is a prediction, the calc. splitting
                       is printed in the obs-calc output field. Furthermore,
                       if the reference transition was given, this frequency
                       is used as an offset to calculate the new frequency.
  
&                      average fit see '#'

Kup   Klo     or  K    K quantum number. This is a good quantum number for 
                       the E-states, especially at low barriers. 
                       For A-states the sign of K indicates the symmetry of
                       the eigenvector.

diff                   splitting input: the number refers to the reference
                       transition. The frequency following '= ' is the 
                       difference between the absolute frequency and the 
                       absolute frequency of the reference transition. 
                       if diff 0 is given and ref '#' has a value the 
                       '#' value will also be taken as the diff value.

GHz   MHz  cm-1        Units for the Frequency

Without Keywords the input is straightforward:

  J K- K+ J K- K+ Frequency Error; all other options must follow with 
                                    keywords

Example (are lines have the same meaning): 
  
   2 1 2  1 0 1  15.785559  
   2 1 2  1 0 1  15785.559  MHz  
   Jup 2  K-up  1  K+up 2   Jlo 1  K-lo 0  K+lo 1  =  15.785559 
   J   2  K-    1  K+   2   J   1  K-   0  K+   1  =  15.785559 
    =  15.785559    J   2  K-    1  K+   2   J   1  K-   0  K+   1  
   J 2 1  K- 1 0  K+ 2 1    =  15.785559   
   J 2 1  t  2 1            =  15.785559   


[1] R.C.Woods, J.Mol.Spectrosc., 21, 4 (1966).
[2] R.C.Woods, J.Mol.Spectrosc., 22, 49 (1967).
[3] J.M.Vacherand, B.P.van Eijck, J.Burie, and J.Demaison, J.Mol.Spectrosc.,
    118, 355 (1986).
[4] J.K.G.Watson, in: Vibrational Spectra and Structure (J.R.Durig, ed.)
    Vol.6, p.39, Elsevier, Amsterdam 1977.
[5] H.P.Benz, A.Bauder, Hs.H.Guenthard, J.Mol.Spectrosc., 21, 156-164 (1966).
[6] H.Hartwig and H.Dreizler, Z.Naturforsch, 51a (1996) 923.


APPENDIX

1. Output file 'eval.out'

This file contains the eigenvalues which are obtained during the calculation.
J : J Quantum no.
S : Symmetry quantum no. S=0: Eigenvalue without any internal rotation
                         S=1: normally A-Transition ... defined in block 5. 
B : 
F : double F-Value (nuclear quadrupole coupling)
T : number of eigenvalues in one matrix (ascending)
t : like T, but ranging only from 1 to 2J+1. t restarts for every v block.
K : assigned K quantum no. 
best K(s)   vec : 
 The K value of the basis function which has the strongest contribution 
 to the eigenvalue is printed. vec is the coefficient of this basis function. 
 If two basis functions have similar contributions both are printed.
V  vec :
 The v value of the eigenvector and its coefficient.
2nd.K vec :
 The second important K Basis function and its coefficient

 Iteration    1

  J  S  B  F    T  t  K  V    Energy/GHz    best K(s)      vec   V   vec 2nd.K vec 2nd.V vec
  1  0  1  2    1  1  0  1    12.63588736      0         1.000   1 0.000   1 0.000   1 0.000
  1  0  1  2    2  2 -1  1    23.97366175   1 -1  0.707 -0.707   1 0.000             1 0.000
  1  0  1  2    3  3  1  1    24.70430165   1 -1 -0.707 -0.707   1 0.000             1 0.000

2. Output control with the aprint and xprint keywords
   aprint
   1-7 AP_PL    ! Parameter List (0..3)
     8 AP_TF    ! List of Transitions zero cycle 
    16 AP_TL    ! List of Transitions 
    32 AP_TE    ! List of Transitions extended
    64 AP_PC    ! additional Parameter information
   128 AP_IO    ! Echo of input transition list 
   256 AP_LT    ! Latex Output (at the end only)
   512 AP_SV    ! SVD-Information
  1024 AP_ST    ! status
  2048 AP_TI    ! Torsional Integrals < v K | v' K' > 
  4096 AP_RM    ! Rotation Matrix D   
  8192 AP_EO    ! Eigen energies of one top operator
 16384 AP_MO    ! Matrix elements of one top operator  
 32768 AP_EH    ! Eigen energies / vectors    
 65536 AP_MH    ! Matrix elements     

   xprint
     1 XP_FI    ! first cycle
     2 XP_LA    ! last cycle
     4 XP_CC    ! every converged cycle
     8 XP_EC    ! every cycle
    16 XP_DE    ! for every dqu derivative


3. Examples of input files: 
-------------------------- snip snip -------------------------------
  A simple rigid rotor prediction        !(1-st line of title )
  A-species of propylenoxid              !(2-nd line of title )
                                         !(one empty line )
 int   2        ! keyword for intensity calculation
 ntop  0        ! no internal rotor 
 freq 12.0 18.0 ! Frequency region 
 limit 0.5      ! lower limit for intensities
 temp 20.0      ! temperature in K for the Boltzmann factor 
                                         !(one empty line ) 
 BJ         6.316693679 
 BK        11.706844345
 B-         0.365382450
 mu_x  1.67                   ! Dipole moment
 mu_y  0.56
 mu_z  0.95
                                         !(one empty line )

 J 5                         
