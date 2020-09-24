# DSToolbox

DSToolbox is a Matlab toolbox for analyzing and modelling of unsteady fluid dynamics experiments.

## Installation

Clone the repository in the folder of choice on your computer, for example `Documents/MATLAB`. Open a Terminal window and type:

```bash
cd Documents/MATLAB
git clone git@github.com:lucasschn/dstoolbox
```
You can now make sure that you have a new folder called `dstoolbox` at the specified location. 

## Organization 

The `src` folder contains the source code and the `script` folder contains scripts and examples that make use of the toolbox functionalities. The `src`folder is further divided in three subfolders: 
- `common` contains the classes definitions for the core-objects of the toolbox, e.g. Airfoil, AirfoilMotion, or SteadyCurv that are used independently of the model. 
- `lib` contains a library of useful functions, not classes, also common to all models.
- `model` contains functions that are model-specific, e.g. only used for Sheng and Expfit models


## Usage

The repository consists in a collection of objects, such as airfoils or typical motions, that can be created and on which functions can be applied. Scripts can then be written where these objects are created, such as in the example below: 

```matlab
airfoil = Airfoil('naca0012',0.5) % creates an Airfoil object with name naca0012 and 0.5m chord length
airfoil.steady = SteadyCurve(alpha,CN,13) % creates a SteadyCurve object
```
The created steady curve is assigned as the steady/static curve to the airfoil, with angle of attack alpha, normal coefficient CN, and static stall angle 13°. Different methods apply to steady curves, for additional computations or for plots. 

```matlab
airfoil.steady.plotCN() % plots the normal coefficient as a function of the AoA
airfoil.steady.fitKirchhoff() % fits a Kirchhoff curve to the static stall curve
```

An other object is needed to represent the dynamic stall experiment. It can be a ramp-up motion with constant pitch rate: 

```matlab
ramp = RampUpMotion('r',0.01,'V',0.5) % creates an ramp-up object with reduced pitch rate 0.01 and incoming flow velocity 0.5m/s.
```

a sinusoidal motion with constant frequency:
```matlab
pitching = PitchingMotion('alpha',alpha,'CN',CN,'k',red_freq) % creates a pitching motion object with angle of attack vector alpha, normal coefficient CN and reduced frequency red_freq.
```

or a general motion with custom angle of attack history:

```matlab
motion = AirfoilMotion('alpha',alpha,'CN',CN)
```
RampUpMotion and PitchingMotion both inherit from AirfoilMotion, meaning that all properties and methods of AirfoilMotion also apply to RampUpMotion and PitchingMotion. Howvever, RampUpMotion and PitchingMotion both individually have properties and methods that AirfoilMotion does not, such as the reduced pitch rate `r` and the reduced frequency `k` respectively. 

All three airfoil motions accept name-value pair arguments when constructed. This means that you can pass any `'name',value` pair as an argument when creating the object to automatically assign the value `value` to the property `name` to the object, as long as the property `name` exists for this object. 

Different methods can be applied to a newly created ramp object, such as `setCL()` for setting the experiment lift coefficient corresponding to this ramp manually. A convenient function `loadRamp(casenumber,filtered)` sets up the experimental data to the ramp automatically from the server data. 

```matlab
ramp = loadRamp(22,false);
ramp.setPitchRate(airfoil);
ramp.findExpOnset()
```

Here the number 22 defines the experimental case number corresponding to the desired experiement. All case numbers are defined in the labbook (`labbook.m`in the repository). `setPitchrate(airfoil)`must be executed independently because it requires an airfoil object as an argument (in order to define the reduced pitch rate, the chord length is required). This will also set the convectime time vector, which allows `findExpOnset()` to be run. It is recommended to take the habbit to declare a ramp using this three methods before any usage. 

Once the airfoil motion has been set up correctly, the corresponding aerodynamic normal coefficients can be predicted using a dynamic stall model. All dynamic stall models are methods that apply to motion objects. The general syntax for models is as follows: 

```matlab
ramp.BeddoesLeishman(airfoil,Tp,Tf,Tv,Tvl,'mode') % computes the aerodynamic loading experienced by an airfoil object describing the motion described by ramp
```
The time constants Tp, Tf, Tv, and Tvl are necessary input arguments to Beddoes-Leishman model. Depending on the selected model the number of time constants can vary from 3 to 4. The 'mode' argument can be either 'experimental' or 'analytical' depending if the user wants numerical or analytical derivatives to be used. 

## App

Before using the app, you have to run the script `setPaths.m`on your machine with the correct path to the folder where you want the produced figures to be saved.

## Troubleshooting

When running a file always make sure that your current folder is the folder containing the file. For example, many scripts won't execute correctly if your current folder is not `dstoolbox/scripts`.


If you see the error message

```
Matlab couldn't read the experimental data. Are you sure you are connected to the server?
```

Make sure you are connected to the raw server. Otherwise, open labbook.m and make sure the path to the smartH folder is correctly set. 


## License
[MIT](https://choosealicense.com/licenses/mit/)
