global inflowlib
global sym_calcoutput
global sym_end
global backup_Vx

path,_ = splitdir(@__FILE__)

ifw_active = true


"""
    ifwinit(inflowlib_filename ;HWindSpeed=6.87,turbsim_filename="path/test.bts")

calls inflow wind init

# Inputs
* `inflowlib_filename::string`: path and name of inflow-wind dynamic library
* `HWindSpeed::float`: optional, backup steady windspeed (m/s)
* `turbsim_filename::string`: path and name of turbsim data e.g. "path/test.bts"

# Outputs:
* `none`:

"""
function ifwinit(inflowlib_filename ;HWindSpeed=6.87,turbsim_filename="$path/test.bts")

    # Where the input is manipulated
    HWindSpeed_str = "       $(round(HWindSpeed,digits=1))   HWindSpeed     - Horizontal windspeed                            (m/s)"
    turbsim_str = "\"$turbsim_filename\"      filename_bts   - name of the full field wind file to use (.bts)"

    # NOTE: that there is an issue when passing in a binary coded string to fortran that an extra line/character is prepended, so we have to removed the first line to compensate
    #\x00------ InflowWind v3.01.* INPUT FILE -------------------------------------------------------------------------
    input_string_array = [codeunits("Steady 15 m/s winds with no shear for IEA 15 MW Offshore Reference Turbine
    \x00--------------------------------------------------------------------------------------------------------------
           \x00false  Echo           - Echo input data to <RootName>.ech (flag)
              \x003   WindType       - switch for wind file type (1=steady; 2=uniform; 3=binary TurbSim FF; 4=binary Bladed-style FF; 5=HAWC format; 6=User defined; 7=native Bladed FF)
              \x000   PropagationDir - Direction of wind propagation (meteoroligical rotation from aligned with X (positive rotates towards -Y) -- degrees)
              \x000   VFlowAng       - Upflow angle (degrees) (not used for native Bladed format WindType=7)
              \x001   NWindVel       - Number of points to output the wind velocity    (0 to 9)
              \x000   WindVxiList    - List of coordinates in the inertial X direction (m)
              \x000   WindVyiList    - List of coordinates in the inertial Y direction (m)
              \x001   WindVziList    - List of coordinates in the inertial Z direction (m)
    \x00================== Parameters for Steady Wind Conditions [used only for WindType = 1] =========================
            \x00$HWindSpeed_str
            \x00150   RefHt          - Reference height for horizontal wind speed      (m)
            \x000.0   PLexp          - Power law exponent                              (-)
    \x00================== Parameters for Uniform wind file   [used only for WindType = 2] ============================
    \x00\"unused\"      FileName_Uni   - Filename of time series data for uniform wind field.      (-)
            \x00150   RefHt_Uni      - Reference height for horizontal wind speed                (m)
         \x00125.88   RefLength      - Reference length for linear horizontal and vertical sheer (-)
    \x00================== Parameters for Binary TurbSim Full-Field files   [used only for WindType = 3] ==============
    \x00$turbsim_str
    \x00================== Parameters for Binary Bladed-style Full-Field files   [used only for WindType = 4] =========
    \x00\"unused\"      FilenameRoot   - Rootname of the full-field wind file to use (.wnd, .sum)
    \x00False         TowerFile      - Have tower file (.twr) (flag)
    \x00================== Parameters for HAWC-format binary files  [Only used with WindType = 5] =====================
    \x00\"unused\"      FileName_u     - name of the file containing the u-component fluctuating wind (.bin)
    \x00\"unused\"      FileName_v     - name of the file containing the v-component fluctuating wind (.bin)
    \x00\"unused\"      FileName_w     - name of the file containing the w-component fluctuating wind (.bin)
             \x0064   nx             - number of grids in the x direction (in the 3 files above) (-)
             \x0032   ny             - number of grids in the y direction (in the 3 files above) (-)
             \x0032   nz             - number of grids in the z direction (in the 3 files above) (-)
             \x0016   dx             - distance (in meters) between points in the x direction    (m)
              \x003   dy             - distance (in meters) between points in the y direction    (m)
              \x003   dz             - distance (in meters) between points in the z direction    (m)
            \x00150   RefHt_HAWC     - reference height; the height (in meters) of the vertical center of the grid (m)
     \x00-------------   Scaling parameters for turbulence   ---------------------------------------------------------
              \x002   ScaleMethod    - Turbulence scaling method   [0 = none, 1 = direct scaling, 2 = calculate scaling factor based on a desired standard deviation]
              \x001   SFx            - Turbulence scaling factor for the x direction (-)   [ScaleMethod=1]
              \x001   SFy            - Turbulence scaling factor for the y direction (-)   [ScaleMethod=1]
              \x001   SFz            - Turbulence scaling factor for the z direction (-)   [ScaleMethod=1]
            \x001.2   SigmaFx        - Turbulence standard deviation to calculate scaling from in x direction (m/s)    [ScaleMethod=2]
            \x000.8   SigmaFy        - Turbulence standard deviation to calculate scaling from in y direction (m/s)    [ScaleMethod=2]
            \x000.2   SigmaFz        - Turbulence standard deviation to calculate scaling from in z direction (m/s)    [ScaleMethod=2]
      \x00-------------   Mean wind profile parameters (added to HAWC-format files)   ---------------------------------
             \x0012   URef           - Mean u-component wind speed at the reference height (m/s)
              \x002   WindProfile    - Wind profile type (0=constant;1=logarithmic,2=power law)
            \x000.2   PLExp_HAWC     - Power law exponent (-) (used for PL wind profile type only)
           \x000.03   Z0             - Surface roughness length (m) (used for LG wind profile type only)
              \x000   XOffset        - Initial offset in +x direction (shift of wind box)
    \x00====================== OUTPUT ==================================================
    \x00False         SumPrint       - Print summary data to <RootName>.IfW.sum (flag)
                  \x00OutList        - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)
    \x00\"Wind1VelX,Wind1VelY,Wind1VelZ\"     - Wind velocity at point WindVxiList(1),WindVyiList(1),WindVziList(1).  X, Y, and Z direction components.
    \x00END of input file (the word \"END\" must appear in the first 3 columns of this last OutList line)
    \x00---------------------------------------------------------------------------------------")]


    input_string_length = [length(input_string_array[1])]
    # Only needed for WindType = 2, can leave it empty if not used, but still need as input
    # ifw_uniform_string_array = [""] # if not used
    # could be an arbitrary number of lines long
    # NOTE: that there is an issue when passing in a binary coded string to fortran that an extra line/character is prepended, so we have to removed the first line to compensate
    # ! OpenFAST InflowWind uniform wind input file for 15 m/s wind.\x00
    uniform_string_array = [codeunits("! Time Wind  Wind  Vert. Horiz. Vert. LinV  Gust   Upflow
    \x00!      Speed Dir   Speed Shear  Shear Shear Speed  Angle
    \x00! (sec) (m/s) (deg) (m/s) (-)    (-)   (-)  (m/s)  (deg)
      \x000.0  15.0  0.0   0.0   5.0    0.0   0.0   0.0    0.0
      \x000.1  15.0  0.0   0.0   0.0    0.0   0.0   0.0    0.0
      \x001.0  15.0  0.0   0.0   0.0    0.0   0.0   0.0    0.0")]

    uniform_string_length = [length(uniform_string_array[1])]
    numWindPts = [1] #[Int(length(positions)/3)]     # total number of wind points requesting velocities for at each time step. must be integer

    # Allocate Outputs
    numChannels = [0]
    channel_names = string(repeat(" ", 20 * 4000))
    channel_units = string(repeat(" ", 20 * 4000))
    error_status = [0]
    error_message = string(repeat(" ", 1025))


    try
        println("Attempting to access inflow wind at: $inflowlib_filename")
        global inflowlib = Libdl.dlopen(inflowlib_filename) # Open the library explicitly.
        sym_init = Libdl.dlsym(inflowlib, :IfW_C_Init)   # Get a symbol for the function to call.
        global sym_calcoutput = Libdl.dlsym(inflowlib, :IfW_C_CalcOutput)   # Get a symbol for the function to call.
        global sym_end = Libdl.dlsym(inflowlib, :IfW_C_End)

        ccall(sym_init,Cint,
            (Ptr{Cchar}, # IN: input_string_array
            Ptr{Cint},     # IN: input file string length
            Ptr{Cchar}, # IN: uniform_string_array
            Ptr{Cint},     # IN: input file string length
            Ptr{Cint}, # IN: numWindPts
            Ptr{Cdouble}, # IN: dt
            Ptr{Cint}, # OUT: numChannels
            Cstring, # OUT: channel_names
            Cstring, # OUT: channel_units
            Ptr{Cint}, # OUT: error_status
            Cstring), # OUT: error_message
            input_string_array, #InputFileStrings_C
            input_string_length,
            uniform_string_array, #InputUniformStrings_C
            uniform_string_length,
            numWindPts, #NumWindPts_C
            [0.0], #DT_C !!!NOT USED!!!
            numChannels,#NumChannels_C
            channel_names, #OutputChannelNames_C #TODO: none of these are coming through properly on the julia or python side
            channel_units, #OutputChannelUnits_C #TODO: none of these are coming through properly on the julia or python side
            error_status, #ErrStat_C
            error_message) #ErrMsg_C #TODO: none of these are coming through properly on the julia or python side
        catch
            global backup_Vx = HWindSpeed
            @warn "inflow wind library is not configured properly, using a steady wind profile [$backup_Vx,0.0,0.0] for Vx, Vy, Vz where Vx is initin[9]"
            global ifw_active = false
        end

end


"""
    ifwcalcoutput(position,time)

calls inflow wind clacoutput

# Inputs
* `position::Array(float)`: x, y, z sample position (m)
* `time::float`: sample time (s)

# Outputs:
* `velocities`: x, y, z velocity at sample position

"""
function ifwcalcoutput(position,time)

    outputChannelValues = zeros(3)

    #Reset error message
    error_message = string(repeat(" ", 1025))
    error_status = [0]
    velocities = zeros(Float32,3) # output velocities (N x 3)

    if ifw_active

        ccall( sym_calcoutput,Cint,
        (Ptr{Cdouble}, #Time_C
        Ptr{Cfloat}, #Positions_C
        Ptr{Cfloat}, #Velocities_C
        Ptr{Cfloat}, #OutputChannelValues_C
        Ptr{Cint}, #ErrStat_C
        Cstring), #ErrMsg_C
        [time],                 # IN: time at which to calculate velocities
        Float32.(position),                      # IN: positions - specified by user
        velocities,                     # OUT: velocities at desired positions
        outputChannelValues,                 # OUT: output channel values as described in input file
        error_status,              # ErrStat_C
        error_message)                     # ErrMsg_C

        return velocities
    else
        # @warn "using 10m/s inflow in the x direction since inflow wind library is misconfigured"
        return [backup_Vx,0.0,0.0]
    end
end

"""
    ifwend()

calls inflow wind end function and cleanup
"""
function ifwend()
    #Reset error message
    error_message = string(repeat(" ", 1025))
    error_status = [0]

    if ifw_active
        ccall( sym_end,Cint,
        (Ptr{Cint}, #ErrStat_C
        Cstring), #ErrMsg_C
        error_status,
        error_message)

        Libdl.dlclose(inflowlib) # Close the library explicitly.

    end
end


# Usage Example
# this was built from https://github.com/nrmendoza/openfast/tree/f/CCT2-IFW_wrapper
# mkdir build
# cd build
# cmake ..
# make ifw_c_binding
# the the library is located in: build/modules/inflowwind/libifw_c_binding.dylib
# and it can be moved anywhere
# import Libdl
# ifwstr = "$path/../../dylibs/libifw_c_binding"
# turbsim_filename = "$path/../test/data/ifw/test.bts"
# ifwinit(ifwstr;turbsim_filename)
# velocity = ifwcalcoutput([0.0,0.0,100.0],0.1)
# println(velocity)
# ifwend()
