# Set project directory one level above of Makefile directory. $(CURDIR) is a GNU make variable containing the path to the current working directory
PROJDIR := $(realpath $(CURDIR)/../)
SOURCEDIR := $(PROJDIR)/src/c
BUILDDIR := $(PROJDIR)/bld

# Name of the final executable
TARGET =$(BUILDDIR)/ForestFV_5

# Decide whether the commands will be shwon or not
VERBOSE = TRUE

# Create the list of directories
DIRS =./
SOURCEDIRS = $(foreach dir, $(DIRS), $(addprefix $(SOURCEDIR)/, $(dir)))
TARGETDIRS = $(foreach dir, $(DIRS), $(addprefix $(BUILDDIR)/, $(dir)))
SOURCEDIRS = $(SOURCEDIR)
TARGETDIRS = $(BUILDDIR)

# Generate the GCC includes parameters by adding -I before each source folder
INCLUDES = $(foreach dir, $(SOURCEDIRS), $(addprefix -I, $(dir)))
INCLUDES= -I/home/ek/Local/petsc_intallation/arch-linux2-c-debug -I/home/ek/Local/petsc_intallation/include/ -I/home/ek/Local/petsc_intallation/arch-linux2-c-debug/lib -L/home/ek/Local/petsc_intallation/arch-linux2-c-debug/lib -I/home/ek/Local/parmets/metis/install/include/ -L /home/ek/Local/parmets/metis/install/lib/  -I/home/ek/Local/parmets/install/include/ -L /home/ek/Local/parmets/install/lib/ 

LIB_DIR=/home/ek/Local/petsc_intallation/arch-linux2-c-debug -I/home/ek/Local/petsc_intallation/include/

# Add this list to VPATH, the place make will look for the source files
VPATH = $(SOURCEDIRS)

# Create a list of *.c sources in DIRS
SOURCES = $(foreach dir,$(SOURCEDIRS),$(wildcard $(dir)/*.c))

# Define objects for all sources
OBJS := $(subst $(SOURCEDIR),$(BUILDDIR),$(SOURCES:.c=.o))

# Define dependencies files for all objects
DEPS = $(OBJS:.o=.d)

# Name the compiler
#CC =/home/ek/Local/petsc_intallation/install/bin/mpicc -Ofast
CC =/home/ek/Local/petsc_intallation/install/bin/mpicc
#CC =/home/ek/Local/petsc_intallation/install/bin/mpicc -fsanitize=address



RM = rm -rf 
RMDIR = rm -rf 
MKDIR = mkdir -p
ERRIGNORE = 2>/dev/null
SEP=/

FLAGS ?= $(INC_FLAGS) -Wl,-rpath,/home/ek/Local/petsc_intallation/arch-linux2-c-debug/lib -Wl,-rpath,/home/ek/Local/petsc_intallation/arch-linux2-c-debug/lib -L/home/ek/Local/petsc_intallation/arch-linux2-c-debug/lib -Wl,-rpath,/home/ek/Local/petsc_intallation/install/lib -L/home/ek/Local/petsc_intallation/install/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lpthread -lscalapack -lsuperlu_dist -lflapack -lfblas -lpthread -lparmetis -lmetis -lm -lX11 -lstdc++ -ldl -lmpifort -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -ldl

# Remove space after separator
PSEP = $(strip $(SEP))

# Hide or not the calls depending of VERBOSE
ifeq ($(VERBOSE),TRUE)
    HIDE =  
else
    HIDE = @
endif

# Define the function that will generate each rule
define generateRules
$(1)/%.o: %.c strdata.h
	@echo Building $$@
	$(HIDE) $(CC) $(FLAGS)  $(INCLUDES) -c -o $$(subst /,$$(PSEP),$$@) $$(subst /,$$(PSEP),$$<) -MMD
endef

.PHONY: all clean directories 

all: directories $(TARGET)

$(TARGET): $(OBJS)
	$(HIDE)echo Linking $@
	$(CC)  $(INCLUDES)  $(OBJS) -o $(TARGET) $(FLAGS)

# Include dependencies
-include $(DEPS)

# Generate rules
$(foreach targetdir, $(TARGETDIRS), $(eval $(call generateRules, $(targetdir))))

directories: 
	$(HIDE)$(MKDIR) $(subst /,$(PSEP),$(TARGETDIRS)) $(ERRIGNORE)


# Remove all objects, dependencies and executable files generated during the build
clean:
	$(HIDE)$(RMDIR) $(subst /,$(PSEP),$(TARGETDIRS)) $(ERRIGNORE)
	$(HIDE)$(RM) $(TARGET) $(ERRIGNORE)
	@echo Cleaning done ! 
