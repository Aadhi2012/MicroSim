#!/bin/sh

##
# @file
# @brief To set all microsim commands to shell scope.
#
# Usage : source setvars.sh $(pwd)	
#	
# @author Chirantandip Mahanta
# @date 2024-09-26
#

# Define MicroSim Paths

MICROSIM_HOME="${1:-${HOME}/MicroSim}" ;
MICROSIM_BIN_DIR="${2:-${HOME}/.local/bin}" ;
mkdir -p $MICROSIM_BIN_DIR ;

# Define functions in shell scope

microsim-export-path(){
    ## Initialize the paths
    if [[ "$1" == "" ]]; then
        echo "- **WARNING** MICROSIM_HOME not provided!"
        echo "- choosing default"
        MICROSIM_HOME=${HOME}/MicroSim
    else
        MICROSIM_HOME=$1
    fi
    echo "- MICROSIM_HOME = $MICROSIM_HOME" ;
    if [[ "$2" == "" ]]; then
        echo "- **WARNING** MICROSIM_BIN_DIR not provided!"
        echo "- choosing default"
        MICROSIM_BIN_DIR=${HOME}/.local/bin/
    else
        MICROSIM_BIN_DIR=$2
    fi
    echo "- MICROSIM_BIN_DIR = $MICROSIM_BIN_DIR" ;
    mkdir -p $MICROSIM_BIN_DIR ;
    ## Check if the paths exist
    if [[ ! -d "$MICROSIM_HOME" ]]; then
        echo "- **WARNING** path $MICROSIM_HOME does not exist." ;
        return 1 ;
    fi
    if [[ ! -d "$MICROSIM_BIN_DIR" ]]; then
        echo "path $MICROSIM_BIN_DIR could not be created." ;
        return 1 ;
    fi
    ## Export the paths
    export MICROSIM_HOME ;
    export MICROSIM_BIN_DIR ;
    ## Export bin path to $PATH
    if [[ ":$PATH:" != *":$MICROSIM_BIN_DIR:"* ]]; then
        export PATH="${PATH}:${MICROSIM_BIN_DIR}" ;
        echo path "${MICROSIM_BIN_DIR}" exported to \$PATH ;
    else
        echo "- $MICROSIM_BIN_DIR is already in the \$PATH" ;
    fi
}

microsim-show-command-paths() {
    microsim-print-path(){
        local filename=$(basename "$1") ;
        local padding=$((20 - ${#filename})) ;
        printf "command: $filename""%s%${padding}s""<-- ""$2\n" ;
    }
    if [[ ! -d "$MICROSIM_BIN_DIR" ]]; then
        echo "- **ERROR** path $MICROSIM_BIN_DIR does not exist." ;
        return 1 ;
    fi
    find "$MICROSIM_BIN_DIR" -type f -executable -exec file {} + | grep -E 'Bourne-Again shell script|bash script' | awk -F: '{print $1}' | while read -r bash_file; do
        local python_path=$(grep -Eo 'python3\s+"[^"]+"' "$bash_file" | sed 's/python3\s*"\([^"]*\)"/\1/') ;
        if [[ -n "$python_path" ]]; then
            microsim-print-path ${bash_file} ${python_path} ;
        fi
    done
    for link_file in "$MICROSIM_BIN_DIR"/*; do
        if [[ -L "$link_file" ]]; then
            local original_path=$(readlink "$link_file") ;
            microsim-print-path ${link_file} ${original_path} ;
        fi
    done
}

microsim-show-env(){
    echo "" ;
    echo "MicroSim Environment :" ;
    echo "----------------------" ;
    echo "MICROSIM_HOME    = $MICROSIM_HOME" ;
    echo "MICROSIM_BIN_DIR = $MICROSIM_BIN_DIR" ;
    microsim-show-command-paths ;
    echo "----------------------" ;
    echo "" ;
}

microsim-init(){

# Define functions in this subshell scope
    microsim-print-path(){
        local filename=$(basename "$1") ;
        local padding=$((20 - ${#1})) ;
        printf "command: $(basename "$1")""%s%${padding}s""<-- ""$2\n" ;
    }
    
    echo "- Initializing paths ..." ;
    microsim-export-path $@
    if [[ $? -eq 1 ]]; then
        echo "- **ERROR** in exporting paths! Aborting initialization."
        return 1
    fi


    echo "- Initializing commands ..." ;


    local MICROSIM_SCRIPT_DIR="${MICROSIM_HOME}/post" ;
    if [[ ! -d "$MICROSIM_SCRIPT_DIR" ]]; then
        echo "- **ERROR** path $MICROSIM_SCRIPT_DIR does not exist." ;
        return 1 ;
    fi

    for script in "$MICROSIM_SCRIPT_DIR"/*.py; do
        local base_name=$(basename "$script" .py) ;
        local executable="$MICROSIM_BIN_DIR/microsim-$base_name" ;
        cat > "$executable" << EOF
#!/bin/bash
python3 "$MICROSIM_SCRIPT_DIR/$base_name.py" "\$@"
EOF
        chmod +x "$executable" ;
        if [ $? -eq 0 ]; then
            microsim-print-path microsim-$base_name $script ;
        else
            echo "- **ERROR** Failed to create executable: $executable" ;
        fi
    done

    find "$MICROSIM_HOME" -type f -name "microsim*" -executable -exec file {} + | grep -E 'ELF' | awk -F: '{print $1}' | while read -r executable; do
        filename=$(basename "$executable") ;
        ln -sf "$executable" "$MICROSIM_BIN_DIR/$filename" ;
        microsim-print-path $filename $executable ;
    done

}

#
microsim-reset(){
    echo "Clearing MicroSim commands..."
    rm $MICROSIM_BIN_DIR/microsim* ;

    echo "Clearing MicroSim environment..."
    microsim-export-path ${HOME}/MicroSim ${HOME}/.local/bin
}

microsim-init $@
