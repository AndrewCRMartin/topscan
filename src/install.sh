BIN=${HOME}/bin
DATA=${HOME}/data

if [ ! -d $BIN ]; then
    echo "Creating binary directory: $BIN"
    mkdir -p $BIN
fi
if [ ! -d $DATA ]; then
    echo "Creating data directory: $DATA"
    mkdir -p $DATA
fi

cp -i topscan mergestride mergepdbsecstr $BIN
cp -i topmat.mat    $DATA
cp -i numtopmat.mat $DATA
