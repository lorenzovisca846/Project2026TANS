Come compilare ed eseguire il codice di simulazione. Dalla cartella home del progetto:

mkdir simulation/build
cd simulation/build
cmake ..
make

A questo punto il programma può essere eseguito dalla cartella build con il comando:

./simulation


Nota: nella cartella simulation è presente il file simConfig.txt contenente le impostazioni della simulazione. In alternativa è possibile utilizzare un file personalizzato, il nome del file in questo caso deve essere passato in fase di esecuzione del file:

./simulation PATH

dove PATH è il path del file rispetto alla cartella config