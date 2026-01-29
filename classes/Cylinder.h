/*
    Classe per la simulazione:
    Data members:
        - double raggio
        - double spessore
        - double lunghezza
        - string materiale
        - double ms coefficient

    Oggetti di tipo cylinder:
        beampipe (Berillio)
        layer1 (Silicio)
        layer2 (Silicio)

    Eventualmente creiamo classi beampipe e layer che ereditano da cylinder (non mi sembra necessario)

    Metodi
        - costruttore
        - getter
        - calcolo del coefficiente di multiple scattering
            NOTA: soprattutto nei layer dei rivelatori la traiettoria non è radiale,
                  approssimiamo comunque la lunghezza del tragitto come lo spessore del layer?
*/