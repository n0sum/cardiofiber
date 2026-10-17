# Dipendenze di terze parti

## happly

Lettore/scrittore PLY a header singolo, di Nicholas Sharp.

- upstream : https://github.com/nmwsharp/happly
- licenza  : MIT
- commit   : 8a606309daaa680eee495c8279feb0b704148f4a

Era un git submodule nel repository precedente. E' stato **vendorizzato**
(copiato in `include/`) perche' il repository compili senza rete e senza
dipendere da un upstream che puo' spostarsi: condizione necessaria perche'
la build sia riproducibile a distanza di tempo.

Per aggiornarlo: scaricare `happly.h` dall'upstream, sostituire il file e
aggiornare il commit qui sopra.
