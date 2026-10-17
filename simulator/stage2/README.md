# Stadio 2: geometrie diverse

Da fare. Due assi geometrici controllati, non rumore casuale.

**Forme.** Istanze dell'SSM ottenute regolando i pesi dei modi, cosi' che
la distanza di forma dalla media sia un parametro e non un caso. La
corrispondenza con il riferimento e' l'identita' sugli indici.

**Risoluzioni.** Remeshing indipendente con lunghezza-target degli spigoli
variabile: vertici e connettivita' nuovi. Qui la corrispondenza non e' per
indice ma per proiezione al piu' vicino sul riferimento. La decimazione
che conserva i vertici non sarebbe un test valido, perche' renderebbe la
corrispondenza banale.

Su ogni geometria le fibre si **trasportano** dal riferimento
(push-forward col Jacobiano del triangolo, rinormalizzato; trasversale
dalla nuova normale), mai si rigenerano: rigenerarle introdurrebbe la
variabilita' delle fibre come confondimento. La LAT invece si **ri-simula**,
perche' dipende dalla geometria e non e' una proprieta' materiale.

Questa cartella esiste in git grazie a questo file: git non versiona le
directory vuote.
