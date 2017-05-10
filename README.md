# group3_2017_healpix

# Scientific context


A gamma-ray transient is an astrophysical source that shows variable behavior, i.e. a significant change in the level of gamma-ray emission (i.e. several times the known quiescent flux level or a significant detection on short timescales) https://youtu.be/u37-D55SZu8 . The transient source could be associated with a known or unknown counterparts and could be of galactic or extragalactic origin, can emit coherent or incoherent radiation and even non electromagnetic signals, and can be the result of thermal runaways, explosions, and particle acceleration. Variability time scales, released energy, the wavelength of the emitted radiation can vary widely depending on the nature of the progenitor and of the physical process. In general, they are associated with catastrophic events involving compact objects, such as white dwarves, neutron stars, and black holes, and as such they offer the possibility to study the most extreme physical conditions in the Universe.
The time variability of these phenomena in the gamma-ray domain span from seconds to days. To capture these phenomena during their evolution and for effective communication to the astrophysical community, the speed is crucial and requires automated data analysis pipelines that will automatically detects and generates science alerts [1], [2]. These automated systems search for a transient event should be performed on multiple time scales (from seconds to hours). These systems will issue science alerts (a significant detection of a transient source that requires a follow-up strategy) immediately upon detection of γ-ray flares or that allow following external triggers in real-time.
The INAF/IASF Bologna institute is involved in different astrophysical projects; between them these proposals are applicable to the following gamma-ray astrophysical projects:
1) Given the flaring and transient nature of many VHE sources, the CTA Observatory [3] (https://youtu.be/1Q9dCqPNAXQ) will couple the unprecedented sensitivity to another performance driver: the fast reaction to VHE events, achieved using the Real-Time Analysis (RTA) science alert system [2], the key CTA system for the fast identification of flaring events. The Real-Time Analysis (RTA) is a data analysis pipeline that will automatically detects and generates science alerts with a maximum latency of 30 seconds with respect to the triggering event collection and ensure fast communication to/from the astrophysics community. 
2) The AGILE space mission (https://youtu.be/KyvOvZc53z8), sensitive in the energy range of 30 MeV-30 GeV, has detected many gamma-ray transients of galactic and extragalactic origins. Data are automatically analyzed at every orbital downlink by an alert pipeline operating on different timescales. As proper flux thresholds are exceeded, alerts are automatically generated and sent as SMS messages to cellular telephones, e-mails, and push notifications of an application for smartphones and tablets. These alerts are cross checked with the results of two pipelines, and a manual analysis is performed. Being a small scientific-class mission, AGILE is characterized by optimization of both scientific analysis and ground-segment resources. The system is capable of generating alerts within two to three hours of a data downlink, an unprecedented reaction time in gamma-ray astrophysics.

# Project

Il contesto astrofisico è quello descritto in precedenza. In questo caso il focus è sulla rappresentazione dell'immagine. Si parte direttamente rappresentando l'immagine astronomica in rappresentazione Healpix (http://healpix.jpl.nasa.gov). Una volta ottenuta questa, si possono poi fare operazioni semplici (e.g. Non-linear stretching e/o smoothing) su una immagine in rappresentazione Healpix. L’input è l'elenco dei fotoni gamma (tempo di arrivo, posizione nel cielo ed energia del fotone) e si costruisce la sky mappa, e poi si può applicare una delle operazioni insegnate al corso. 

Ci sono dei visualizzatori che permettono di vedere il risultato. 
Nota: lo smoothing richiede di conoscere l'adiacenza dei pixel, ma le librerie dovrebbero già fornire tutto quanto. È possibile fermarsi al non-linear stretching: in questo caso si potrebbe 'usare' l'energia dei fotoni per fare uno stretching più complesso (ogni pixel contiene il numero di fotoni presenti in quel pixel, cioè in un'area di cielo, senza considerare la loro energia).

# Attività Progettuale
Esperienza sul deep learning. Si potrebbero generare una serie di immagini astronomiche simulate, alcune con solo background e altre con background più un segnale di un certo livello da una sorgente astrofisica (si potrebbero creare classi diverse in base all'intensità del segnale dalla sorgente astrofisica) e poi provare il Deep Learning su questi set di immagini, facendo misure di efficienza e di prestazioni. Il team AFILE ha anche un simulatore di immagini astronomiche in grado di produrre tutte le immagini necessarie (ed etichettate) i parametri che vogliamo.
L’attivitá andrebbe svolta su macchina Minsky (Power 8+ con 2 K100 su nvlink).

# Reference material
- formato di memorizzazione delle mappe: https://fits.gsfc.nasa.gov
- libreria C per aprire le mappe: https://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html. Esiste anche la versione C++: https://heasarc.gsfc.nasa.gov/docs/software/fitsio/ccfits/ ma non mi è chiaro quanto è mantenuta. Dateci un occhio.
- ds9: un tool grafico per aprire le mappe. http://ds9.si.edu/site/Home.html

- Reference systems and projections: http://pro.arcgis.com/en/pro-app/help/mapping/properties/coordinate-systems-and-projections.htm
un minimo di info sui sistemi di riferimento celesti: https://en.wikipedia.org/wiki/Celestial_coordinate_system. Noi usiamo il Galattico, ma l’Equatoriale è forse il più diffuso
- Map projections systems: aitoff projection: https://en.wikipedia.org/wiki/Aitoff_projection, arc projection
- Healpix (http://healpix.jpl.nasa.gov),

# References
[1] A. Bulgarelli, et al., doi:10.1088/0004-637X/781/1/19 or https://arxiv.org/abs/1401.3573
[2] A. Bulgarelli, at al., https://arxiv.org/abs/1509.01963
[3] B. S. Acharya, et al., Introducing the CTA concept, Astr. Phys., 43, 3-18, 2013.



