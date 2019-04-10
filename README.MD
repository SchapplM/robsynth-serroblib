# Bibliothek serieller Roboter

Dieses Repo enthält eine Sammlung von Modellen serieller Roboter. Die Modelle beinhalten symbolisch generierte Kinematik und Dynamik der Roboter. Außerdem stehen Parameter-Datensätze für existierende Roboter zur Verfügung.

Moritz Schappler, schappler@imes.uni-hannover.de, 2018-08  
(C) Institut für mechatronische Systeme, Universität Hannover

## Benutzung

Für die Benutzung des Repos muss mit dem Skript `serroblib_path_init.m` aus dem Hauptordner die Pfadinitialisierung durchgeführt werden.

### Abhängigkeiten

Die Robotermodelle setzen die Toolbox imes-matlab und imes-matlab-robotics mit initialisiertem Pfad voraus.

## Inhalt

Die Robotermodelle sind nach der Anzahl der Gelenk-FG der Roboter aufgeteilt und in csv-Dateien tabellarisch aufgelistet.
Zum Eintragen, Auslesen und Filtern von Robotermodellen existieren Funktionen.

Im Folgenden wird unterschieden zwischen

* Typ: Gekennzeichnet durch Anzahl und Art der Gelenke. Z.B. `S4RRPR` (S=Seriell, 4=4FG, RRPR=Reihenfolge der Gelenke)
* Modell: Genaue Anordnung der Gelenke (Senkrecht, Parallel, Verschoben, in einem Punkt). Z.B. `S4RRPR1` als eine Möglichkeit, die Gelenke R,R,P,R in dieser Reihenfolge anzuordnen.
* Roboter: Parametriertes Robotermodell. Z.B. "Kuka KR 5 scara R350 Z320" für Kuka-Roboter mit Kinematik vom Typ RRPR1

### Funktionen

Die genauen Übergabeargumente und Funktionsweisen sind dem Quelltext der Funktionen und Beispielaufrufen zu entnehmen.

Datenbank aufbauen (z.B. mit Struktursynthese-Ergebnissen):
* `serroblib_add_robot`: Hinzufügen eines neuen Roboters zur Datenbank. Der Roboter wird mit modifizierten DH-Parametern nach [KhalilKle1986] beschrieben.

Datenbank sortieren:
* `serroblib_gen_bitarrays`: Robotermodelle aus den csv-Tabellen effizient als Bit-Array in Matlab-Datendateien speichern. Muss durchgeführt werden, bevor die Datenbank ausgelesen wird.
* `serroblib_create_robot_csv_all`: Skript zum Schreiben von csv-Tabellen zur Beschreibung der einzelnen Roboter zu den Modellen.
* `serroblib_remove_robot`: Robotermodell aus Datenbank entfernen.

Datenbank auslesen und benutzen:
* `serroblib_filter_robots`: Liste von Robotern nach Filterkriterien (wie Anzahl und Art der FG) ausgeben
* `serroblib_create_robot_class`: Instanz einer Matlab-Klasse für einen Roboter erstellen, damit sind dann weitere Berechnungen mit dem Modell des Roboters möglich
* `serroblib_generate_mapleinput`: Eingabedaten zur Code-Generierung für diesen Roboter erstellen
* `serroblib_generate_code`: Kinematik- und Dynamikfunktionen symbolisch mit der HybrDyn-Toolbox berechnen lassen.
* `serroblib_addtopath`: Matlab-Funktionen dieses Roboters zum Matlab-Pfad hinzufügen, damit die Funktionsaufrufe funktionieren.

Zusätzliche Eigenschaften des Roboters ermitteln und ebenfalls eintragen (ist notwendig, da manche Eigenschaften erst mit den generierten Matlab-Funktionen berechnet werden können): 
* `scripts/correct_phi_N_E.m`
* `scripts/determine_jointnumber_influence_ee_position.m`

Datenbank prüfen:
* Mit `scripts/check_model_enumeration.m` lässt sich prüfen, ob nach Lösung von Kinematiken Lücken in der Nummerierung auftreten.

Datenbank modifizieren: Wenn nachträglich das Format der Datenbank geändert wird, geschieht dies automatisiert über Skripte (Ordner `scripts`). Die folgenden Skripte können als Beispiele für weitere Modifikationen dienen:
* `add_csv_column_ee_phi.m` und `correct_phi_N_E.m`: Hinzufügen und Belegen von Spalten für zusätzliche EE-Drehung
* `add_csv_column_pos_infl_joint.m` und `determine_jointnumber_influence_ee_position.m`: Zusätzliche Spalte für Nummer des Gelenks, das zuletzt die EE-Position beeinflusst.

### Übersichtsdateien

Alle Informationen der Roboterbibliothek liegen in csv-Dateien (Trennzeichen: Semikolon) vor. Dadurch ist eine automatische Verarbeitung und eine Versionsverwaltung möglich.

#### Liste der Robotermodelle für eine Gelenkfolge

Die Datei `mdl_xdof/SxRRRRRR.csv` enthält die MDH-Parameter aller Robotermodelle mit der entsprechenden Gelenkfolge. Dabei ist

* S die Kennzeichnung für serielle Roboter,
* x die Anzahl der Gelenk-FG,
* R/P die Kennzeichnung für Dreh- oder Schubgelenk.

Die Tabelle enthält die MDH-Parameter in symbolischer Form und zeigt die kinematische Struktur anhand von 0-Einträgen und 90°-Einträgen in den Parametern.

#### Liste der Roboterparameter für ein Robotermodell

Die Datei `mdl_xdof/SxRRRRRRyyyyy/models.csv` enthält mögliche Roboterparameter für das entsprechende Robotermodell. Dabei ist

* die Bezeichnung des Modells wie oben,
* y die laufende Nummer des Robotermodells.

Die Tabelle enthält Zahlenwerte für die einzelnen MDH-Parameter.

#### Matlab-Datei mit Modellübersicht

Die csv-Dateien enhalten alle Informationen der Roboterbibliothek. Die Daten werden zusätzlich in komprimierter Form gespeichert. In `mdl_xdof/Sx_list.mat` liegen die gleichen Daten wie in den csv-Tabellen zur einfacheren Verarbeitung in Matlab.

### Kinematik- und Dynamikfunktionen

Die Modelle für Kinematik und Dynamik, die den Kern dieser Sammlung darstellt, sind in den jeweiligen Modellordnern `mdl_xdof/SxRRRRRRyyyyy` im Unterordner `hd` (Abkürzung für "HybrDyn") als Matlab-Funktionen abgelegt.
Die Bedeutung der Funktionsnamen ist in der README von HybrDyn dargestellt. Es werden Abkürzungen für Funktionsart, Art der Herleitung und verwendete Parameter verwendet.

Zusätzlich werden hier auch Simulink-Bibliotheken abgelegt (mdl-Format, beginnend mit dem Präfix "lib"). Die Simulink-Bibliotheken können für Simulationen oder auch direkt für Regler auf Echtzeitsystemen verwendet werden.

### CAD-Modelle

Zu den Parametrierten Robotermodellen können in der Datenbank auch CAD-Modelle abgelegt werden. Diese können dann mit Hilfe der Matlab-Klasse `SerRob` zur Visualisierung verwendet werden.
Die CAD-Modelle werden im Unterordner `CAD_Robotername` von `mdl_xdof/Typbezeichnung` abgelegt.  
In jedem CAD-Ordner liegen die CAD-Modelle der einzelnen Roboterkörper im STL-Format (binär, nicht ASCI) und eine Initialisierungsdatei `Robotername_init_CAD.m`. Dort werden z.B. die Koordinatentransformationen zwischen den Körper-Koordinatensystemen und dem Ursprung und der Ausrichtung in den STL-Dateien definiert.

### Tests

Im Ordner `examples_tests` befinden sich Testskripte, mit denen sich die Funktionen aus dieser Bibliothek ausführen lassen.

Weitere Beispiele befinden sich im Robotik-Repo für die Roboterklasse `SerRob` unter Benutzung der Modelle aus dieser Toolbox.

## Quellen

* [KhalilKle1986] Khalil, W. and Kleinfinger, J.-F.: "A new geometric notation for open and closed-loop robots" (1986) <a name="KhalilKle1986"></a> 
