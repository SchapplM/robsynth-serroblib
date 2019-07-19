% Diese Datei dient zum Auffinden dieses Ordners mit dem Befehl
% `which('serroblib_path_init.m')` und zum Initialisieren einiger Dateien
% Die Initialisierung der Dateien muss bei jeder Änderung im Repo gemacht
% werden und sollte daher in die startup.m aufgenommen werden.

this_tb_path = fileparts( mfilename('fullpath') );
addpath(this_tb_path);

% Alle csv-Tabellen (versionsverwaltet) nach .mat konvertieren (nicht
% versionsverwaltet, da binär-Format)
serroblib_gen_bitarrays

% Alle Vorlagen-Funktionen erstellen. Diese Funktionen sind nicht im Repo,
% da sie im wesentlichen identisch sind bis auf wenige Zeilen.
serroblib_create_template_functions