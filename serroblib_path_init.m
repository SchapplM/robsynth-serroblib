% Diese Datei dient nur zum Auffinden dieses Ordners mit dem Befehl
% `which('serroblib_path_init.m')` und zum Hinzufügen des 
% Hauptverzeichnisses des Repos

this_tb_path = fileparts( mfilename('fullpath') );
addpath(this_tb_path);

% Alle csv-Tabellen (versionsverwaltet) nach .mat konvertieren (nicht
% versionsverwaltet, da binär-Format)
serroblib_gen_bitarrays