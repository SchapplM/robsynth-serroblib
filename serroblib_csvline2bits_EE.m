% Wandle eine Roboterstruktur aus gegebener CSV-Zeile in Bit-Array um
% Die csv-Zeile ist besser lesbar, das Bit-Array ist schneller
% Maschinen-Verarbeitbar
% 
% Eingabe:
% csvline
%   Cell-Array mit Spaltenweise EE-FG ('1' oder '0')
%   Oder: Vektor mit Spaltenweise 1,0
% 
% Ausgabe:
% BAE [1x1] uint16
%   Bit-Array zur Kennzeichnung der EE-FG
%   Bits:
%   01 (LSB): vx0
%   ...
%   06: wz0
% 
% Siehe auch: serroblib_bits2csvline

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-08
% (C) Institut für mechatronische Systeme, Universität Hannover

function BAE = serroblib_csvline2bits_EE(csvline)

%% Initialisierung
% Prüfe Eingabe
if length(csvline) ~= 6
  error('falsche Anzahl Einträge in csvline');
end

% Ausgabevariable vorbelegen
BAE = uint16(0);

%% Bit-Vektor für EE-FG aus csv-Zeile gewinnen
b = 0; % Bit-Offset zur Verschiebung der Parameter-Bits in der Gesamtvariable
c = 0;
for i = 1:6
  c = c+1;
  if strcmp(csvline{c}, '0')
    BAE = bitor( BAE, bitshift(0,b)); b = b+1;
  else
    BAE = bitor( BAE, bitshift(1,b)); b = b+1;
  end
end
