% Wandle eine binär kodierte Zeile der EE-FG in eine CSV-Zeile um.
% Damit lässt sich aus den Bit-Arrays gut die vollständige CSV-Zeile
% (Kinematik und EE-FG) rekonstruieren
% 
% Eingabe:
% BAE [1x1] uint16
%   Bit-Array zur Kennzeichnung der EE-FG
%   Bits:
%   01 (LSB): vx0
%   ...
%   06: wz0
% 
% Ausgabe:
% csvline
%   Cell-Array mit Spaltenweise EE-FG ('1' oder '0')

% Siehe auch: serroblib_csvline2bits_EE

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-02
% (C) Institut für Mechatronische Systeme, Universität Hannover

function csvline_EE = serroblib_bits2csvline_EE(BAE)

csvline_EE = cell(1,6);
c=0;
for i = 1:6
  c = c+1; csvline_EE{c} = sprintf('%d', bitget(BAE,i));
end