% Umwandlung als numerische Werte gegebene MDH-Parameter in eine Index-Struktur
% Damit können schnell aus der Modifikation von Instanzen der SerRob-Klasse
% neue Einträge (v.a. Varianten von Modellen) für die Modellbibliothek
% generiert werden.
% 
% Eingabe:
% MDH_struct_num
%   Struktur mit numerischen Werten für alle MDH-Parameter
%   (a,alpha,d,theta, ...), wie sie in SerRob verwendet wird.
% 
% Ausgabe:
% MDH_struct_idx
%   Struktur mit Index-Werten für die MDH-Parameter, die sich auf die
%   mit Bits kodierte Reihenfolge möglicher Werte in der Datenbank bezieht
%   (insbesondere in serroblib_gen_bitarrays, serroblib_bits2csvline

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-04
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MDH_struct_idx = serroblib_mdh_numparam2indexstruct(MDH_struct_num)
N = length(MDH_struct_num.d);
MDH_struct_idx = struct(          'type', 255*uint8(ones(1,N)), ...
  'beta',   255*uint8(ones(1,N)), 'b', 255*uint8(ones(1,N)), ...
  'alpha',  255*uint8(ones(1,N)), 'a', 255*uint8(ones(1,N)), ...
  'theta',  255*uint8(ones(1,N)), 'd', 255*uint8(ones(1,N)), ...
  'offset', 255*uint8(ones(1,N)));

% Die Indizes verweisen auf die Reihenfolge möglicher Zustände:
values_angles = [0, pi/2, pi, -pi/2, NaN];

% Eingabe-Struktur durchgehen und Werte für Ausgabe-Struktur belegen
for i = 1:N
  if MDH_struct_num.sigma(i) == 0, MDH_struct_idx.type(i) = 0;
  else,                            MDH_struct_idx.type(i) = 1; end
  
  idx_beta = find(values_angles == MDH_struct_num.beta(i))-1;
  if isnan(MDH_struct_num.beta(i)), idx_beta = 5-1; end
  MDH_struct_idx.beta(i) = idx_beta;
  
  idx_alpha = find(values_angles == MDH_struct_num.alpha(i))-1;
  if isnan(MDH_struct_num.alpha(i)), idx_alpha = 5-1; end
  MDH_struct_idx.alpha(i) = idx_alpha;
  
  idx_theta = find(values_angles == MDH_struct_num.theta(i))-1;
  if isnan(MDH_struct_num.theta(i)), idx_theta = 5-1; end
  MDH_struct_idx.theta(i) = idx_theta;
  
  if isnan(MDH_struct_num.b(i)), MDH_struct_idx.b(i) = 1;
  else,                          MDH_struct_idx.b(i) = 0; end
  
  if isnan(MDH_struct_num.a(i)), MDH_struct_idx.a(i) = 1;
  else,                          MDH_struct_idx.a(i) = 0; end  
  
  if isnan(MDH_struct_num.d(i)), MDH_struct_idx.d(i) = 1;
  else,                          MDH_struct_idx.d(i) = 0; end
  
  idx_off = find(values_angles == MDH_struct_num.offset(i))-1;
  MDH_struct_idx.offset(i) = idx_off;
end
