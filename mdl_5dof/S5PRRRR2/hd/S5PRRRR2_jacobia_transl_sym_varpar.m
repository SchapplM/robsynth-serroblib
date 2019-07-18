% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRR2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:30
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRRR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_jacobia_transl_sym_varpar: pkin has to be [6x1] (double)');
%% Function calls
if link_index == 0
	Ja_transl=S5PRRRR2_jacobia_transl_0_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 1
	Ja_transl=S5PRRRR2_jacobia_transl_1_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 2
	Ja_transl=S5PRRRR2_jacobia_transl_2_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 3
	Ja_transl=S5PRRRR2_jacobia_transl_3_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 4
	Ja_transl=S5PRRRR2_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 5
	Ja_transl=S5PRRRR2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, pkin);
else
	Ja_transl=NaN(3,5);
end