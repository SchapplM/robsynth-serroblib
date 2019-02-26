% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4RPRP2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% Ja_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S4RPRP2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPRP2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPRP2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_jacobia_transl_sym_varpar: pkin has to be [5x1] (double)');
%% Function calls
if link_index == 0
	Ja_transl=S4RPRP2_jacobia_transl_0_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 1
	Ja_transl=S4RPRP2_jacobia_transl_1_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 2
	Ja_transl=S4RPRP2_jacobia_transl_2_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 3
	Ja_transl=S4RPRP2_jacobia_transl_3_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 4
	Ja_transl=S4RPRP2_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, pkin);
else
	Ja_transl=NaN(3,4);
end