% Analytische Jacobi-Matrix für beliebiges Segment von
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RPPRRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% Ja [6x6]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja = S6RPPRRP5_jacobia_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_jacobia_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP5_jacobia_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP5_jacobia_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_jacobia_sym_varpar: pkin has to be [8x1] (double)');
%% Function calls
if link_index == 0
	Ja=S6RPPRRP5_jacobia_0_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 1
	Ja=S6RPPRRP5_jacobia_1_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 2
	Ja=S6RPPRRP5_jacobia_2_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 3
	Ja=S6RPPRRP5_jacobia_3_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 4
	Ja=S6RPPRRP5_jacobia_4_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 5
	Ja=S6RPPRRP5_jacobia_5_sym_varpar(qJ, r_i_i_C, pkin);
elseif link_index == 6
	Ja=S6RPPRRP5_jacobia_6_sym_varpar(qJ, r_i_i_C, pkin);
else
	Ja=NaN(6,6);
end