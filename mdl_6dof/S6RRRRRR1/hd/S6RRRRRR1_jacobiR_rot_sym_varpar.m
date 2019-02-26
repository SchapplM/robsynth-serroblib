% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
%% Function calls
if link_index == 0
	JR_rot=S6RRRRRR1_jacobiR_rot_0_sym_varpar(qJ, pkin);
elseif link_index == 1
	JR_rot=S6RRRRRR1_jacobiR_rot_1_sym_varpar(qJ, pkin);
elseif link_index == 2
	JR_rot=S6RRRRRR1_jacobiR_rot_2_sym_varpar(qJ, pkin);
elseif link_index == 3
	JR_rot=S6RRRRRR1_jacobiR_rot_3_sym_varpar(qJ, pkin);
elseif link_index == 4
	JR_rot=S6RRRRRR1_jacobiR_rot_4_sym_varpar(qJ, pkin);
elseif link_index == 5
	JR_rot=S6RRRRRR1_jacobiR_rot_5_sym_varpar(qJ, pkin);
elseif link_index == 6
	JR_rot=S6RRRRRR1_jacobiR_rot_6_sym_varpar(qJ, pkin);
else
	JR_rot=NaN(9,6);
end