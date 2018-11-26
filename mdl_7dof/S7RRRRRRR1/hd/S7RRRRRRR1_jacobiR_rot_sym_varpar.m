% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% 
% Output:
% JR_rot [9x7]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 19:01
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JR_rot = S7RRRRRRR1_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),uint8(0),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiR_rot_sym_varpar: qJ has to be [7x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S7RRRRRRR1_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiR_rot_sym_varpar: pkin has to be [4x1] (double)');
%% Function calls
if link_index == 0
	JR_rot=S7RRRRRRR1_jacobiR_rot_0_sym_varpar(qJ, pkin);
elseif link_index == 1
	JR_rot=S7RRRRRRR1_jacobiR_rot_1_sym_varpar(qJ, pkin);
elseif link_index == 2
	JR_rot=S7RRRRRRR1_jacobiR_rot_2_sym_varpar(qJ, pkin);
elseif link_index == 3
	JR_rot=S7RRRRRRR1_jacobiR_rot_3_sym_varpar(qJ, pkin);
elseif link_index == 4
	JR_rot=S7RRRRRRR1_jacobiR_rot_4_sym_varpar(qJ, pkin);
elseif link_index == 5
	JR_rot=S7RRRRRRR1_jacobiR_rot_5_sym_varpar(qJ, pkin);
elseif link_index == 6
	JR_rot=S7RRRRRRR1_jacobiR_rot_6_sym_varpar(qJ, pkin);
elseif link_index == 7
	JR_rot=S7RRRRRRR1_jacobiR_rot_7_sym_varpar(qJ, pkin);
else
	JR_rot=NaN(9,7);
end