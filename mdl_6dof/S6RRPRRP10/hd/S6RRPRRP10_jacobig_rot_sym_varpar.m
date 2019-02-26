% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
% 
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRP10_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP10_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
%% Function calls
if link_index == 0
	Jg_rot=S6RRPRRP10_jacobig_rot_0_sym_varpar(qJ, pkin);
elseif link_index == 1
	Jg_rot=S6RRPRRP10_jacobig_rot_1_sym_varpar(qJ, pkin);
elseif link_index == 2
	Jg_rot=S6RRPRRP10_jacobig_rot_2_sym_varpar(qJ, pkin);
elseif link_index == 3
	Jg_rot=S6RRPRRP10_jacobig_rot_3_sym_varpar(qJ, pkin);
elseif link_index == 4
	Jg_rot=S6RRPRRP10_jacobig_rot_4_sym_varpar(qJ, pkin);
elseif link_index == 5
	Jg_rot=S6RRPRRP10_jacobig_rot_5_sym_varpar(qJ, pkin);
elseif link_index == 6
	Jg_rot=S6RRPRRP10_jacobig_rot_6_sym_varpar(qJ, pkin);
else
	Jg_rot=NaN(3,6);
end