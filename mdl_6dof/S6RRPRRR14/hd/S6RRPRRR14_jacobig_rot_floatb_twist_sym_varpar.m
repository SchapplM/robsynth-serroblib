% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR14
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:26
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR14_jacobig_rot_floatb_twist_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobig_rot_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR14_jacobig_rot_floatb_twist_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobig_rot_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');
%% Function calls
if link_index == 0
	Jg_rot=S6RRPRRR14_jacobig_rot_0_floatb_twist_sym_varpar(qJ, pkin);
elseif link_index == 1
	Jg_rot=S6RRPRRR14_jacobig_rot_1_floatb_twist_sym_varpar(qJ, pkin);
elseif link_index == 2
	Jg_rot=S6RRPRRR14_jacobig_rot_2_floatb_twist_sym_varpar(qJ, pkin);
elseif link_index == 3
	Jg_rot=S6RRPRRR14_jacobig_rot_3_floatb_twist_sym_varpar(qJ, pkin);
elseif link_index == 4
	Jg_rot=S6RRPRRR14_jacobig_rot_4_floatb_twist_sym_varpar(qJ, pkin);
elseif link_index == 5
	Jg_rot=S6RRPRRR14_jacobig_rot_5_floatb_twist_sym_varpar(qJ, pkin);
elseif link_index == 6
	Jg_rot=S6RRPRRR14_jacobig_rot_6_floatb_twist_sym_varpar(qJ, pkin);
else
	Jg_rot=NaN(3,6);
end