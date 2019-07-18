% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S4RRPR2
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% Jg_rot [3x4]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S4RRPR2_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_jacobig_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR2_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_jacobig_rot_sym_varpar: pkin has to be [5x1] (double)');
%% Function calls
if link_index == 0
	Jg_rot=S4RRPR2_jacobig_rot_0_sym_varpar(qJ, pkin);
elseif link_index == 1
	Jg_rot=S4RRPR2_jacobig_rot_1_sym_varpar(qJ, pkin);
elseif link_index == 2
	Jg_rot=S4RRPR2_jacobig_rot_2_sym_varpar(qJ, pkin);
elseif link_index == 3
	Jg_rot=S4RRPR2_jacobig_rot_3_sym_varpar(qJ, pkin);
elseif link_index == 4
	Jg_rot=S4RRPR2_jacobig_rot_4_sym_varpar(qJ, pkin);
else
	Jg_rot=NaN(3,4);
end