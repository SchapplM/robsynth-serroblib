% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: S5RRRRR15_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% Jg_rot [3x5]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S5RRRRR15_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_jacobig_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR15_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
Jg_rot=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(5));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0; 1, cos(pkin(5)), 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (2->2), mult. (4->2), div. (0->0), fcn. (10->4), ass. (0->5)
	t70 = sin(pkin(5));
	t73 = cos(qJ(1)) * t70;
	t71 = cos(pkin(5));
	t69 = sin(qJ(1)) * t70;
	t1 = [0, t69, t69, 0, 0; 0, -t73, -t73, 0, 0; 1, t71, t71, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (3->3), mult. (6->2), div. (0->0), fcn. (15->4), ass. (0->5)
	t79 = sin(pkin(5));
	t82 = cos(qJ(1)) * t79;
	t80 = cos(pkin(5));
	t78 = sin(qJ(1)) * t79;
	t1 = [0, t78, t78, t78, 0; 0, -t82, -t82, -t82, 0; 1, t80, t80, t80, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:08
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (18->9), mult. (23->14), div. (0->0), fcn. (40->8), ass. (0->13)
	t153 = qJ(2) + qJ(3) + qJ(4);
	t152 = cos(t153);
	t157 = cos(pkin(5));
	t161 = t152 * t157;
	t155 = sin(pkin(5));
	t158 = sin(qJ(1));
	t150 = t158 * t155;
	t159 = cos(qJ(1));
	t160 = t159 * t155;
	t156 = cos(pkin(6));
	t154 = sin(pkin(6));
	t151 = sin(t153);
	t1 = [0, t150, t150, t150, t156 * t150 + (t151 * t159 + t158 * t161) * t154; 0, -t160, -t160, -t160, -t156 * t160 + (t151 * t158 - t159 * t161) * t154; 1, t157, t157, t157, -t155 * t152 * t154 + t157 * t156;];
	Jg_rot = t1;
end