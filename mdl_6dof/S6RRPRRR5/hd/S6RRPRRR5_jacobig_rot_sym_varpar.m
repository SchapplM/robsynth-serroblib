% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR5
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:57
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR5_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR5_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t59 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, -cos(qJ(1)) * t59, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->5), mult. (22->12), div. (0->0), fcn. (35->8), ass. (0->12)
	t91 = sin(pkin(12));
	t93 = cos(pkin(12));
	t95 = sin(qJ(2));
	t97 = cos(qJ(2));
	t99 = t91 * t95 - t93 * t97;
	t98 = cos(qJ(1));
	t96 = sin(qJ(1));
	t94 = cos(pkin(6));
	t92 = sin(pkin(6));
	t90 = -t97 * t91 - t95 * t93;
	t89 = t99 * t94;
	t1 = [0, t96 * t92, 0, -t96 * t89 - t98 * t90, 0, 0; 0, -t98 * t92, 0, t98 * t89 - t96 * t90, 0, 0; 1, t94, 0, t99 * t92, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (18->10), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->17)
	t131 = sin(pkin(6));
	t136 = sin(qJ(1));
	t143 = t136 * t131;
	t139 = cos(qJ(1));
	t142 = t139 * t131;
	t130 = sin(pkin(12));
	t132 = cos(pkin(12));
	t135 = sin(qJ(2));
	t138 = cos(qJ(2));
	t141 = t138 * t130 + t135 * t132;
	t140 = t135 * t130 - t138 * t132;
	t137 = cos(qJ(4));
	t134 = sin(qJ(4));
	t133 = cos(pkin(6));
	t127 = t141 * t133;
	t126 = t140 * t133;
	t1 = [0, t143, 0, -t136 * t126 + t139 * t141, (-t136 * t127 - t139 * t140) * t134 - t137 * t143, 0; 0, -t142, 0, t139 * t126 + t136 * t141, (t139 * t127 - t136 * t140) * t134 + t137 * t142, 0; 1, t133, 0, t140 * t131, t141 * t134 * t131 - t133 * t137, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:45
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (28->10), mult. (78->24), div. (0->0), fcn. (117->10), ass. (0->20)
	t146 = sin(pkin(6));
	t151 = sin(qJ(1));
	t158 = t151 * t146;
	t154 = cos(qJ(1));
	t157 = t154 * t146;
	t145 = sin(pkin(12));
	t147 = cos(pkin(12));
	t150 = sin(qJ(2));
	t153 = cos(qJ(2));
	t156 = t153 * t145 + t150 * t147;
	t155 = t150 * t145 - t153 * t147;
	t152 = cos(qJ(4));
	t149 = sin(qJ(4));
	t148 = cos(pkin(6));
	t142 = t156 * t148;
	t141 = t155 * t148;
	t140 = t156 * t149 * t146 - t148 * t152;
	t139 = (-t151 * t142 - t154 * t155) * t149 - t152 * t158;
	t138 = (t154 * t142 - t151 * t155) * t149 + t152 * t157;
	t1 = [0, t158, 0, -t151 * t141 + t154 * t156, t139, t139; 0, -t157, 0, t154 * t141 + t151 * t156, t138, t138; 1, t148, 0, t155 * t146, t140, t140;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end