% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRPR5_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR5_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t59 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, -cos(qJ(1)) * t59, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->5), mult. (22->12), div. (0->0), fcn. (35->8), ass. (0->12)
	t91 = sin(pkin(11));
	t93 = cos(pkin(11));
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
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (8->5), mult. (22->12), div. (0->0), fcn. (35->8), ass. (0->12)
	t118 = sin(pkin(11));
	t120 = cos(pkin(11));
	t122 = sin(qJ(2));
	t124 = cos(qJ(2));
	t126 = t118 * t122 - t120 * t124;
	t125 = cos(qJ(1));
	t123 = sin(qJ(1));
	t121 = cos(pkin(6));
	t119 = sin(pkin(6));
	t117 = -t124 * t118 - t122 * t120;
	t116 = t126 * t121;
	t1 = [0, t123 * t119, 0, -t123 * t116 - t125 * t117, 0, 0; 0, -t125 * t119, 0, t125 * t116 - t123 * t117, 0, 0; 1, t121, 0, t126 * t119, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:38
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (18->10), mult. (50->24), div. (0->0), fcn. (76->10), ass. (0->17)
	t132 = sin(pkin(6));
	t137 = sin(qJ(1));
	t144 = t137 * t132;
	t140 = cos(qJ(1));
	t143 = t140 * t132;
	t131 = sin(pkin(11));
	t133 = cos(pkin(11));
	t136 = sin(qJ(2));
	t139 = cos(qJ(2));
	t142 = t139 * t131 + t136 * t133;
	t141 = t136 * t131 - t139 * t133;
	t138 = cos(qJ(4));
	t135 = sin(qJ(4));
	t134 = cos(pkin(6));
	t128 = t142 * t134;
	t127 = t141 * t134;
	t1 = [0, t144, 0, -t137 * t127 + t140 * t142, 0, (-t137 * t128 - t140 * t141) * t135 - t138 * t144; 0, -t143, 0, t140 * t127 + t137 * t142, 0, (t140 * t128 - t137 * t141) * t135 + t138 * t143; 1, t134, 0, t141 * t132, 0, t142 * t135 * t132 - t134 * t138;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end