% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:26
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR13_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR13_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t50 = sin(qJ(2));
	t51 = sin(qJ(1));
	t57 = t51 * t50;
	t52 = cos(qJ(2));
	t56 = t51 * t52;
	t53 = cos(qJ(1));
	t55 = t53 * t50;
	t54 = t53 * t52;
	t49 = cos(pkin(6));
	t48 = sin(pkin(6));
	t47 = -t49 * t57 + t54;
	t46 = -t49 * t56 - t55;
	t45 = -t49 * t55 - t56;
	t44 = -t49 * t54 + t57;
	t1 = [t45, t46, 0, 0, 0, 0; t47, -t44, 0, 0, 0, 0; 0, t48 * t52, 0, 0, 0, 0; t44, -t47, 0, 0, 0, 0; t46, t45, 0, 0, 0, 0; 0, -t48 * t50, 0, 0, 0, 0; t53 * t48, 0, 0, 0, 0, 0; t51 * t48, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (9->7), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t62 = sin(qJ(2));
	t63 = sin(qJ(1));
	t69 = t63 * t62;
	t64 = cos(qJ(2));
	t68 = t63 * t64;
	t65 = cos(qJ(1));
	t67 = t65 * t62;
	t66 = t65 * t64;
	t61 = cos(pkin(6));
	t60 = sin(pkin(6));
	t59 = -t61 * t69 + t66;
	t58 = t61 * t68 + t67;
	t57 = t61 * t67 + t68;
	t56 = -t61 * t66 + t69;
	t1 = [t65 * t60, 0, 0, 0, 0, 0; t63 * t60, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t57, t58, 0, 0, 0, 0; -t59, t56, 0, 0, 0, 0; 0, -t60 * t64, 0, 0, 0, 0; -t56, t59, 0, 0, 0, 0; t58, t57, 0, 0, 0, 0; 0, t60 * t62, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (26->15), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->25)
	t85 = sin(pkin(6));
	t87 = sin(qJ(4));
	t102 = t85 * t87;
	t90 = cos(qJ(4));
	t101 = t85 * t90;
	t91 = cos(qJ(2));
	t100 = t85 * t91;
	t92 = cos(qJ(1));
	t99 = t85 * t92;
	t88 = sin(qJ(2));
	t89 = sin(qJ(1));
	t98 = t89 * t88;
	t97 = t89 * t91;
	t96 = t92 * t88;
	t95 = t92 * t91;
	t86 = cos(pkin(6));
	t79 = -t86 * t95 + t98;
	t94 = -t79 * t87 + t90 * t99;
	t93 = t79 * t90 + t87 * t99;
	t82 = -t86 * t98 + t95;
	t81 = t86 * t97 + t96;
	t80 = t86 * t96 + t97;
	t78 = t101 * t89 + t81 * t87;
	t77 = -t102 * t89 + t81 * t90;
	t1 = [t94, t82 * t87, 0, t77, 0, 0; t78, t80 * t87, 0, t93, 0, 0; 0, t88 * t102, 0, -t100 * t90 - t86 * t87, 0, 0; -t93, t82 * t90, 0, -t78, 0, 0; t77, t80 * t90, 0, t94, 0, 0; 0, t88 * t101, 0, t100 * t87 - t86 * t90, 0, 0; -t80, -t81, 0, 0, 0, 0; t82, -t79, 0, 0, 0, 0; 0, t100, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (57->27), mult. (163->58), div. (0->0), fcn. (238->10), ass. (0->30)
	t112 = sin(pkin(11));
	t116 = sin(qJ(4));
	t131 = t112 * t116;
	t113 = sin(pkin(6));
	t119 = cos(qJ(4));
	t130 = t113 * t119;
	t120 = cos(qJ(2));
	t129 = t113 * t120;
	t121 = cos(qJ(1));
	t128 = t113 * t121;
	t114 = cos(pkin(11));
	t127 = t114 * t116;
	t117 = sin(qJ(2));
	t126 = t116 * t117;
	t118 = sin(qJ(1));
	t125 = t118 * t117;
	t124 = t118 * t120;
	t123 = t121 * t117;
	t122 = t121 * t120;
	t115 = cos(pkin(6));
	t106 = -t115 * t122 + t125;
	t102 = t106 * t119 + t116 * t128;
	t103 = -t106 * t116 + t119 * t128;
	t109 = -t115 * t125 + t122;
	t108 = t115 * t124 + t123;
	t107 = t115 * t123 + t124;
	t105 = -t115 * t116 - t119 * t129;
	t101 = t108 * t116 + t118 * t130;
	t100 = t118 * t113 * t116 - t108 * t119;
	t1 = [t103 * t114 - t107 * t112, -t108 * t112 + t109 * t127, 0, -t100 * t114, 0, 0; t101 * t114 + t109 * t112, -t106 * t112 + t107 * t127, 0, t102 * t114, 0, 0; 0, (t112 * t120 + t114 * t126) * t113, 0, t105 * t114, 0, 0; -t103 * t112 - t107 * t114, -t108 * t114 - t109 * t131, 0, t100 * t112, 0, 0; -t101 * t112 + t109 * t114, -t106 * t114 - t107 * t131, 0, -t102 * t112, 0, 0; 0, (-t112 * t126 + t114 * t120) * t113, 0, -t105 * t112, 0, 0; t102, -t109 * t119, 0, t101, 0, 0; t100, -t107 * t119, 0, -t103, 0, 0; 0, -t117 * t130, 0, t115 * t119 - t116 * t129, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:26:24
	% EndTime: 2019-10-10 10:26:24
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (115->32), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
	t137 = cos(pkin(6));
	t142 = cos(qJ(2));
	t143 = cos(qJ(1));
	t144 = t143 * t142;
	t139 = sin(qJ(2));
	t140 = sin(qJ(1));
	t147 = t140 * t139;
	t127 = -t137 * t144 + t147;
	t138 = sin(qJ(4));
	t141 = cos(qJ(4));
	t136 = sin(pkin(6));
	t149 = t136 * t143;
	t122 = -t127 * t138 + t141 * t149;
	t145 = t143 * t139;
	t146 = t140 * t142;
	t128 = t137 * t145 + t146;
	t135 = pkin(11) + qJ(6);
	t133 = sin(t135);
	t134 = cos(t135);
	t158 = t122 * t133 + t128 * t134;
	t157 = t122 * t134 - t128 * t133;
	t154 = t133 * t138;
	t153 = t134 * t138;
	t152 = t136 * t139;
	t151 = t136 * t141;
	t150 = t136 * t142;
	t148 = t138 * t139;
	t121 = t127 * t141 + t138 * t149;
	t130 = -t137 * t147 + t144;
	t129 = t137 * t146 + t145;
	t126 = t137 * t141 - t138 * t150;
	t125 = -t137 * t138 - t141 * t150;
	t120 = t129 * t138 + t140 * t151;
	t119 = t140 * t136 * t138 - t129 * t141;
	t118 = t120 * t134 + t130 * t133;
	t117 = -t120 * t133 + t130 * t134;
	t1 = [t157, -t129 * t133 + t130 * t153, 0, -t119 * t134, 0, t117; t118, -t127 * t133 + t128 * t153, 0, t121 * t134, 0, t158; 0, (t133 * t142 + t134 * t148) * t136, 0, t125 * t134, 0, -t126 * t133 + t134 * t152; -t158, -t129 * t134 - t130 * t154, 0, t119 * t133, 0, -t118; t117, -t127 * t134 - t128 * t154, 0, -t121 * t133, 0, t157; 0, (-t133 * t148 + t134 * t142) * t136, 0, -t125 * t133, 0, -t126 * t134 - t133 * t152; t121, -t130 * t141, 0, t120, 0, 0; t119, -t128 * t141, 0, -t122, 0, 0; 0, -t139 * t151, 0, t126, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end