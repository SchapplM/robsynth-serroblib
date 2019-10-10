% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:20
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR10_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR10_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
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
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
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
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (19->11), mult. (57->26), div. (0->0), fcn. (88->8), ass. (0->20)
	t66 = sin(pkin(6));
	t70 = sin(qJ(1));
	t79 = t66 * t70;
	t71 = cos(qJ(2));
	t78 = t66 * t71;
	t72 = cos(qJ(1));
	t77 = t66 * t72;
	t69 = sin(qJ(2));
	t76 = t70 * t69;
	t75 = t70 * t71;
	t74 = t72 * t69;
	t73 = t72 * t71;
	t68 = cos(pkin(6));
	t67 = cos(pkin(11));
	t65 = sin(pkin(11));
	t64 = -t68 * t76 + t73;
	t63 = t68 * t75 + t74;
	t62 = t68 * t74 + t75;
	t61 = t68 * t73 - t76;
	t1 = [-t62 * t67 + t65 * t77, -t63 * t67, 0, 0, 0, 0; t64 * t67 + t65 * t79, t61 * t67, 0, 0, 0, 0; 0, t67 * t78, 0, 0, 0, 0; t62 * t65 + t67 * t77, t63 * t65, 0, 0, 0, 0; -t64 * t65 + t67 * t79, -t61 * t65, 0, 0, 0, 0; 0, -t65 * t78, 0, 0, 0, 0; t61, t64, 0, 0, 0, 0; t63, t62, 0, 0, 0, 0; 0, t66 * t69, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (55->16), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->26)
	t90 = sin(pkin(6));
	t92 = sin(qJ(2));
	t105 = t90 * t92;
	t93 = sin(qJ(1));
	t104 = t90 * t93;
	t94 = cos(qJ(2));
	t103 = t90 * t94;
	t95 = cos(qJ(1));
	t102 = t90 * t95;
	t101 = t93 * t92;
	t100 = t93 * t94;
	t99 = t95 * t92;
	t98 = t95 * t94;
	t91 = cos(pkin(6));
	t83 = t91 * t99 + t100;
	t89 = pkin(11) + qJ(4);
	t87 = sin(t89);
	t88 = cos(t89);
	t97 = t87 * t102 - t83 * t88;
	t96 = t88 * t102 + t83 * t87;
	t85 = -t91 * t101 + t98;
	t84 = t91 * t100 + t99;
	t82 = t91 * t98 - t101;
	t81 = t87 * t104 + t85 * t88;
	t80 = t88 * t104 - t85 * t87;
	t1 = [t97, -t84 * t88, 0, t80, 0, 0; t81, t82 * t88, 0, -t96, 0, 0; 0, t88 * t103, 0, -t87 * t105 + t91 * t88, 0, 0; t96, t84 * t87, 0, -t81, 0, 0; t80, -t82 * t87, 0, t97, 0, 0; 0, -t87 * t103, 0, -t88 * t105 - t91 * t87, 0, 0; t82, t85, 0, 0, 0, 0; t84, t83, 0, 0, 0, 0; 0, t105, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (55->16), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->26)
	t106 = sin(pkin(6));
	t108 = sin(qJ(2));
	t121 = t106 * t108;
	t109 = sin(qJ(1));
	t120 = t106 * t109;
	t110 = cos(qJ(2));
	t119 = t106 * t110;
	t111 = cos(qJ(1));
	t118 = t106 * t111;
	t117 = t109 * t108;
	t116 = t109 * t110;
	t115 = t111 * t108;
	t114 = t111 * t110;
	t107 = cos(pkin(6));
	t100 = t107 * t115 + t116;
	t105 = pkin(11) + qJ(4);
	t103 = sin(t105);
	t104 = cos(t105);
	t113 = t100 * t103 + t104 * t118;
	t112 = t100 * t104 - t103 * t118;
	t102 = -t107 * t117 + t114;
	t101 = t107 * t116 + t115;
	t99 = t107 * t114 - t117;
	t98 = t102 * t104 + t103 * t120;
	t97 = t102 * t103 - t104 * t120;
	t1 = [t99, t102, 0, 0, 0, 0; t101, t100, 0, 0, 0, 0; 0, t121, 0, 0, 0, 0; t112, t101 * t104, 0, t97, 0, 0; -t98, -t99 * t104, 0, t113, 0, 0; 0, -t104 * t119, 0, t103 * t121 - t107 * t104, 0, 0; -t113, -t101 * t103, 0, t98, 0, 0; t97, t99 * t103, 0, t112, 0, 0; 0, t103 * t119, 0, t107 * t103 + t104 * t121, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:53
	% EndTime: 2019-10-10 10:20:53
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (122->32), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
	t137 = cos(pkin(6));
	t139 = sin(qJ(2));
	t143 = cos(qJ(1));
	t146 = t143 * t139;
	t140 = sin(qJ(1));
	t142 = cos(qJ(2));
	t148 = t140 * t142;
	t129 = t137 * t146 + t148;
	t135 = pkin(11) + qJ(4);
	t133 = sin(t135);
	t134 = cos(t135);
	t136 = sin(pkin(6));
	t151 = t136 * t143;
	t121 = t129 * t133 + t134 * t151;
	t145 = t143 * t142;
	t149 = t140 * t139;
	t128 = -t137 * t145 + t149;
	t138 = sin(qJ(6));
	t141 = cos(qJ(6));
	t159 = -t121 * t138 - t128 * t141;
	t158 = t121 * t141 - t128 * t138;
	t155 = t133 * t138;
	t154 = t133 * t141;
	t153 = t136 * t139;
	t152 = t136 * t140;
	t150 = t138 * t142;
	t147 = t141 * t142;
	t144 = -t129 * t134 + t133 * t151;
	t131 = -t137 * t149 + t145;
	t130 = t137 * t148 + t146;
	t127 = t137 * t133 + t134 * t153;
	t126 = t133 * t153 - t137 * t134;
	t125 = t131 * t134 + t133 * t152;
	t124 = t131 * t133 - t134 * t152;
	t120 = t124 * t138 + t130 * t141;
	t119 = t124 * t141 - t130 * t138;
	t1 = [t159, -t130 * t155 + t131 * t141, 0, t125 * t138, 0, t119; t120, -t128 * t155 + t129 * t141, 0, -t144 * t138, 0, t158; 0, (t133 * t150 + t139 * t141) * t136, 0, t127 * t138, 0, t126 * t141 + t136 * t150; -t158, -t130 * t154 - t131 * t138, 0, t125 * t141, 0, -t120; t119, -t128 * t154 - t129 * t138, 0, -t144 * t141, 0, t159; 0, (t133 * t147 - t138 * t139) * t136, 0, t127 * t141, 0, -t126 * t138 + t136 * t147; t144, -t130 * t134, 0, -t124, 0, 0; t125, -t128 * t134, 0, -t121, 0, 0; 0, t136 * t142 * t134, 0, -t126, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end