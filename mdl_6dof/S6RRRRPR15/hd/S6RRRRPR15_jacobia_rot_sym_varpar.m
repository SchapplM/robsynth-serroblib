% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RRRRPR15_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:56
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR15_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR15_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (35->13), mult. (116->34), div. (25->9), fcn. (175->9), ass. (0->25)
	t30 = cos(pkin(6));
	t29 = sin(pkin(6));
	t34 = cos(qJ(1));
	t38 = t34 * t29;
	t26 = atan2(t38, t30);
	t23 = sin(t26);
	t24 = cos(t26);
	t18 = t23 * t38 + t24 * t30;
	t32 = sin(qJ(1));
	t42 = 0.1e1 / t18 ^ 2 * t32 ^ 2;
	t27 = t29 ^ 2;
	t25 = 0.1e1 / (0.1e1 + t34 ^ 2 * t27 / t30 ^ 2);
	t41 = t25 / t30;
	t31 = sin(qJ(2));
	t40 = t32 * t31;
	t33 = cos(qJ(2));
	t39 = t32 * t33;
	t37 = t34 * t31;
	t36 = t34 * t33;
	t22 = -t30 * t40 + t36;
	t20 = 0.1e1 / t22 ^ 2;
	t21 = t30 * t39 + t37;
	t35 = t21 ^ 2 * t20 + 0.1e1;
	t19 = 0.1e1 / t35;
	t1 = [-t32 * t29 * t41, 0, 0, 0, 0, 0; (0.1e1 / t18 * t38 - (-t24 * t27 * t34 * t41 + (t25 - 0.1e1) * t29 * t23) * t29 * t42) / (t27 * t42 + 0.1e1), 0, 0, 0, 0, 0; ((t30 * t36 - t40) / t22 - (-t30 * t37 - t39) * t21 * t20) * t19, t35 * t19, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (304->29), mult. (885->75), div. (55->9), fcn. (1254->13), ass. (0->49)
	t69 = cos(pkin(6));
	t74 = cos(qJ(2));
	t75 = cos(qJ(1));
	t79 = t74 * t75;
	t71 = sin(qJ(2));
	t72 = sin(qJ(1));
	t81 = t72 * t71;
	t60 = -t69 * t79 + t81;
	t66 = sin(pkin(7));
	t68 = cos(pkin(7));
	t67 = sin(pkin(6));
	t83 = t67 * t75;
	t54 = -t60 * t66 + t68 * t83;
	t59 = -t66 * t67 * t74 + t68 * t69;
	t53 = atan2(t54, t59);
	t50 = sin(t53);
	t51 = cos(t53);
	t44 = t50 * t54 + t51 * t59;
	t43 = 0.1e1 / t44 ^ 2;
	t80 = t72 * t74;
	t82 = t71 * t75;
	t62 = -t69 * t80 - t82;
	t84 = t67 * t72;
	t55 = t62 * t66 - t68 * t84;
	t90 = t43 * t55 ^ 2;
	t70 = sin(qJ(3));
	t76 = t62 * t68 + t66 * t84;
	t63 = -t69 * t81 + t79;
	t73 = cos(qJ(3));
	t86 = t63 * t73;
	t49 = t76 * t70 + t86;
	t47 = 0.1e1 / t49 ^ 2;
	t87 = t63 * t70;
	t48 = -t76 * t73 + t87;
	t89 = t47 * t48;
	t88 = t51 * t54;
	t85 = t67 * t71;
	t78 = t47 * t48 ^ 2 + 0.1e1;
	t77 = t60 * t68 + t66 * t83;
	t61 = -t69 * t82 - t80;
	t58 = 0.1e1 / t59 ^ 2;
	t57 = 0.1e1 / t59;
	t52 = 0.1e1 / (t54 ^ 2 * t58 + 0.1e1);
	t46 = 0.1e1 / t49;
	t45 = 0.1e1 / t78;
	t42 = 0.1e1 / t44;
	t41 = 0.1e1 / (0.1e1 + t90);
	t40 = (-t54 * t58 * t85 + t57 * t61) * t66 * t52;
	t1 = [t55 * t57 * t52, t40, 0, 0, 0, 0; (t54 * t42 + (t50 + (t57 * t88 - t50) * t52) * t90) * t41, (t63 * t66 * t42 + ((t50 * t61 + t51 * t85) * t66 + (-t50 * t59 + t88) * t40) * t55 * t43) * t41, 0, 0, 0, 0; ((t61 * t70 - t77 * t73) * t46 - (t61 * t73 + t77 * t70) * t89) * t45, ((t62 * t70 + t68 * t86) * t46 - (t62 * t73 - t68 * t87) * t89) * t45, t78 * t45, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (833->51), mult. (2499->114), div. (85->9), fcn. (3431->15), ass. (0->68)
	t136 = sin(qJ(1));
	t109 = cos(qJ(3));
	t101 = sin(pkin(7));
	t102 = sin(pkin(6));
	t111 = cos(qJ(1));
	t127 = t102 * t111;
	t119 = t101 * t127;
	t103 = cos(pkin(7));
	t125 = t103 * t109;
	t106 = sin(qJ(3));
	t104 = cos(pkin(6));
	t110 = cos(qJ(2));
	t116 = t136 * t110;
	t107 = sin(qJ(2));
	t122 = t111 * t107;
	t96 = t104 * t122 + t116;
	t129 = t96 * t106;
	t117 = t136 * t107;
	t121 = t111 * t110;
	t95 = -t104 * t121 + t117;
	t78 = t109 * t119 + t95 * t125 + t129;
	t128 = t101 * t104;
	t88 = -t109 * t128 + (t106 * t107 - t110 * t125) * t102;
	t77 = atan2(-t78, t88);
	t74 = sin(t77);
	t75 = cos(t77);
	t68 = -t74 * t78 + t75 * t88;
	t67 = 0.1e1 / t68 ^ 2;
	t113 = t104 * t116 + t122;
	t112 = t113 * t109;
	t118 = t102 * t136;
	t115 = t101 * t118;
	t97 = -t104 * t117 + t121;
	t82 = t103 * t112 + t97 * t106 - t109 * t115;
	t135 = t67 * t82;
	t105 = sin(qJ(4));
	t108 = cos(qJ(4));
	t83 = t97 * t109 + (-t113 * t103 + t115) * t106;
	t91 = t113 * t101 + t103 * t118;
	t73 = t91 * t105 + t83 * t108;
	t71 = 0.1e1 / t73 ^ 2;
	t72 = t83 * t105 - t91 * t108;
	t134 = t71 * t72;
	t133 = t75 * t78;
	t87 = 0.1e1 / t88 ^ 2;
	t132 = t78 * t87;
	t131 = t82 ^ 2 * t67;
	t130 = t101 * t97;
	t126 = t103 * t106;
	t124 = t106 * t110;
	t123 = t107 * t109;
	t120 = t72 ^ 2 * t71 + 0.1e1;
	t114 = -t74 * t88 - t133;
	t81 = t106 * t119 - t96 * t109 + t95 * t126;
	t94 = (t103 * t123 + t124) * t102;
	t90 = -t95 * t101 + t103 * t127;
	t89 = t106 * t128 + (t103 * t124 + t123) * t102;
	t86 = 0.1e1 / t88;
	t85 = -t97 * t126 - t112;
	t84 = -t95 * t106 + t96 * t125;
	t76 = 0.1e1 / (t78 ^ 2 * t87 + 0.1e1);
	t70 = 0.1e1 / t73;
	t69 = 0.1e1 / t120;
	t66 = 0.1e1 / t68;
	t65 = 0.1e1 / (0.1e1 + t131);
	t64 = (t94 * t132 - t84 * t86) * t76;
	t63 = (t89 * t132 + t81 * t86) * t76;
	t1 = [-t82 * t86 * t76, t64, t63, 0, 0, 0; ((-t129 + (-t103 * t95 - t119) * t109) * t66 - (-t74 + (t86 * t133 + t74) * t76) * t131) * t65, ((-t113 * t106 + t97 * t125) * t66 - (t114 * t64 - t74 * t84 + t75 * t94) * t135) * t65, (t83 * t66 - (t114 * t63 + t74 * t81 + t75 * t89) * t135) * t65, 0, 0, 0; ((t81 * t105 - t90 * t108) * t70 - (t90 * t105 + t81 * t108) * t134) * t69, ((t85 * t105 - t108 * t130) * t70 - (t105 * t130 + t85 * t108) * t134) * t69, (-t105 * t70 + t108 * t134) * t82 * t69, t120 * t69, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:27
	% DurationCPUTime: 0.66s
	% Computational Cost: add. (1577->57), mult. (4547->136), div. (107->9), fcn. (6260->15), ass. (0->71)
	t127 = cos(pkin(6));
	t133 = cos(qJ(2));
	t159 = sin(qJ(1));
	t139 = t159 * t133;
	t130 = sin(qJ(2));
	t134 = cos(qJ(1));
	t143 = t134 * t130;
	t119 = t127 * t143 + t139;
	t129 = sin(qJ(3));
	t132 = cos(qJ(3));
	t140 = t159 * t130;
	t142 = t134 * t133;
	t118 = -t127 * t142 + t140;
	t124 = sin(pkin(7));
	t126 = cos(pkin(7));
	t125 = sin(pkin(6));
	t148 = t125 * t134;
	t136 = t118 * t126 + t124 * t148;
	t106 = -t119 * t132 + t136 * t129;
	t115 = -t118 * t124 + t126 * t148;
	t128 = sin(qJ(4));
	t131 = cos(qJ(4));
	t161 = t106 * t128 - t115 * t131;
	t160 = t106 * t131 + t115 * t128;
	t147 = t126 * t129;
	t150 = t124 * t127;
	t114 = t129 * t150 + (t130 * t132 + t133 * t147) * t125;
	t117 = -t125 * t133 * t124 + t127 * t126;
	t100 = t114 * t128 - t117 * t131;
	t90 = atan2(t161, t100);
	t87 = sin(t90);
	t88 = cos(t90);
	t86 = t100 * t88 + t161 * t87;
	t85 = 0.1e1 / t86 ^ 2;
	t120 = -t127 * t139 - t143;
	t121 = -t127 * t140 + t142;
	t141 = t125 * t159;
	t138 = t124 * t141;
	t108 = t121 * t132 + (t120 * t126 + t138) * t129;
	t135 = -t120 * t124 + t126 * t141;
	t95 = t108 * t128 - t135 * t131;
	t158 = t85 * t95;
	t157 = t85 * t95 ^ 2;
	t156 = t88 * t161;
	t99 = 0.1e1 / t100 ^ 2;
	t155 = t161 * t99;
	t146 = t126 * t132;
	t107 = -t120 * t146 + t121 * t129 - t132 * t138;
	t103 = 0.1e1 / t107 ^ 2;
	t96 = t108 * t131 + t135 * t128;
	t154 = t103 * t96;
	t149 = t124 * t131;
	t145 = t129 * t130;
	t144 = t132 * t133;
	t137 = -t100 * t87 + t156;
	t104 = -t119 * t129 - t136 * t132;
	t113 = t132 * t150 + (t126 * t144 - t145) * t125;
	t110 = ((-t126 * t145 + t144) * t128 - t130 * t149) * t125;
	t109 = t120 * t132 - t121 * t147;
	t102 = 0.1e1 / t107;
	t101 = t114 * t131 + t117 * t128;
	t98 = 0.1e1 / t100;
	t97 = (-t118 * t132 - t119 * t147) * t128 - t119 * t149;
	t91 = 0.1e1 / (t103 * t96 ^ 2 + 0.1e1);
	t89 = 0.1e1 / (t161 ^ 2 * t99 + 0.1e1);
	t84 = 0.1e1 / t86;
	t83 = 0.1e1 / (0.1e1 + t157);
	t82 = (-t104 * t98 - t113 * t155) * t89 * t128;
	t81 = (-t110 * t155 - t97 * t98) * t89;
	t80 = (-t101 * t155 + t160 * t98) * t89;
	t1 = [-t95 * t98 * t89, t81, t82, t80, 0, 0; (t161 * t84 - (-t87 + (-t98 * t156 + t87) * t89) * t157) * t83, ((t109 * t128 - t121 * t149) * t84 - (t110 * t88 + t137 * t81 - t87 * t97) * t158) * t83, (-t107 * t128 * t84 - (t137 * t82 + (-t104 * t87 + t113 * t88) * t128) * t158) * t83, (t96 * t84 - (t101 * t88 + t137 * t80 + t160 * t87) * t158) * t83, 0, 0; (t160 * t102 - t104 * t154) * t91, ((t121 * t124 * t128 + t109 * t131) * t102 - (t120 * t129 + t121 * t146) * t154) * t91, (-t102 * t107 * t131 - t108 * t154) * t91, -t95 * t102 * t91, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:27
	% DurationCPUTime: 0.82s
	% Computational Cost: add. (1916->67), mult. (5494->159), div. (115->9), fcn. (7543->17), ass. (0->76)
	t144 = sin(qJ(2));
	t145 = sin(qJ(1));
	t149 = cos(qJ(2));
	t150 = cos(qJ(1));
	t169 = cos(pkin(6));
	t154 = t150 * t169;
	t132 = t144 * t154 + t145 * t149;
	t143 = sin(qJ(3));
	t148 = cos(qJ(3));
	t131 = t145 * t144 - t149 * t154;
	t138 = sin(pkin(7));
	t140 = cos(pkin(7));
	t139 = sin(pkin(6));
	t162 = t139 * t150;
	t152 = t131 * t140 + t138 * t162;
	t118 = -t132 * t148 + t152 * t143;
	t126 = -t131 * t138 + t140 * t162;
	t142 = sin(qJ(4));
	t147 = cos(qJ(4));
	t106 = t118 * t142 - t126 * t147;
	t175 = t118 * t147 + t126 * t142;
	t155 = t145 * t169;
	t133 = -t150 * t144 - t149 * t155;
	t134 = -t144 * t155 + t150 * t149;
	t163 = t139 * t145;
	t156 = t138 * t163;
	t120 = t134 * t148 + (t133 * t140 + t156) * t143;
	t128 = -t133 * t138 + t140 * t163;
	t107 = t120 * t142 - t128 * t147;
	t141 = sin(qJ(6));
	t160 = t140 * t148;
	t119 = -t133 * t160 + t134 * t143 - t148 * t156;
	t146 = cos(qJ(6));
	t166 = t119 * t146;
	t98 = t107 * t141 + t166;
	t96 = 0.1e1 / t98 ^ 2;
	t167 = t119 * t141;
	t97 = -t107 * t146 + t167;
	t172 = t96 * t97;
	t108 = t120 * t147 + t128 * t142;
	t153 = t169 * t138;
	t161 = t140 * t143;
	t125 = t143 * t153 + (t144 * t148 + t149 * t161) * t139;
	t130 = -t139 * t149 * t138 + t169 * t140;
	t114 = t125 * t147 + t130 * t142;
	t102 = atan2(t175, t114);
	t100 = cos(t102);
	t99 = sin(t102);
	t93 = t100 * t114 + t175 * t99;
	t92 = 0.1e1 / t93 ^ 2;
	t171 = t108 * t92;
	t170 = t108 ^ 2 * t92;
	t112 = 0.1e1 / t114 ^ 2;
	t168 = t175 * t112;
	t164 = t138 * t142;
	t159 = t143 * t144;
	t158 = t148 * t149;
	t157 = t97 ^ 2 * t96 + 0.1e1;
	t151 = -t132 * t143 - t152 * t148;
	t124 = t148 * t153 + (t140 * t158 - t159) * t139;
	t123 = ((-t140 * t159 + t158) * t147 + t144 * t164) * t139;
	t122 = t133 * t148 - t134 * t161;
	t121 = t133 * t143 + t134 * t160;
	t113 = -t125 * t142 + t130 * t147;
	t111 = 0.1e1 / t114;
	t110 = -t134 * t138 * t147 + t122 * t142;
	t109 = (-t131 * t148 - t132 * t161) * t147 + t132 * t164;
	t101 = 0.1e1 / (t112 * t175 ^ 2 + 0.1e1);
	t95 = 0.1e1 / t98;
	t94 = 0.1e1 / t157;
	t91 = 0.1e1 / t93;
	t90 = 0.1e1 / (0.1e1 + t170);
	t89 = (-t111 * t151 - t124 * t168) * t147 * t101;
	t88 = (-t109 * t111 - t123 * t168) * t101;
	t87 = (-t106 * t111 - t113 * t168) * t101;
	t1 = [-t108 * t111 * t101, t88, t89, t87, 0, 0; (t175 * t91 - (-t99 + (-t100 * t111 * t175 + t99) * t101) * t170) * t90, ((t122 * t147 + t134 * t164) * t91 - ((-t114 * t88 - t109) * t99 + (t175 * t88 + t123) * t100) * t171) * t90, (-t119 * t147 * t91 - ((-t114 * t89 - t147 * t151) * t99 + (t124 * t147 + t175 * t89) * t100) * t171) * t90, (-t107 * t91 - ((-t114 * t87 - t106) * t99 + (t175 * t87 + t113) * t100) * t171) * t90, 0, 0; ((-t106 * t146 + t141 * t151) * t95 - (t106 * t141 + t146 * t151) * t172) * t94, ((-t110 * t146 + t121 * t141) * t95 - (t110 * t141 + t121 * t146) * t172) * t94, ((t120 * t141 + t142 * t166) * t95 - (t120 * t146 - t142 * t167) * t172) * t94, (-t141 * t172 - t146 * t95) * t94 * t108, 0, t157 * t94;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end