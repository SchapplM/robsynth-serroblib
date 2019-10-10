% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR13
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
%   Wie in S6RRRPRR13_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:14
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR13_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR13_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
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
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:18
	% DurationCPUTime: 0.23s
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
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:18
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (781->51), mult. (2353->113), div. (80->9), fcn. (3235->15), ass. (0->67)
	t132 = sin(qJ(1));
	t106 = cos(qJ(3));
	t100 = sin(pkin(6));
	t108 = cos(qJ(1));
	t123 = t100 * t108;
	t99 = sin(pkin(7));
	t116 = t99 * t123;
	t102 = cos(pkin(7));
	t121 = t102 * t106;
	t104 = sin(qJ(3));
	t103 = cos(pkin(6));
	t107 = cos(qJ(2));
	t113 = t132 * t107;
	t105 = sin(qJ(2));
	t118 = t108 * t105;
	t93 = t103 * t118 + t113;
	t124 = t93 * t104;
	t114 = t132 * t105;
	t117 = t108 * t107;
	t92 = -t103 * t117 + t114;
	t75 = t106 * t116 + t92 * t121 + t124;
	t125 = t103 * t99;
	t85 = -t106 * t125 + (t104 * t105 - t107 * t121) * t100;
	t74 = atan2(-t75, t85);
	t71 = sin(t74);
	t72 = cos(t74);
	t65 = -t71 * t75 + t72 * t85;
	t64 = 0.1e1 / t65 ^ 2;
	t110 = t103 * t113 + t118;
	t109 = t110 * t106;
	t115 = t100 * t132;
	t112 = t99 * t115;
	t94 = -t103 * t114 + t117;
	t79 = t102 * t109 + t94 * t104 - t106 * t112;
	t131 = t64 * t79;
	t130 = t64 * t79 ^ 2;
	t101 = cos(pkin(13));
	t80 = t94 * t106 + (-t110 * t102 + t112) * t104;
	t88 = t102 * t115 + t110 * t99;
	t98 = sin(pkin(13));
	t70 = t80 * t101 + t88 * t98;
	t68 = 0.1e1 / t70 ^ 2;
	t69 = -t88 * t101 + t80 * t98;
	t129 = t68 * t69;
	t128 = t72 * t75;
	t84 = 0.1e1 / t85 ^ 2;
	t127 = t75 * t84;
	t126 = t94 * t99;
	t122 = t102 * t104;
	t120 = t104 * t107;
	t119 = t105 * t106;
	t111 = -t71 * t85 - t128;
	t78 = t104 * t116 - t93 * t106 + t92 * t122;
	t91 = (t102 * t119 + t120) * t100;
	t87 = t102 * t123 - t92 * t99;
	t86 = t104 * t125 + (t102 * t120 + t119) * t100;
	t83 = 0.1e1 / t85;
	t82 = -t94 * t122 - t109;
	t81 = -t92 * t104 + t93 * t121;
	t73 = 0.1e1 / (t75 ^ 2 * t84 + 0.1e1);
	t67 = 0.1e1 / t70;
	t66 = 0.1e1 / (t68 * t69 ^ 2 + 0.1e1);
	t63 = 0.1e1 / t65;
	t62 = 0.1e1 / (0.1e1 + t130);
	t61 = (t91 * t127 - t81 * t83) * t73;
	t60 = (t86 * t127 + t78 * t83) * t73;
	t1 = [-t79 * t83 * t73, t61, t60, 0, 0, 0; ((-t124 + (-t102 * t92 - t116) * t106) * t63 - (-t71 + (t83 * t128 + t71) * t73) * t130) * t62, ((-t110 * t104 + t94 * t121) * t63 - (t111 * t61 - t71 * t81 + t72 * t91) * t131) * t62, (t80 * t63 - (t111 * t60 + t71 * t78 + t72 * t86) * t131) * t62, 0, 0, 0; ((-t87 * t101 + t78 * t98) * t67 - (t78 * t101 + t87 * t98) * t129) * t66, ((-t101 * t126 + t82 * t98) * t67 - (t82 * t101 + t98 * t126) * t129) * t66, (t101 * t129 - t98 * t67) * t79 * t66, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:18
	% EndTime: 2019-10-10 12:14:18
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (899->52), mult. (2499->114), div. (85->9), fcn. (3431->15), ass. (0->69)
	t143 = sin(qJ(1));
	t113 = cos(pkin(6));
	t115 = sin(qJ(2));
	t124 = t143 * t115;
	t117 = cos(qJ(2));
	t118 = cos(qJ(1));
	t128 = t118 * t117;
	t101 = -t113 * t128 + t124;
	t116 = cos(qJ(3));
	t110 = sin(pkin(7));
	t111 = sin(pkin(6));
	t134 = t111 * t118;
	t126 = t110 * t134;
	t112 = cos(pkin(7));
	t132 = t112 * t116;
	t123 = t143 * t117;
	t129 = t118 * t115;
	t102 = t113 * t129 + t123;
	t114 = sin(qJ(3));
	t137 = t102 * t114;
	t84 = t101 * t132 + t116 * t126 + t137;
	t135 = t110 * t113;
	t94 = -t116 * t135 + (t114 * t115 - t117 * t132) * t111;
	t83 = atan2(-t84, t94);
	t80 = sin(t83);
	t81 = cos(t83);
	t74 = -t80 * t84 + t81 * t94;
	t73 = 0.1e1 / t74 ^ 2;
	t103 = -t113 * t124 + t128;
	t120 = t113 * t123 + t129;
	t119 = t120 * t116;
	t125 = t111 * t143;
	t122 = t110 * t125;
	t88 = t103 * t114 + t112 * t119 - t116 * t122;
	t142 = t73 * t88;
	t109 = pkin(13) + qJ(5);
	t107 = sin(t109);
	t108 = cos(t109);
	t89 = t103 * t116 + (-t120 * t112 + t122) * t114;
	t97 = t120 * t110 + t112 * t125;
	t79 = t97 * t107 + t89 * t108;
	t77 = 0.1e1 / t79 ^ 2;
	t78 = t89 * t107 - t97 * t108;
	t141 = t77 * t78;
	t140 = t81 * t84;
	t93 = 0.1e1 / t94 ^ 2;
	t139 = t84 * t93;
	t138 = t88 ^ 2 * t73;
	t136 = t103 * t110;
	t133 = t112 * t114;
	t131 = t114 * t117;
	t130 = t115 * t116;
	t127 = t78 ^ 2 * t77 + 0.1e1;
	t121 = -t80 * t94 - t140;
	t87 = t101 * t133 - t102 * t116 + t114 * t126;
	t100 = (t112 * t130 + t131) * t111;
	t96 = -t101 * t110 + t112 * t134;
	t95 = t114 * t135 + (t112 * t131 + t130) * t111;
	t92 = 0.1e1 / t94;
	t91 = -t103 * t133 - t119;
	t90 = -t101 * t114 + t102 * t132;
	t82 = 0.1e1 / (t84 ^ 2 * t93 + 0.1e1);
	t76 = 0.1e1 / t79;
	t75 = 0.1e1 / t127;
	t72 = 0.1e1 / t74;
	t71 = 0.1e1 / (0.1e1 + t138);
	t70 = (t100 * t139 - t90 * t92) * t82;
	t69 = (t95 * t139 + t87 * t92) * t82;
	t1 = [-t88 * t92 * t82, t70, t69, 0, 0, 0; ((-t137 + (-t101 * t112 - t126) * t116) * t72 - (-t80 + (t92 * t140 + t80) * t82) * t138) * t71, ((t103 * t132 - t120 * t114) * t72 - (t81 * t100 + t121 * t70 - t80 * t90) * t142) * t71, (t89 * t72 - (t121 * t69 + t80 * t87 + t81 * t95) * t142) * t71, 0, 0, 0; ((t87 * t107 - t96 * t108) * t76 - (t96 * t107 + t87 * t108) * t141) * t75, ((t91 * t107 - t108 * t136) * t76 - (t107 * t136 + t91 * t108) * t141) * t75, (-t107 * t76 + t108 * t141) * t88 * t75, 0, t127 * t75, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:18
	% EndTime: 2019-10-10 12:14:19
	% DurationCPUTime: 0.89s
	% Computational Cost: add. (2487->68), mult. (5494->156), div. (115->9), fcn. (7543->17), ass. (0->80)
	t150 = cos(pkin(6));
	t156 = cos(qJ(2));
	t186 = sin(qJ(1));
	t165 = t186 * t156;
	t153 = sin(qJ(2));
	t157 = cos(qJ(1));
	t169 = t157 * t153;
	t140 = t150 * t169 + t165;
	t152 = sin(qJ(3));
	t155 = cos(qJ(3));
	t166 = t186 * t153;
	t168 = t157 * t156;
	t139 = -t150 * t168 + t166;
	t147 = sin(pkin(7));
	t149 = cos(pkin(7));
	t148 = sin(pkin(6));
	t173 = t148 * t157;
	t162 = t139 * t149 + t147 * t173;
	t126 = -t140 * t155 + t162 * t152;
	t136 = -t139 * t147 + t149 * t173;
	t146 = pkin(13) + qJ(5);
	t144 = sin(t146);
	t145 = cos(t146);
	t187 = t126 * t144 - t136 * t145;
	t114 = t126 * t145 + t136 * t144;
	t172 = t149 * t152;
	t174 = t147 * t150;
	t135 = t152 * t174 + (t153 * t155 + t156 * t172) * t148;
	t138 = -t148 * t156 * t147 + t150 * t149;
	t121 = t135 * t144 - t138 * t145;
	t108 = atan2(t187, t121);
	t103 = sin(t108);
	t104 = cos(t108);
	t101 = t103 * t187 + t104 * t121;
	t100 = 0.1e1 / t101 ^ 2;
	t161 = t150 * t165 + t169;
	t167 = t148 * t186;
	t163 = t147 * t167;
	t141 = -t150 * t166 + t168;
	t176 = t141 * t155;
	t128 = t176 + (-t161 * t149 + t163) * t152;
	t158 = t161 * t147 + t149 * t167;
	t115 = t128 * t144 - t158 * t145;
	t185 = t100 * t115;
	t116 = t128 * t145 + t158 * t144;
	t154 = cos(qJ(6));
	t160 = t161 * t155;
	t127 = t141 * t152 + t149 * t160 - t155 * t163;
	t151 = sin(qJ(6));
	t181 = t127 * t151;
	t110 = t116 * t154 + t181;
	t107 = 0.1e1 / t110 ^ 2;
	t180 = t127 * t154;
	t109 = t116 * t151 - t180;
	t184 = t107 * t109;
	t120 = 0.1e1 / t121 ^ 2;
	t183 = t187 * t120;
	t182 = t115 ^ 2 * t100;
	t175 = t147 * t145;
	t171 = t152 * t153;
	t170 = t155 * t156;
	t164 = t109 ^ 2 * t107 + 0.1e1;
	t159 = -t140 * t152 - t162 * t155;
	t134 = t155 * t174 + (t149 * t170 - t171) * t148;
	t131 = -t141 * t172 - t160;
	t130 = t149 * t176 - t161 * t152;
	t129 = ((-t149 * t171 + t170) * t144 - t153 * t175) * t148;
	t122 = t135 * t145 + t138 * t144;
	t119 = 0.1e1 / t121;
	t118 = t141 * t147 * t144 + t131 * t145;
	t117 = (-t139 * t155 - t140 * t172) * t144 - t140 * t175;
	t106 = 0.1e1 / t110;
	t105 = 0.1e1 / (t120 * t187 ^ 2 + 0.1e1);
	t102 = 0.1e1 / t164;
	t99 = 0.1e1 / t101;
	t98 = 0.1e1 / (0.1e1 + t182);
	t97 = (-t119 * t159 - t134 * t183) * t144 * t105;
	t96 = (-t117 * t119 - t129 * t183) * t105;
	t95 = (t114 * t119 - t122 * t183) * t105;
	t1 = [-t115 * t119 * t105, t96, t97, 0, t95, 0; (t187 * t99 - (-t103 + (-t104 * t119 * t187 + t103) * t105) * t182) * t98, ((t131 * t144 - t141 * t175) * t99 - ((t187 * t96 + t129) * t104 + (-t121 * t96 - t117) * t103) * t185) * t98, (-t127 * t144 * t99 - ((t134 * t144 + t187 * t97) * t104 + (-t121 * t97 - t144 * t159) * t103) * t185) * t98, 0, (t116 * t99 - ((t187 * t95 + t122) * t104 + (-t121 * t95 + t114) * t103) * t185) * t98, 0; ((t114 * t151 - t154 * t159) * t106 - (t114 * t154 + t151 * t159) * t184) * t102, ((t118 * t151 - t130 * t154) * t106 - (t118 * t154 + t130 * t151) * t184) * t102, ((-t128 * t154 - t145 * t181) * t106 - (t128 * t151 - t145 * t180) * t184) * t102, 0, (-t151 * t106 + t154 * t184) * t115 * t102, t164 * t102;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end