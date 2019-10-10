% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR12
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
%   Wie in S6RRRRPR12_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:50
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR12_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR12_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
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
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
	% DurationCPUTime: 0.22s
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
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:59
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
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:59
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (899->52), mult. (2499->114), div. (85->9), fcn. (3431->15), ass. (0->69)
	t144 = sin(qJ(1));
	t114 = cos(pkin(6));
	t116 = sin(qJ(2));
	t125 = t144 * t116;
	t118 = cos(qJ(2));
	t119 = cos(qJ(1));
	t129 = t119 * t118;
	t102 = -t114 * t129 + t125;
	t117 = cos(qJ(3));
	t111 = sin(pkin(7));
	t112 = sin(pkin(6));
	t135 = t112 * t119;
	t127 = t111 * t135;
	t113 = cos(pkin(7));
	t133 = t113 * t117;
	t124 = t144 * t118;
	t130 = t119 * t116;
	t103 = t114 * t130 + t124;
	t115 = sin(qJ(3));
	t138 = t103 * t115;
	t85 = t102 * t133 + t117 * t127 + t138;
	t136 = t111 * t114;
	t95 = -t117 * t136 + (t115 * t116 - t118 * t133) * t112;
	t84 = atan2(-t85, t95);
	t81 = sin(t84);
	t82 = cos(t84);
	t75 = -t81 * t85 + t82 * t95;
	t74 = 0.1e1 / t75 ^ 2;
	t104 = -t114 * t125 + t129;
	t121 = t114 * t124 + t130;
	t120 = t121 * t117;
	t126 = t112 * t144;
	t123 = t111 * t126;
	t89 = t104 * t115 + t113 * t120 - t117 * t123;
	t143 = t74 * t89;
	t110 = qJ(4) + pkin(13);
	t108 = sin(t110);
	t109 = cos(t110);
	t90 = t104 * t117 + (-t121 * t113 + t123) * t115;
	t98 = t121 * t111 + t113 * t126;
	t80 = t98 * t108 + t90 * t109;
	t78 = 0.1e1 / t80 ^ 2;
	t79 = t90 * t108 - t98 * t109;
	t142 = t78 * t79;
	t141 = t82 * t85;
	t94 = 0.1e1 / t95 ^ 2;
	t140 = t85 * t94;
	t139 = t89 ^ 2 * t74;
	t137 = t104 * t111;
	t134 = t113 * t115;
	t132 = t115 * t118;
	t131 = t116 * t117;
	t128 = t79 ^ 2 * t78 + 0.1e1;
	t122 = -t81 * t95 - t141;
	t88 = t102 * t134 - t103 * t117 + t115 * t127;
	t101 = (t113 * t131 + t132) * t112;
	t97 = -t102 * t111 + t113 * t135;
	t96 = t115 * t136 + (t113 * t132 + t131) * t112;
	t93 = 0.1e1 / t95;
	t92 = -t104 * t134 - t120;
	t91 = -t102 * t115 + t103 * t133;
	t83 = 0.1e1 / (t85 ^ 2 * t94 + 0.1e1);
	t77 = 0.1e1 / t80;
	t76 = 0.1e1 / t128;
	t73 = 0.1e1 / t75;
	t72 = 0.1e1 / (0.1e1 + t139);
	t71 = (t101 * t140 - t91 * t93) * t83;
	t70 = (t96 * t140 + t88 * t93) * t83;
	t1 = [-t89 * t93 * t83, t71, t70, 0, 0, 0; ((-t138 + (-t102 * t113 - t127) * t117) * t73 - (-t81 + (t93 * t141 + t81) * t83) * t139) * t72, ((t104 * t133 - t121 * t115) * t73 - (t82 * t101 + t122 * t71 - t81 * t91) * t143) * t72, (t90 * t73 - (t122 * t70 + t81 * t88 + t82 * t96) * t143) * t72, 0, 0, 0; ((t88 * t108 - t97 * t109) * t77 - (t97 * t108 + t88 * t109) * t142) * t76, ((t92 * t108 - t109 * t137) * t77 - (t108 * t137 + t92 * t109) * t142) * t76, (-t108 * t77 + t109 * t142) * t89 * t76, t128 * t76, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:59
	% EndTime: 2019-10-10 12:50:00
	% DurationCPUTime: 0.86s
	% Computational Cost: add. (2487->68), mult. (5494->156), div. (115->9), fcn. (7543->17), ass. (0->80)
	t151 = cos(pkin(6));
	t157 = cos(qJ(2));
	t187 = sin(qJ(1));
	t166 = t187 * t157;
	t154 = sin(qJ(2));
	t158 = cos(qJ(1));
	t169 = t158 * t154;
	t141 = t151 * t169 + t166;
	t153 = sin(qJ(3));
	t156 = cos(qJ(3));
	t167 = t187 * t154;
	t170 = t157 * t158;
	t140 = -t151 * t170 + t167;
	t148 = sin(pkin(7));
	t150 = cos(pkin(7));
	t149 = sin(pkin(6));
	t174 = t149 * t158;
	t163 = t140 * t150 + t148 * t174;
	t127 = -t141 * t156 + t163 * t153;
	t137 = -t140 * t148 + t150 * t174;
	t147 = qJ(4) + pkin(13);
	t145 = sin(t147);
	t146 = cos(t147);
	t188 = t127 * t145 - t137 * t146;
	t115 = t127 * t146 + t137 * t145;
	t173 = t150 * t153;
	t175 = t148 * t151;
	t136 = t153 * t175 + (t154 * t156 + t157 * t173) * t149;
	t139 = -t148 * t149 * t157 + t150 * t151;
	t122 = t136 * t145 - t139 * t146;
	t109 = atan2(t188, t122);
	t104 = sin(t109);
	t105 = cos(t109);
	t102 = t104 * t188 + t105 * t122;
	t101 = 0.1e1 / t102 ^ 2;
	t162 = t151 * t166 + t169;
	t168 = t149 * t187;
	t164 = t148 * t168;
	t142 = -t151 * t167 + t170;
	t177 = t142 * t156;
	t129 = t177 + (-t162 * t150 + t164) * t153;
	t159 = t162 * t148 + t150 * t168;
	t116 = t129 * t145 - t159 * t146;
	t186 = t101 * t116;
	t185 = t101 * t116 ^ 2;
	t117 = t129 * t146 + t159 * t145;
	t155 = cos(qJ(6));
	t161 = t162 * t156;
	t128 = t142 * t153 + t150 * t161 - t156 * t164;
	t152 = sin(qJ(6));
	t182 = t128 * t152;
	t111 = t117 * t155 + t182;
	t108 = 0.1e1 / t111 ^ 2;
	t181 = t128 * t155;
	t110 = t117 * t152 - t181;
	t184 = t108 * t110;
	t121 = 0.1e1 / t122 ^ 2;
	t183 = t188 * t121;
	t176 = t146 * t148;
	t172 = t153 * t154;
	t171 = t156 * t157;
	t165 = t108 * t110 ^ 2 + 0.1e1;
	t160 = -t141 * t153 - t163 * t156;
	t135 = t156 * t175 + (t150 * t171 - t172) * t149;
	t132 = -t142 * t173 - t161;
	t131 = t150 * t177 - t162 * t153;
	t130 = ((-t150 * t172 + t171) * t145 - t154 * t176) * t149;
	t123 = t136 * t146 + t139 * t145;
	t120 = 0.1e1 / t122;
	t119 = t142 * t145 * t148 + t132 * t146;
	t118 = (-t140 * t156 - t141 * t173) * t145 - t141 * t176;
	t107 = 0.1e1 / t111;
	t106 = 0.1e1 / (t121 * t188 ^ 2 + 0.1e1);
	t103 = 0.1e1 / t165;
	t100 = 0.1e1 / t102;
	t99 = 0.1e1 / (0.1e1 + t185);
	t98 = (-t120 * t160 - t135 * t183) * t145 * t106;
	t97 = (-t118 * t120 - t130 * t183) * t106;
	t96 = (t115 * t120 - t123 * t183) * t106;
	t1 = [-t116 * t120 * t106, t97, t98, t96, 0, 0; (t188 * t100 - (-t104 + (-t105 * t120 * t188 + t104) * t106) * t185) * t99, ((t132 * t145 - t142 * t176) * t100 - ((t188 * t97 + t130) * t105 + (-t122 * t97 - t118) * t104) * t186) * t99, (-t128 * t145 * t100 - ((t135 * t145 + t188 * t98) * t105 + (-t122 * t98 - t145 * t160) * t104) * t186) * t99, (t117 * t100 - ((t188 * t96 + t123) * t105 + (-t122 * t96 + t115) * t104) * t186) * t99, 0, 0; ((t115 * t152 - t155 * t160) * t107 - (t115 * t155 + t152 * t160) * t184) * t103, ((t119 * t152 - t131 * t155) * t107 - (t119 * t155 + t131 * t152) * t184) * t103, ((-t129 * t155 - t146 * t182) * t107 - (t129 * t152 - t146 * t181) * t184) * t103, (-t107 * t152 + t155 * t184) * t116 * t103, 0, t165 * t103;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end