% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR8
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
%   Wie in S6RRRRRR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:29
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
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
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:31
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
	% StartTime: 2019-10-10 13:29:31
	% EndTime: 2019-10-10 13:29:31
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
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:31
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (965->52), mult. (2645->114), div. (90->9), fcn. (3627->15), ass. (0->70)
	t151 = sin(qJ(1));
	t119 = sin(pkin(6));
	t122 = sin(qJ(3));
	t123 = sin(qJ(2));
	t124 = cos(qJ(3));
	t125 = cos(qJ(2));
	t120 = cos(pkin(7));
	t140 = t120 * t124;
	t118 = sin(pkin(7));
	t121 = cos(pkin(6));
	t143 = t118 * t121;
	t102 = -t124 * t143 + (t122 * t123 - t125 * t140) * t119;
	t132 = t151 * t123;
	t126 = cos(qJ(1));
	t136 = t126 * t125;
	t109 = -t121 * t136 + t132;
	t142 = t119 * t126;
	t134 = t118 * t142;
	t131 = t151 * t125;
	t137 = t126 * t123;
	t110 = t121 * t137 + t131;
	t145 = t110 * t122;
	t92 = t109 * t140 + t124 * t134 + t145;
	t91 = atan2(-t92, t102);
	t88 = sin(t91);
	t89 = cos(t91);
	t82 = t89 * t102 - t88 * t92;
	t81 = 0.1e1 / t82 ^ 2;
	t111 = -t121 * t132 + t136;
	t128 = t121 * t131 + t137;
	t127 = t128 * t124;
	t133 = t119 * t151;
	t130 = t118 * t133;
	t96 = t111 * t122 + t120 * t127 - t124 * t130;
	t150 = t81 * t96;
	t105 = t128 * t118 + t120 * t133;
	t117 = qJ(4) + qJ(5);
	t115 = sin(t117);
	t116 = cos(t117);
	t97 = t111 * t124 + (-t128 * t120 + t130) * t122;
	t87 = t105 * t115 + t97 * t116;
	t85 = 0.1e1 / t87 ^ 2;
	t86 = -t105 * t116 + t97 * t115;
	t149 = t85 * t86;
	t148 = t89 * t92;
	t147 = t96 ^ 2 * t81;
	t101 = 0.1e1 / t102 ^ 2;
	t146 = t101 * t92;
	t144 = t111 * t118;
	t141 = t120 * t122;
	t139 = t122 * t125;
	t138 = t123 * t124;
	t135 = t86 ^ 2 * t85 + 0.1e1;
	t129 = -t102 * t88 - t148;
	t95 = t109 * t141 - t110 * t124 + t122 * t134;
	t108 = (t120 * t138 + t139) * t119;
	t104 = -t109 * t118 + t120 * t142;
	t103 = t122 * t143 + (t120 * t139 + t138) * t119;
	t100 = 0.1e1 / t102;
	t99 = -t111 * t141 - t127;
	t98 = -t109 * t122 + t110 * t140;
	t90 = 0.1e1 / (t92 ^ 2 * t101 + 0.1e1);
	t84 = 0.1e1 / t87;
	t83 = 0.1e1 / t135;
	t80 = 0.1e1 / t82;
	t79 = 0.1e1 / (0.1e1 + t147);
	t78 = (-t100 * t98 + t108 * t146) * t90;
	t77 = (t100 * t95 + t103 * t146) * t90;
	t76 = t135 * t83;
	t1 = [-t96 * t100 * t90, t78, t77, 0, 0, 0; ((-t145 + (-t109 * t120 - t134) * t124) * t80 - (-t88 + (t100 * t148 + t88) * t90) * t147) * t79, ((t111 * t140 - t128 * t122) * t80 - (t89 * t108 + t129 * t78 - t88 * t98) * t150) * t79, (t97 * t80 - (t89 * t103 + t129 * t77 + t88 * t95) * t150) * t79, 0, 0, 0; ((-t104 * t116 + t95 * t115) * t84 - (t104 * t115 + t95 * t116) * t149) * t83, ((t99 * t115 - t116 * t144) * t84 - (t115 * t144 + t99 * t116) * t149) * t83, (-t115 * t84 + t116 * t149) * t96 * t83, t76, t76, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:31
	% EndTime: 2019-10-10 13:29:32
	% DurationCPUTime: 0.96s
	% Computational Cost: add. (3129->66), mult. (6894->153), div. (145->9), fcn. (9467->17), ass. (0->84)
	t167 = cos(pkin(6));
	t173 = cos(qJ(2));
	t205 = sin(qJ(1));
	t183 = t205 * t173;
	t170 = sin(qJ(2));
	t174 = cos(qJ(1));
	t187 = t174 * t170;
	t157 = t167 * t187 + t183;
	t169 = sin(qJ(3));
	t172 = cos(qJ(3));
	t184 = t205 * t170;
	t186 = t174 * t173;
	t156 = -t167 * t186 + t184;
	t164 = sin(pkin(7));
	t166 = cos(pkin(7));
	t165 = sin(pkin(6));
	t191 = t165 * t174;
	t179 = t156 * t166 + t164 * t191;
	t143 = -t157 * t172 + t179 * t169;
	t153 = -t156 * t164 + t166 * t191;
	t163 = qJ(4) + qJ(5);
	t161 = sin(t163);
	t162 = cos(t163);
	t206 = t143 * t161 - t153 * t162;
	t131 = t143 * t162 + t153 * t161;
	t190 = t166 * t169;
	t192 = t164 * t167;
	t152 = t169 * t192 + (t170 * t172 + t173 * t190) * t165;
	t155 = -t165 * t173 * t164 + t167 * t166;
	t138 = t152 * t161 - t155 * t162;
	t127 = atan2(t206, t138);
	t120 = sin(t127);
	t121 = cos(t127);
	t118 = t120 * t206 + t121 * t138;
	t117 = 0.1e1 / t118 ^ 2;
	t178 = t167 * t183 + t187;
	t185 = t165 * t205;
	t181 = t164 * t185;
	t158 = -t167 * t184 + t186;
	t194 = t158 * t172;
	t145 = t194 + (-t178 * t166 + t181) * t169;
	t175 = t178 * t164 + t166 * t185;
	t132 = t145 * t161 - t175 * t162;
	t204 = t117 * t132;
	t203 = t121 * t206;
	t133 = t145 * t162 + t175 * t161;
	t171 = cos(qJ(6));
	t177 = t178 * t172;
	t144 = t158 * t169 + t166 * t177 - t172 * t181;
	t168 = sin(qJ(6));
	t199 = t144 * t168;
	t126 = t133 * t171 + t199;
	t123 = 0.1e1 / t126 ^ 2;
	t198 = t144 * t171;
	t125 = t133 * t168 - t198;
	t202 = t123 * t125;
	t137 = 0.1e1 / t138 ^ 2;
	t201 = t206 * t137;
	t200 = t132 ^ 2 * t117;
	t193 = t164 * t162;
	t189 = t169 * t170;
	t188 = t172 * t173;
	t182 = t125 ^ 2 * t123 + 0.1e1;
	t180 = -t120 * t138 + t203;
	t176 = -t157 * t169 - t179 * t172;
	t151 = t172 * t192 + (t166 * t188 - t189) * t165;
	t148 = ((-t166 * t189 + t188) * t161 - t170 * t193) * t165;
	t147 = -t158 * t190 - t177;
	t146 = t166 * t194 - t178 * t169;
	t139 = t152 * t162 + t155 * t161;
	t136 = 0.1e1 / t138;
	t135 = t158 * t164 * t161 + t147 * t162;
	t134 = (-t156 * t172 - t157 * t190) * t161 - t157 * t193;
	t124 = 0.1e1 / (t137 * t206 ^ 2 + 0.1e1);
	t122 = 0.1e1 / t126;
	t119 = 0.1e1 / t182;
	t116 = 0.1e1 / t118;
	t115 = 0.1e1 / (0.1e1 + t200);
	t114 = (-t136 * t176 - t151 * t201) * t161 * t124;
	t113 = (-t134 * t136 - t148 * t201) * t124;
	t112 = (t131 * t136 - t139 * t201) * t124;
	t111 = (-t168 * t122 + t171 * t202) * t132 * t119;
	t110 = (t133 * t116 - (t180 * t112 + t120 * t131 + t121 * t139) * t204) * t115;
	t1 = [-t132 * t136 * t124, t113, t114, t112, t112, 0; (t206 * t116 - (-t120 + (-t136 * t203 + t120) * t124) * t200) * t115, ((t147 * t161 - t158 * t193) * t116 - (t180 * t113 - t120 * t134 + t121 * t148) * t204) * t115, (-t144 * t161 * t116 - ((-t120 * t176 + t121 * t151) * t161 + t180 * t114) * t204) * t115, t110, t110, 0; ((t131 * t168 - t171 * t176) * t122 - (t131 * t171 + t168 * t176) * t202) * t119, ((t135 * t168 - t146 * t171) * t122 - (t135 * t171 + t146 * t168) * t202) * t119, ((-t145 * t171 - t162 * t199) * t122 - (t145 * t168 - t162 * t198) * t202) * t119, t111, t111, t182 * t119;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end