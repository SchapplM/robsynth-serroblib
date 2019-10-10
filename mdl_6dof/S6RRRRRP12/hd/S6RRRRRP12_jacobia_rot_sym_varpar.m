% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRP12
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
%   Wie in S6RRRRRP12_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:16
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRP12_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP12_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:15:47
	% EndTime: 2019-10-10 13:15:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:15:47
	% EndTime: 2019-10-10 13:15:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:15:47
	% EndTime: 2019-10-10 13:15:47
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
	% StartTime: 2019-10-10 13:15:47
	% EndTime: 2019-10-10 13:15:47
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
	% StartTime: 2019-10-10 13:15:47
	% EndTime: 2019-10-10 13:15:48
	% DurationCPUTime: 0.43s
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
	% StartTime: 2019-10-10 13:15:48
	% EndTime: 2019-10-10 13:15:48
	% DurationCPUTime: 0.86s
	% Computational Cost: add. (1916->67), mult. (5494->157), div. (115->9), fcn. (7543->17), ass. (0->78)
	t141 = cos(pkin(6));
	t149 = cos(qJ(2));
	t178 = sin(qJ(1));
	t157 = t178 * t149;
	t145 = sin(qJ(2));
	t150 = cos(qJ(1));
	t162 = t150 * t145;
	t134 = t141 * t162 + t157;
	t144 = sin(qJ(3));
	t148 = cos(qJ(3));
	t158 = t178 * t145;
	t161 = t150 * t149;
	t133 = -t141 * t161 + t158;
	t138 = sin(pkin(7));
	t140 = cos(pkin(7));
	t139 = sin(pkin(6));
	t166 = t139 * t150;
	t155 = t133 * t140 + t138 * t166;
	t120 = -t134 * t148 + t155 * t144;
	t130 = -t133 * t138 + t140 * t166;
	t143 = sin(qJ(4));
	t147 = cos(qJ(4));
	t179 = t120 * t143 - t130 * t147;
	t108 = t120 * t147 + t130 * t143;
	t154 = t141 * t157 + t162;
	t159 = t139 * t178;
	t156 = t138 * t159;
	t135 = -t141 * t158 + t161;
	t169 = t135 * t148;
	t122 = t169 + (-t154 * t140 + t156) * t144;
	t151 = t154 * t138 + t140 * t159;
	t110 = t122 * t147 + t151 * t143;
	t153 = t154 * t148;
	t121 = t135 * t144 + t140 * t153 - t148 * t156;
	t142 = sin(qJ(5));
	t146 = cos(qJ(5));
	t100 = t110 * t146 + t121 * t142;
	t98 = 0.1e1 / t100 ^ 2;
	t99 = t110 * t142 - t121 * t146;
	t177 = t98 * t99;
	t109 = t122 * t143 - t151 * t147;
	t165 = t140 * t144;
	t168 = t138 * t141;
	t129 = t144 * t168 + (t145 * t148 + t149 * t165) * t139;
	t132 = -t139 * t149 * t138 + t141 * t140;
	t115 = t129 * t143 - t132 * t147;
	t104 = atan2(t179, t115);
	t101 = sin(t104);
	t102 = cos(t104);
	t95 = t101 * t179 + t102 * t115;
	t94 = 0.1e1 / t95 ^ 2;
	t176 = t109 * t94;
	t175 = t109 ^ 2 * t94;
	t114 = 0.1e1 / t115 ^ 2;
	t174 = t179 * t114;
	t173 = t121 * t147;
	t167 = t138 * t147;
	t164 = t144 * t145;
	t163 = t148 * t149;
	t160 = t99 ^ 2 * t98 + 0.1e1;
	t152 = -t134 * t144 - t155 * t148;
	t128 = t148 * t168 + (t140 * t163 - t164) * t139;
	t125 = ((-t140 * t164 + t163) * t143 - t145 * t167) * t139;
	t124 = -t135 * t165 - t153;
	t123 = t140 * t169 - t154 * t144;
	t116 = t129 * t147 + t132 * t143;
	t113 = 0.1e1 / t115;
	t112 = t135 * t138 * t143 + t124 * t147;
	t111 = (-t133 * t148 - t134 * t165) * t143 - t134 * t167;
	t103 = 0.1e1 / (t114 * t179 ^ 2 + 0.1e1);
	t97 = 0.1e1 / t100;
	t96 = 0.1e1 / t160;
	t93 = 0.1e1 / t95;
	t92 = 0.1e1 / (0.1e1 + t175);
	t91 = (-t113 * t152 - t128 * t174) * t143 * t103;
	t90 = (-t111 * t113 - t125 * t174) * t103;
	t89 = (t108 * t113 - t116 * t174) * t103;
	t1 = [-t109 * t113 * t103, t90, t91, t89, 0, 0; (t179 * t93 - (-t101 + (-t102 * t113 * t179 + t101) * t103) * t175) * t92, ((t124 * t143 - t135 * t167) * t93 - ((t179 * t90 + t125) * t102 + (-t115 * t90 - t111) * t101) * t176) * t92, (-t121 * t143 * t93 - ((t128 * t143 + t179 * t91) * t102 + (-t115 * t91 - t143 * t152) * t101) * t176) * t92, (t110 * t93 - ((t179 * t89 + t116) * t102 + (-t115 * t89 + t108) * t101) * t176) * t92, 0, 0; ((t108 * t142 - t146 * t152) * t97 - (t108 * t146 + t142 * t152) * t177) * t96, ((t112 * t142 - t123 * t146) * t97 - (t112 * t146 + t123 * t142) * t177) * t96, ((-t122 * t146 - t142 * t173) * t97 - (t122 * t142 - t146 * t173) * t177) * t96, (-t142 * t97 + t146 * t177) * t96 * t109, t160 * t96, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:15:48
	% EndTime: 2019-10-10 13:15:49
	% DurationCPUTime: 1.31s
	% Computational Cost: add. (3414->77), mult. (9812->180), div. (137->9), fcn. (13393->17), ass. (0->89)
	t174 = cos(pkin(6));
	t178 = sin(qJ(2));
	t184 = cos(qJ(1));
	t190 = t184 * t178;
	t179 = sin(qJ(1));
	t183 = cos(qJ(2));
	t192 = t179 * t183;
	t165 = t174 * t190 + t192;
	t177 = sin(qJ(3));
	t182 = cos(qJ(3));
	t189 = t184 * t183;
	t193 = t179 * t178;
	t164 = -t174 * t189 + t193;
	t171 = sin(pkin(7));
	t173 = cos(pkin(7));
	t172 = sin(pkin(6));
	t200 = t172 * t184;
	t186 = t164 * t173 + t171 * t200;
	t152 = -t165 * t182 + t186 * t177;
	t159 = -t164 * t171 + t173 * t200;
	t176 = sin(qJ(4));
	t181 = cos(qJ(4));
	t140 = t152 * t181 + t159 * t176;
	t175 = sin(qJ(5));
	t180 = cos(qJ(5));
	t185 = t165 * t177 + t186 * t182;
	t218 = t140 * t175 + t185 * t180;
	t217 = t140 * t180 - t185 * t175;
	t138 = t152 * t176 - t159 * t181;
	t194 = t178 * t182;
	t195 = t177 * t183;
	t203 = t171 * t174;
	t158 = t177 * t203 + (t173 * t195 + t194) * t172;
	t163 = -t172 * t183 * t171 + t174 * t173;
	t149 = t158 * t181 + t163 * t176;
	t191 = t182 * t183;
	t196 = t177 * t178;
	t157 = -t182 * t203 + (-t173 * t191 + t196) * t172;
	t135 = t149 * t175 - t157 * t180;
	t122 = atan2(t218, t135);
	t119 = sin(t122);
	t120 = cos(t122);
	t118 = t119 * t218 + t120 * t135;
	t117 = 0.1e1 / t118 ^ 2;
	t166 = -t174 * t192 - t190;
	t167 = -t174 * t193 + t189;
	t201 = t172 * t179;
	t188 = t171 * t201;
	t154 = t167 * t182 + (t166 * t173 + t188) * t177;
	t161 = -t166 * t171 + t173 * t201;
	t142 = t154 * t181 + t161 * t176;
	t198 = t173 * t182;
	t153 = -t166 * t198 + t167 * t177 - t182 * t188;
	t206 = t153 * t180;
	t129 = t142 * t175 - t206;
	t212 = t117 * t129;
	t211 = t120 * t218;
	t130 = t142 * t180 + t153 * t175;
	t125 = 0.1e1 / t130 ^ 2;
	t141 = -t154 * t176 + t161 * t181;
	t210 = t125 * t141;
	t134 = 0.1e1 / t135 ^ 2;
	t209 = t218 * t134;
	t208 = t129 ^ 2 * t117;
	t207 = t141 ^ 2 * t125;
	t202 = t171 * t176;
	t199 = t173 * t177;
	t197 = t175 * t181;
	t187 = -t119 * t135 + t211;
	t156 = t166 * t182 - t167 * t199;
	t155 = t166 * t177 + t167 * t198;
	t148 = -t158 * t176 + t163 * t181;
	t145 = (((-t173 * t196 + t191) * t181 + t178 * t202) * t175 - (t173 * t194 + t195) * t180) * t172;
	t144 = t156 * t181 + t167 * t202;
	t143 = -t157 * t197 - t158 * t180;
	t136 = t149 * t180 + t157 * t175;
	t133 = 0.1e1 / t135;
	t132 = ((-t164 * t182 - t165 * t199) * t181 + t165 * t202) * t175 - (-t164 * t177 + t165 * t198) * t180;
	t131 = t152 * t180 - t185 * t197;
	t124 = 0.1e1 / t130;
	t123 = 0.1e1 / (0.1e1 + t207);
	t121 = 0.1e1 / (t134 * t218 ^ 2 + 0.1e1);
	t116 = 0.1e1 / t118;
	t115 = 0.1e1 / (0.1e1 + t208);
	t114 = (-t133 * t138 - t148 * t209) * t175 * t121;
	t113 = (-t132 * t133 - t145 * t209) * t121;
	t112 = (-t131 * t133 - t143 * t209) * t121;
	t111 = (t133 * t217 - t136 * t209) * t121;
	t1 = [-t129 * t133 * t121, t113, t112, t114, t111, 0; (t218 * t116 - (-t119 + (-t133 * t211 + t119) * t121) * t208) * t115, ((t144 * t175 - t155 * t180) * t116 - (t187 * t113 - t119 * t132 + t120 * t145) * t212) * t115, ((-t153 * t197 - t154 * t180) * t116 - (t187 * t112 - t119 * t131 + t120 * t143) * t212) * t115, (t141 * t175 * t116 - ((-t119 * t138 + t120 * t148) * t175 + t187 * t114) * t212) * t115, (t130 * t116 - (t187 * t111 + t119 * t217 + t120 * t136) * t212) * t115, 0; (-t138 * t124 - t217 * t210) * t123, ((t167 * t171 * t181 - t156 * t176) * t124 - (t144 * t180 + t155 * t175) * t210) * t123, (t153 * t176 * t124 - (t154 * t175 - t181 * t206) * t210) * t123, (-t124 * t142 - t180 * t207) * t123, t129 * t123 * t210, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end