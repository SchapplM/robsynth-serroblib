% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR9
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
%   Wie in S6RRRRRR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:32
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR9_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR9_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:32:11
	% EndTime: 2019-10-10 13:32:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:32:11
	% EndTime: 2019-10-10 13:32:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:32:11
	% EndTime: 2019-10-10 13:32:11
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
	% StartTime: 2019-10-10 13:32:11
	% EndTime: 2019-10-10 13:32:11
	% DurationCPUTime: 0.28s
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
	% StartTime: 2019-10-10 13:32:11
	% EndTime: 2019-10-10 13:32:11
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
	% StartTime: 2019-10-10 13:32:11
	% EndTime: 2019-10-10 13:32:12
	% DurationCPUTime: 0.87s
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
	% StartTime: 2019-10-10 13:32:11
	% EndTime: 2019-10-10 13:32:12
	% DurationCPUTime: 0.83s
	% Computational Cost: add. (2101->66), mult. (5738->154), div. (120->9), fcn. (7872->17), ass. (0->82)
	t156 = cos(pkin(6));
	t162 = cos(qJ(2));
	t193 = sin(qJ(1));
	t172 = t193 * t162;
	t159 = sin(qJ(2));
	t163 = cos(qJ(1));
	t176 = t163 * t159;
	t146 = t156 * t176 + t172;
	t158 = sin(qJ(3));
	t161 = cos(qJ(3));
	t173 = t193 * t159;
	t175 = t163 * t162;
	t145 = -t156 * t175 + t173;
	t153 = sin(pkin(7));
	t155 = cos(pkin(7));
	t154 = sin(pkin(6));
	t180 = t154 * t163;
	t168 = t145 * t155 + t153 * t180;
	t132 = -t146 * t161 + t168 * t158;
	t142 = -t145 * t153 + t155 * t180;
	t157 = sin(qJ(4));
	t160 = cos(qJ(4));
	t194 = t132 * t157 - t142 * t160;
	t120 = t132 * t160 + t142 * t157;
	t179 = t155 * t158;
	t182 = t153 * t156;
	t141 = t158 * t182 + (t159 * t161 + t162 * t179) * t154;
	t144 = -t154 * t162 * t153 + t156 * t155;
	t127 = t141 * t157 - t144 * t160;
	t116 = atan2(t194, t127);
	t113 = sin(t116);
	t114 = cos(t116);
	t107 = t113 * t194 + t114 * t127;
	t106 = 0.1e1 / t107 ^ 2;
	t167 = t156 * t172 + t176;
	t174 = t154 * t193;
	t170 = t153 * t174;
	t147 = -t156 * t173 + t175;
	t183 = t147 * t161;
	t134 = t183 + (-t167 * t155 + t170) * t158;
	t164 = t167 * t153 + t155 * t174;
	t121 = t134 * t157 - t164 * t160;
	t192 = t106 * t121;
	t122 = t134 * t160 + t164 * t157;
	t166 = t167 * t161;
	t133 = t147 * t158 + t155 * t166 - t161 * t170;
	t152 = qJ(5) + qJ(6);
	t150 = sin(t152);
	t151 = cos(t152);
	t112 = t122 * t151 + t133 * t150;
	t110 = 0.1e1 / t112 ^ 2;
	t111 = t122 * t150 - t133 * t151;
	t191 = t110 * t111;
	t190 = t114 * t194;
	t126 = 0.1e1 / t127 ^ 2;
	t189 = t194 * t126;
	t188 = t121 ^ 2 * t106;
	t187 = t133 * t160;
	t181 = t153 * t160;
	t178 = t158 * t159;
	t177 = t161 * t162;
	t171 = t111 ^ 2 * t110 + 0.1e1;
	t169 = -t113 * t127 + t190;
	t165 = -t146 * t158 - t168 * t161;
	t140 = t161 * t182 + (t155 * t177 - t178) * t154;
	t137 = ((-t155 * t178 + t177) * t157 - t159 * t181) * t154;
	t136 = -t147 * t179 - t166;
	t135 = t155 * t183 - t167 * t158;
	t128 = t141 * t160 + t144 * t157;
	t125 = 0.1e1 / t127;
	t124 = t147 * t153 * t157 + t136 * t160;
	t123 = (-t145 * t161 - t146 * t179) * t157 - t146 * t181;
	t115 = 0.1e1 / (t126 * t194 ^ 2 + 0.1e1);
	t109 = 0.1e1 / t112;
	t108 = 0.1e1 / t171;
	t105 = 0.1e1 / t107;
	t104 = 0.1e1 / (0.1e1 + t188);
	t103 = (-t125 * t165 - t140 * t189) * t157 * t115;
	t102 = (-t123 * t125 - t137 * t189) * t115;
	t101 = (t120 * t125 - t128 * t189) * t115;
	t100 = t171 * t108;
	t1 = [-t121 * t125 * t115, t102, t103, t101, 0, 0; (t194 * t105 - (-t113 + (-t125 * t190 + t113) * t115) * t188) * t104, ((t136 * t157 - t147 * t181) * t105 - (t169 * t102 - t113 * t123 + t114 * t137) * t192) * t104, (-t133 * t157 * t105 - ((-t113 * t165 + t114 * t140) * t157 + t169 * t103) * t192) * t104, (t122 * t105 - (t169 * t101 + t113 * t120 + t114 * t128) * t192) * t104, 0, 0; ((t120 * t150 - t151 * t165) * t109 - (t120 * t151 + t150 * t165) * t191) * t108, ((t124 * t150 - t135 * t151) * t109 - (t124 * t151 + t135 * t150) * t191) * t108, ((-t134 * t151 - t150 * t187) * t109 - (t134 * t150 - t151 * t187) * t191) * t108, (-t150 * t109 + t151 * t191) * t121 * t108, t100, t100;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end