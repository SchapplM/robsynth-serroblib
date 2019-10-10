% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR15
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
%   Wie in S6RRRPRR15_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR15_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR15_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
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
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
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
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (662->45), mult. (2020->99), div. (77->9), fcn. (2784->13), ass. (0->59)
	t90 = sin(pkin(6));
	t98 = cos(qJ(1));
	t111 = t90 * t98;
	t89 = sin(pkin(7));
	t101 = t89 * t111;
	t91 = cos(pkin(7));
	t96 = cos(qJ(3));
	t109 = t91 * t96;
	t94 = sin(qJ(2));
	t104 = t98 * t94;
	t95 = sin(qJ(1));
	t97 = cos(qJ(2));
	t105 = t95 * t97;
	t92 = cos(pkin(6));
	t83 = t92 * t104 + t105;
	t93 = sin(qJ(3));
	t114 = t83 * t93;
	t103 = t98 * t97;
	t106 = t95 * t94;
	t82 = -t92 * t103 + t106;
	t66 = t96 * t101 + t82 * t109 + t114;
	t113 = t89 * t92;
	t74 = -t96 * t113 + (-t97 * t109 + t93 * t94) * t90;
	t64 = atan2(-t66, t74);
	t61 = sin(t64);
	t62 = cos(t64);
	t60 = -t61 * t66 + t62 * t74;
	t59 = 0.1e1 / t60 ^ 2;
	t112 = t90 * t95;
	t102 = t89 * t112;
	t84 = -t92 * t105 - t104;
	t85 = -t92 * t106 + t103;
	t69 = -t96 * t102 - t84 * t109 + t85 * t93;
	t119 = t59 * t69;
	t118 = t62 * t66;
	t73 = 0.1e1 / t74 ^ 2;
	t117 = t66 * t73;
	t116 = t69 ^ 2 * t59;
	t70 = t85 * t96 + (t84 * t91 + t102) * t93;
	t78 = t91 * t112 - t84 * t89;
	t77 = 0.1e1 / t78 ^ 2;
	t115 = t70 * t77;
	t110 = t91 * t93;
	t108 = t93 * t97;
	t107 = t94 * t96;
	t100 = -t61 * t74 - t118;
	t99 = t93 * t101 + t82 * t110 - t83 * t96;
	t81 = (t91 * t107 + t108) * t90;
	t76 = 0.1e1 / t78;
	t75 = t93 * t113 + (t91 * t108 + t107) * t90;
	t72 = 0.1e1 / t74;
	t71 = t83 * t109 - t82 * t93;
	t65 = 0.1e1 / (t70 ^ 2 * t77 + 0.1e1);
	t63 = 0.1e1 / (t66 ^ 2 * t73 + 0.1e1);
	t58 = 0.1e1 / t60;
	t57 = 0.1e1 / (0.1e1 + t116);
	t56 = (t81 * t117 - t71 * t72) * t63;
	t55 = (t75 * t117 + t72 * t99) * t63;
	t1 = [-t69 * t72 * t63, t56, t55, 0, 0, 0; ((-t114 + (-t82 * t91 - t101) * t96) * t58 - (-t61 + (t72 * t118 + t61) * t63) * t116) * t57, ((t85 * t109 + t84 * t93) * t58 - (t100 * t56 - t61 * t71 + t62 * t81) * t119) * t57, (t70 * t58 - (t100 * t55 + t61 * t99 + t62 * t75) * t119) * t57, 0, 0, 0; (t99 * t76 - (t91 * t111 - t82 * t89) * t115) * t65, ((-t85 * t110 + t84 * t96) * t76 - t85 * t89 * t115) * t65, -t69 * t76 * t65, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (833->47), mult. (2499->112), div. (85->9), fcn. (3431->15), ass. (0->64)
	t105 = cos(pkin(7));
	t107 = sin(qJ(3));
	t111 = cos(qJ(3));
	t103 = sin(pkin(7));
	t104 = sin(pkin(6));
	t113 = cos(qJ(1));
	t126 = t104 * t113;
	t121 = t103 * t126;
	t108 = sin(qJ(2));
	t109 = sin(qJ(1));
	t112 = cos(qJ(2));
	t128 = cos(pkin(6));
	t119 = t113 * t128;
	t96 = t109 * t108 - t112 * t119;
	t97 = t108 * t119 + t109 * t112;
	t81 = (t105 * t96 + t121) * t111 + t97 * t107;
	t125 = t105 * t107;
	t114 = t107 * t121 - t97 * t111 + t96 * t125;
	t118 = t128 * t103;
	t92 = t107 * t118 + (t108 * t111 + t112 * t125) * t104;
	t80 = atan2(t114, t92);
	t77 = sin(t80);
	t78 = cos(t80);
	t71 = t114 * t77 + t78 * t92;
	t70 = 0.1e1 / t71 ^ 2;
	t127 = t104 * t109;
	t120 = t109 * t128;
	t98 = -t113 * t108 - t112 * t120;
	t116 = t103 * t127 + t105 * t98;
	t99 = -t108 * t120 + t113 * t112;
	t129 = t99 * t111;
	t86 = t116 * t107 + t129;
	t136 = t70 * t86;
	t106 = sin(qJ(5));
	t110 = cos(qJ(5));
	t85 = t99 * t107 - t116 * t111;
	t94 = -t98 * t103 + t105 * t127;
	t76 = t85 * t106 + t94 * t110;
	t74 = 0.1e1 / t76 ^ 2;
	t75 = t94 * t106 - t85 * t110;
	t135 = t74 * t75;
	t134 = t78 * t114;
	t90 = 0.1e1 / t92 ^ 2;
	t133 = t114 * t90;
	t132 = t86 ^ 2 * t70;
	t131 = t103 * t99;
	t124 = t107 * t108;
	t123 = t111 * t112;
	t122 = t75 ^ 2 * t74 + 0.1e1;
	t117 = -t77 * t92 + t134;
	t95 = (-t105 * t124 + t123) * t104;
	t93 = -t96 * t103 + t105 * t126;
	t91 = t111 * t118 + (t105 * t123 - t124) * t104;
	t89 = 0.1e1 / t92;
	t88 = t105 * t129 + t98 * t107;
	t87 = -t96 * t111 - t97 * t125;
	t79 = 0.1e1 / (t114 ^ 2 * t90 + 0.1e1);
	t73 = 0.1e1 / t76;
	t72 = 0.1e1 / t122;
	t69 = 0.1e1 / t71;
	t68 = 0.1e1 / (0.1e1 + t132);
	t67 = (-t95 * t133 - t87 * t89) * t79;
	t66 = (-t91 * t133 + t81 * t89) * t79;
	t1 = [-t86 * t89 * t79, t67, t66, 0, 0, 0; (t114 * t69 - (-t77 + (-t89 * t134 + t77) * t79) * t132) * t68, ((t98 * t111 - t99 * t125) * t69 - (t117 * t67 - t77 * t87 + t78 * t95) * t136) * t68, (-t85 * t69 - (t117 * t66 + t77 * t81 + t78 * t91) * t136) * t68, 0, 0, 0; ((t93 * t106 + t110 * t81) * t73 - (-t106 * t81 + t93 * t110) * t135) * t72, ((t106 * t131 - t88 * t110) * t73 - (t88 * t106 + t110 * t131) * t135) * t72, (-t106 * t135 - t110 * t73) * t86 * t72, 0, t122 * t72, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:24
	% DurationCPUTime: 0.85s
	% Computational Cost: add. (1916->67), mult. (5494->156), div. (115->9), fcn. (7543->17), ass. (0->78)
	t142 = cos(pkin(6));
	t151 = cos(qJ(2));
	t152 = cos(qJ(1));
	t158 = t152 * t151;
	t146 = sin(qJ(2));
	t147 = sin(qJ(1));
	t161 = t147 * t146;
	t133 = -t142 * t158 + t161;
	t139 = sin(pkin(7));
	t141 = cos(pkin(7));
	t140 = sin(pkin(6));
	t165 = t140 * t152;
	t128 = -t133 * t139 + t141 * t165;
	t144 = sin(qJ(5));
	t149 = cos(qJ(5));
	t159 = t152 * t146;
	t160 = t147 * t151;
	t134 = t142 * t159 + t160;
	t145 = sin(qJ(3));
	t150 = cos(qJ(3));
	t154 = t133 * t141 + t139 * t165;
	t153 = t134 * t145 + t154 * t150;
	t108 = -t128 * t149 + t153 * t144;
	t180 = t128 * t144 + t153 * t149;
	t121 = -t134 * t150 + t154 * t145;
	t135 = -t142 * t160 - t159;
	t166 = t140 * t147;
	t130 = -t135 * t139 + t141 * t166;
	t156 = t139 * t166;
	t164 = t141 * t150;
	t136 = -t142 * t161 + t158;
	t169 = t136 * t145;
	t155 = -t135 * t164 - t150 * t156 + t169;
	t110 = t130 * t149 + t155 * t144;
	t148 = cos(qJ(6));
	t122 = t136 * t150 + (t135 * t141 + t156) * t145;
	t143 = sin(qJ(6));
	t173 = t122 * t143;
	t100 = t110 * t148 + t173;
	t98 = 0.1e1 / t100 ^ 2;
	t172 = t122 * t148;
	t99 = t110 * t143 - t172;
	t177 = t98 * t99;
	t109 = t130 * t144 - t155 * t149;
	t168 = t139 * t142;
	t126 = -t150 * t168 + (t145 * t146 - t151 * t164) * t140;
	t132 = -t140 * t151 * t139 + t142 * t141;
	t117 = -t126 * t149 + t132 * t144;
	t104 = atan2(t180, t117);
	t101 = sin(t104);
	t102 = cos(t104);
	t95 = t101 * t180 + t102 * t117;
	t94 = 0.1e1 / t95 ^ 2;
	t176 = t109 * t94;
	t175 = t109 ^ 2 * t94;
	t116 = 0.1e1 / t117 ^ 2;
	t174 = t180 * t116;
	t167 = t139 * t144;
	t163 = t145 * t151;
	t162 = t146 * t150;
	t157 = t99 ^ 2 * t98 + 0.1e1;
	t127 = t145 * t168 + (t141 * t163 + t162) * t140;
	t125 = (t146 * t167 - (t141 * t162 + t163) * t149) * t140;
	t124 = t135 * t150 - t141 * t169;
	t123 = t135 * t145 + t136 * t164;
	t118 = t126 * t144 + t132 * t149;
	t115 = 0.1e1 / t117;
	t112 = t136 * t139 * t149 + t123 * t144;
	t111 = -t134 * t167 + (-t133 * t145 + t134 * t164) * t149;
	t103 = 0.1e1 / (t116 * t180 ^ 2 + 0.1e1);
	t97 = 0.1e1 / t100;
	t96 = 0.1e1 / t157;
	t93 = 0.1e1 / t95;
	t92 = 0.1e1 / (0.1e1 + t175);
	t91 = (-t115 * t121 + t127 * t174) * t149 * t103;
	t90 = (t111 * t115 - t125 * t174) * t103;
	t89 = (-t108 * t115 - t118 * t174) * t103;
	t1 = [-t109 * t115 * t103, t90, t91, 0, t89, 0; (t180 * t93 - (-t101 + (-t102 * t115 * t180 + t101) * t103) * t175) * t92, ((-t123 * t149 + t136 * t167) * t93 - ((t180 * t90 + t125) * t102 + (-t117 * t90 + t111) * t101) * t176) * t92, (-t122 * t149 * t93 - ((-t127 * t149 + t180 * t91) * t102 + (-t117 * t91 - t121 * t149) * t101) * t176) * t92, 0, (t110 * t93 - ((t180 * t89 + t118) * t102 + (-t117 * t89 - t108) * t101) * t176) * t92, 0; ((-t108 * t143 - t121 * t148) * t97 - (-t108 * t148 + t121 * t143) * t177) * t96, ((t112 * t143 - t124 * t148) * t97 - (t112 * t148 + t124 * t143) * t177) * t96, ((t144 * t173 + t155 * t148) * t97 - (-t155 * t143 + t144 * t172) * t177) * t96, 0, (-t143 * t97 + t148 * t177) * t96 * t109, t157 * t96;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end