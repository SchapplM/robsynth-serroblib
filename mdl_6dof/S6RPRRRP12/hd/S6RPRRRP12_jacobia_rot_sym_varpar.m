% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP12
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
%   Wie in S6RPRRRP12_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:58
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRP12_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP12_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:24
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (25->13), mult. (89->32), div. (20->9), fcn. (140->9), ass. (0->23)
	t28 = cos(pkin(6));
	t26 = sin(pkin(6));
	t30 = cos(qJ(1));
	t32 = t30 * t26;
	t22 = atan2(t32, t28);
	t19 = sin(t22);
	t20 = cos(t22);
	t14 = t19 * t32 + t20 * t28;
	t29 = sin(qJ(1));
	t37 = 0.1e1 / t14 ^ 2 * t29 ^ 2;
	t23 = t26 ^ 2;
	t21 = 0.1e1 / (0.1e1 + t30 ^ 2 * t23 / t28 ^ 2);
	t36 = t21 / t28;
	t25 = sin(pkin(12));
	t35 = t29 * t25;
	t27 = cos(pkin(12));
	t34 = t29 * t27;
	t33 = t30 * t25;
	t31 = t30 * t27;
	t18 = -t28 * t35 + t31;
	t17 = t28 * t34 + t33;
	t16 = 0.1e1 / t18 ^ 2;
	t1 = [-t29 * t26 * t36, 0, 0, 0, 0, 0; (0.1e1 / t14 * t32 - (-t20 * t23 * t30 * t36 + (t21 - 0.1e1) * t26 * t19) * t26 * t37) / (t23 * t37 + 0.1e1), 0, 0, 0, 0, 0; ((t28 * t31 - t35) / t18 - (-t28 * t33 - t34) * t17 * t16) / (t17 ^ 2 * t16 + 0.1e1), 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:24
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (156->21), mult. (450->48), div. (25->9), fcn. (637->13), ass. (0->38)
	t61 = cos(pkin(6));
	t59 = cos(pkin(12));
	t65 = cos(qJ(1));
	t69 = t65 * t59;
	t56 = sin(pkin(12));
	t63 = sin(qJ(1));
	t72 = t63 * t56;
	t51 = -t61 * t69 + t72;
	t57 = sin(pkin(7));
	t60 = cos(pkin(7));
	t58 = sin(pkin(6));
	t73 = t58 * t65;
	t46 = -t51 * t57 + t60 * t73;
	t50 = -t58 * t59 * t57 + t61 * t60;
	t45 = atan2(t46, t50);
	t42 = sin(t45);
	t43 = cos(t45);
	t37 = t42 * t46 + t43 * t50;
	t70 = t65 * t56;
	t71 = t63 * t59;
	t53 = -t61 * t71 - t70;
	t74 = t58 * t63;
	t47 = t53 * t57 - t60 * t74;
	t75 = t47 ^ 2 / t37 ^ 2;
	t54 = -t61 * t72 + t69;
	t62 = sin(qJ(3));
	t64 = cos(qJ(3));
	t66 = t53 * t60 + t57 * t74;
	t41 = t54 * t64 + t66 * t62;
	t39 = 0.1e1 / t41 ^ 2;
	t40 = t54 * t62 - t66 * t64;
	t68 = t40 ^ 2 * t39 + 0.1e1;
	t67 = t51 * t60 + t57 * t73;
	t52 = -t61 * t70 - t71;
	t49 = 0.1e1 / t50;
	t44 = 0.1e1 / (0.1e1 + t46 ^ 2 / t50 ^ 2);
	t38 = 0.1e1 / t68;
	t1 = [t47 * t49 * t44, 0, 0, 0, 0, 0; (t46 / t37 + (t42 + (t43 * t46 * t49 - t42) * t44) * t75) / (0.1e1 + t75), 0, 0, 0, 0, 0; ((t52 * t62 - t67 * t64) / t41 - (t52 * t64 + t67 * t62) * t40 * t39) * t38, 0, t68 * t38, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:24
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (559->39), mult. (1669->87), div. (55->9), fcn. (2293->15), ass. (0->57)
	t115 = sin(qJ(1));
	t88 = sin(pkin(6));
	t102 = t88 * t115;
	t87 = sin(pkin(7));
	t90 = cos(pkin(7));
	t89 = cos(pkin(12));
	t100 = t115 * t89;
	t86 = sin(pkin(12));
	t96 = cos(qJ(1));
	t106 = t96 * t86;
	t91 = cos(pkin(6));
	t98 = t91 * t100 + t106;
	t116 = -t87 * t102 + t98 * t90;
	t101 = t115 * t86;
	t105 = t96 * t89;
	t83 = -t91 * t101 + t105;
	t93 = sin(qJ(3));
	t95 = cos(qJ(3));
	t72 = -t116 * t93 + t83 * t95;
	t78 = t90 * t102 + t98 * t87;
	t92 = sin(qJ(4));
	t94 = cos(qJ(4));
	t62 = t72 * t94 + t78 * t92;
	t60 = 0.1e1 / t62 ^ 2;
	t61 = t72 * t92 - t78 * t94;
	t114 = t60 * t61;
	t109 = t88 * t96;
	t104 = t87 * t109;
	t107 = t90 * t95;
	t82 = t91 * t106 + t100;
	t111 = t82 * t93;
	t81 = -t91 * t105 + t101;
	t67 = t95 * t104 + t81 * t107 + t111;
	t110 = t87 * t91;
	t75 = -t95 * t110 + (-t89 * t107 + t86 * t93) * t88;
	t66 = atan2(-t67, t75);
	t64 = cos(t66);
	t113 = t64 * t67;
	t63 = sin(t66);
	t57 = -t63 * t67 + t64 * t75;
	t56 = 0.1e1 / t57 ^ 2;
	t71 = t116 * t95 + t83 * t93;
	t112 = t71 ^ 2 * t56;
	t108 = t90 * t93;
	t103 = t61 ^ 2 * t60 + 0.1e1;
	t70 = t93 * t104 + t81 * t108 - t82 * t95;
	t77 = t90 * t109 - t81 * t87;
	t76 = t93 * t110 + (t89 * t108 + t86 * t95) * t88;
	t74 = 0.1e1 / t75 ^ 2;
	t73 = 0.1e1 / t75;
	t65 = 0.1e1 / (t67 ^ 2 * t74 + 0.1e1);
	t59 = 0.1e1 / t62;
	t58 = 0.1e1 / t103;
	t55 = 0.1e1 / t57;
	t54 = 0.1e1 / (0.1e1 + t112);
	t53 = (t67 * t74 * t76 + t70 * t73) * t65;
	t1 = [-t71 * t73 * t65, 0, t53, 0, 0, 0; ((-t111 + (-t81 * t90 - t104) * t95) * t55 - (-t63 + (t73 * t113 + t63) * t65) * t112) * t54, 0, (t72 * t55 - (t63 * t70 + t64 * t76 + (-t63 * t75 - t113) * t53) * t71 * t56) * t54, 0, 0, 0; ((t70 * t92 - t77 * t94) * t59 - (t70 * t94 + t77 * t92) * t114) * t58, 0, (t94 * t114 - t92 * t59) * t71 * t58, t103 * t58, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:25
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (1440->49), mult. (4121->118), div. (85->9), fcn. (5660->17), ass. (0->69)
	t119 = cos(pkin(6));
	t117 = cos(pkin(12));
	t152 = sin(qJ(1));
	t134 = t152 * t117;
	t114 = sin(pkin(12));
	t126 = cos(qJ(1));
	t139 = t126 * t114;
	t111 = t119 * t139 + t134;
	t122 = sin(qJ(3));
	t125 = cos(qJ(3));
	t135 = t152 * t114;
	t138 = t126 * t117;
	t110 = -t119 * t138 + t135;
	t115 = sin(pkin(7));
	t118 = cos(pkin(7));
	t116 = sin(pkin(6));
	t141 = t116 * t126;
	t131 = t110 * t118 + t115 * t141;
	t100 = -t111 * t125 + t122 * t131;
	t107 = -t110 * t115 + t118 * t141;
	t121 = sin(qJ(4));
	t124 = cos(qJ(4));
	t154 = t100 * t121 - t107 * t124;
	t90 = t100 * t124 + t107 * t121;
	t130 = t119 * t134 + t139;
	t136 = t116 * t152;
	t153 = -t115 * t136 + t130 * t118;
	t140 = t117 * t118;
	t142 = t115 * t119;
	t106 = t122 * t142 + (t114 * t125 + t122 * t140) * t116;
	t109 = -t116 * t117 * t115 + t119 * t118;
	t95 = t106 * t121 - t109 * t124;
	t86 = atan2(t154, t95);
	t83 = sin(t86);
	t84 = cos(t86);
	t77 = t154 * t83 + t84 * t95;
	t76 = 0.1e1 / t77 ^ 2;
	t112 = -t119 * t135 + t138;
	t102 = t112 * t125 - t153 * t122;
	t127 = t115 * t130 + t118 * t136;
	t91 = t102 * t121 - t124 * t127;
	t151 = t76 * t91;
	t101 = t112 * t122 + t153 * t125;
	t120 = sin(qJ(5));
	t123 = cos(qJ(5));
	t92 = t102 * t124 + t121 * t127;
	t82 = t101 * t120 + t92 * t123;
	t80 = 0.1e1 / t82 ^ 2;
	t81 = -t101 * t123 + t92 * t120;
	t150 = t80 * t81;
	t149 = t84 * t154;
	t94 = 0.1e1 / t95 ^ 2;
	t148 = t154 * t94;
	t147 = t91 ^ 2 * t76;
	t146 = t101 * t124;
	t137 = t81 ^ 2 * t80 + 0.1e1;
	t132 = -t83 * t95 + t149;
	t128 = -t111 * t122 - t125 * t131;
	t105 = t125 * t142 + (-t114 * t122 + t125 * t140) * t116;
	t96 = t106 * t124 + t109 * t121;
	t93 = 0.1e1 / t95;
	t85 = 0.1e1 / (t154 ^ 2 * t94 + 0.1e1);
	t79 = 0.1e1 / t82;
	t78 = 0.1e1 / t137;
	t75 = 0.1e1 / t77;
	t74 = 0.1e1 / (0.1e1 + t147);
	t73 = (-t105 * t148 - t128 * t93) * t85 * t121;
	t72 = (-t148 * t96 + t90 * t93) * t85;
	t1 = [-t91 * t93 * t85, 0, t73, t72, 0, 0; (t154 * t75 - (-t83 + (-t149 * t93 + t83) * t85) * t147) * t74, 0, (-t101 * t121 * t75 - (t132 * t73 + (t105 * t84 - t128 * t83) * t121) * t151) * t74, (t92 * t75 - (t132 * t72 + t83 * t90 + t84 * t96) * t151) * t74, 0, 0; ((t90 * t120 - t123 * t128) * t79 - (t120 * t128 + t90 * t123) * t150) * t78, 0, ((-t102 * t123 - t120 * t146) * t79 - (t102 * t120 - t123 * t146) * t150) * t78, (-t120 * t79 + t123 * t150) * t91 * t78, t137 * t78, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:25
	% DurationCPUTime: 0.95s
	% Computational Cost: add. (2697->57), mult. (7742->135), div. (107->9), fcn. (10572->17), ass. (0->77)
	t145 = cos(pkin(6));
	t140 = sin(pkin(12));
	t153 = cos(qJ(1));
	t161 = t153 * t140;
	t143 = cos(pkin(12));
	t149 = sin(qJ(1));
	t162 = t149 * t143;
	t137 = t145 * t161 + t162;
	t148 = sin(qJ(3));
	t152 = cos(qJ(3));
	t160 = t153 * t143;
	t163 = t149 * t140;
	t136 = -t145 * t160 + t163;
	t141 = sin(pkin(7));
	t144 = cos(pkin(7));
	t142 = sin(pkin(6));
	t166 = t142 * t153;
	t157 = t136 * t144 + t141 * t166;
	t126 = -t137 * t152 + t157 * t148;
	t131 = -t136 * t141 + t144 * t166;
	t147 = sin(qJ(4));
	t151 = cos(qJ(4));
	t116 = t126 * t151 + t131 * t147;
	t146 = sin(qJ(5));
	t150 = cos(qJ(5));
	t154 = t137 * t148 + t157 * t152;
	t184 = t116 * t146 + t154 * t150;
	t183 = t116 * t150 - t154 * t146;
	t114 = t126 * t147 - t131 * t151;
	t156 = t145 * t162 + t161;
	t167 = t142 * t149;
	t178 = -t141 * t167 + t156 * t144;
	t165 = t143 * t144;
	t168 = t141 * t145;
	t130 = t148 * t168 + (t140 * t152 + t148 * t165) * t142;
	t135 = -t142 * t143 * t141 + t145 * t144;
	t123 = t130 * t151 + t135 * t147;
	t129 = -t152 * t168 + (t140 * t148 - t152 * t165) * t142;
	t111 = t123 * t146 - t129 * t150;
	t99 = atan2(t184, t111);
	t97 = cos(t99);
	t177 = t184 * t97;
	t138 = -t145 * t163 + t160;
	t128 = t138 * t152 - t178 * t148;
	t133 = t156 * t141 + t144 * t167;
	t118 = t128 * t151 + t133 * t147;
	t127 = t138 * t148 + t178 * t152;
	t171 = t127 * t150;
	t106 = t118 * t146 - t171;
	t96 = sin(t99);
	t95 = t97 * t111 + t184 * t96;
	t94 = 0.1e1 / t95 ^ 2;
	t176 = t106 * t94;
	t175 = t106 ^ 2 * t94;
	t107 = t118 * t150 + t127 * t146;
	t102 = 0.1e1 / t107 ^ 2;
	t117 = -t128 * t147 + t133 * t151;
	t174 = t102 * t117;
	t110 = 0.1e1 / t111 ^ 2;
	t173 = t184 * t110;
	t172 = t117 ^ 2 * t102;
	t164 = t146 * t151;
	t158 = -t111 * t96 + t177;
	t122 = -t130 * t147 + t135 * t151;
	t119 = -t129 * t164 - t130 * t150;
	t112 = t123 * t150 + t129 * t146;
	t109 = 0.1e1 / t111;
	t108 = t126 * t150 - t154 * t164;
	t101 = 0.1e1 / t107;
	t100 = 0.1e1 / (0.1e1 + t172);
	t98 = 0.1e1 / (t110 * t184 ^ 2 + 0.1e1);
	t93 = 0.1e1 / t95;
	t92 = 0.1e1 / (0.1e1 + t175);
	t91 = (-t109 * t114 - t122 * t173) * t98 * t146;
	t90 = (-t108 * t109 - t119 * t173) * t98;
	t89 = (t109 * t183 - t112 * t173) * t98;
	t1 = [-t106 * t109 * t98, 0, t90, t91, t89, 0; (t184 * t93 - (-t96 + (-t109 * t177 + t96) * t98) * t175) * t92, 0, ((-t127 * t164 - t128 * t150) * t93 - (-t96 * t108 + t97 * t119 + t158 * t90) * t176) * t92, (t117 * t146 * t93 - (t158 * t91 + (-t114 * t96 + t122 * t97) * t146) * t176) * t92, (t107 * t93 - (t97 * t112 + t158 * t89 + t183 * t96) * t176) * t92, 0; (-t114 * t101 - t183 * t174) * t100, 0, (t127 * t147 * t101 - (t128 * t146 - t151 * t171) * t174) * t100, (-t101 * t118 - t150 * t172) * t100, t106 * t100 * t174, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end