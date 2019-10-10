% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR14
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
%   Wie in S6RRRRPR14_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:54
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR14_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR14_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
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
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:13
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
	% StartTime: 2019-10-10 12:54:13
	% EndTime: 2019-10-10 12:54:13
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
	% StartTime: 2019-10-10 12:54:13
	% EndTime: 2019-10-10 12:54:14
	% DurationCPUTime: 0.79s
	% Computational Cost: add. (1829->65), mult. (5250->153), div. (110->9), fcn. (7214->17), ass. (0->79)
	t138 = cos(pkin(6));
	t144 = cos(qJ(2));
	t174 = sin(qJ(1));
	t153 = t174 * t144;
	t141 = sin(qJ(2));
	t145 = cos(qJ(1));
	t157 = t145 * t141;
	t129 = t138 * t157 + t153;
	t140 = sin(qJ(3));
	t143 = cos(qJ(3));
	t154 = t174 * t141;
	t156 = t145 * t144;
	t128 = -t138 * t156 + t154;
	t134 = sin(pkin(7));
	t137 = cos(pkin(7));
	t135 = sin(pkin(6));
	t161 = t135 * t145;
	t150 = t128 * t137 + t134 * t161;
	t115 = -t129 * t143 + t150 * t140;
	t125 = -t128 * t134 + t137 * t161;
	t139 = sin(qJ(4));
	t142 = cos(qJ(4));
	t175 = t115 * t139 - t125 * t142;
	t103 = t115 * t142 + t125 * t139;
	t149 = t138 * t153 + t157;
	t155 = t135 * t174;
	t152 = t134 * t155;
	t130 = -t138 * t154 + t156;
	t164 = t130 * t143;
	t117 = t164 + (-t149 * t137 + t152) * t140;
	t146 = t149 * t134 + t137 * t155;
	t105 = t117 * t142 + t146 * t139;
	t148 = t149 * t143;
	t116 = t130 * t140 + t137 * t148 - t143 * t152;
	t133 = sin(pkin(13));
	t136 = cos(pkin(13));
	t95 = t105 * t136 + t116 * t133;
	t93 = 0.1e1 / t95 ^ 2;
	t94 = t105 * t133 - t116 * t136;
	t173 = t93 * t94;
	t160 = t137 * t140;
	t163 = t134 * t138;
	t124 = t140 * t163 + (t141 * t143 + t144 * t160) * t135;
	t127 = -t135 * t144 * t134 + t138 * t137;
	t110 = t124 * t139 - t127 * t142;
	t99 = atan2(t175, t110);
	t97 = cos(t99);
	t172 = t175 * t97;
	t104 = t117 * t139 - t146 * t142;
	t96 = sin(t99);
	t90 = t97 * t110 + t175 * t96;
	t89 = 0.1e1 / t90 ^ 2;
	t171 = t104 * t89;
	t170 = t104 ^ 2 * t89;
	t109 = 0.1e1 / t110 ^ 2;
	t169 = t175 * t109;
	t168 = t116 * t142;
	t162 = t134 * t142;
	t159 = t140 * t141;
	t158 = t143 * t144;
	t151 = -t110 * t96 + t172;
	t147 = -t129 * t140 - t150 * t143;
	t123 = t143 * t163 + (t137 * t158 - t159) * t135;
	t120 = ((-t137 * t159 + t158) * t139 - t141 * t162) * t135;
	t119 = -t130 * t160 - t148;
	t118 = t137 * t164 - t149 * t140;
	t111 = t124 * t142 + t127 * t139;
	t108 = 0.1e1 / t110;
	t107 = t130 * t134 * t139 + t119 * t142;
	t106 = (-t128 * t143 - t129 * t160) * t139 - t129 * t162;
	t98 = 0.1e1 / (t109 * t175 ^ 2 + 0.1e1);
	t92 = 0.1e1 / t95;
	t91 = 0.1e1 / (t94 ^ 2 * t93 + 0.1e1);
	t88 = 0.1e1 / t90;
	t87 = 0.1e1 / (0.1e1 + t170);
	t86 = (-t108 * t147 - t123 * t169) * t98 * t139;
	t85 = (-t106 * t108 - t120 * t169) * t98;
	t84 = (t103 * t108 - t111 * t169) * t98;
	t1 = [-t104 * t108 * t98, t85, t86, t84, 0, 0; (t175 * t88 - (-t96 + (-t108 * t172 + t96) * t98) * t170) * t87, ((t119 * t139 - t130 * t162) * t88 - (-t96 * t106 + t97 * t120 + t151 * t85) * t171) * t87, (-t116 * t139 * t88 - (t151 * t86 + (t123 * t97 - t147 * t96) * t139) * t171) * t87, (t105 * t88 - (t103 * t96 + t97 * t111 + t151 * t84) * t171) * t87, 0, 0; ((t103 * t133 - t136 * t147) * t92 - (t103 * t136 + t133 * t147) * t173) * t91, ((t107 * t133 - t118 * t136) * t92 - (t107 * t136 + t118 * t133) * t173) * t91, ((-t117 * t136 - t133 * t168) * t92 - (t117 * t133 - t136 * t168) * t173) * t91, (-t133 * t92 + t136 * t173) * t91 * t104, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:13
	% EndTime: 2019-10-10 12:54:14
	% DurationCPUTime: 0.84s
	% Computational Cost: add. (2000->68), mult. (5494->157), div. (115->9), fcn. (7543->17), ass. (0->79)
	t145 = cos(pkin(6));
	t151 = cos(qJ(2));
	t180 = sin(qJ(1));
	t159 = t180 * t151;
	t148 = sin(qJ(2));
	t152 = cos(qJ(1));
	t164 = t152 * t148;
	t135 = t145 * t164 + t159;
	t147 = sin(qJ(3));
	t150 = cos(qJ(3));
	t160 = t180 * t148;
	t163 = t152 * t151;
	t134 = -t145 * t163 + t160;
	t142 = sin(pkin(7));
	t144 = cos(pkin(7));
	t143 = sin(pkin(6));
	t168 = t143 * t152;
	t157 = t134 * t144 + t142 * t168;
	t121 = -t135 * t150 + t157 * t147;
	t131 = -t134 * t142 + t144 * t168;
	t146 = sin(qJ(4));
	t149 = cos(qJ(4));
	t181 = t121 * t146 - t131 * t149;
	t109 = t121 * t149 + t131 * t146;
	t156 = t145 * t159 + t164;
	t161 = t143 * t180;
	t158 = t142 * t161;
	t136 = -t145 * t160 + t163;
	t171 = t136 * t150;
	t123 = t171 + (-t156 * t144 + t158) * t147;
	t153 = t156 * t142 + t144 * t161;
	t111 = t123 * t149 + t153 * t146;
	t155 = t156 * t150;
	t122 = t136 * t147 + t144 * t155 - t150 * t158;
	t141 = pkin(13) + qJ(6);
	t139 = sin(t141);
	t140 = cos(t141);
	t100 = t111 * t139 - t122 * t140;
	t101 = t111 * t140 + t122 * t139;
	t99 = 0.1e1 / t101 ^ 2;
	t179 = t100 * t99;
	t110 = t123 * t146 - t153 * t149;
	t167 = t144 * t147;
	t170 = t142 * t145;
	t130 = t147 * t170 + (t148 * t150 + t151 * t167) * t143;
	t133 = -t143 * t151 * t142 + t145 * t144;
	t116 = t130 * t146 - t133 * t149;
	t105 = atan2(t181, t116);
	t102 = sin(t105);
	t103 = cos(t105);
	t96 = t102 * t181 + t103 * t116;
	t95 = 0.1e1 / t96 ^ 2;
	t178 = t110 * t95;
	t177 = t110 ^ 2 * t95;
	t115 = 0.1e1 / t116 ^ 2;
	t176 = t181 * t115;
	t175 = t122 * t149;
	t169 = t142 * t149;
	t166 = t147 * t148;
	t165 = t150 * t151;
	t162 = t100 ^ 2 * t99 + 0.1e1;
	t154 = -t135 * t147 - t157 * t150;
	t129 = t150 * t170 + (t144 * t165 - t166) * t143;
	t126 = ((-t144 * t166 + t165) * t146 - t148 * t169) * t143;
	t125 = -t136 * t167 - t155;
	t124 = t144 * t171 - t156 * t147;
	t117 = t130 * t149 + t133 * t146;
	t114 = 0.1e1 / t116;
	t113 = t136 * t142 * t146 + t125 * t149;
	t112 = (-t134 * t150 - t135 * t167) * t146 - t135 * t169;
	t104 = 0.1e1 / (t115 * t181 ^ 2 + 0.1e1);
	t98 = 0.1e1 / t101;
	t97 = 0.1e1 / t162;
	t94 = 0.1e1 / t96;
	t93 = 0.1e1 / (0.1e1 + t177);
	t92 = (-t114 * t154 - t129 * t176) * t146 * t104;
	t91 = (-t112 * t114 - t126 * t176) * t104;
	t90 = (t109 * t114 - t117 * t176) * t104;
	t1 = [-t110 * t114 * t104, t91, t92, t90, 0, 0; (t181 * t94 - (-t102 + (-t103 * t114 * t181 + t102) * t104) * t177) * t93, ((t125 * t146 - t136 * t169) * t94 - ((t181 * t91 + t126) * t103 + (-t116 * t91 - t112) * t102) * t178) * t93, (-t122 * t146 * t94 - ((t129 * t146 + t181 * t92) * t103 + (-t116 * t92 - t146 * t154) * t102) * t178) * t93, (t111 * t94 - ((t181 * t90 + t117) * t103 + (-t116 * t90 + t109) * t102) * t178) * t93, 0, 0; ((t109 * t139 - t140 * t154) * t98 - (t109 * t140 + t139 * t154) * t179) * t97, ((t113 * t139 - t124 * t140) * t98 - (t113 * t140 + t124 * t139) * t179) * t97, ((-t123 * t140 - t139 * t175) * t98 - (t123 * t139 - t140 * t175) * t179) * t97, (-t139 * t98 + t140 * t179) * t97 * t110, 0, t162 * t97;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end