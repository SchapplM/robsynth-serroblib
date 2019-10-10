% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR9
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
%   Wie in S6RRRPRR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:08
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR9_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR9_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
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
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
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
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (391->34), mult. (1119->83), div. (56->9), fcn. (1577->15), ass. (0->53)
	t79 = sin(pkin(7));
	t78 = sin(pkin(13));
	t81 = cos(pkin(13));
	t84 = sin(qJ(3));
	t87 = cos(qJ(3));
	t90 = t87 * t78 + t84 * t81;
	t63 = t90 * t79;
	t82 = cos(pkin(7));
	t65 = t90 * t82;
	t83 = cos(pkin(6));
	t85 = sin(qJ(2));
	t89 = cos(qJ(1));
	t92 = t89 * t85;
	t86 = sin(qJ(1));
	t88 = cos(qJ(2));
	t93 = t86 * t88;
	t71 = -t83 * t93 - t92;
	t91 = t89 * t88;
	t94 = t86 * t85;
	t72 = -t83 * t94 + t91;
	t73 = t84 * t78 - t87 * t81;
	t80 = sin(pkin(6));
	t96 = t80 * t86;
	t54 = t63 * t96 + t71 * t65 - t72 * t73;
	t51 = 0.1e1 / t54 ^ 2;
	t62 = t73 * t79;
	t64 = t73 * t82;
	t52 = -t62 * t96 - t71 * t64 - t72 * t90;
	t101 = t51 * t52;
	t100 = t52 ^ 2 * t51;
	t69 = -t83 * t91 + t94;
	t95 = t80 * t89;
	t59 = -t69 * t79 + t82 * t95;
	t68 = -t80 * t88 * t79 + t83 * t82;
	t58 = atan2(t59, t68);
	t56 = cos(t58);
	t99 = t56 * t59;
	t55 = sin(t58);
	t49 = t55 * t59 + t56 * t68;
	t48 = 0.1e1 / t49 ^ 2;
	t60 = t71 * t79 - t82 * t96;
	t98 = t60 ^ 2 * t48;
	t97 = t80 * t85;
	t70 = -t83 * t92 - t93;
	t67 = 0.1e1 / t68 ^ 2;
	t66 = 0.1e1 / t68;
	t57 = 0.1e1 / (t59 ^ 2 * t67 + 0.1e1);
	t50 = 0.1e1 / t54;
	t47 = 0.1e1 / t49;
	t46 = 0.1e1 / (0.1e1 + t98);
	t45 = 0.1e1 / (0.1e1 + t100);
	t44 = (-t59 * t67 * t97 + t66 * t70) * t79 * t57;
	t1 = [t60 * t66 * t57, t44, 0, 0, 0, 0; (t59 * t47 + (t55 + (t66 * t99 - t55) * t57) * t98) * t46, (t72 * t79 * t47 + ((t55 * t70 + t56 * t97) * t79 + (-t55 * t68 + t99) * t44) * t60 * t48) * t46, 0, 0, 0, 0; ((t62 * t95 + t69 * t64 + t70 * t90) * t50 + (t63 * t95 + t69 * t65 - t70 * t73) * t101) * t45, ((-t72 * t64 + t71 * t90) * t50 + (-t72 * t65 - t71 * t73) * t101) * t45, (t54 * t50 + t100) * t45, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:08
	% DurationCPUTime: 0.48s
	% Computational Cost: add. (1467->49), mult. (4186->118), div. (85->9), fcn. (5752->17), ass. (0->66)
	t114 = sin(pkin(13));
	t117 = cos(pkin(13));
	t120 = sin(qJ(3));
	t124 = cos(qJ(3));
	t109 = t120 * t114 - t124 * t117;
	t118 = cos(pkin(7));
	t103 = t109 * t118;
	t121 = sin(qJ(2));
	t122 = sin(qJ(1));
	t125 = cos(qJ(2));
	t126 = cos(qJ(1));
	t138 = cos(pkin(6));
	t132 = t126 * t138;
	t105 = t122 * t121 - t125 * t132;
	t106 = t121 * t132 + t122 * t125;
	t116 = sin(pkin(6));
	t115 = sin(pkin(7));
	t129 = t109 * t115;
	t128 = t116 * t129;
	t130 = t124 * t114 + t120 * t117;
	t87 = t105 * t103 - t106 * t130 + t126 * t128;
	t95 = (t125 * t103 + t121 * t130) * t116 + t138 * t129;
	t81 = atan2(t87, t95);
	t78 = sin(t81);
	t79 = cos(t81);
	t76 = t78 * t87 + t79 * t95;
	t75 = 0.1e1 / t76 ^ 2;
	t133 = t122 * t138;
	t107 = -t126 * t121 - t125 * t133;
	t108 = -t121 * t133 + t126 * t125;
	t89 = -t107 * t103 - t108 * t130 - t122 * t128;
	t143 = t75 * t89;
	t142 = t79 * t87;
	t135 = t122 * t116;
	t100 = -t107 * t115 + t118 * t135;
	t119 = sin(qJ(5));
	t123 = cos(qJ(5));
	t102 = t130 * t115;
	t104 = t130 * t118;
	t91 = t102 * t135 + t107 * t104 - t108 * t109;
	t85 = t100 * t119 + t91 * t123;
	t83 = 0.1e1 / t85 ^ 2;
	t84 = -t100 * t123 + t91 * t119;
	t141 = t83 * t84;
	t93 = 0.1e1 / t95 ^ 2;
	t140 = t87 * t93;
	t139 = t89 ^ 2 * t75;
	t137 = t108 * t115;
	t136 = t116 * t126;
	t134 = t84 ^ 2 * t83 + 0.1e1;
	t131 = -t78 * t95 + t142;
	t127 = t102 * t136 + t105 * t104 + t106 * t109;
	t99 = -t105 * t115 + t118 * t136;
	t98 = (-t103 * t121 + t125 * t130) * t116;
	t97 = -t108 * t104 - t107 * t109;
	t96 = t106 * t103 + t105 * t130;
	t94 = t138 * t102 + (t104 * t125 - t109 * t121) * t116;
	t92 = 0.1e1 / t95;
	t82 = 0.1e1 / t85;
	t80 = 0.1e1 / (t87 ^ 2 * t93 + 0.1e1);
	t77 = 0.1e1 / t134;
	t74 = 0.1e1 / t76;
	t73 = 0.1e1 / (0.1e1 + t139);
	t72 = (-t98 * t140 + t92 * t96) * t80;
	t71 = (t127 * t92 - t94 * t140) * t80;
	t1 = [t89 * t92 * t80, t72, t71, 0, 0, 0; (t87 * t74 + (t78 + (t92 * t142 - t78) * t80) * t139) * t73, ((-t108 * t103 + t107 * t130) * t74 + (t131 * t72 + t78 * t96 + t79 * t98) * t143) * t73, (t91 * t74 + (t127 * t78 + t131 * t71 + t79 * t94) * t143) * t73, 0, 0, 0; ((t119 * t127 - t99 * t123) * t82 - (t99 * t119 + t123 * t127) * t141) * t77, ((t97 * t119 - t123 * t137) * t82 - (t119 * t137 + t97 * t123) * t141) * t77, (t119 * t82 - t123 * t141) * t89 * t77, 0, t134 * t77, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:08
	% EndTime: 2019-10-10 12:08:09
	% DurationCPUTime: 0.94s
	% Computational Cost: add. (2899->70), mult. (8111->162), div. (115->9), fcn. (11143->19), ass. (0->78)
	t153 = sin(pkin(13));
	t160 = sin(qJ(3));
	t186 = cos(pkin(13));
	t187 = cos(qJ(3));
	t150 = -t187 * t153 - t160 * t186;
	t154 = sin(pkin(7));
	t141 = t150 * t154;
	t156 = cos(pkin(7));
	t143 = t150 * t156;
	t157 = cos(pkin(6));
	t165 = cos(qJ(2));
	t166 = cos(qJ(1));
	t172 = t166 * t165;
	t161 = sin(qJ(2));
	t162 = sin(qJ(1));
	t175 = t162 * t161;
	t145 = -t157 * t172 + t175;
	t173 = t166 * t161;
	t174 = t162 * t165;
	t146 = t157 * t173 + t174;
	t169 = -t160 * t153 + t187 * t186;
	t155 = sin(pkin(6));
	t176 = t155 * t166;
	t128 = -t141 * t176 - t145 * t143 - t146 * t169;
	t139 = -t145 * t154 + t156 * t176;
	t159 = sin(qJ(5));
	t164 = cos(qJ(5));
	t188 = t128 * t159 - t139 * t164;
	t116 = t128 * t164 + t139 * t159;
	t133 = -t157 * t141 + (-t143 * t165 + t161 * t169) * t155;
	t144 = -t155 * t165 * t154 + t157 * t156;
	t121 = t133 * t159 - t144 * t164;
	t112 = atan2(t188, t121);
	t109 = sin(t112);
	t110 = cos(t112);
	t103 = t109 * t188 + t110 * t121;
	t102 = 0.1e1 / t103 ^ 2;
	t147 = -t157 * t174 - t173;
	t148 = -t157 * t175 + t172;
	t177 = t155 * t162;
	t167 = -t141 * t177 - t147 * t143 + t148 * t169;
	t170 = -t147 * t154 + t156 * t177;
	t117 = t159 * t167 - t170 * t164;
	t185 = t102 * t117;
	t118 = t170 * t159 + t164 * t167;
	t140 = t169 * t154;
	t142 = t169 * t156;
	t130 = t140 * t177 + t147 * t142 + t148 * t150;
	t158 = sin(qJ(6));
	t163 = cos(qJ(6));
	t108 = t118 * t163 - t130 * t158;
	t106 = 0.1e1 / t108 ^ 2;
	t107 = t118 * t158 + t130 * t163;
	t184 = t106 * t107;
	t120 = 0.1e1 / t121 ^ 2;
	t183 = t188 * t120;
	t182 = t117 ^ 2 * t102;
	t181 = t130 * t164;
	t178 = t154 * t164;
	t171 = t107 ^ 2 * t106 + 0.1e1;
	t168 = -t140 * t176 - t145 * t142 + t146 * t150;
	t136 = ((t143 * t161 + t165 * t169) * t159 - t161 * t178) * t155;
	t135 = t148 * t143 + t147 * t169;
	t134 = t148 * t142 - t147 * t150;
	t132 = t157 * t140 + (t142 * t165 + t150 * t161) * t155;
	t124 = t148 * t154 * t159 + t135 * t164;
	t123 = (t146 * t143 - t145 * t169) * t159 - t146 * t178;
	t122 = t133 * t164 + t144 * t159;
	t119 = 0.1e1 / t121;
	t111 = 0.1e1 / (t120 * t188 ^ 2 + 0.1e1);
	t105 = 0.1e1 / t108;
	t104 = 0.1e1 / t171;
	t101 = 0.1e1 / t103;
	t100 = 0.1e1 / (0.1e1 + t182);
	t99 = (-t119 * t123 - t136 * t183) * t111;
	t98 = (-t119 * t168 - t132 * t183) * t159 * t111;
	t97 = (t116 * t119 - t122 * t183) * t111;
	t1 = [-t117 * t119 * t111, t99, t98, 0, t97, 0; (t188 * t101 - (-t109 + (-t110 * t119 * t188 + t109) * t111) * t182) * t100, ((t135 * t159 - t148 * t178) * t101 - ((t188 * t99 + t136) * t110 + (-t121 * t99 - t123) * t109) * t185) * t100, (t130 * t159 * t101 - ((t132 * t159 + t188 * t98) * t110 + (-t121 * t98 - t159 * t168) * t109) * t185) * t100, 0, (t118 * t101 - ((t188 * t97 + t122) * t110 + (-t121 * t97 + t116) * t109) * t185) * t100, 0; ((t116 * t158 - t163 * t168) * t105 - (t116 * t163 + t158 * t168) * t184) * t104, ((t124 * t158 - t134 * t163) * t105 - (t124 * t163 + t134 * t158) * t184) * t104, ((t158 * t181 - t163 * t167) * t105 - (t158 * t167 + t163 * t181) * t184) * t104, 0, (-t158 * t105 + t163 * t184) * t117 * t104, t171 * t104;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end