% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR11
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
%   Wie in S6RPRRPR11_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR11_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR11_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
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
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.12s
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
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
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
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:21
	% DurationCPUTime: 0.69s
	% Computational Cost: add. (1353->49), mult. (3877->116), div. (80->9), fcn. (5331->17), ass. (0->69)
	t114 = sin(qJ(4));
	t116 = cos(qJ(4));
	t113 = cos(pkin(6));
	t111 = cos(pkin(12));
	t144 = sin(qJ(1));
	t126 = t144 * t111;
	t107 = sin(pkin(12));
	t118 = cos(qJ(1));
	t130 = t118 * t107;
	t103 = t113 * t130 + t126;
	t115 = sin(qJ(3));
	t117 = cos(qJ(3));
	t127 = t144 * t107;
	t129 = t118 * t111;
	t102 = -t113 * t129 + t127;
	t108 = sin(pkin(7));
	t112 = cos(pkin(7));
	t109 = sin(pkin(6));
	t132 = t109 * t118;
	t123 = t102 * t112 + t108 * t132;
	t92 = -t103 * t117 + t123 * t115;
	t99 = -t102 * t108 + t112 * t132;
	t146 = t92 * t114 - t99 * t116;
	t82 = t99 * t114 + t92 * t116;
	t122 = t113 * t126 + t130;
	t128 = t109 * t144;
	t145 = -t108 * t128 + t122 * t112;
	t101 = -t109 * t111 * t108 + t113 * t112;
	t131 = t111 * t112;
	t133 = t108 * t113;
	t98 = t115 * t133 + (t107 * t117 + t115 * t131) * t109;
	t87 = -t101 * t116 + t98 * t114;
	t78 = atan2(t146, t87);
	t75 = sin(t78);
	t76 = cos(t78);
	t69 = t146 * t75 + t76 * t87;
	t68 = 0.1e1 / t69 ^ 2;
	t119 = t122 * t108 + t112 * t128;
	t104 = -t113 * t127 + t129;
	t94 = t104 * t117 - t145 * t115;
	t83 = t94 * t114 - t119 * t116;
	t143 = t68 * t83;
	t110 = cos(pkin(13));
	t106 = sin(pkin(13));
	t93 = t104 * t115 + t145 * t117;
	t138 = t93 * t106;
	t84 = t119 * t114 + t94 * t116;
	t74 = t84 * t110 + t138;
	t72 = 0.1e1 / t74 ^ 2;
	t137 = t93 * t110;
	t73 = t84 * t106 - t137;
	t142 = t72 * t73;
	t141 = t76 * t146;
	t86 = 0.1e1 / t87 ^ 2;
	t140 = t146 * t86;
	t139 = t83 ^ 2 * t68;
	t124 = -t75 * t87 + t141;
	t120 = -t103 * t115 - t123 * t117;
	t97 = t117 * t133 + (-t107 * t115 + t117 * t131) * t109;
	t88 = t101 * t114 + t98 * t116;
	t85 = 0.1e1 / t87;
	t77 = 0.1e1 / (t146 ^ 2 * t86 + 0.1e1);
	t71 = 0.1e1 / t74;
	t70 = 0.1e1 / (t73 ^ 2 * t72 + 0.1e1);
	t67 = 0.1e1 / t69;
	t66 = 0.1e1 / (0.1e1 + t139);
	t65 = (-t120 * t85 - t97 * t140) * t77 * t114;
	t64 = (-t88 * t140 + t82 * t85) * t77;
	t1 = [-t83 * t85 * t77, 0, t65, t64, 0, 0; (t146 * t67 - (-t75 + (-t85 * t141 + t75) * t77) * t139) * t66, 0, (-t93 * t114 * t67 - (t124 * t65 + (-t120 * t75 + t76 * t97) * t114) * t143) * t66, (t84 * t67 - (t124 * t64 + t75 * t82 + t76 * t88) * t143) * t66, 0, 0; ((t82 * t106 - t110 * t120) * t71 - (t106 * t120 + t82 * t110) * t142) * t70, 0, ((-t94 * t110 - t116 * t138) * t71 - (t94 * t106 - t116 * t137) * t142) * t70, (-t106 * t71 + t110 * t142) * t83 * t70, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:21
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (1506->50), mult. (4121->118), div. (85->9), fcn. (5660->17), ass. (0->70)
	t133 = cos(pkin(6));
	t131 = cos(pkin(12));
	t164 = sin(qJ(1));
	t146 = t164 * t131;
	t128 = sin(pkin(12));
	t138 = cos(qJ(1));
	t151 = t138 * t128;
	t122 = t133 * t151 + t146;
	t135 = sin(qJ(3));
	t137 = cos(qJ(3));
	t147 = t164 * t128;
	t150 = t138 * t131;
	t121 = -t133 * t150 + t147;
	t129 = sin(pkin(7));
	t132 = cos(pkin(7));
	t130 = sin(pkin(6));
	t153 = t130 * t138;
	t143 = t121 * t132 + t129 * t153;
	t111 = -t122 * t137 + t143 * t135;
	t118 = -t121 * t129 + t132 * t153;
	t134 = sin(qJ(4));
	t136 = cos(qJ(4));
	t166 = t111 * t134 - t118 * t136;
	t101 = t111 * t136 + t118 * t134;
	t142 = t133 * t146 + t151;
	t148 = t130 * t164;
	t165 = -t129 * t148 + t142 * t132;
	t123 = -t133 * t147 + t150;
	t113 = t123 * t137 - t165 * t135;
	t139 = t142 * t129 + t132 * t148;
	t103 = t113 * t136 + t139 * t134;
	t112 = t123 * t135 + t165 * t137;
	t127 = pkin(13) + qJ(6);
	t125 = sin(t127);
	t126 = cos(t127);
	t93 = t103 * t126 + t112 * t125;
	t91 = 0.1e1 / t93 ^ 2;
	t92 = t103 * t125 - t112 * t126;
	t163 = t91 * t92;
	t152 = t131 * t132;
	t154 = t129 * t133;
	t117 = t135 * t154 + (t128 * t137 + t135 * t152) * t130;
	t120 = -t129 * t130 * t131 + t132 * t133;
	t106 = t117 * t134 - t120 * t136;
	t97 = atan2(t166, t106);
	t95 = cos(t97);
	t162 = t95 * t166;
	t102 = t113 * t134 - t139 * t136;
	t94 = sin(t97);
	t88 = t106 * t95 + t166 * t94;
	t87 = 0.1e1 / t88 ^ 2;
	t161 = t102 * t87;
	t160 = t102 ^ 2 * t87;
	t105 = 0.1e1 / t106 ^ 2;
	t159 = t105 * t166;
	t158 = t112 * t136;
	t149 = t91 * t92 ^ 2 + 0.1e1;
	t144 = -t106 * t94 + t162;
	t140 = -t122 * t135 - t143 * t137;
	t116 = t137 * t154 + (-t128 * t135 + t137 * t152) * t130;
	t107 = t117 * t136 + t120 * t134;
	t104 = 0.1e1 / t106;
	t96 = 0.1e1 / (t105 * t166 ^ 2 + 0.1e1);
	t90 = 0.1e1 / t93;
	t89 = 0.1e1 / t149;
	t86 = 0.1e1 / t88;
	t85 = 0.1e1 / (0.1e1 + t160);
	t84 = (-t104 * t140 - t116 * t159) * t96 * t134;
	t83 = (t101 * t104 - t107 * t159) * t96;
	t1 = [-t102 * t104 * t96, 0, t84, t83, 0, 0; (t166 * t86 - (-t94 + (-t104 * t162 + t94) * t96) * t160) * t85, 0, (-t112 * t134 * t86 - (t144 * t84 + (t116 * t95 - t140 * t94) * t134) * t161) * t85, (t103 * t86 - (t101 * t94 + t95 * t107 + t144 * t83) * t161) * t85, 0, 0; ((t101 * t125 - t126 * t140) * t90 - (t101 * t126 + t125 * t140) * t163) * t89, 0, ((-t113 * t126 - t125 * t158) * t90 - (t113 * t125 - t126 * t158) * t163) * t89, (-t125 * t90 + t126 * t163) * t89 * t102, 0, t149 * t89;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end