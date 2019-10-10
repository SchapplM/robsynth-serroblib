% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR9
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
%   Wie in S6RPRRPR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR9_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR9_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
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
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:32
	% DurationCPUTime: 0.13s
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
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:32
	% DurationCPUTime: 0.35s
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
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:32
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (607->36), mult. (1669->81), div. (55->9), fcn. (2293->15), ass. (0->56)
	t101 = cos(pkin(7));
	t105 = cos(qJ(1));
	t99 = sin(pkin(6));
	t117 = t105 * t99;
	t102 = cos(pkin(6));
	t124 = sin(qJ(1));
	t97 = sin(pkin(12));
	t111 = t124 * t97;
	t100 = cos(pkin(12));
	t115 = t105 * t100;
	t89 = -t102 * t115 + t111;
	t98 = sin(pkin(7));
	t125 = t101 * t89 + t98 * t117;
	t103 = sin(qJ(3));
	t104 = cos(qJ(3));
	t110 = t124 * t100;
	t118 = t105 * t97;
	t90 = t102 * t118 + t110;
	t75 = t90 * t103 + t104 * t125;
	t107 = t102 * t110 + t118;
	t112 = t99 * t124;
	t126 = t107 * t101 - t98 * t112;
	t91 = -t102 * t111 + t115;
	t80 = -t103 * t126 + t91 * t104;
	t86 = t101 * t112 + t107 * t98;
	t96 = qJ(4) + pkin(13);
	t94 = sin(t96);
	t95 = cos(t96);
	t70 = t80 * t95 + t86 * t94;
	t68 = 0.1e1 / t70 ^ 2;
	t69 = t80 * t94 - t86 * t95;
	t123 = t68 * t69;
	t108 = t100 * t101 * t99 + t102 * t98;
	t120 = t97 * t99;
	t83 = t103 * t120 - t104 * t108;
	t74 = atan2(-t75, t83);
	t72 = cos(t74);
	t122 = t72 * t75;
	t71 = sin(t74);
	t65 = -t71 * t75 + t72 * t83;
	t64 = 0.1e1 / t65 ^ 2;
	t79 = t91 * t103 + t104 * t126;
	t121 = t79 ^ 2 * t64;
	t113 = t68 * t69 ^ 2 + 0.1e1;
	t78 = t103 * t125 - t90 * t104;
	t85 = t101 * t117 - t89 * t98;
	t84 = t103 * t108 + t104 * t120;
	t82 = 0.1e1 / t83 ^ 2;
	t81 = 0.1e1 / t83;
	t73 = 0.1e1 / (t75 ^ 2 * t82 + 0.1e1);
	t67 = 0.1e1 / t70;
	t66 = 0.1e1 / t113;
	t63 = 0.1e1 / t65;
	t62 = 0.1e1 / (0.1e1 + t121);
	t61 = (t75 * t82 * t84 + t78 * t81) * t73;
	t1 = [-t79 * t81 * t73, 0, t61, 0, 0, 0; (-t75 * t63 - (-t71 + (t122 * t81 + t71) * t73) * t121) * t62, 0, (t80 * t63 - (t71 * t78 + t72 * t84 + (-t71 * t83 - t122) * t61) * t79 * t64) * t62, 0, 0, 0; ((t78 * t94 - t85 * t95) * t67 - (t78 * t95 + t85 * t94) * t123) * t66, 0, (t123 * t95 - t67 * t94) * t79 * t66, t113 * t66, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:32
	% DurationCPUTime: 0.67s
	% Computational Cost: add. (1859->50), mult. (4121->117), div. (85->9), fcn. (5660->17), ass. (0->71)
	t127 = cos(pkin(6));
	t125 = cos(pkin(12));
	t159 = sin(qJ(1));
	t140 = t159 * t125;
	t122 = sin(pkin(12));
	t132 = cos(qJ(1));
	t145 = t132 * t122;
	t116 = t127 * t145 + t140;
	t129 = sin(qJ(3));
	t131 = cos(qJ(3));
	t141 = t159 * t122;
	t144 = t132 * t125;
	t115 = -t127 * t144 + t141;
	t123 = sin(pkin(7));
	t126 = cos(pkin(7));
	t124 = sin(pkin(6));
	t147 = t124 * t132;
	t137 = t115 * t126 + t123 * t147;
	t105 = -t116 * t131 + t137 * t129;
	t112 = -t115 * t123 + t126 * t147;
	t121 = qJ(4) + pkin(13);
	t119 = sin(t121);
	t120 = cos(t121);
	t161 = t105 * t119 - t112 * t120;
	t95 = t105 * t120 + t112 * t119;
	t136 = t127 * t140 + t145;
	t142 = t124 * t159;
	t160 = -t123 * t142 + t136 * t126;
	t146 = t125 * t126;
	t148 = t123 * t127;
	t111 = t129 * t148 + (t122 * t131 + t129 * t146) * t124;
	t114 = -t124 * t125 * t123 + t127 * t126;
	t100 = t111 * t119 - t114 * t120;
	t89 = atan2(t161, t100);
	t84 = sin(t89);
	t85 = cos(t89);
	t82 = t85 * t100 + t161 * t84;
	t81 = 0.1e1 / t82 ^ 2;
	t117 = -t127 * t141 + t144;
	t107 = t117 * t131 - t160 * t129;
	t133 = t136 * t123 + t126 * t142;
	t96 = t107 * t119 - t133 * t120;
	t158 = t81 * t96;
	t157 = t85 * t161;
	t130 = cos(qJ(6));
	t106 = t117 * t129 + t160 * t131;
	t128 = sin(qJ(6));
	t153 = t106 * t128;
	t97 = t107 * t120 + t133 * t119;
	t91 = t97 * t130 + t153;
	t88 = 0.1e1 / t91 ^ 2;
	t152 = t106 * t130;
	t90 = t97 * t128 - t152;
	t156 = t88 * t90;
	t99 = 0.1e1 / t100 ^ 2;
	t155 = t161 * t99;
	t154 = t96 ^ 2 * t81;
	t143 = t90 ^ 2 * t88 + 0.1e1;
	t138 = -t100 * t84 + t157;
	t134 = -t116 * t129 - t137 * t131;
	t110 = t131 * t148 + (-t122 * t129 + t131 * t146) * t124;
	t101 = t111 * t120 + t114 * t119;
	t98 = 0.1e1 / t100;
	t87 = 0.1e1 / t91;
	t86 = 0.1e1 / (t161 ^ 2 * t99 + 0.1e1);
	t83 = 0.1e1 / t143;
	t80 = 0.1e1 / t82;
	t79 = 0.1e1 / (0.1e1 + t154);
	t78 = (-t110 * t155 - t134 * t98) * t86 * t119;
	t77 = (-t101 * t155 + t95 * t98) * t86;
	t1 = [-t96 * t98 * t86, 0, t78, t77, 0, 0; (t161 * t80 - (-t84 + (-t98 * t157 + t84) * t86) * t154) * t79, 0, (-t106 * t119 * t80 - (t138 * t78 + (t110 * t85 - t134 * t84) * t119) * t158) * t79, (t97 * t80 - (t85 * t101 + t138 * t77 + t84 * t95) * t158) * t79, 0, 0; ((t95 * t128 - t130 * t134) * t87 - (t128 * t134 + t95 * t130) * t156) * t83, 0, ((-t107 * t130 - t120 * t153) * t87 - (t107 * t128 - t120 * t152) * t156) * t83, (-t128 * t87 + t130 * t156) * t96 * t83, 0, t143 * t83;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end