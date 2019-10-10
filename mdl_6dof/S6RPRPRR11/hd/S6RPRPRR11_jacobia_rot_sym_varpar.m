% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR11
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
%   Wie in S6RPRPRR11_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR11_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR11_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
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
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
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
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (507->39), mult. (1523->86), div. (50->9), fcn. (2097->15), ass. (0->56)
	t84 = sin(pkin(7));
	t88 = cos(pkin(7));
	t83 = sin(pkin(12));
	t92 = cos(qJ(1));
	t101 = t92 * t83;
	t89 = cos(pkin(6));
	t110 = sin(qJ(1));
	t87 = cos(pkin(12));
	t96 = t110 * t87;
	t94 = t89 * t96 + t101;
	t85 = sin(pkin(6));
	t98 = t85 * t110;
	t111 = -t84 * t98 + t94 * t88;
	t100 = t92 * t87;
	t97 = t110 * t83;
	t79 = -t89 * t97 + t100;
	t90 = sin(qJ(3));
	t91 = cos(qJ(3));
	t68 = -t111 * t90 + t79 * t91;
	t74 = t94 * t84 + t88 * t98;
	t82 = sin(pkin(13));
	t86 = cos(pkin(13));
	t58 = t68 * t86 + t74 * t82;
	t56 = 0.1e1 / t58 ^ 2;
	t57 = t68 * t82 - t74 * t86;
	t109 = t56 * t57;
	t102 = t88 * t91;
	t78 = t89 * t101 + t96;
	t106 = t78 * t90;
	t77 = -t89 * t100 + t97;
	t104 = t85 * t92;
	t99 = t84 * t104;
	t63 = t77 * t102 + t91 * t99 + t106;
	t105 = t84 * t89;
	t71 = -t91 * t105 + (-t87 * t102 + t83 * t90) * t85;
	t62 = atan2(-t63, t71);
	t60 = cos(t62);
	t108 = t60 * t63;
	t59 = sin(t62);
	t53 = -t59 * t63 + t60 * t71;
	t52 = 0.1e1 / t53 ^ 2;
	t67 = t111 * t91 + t79 * t90;
	t107 = t67 ^ 2 * t52;
	t103 = t88 * t90;
	t66 = t77 * t103 - t78 * t91 + t90 * t99;
	t73 = t88 * t104 - t77 * t84;
	t72 = t90 * t105 + (t87 * t103 + t83 * t91) * t85;
	t70 = 0.1e1 / t71 ^ 2;
	t69 = 0.1e1 / t71;
	t61 = 0.1e1 / (t63 ^ 2 * t70 + 0.1e1);
	t55 = 0.1e1 / t58;
	t54 = 0.1e1 / (t57 ^ 2 * t56 + 0.1e1);
	t51 = 0.1e1 / t53;
	t50 = 0.1e1 / (0.1e1 + t107);
	t49 = (t63 * t70 * t72 + t66 * t69) * t61;
	t1 = [-t67 * t69 * t61, 0, t49, 0, 0, 0; ((-t106 + (-t77 * t88 - t99) * t91) * t51 - (-t59 + (t69 * t108 + t59) * t61) * t107) * t50, 0, (t68 * t51 - (t59 * t66 + t60 * t72 + (-t59 * t71 - t108) * t49) * t67 * t52) * t50, 0, 0, 0; ((t66 * t82 - t73 * t86) * t55 - (t66 * t86 + t73 * t82) * t109) * t54, 0, (t86 * t109 - t82 * t55) * t67 * t54, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (607->36), mult. (1669->83), div. (55->9), fcn. (2293->15), ass. (0->54)
	t100 = cos(pkin(7));
	t104 = cos(qJ(1));
	t98 = sin(pkin(6));
	t115 = t104 * t98;
	t101 = cos(pkin(6));
	t113 = t101 * t104;
	t121 = sin(qJ(1));
	t96 = sin(pkin(12));
	t99 = cos(pkin(12));
	t88 = -t99 * t113 + t121 * t96;
	t97 = sin(pkin(7));
	t122 = t100 * t88 + t97 * t115;
	t102 = sin(qJ(3));
	t103 = cos(qJ(3));
	t89 = t96 * t113 + t121 * t99;
	t74 = t89 * t102 + t122 * t103;
	t109 = t101 * t121;
	t106 = t104 * t96 + t99 * t109;
	t110 = t98 * t121;
	t123 = t106 * t100 - t97 * t110;
	t90 = t104 * t99 - t96 * t109;
	t79 = -t123 * t102 + t90 * t103;
	t85 = t100 * t110 + t106 * t97;
	t95 = pkin(13) + qJ(5);
	t93 = sin(t95);
	t94 = cos(t95);
	t69 = t79 * t94 + t85 * t93;
	t67 = 0.1e1 / t69 ^ 2;
	t68 = t79 * t93 - t85 * t94;
	t120 = t67 * t68;
	t107 = t100 * t98 * t99 + t101 * t97;
	t117 = t96 * t98;
	t82 = t102 * t117 - t107 * t103;
	t73 = atan2(-t74, t82);
	t71 = cos(t73);
	t119 = t71 * t74;
	t70 = sin(t73);
	t64 = -t70 * t74 + t71 * t82;
	t63 = 0.1e1 / t64 ^ 2;
	t78 = t90 * t102 + t123 * t103;
	t118 = t78 ^ 2 * t63;
	t111 = t68 ^ 2 * t67 + 0.1e1;
	t77 = t122 * t102 - t89 * t103;
	t84 = t100 * t115 - t88 * t97;
	t83 = t107 * t102 + t103 * t117;
	t81 = 0.1e1 / t82 ^ 2;
	t80 = 0.1e1 / t82;
	t72 = 0.1e1 / (t74 ^ 2 * t81 + 0.1e1);
	t66 = 0.1e1 / t69;
	t65 = 0.1e1 / t111;
	t62 = 0.1e1 / t64;
	t61 = 0.1e1 / (0.1e1 + t118);
	t60 = (t74 * t81 * t83 + t77 * t80) * t72;
	t1 = [-t78 * t80 * t72, 0, t60, 0, 0, 0; (-t74 * t62 - (-t70 + (t80 * t119 + t70) * t72) * t118) * t61, 0, (t79 * t62 - (t70 * t77 + t71 * t83 + (-t70 * t82 - t119) * t60) * t78 * t63) * t61, 0, 0, 0; ((t77 * t93 - t84 * t94) * t66 - (t77 * t94 + t84 * t93) * t120) * t65, 0, (t94 * t120 - t93 * t66) * t78 * t65, 0, t111 * t65, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:01
	% DurationCPUTime: 0.67s
	% Computational Cost: add. (1859->50), mult. (4121->117), div. (85->9), fcn. (5660->17), ass. (0->71)
	t128 = cos(pkin(6));
	t126 = cos(pkin(12));
	t160 = sin(qJ(1));
	t141 = t160 * t126;
	t123 = sin(pkin(12));
	t133 = cos(qJ(1));
	t146 = t133 * t123;
	t117 = t128 * t146 + t141;
	t130 = sin(qJ(3));
	t132 = cos(qJ(3));
	t142 = t160 * t123;
	t145 = t133 * t126;
	t116 = -t128 * t145 + t142;
	t124 = sin(pkin(7));
	t127 = cos(pkin(7));
	t125 = sin(pkin(6));
	t148 = t125 * t133;
	t138 = t116 * t127 + t124 * t148;
	t106 = -t117 * t132 + t138 * t130;
	t113 = -t116 * t124 + t127 * t148;
	t122 = pkin(13) + qJ(5);
	t120 = sin(t122);
	t121 = cos(t122);
	t162 = t106 * t120 - t113 * t121;
	t96 = t106 * t121 + t113 * t120;
	t137 = t128 * t141 + t146;
	t143 = t125 * t160;
	t161 = -t124 * t143 + t137 * t127;
	t147 = t126 * t127;
	t149 = t124 * t128;
	t112 = t130 * t149 + (t123 * t132 + t130 * t147) * t125;
	t115 = -t125 * t126 * t124 + t128 * t127;
	t101 = t112 * t120 - t115 * t121;
	t90 = atan2(t162, t101);
	t85 = sin(t90);
	t86 = cos(t90);
	t83 = t86 * t101 + t162 * t85;
	t82 = 0.1e1 / t83 ^ 2;
	t118 = -t128 * t142 + t145;
	t108 = t118 * t132 - t161 * t130;
	t134 = t137 * t124 + t127 * t143;
	t97 = t108 * t120 - t134 * t121;
	t159 = t82 * t97;
	t158 = t86 * t162;
	t131 = cos(qJ(6));
	t107 = t118 * t130 + t161 * t132;
	t129 = sin(qJ(6));
	t154 = t107 * t129;
	t98 = t108 * t121 + t134 * t120;
	t92 = t98 * t131 + t154;
	t89 = 0.1e1 / t92 ^ 2;
	t153 = t107 * t131;
	t91 = t98 * t129 - t153;
	t157 = t89 * t91;
	t156 = t97 ^ 2 * t82;
	t100 = 0.1e1 / t101 ^ 2;
	t155 = t100 * t162;
	t144 = t91 ^ 2 * t89 + 0.1e1;
	t139 = -t101 * t85 + t158;
	t135 = -t117 * t130 - t138 * t132;
	t111 = t132 * t149 + (-t123 * t130 + t132 * t147) * t125;
	t102 = t112 * t121 + t115 * t120;
	t99 = 0.1e1 / t101;
	t88 = 0.1e1 / t92;
	t87 = 0.1e1 / (t100 * t162 ^ 2 + 0.1e1);
	t84 = 0.1e1 / t144;
	t81 = 0.1e1 / t83;
	t80 = 0.1e1 / (0.1e1 + t156);
	t79 = (-t111 * t155 - t135 * t99) * t87 * t120;
	t78 = (-t102 * t155 + t96 * t99) * t87;
	t1 = [-t97 * t99 * t87, 0, t79, 0, t78, 0; (t162 * t81 - (-t85 + (-t99 * t158 + t85) * t87) * t156) * t80, 0, (-t107 * t120 * t81 - (t139 * t79 + (t111 * t86 - t135 * t85) * t120) * t159) * t80, 0, (t98 * t81 - (t86 * t102 + t139 * t78 + t85 * t96) * t159) * t80, 0; ((t96 * t129 - t131 * t135) * t88 - (t129 * t135 + t96 * t131) * t157) * t84, 0, ((-t108 * t131 - t121 * t154) * t88 - (t108 * t129 - t121 * t153) * t157) * t84, 0, (-t129 * t88 + t131 * t157) * t97 * t84, t144 * t84;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end