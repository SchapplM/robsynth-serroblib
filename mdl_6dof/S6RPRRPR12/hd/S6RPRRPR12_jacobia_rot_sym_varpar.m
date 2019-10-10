% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR12
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
%   Wie in S6RPRRPR12_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:43
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR12_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR12_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
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
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.22s
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
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:26
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
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:26
	% DurationCPUTime: 0.53s
	% Computational Cost: add. (1158->43), mult. (3332->103), div. (77->9), fcn. (4592->15), ass. (0->63)
	t106 = sin(qJ(4));
	t109 = cos(qJ(4));
	t107 = sin(qJ(3));
	t110 = cos(qJ(3));
	t101 = sin(pkin(7));
	t104 = cos(pkin(7));
	t102 = sin(pkin(6));
	t111 = cos(qJ(1));
	t123 = t102 * t111;
	t105 = cos(pkin(6));
	t103 = cos(pkin(12));
	t118 = t111 * t103;
	t100 = sin(pkin(12));
	t108 = sin(qJ(1));
	t121 = t108 * t100;
	t96 = -t105 * t118 + t121;
	t115 = t101 * t123 + t104 * t96;
	t119 = t111 * t100;
	t120 = t108 * t103;
	t97 = t105 * t119 + t120;
	t86 = t115 * t107 - t97 * t110;
	t93 = -t96 * t101 + t104 * t123;
	t136 = t86 * t106 - t93 * t109;
	t135 = t93 * t106 + t86 * t109;
	t114 = t105 * t120 + t119;
	t124 = t102 * t108;
	t134 = -t101 * t124 + t114 * t104;
	t122 = t103 * t104;
	t125 = t101 * t105;
	t92 = t107 * t125 + (t100 * t110 + t107 * t122) * t102;
	t95 = -t102 * t103 * t101 + t105 * t104;
	t80 = t92 * t106 - t95 * t109;
	t71 = atan2(t136, t80);
	t68 = sin(t71);
	t69 = cos(t71);
	t67 = t136 * t68 + t69 * t80;
	t66 = 0.1e1 / t67 ^ 2;
	t112 = t114 * t101 + t104 * t124;
	t98 = -t105 * t121 + t118;
	t88 = -t134 * t107 + t98 * t110;
	t76 = t88 * t106 - t112 * t109;
	t133 = t66 * t76;
	t132 = t69 * t136;
	t79 = 0.1e1 / t80 ^ 2;
	t131 = t136 * t79;
	t130 = t76 ^ 2 * t66;
	t77 = t112 * t106 + t88 * t109;
	t87 = t98 * t107 + t134 * t110;
	t83 = 0.1e1 / t87 ^ 2;
	t129 = t77 * t83;
	t116 = -t68 * t80 + t132;
	t84 = -t97 * t107 - t115 * t110;
	t91 = t110 * t125 + (-t100 * t107 + t110 * t122) * t102;
	t82 = 0.1e1 / t87;
	t81 = t95 * t106 + t92 * t109;
	t78 = 0.1e1 / t80;
	t72 = 0.1e1 / (t77 ^ 2 * t83 + 0.1e1);
	t70 = 0.1e1 / (t136 ^ 2 * t79 + 0.1e1);
	t65 = 0.1e1 / t67;
	t64 = 0.1e1 / (0.1e1 + t130);
	t63 = (-t91 * t131 - t78 * t84) * t70 * t106;
	t62 = (-t81 * t131 + t135 * t78) * t70;
	t1 = [-t76 * t78 * t70, 0, t63, t62, 0, 0; (t136 * t65 - (-t68 + (-t78 * t132 + t68) * t70) * t130) * t64, 0, (-t87 * t106 * t65 - (t116 * t63 + (-t68 * t84 + t69 * t91) * t106) * t133) * t64, (t77 * t65 - (t116 * t62 + t135 * t68 + t69 * t81) * t133) * t64, 0, 0; (-t84 * t129 + t135 * t82) * t72, 0, (-t109 * t82 * t87 - t88 * t129) * t72, -t76 * t82 * t72, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:26
	% DurationCPUTime: 0.60s
	% Computational Cost: add. (1440->49), mult. (4121->117), div. (85->9), fcn. (5660->17), ass. (0->70)
	t127 = cos(pkin(6));
	t122 = sin(pkin(12));
	t135 = cos(qJ(1));
	t144 = t135 * t122;
	t125 = cos(pkin(12));
	t131 = sin(qJ(1));
	t145 = t131 * t125;
	t118 = t127 * t144 + t145;
	t130 = sin(qJ(3));
	t134 = cos(qJ(3));
	t143 = t135 * t125;
	t146 = t131 * t122;
	t117 = -t127 * t143 + t146;
	t123 = sin(pkin(7));
	t126 = cos(pkin(7));
	t124 = sin(pkin(6));
	t148 = t124 * t135;
	t139 = t117 * t126 + t123 * t148;
	t107 = -t118 * t134 + t139 * t130;
	t112 = -t117 * t123 + t126 * t148;
	t129 = sin(qJ(4));
	t133 = cos(qJ(4));
	t97 = t107 * t129 - t112 * t133;
	t162 = t107 * t133 + t112 * t129;
	t138 = t127 * t145 + t144;
	t149 = t124 * t131;
	t159 = -t123 * t149 + t138 * t126;
	t147 = t125 * t126;
	t150 = t123 * t127;
	t111 = t130 * t150 + (t122 * t134 + t130 * t147) * t124;
	t116 = -t124 * t125 * t123 + t127 * t126;
	t103 = t111 * t133 + t116 * t129;
	t93 = atan2(t162, t103);
	t90 = sin(t93);
	t91 = cos(t93);
	t84 = t91 * t103 + t162 * t90;
	t83 = 0.1e1 / t84 ^ 2;
	t119 = -t127 * t146 + t143;
	t109 = t119 * t134 - t159 * t130;
	t114 = t138 * t123 + t126 * t149;
	t99 = t109 * t133 + t114 * t129;
	t158 = t83 * t99;
	t128 = sin(qJ(6));
	t108 = t119 * t130 + t159 * t134;
	t132 = cos(qJ(6));
	t152 = t108 * t132;
	t98 = t109 * t129 - t114 * t133;
	t89 = t98 * t128 + t152;
	t87 = 0.1e1 / t89 ^ 2;
	t153 = t108 * t128;
	t88 = -t98 * t132 + t153;
	t157 = t87 * t88;
	t156 = t91 * t162;
	t155 = t99 ^ 2 * t83;
	t101 = 0.1e1 / t103 ^ 2;
	t154 = t101 * t162;
	t142 = t88 ^ 2 * t87 + 0.1e1;
	t140 = -t103 * t90 + t156;
	t136 = -t118 * t130 - t139 * t134;
	t110 = t134 * t150 + (-t122 * t130 + t134 * t147) * t124;
	t102 = -t111 * t129 + t116 * t133;
	t100 = 0.1e1 / t103;
	t92 = 0.1e1 / (t101 * t162 ^ 2 + 0.1e1);
	t86 = 0.1e1 / t89;
	t85 = 0.1e1 / t142;
	t82 = 0.1e1 / t84;
	t81 = 0.1e1 / (0.1e1 + t155);
	t80 = (-t100 * t136 - t110 * t154) * t92 * t133;
	t79 = (-t100 * t97 - t102 * t154) * t92;
	t1 = [-t99 * t100 * t92, 0, t80, t79, 0, 0; (t162 * t82 - (-t90 + (-t100 * t156 + t90) * t92) * t155) * t81, 0, (-t108 * t133 * t82 - (t140 * t80 + (t110 * t91 - t136 * t90) * t133) * t158) * t81, (-t98 * t82 - (t91 * t102 + t140 * t79 - t90 * t97) * t158) * t81, 0, 0; ((t128 * t136 - t97 * t132) * t86 - (t97 * t128 + t132 * t136) * t157) * t85, 0, ((t109 * t128 + t129 * t152) * t86 - (t109 * t132 - t129 * t153) * t157) * t85, (-t128 * t157 - t132 * t86) * t99 * t85, 0, t142 * t85;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end