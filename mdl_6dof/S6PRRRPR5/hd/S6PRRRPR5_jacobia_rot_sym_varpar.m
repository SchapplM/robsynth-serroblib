% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR5
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
%   Wie in S6PRRRPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:54
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:08
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (179->22), mult. (525->56), div. (35->9), fcn. (736->13), ass. (0->37)
	t45 = sin(pkin(12));
	t48 = cos(pkin(12));
	t54 = cos(qJ(2));
	t50 = cos(pkin(6));
	t52 = sin(qJ(2));
	t58 = t50 * t52;
	t43 = -t45 * t58 + t48 * t54;
	t51 = sin(qJ(3));
	t63 = t43 * t51;
	t53 = cos(qJ(3));
	t62 = t43 * t53;
	t46 = sin(pkin(7));
	t47 = sin(pkin(6));
	t61 = t46 * t47;
	t49 = cos(pkin(7));
	t60 = t47 * t49;
	t59 = t47 * t52;
	t57 = t50 * t54;
	t42 = -t45 * t57 - t48 * t52;
	t55 = t42 * t49 + t45 * t61;
	t32 = t51 * t55 + t62;
	t30 = 0.1e1 / t32 ^ 2;
	t31 = -t53 * t55 + t63;
	t56 = t31 ^ 2 * t30 + 0.1e1;
	t41 = -t45 * t54 - t48 * t58;
	t40 = t50 * t49 - t54 * t61;
	t39 = 0.1e1 / t40 ^ 2;
	t38 = -t42 * t46 + t45 * t60;
	t37 = (-t45 * t52 + t48 * t57) * t46 + t48 * t60;
	t36 = atan2(t37, t40);
	t34 = cos(t36);
	t33 = sin(t36);
	t29 = 0.1e1 / t56;
	t28 = t33 * t37 + t34 * t40;
	t27 = 0.1e1 / t28 ^ 2;
	t25 = (t41 / t40 - t37 * t39 * t59) * t46 / (t37 ^ 2 * t39 + 0.1e1);
	t1 = [0, t25, 0, 0, 0, 0; 0, (t43 * t46 / t28 - ((t33 * t41 + t34 * t59) * t46 + (-t33 * t40 + t34 * t37) * t25) * t38 * t27) / (t38 ^ 2 * t27 + 0.1e1), 0, 0, 0, 0; 0, ((t42 * t51 + t49 * t62) / t32 - (t42 * t53 - t49 * t63) * t31 * t30) * t29, t56 * t29, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:08
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (607->40), mult. (1825->97), div. (65->9), fcn. (2498->15), ass. (0->60)
	t87 = sin(pkin(7));
	t88 = sin(pkin(6));
	t109 = t88 * t87;
	t89 = cos(pkin(12));
	t90 = cos(pkin(7));
	t91 = cos(pkin(6));
	t97 = cos(qJ(2));
	t106 = t91 * t97;
	t86 = sin(pkin(12));
	t94 = sin(qJ(2));
	t99 = t89 * t106 - t86 * t94;
	t116 = -t89 * t109 + t99 * t90;
	t107 = t91 * t94;
	t81 = t89 * t107 + t86 * t97;
	t93 = sin(qJ(3));
	t96 = cos(qJ(3));
	t66 = -t116 * t96 + t81 * t93;
	t108 = t90 * t96;
	t110 = t87 * t91;
	t75 = -t96 * t110 + (-t97 * t108 + t93 * t94) * t88;
	t65 = atan2(-t66, t75);
	t62 = sin(t65);
	t63 = cos(t65);
	t56 = -t62 * t66 + t63 * t75;
	t55 = 0.1e1 / t56 ^ 2;
	t103 = t86 * t109;
	t83 = -t86 * t107 + t89 * t97;
	t111 = t83 * t93;
	t82 = -t86 * t106 - t89 * t94;
	t69 = -t96 * t103 - t82 * t108 + t111;
	t115 = t55 * t69;
	t70 = t83 * t96 + (t82 * t90 + t103) * t93;
	t77 = t86 * t88 * t90 - t82 * t87;
	t92 = sin(qJ(4));
	t95 = cos(qJ(4));
	t61 = t70 * t95 + t77 * t92;
	t59 = 0.1e1 / t61 ^ 2;
	t60 = t70 * t92 - t77 * t95;
	t114 = t59 * t60;
	t74 = 0.1e1 / t75 ^ 2;
	t113 = t66 * t74;
	t112 = t83 * t87;
	t105 = t93 * t97;
	t104 = t94 * t96;
	t101 = t60 ^ 2 * t59 + 0.1e1;
	t100 = -t62 * t75 - t63 * t66;
	t80 = (t90 * t104 + t105) * t88;
	t76 = t93 * t110 + (t90 * t105 + t104) * t88;
	t73 = 0.1e1 / t75;
	t72 = -t90 * t111 + t82 * t96;
	t71 = t81 * t108 + t99 * t93;
	t68 = t116 * t93 + t81 * t96;
	t64 = 0.1e1 / (t66 ^ 2 * t74 + 0.1e1);
	t58 = 0.1e1 / t61;
	t57 = 0.1e1 / t101;
	t54 = 0.1e1 / t56;
	t53 = 0.1e1 / (t69 ^ 2 * t55 + 0.1e1);
	t52 = (t80 * t113 - t71 * t73) * t64;
	t51 = (t76 * t113 - t68 * t73) * t64;
	t1 = [0, t52, t51, 0, 0, 0; 0, ((t83 * t108 + t82 * t93) * t54 - (t100 * t52 - t62 * t71 + t63 * t80) * t115) * t53, (t70 * t54 - (t100 * t51 - t62 * t68 + t63 * t76) * t115) * t53, 0, 0, 0; 0, ((-t95 * t112 + t72 * t92) * t58 - (t92 * t112 + t72 * t95) * t114) * t57, (t95 * t114 - t92 * t58) * t69 * t57, t101 * t57, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:08
	% EndTime: 2019-10-09 22:54:08
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (655->41), mult. (1825->97), div. (65->9), fcn. (2498->15), ass. (0->61)
	t101 = cos(pkin(12));
	t102 = cos(pkin(7));
	t105 = sin(qJ(2));
	t103 = cos(pkin(6));
	t107 = cos(qJ(2));
	t116 = t103 * t107;
	t98 = sin(pkin(12));
	t109 = t101 * t116 - t98 * t105;
	t100 = sin(pkin(6));
	t99 = sin(pkin(7));
	t121 = t100 * t99;
	t126 = -t101 * t121 + t109 * t102;
	t104 = sin(qJ(3));
	t106 = cos(qJ(3));
	t117 = t103 * t105;
	t90 = t101 * t117 + t98 * t107;
	t75 = t90 * t104 - t126 * t106;
	t118 = t102 * t106;
	t120 = t103 * t99;
	t84 = -t106 * t120 + (t104 * t105 - t107 * t118) * t100;
	t74 = atan2(-t75, t84);
	t71 = sin(t74);
	t72 = cos(t74);
	t65 = -t71 * t75 + t72 * t84;
	t64 = 0.1e1 / t65 ^ 2;
	t113 = t98 * t121;
	t92 = t101 * t107 - t98 * t117;
	t119 = t92 * t104;
	t91 = -t101 * t105 - t98 * t116;
	t78 = -t106 * t113 - t91 * t118 + t119;
	t125 = t64 * t78;
	t79 = t92 * t106 + (t102 * t91 + t113) * t104;
	t86 = t98 * t100 * t102 - t91 * t99;
	t97 = qJ(4) + pkin(13);
	t95 = sin(t97);
	t96 = cos(t97);
	t70 = t79 * t96 + t86 * t95;
	t68 = 0.1e1 / t70 ^ 2;
	t69 = t79 * t95 - t86 * t96;
	t124 = t68 * t69;
	t83 = 0.1e1 / t84 ^ 2;
	t123 = t75 * t83;
	t122 = t92 * t99;
	t115 = t104 * t107;
	t114 = t105 * t106;
	t111 = t69 ^ 2 * t68 + 0.1e1;
	t110 = -t71 * t84 - t72 * t75;
	t89 = (t102 * t114 + t115) * t100;
	t85 = t104 * t120 + (t102 * t115 + t114) * t100;
	t82 = 0.1e1 / t84;
	t81 = -t102 * t119 + t91 * t106;
	t80 = t109 * t104 + t90 * t118;
	t77 = t126 * t104 + t90 * t106;
	t73 = 0.1e1 / (t75 ^ 2 * t83 + 0.1e1);
	t67 = 0.1e1 / t70;
	t66 = 0.1e1 / t111;
	t63 = 0.1e1 / t65;
	t62 = 0.1e1 / (t78 ^ 2 * t64 + 0.1e1);
	t61 = (t89 * t123 - t80 * t82) * t73;
	t60 = (t85 * t123 - t77 * t82) * t73;
	t1 = [0, t61, t60, 0, 0, 0; 0, ((t91 * t104 + t92 * t118) * t63 - (t110 * t61 - t71 * t80 + t72 * t89) * t125) * t62, (t79 * t63 - (t110 * t60 - t71 * t77 + t72 * t85) * t125) * t62, 0, 0, 0; 0, ((-t96 * t122 + t81 * t95) * t67 - (t95 * t122 + t81 * t96) * t124) * t66, (t96 * t124 - t95 * t67) * t78 * t66, t111 * t66, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:08
	% EndTime: 2019-10-09 22:54:08
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (1983->60), mult. (4378->142), div. (95->9), fcn. (6002->17), ass. (0->77)
	t130 = sin(pkin(6));
	t135 = sin(qJ(3));
	t136 = sin(qJ(2));
	t138 = cos(qJ(3));
	t139 = cos(qJ(2));
	t132 = cos(pkin(7));
	t152 = t132 * t135;
	t129 = sin(pkin(7));
	t133 = cos(pkin(6));
	t155 = t129 * t133;
	t117 = t135 * t155 + (t136 * t138 + t139 * t152) * t130;
	t154 = t130 * t129;
	t119 = t133 * t132 - t139 * t154;
	t127 = qJ(4) + pkin(13);
	t125 = sin(t127);
	t126 = cos(t127);
	t105 = t117 * t125 - t119 * t126;
	t128 = sin(pkin(12));
	t131 = cos(pkin(12));
	t150 = t133 * t136;
	t120 = t128 * t139 + t131 * t150;
	t149 = t133 * t139;
	t142 = -t128 * t136 + t131 * t149;
	t140 = -t131 * t154 + t142 * t132;
	t108 = t120 * t138 + t140 * t135;
	t153 = t130 * t132;
	t141 = -t142 * t129 - t131 * t153;
	t96 = t108 * t125 - t141 * t126;
	t95 = atan2(-t96, t105);
	t90 = sin(t95);
	t91 = cos(t95);
	t86 = t91 * t105 - t90 * t96;
	t85 = 0.1e1 / t86 ^ 2;
	t121 = -t128 * t149 - t131 * t136;
	t122 = -t128 * t150 + t131 * t139;
	t145 = t128 * t154;
	t110 = t122 * t138 + (t121 * t132 + t145) * t135;
	t143 = -t121 * t129 + t128 * t153;
	t99 = t110 * t125 - t143 * t126;
	t161 = t85 * t99;
	t100 = t110 * t126 + t143 * t125;
	t137 = cos(qJ(6));
	t151 = t132 * t138;
	t109 = -t121 * t151 + t122 * t135 - t138 * t145;
	t134 = sin(qJ(6));
	t158 = t109 * t134;
	t93 = t100 * t137 + t158;
	t89 = 0.1e1 / t93 ^ 2;
	t157 = t109 * t137;
	t92 = t100 * t134 - t157;
	t160 = t89 * t92;
	t104 = 0.1e1 / t105 ^ 2;
	t159 = t104 * t96;
	t156 = t126 * t129;
	t148 = t135 * t136;
	t147 = t138 * t139;
	t146 = t92 ^ 2 * t89 + 0.1e1;
	t144 = -t105 * t90 - t91 * t96;
	t116 = t138 * t155 + (t132 * t147 - t148) * t130;
	t113 = ((-t132 * t148 + t147) * t125 - t136 * t156) * t130;
	t112 = t121 * t138 - t122 * t152;
	t111 = t121 * t135 + t122 * t151;
	t107 = -t120 * t135 + t140 * t138;
	t106 = t117 * t126 + t119 * t125;
	t103 = 0.1e1 / t105;
	t102 = t122 * t129 * t125 + t112 * t126;
	t101 = (-t120 * t152 + t142 * t138) * t125 - t120 * t156;
	t98 = t108 * t126 + t141 * t125;
	t94 = 0.1e1 / (t96 ^ 2 * t104 + 0.1e1);
	t88 = 0.1e1 / t93;
	t87 = 0.1e1 / t146;
	t84 = 0.1e1 / t86;
	t83 = 0.1e1 / (t99 ^ 2 * t85 + 0.1e1);
	t82 = (-t103 * t107 + t116 * t159) * t94 * t125;
	t81 = (-t101 * t103 + t113 * t159) * t94;
	t80 = (-t103 * t98 + t106 * t159) * t94;
	t1 = [0, t81, t82, t80, 0, 0; 0, ((t112 * t125 - t122 * t156) * t84 - (-t90 * t101 + t91 * t113 + t144 * t81) * t161) * t83, (-t109 * t125 * t84 - (t144 * t82 + (-t107 * t90 + t116 * t91) * t125) * t161) * t83, (t100 * t84 - (t91 * t106 + t144 * t80 - t90 * t98) * t161) * t83, 0, 0; 0, ((t102 * t134 - t111 * t137) * t88 - (t102 * t137 + t111 * t134) * t160) * t87, ((-t110 * t137 - t126 * t158) * t88 - (t110 * t134 - t126 * t157) * t160) * t87, (-t134 * t88 + t137 * t160) * t99 * t87, 0, t146 * t87;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end