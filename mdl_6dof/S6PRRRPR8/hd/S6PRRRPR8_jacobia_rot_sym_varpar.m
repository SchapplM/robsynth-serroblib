% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR8
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
%   Wie in S6PRRRPR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:00
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:05
	% DurationCPUTime: 0.17s
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
	t32 = t55 * t51 + t62;
	t30 = 0.1e1 / t32 ^ 2;
	t31 = -t55 * t53 + t63;
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
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:05
	% DurationCPUTime: 0.35s
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
	% StartTime: 2019-10-09 23:00:05
	% EndTime: 2019-10-09 23:00:05
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (1245->53), mult. (3597->128), div. (87->9), fcn. (4945->15), ass. (0->67)
	t114 = sin(qJ(4));
	t117 = cos(qJ(4));
	t109 = sin(pkin(7));
	t111 = cos(pkin(12));
	t108 = sin(pkin(12));
	t116 = sin(qJ(2));
	t113 = cos(pkin(6));
	t119 = cos(qJ(2));
	t127 = t113 * t119;
	t122 = -t108 * t116 + t111 * t127;
	t110 = sin(pkin(6));
	t112 = cos(pkin(7));
	t131 = t110 * t112;
	t121 = -t122 * t109 - t111 * t131;
	t128 = t113 * t116;
	t103 = t108 * t119 + t111 * t128;
	t115 = sin(qJ(3));
	t118 = cos(qJ(3));
	t132 = t110 * t109;
	t120 = -t111 * t132 + t122 * t112;
	t90 = t103 * t118 + t120 * t115;
	t79 = t90 * t114 - t121 * t117;
	t102 = t113 * t112 - t119 * t132;
	t130 = t112 * t115;
	t134 = t109 * t113;
	t99 = t115 * t134 + (t116 * t118 + t119 * t130) * t110;
	t93 = -t102 * t117 + t99 * t114;
	t78 = atan2(-t79, t93);
	t74 = sin(t78);
	t75 = cos(t78);
	t73 = -t74 * t79 + t75 * t93;
	t72 = 0.1e1 / t73 ^ 2;
	t104 = -t108 * t127 - t111 * t116;
	t100 = -t104 * t109 + t108 * t131;
	t105 = -t108 * t128 + t111 * t119;
	t124 = t108 * t132;
	t92 = t105 * t118 + (t104 * t112 + t124) * t115;
	t82 = -t100 * t117 + t92 * t114;
	t137 = t72 * t82;
	t88 = 0.1e1 / t93 ^ 2;
	t136 = t79 * t88;
	t83 = t100 * t114 + t92 * t117;
	t129 = t112 * t118;
	t91 = -t104 * t129 + t105 * t115 - t118 * t124;
	t86 = 0.1e1 / t91 ^ 2;
	t135 = t83 * t86;
	t133 = t109 * t117;
	t126 = t115 * t116;
	t125 = t118 * t119;
	t123 = -t74 * t93 - t75 * t79;
	t98 = t118 * t134 + (t112 * t125 - t126) * t110;
	t96 = ((-t112 * t126 + t125) * t114 - t116 * t133) * t110;
	t95 = t104 * t118 - t105 * t130;
	t94 = t102 * t114 + t99 * t117;
	t89 = -t103 * t115 + t120 * t118;
	t87 = 0.1e1 / t93;
	t85 = 0.1e1 / t91;
	t84 = (-t103 * t130 + t122 * t118) * t114 - t103 * t133;
	t81 = t121 * t114 + t90 * t117;
	t77 = 0.1e1 / (t79 ^ 2 * t88 + 0.1e1);
	t76 = 0.1e1 / (t83 ^ 2 * t86 + 0.1e1);
	t71 = 0.1e1 / t73;
	t70 = 0.1e1 / (t82 ^ 2 * t72 + 0.1e1);
	t69 = (t98 * t136 - t87 * t89) * t77 * t114;
	t68 = (t96 * t136 - t84 * t87) * t77;
	t67 = (t94 * t136 - t81 * t87) * t77;
	t1 = [0, t68, t69, t67, 0, 0; 0, ((-t105 * t133 + t95 * t114) * t71 - (t123 * t68 - t74 * t84 + t75 * t96) * t137) * t70, (-t91 * t114 * t71 - (t123 * t69 + (-t74 * t89 + t75 * t98) * t114) * t137) * t70, (t83 * t71 - (t123 * t67 - t74 * t81 + t75 * t94) * t137) * t70, 0, 0; 0, ((t105 * t109 * t114 + t95 * t117) * t85 - (t104 * t115 + t105 * t129) * t135) * t76, (-t117 * t85 * t91 - t92 * t135) * t76, -t82 * t85 * t76, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:05
	% EndTime: 2019-10-09 23:00:05
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (1524->59), mult. (4378->142), div. (95->9), fcn. (6002->17), ass. (0->76)
	t126 = sin(pkin(6));
	t132 = sin(qJ(3));
	t133 = sin(qJ(2));
	t136 = cos(qJ(3));
	t137 = cos(qJ(2));
	t128 = cos(pkin(7));
	t147 = t128 * t132;
	t125 = sin(pkin(7));
	t129 = cos(pkin(6));
	t150 = t125 * t129;
	t113 = t132 * t150 + (t133 * t136 + t137 * t147) * t126;
	t151 = t125 * t126;
	t117 = t129 * t128 - t137 * t151;
	t131 = sin(qJ(4));
	t135 = cos(qJ(4));
	t108 = t113 * t135 + t117 * t131;
	t124 = sin(pkin(12));
	t127 = cos(pkin(12));
	t145 = t129 * t133;
	t119 = t124 * t137 + t127 * t145;
	t144 = t129 * t137;
	t118 = -t124 * t133 + t127 * t144;
	t138 = t118 * t128 - t127 * t151;
	t104 = t119 * t136 + t138 * t132;
	t148 = t126 * t128;
	t114 = -t118 * t125 - t127 * t148;
	t95 = t104 * t135 + t114 * t131;
	t93 = atan2(-t95, t108);
	t90 = sin(t93);
	t91 = cos(t93);
	t84 = t91 * t108 - t90 * t95;
	t83 = 0.1e1 / t84 ^ 2;
	t120 = -t124 * t144 - t127 * t133;
	t121 = -t124 * t145 + t127 * t137;
	t140 = t124 * t151;
	t106 = t121 * t136 + (t120 * t128 + t140) * t132;
	t115 = -t120 * t125 + t124 * t148;
	t98 = t106 * t135 + t115 * t131;
	t156 = t83 * t98;
	t130 = sin(qJ(6));
	t146 = t128 * t136;
	t105 = -t120 * t146 + t121 * t132 - t136 * t140;
	t134 = cos(qJ(6));
	t152 = t105 * t134;
	t97 = t106 * t131 - t115 * t135;
	t89 = t97 * t130 + t152;
	t87 = 0.1e1 / t89 ^ 2;
	t153 = t105 * t130;
	t88 = -t97 * t134 + t153;
	t155 = t87 * t88;
	t102 = 0.1e1 / t108 ^ 2;
	t154 = t102 * t95;
	t149 = t125 * t131;
	t143 = t132 * t133;
	t142 = t136 * t137;
	t141 = t88 ^ 2 * t87 + 0.1e1;
	t139 = -t108 * t90 - t91 * t95;
	t112 = t136 * t150 + (t128 * t142 - t143) * t126;
	t111 = ((-t128 * t143 + t142) * t135 + t133 * t149) * t126;
	t110 = t120 * t136 - t121 * t147;
	t109 = t120 * t132 + t121 * t146;
	t107 = -t113 * t131 + t117 * t135;
	t103 = -t119 * t132 + t138 * t136;
	t101 = 0.1e1 / t108;
	t100 = -t121 * t125 * t135 + t110 * t131;
	t99 = (t118 * t136 - t119 * t147) * t135 + t119 * t149;
	t94 = t104 * t131 - t114 * t135;
	t92 = 0.1e1 / (t95 ^ 2 * t102 + 0.1e1);
	t86 = 0.1e1 / t89;
	t85 = 0.1e1 / t141;
	t82 = 0.1e1 / t84;
	t81 = 0.1e1 / (t98 ^ 2 * t83 + 0.1e1);
	t80 = (-t101 * t103 + t112 * t154) * t92 * t135;
	t79 = (-t101 * t99 + t111 * t154) * t92;
	t78 = (t101 * t94 + t107 * t154) * t92;
	t1 = [0, t79, t80, t78, 0, 0; 0, ((t110 * t135 + t121 * t149) * t82 - (t91 * t111 + t139 * t79 - t90 * t99) * t156) * t81, (-t105 * t135 * t82 - (t139 * t80 + (-t103 * t90 + t112 * t91) * t135) * t156) * t81, (-t97 * t82 - (t91 * t107 + t139 * t78 + t90 * t94) * t156) * t81, 0, 0; 0, ((-t100 * t134 + t109 * t130) * t86 - (t100 * t130 + t109 * t134) * t155) * t85, ((t106 * t130 + t131 * t152) * t86 - (t106 * t134 - t131 * t153) * t155) * t85, (-t130 * t155 - t134 * t86) * t98 * t85, 0, t141 * t85;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end