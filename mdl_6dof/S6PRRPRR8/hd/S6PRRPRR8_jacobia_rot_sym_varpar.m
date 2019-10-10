% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR8
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
%   Wie in S6PRRPRR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
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
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:10
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (473->35), mult. (1448->85), div. (57->9), fcn. (1990->13), ass. (0->52)
	t76 = cos(pkin(12));
	t77 = cos(pkin(7));
	t73 = sin(pkin(12));
	t80 = sin(qJ(2));
	t78 = cos(pkin(6));
	t82 = cos(qJ(2));
	t90 = t78 * t82;
	t84 = -t73 * t80 + t76 * t90;
	t74 = sin(pkin(7));
	t75 = sin(pkin(6));
	t93 = t75 * t74;
	t98 = -t76 * t93 + t84 * t77;
	t91 = t78 * t80;
	t69 = t73 * t82 + t76 * t91;
	t79 = sin(qJ(3));
	t81 = cos(qJ(3));
	t54 = t69 * t79 - t98 * t81;
	t92 = t77 * t81;
	t94 = t74 * t78;
	t64 = -t81 * t94 + (t79 * t80 - t82 * t92) * t75;
	t53 = atan2(-t54, t64);
	t49 = sin(t53);
	t50 = cos(t53);
	t48 = -t49 * t54 + t50 * t64;
	t47 = 0.1e1 / t48 ^ 2;
	t70 = -t73 * t90 - t76 * t80;
	t85 = t70 * t77 + t73 * t93;
	t71 = -t73 * t91 + t76 * t82;
	t95 = t71 * t79;
	t57 = -t85 * t81 + t95;
	t97 = t47 * t57;
	t61 = 0.1e1 / t64 ^ 2;
	t96 = t54 * t61;
	t89 = t79 * t82;
	t88 = t80 * t81;
	t86 = -t49 * t64 - t50 * t54;
	t68 = (t77 * t88 + t89) * t75;
	t66 = t73 * t75 * t77 - t70 * t74;
	t65 = t79 * t94 + (t77 * t89 + t88) * t75;
	t63 = 0.1e1 / t66 ^ 2;
	t62 = 0.1e1 / t66;
	t60 = 0.1e1 / t64;
	t59 = t69 * t92 + t84 * t79;
	t58 = t71 * t81 + t85 * t79;
	t56 = t69 * t81 + t98 * t79;
	t52 = 0.1e1 / (t58 ^ 2 * t63 + 0.1e1);
	t51 = 0.1e1 / (t54 ^ 2 * t61 + 0.1e1);
	t46 = 0.1e1 / t48;
	t45 = 0.1e1 / (t57 ^ 2 * t47 + 0.1e1);
	t44 = (-t59 * t60 + t68 * t96) * t51;
	t43 = (-t56 * t60 + t65 * t96) * t51;
	t1 = [0, t44, t43, 0, 0, 0; 0, ((t70 * t79 + t71 * t92) * t46 - (t86 * t44 - t49 * t59 + t50 * t68) * t97) * t45, (t58 * t46 - (t86 * t43 - t49 * t56 + t50 * t65) * t97) * t45, 0, 0, 0; 0, ((t70 * t81 - t77 * t95) * t62 - t71 * t74 * t58 * t63) * t52, -t57 * t62 * t52, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:10
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (607->40), mult. (1825->96), div. (65->9), fcn. (2498->15), ass. (0->61)
	t90 = sin(pkin(7));
	t91 = sin(pkin(6));
	t112 = t90 * t91;
	t92 = cos(pkin(12));
	t104 = t92 * t112;
	t93 = cos(pkin(7));
	t96 = sin(qJ(3));
	t110 = t93 * t96;
	t100 = cos(qJ(2));
	t105 = t92 * t100;
	t89 = sin(pkin(12));
	t94 = cos(pkin(6));
	t97 = sin(qJ(2));
	t83 = t94 * t105 - t89 * t97;
	t106 = t89 * t100;
	t109 = t94 * t97;
	t84 = t92 * t109 + t106;
	t99 = cos(qJ(3));
	t71 = -t96 * t104 + t83 * t110 + t84 * t99;
	t111 = t90 * t94;
	t80 = t96 * t111 + (t100 * t110 + t97 * t99) * t91;
	t69 = atan2(-t71, t80);
	t66 = sin(t69);
	t67 = cos(t69);
	t60 = -t66 * t71 + t67 * t80;
	t59 = 0.1e1 / t60 ^ 2;
	t85 = -t94 * t106 - t92 * t97;
	t101 = t89 * t112 + t85 * t93;
	t86 = -t89 * t109 + t105;
	t113 = t86 * t99;
	t74 = t101 * t96 + t113;
	t117 = t59 * t74;
	t73 = -t101 * t99 + t86 * t96;
	t81 = t89 * t91 * t93 - t85 * t90;
	t95 = sin(qJ(5));
	t98 = cos(qJ(5));
	t65 = t73 * t95 + t81 * t98;
	t63 = 0.1e1 / t65 ^ 2;
	t64 = -t73 * t98 + t81 * t95;
	t116 = t63 * t64;
	t78 = 0.1e1 / t80 ^ 2;
	t115 = t71 * t78;
	t114 = t86 * t90;
	t108 = t96 * t97;
	t107 = t100 * t99;
	t103 = t64 ^ 2 * t63 + 0.1e1;
	t102 = -t66 * t80 - t67 * t71;
	t82 = (-t93 * t108 + t107) * t91;
	t79 = t99 * t111 + (t93 * t107 - t108) * t91;
	t77 = 0.1e1 / t80;
	t76 = t93 * t113 + t85 * t96;
	t75 = -t84 * t110 + t83 * t99;
	t70 = t84 * t96 + (-t83 * t93 + t104) * t99;
	t68 = 0.1e1 / (t71 ^ 2 * t78 + 0.1e1);
	t62 = 0.1e1 / t65;
	t61 = 0.1e1 / t103;
	t58 = 0.1e1 / t60;
	t57 = 0.1e1 / (t74 ^ 2 * t59 + 0.1e1);
	t56 = (t82 * t115 - t75 * t77) * t68;
	t55 = (t79 * t115 + t70 * t77) * t68;
	t1 = [0, t56, t55, 0, 0, 0; 0, ((-t86 * t110 + t85 * t99) * t58 - (t102 * t56 - t66 * t75 + t67 * t82) * t117) * t57, (-t73 * t58 - (t102 * t55 + t66 * t70 + t67 * t79) * t117) * t57, 0, 0, 0; 0, ((t95 * t114 - t76 * t98) * t62 - (t98 * t114 + t76 * t95) * t116) * t61, (-t95 * t116 - t98 * t62) * t74 * t61, 0, t103 * t61, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:10
	% EndTime: 2019-10-09 22:39:10
	% DurationCPUTime: 0.57s
	% Computational Cost: add. (1524->59), mult. (4378->141), div. (95->9), fcn. (6002->17), ass. (0->76)
	t121 = sin(pkin(6));
	t127 = sin(qJ(3));
	t128 = sin(qJ(2));
	t131 = cos(qJ(3));
	t132 = cos(qJ(2));
	t123 = cos(pkin(7));
	t144 = t123 * t131;
	t120 = sin(pkin(7));
	t124 = cos(pkin(6));
	t148 = t120 * t124;
	t109 = -t131 * t148 + (t128 * t127 - t132 * t144) * t121;
	t146 = t121 * t120;
	t114 = t124 * t123 - t132 * t146;
	t126 = sin(qJ(5));
	t130 = cos(qJ(5));
	t104 = -t109 * t130 + t114 * t126;
	t122 = cos(pkin(12));
	t119 = sin(pkin(12));
	t142 = t124 * t132;
	t136 = -t119 * t128 + t122 * t142;
	t145 = t121 * t123;
	t111 = -t136 * t120 - t122 * t145;
	t134 = -t122 * t146 + t136 * t123;
	t143 = t124 * t128;
	t135 = t119 * t132 + t122 * t143;
	t133 = t135 * t127 - t134 * t131;
	t90 = t111 * t126 - t133 * t130;
	t89 = atan2(-t90, t104);
	t86 = sin(t89);
	t87 = cos(t89);
	t80 = t87 * t104 - t86 * t90;
	t79 = 0.1e1 / t80 ^ 2;
	t115 = -t119 * t142 - t122 * t128;
	t138 = t119 * t146;
	t116 = -t119 * t143 + t122 * t132;
	t149 = t116 * t127;
	t102 = -t115 * t144 - t131 * t138 + t149;
	t112 = -t115 * t120 + t119 * t145;
	t93 = -t102 * t130 + t112 * t126;
	t154 = t79 * t93;
	t129 = cos(qJ(6));
	t103 = t116 * t131 + (t115 * t123 + t138) * t127;
	t125 = sin(qJ(6));
	t151 = t103 * t125;
	t94 = t102 * t126 + t112 * t130;
	t85 = t94 * t129 + t151;
	t83 = 0.1e1 / t85 ^ 2;
	t150 = t103 * t129;
	t84 = t94 * t125 - t150;
	t153 = t83 * t84;
	t100 = 0.1e1 / t104 ^ 2;
	t152 = t100 * t90;
	t147 = t120 * t126;
	t141 = t128 * t131;
	t140 = t132 * t127;
	t139 = t84 ^ 2 * t83 + 0.1e1;
	t137 = -t104 * t86 - t87 * t90;
	t110 = t127 * t148 + (t123 * t140 + t141) * t121;
	t108 = (t128 * t147 - (t123 * t141 + t140) * t130) * t121;
	t107 = t115 * t131 - t123 * t149;
	t106 = t115 * t127 + t116 * t144;
	t105 = t109 * t126 + t114 * t130;
	t101 = t134 * t127 + t135 * t131;
	t99 = 0.1e1 / t104;
	t96 = t116 * t120 * t130 + t106 * t126;
	t95 = -t135 * t147 + (t136 * t127 + t135 * t144) * t130;
	t92 = t111 * t130 + t133 * t126;
	t88 = 0.1e1 / (t90 ^ 2 * t100 + 0.1e1);
	t82 = 0.1e1 / t85;
	t81 = 0.1e1 / t139;
	t78 = 0.1e1 / t80;
	t77 = 0.1e1 / (t93 ^ 2 * t79 + 0.1e1);
	t76 = (t101 * t99 - t110 * t152) * t88 * t130;
	t75 = (t108 * t152 + t95 * t99) * t88;
	t74 = (t105 * t152 - t92 * t99) * t88;
	t1 = [0, t75, t76, 0, t74, 0; 0, ((-t106 * t130 + t116 * t147) * t78 - (t87 * t108 + t137 * t75 + t86 * t95) * t154) * t77, (-t103 * t130 * t78 - (t137 * t76 + (t101 * t86 - t110 * t87) * t130) * t154) * t77, 0, (t94 * t78 - (t87 * t105 + t137 * t74 - t86 * t92) * t154) * t77, 0; 0, ((-t107 * t129 + t96 * t125) * t82 - (t107 * t125 + t96 * t129) * t153) * t81, ((t102 * t129 + t126 * t151) * t82 - (-t102 * t125 + t126 * t150) * t153) * t81, 0, (-t125 * t82 + t129 * t153) * t93 * t81, t139 * t81;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end