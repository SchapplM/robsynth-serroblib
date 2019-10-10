% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR6
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
%   Wie in S6PRRPRR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
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
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:20
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (555->40), mult. (1679->96), div. (60->9), fcn. (2302->15), ass. (0->59)
	t94 = cos(pkin(6));
	t98 = cos(qJ(2));
	t106 = t94 * t98;
	t88 = sin(pkin(12));
	t92 = cos(pkin(12));
	t96 = sin(qJ(2));
	t100 = t92 * t106 - t88 * t96;
	t89 = sin(pkin(7));
	t90 = sin(pkin(6));
	t109 = t90 * t89;
	t93 = cos(pkin(7));
	t116 = t100 * t93 - t92 * t109;
	t107 = t94 * t96;
	t82 = t92 * t107 + t88 * t98;
	t95 = sin(qJ(3));
	t97 = cos(qJ(3));
	t67 = -t116 * t97 + t82 * t95;
	t108 = t93 * t97;
	t110 = t89 * t94;
	t76 = -t97 * t110 + (-t98 * t108 + t95 * t96) * t90;
	t66 = atan2(-t67, t76);
	t63 = sin(t66);
	t64 = cos(t66);
	t57 = -t63 * t67 + t64 * t76;
	t56 = 0.1e1 / t57 ^ 2;
	t103 = t88 * t109;
	t84 = -t88 * t107 + t92 * t98;
	t111 = t84 * t95;
	t83 = -t88 * t106 - t92 * t96;
	t70 = -t97 * t103 - t83 * t108 + t111;
	t115 = t56 * t70;
	t71 = t84 * t97 + (t83 * t93 + t103) * t95;
	t78 = t88 * t90 * t93 - t83 * t89;
	t87 = sin(pkin(13));
	t91 = cos(pkin(13));
	t62 = t71 * t91 + t78 * t87;
	t60 = 0.1e1 / t62 ^ 2;
	t61 = t71 * t87 - t78 * t91;
	t114 = t60 * t61;
	t75 = 0.1e1 / t76 ^ 2;
	t113 = t67 * t75;
	t112 = t84 * t89;
	t105 = t95 * t98;
	t104 = t96 * t97;
	t101 = -t63 * t76 - t64 * t67;
	t81 = (t93 * t104 + t105) * t90;
	t77 = t95 * t110 + (t93 * t105 + t104) * t90;
	t74 = 0.1e1 / t76;
	t73 = -t93 * t111 + t83 * t97;
	t72 = t100 * t95 + t82 * t108;
	t69 = t116 * t95 + t82 * t97;
	t65 = 0.1e1 / (t67 ^ 2 * t75 + 0.1e1);
	t59 = 0.1e1 / t62;
	t58 = 0.1e1 / (t61 ^ 2 * t60 + 0.1e1);
	t55 = 0.1e1 / t57;
	t54 = 0.1e1 / (t70 ^ 2 * t56 + 0.1e1);
	t53 = (t81 * t113 - t72 * t74) * t65;
	t52 = (t77 * t113 - t69 * t74) * t65;
	t1 = [0, t53, t52, 0, 0, 0; 0, ((t84 * t108 + t83 * t95) * t55 - (t101 * t53 - t63 * t72 + t64 * t81) * t115) * t54, (t71 * t55 - (t101 * t52 - t63 * t69 + t64 * t77) * t115) * t54, 0, 0, 0; 0, ((-t91 * t112 + t73 * t87) * t59 - (t87 * t112 + t73 * t91) * t114) * t58, (t91 * t114 - t87 * t59) * t70 * t58, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:20
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (655->41), mult. (1825->98), div. (65->9), fcn. (2498->15), ass. (0->61)
	t100 = cos(pkin(12));
	t103 = sin(qJ(3));
	t105 = cos(qJ(3));
	t101 = cos(pkin(7));
	t104 = sin(qJ(2));
	t102 = cos(pkin(6));
	t106 = cos(qJ(2));
	t114 = t102 * t106;
	t97 = sin(pkin(12));
	t108 = t100 * t114 - t97 * t104;
	t107 = t108 * t101;
	t98 = sin(pkin(7));
	t99 = sin(pkin(6));
	t119 = t98 * t99;
	t112 = t105 * t119;
	t115 = t102 * t104;
	t89 = t100 * t115 + t97 * t106;
	t74 = t100 * t112 + t89 * t103 - t105 * t107;
	t118 = t101 * t99;
	t109 = t102 * t98 + t106 * t118;
	t83 = t99 * t104 * t103 - t109 * t105;
	t73 = atan2(-t74, t83);
	t70 = sin(t73);
	t71 = cos(t73);
	t64 = -t70 * t74 + t71 * t83;
	t63 = 0.1e1 / t64 ^ 2;
	t116 = t101 * t105;
	t91 = t100 * t106 - t97 * t115;
	t117 = t91 * t103;
	t90 = -t100 * t104 - t97 * t114;
	t77 = -t97 * t112 - t90 * t116 + t117;
	t123 = t63 * t77;
	t78 = t91 * t105 + (t101 * t90 + t97 * t119) * t103;
	t85 = t97 * t118 - t90 * t98;
	t96 = pkin(13) + qJ(5);
	t94 = sin(t96);
	t95 = cos(t96);
	t69 = t78 * t95 + t85 * t94;
	t67 = 0.1e1 / t69 ^ 2;
	t68 = t78 * t94 - t85 * t95;
	t122 = t67 * t68;
	t82 = 0.1e1 / t83 ^ 2;
	t121 = t74 * t82;
	t120 = t91 * t98;
	t113 = t104 * t105;
	t111 = t68 ^ 2 * t67 + 0.1e1;
	t110 = -t70 * t83 - t71 * t74;
	t88 = (t101 * t113 + t103 * t106) * t99;
	t84 = t109 * t103 + t99 * t113;
	t81 = 0.1e1 / t83;
	t80 = -t101 * t117 + t90 * t105;
	t79 = t108 * t103 + t89 * t116;
	t76 = t89 * t105 + (-t100 * t119 + t107) * t103;
	t72 = 0.1e1 / (t74 ^ 2 * t82 + 0.1e1);
	t66 = 0.1e1 / t69;
	t65 = 0.1e1 / t111;
	t62 = 0.1e1 / t64;
	t61 = 0.1e1 / (t77 ^ 2 * t63 + 0.1e1);
	t60 = (t88 * t121 - t79 * t81) * t72;
	t59 = (t84 * t121 - t76 * t81) * t72;
	t1 = [0, t60, t59, 0, 0, 0; 0, ((t90 * t103 + t91 * t116) * t62 - (t110 * t60 - t70 * t79 + t71 * t88) * t123) * t61, (t78 * t62 - (t110 * t59 - t70 * t76 + t71 * t84) * t123) * t61, 0, 0, 0; 0, ((-t95 * t120 + t80 * t94) * t66 - (t94 * t120 + t80 * t95) * t122) * t65, (t95 * t122 - t94 * t66) * t77 * t65, 0, t111 * t65, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:20
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (1983->60), mult. (4378->142), div. (95->9), fcn. (6002->17), ass. (0->77)
	t129 = sin(pkin(6));
	t134 = sin(qJ(3));
	t135 = sin(qJ(2));
	t137 = cos(qJ(3));
	t138 = cos(qJ(2));
	t131 = cos(pkin(7));
	t151 = t131 * t134;
	t128 = sin(pkin(7));
	t132 = cos(pkin(6));
	t154 = t128 * t132;
	t116 = t134 * t154 + (t135 * t137 + t138 * t151) * t129;
	t153 = t129 * t128;
	t118 = t132 * t131 - t138 * t153;
	t126 = pkin(13) + qJ(5);
	t124 = sin(t126);
	t125 = cos(t126);
	t104 = t116 * t124 - t118 * t125;
	t127 = sin(pkin(12));
	t130 = cos(pkin(12));
	t149 = t132 * t135;
	t119 = t127 * t138 + t130 * t149;
	t148 = t132 * t138;
	t141 = -t127 * t135 + t130 * t148;
	t139 = -t130 * t153 + t141 * t131;
	t107 = t119 * t137 + t139 * t134;
	t152 = t129 * t131;
	t140 = -t141 * t128 - t130 * t152;
	t95 = t107 * t124 - t140 * t125;
	t94 = atan2(-t95, t104);
	t89 = sin(t94);
	t90 = cos(t94);
	t85 = t90 * t104 - t89 * t95;
	t84 = 0.1e1 / t85 ^ 2;
	t120 = -t127 * t148 - t130 * t135;
	t121 = -t127 * t149 + t130 * t138;
	t144 = t127 * t153;
	t109 = t121 * t137 + (t120 * t131 + t144) * t134;
	t142 = -t120 * t128 + t127 * t152;
	t98 = t109 * t124 - t142 * t125;
	t160 = t84 * t98;
	t136 = cos(qJ(6));
	t150 = t131 * t137;
	t108 = -t120 * t150 + t121 * t134 - t137 * t144;
	t133 = sin(qJ(6));
	t157 = t108 * t133;
	t99 = t109 * t125 + t142 * t124;
	t92 = t99 * t136 + t157;
	t88 = 0.1e1 / t92 ^ 2;
	t156 = t108 * t136;
	t91 = t99 * t133 - t156;
	t159 = t88 * t91;
	t103 = 0.1e1 / t104 ^ 2;
	t158 = t103 * t95;
	t155 = t128 * t125;
	t147 = t134 * t135;
	t146 = t137 * t138;
	t145 = t91 ^ 2 * t88 + 0.1e1;
	t143 = -t104 * t89 - t90 * t95;
	t115 = t137 * t154 + (t131 * t146 - t147) * t129;
	t112 = ((-t131 * t147 + t146) * t124 - t135 * t155) * t129;
	t111 = t120 * t137 - t121 * t151;
	t110 = t120 * t134 + t121 * t150;
	t106 = -t119 * t134 + t139 * t137;
	t105 = t116 * t125 + t118 * t124;
	t102 = 0.1e1 / t104;
	t101 = t121 * t128 * t124 + t111 * t125;
	t100 = (-t119 * t151 + t141 * t137) * t124 - t119 * t155;
	t97 = t107 * t125 + t140 * t124;
	t93 = 0.1e1 / (t95 ^ 2 * t103 + 0.1e1);
	t87 = 0.1e1 / t92;
	t86 = 0.1e1 / t145;
	t83 = 0.1e1 / t85;
	t82 = 0.1e1 / (t98 ^ 2 * t84 + 0.1e1);
	t81 = (-t102 * t106 + t115 * t158) * t93 * t124;
	t80 = (-t100 * t102 + t112 * t158) * t93;
	t79 = (-t102 * t97 + t105 * t158) * t93;
	t1 = [0, t80, t81, 0, t79, 0; 0, ((t111 * t124 - t121 * t155) * t83 - (-t89 * t100 + t90 * t112 + t143 * t80) * t160) * t82, (-t108 * t124 * t83 - (t143 * t81 + (-t106 * t89 + t115 * t90) * t124) * t160) * t82, 0, (t99 * t83 - (t90 * t105 + t143 * t79 - t89 * t97) * t160) * t82, 0; 0, ((t101 * t133 - t110 * t136) * t87 - (t101 * t136 + t110 * t133) * t159) * t86, ((-t109 * t136 - t125 * t157) * t87 - (t109 * t133 - t125 * t156) * t159) * t86, 0, (-t133 * t87 + t136 * t159) * t98 * t86, t145 * t86;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end