% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR7
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
%   Wie in S6PRRRPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:58
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobia_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
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
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.40s
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
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (1437->59), mult. (4134->141), div. (90->9), fcn. (5673->17), ass. (0->75)
	t120 = sin(pkin(6));
	t126 = sin(qJ(3));
	t127 = sin(qJ(2));
	t129 = cos(qJ(3));
	t130 = cos(qJ(2));
	t123 = cos(pkin(7));
	t142 = t123 * t126;
	t119 = sin(pkin(7));
	t124 = cos(pkin(6));
	t146 = t119 * t124;
	t109 = t126 * t146 + (t127 * t129 + t130 * t142) * t120;
	t144 = t120 * t119;
	t111 = t124 * t123 - t130 * t144;
	t125 = sin(qJ(4));
	t128 = cos(qJ(4));
	t101 = t109 * t125 - t111 * t128;
	t122 = cos(pkin(12));
	t118 = sin(pkin(12));
	t139 = t124 * t130;
	t133 = -t118 * t127 + t122 * t139;
	t143 = t120 * t123;
	t132 = -t133 * t119 - t122 * t143;
	t140 = t124 * t127;
	t112 = t118 * t130 + t122 * t140;
	t131 = -t122 * t144 + t133 * t123;
	t98 = t112 * t129 + t131 * t126;
	t88 = t98 * t125 - t132 * t128;
	t87 = atan2(-t88, t101);
	t84 = sin(t87);
	t85 = cos(t87);
	t78 = t85 * t101 - t84 * t88;
	t77 = 0.1e1 / t78 ^ 2;
	t113 = -t118 * t139 - t122 * t127;
	t114 = -t118 * t140 + t122 * t130;
	t136 = t118 * t144;
	t100 = t114 * t129 + (t113 * t123 + t136) * t126;
	t134 = -t113 * t119 + t118 * t143;
	t91 = t100 * t125 - t134 * t128;
	t151 = t77 * t91;
	t121 = cos(pkin(13));
	t117 = sin(pkin(13));
	t141 = t123 * t129;
	t99 = -t113 * t141 + t114 * t126 - t129 * t136;
	t148 = t99 * t117;
	t92 = t100 * t128 + t134 * t125;
	t83 = t92 * t121 + t148;
	t81 = 0.1e1 / t83 ^ 2;
	t147 = t99 * t121;
	t82 = t92 * t117 - t147;
	t150 = t81 * t82;
	t96 = 0.1e1 / t101 ^ 2;
	t149 = t88 * t96;
	t145 = t119 * t128;
	t138 = t126 * t127;
	t137 = t129 * t130;
	t135 = -t101 * t84 - t85 * t88;
	t108 = t129 * t146 + (t123 * t137 - t138) * t120;
	t105 = ((-t123 * t138 + t137) * t125 - t127 * t145) * t120;
	t104 = t113 * t129 - t114 * t142;
	t103 = t113 * t126 + t114 * t141;
	t102 = t109 * t128 + t111 * t125;
	t97 = -t112 * t126 + t131 * t129;
	t95 = 0.1e1 / t101;
	t94 = t114 * t119 * t125 + t104 * t128;
	t93 = (-t112 * t142 + t133 * t129) * t125 - t112 * t145;
	t90 = t132 * t125 + t98 * t128;
	t86 = 0.1e1 / (t88 ^ 2 * t96 + 0.1e1);
	t80 = 0.1e1 / t83;
	t79 = 0.1e1 / (t82 ^ 2 * t81 + 0.1e1);
	t76 = 0.1e1 / t78;
	t75 = 0.1e1 / (t91 ^ 2 * t77 + 0.1e1);
	t74 = (t108 * t149 - t95 * t97) * t86 * t125;
	t73 = (t105 * t149 - t93 * t95) * t86;
	t72 = (t102 * t149 - t90 * t95) * t86;
	t1 = [0, t73, t74, t72, 0, 0; 0, ((t104 * t125 - t114 * t145) * t76 - (t85 * t105 + t135 * t73 - t84 * t93) * t151) * t75, (-t99 * t125 * t76 - (t135 * t74 + (t108 * t85 - t84 * t97) * t125) * t151) * t75, (t92 * t76 - (t85 * t102 + t135 * t72 - t84 * t90) * t151) * t75, 0, 0; 0, ((-t103 * t121 + t94 * t117) * t80 - (t103 * t117 + t94 * t121) * t150) * t79, ((-t100 * t121 - t128 * t148) * t80 - (t100 * t117 - t128 * t147) * t150) * t79, (-t117 * t80 + t121 * t150) * t91 * t79, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:05
	% DurationCPUTime: 0.53s
	% Computational Cost: add. (1590->60), mult. (4378->143), div. (95->9), fcn. (6002->17), ass. (0->76)
	t130 = sin(pkin(12));
	t133 = cos(pkin(12));
	t138 = sin(qJ(2));
	t135 = cos(pkin(6));
	t141 = cos(qJ(2));
	t151 = t135 * t141;
	t123 = -t130 * t151 - t133 * t138;
	t152 = t135 * t138;
	t124 = -t130 * t152 + t133 * t141;
	t134 = cos(pkin(7));
	t137 = sin(qJ(3));
	t140 = cos(qJ(3));
	t131 = sin(pkin(7));
	t132 = sin(pkin(6));
	t156 = t132 * t131;
	t147 = t130 * t156;
	t110 = t124 * t140 + (t123 * t134 + t147) * t137;
	t136 = sin(qJ(4));
	t139 = cos(qJ(4));
	t155 = t132 * t134;
	t145 = -t123 * t131 + t130 * t155;
	t102 = t110 * t139 + t145 * t136;
	t153 = t134 * t140;
	t109 = -t123 * t153 + t124 * t137 - t140 * t147;
	t129 = pkin(13) + qJ(6);
	t127 = sin(t129);
	t128 = cos(t129);
	t93 = t102 * t128 + t109 * t127;
	t91 = 0.1e1 / t93 ^ 2;
	t92 = t102 * t127 - t109 * t128;
	t162 = t91 * t92;
	t101 = t110 * t136 - t145 * t139;
	t154 = t134 * t137;
	t158 = t131 * t135;
	t119 = t137 * t158 + (t138 * t140 + t141 * t154) * t132;
	t121 = t135 * t134 - t141 * t156;
	t111 = t119 * t136 - t121 * t139;
	t122 = t130 * t141 + t133 * t152;
	t144 = -t130 * t138 + t133 * t151;
	t142 = -t133 * t156 + t144 * t134;
	t108 = t122 * t140 + t142 * t137;
	t143 = -t144 * t131 - t133 * t155;
	t98 = t108 * t136 - t143 * t139;
	t97 = atan2(-t98, t111);
	t94 = sin(t97);
	t95 = cos(t97);
	t88 = t111 * t95 - t94 * t98;
	t87 = 0.1e1 / t88 ^ 2;
	t161 = t101 * t87;
	t106 = 0.1e1 / t111 ^ 2;
	t160 = t106 * t98;
	t159 = t109 * t139;
	t157 = t131 * t139;
	t150 = t137 * t138;
	t149 = t140 * t141;
	t148 = t91 * t92 ^ 2 + 0.1e1;
	t146 = -t111 * t94 - t95 * t98;
	t118 = t140 * t158 + (t134 * t149 - t150) * t132;
	t115 = ((-t134 * t150 + t149) * t136 - t138 * t157) * t132;
	t114 = t123 * t140 - t124 * t154;
	t113 = t123 * t137 + t124 * t153;
	t112 = t119 * t139 + t121 * t136;
	t107 = -t122 * t137 + t142 * t140;
	t105 = 0.1e1 / t111;
	t104 = t124 * t131 * t136 + t114 * t139;
	t103 = (-t122 * t154 + t144 * t140) * t136 - t122 * t157;
	t100 = t108 * t139 + t143 * t136;
	t96 = 0.1e1 / (t106 * t98 ^ 2 + 0.1e1);
	t90 = 0.1e1 / t93;
	t89 = 0.1e1 / t148;
	t86 = 0.1e1 / t88;
	t85 = 0.1e1 / (t101 ^ 2 * t87 + 0.1e1);
	t84 = (-t105 * t107 + t118 * t160) * t96 * t136;
	t83 = (-t103 * t105 + t115 * t160) * t96;
	t82 = (-t100 * t105 + t112 * t160) * t96;
	t1 = [0, t83, t84, t82, 0, 0; 0, ((t114 * t136 - t124 * t157) * t86 - (-t103 * t94 + t115 * t95 + t146 * t83) * t161) * t85, (-t109 * t136 * t86 - (t146 * t84 + (-t107 * t94 + t118 * t95) * t136) * t161) * t85, (t102 * t86 - (-t100 * t94 + t112 * t95 + t146 * t82) * t161) * t85, 0, 0; 0, ((t104 * t127 - t113 * t128) * t90 - (t104 * t128 + t113 * t127) * t162) * t89, ((-t110 * t128 - t127 * t159) * t90 - (t110 * t127 - t128 * t159) * t162) * t89, (-t127 * t90 + t128 * t162) * t89 * t101, 0, t148 * t89;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end