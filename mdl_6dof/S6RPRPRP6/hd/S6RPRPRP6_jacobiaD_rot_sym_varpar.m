% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRP6
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
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RPRPRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:42
	% DurationCPUTime: 0.78s
	% Computational Cost: add. (1883->70), mult. (1839->155), div. (436->14), fcn. (2165->7), ass. (0->75)
	t92 = sin(qJ(1));
	t142 = 0.2e1 * t92;
	t87 = 0.1e1 / t92;
	t93 = cos(qJ(1));
	t126 = t87 * t93;
	t120 = qJD(3) * t92;
	t122 = qJD(1) * t93;
	t85 = pkin(9) + qJ(3);
	t83 = sin(t85);
	t112 = t83 * t122;
	t78 = t83 ^ 2;
	t84 = cos(t85);
	t81 = 0.1e1 / t84 ^ 2;
	t129 = t78 * t81;
	t86 = t92 ^ 2;
	t73 = t86 * t129 + 0.1e1;
	t71 = 0.1e1 / t73;
	t80 = 0.1e1 / t84;
	t57 = (-(-t84 * t120 - t112) * t80 + t120 * t129) * t71;
	t141 = t57 - t120;
	t91 = t93 ^ 2;
	t140 = qJD(1) * (0.1e1 / t86 * t91 + 0.1e1) * t126;
	t124 = t92 * t83;
	t70 = atan2(-t124, -t84);
	t68 = sin(t70);
	t114 = t68 * t124;
	t69 = cos(t70);
	t64 = -t69 * t84 - t114;
	t61 = 0.1e1 / t64;
	t62 = 0.1e1 / t64 ^ 2;
	t139 = -0.2e1 * t83;
	t138 = t71 - 0.1e1;
	t131 = t69 * t83;
	t53 = (-t57 * t92 + qJD(3)) * t131 + (t141 * t84 - t112) * t68;
	t137 = t53 * t61 * t62;
	t136 = t57 * t83;
	t135 = t62 * t83;
	t134 = t62 * t93;
	t127 = t80 * t83;
	t77 = t83 * t78;
	t79 = t84 ^ 2;
	t100 = qJD(3) * (t77 * t80 / t79 + t127);
	t105 = t78 * t92 * t122;
	t133 = (t86 * t100 + t81 * t105) / t73 ^ 2;
	t132 = t68 * t92;
	t130 = t78 * t80;
	t128 = t78 * t91;
	t88 = 0.1e1 / t92 ^ 2;
	t125 = t88 * t91;
	t123 = qJD(1) * t92;
	t121 = qJD(3) * t84;
	t119 = qJD(3) * t93;
	t104 = t83 * t91 * t121;
	t60 = t62 * t128 + 0.1e1;
	t118 = 0.2e1 * (-t128 * t137 + (t104 - t105) * t62) / t60 ^ 2;
	t117 = 0.2e1 * t137;
	t76 = t79 * t125 + 0.1e1;
	t116 = 0.2e1 * (-t88 * t104 - t79 * t140) / t76 ^ 2;
	t115 = t83 * t134;
	t113 = t71 * t130;
	t111 = 0.1e1 + t129;
	t110 = 0.1e1 + t125;
	t109 = t83 * t118;
	t108 = t133 * t139;
	t107 = t133 * t142;
	t106 = t92 * t113;
	t103 = t111 * t93;
	t102 = t110 * t83;
	t74 = 0.1e1 / t76;
	t66 = t111 * t92 * t71;
	t58 = 0.1e1 / t60;
	t56 = (t138 * t83 * t68 - t69 * t106) * t93;
	t55 = -t84 * t132 + t131 + (-t69 * t124 + t68 * t84) * t66;
	t54 = -t111 * t107 + (qJD(1) * t103 + t100 * t142) * t71;
	t1 = [t93 * t80 * t108 + (qJD(3) * t103 - t123 * t127) * t71, 0, t54, 0, 0, 0; (t61 * t109 + (-t61 * t121 + (qJD(1) * t56 + t53) * t135) * t58) * t92 + (t62 * t109 * t56 + (-((t57 * t106 + t138 * t121 + t108) * t68 + (t107 * t130 - t136 + (t136 + (-t77 * t81 + t139) * t120) * t71) * t69) * t115 + (t83 * t117 - t62 * t121) * t56 + (-t61 + ((-t86 + t91) * t69 * t113 + t138 * t114) * t62) * t83 * qJD(1)) * t58) * t93, 0, (t55 * t135 - t61 * t84) * t93 * t118 + ((-t61 * t123 + (-qJD(3) * t55 - t53) * t134) * t84 + (-t61 * t119 - (-t54 * t69 * t92 - t141 * t68 + (-qJD(3) * t68 - t122 * t69 + t132 * t57) * t66) * t115 + (t93 * t117 + t62 * t123) * t55 - ((t54 - t122) * t68 + ((-t66 * t92 + 0.1e1) * qJD(3) + (t66 - t92) * t57) * t69) * t84 * t134) * t83) * t58, 0, 0, 0; t110 * t84 * t116 + (qJD(3) * t102 + 0.2e1 * t84 * t140) * t74, 0, t83 * t116 * t126 + (-t84 * t87 * t119 + qJD(1) * t102) * t74, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:42
	% DurationCPUTime: 1.02s
	% Computational Cost: add. (2081->92), mult. (2519->204), div. (480->12), fcn. (2968->9), ass. (0->90)
	t135 = sin(qJ(1));
	t137 = cos(qJ(1));
	t131 = pkin(9) + qJ(3);
	t129 = sin(t131);
	t152 = qJD(1) * t129 + qJD(5);
	t130 = cos(t131);
	t172 = qJD(3) * t130;
	t197 = t152 * t135 - t137 * t172;
	t178 = t135 * t130;
	t119 = atan2(-t178, t129);
	t118 = cos(t119);
	t117 = sin(t119);
	t165 = t117 * t178;
	t103 = t118 * t129 - t165;
	t100 = 0.1e1 / t103;
	t136 = cos(qJ(5));
	t177 = t135 * t136;
	t134 = sin(qJ(5));
	t179 = t134 * t137;
	t114 = t129 * t179 + t177;
	t110 = 0.1e1 / t114;
	t124 = 0.1e1 / t129;
	t101 = 0.1e1 / t103 ^ 2;
	t111 = 0.1e1 / t114 ^ 2;
	t125 = 0.1e1 / t129 ^ 2;
	t132 = t135 ^ 2;
	t128 = t130 ^ 2;
	t183 = t125 * t128;
	t122 = t132 * t183 + 0.1e1;
	t120 = 0.1e1 / t122;
	t196 = t120 - 0.1e1;
	t133 = t137 ^ 2;
	t182 = t128 * t133;
	t97 = t101 * t182 + 0.1e1;
	t95 = 0.1e1 / t97;
	t195 = t101 * t95;
	t171 = qJD(3) * t135;
	t162 = t125 * t171;
	t174 = qJD(1) * t137;
	t163 = t130 * t174;
	t94 = ((t129 * t171 - t163) * t124 + t128 * t162) * t120;
	t154 = -t94 + t171;
	t155 = -t135 * t94 + qJD(3);
	t186 = t118 * t130;
	t89 = t155 * t186 + (t154 * t129 - t163) * t117;
	t194 = t100 * t101 * t89;
	t176 = t136 * t137;
	t180 = t134 * t135;
	t113 = -t129 * t176 + t180;
	t109 = t113 ^ 2;
	t108 = t109 * t111 + 0.1e1;
	t189 = t111 * t113;
	t153 = qJD(5) * t129 + qJD(1);
	t148 = t153 * t137;
	t99 = -t197 * t134 + t136 * t148;
	t192 = t110 * t111 * t99;
	t98 = t134 * t148 + t197 * t136;
	t193 = 0.1e1 / t108 ^ 2 * (-t109 * t192 + t98 * t189);
	t188 = t113 * t134;
	t187 = t117 * t129;
	t185 = t124 * t128;
	t184 = t124 * t130;
	t175 = qJD(1) * t135;
	t173 = qJD(3) * t129;
	t150 = t128 * t135 * t174;
	t170 = 0.2e1 * (-t182 * t194 + (-t129 * t133 * t172 - t150) * t101) / t97 ^ 2;
	t169 = 0.2e1 * t194;
	t168 = 0.2e1 * t193;
	t127 = t130 * t128;
	t146 = qJD(3) * (-t124 * t125 * t127 - t184);
	t167 = 0.2e1 * (t125 * t150 + t132 * t146) / t122 ^ 2;
	t166 = t95 * t173;
	t164 = t120 * t185;
	t160 = 0.1e1 + t183;
	t159 = t100 * t170;
	t158 = 0.2e1 * t113 * t192;
	t157 = t130 * t167;
	t156 = t135 * t167;
	t151 = t135 * t164;
	t149 = t160 * t137;
	t147 = t110 * t136 + t111 * t188;
	t145 = t147 * t137;
	t116 = -t129 * t180 + t176;
	t115 = t129 * t177 + t179;
	t106 = 0.1e1 / t108;
	t105 = t160 * t135 * t120;
	t93 = (t196 * t130 * t117 + t118 * t151) * t137;
	t91 = t135 * t187 + t186 + (-t118 * t178 - t187) * t105;
	t90 = -t160 * t156 + (qJD(1) * t149 + 0.2e1 * t135 * t146) * t120;
	t1 = [t124 * t137 * t157 + (qJD(3) * t149 + t175 * t184) * t120, 0, t90, 0, 0, 0; (t100 * t166 + (t159 + (qJD(1) * t93 + t89) * t195) * t130) * t135 + ((t93 * t166 + (t93 * t170 + ((t94 * t151 + t196 * t173 + t157) * t117 + (t156 * t185 + t130 * t94 + (t127 * t162 - (t94 - 0.2e1 * t171) * t130) * t120) * t118) * t95 * t137) * t130) * t101 + (t93 * t169 + (-t100 + ((t132 - t133) * t118 * t164 + t196 * t165) * t101) * qJD(1)) * t130 * t95) * t137, 0, (t100 * t95 * t175 + (t159 + (qJD(3) * t91 + t89) * t195) * t137) * t129 + (((-qJD(3) * t100 + t91 * t169) * t137 + (t91 * t175 + (-(-t105 * t174 - t135 * t90) * t118 - ((t105 * t135 - 0.1e1) * t94 + (-t105 + t135) * qJD(3)) * t117) * t130 * t137) * t101) * t95 + (t91 * t170 - ((-t90 + t174) * t117 + (t154 * t105 - t155) * t118) * t95 * t129) * t101 * t137) * t130, 0, 0, 0; (-t110 * t115 + t116 * t189) * t168 + (t116 * t158 + (-t115 * t99 - t116 * t98 + t153 * t113 * t177 - (-t130 * t171 - t152 * t137) * t188) * t111 + (t152 * t176 + (-t153 * t134 + t136 * t172) * t135) * t110) * t106, 0, t130 * t145 * t168 + (t145 * t173 + (t147 * t175 + ((qJD(5) * t110 + t158) * t134 + (-t134 * t98 + (-qJD(5) * t113 + t99) * t136) * t111) * t137) * t130) * t106, 0, -0.2e1 * t193 + 0.2e1 * (t106 * t111 * t98 + (-t106 * t192 - t111 * t193) * t113) * t113, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:42
	% DurationCPUTime: 0.99s
	% Computational Cost: add. (2081->92), mult. (2519->204), div. (480->12), fcn. (2968->9), ass. (0->90)
	t139 = sin(qJ(1));
	t141 = cos(qJ(1));
	t135 = pkin(9) + qJ(3);
	t133 = sin(t135);
	t156 = qJD(1) * t133 + qJD(5);
	t134 = cos(t135);
	t176 = qJD(3) * t134;
	t201 = t156 * t139 - t141 * t176;
	t183 = t139 * t134;
	t123 = atan2(-t183, t133);
	t122 = cos(t123);
	t121 = sin(t123);
	t169 = t121 * t183;
	t107 = t122 * t133 - t169;
	t104 = 0.1e1 / t107;
	t138 = sin(qJ(5));
	t181 = t141 * t138;
	t140 = cos(qJ(5));
	t182 = t139 * t140;
	t118 = t133 * t181 + t182;
	t114 = 0.1e1 / t118;
	t128 = 0.1e1 / t133;
	t105 = 0.1e1 / t107 ^ 2;
	t115 = 0.1e1 / t118 ^ 2;
	t129 = 0.1e1 / t133 ^ 2;
	t136 = t139 ^ 2;
	t132 = t134 ^ 2;
	t187 = t129 * t132;
	t126 = t136 * t187 + 0.1e1;
	t124 = 0.1e1 / t126;
	t200 = t124 - 0.1e1;
	t137 = t141 ^ 2;
	t186 = t132 * t137;
	t101 = t105 * t186 + 0.1e1;
	t99 = 0.1e1 / t101;
	t199 = t105 * t99;
	t175 = qJD(3) * t139;
	t166 = t129 * t175;
	t178 = qJD(1) * t141;
	t167 = t134 * t178;
	t98 = ((t133 * t175 - t167) * t128 + t132 * t166) * t124;
	t158 = -t98 + t175;
	t159 = -t139 * t98 + qJD(3);
	t190 = t122 * t134;
	t93 = t159 * t190 + (t158 * t133 - t167) * t121;
	t198 = t104 * t105 * t93;
	t157 = qJD(5) * t133 + qJD(1);
	t152 = t157 * t141;
	t102 = t138 * t152 + t201 * t140;
	t180 = t141 * t140;
	t184 = t138 * t139;
	t117 = -t133 * t180 + t184;
	t113 = t117 ^ 2;
	t112 = t113 * t115 + 0.1e1;
	t193 = t115 * t117;
	t103 = -t201 * t138 + t140 * t152;
	t195 = t103 * t114 * t115;
	t197 = 0.1e1 / t112 ^ 2 * (t102 * t193 - t113 * t195);
	t192 = t117 * t138;
	t191 = t121 * t133;
	t189 = t128 * t132;
	t188 = t128 * t134;
	t179 = qJD(1) * t139;
	t177 = qJD(3) * t133;
	t154 = t132 * t139 * t178;
	t174 = 0.2e1 / t101 ^ 2 * (-t186 * t198 + (-t133 * t137 * t176 - t154) * t105);
	t173 = 0.2e1 * t198;
	t172 = 0.2e1 * t197;
	t131 = t134 * t132;
	t150 = qJD(3) * (-t128 * t129 * t131 - t188);
	t171 = 0.2e1 * (t129 * t154 + t136 * t150) / t126 ^ 2;
	t170 = t99 * t177;
	t168 = t124 * t189;
	t164 = 0.1e1 + t187;
	t163 = t104 * t174;
	t162 = t134 * t171;
	t161 = t139 * t171;
	t160 = 0.2e1 * t117 * t195;
	t155 = t139 * t168;
	t153 = t164 * t141;
	t151 = t114 * t140 + t115 * t192;
	t149 = t151 * t141;
	t120 = -t133 * t184 + t180;
	t119 = t133 * t182 + t181;
	t110 = 0.1e1 / t112;
	t109 = t164 * t139 * t124;
	t97 = (t200 * t134 * t121 + t122 * t155) * t141;
	t95 = t139 * t191 + t190 + (-t122 * t183 - t191) * t109;
	t94 = -t164 * t161 + (qJD(1) * t153 + 0.2e1 * t139 * t150) * t124;
	t1 = [t128 * t141 * t162 + (qJD(3) * t153 + t179 * t188) * t124, 0, t94, 0, 0, 0; (t104 * t170 + (t163 + (qJD(1) * t97 + t93) * t199) * t134) * t139 + ((t97 * t170 + (t97 * t174 + ((t98 * t155 + t200 * t177 + t162) * t121 + (t161 * t189 + t134 * t98 + (t131 * t166 - (t98 - 0.2e1 * t175) * t134) * t124) * t122) * t99 * t141) * t134) * t105 + (t97 * t173 + (-t104 + ((t136 - t137) * t122 * t168 + t200 * t169) * t105) * qJD(1)) * t134 * t99) * t141, 0, (t104 * t99 * t179 + (t163 + (qJD(3) * t95 + t93) * t199) * t141) * t133 + (((-qJD(3) * t104 + t95 * t173) * t141 + (t95 * t179 + (-(-t109 * t178 - t139 * t94) * t122 - ((t109 * t139 - 0.1e1) * t98 + (-t109 + t139) * qJD(3)) * t121) * t134 * t141) * t105) * t99 + (t95 * t174 - ((-t94 + t178) * t121 + (t158 * t109 - t159) * t122) * t99 * t133) * t105 * t141) * t134, 0, 0, 0; (-t114 * t119 + t120 * t193) * t172 + (t120 * t160 + (-t120 * t102 - t119 * t103 + t157 * t117 * t182 - (-t134 * t175 - t156 * t141) * t192) * t115 + (t156 * t180 + (-t157 * t138 + t140 * t176) * t139) * t114) * t110, 0, t134 * t149 * t172 + (t149 * t177 + (t151 * t179 + ((qJD(5) * t114 + t160) * t138 + (-t102 * t138 + (-qJD(5) * t117 + t103) * t140) * t115) * t141) * t134) * t110, 0, -0.2e1 * t197 + 0.2e1 * (t102 * t110 * t115 + (-t110 * t195 - t115 * t197) * t117) * t117, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end