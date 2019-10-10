% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPPR5
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
%   Wie in S6RPRPPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:22
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPPR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:18
	% EndTime: 2019-10-10 00:22:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:18
	% DurationCPUTime: 0.72s
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
	% StartTime: 2019-10-10 00:22:17
	% EndTime: 2019-10-10 00:22:19
	% DurationCPUTime: 1.01s
	% Computational Cost: add. (1897->83), mult. (2191->183), div. (456->12), fcn. (2616->9), ass. (0->85)
	t104 = pkin(9) + qJ(3);
	t103 = cos(t104);
	t101 = t103 ^ 2;
	t102 = sin(t104);
	t154 = t101 / t102 ^ 2;
	t109 = sin(qJ(1));
	t130 = 0.1e1 + t154;
	t105 = t109 ^ 2;
	t96 = t105 * t154 + 0.1e1;
	t94 = 0.1e1 / t96;
	t121 = t130 * t94;
	t77 = t109 * t121;
	t166 = t109 * t77 - 0.1e1;
	t163 = t103 * t154;
	t97 = 0.1e1 / t102;
	t119 = qJD(3) * (-t103 - t163) * t97;
	t110 = cos(qJ(1));
	t142 = qJD(1) * t110;
	t129 = t109 * t142;
	t165 = (t105 * t119 + t129 * t154) / t96 ^ 2;
	t107 = sin(pkin(10));
	t139 = qJD(3) * t110;
	t127 = t103 * t139;
	t146 = t109 * t107;
	t108 = cos(pkin(10));
	t148 = t108 * t110;
	t90 = -t102 * t146 + t148;
	t82 = t90 * qJD(1) + t107 * t127;
	t145 = t109 * t108;
	t149 = t107 * t110;
	t88 = t102 * t149 + t145;
	t85 = 0.1e1 / t88 ^ 2;
	t164 = t82 * t85;
	t147 = t109 * t103;
	t93 = atan2(-t147, t102);
	t91 = sin(t93);
	t133 = t91 * t147;
	t92 = cos(t93);
	t75 = t102 * t92 - t133;
	t72 = 0.1e1 / t75;
	t84 = 0.1e1 / t88;
	t73 = 0.1e1 / t75 ^ 2;
	t162 = t94 - 0.1e1;
	t140 = qJD(3) * t109;
	t153 = t102 * t91;
	t68 = ((t102 * t140 - t103 * t142) * t97 + t140 * t154) * t94;
	t63 = (-t68 + t140) * t153 + (-t91 * t142 + (-t109 * t68 + qJD(3)) * t92) * t103;
	t161 = t63 * t72 * t73;
	t106 = t110 ^ 2;
	t71 = t101 * t106 * t73 + 0.1e1;
	t69 = 0.1e1 / t71;
	t159 = t69 * t73;
	t158 = t72 * t69;
	t157 = t84 * t164;
	t120 = t102 * t148 - t146;
	t156 = t85 * t120;
	t155 = t101 * t97;
	t151 = t110 * t73;
	t150 = qJD(3) * t77;
	t144 = qJD(1) * t103;
	t143 = qJD(1) * t109;
	t141 = qJD(3) * t102;
	t131 = t73 * t141;
	t138 = 0.2e1 * (-t106 * t103 * t131 + (-t106 * t161 - t73 * t129) * t101) / t71 ^ 2;
	t137 = 0.2e1 * t161;
	t83 = t120 ^ 2;
	t80 = t83 * t85 + 0.1e1;
	t89 = t102 * t145 + t149;
	t81 = t89 * qJD(1) - t108 * t127;
	t136 = 0.2e1 * (-t81 * t156 - t83 * t157) / t80 ^ 2;
	t135 = 0.2e1 * t165;
	t134 = t109 * t94 * t97;
	t132 = t103 * t162;
	t128 = t103 * t140;
	t126 = t72 * t138;
	t125 = t73 * t138;
	t124 = -0.2e1 * t120 * t157;
	t123 = t103 * t135;
	t122 = t101 * t134;
	t78 = 0.1e1 / t80;
	t118 = (-t107 * t156 + t108 * t84) * t78;
	t67 = (t92 * t122 + t91 * t132) * t110;
	t65 = -t166 * t92 * t103 + (t109 - t77) * t153;
	t64 = t121 * t142 + 0.2e1 * (t119 * t94 - t130 * t165) * t109;
	t1 = [t134 * t144 + (qJD(3) * t121 + t97 * t123) * t110, 0, t64, 0, 0, 0; (t141 * t158 + (t126 + (qJD(1) * t67 + t63) * t159) * t103) * t109 + (t67 * t125 * t103 + (t67 * t131 + (t67 * t137 + ((t68 * t122 + t162 * t141 + t123) * t91 + (-t68 * t132 + (t135 * t155 + (0.2e1 * t103 + t163) * t94 * qJD(3)) * t109) * t92) * t151) * t103 + (-t72 + (t162 * t133 - (-t105 + t106) * t94 * t92 * t155) * t73) * t144) * t69) * t110, 0, (t143 * t158 + (t126 + (qJD(3) * t65 + t63) * t159) * t110) * t102 + (t65 * t110 * t125 + (-t72 * t139 + (t110 * t137 + t73 * t143) * t65 + (-((-t109 * t64 - t142 * t77) * t92 + (t166 * t68 + t140 - t150) * t91) * t103 - ((-t64 + t142) * t91 + (-t68 * t77 - qJD(3) + (t68 + t150) * t109) * t92) * t102) * t151) * t69) * t103, 0, 0, 0; (-t90 * t156 - t84 * t89) * t136 + ((t120 * qJD(1) + t108 * t128) * t84 + t90 * t124 + (-t89 * t82 + (-t88 * qJD(1) - t107 * t128) * t120 - t90 * t81) * t85) * t78, 0, t102 * t118 * t139 + (t118 * t143 + ((t84 * t136 + t78 * t164) * t108 + (-t136 * t156 + (-t81 * t85 + t124) * t78) * t107) * t110) * t103, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:22:18
	% EndTime: 2019-10-10 00:22:19
	% DurationCPUTime: 1.10s
	% Computational Cost: add. (2429->94), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->95)
	t143 = pkin(9) + qJ(3);
	t139 = sin(t143);
	t134 = 0.1e1 / t139 ^ 2;
	t141 = cos(t143);
	t137 = t141 ^ 2;
	t191 = t134 * t137;
	t210 = t141 * t191;
	t146 = sin(qJ(1));
	t167 = 0.1e1 + t191;
	t209 = t146 * t167;
	t163 = qJD(1) * t139 + qJD(6);
	t147 = cos(qJ(1));
	t179 = qJD(3) * t147;
	t208 = -t141 * t179 + t163 * t146;
	t180 = qJD(3) * t146;
	t207 = t141 * t180 + t163 * t147;
	t184 = t146 * t141;
	t128 = atan2(-t184, t139);
	t127 = cos(t128);
	t126 = sin(t128);
	t172 = t126 * t184;
	t112 = t127 * t139 - t172;
	t109 = 0.1e1 / t112;
	t142 = pkin(10) + qJ(6);
	t138 = sin(t142);
	t140 = cos(t142);
	t185 = t146 * t140;
	t187 = t139 * t147;
	t123 = t138 * t187 + t185;
	t119 = 0.1e1 / t123;
	t133 = 0.1e1 / t139;
	t110 = 0.1e1 / t112 ^ 2;
	t120 = 0.1e1 / t123 ^ 2;
	t144 = t146 ^ 2;
	t131 = t144 * t191 + 0.1e1;
	t129 = 0.1e1 / t131;
	t206 = t129 - 0.1e1;
	t182 = qJD(1) * t147;
	t170 = t141 * t182;
	t103 = ((t139 * t180 - t170) * t133 + t180 * t191) * t129;
	t193 = t127 * t141;
	t98 = (-t103 * t146 + qJD(3)) * t193 + (-t170 + (-t103 + t180) * t139) * t126;
	t205 = t109 * t110 * t98;
	t164 = qJD(6) * t139 + qJD(1);
	t158 = t164 * t147;
	t104 = t138 * t158 + t208 * t140;
	t186 = t140 * t147;
	t122 = t138 * t146 - t139 * t186;
	t118 = t122 ^ 2;
	t117 = t118 * t120 + 0.1e1;
	t195 = t120 * t122;
	t105 = -t208 * t138 + t140 * t158;
	t201 = t105 * t119 * t120;
	t204 = (t104 * t195 - t118 * t201) / t117 ^ 2;
	t203 = t103 * t126;
	t202 = t103 * t141;
	t200 = t110 * t141;
	t199 = t110 * t147;
	t192 = t133 * t141;
	t156 = qJD(3) * (-t133 * t210 - t192);
	t189 = t137 * t146;
	t161 = t182 * t189;
	t198 = (t134 * t161 + t144 * t156) / t131 ^ 2;
	t116 = t129 * t209;
	t197 = t116 * t146;
	t196 = t119 * t140;
	t194 = t122 * t138;
	t145 = t147 ^ 2;
	t190 = t137 * t145;
	t188 = t139 * t146;
	t183 = qJD(1) * t146;
	t181 = qJD(3) * t139;
	t108 = t110 * t190 + 0.1e1;
	t178 = 0.2e1 / t108 ^ 2 * (-t190 * t205 + (-t141 * t145 * t181 - t161) * t110);
	t177 = 0.2e1 * t205;
	t176 = 0.2e1 * t204;
	t175 = -0.2e1 * t198;
	t174 = t141 * t199;
	t173 = t141 * t198;
	t171 = t133 * t189;
	t166 = t141 * t178;
	t165 = 0.2e1 * t122 * t201;
	t162 = t129 * t171;
	t160 = t167 * t147;
	t159 = t164 * t146;
	t157 = t120 * t194 + t196;
	t155 = t157 * t147;
	t125 = -t138 * t188 + t186;
	t124 = t138 * t147 + t139 * t185;
	t114 = 0.1e1 / t117;
	t106 = 0.1e1 / t108;
	t102 = (t206 * t141 * t126 + t127 * t162) * t147;
	t101 = t126 * t188 + t193 + (-t126 * t139 - t127 * t184) * t116;
	t99 = t175 * t209 + (qJD(1) * t160 + 0.2e1 * t146 * t156) * t129;
	t1 = [0.2e1 * t133 * t147 * t173 + (qJD(3) * t160 + t183 * t192) * t129, 0, t99, 0, 0, 0; (t109 * t166 + (t109 * t181 + (qJD(1) * t102 + t98) * t200) * t106) * t146 + (t110 * t166 * t102 + (-((-t103 * t162 - t206 * t181 - 0.2e1 * t173) * t126 + (t171 * t175 - t202 + (t202 + (-0.2e1 * t141 - t210) * t180) * t129) * t127) * t174 + (t110 * t181 + t141 * t177) * t102 + (-t109 + ((t144 - t145) * t137 * t133 * t129 * t127 + t206 * t172) * t110) * t141 * qJD(1)) * t106) * t147, 0, (t101 * t200 + t109 * t139) * t147 * t178 + ((t109 * t183 + (qJD(3) * t101 + t98) * t199) * t139 + (-t109 * t179 - (-t127 * t146 * t99 + t126 * t180 + t197 * t203 - t203 + (-qJD(3) * t126 - t127 * t182) * t116) * t174 + (t110 * t183 + t147 * t177) * t101 - ((-t99 + t182) * t126 + ((-0.1e1 + t197) * qJD(3) + (-t116 + t146) * t103) * t127) * t110 * t187) * t141) * t106, 0, 0, 0; (-t119 * t124 + t125 * t195) * t176 + (t125 * t165 - t119 * t138 * t159 + t207 * t196 + (t122 * t140 * t159 - t125 * t104 - t124 * t105 + t207 * t194) * t120) * t114, 0, t141 * t155 * t176 + (t155 * t181 + (t157 * t183 + ((qJD(6) * t119 + t165) * t138 + (-t104 * t138 + (-qJD(6) * t122 + t105) * t140) * t120) * t147) * t141) * t114, 0, 0, -0.2e1 * t204 + 0.2e1 * (t104 * t114 * t120 + (-t114 * t201 - t120 * t204) * t122) * t122;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end