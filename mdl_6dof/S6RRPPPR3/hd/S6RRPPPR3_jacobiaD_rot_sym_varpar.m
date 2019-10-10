% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR3
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
%   Wie in S6RRPPPR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:21
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:36
	% EndTime: 2019-10-10 09:21:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:36
	% EndTime: 2019-10-10 09:21:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:36
	% EndTime: 2019-10-10 09:21:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:36
	% EndTime: 2019-10-10 09:21:37
	% DurationCPUTime: 0.80s
	% Computational Cost: add. (776->71), mult. (1835->159), div. (470->13), fcn. (2177->7), ass. (0->71)
	t82 = sin(qJ(1));
	t113 = qJD(1) * t82;
	t133 = 0.2e1 * t82;
	t73 = t82 ^ 2;
	t84 = cos(qJ(1));
	t77 = t84 ^ 2;
	t78 = 0.1e1 / t84;
	t131 = (t73 / t77 + 0.1e1) * t78 * t113;
	t81 = sin(qJ(2));
	t114 = t82 * t81;
	t83 = cos(qJ(2));
	t63 = atan2(-t114, -t83);
	t61 = sin(t63);
	t105 = t61 * t114;
	t62 = cos(t63);
	t57 = -t62 * t83 - t105;
	t54 = 0.1e1 / t57;
	t74 = 0.1e1 / t83;
	t55 = 0.1e1 / t57 ^ 2;
	t75 = 0.1e1 / t83 ^ 2;
	t130 = -0.2e1 * t81;
	t71 = t81 ^ 2;
	t118 = t71 * t75;
	t68 = t73 * t118 + 0.1e1;
	t64 = 0.1e1 / t68;
	t129 = t64 - 0.1e1;
	t112 = qJD(1) * t84;
	t103 = t81 * t112;
	t111 = qJD(2) * t82;
	t120 = t62 * t81;
	t110 = qJD(2) * t83;
	t50 = (-(-t82 * t110 - t103) * t74 + t111 * t118) * t64;
	t46 = (-t50 * t82 + qJD(2)) * t120 + (-t103 + (t50 - t111) * t83) * t61;
	t128 = t46 * t54 * t55;
	t127 = t50 * t61;
	t126 = t50 * t81;
	t125 = t55 * t81;
	t124 = t55 * t84;
	t115 = t74 * t81;
	t70 = t81 * t71;
	t76 = t74 * t75;
	t92 = qJD(2) * (t70 * t76 + t115);
	t96 = t71 * t82 * t112;
	t123 = (t73 * t92 + t75 * t96) / t68 ^ 2;
	t102 = 0.1e1 + t118;
	t60 = t102 * t82 * t64;
	t122 = t60 * t82;
	t121 = t61 * t83;
	t119 = t71 * t74;
	t117 = t71 * t77;
	t116 = t73 / t84 ^ 2;
	t53 = t55 * t117 + 0.1e1;
	t109 = 0.2e1 * (-t117 * t128 + (t77 * t81 * t110 - t96) * t55) / t53 ^ 2;
	t108 = 0.2e1 * t128;
	t69 = t75 * t116 + 0.1e1;
	t107 = 0.2e1 * (t76 * qJD(2) * t81 * t116 + t75 * t131) / t69 ^ 2;
	t106 = t81 * t124;
	t104 = t64 * t119;
	t101 = 0.1e1 + t116;
	t100 = t81 * t109;
	t99 = t123 * t130;
	t98 = t123 * t133;
	t97 = t82 * t104;
	t95 = t102 * t84;
	t93 = t101 * t81 * t75;
	t66 = 0.1e1 / t69;
	t51 = 0.1e1 / t53;
	t49 = (t129 * t81 * t61 - t62 * t97) * t84;
	t48 = -t82 * t121 + t120 + (-t62 * t114 + t121) * t60;
	t47 = -t102 * t98 + (qJD(1) * t95 + t92 * t133) * t64;
	t1 = [t84 * t74 * t99 + (qJD(2) * t95 - t113 * t115) * t64, t47, 0, 0, 0, 0; (t54 * t100 + (-t54 * t110 + (qJD(1) * t49 + t46) * t125) * t51) * t82 + (t55 * t100 * t49 + (-((t129 * t110 + t50 * t97 + t99) * t61 + (t98 * t119 - t126 + (t126 + (-t70 * t75 + t130) * t111) * t64) * t62) * t106 + (t81 * t108 - t55 * t110) * t49 + (-t54 + ((-t73 + t77) * t62 * t104 + t129 * t105) * t55) * t81 * qJD(1)) * t51) * t84, (t48 * t125 - t54 * t83) * t84 * t109 + ((-t54 * t113 + (-qJD(2) * t48 - t46) * t124) * t83 + (-t84 * qJD(2) * t54 - (-t47 * t62 * t82 + t61 * t111 + t122 * t127 - t127 + (-qJD(2) * t61 - t112 * t62) * t60) * t106 + (t84 * t108 + t55 * t113) * t48 - ((t47 - t112) * t61 + ((0.1e1 - t122) * qJD(2) + (t60 - t82) * t50) * t62) * t83 * t124) * t81) * t51, 0, 0, 0, 0; t101 * t74 * t107 + (-qJD(2) * t93 - 0.2e1 * t74 * t131) * t66, t78 * t75 * t107 * t114 + ((-0.2e1 * t71 * t76 - t74) * t78 * t111 - qJD(1) * t93) * t66, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:36
	% EndTime: 2019-10-10 09:21:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:37
	% EndTime: 2019-10-10 09:21:37
	% DurationCPUTime: 0.94s
	% Computational Cost: add. (703->82), mult. (2191->188), div. (456->12), fcn. (2616->9), ass. (0->87)
	t101 = sin(qJ(2));
	t93 = 0.1e1 / t101 ^ 2;
	t103 = cos(qJ(2));
	t97 = t103 ^ 2;
	t149 = t93 * t97;
	t102 = sin(qJ(1));
	t125 = 0.1e1 + t149;
	t95 = t102 ^ 2;
	t91 = t95 * t149 + 0.1e1;
	t89 = 0.1e1 / t91;
	t114 = t125 * t89;
	t75 = t102 * t114;
	t161 = t102 * t75 - 0.1e1;
	t158 = t103 * t149;
	t92 = 0.1e1 / t101;
	t113 = qJD(2) * (-t103 - t158) * t92;
	t104 = cos(qJ(1));
	t136 = qJD(1) * t104;
	t115 = t102 * t97 * t136;
	t160 = (t95 * t113 + t93 * t115) / t91 ^ 2;
	t100 = cos(pkin(9));
	t132 = qJD(2) * t104;
	t121 = t103 * t132;
	t141 = t102 * t100;
	t99 = sin(pkin(9));
	t83 = -t101 * t141 - t104 * t99;
	t77 = t83 * qJD(1) + t100 * t121;
	t139 = t104 * t100;
	t145 = t102 * t99;
	t85 = t101 * t139 - t145;
	t80 = 0.1e1 / t85 ^ 2;
	t159 = t77 * t80;
	t140 = t102 * t103;
	t88 = atan2(-t140, t101);
	t86 = sin(t88);
	t126 = t86 * t140;
	t87 = cos(t88);
	t70 = t87 * t101 - t126;
	t67 = 0.1e1 / t70;
	t79 = 0.1e1 / t85;
	t68 = 0.1e1 / t70 ^ 2;
	t157 = t89 - 0.1e1;
	t134 = qJD(2) * t102;
	t147 = t101 * t86;
	t63 = ((t101 * t134 - t103 * t136) * t92 + t134 * t149) * t89;
	t58 = (-t63 + t134) * t147 + (-t86 * t136 + (-t102 * t63 + qJD(2)) * t87) * t103;
	t156 = t58 * t67 * t68;
	t142 = t101 * t104;
	t84 = t99 * t142 + t141;
	t151 = t80 * t84;
	t152 = t79 * t159;
	t78 = t84 ^ 2;
	t73 = t78 * t80 + 0.1e1;
	t82 = -t101 * t145 + t139;
	t76 = t82 * qJD(1) + t99 * t121;
	t155 = (t76 * t151 - t78 * t152) / t73 ^ 2;
	t98 = t104 ^ 2;
	t148 = t97 * t98;
	t66 = t68 * t148 + 0.1e1;
	t64 = 0.1e1 / t66;
	t153 = t64 * t68;
	t150 = t92 * t97;
	t144 = t104 * t68;
	t143 = qJD(2) * t75;
	t138 = qJD(1) * t102;
	t137 = qJD(1) * t103;
	t135 = qJD(2) * t101;
	t133 = qJD(2) * t103;
	t131 = 0.2e1 * (-t148 * t156 + (-t101 * t98 * t133 - t115) * t68) / t66 ^ 2;
	t130 = 0.2e1 * t156;
	t129 = 0.2e1 * t155;
	t128 = 0.2e1 * t160;
	t127 = t102 * t89 * t92;
	t124 = t157 * t103;
	t123 = t64 * t135;
	t122 = t102 * t133;
	t120 = t67 * t131;
	t119 = t68 * t131;
	t118 = 0.2e1 * t84 * t152;
	t117 = t103 * t128;
	t116 = t97 * t127;
	t71 = 0.1e1 / t73;
	t112 = (t100 * t151 - t79 * t99) * t71;
	t62 = (t87 * t116 + t86 * t124) * t104;
	t60 = -t161 * t87 * t103 + (t102 - t75) * t147;
	t59 = t114 * t136 + 0.2e1 * (t113 * t89 - t125 * t160) * t102;
	t1 = [t127 * t137 + (qJD(2) * t114 + t92 * t117) * t104, t59, 0, 0, 0, 0; (t67 * t123 + (t120 + (qJD(1) * t62 + t58) * t153) * t103) * t102 + (t62 * t68 * t123 + (t62 * t119 + (t62 * t130 + ((t63 * t116 + t157 * t135 + t117) * t86 + (-t63 * t124 + (t128 * t150 + (0.2e1 * t103 + t158) * t89 * qJD(2)) * t102) * t87) * t144) * t64) * t103 + (-t67 + (-(-t95 + t98) * t89 * t87 * t150 + t157 * t126) * t68) * t64 * t137) * t104, (t67 * t64 * t138 + (t120 + (qJD(2) * t60 + t58) * t153) * t104) * t101 + (t60 * t104 * t119 + (-t67 * t132 - ((-t102 * t59 - t136 * t75) * t87 + (t161 * t63 + t134 - t143) * t86) * t103 * t144 + (t104 * t130 + t68 * t138) * t60) * t64 - ((-t59 + t136) * t86 + (-t63 * t75 - qJD(2) + (t63 + t143) * t102) * t87) * t142 * t153) * t103, 0, 0, 0, 0; (t83 * t151 - t79 * t82) * t129 + ((-t84 * qJD(1) - t99 * t122) * t79 + t83 * t118 + (-t82 * t77 - (-t85 * qJD(1) - t100 * t122) * t84 - t83 * t76) * t80) * t71, t101 * t112 * t132 + (t112 * t138 + ((-0.2e1 * t79 * t155 - t71 * t159) * t99 + (t129 * t151 + (-t76 * t80 + t118) * t71) * t100) * t104) * t103, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:21:37
	% EndTime: 2019-10-10 09:21:38
	% DurationCPUTime: 0.95s
	% Computational Cost: add. (1161->93), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->93)
	t135 = sin(qJ(2));
	t136 = sin(qJ(1));
	t137 = cos(qJ(2));
	t181 = t136 * t137;
	t120 = atan2(-t181, t135);
	t117 = cos(t120);
	t116 = sin(t120);
	t167 = t116 * t181;
	t105 = t117 * t135 - t167;
	t102 = 0.1e1 / t105;
	t127 = pkin(9) + qJ(6);
	t126 = cos(t127);
	t138 = cos(qJ(1));
	t184 = t135 * t138;
	t165 = t126 * t184;
	t125 = sin(t127);
	t183 = t136 * t125;
	t115 = t165 - t183;
	t109 = 0.1e1 / t115;
	t128 = 0.1e1 / t135;
	t103 = 0.1e1 / t105 ^ 2;
	t110 = 0.1e1 / t115 ^ 2;
	t129 = 0.1e1 / t135 ^ 2;
	t131 = t136 ^ 2;
	t133 = t137 ^ 2;
	t186 = t129 * t133;
	t123 = t131 * t186 + 0.1e1;
	t121 = 0.1e1 / t123;
	t198 = t121 - 0.1e1;
	t182 = t136 * t126;
	t114 = t125 * t184 + t182;
	t108 = t114 ^ 2;
	t101 = t108 * t110 + 0.1e1;
	t192 = t110 * t114;
	t153 = qJD(1) * t135 + qJD(6);
	t154 = qJD(6) * t135 + qJD(1);
	t173 = qJD(2) * t138;
	t162 = t137 * t173;
	t188 = t125 * t138;
	t94 = -t154 * t188 + (-t153 * t136 + t162) * t126;
	t194 = t109 * t110 * t94;
	t177 = qJD(1) * t138;
	t93 = -qJD(6) * t165 - t125 * t162 - t126 * t177 + t153 * t183;
	t197 = (-t108 * t194 - t93 * t192) / t101 ^ 2;
	t134 = t138 ^ 2;
	t185 = t133 * t134;
	t100 = t103 * t185 + 0.1e1;
	t96 = 0.1e1 / t100;
	t196 = t103 * t96;
	t175 = qJD(2) * t136;
	t163 = t129 * t175;
	t164 = t137 * t177;
	t95 = ((t135 * t175 - t164) * t128 + t133 * t163) * t121;
	t155 = -t95 + t175;
	t156 = -t136 * t95 + qJD(2);
	t189 = t117 * t137;
	t89 = t156 * t189 + (t155 * t135 - t164) * t116;
	t195 = t102 * t103 * t89;
	t193 = t109 * t125;
	t191 = t114 * t126;
	t190 = t116 * t135;
	t187 = t128 * t133;
	t179 = qJD(1) * t136;
	t178 = qJD(1) * t137;
	t176 = qJD(2) * t135;
	t174 = qJD(2) * t137;
	t151 = t133 * t136 * t177;
	t172 = 0.2e1 * (-t185 * t195 + (-t134 * t135 * t174 - t151) * t103) / t100 ^ 2;
	t171 = 0.2e1 * t197;
	t170 = 0.2e1 * t195;
	t132 = t137 * t133;
	t148 = qJD(2) * (-t129 * t132 - t137) * t128;
	t169 = 0.2e1 * (t129 * t151 + t131 * t148) / t123 ^ 2;
	t168 = t96 * t176;
	t166 = t121 * t187;
	t161 = 0.1e1 + t186;
	t160 = t102 * t172;
	t159 = 0.2e1 * t114 * t194;
	t158 = t136 * t169;
	t157 = t137 * t169;
	t152 = t136 * t166;
	t150 = t161 * t138;
	t149 = t110 * t191 - t193;
	t98 = 0.1e1 / t101;
	t147 = t149 * t98;
	t146 = -t136 * t174 - t153 * t138;
	t113 = -t135 * t182 - t188;
	t112 = t126 * t138 - t135 * t183;
	t107 = t161 * t136 * t121;
	t92 = (t198 * t137 * t116 + t117 * t152) * t138;
	t91 = t136 * t190 + t189 + (-t117 * t181 - t190) * t107;
	t90 = -t161 * t158 + (qJD(1) * t150 + 0.2e1 * t136 * t148) * t121;
	t1 = [t128 * t138 * t157 + (t128 * t136 * t178 + qJD(2) * t150) * t121, t90, 0, 0, 0, 0; (t102 * t168 + (t160 + (qJD(1) * t92 + t89) * t196) * t137) * t136 + (t92 * t137 * t96 * t170 + (t92 * t168 + (t92 * t172 + ((t95 * t152 + t198 * t176 + t157) * t116 + (t158 * t187 + t137 * t95 + (t132 * t163 - (t95 - 0.2e1 * t175) * t137) * t121) * t117) * t96 * t138) * t137) * t103 + (-t102 + ((t131 - t134) * t117 * t166 + t198 * t167) * t103) * t96 * t178) * t138, (t102 * t96 * t179 + (t160 + (qJD(2) * t91 + t89) * t196) * t138) * t135 + (t91 * t138 * t103 * t172 + ((-qJD(2) * t102 + t91 * t170) * t138 + (t91 * t179 + (-(-t107 * t177 - t136 * t90) * t117 - ((t107 * t136 - 0.1e1) * t95 + (-t107 + t136) * qJD(2)) * t116) * t137 * t138) * t103) * t96 - ((-t90 + t177) * t116 + (t155 * t107 - t156) * t117) * t184 * t196) * t137, 0, 0, 0, 0; (-t109 * t112 + t113 * t192) * t171 + (t113 * t159 - t154 * t109 * t182 + t146 * t193 + (-t154 * t114 * t183 - t112 * t94 + t113 * t93 - t146 * t191) * t110) * t98, t135 * t147 * t173 + (t147 * t179 + (t149 * t171 + ((qJD(6) * t109 + t159) * t126 + (t126 * t93 + (qJD(6) * t114 - t94) * t125) * t110) * t98) * t138) * t137, 0, 0, 0, -0.2e1 * t197 + 0.2e1 * (-t110 * t93 * t98 + (-t110 * t197 - t98 * t194) * t114) * t114;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end