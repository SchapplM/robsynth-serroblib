% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR10
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
%   Wie in S6RRPPRR10_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:52
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR10_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR10_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:52:05
	% EndTime: 2019-10-10 09:52:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:52:05
	% EndTime: 2019-10-10 09:52:05
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:52:05
	% EndTime: 2019-10-10 09:52:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:52:05
	% EndTime: 2019-10-10 09:52:06
	% DurationCPUTime: 0.78s
	% Computational Cost: add. (774->70), mult. (1839->155), div. (436->14), fcn. (2165->7), ass. (0->75)
	t87 = cos(qJ(1));
	t116 = qJD(1) * t87;
	t85 = sin(qJ(1));
	t138 = 0.2e1 * t85;
	t74 = t85 ^ 2;
	t75 = 0.1e1 / t85;
	t83 = t87 ^ 2;
	t136 = (0.1e1 + 0.1e1 / t74 * t83) * t75 * t116;
	t84 = sin(qJ(2));
	t118 = t85 * t84;
	t86 = cos(qJ(2));
	t65 = atan2(-t118, -t86);
	t63 = sin(t65);
	t108 = t63 * t118;
	t64 = cos(t65);
	t59 = -t64 * t86 - t108;
	t56 = 0.1e1 / t59;
	t79 = 0.1e1 / t86;
	t57 = 0.1e1 / t59 ^ 2;
	t135 = -0.2e1 * t84;
	t73 = t84 ^ 2;
	t80 = 0.1e1 / t86 ^ 2;
	t123 = t73 * t80;
	t70 = t74 * t123 + 0.1e1;
	t66 = 0.1e1 / t70;
	t134 = t66 - 0.1e1;
	t106 = t84 * t116;
	t115 = qJD(2) * t85;
	t125 = t64 * t84;
	t114 = qJD(2) * t86;
	t52 = (-(-t85 * t114 - t106) * t79 + t115 * t123) * t66;
	t48 = (-t52 * t85 + qJD(2)) * t125 + (-t106 + (t52 - t115) * t86) * t63;
	t133 = t48 * t56 * t57;
	t132 = t52 * t63;
	t131 = t52 * t84;
	t130 = t57 * t84;
	t129 = t57 * t87;
	t120 = t79 * t84;
	t72 = t84 * t73;
	t78 = t86 ^ 2;
	t94 = qJD(2) * (t72 * t79 / t78 + t120);
	t99 = t73 * t85 * t116;
	t128 = (t74 * t94 + t80 * t99) / t70 ^ 2;
	t105 = 0.1e1 + t123;
	t62 = t105 * t85 * t66;
	t127 = t62 * t85;
	t126 = t63 * t86;
	t124 = t73 * t79;
	t122 = t73 * t83;
	t76 = 0.1e1 / t85 ^ 2;
	t121 = t76 * t83;
	t119 = t84 * t87;
	t117 = qJD(1) * t85;
	t113 = qJD(2) * t87;
	t55 = t57 * t122 + 0.1e1;
	t98 = t83 * t84 * t114;
	t112 = 0.2e1 * (-t122 * t133 + (t98 - t99) * t57) / t55 ^ 2;
	t111 = 0.2e1 * t133;
	t71 = t78 * t121 + 0.1e1;
	t110 = 0.2e1 * (-t78 * t136 - t76 * t98) / t71 ^ 2;
	t109 = t57 * t119;
	t107 = t66 * t124;
	t104 = 0.1e1 + t121;
	t103 = t84 * t112;
	t102 = t128 * t135;
	t101 = t128 * t138;
	t100 = t85 * t107;
	t97 = t105 * t87;
	t96 = t104 * t84;
	t68 = 0.1e1 / t71;
	t53 = 0.1e1 / t55;
	t51 = (t134 * t84 * t63 - t64 * t100) * t87;
	t50 = -t85 * t126 + t125 + (-t64 * t118 + t126) * t62;
	t49 = -t105 * t101 + (qJD(1) * t97 + t94 * t138) * t66;
	t1 = [t87 * t79 * t102 + (qJD(2) * t97 - t117 * t120) * t66, t49, 0, 0, 0, 0; (t56 * t103 + (-t56 * t114 + (qJD(1) * t51 + t48) * t130) * t53) * t85 + (t57 * t103 * t51 + (-((t52 * t100 + t134 * t114 + t102) * t63 + (t101 * t124 - t131 + (t131 + (-t72 * t80 + t135) * t115) * t66) * t64) * t109 + (t84 * t111 - t57 * t114) * t51 + (-t56 + ((-t74 + t83) * t64 * t107 + t134 * t108) * t57) * t84 * qJD(1)) * t53) * t87, (t50 * t130 - t56 * t86) * t87 * t112 + ((-t56 * t117 + (-qJD(2) * t50 - t48) * t129) * t86 + (-t56 * t113 - (-t49 * t64 * t85 + t63 * t115 + t127 * t132 - t132 + (-qJD(2) * t63 - t116 * t64) * t62) * t109 + (t87 * t111 + t57 * t117) * t50 - ((t49 - t116) * t63 + ((0.1e1 - t127) * qJD(2) + (t62 - t85) * t52) * t64) * t86 * t129) * t84) * t53, 0, 0, 0, 0; t104 * t86 * t110 + (qJD(2) * t96 + 0.2e1 * t86 * t136) * t68, t75 * t110 * t119 + (-t75 * t86 * t113 + qJD(1) * t96) * t68, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:52:05
	% EndTime: 2019-10-10 09:52:06
	% DurationCPUTime: 0.91s
	% Computational Cost: add. (703->81), mult. (2191->188), div. (456->12), fcn. (2616->9), ass. (0->87)
	t101 = sin(qJ(2));
	t93 = 0.1e1 / t101 ^ 2;
	t103 = cos(qJ(2));
	t97 = t103 ^ 2;
	t151 = t93 * t97;
	t102 = sin(qJ(1));
	t126 = 0.1e1 + t151;
	t95 = t102 ^ 2;
	t91 = t95 * t151 + 0.1e1;
	t89 = 0.1e1 / t91;
	t116 = t126 * t89;
	t75 = t102 * t116;
	t161 = t102 * t75 - 0.1e1;
	t159 = t103 * t151;
	t92 = 0.1e1 / t101;
	t113 = qJD(2) * (-t103 - t159) * t92;
	t104 = cos(qJ(1));
	t138 = qJD(1) * t104;
	t117 = t102 * t97 * t138;
	t160 = (t95 * t113 + t93 * t117) / t91 ^ 2;
	t142 = t102 * t103;
	t88 = atan2(-t142, t101);
	t86 = sin(t88);
	t127 = t86 * t142;
	t87 = cos(t88);
	t70 = t87 * t101 - t127;
	t67 = 0.1e1 / t70;
	t100 = cos(pkin(10));
	t143 = t102 * t100;
	t144 = t101 * t104;
	t99 = sin(pkin(10));
	t83 = t99 * t144 + t143;
	t79 = 0.1e1 / t83;
	t68 = 0.1e1 / t70 ^ 2;
	t80 = 0.1e1 / t83 ^ 2;
	t158 = t89 - 0.1e1;
	t136 = qJD(2) * t102;
	t149 = t101 * t86;
	t63 = ((t101 * t136 - t103 * t138) * t92 + t136 * t151) * t89;
	t58 = (-t63 + t136) * t149 + (-t86 * t138 + (-t102 * t63 + qJD(2)) * t87) * t103;
	t157 = t58 * t67 * t68;
	t98 = t104 ^ 2;
	t150 = t97 * t98;
	t66 = t68 * t150 + 0.1e1;
	t64 = 0.1e1 / t66;
	t155 = t64 * t68;
	t141 = t104 * t100;
	t147 = t102 * t99;
	t114 = t101 * t141 - t147;
	t154 = t80 * t114;
	t153 = t80 * t99;
	t152 = t92 * t97;
	t146 = t104 * t68;
	t145 = qJD(2) * t75;
	t140 = qJD(1) * t102;
	t139 = qJD(1) * t103;
	t137 = qJD(2) * t101;
	t135 = qJD(2) * t103;
	t134 = qJD(2) * t104;
	t133 = 0.2e1 * (-t150 * t157 + (-t101 * t98 * t135 - t117) * t68) / t66 ^ 2;
	t132 = 0.2e1 * t157;
	t78 = t114 ^ 2;
	t73 = t78 * t80 + 0.1e1;
	t122 = t103 * t134;
	t84 = t101 * t143 + t104 * t99;
	t76 = t84 * qJD(1) - t100 * t122;
	t85 = -t101 * t147 + t141;
	t77 = t85 * qJD(1) + t99 * t122;
	t81 = t79 * t80;
	t131 = 0.2e1 * (-t78 * t81 * t77 - t76 * t154) / t73 ^ 2;
	t130 = 0.2e1 * t160;
	t129 = -0.2e1 * t81 * t114;
	t128 = t102 * t89 * t92;
	t125 = t103 * t158;
	t124 = t64 * t137;
	t123 = t102 * t135;
	t121 = t67 * t133;
	t120 = t68 * t133;
	t119 = t103 * t130;
	t118 = t97 * t128;
	t115 = t100 * t79 - t114 * t153;
	t71 = 0.1e1 / t73;
	t112 = t71 * t115;
	t62 = (t87 * t118 + t86 * t125) * t104;
	t60 = -t161 * t87 * t103 + (t102 - t75) * t149;
	t59 = t116 * t138 + 0.2e1 * (t113 * t89 - t126 * t160) * t102;
	t1 = [t128 * t139 + (qJD(2) * t116 + t92 * t119) * t104, t59, 0, 0, 0, 0; (t67 * t124 + (t121 + (qJD(1) * t62 + t58) * t155) * t103) * t102 + (t62 * t68 * t124 + (t62 * t120 + (t62 * t132 + ((t63 * t118 + t158 * t137 + t119) * t86 + (-t63 * t125 + (t130 * t152 + (0.2e1 * t103 + t159) * t89 * qJD(2)) * t102) * t87) * t146) * t64) * t103 + (-t67 + (-(-t95 + t98) * t89 * t87 * t152 + t158 * t127) * t68) * t64 * t139) * t104, (t67 * t64 * t140 + (t121 + (qJD(2) * t60 + t58) * t155) * t104) * t101 + (t60 * t104 * t120 + (-t67 * t134 - ((-t102 * t59 - t138 * t75) * t87 + (t161 * t63 + t136 - t145) * t86) * t103 * t146 + (t104 * t132 + t68 * t140) * t60) * t64 - ((-t59 + t138) * t86 + (-t63 * t75 - qJD(2) + (t63 + t145) * t102) * t87) * t144 * t155) * t103, 0, 0, 0, 0; (-t85 * t154 - t79 * t84) * t131 + ((t114 * qJD(1) + t100 * t123) * t79 + t85 * t77 * t129 + (-t84 * t77 + (-t83 * qJD(1) - t99 * t123) * t114 - t85 * t76) * t80) * t71, t101 * t112 * t134 + (t112 * t140 + (t115 * t131 + (-t76 * t153 + (t100 * t80 + t99 * t129) * t77) * t71) * t104) * t103, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:52:05
	% EndTime: 2019-10-10 09:52:06
	% DurationCPUTime: 0.95s
	% Computational Cost: add. (1161->91), mult. (2519->208), div. (480->12), fcn. (2968->9), ass. (0->92)
	t134 = sin(qJ(1));
	t135 = cos(qJ(2));
	t133 = sin(qJ(2));
	t152 = qJD(1) * t133 + qJD(5);
	t136 = cos(qJ(1));
	t172 = qJD(2) * t136;
	t197 = t152 * t134 - t135 * t172;
	t173 = qJD(2) * t135;
	t196 = t134 * t173 + t152 * t136;
	t180 = t134 * t135;
	t118 = atan2(-t180, t133);
	t117 = cos(t118);
	t116 = sin(t118);
	t166 = t116 * t180;
	t105 = t117 * t133 - t166;
	t102 = 0.1e1 / t105;
	t125 = pkin(10) + qJ(5);
	t123 = sin(t125);
	t124 = cos(t125);
	t181 = t134 * t124;
	t182 = t133 * t136;
	t113 = t123 * t182 + t181;
	t109 = 0.1e1 / t113;
	t126 = 0.1e1 / t133;
	t103 = 0.1e1 / t105 ^ 2;
	t110 = 0.1e1 / t113 ^ 2;
	t127 = 0.1e1 / t133 ^ 2;
	t129 = t134 ^ 2;
	t131 = t135 ^ 2;
	t185 = t127 * t131;
	t122 = t129 * t185 + 0.1e1;
	t120 = 0.1e1 / t122;
	t195 = t120 - 0.1e1;
	t112 = t123 * t134 - t124 * t182;
	t108 = t112 ^ 2;
	t101 = t108 * t110 + 0.1e1;
	t189 = t110 * t112;
	t153 = qJD(5) * t133 + qJD(1);
	t147 = t153 * t136;
	t94 = -t197 * t123 + t124 * t147;
	t191 = t109 * t110 * t94;
	t93 = t123 * t147 + t197 * t124;
	t194 = (-t108 * t191 + t93 * t189) / t101 ^ 2;
	t132 = t136 ^ 2;
	t184 = t131 * t132;
	t100 = t103 * t184 + 0.1e1;
	t96 = 0.1e1 / t100;
	t193 = t103 * t96;
	t174 = qJD(2) * t134;
	t163 = t127 * t174;
	t176 = qJD(1) * t136;
	t164 = t135 * t176;
	t95 = ((t133 * t174 - t164) * t126 + t131 * t163) * t120;
	t154 = -t95 + t174;
	t155 = -t134 * t95 + qJD(2);
	t187 = t117 * t135;
	t89 = t155 * t187 + (t154 * t133 - t164) * t116;
	t192 = t102 * t103 * t89;
	t190 = t109 * t124;
	t188 = t112 * t123;
	t186 = t126 * t131;
	t183 = t133 * t134;
	t178 = qJD(1) * t134;
	t177 = qJD(1) * t135;
	t175 = qJD(2) * t133;
	t150 = t131 * t134 * t176;
	t171 = 0.2e1 * (-t184 * t192 + (-t132 * t133 * t173 - t150) * t103) / t100 ^ 2;
	t170 = 0.2e1 * t194;
	t169 = 0.2e1 * t192;
	t130 = t135 * t131;
	t145 = qJD(2) * (-t127 * t130 - t135) * t126;
	t168 = 0.2e1 * (t127 * t150 + t129 * t145) / t122 ^ 2;
	t167 = t96 * t175;
	t165 = t120 * t186;
	t160 = 0.1e1 + t185;
	t159 = t102 * t171;
	t158 = 0.2e1 * t112 * t191;
	t157 = t134 * t168;
	t156 = t135 * t168;
	t151 = t134 * t165;
	t149 = t160 * t136;
	t148 = t134 * t153;
	t146 = t110 * t188 + t190;
	t98 = 0.1e1 / t101;
	t144 = t146 * t98;
	t115 = -t123 * t183 + t124 * t136;
	t114 = t123 * t136 + t133 * t181;
	t107 = t160 * t134 * t120;
	t92 = (t195 * t135 * t116 + t117 * t151) * t136;
	t91 = t116 * t183 + t187 + (-t116 * t133 - t117 * t180) * t107;
	t90 = -t160 * t157 + (qJD(1) * t149 + 0.2e1 * t134 * t145) * t120;
	t1 = [t126 * t136 * t156 + (t126 * t134 * t177 + qJD(2) * t149) * t120, t90, 0, 0, 0, 0; (t102 * t167 + (t159 + (qJD(1) * t92 + t89) * t193) * t135) * t134 + (t92 * t135 * t96 * t169 + (t92 * t167 + (t92 * t171 + ((t95 * t151 + t195 * t175 + t156) * t116 + (t157 * t186 + t135 * t95 + (t130 * t163 - (t95 - 0.2e1 * t174) * t135) * t120) * t117) * t96 * t136) * t135) * t103 + (-t102 + ((t129 - t132) * t117 * t165 + t195 * t166) * t103) * t96 * t177) * t136, (t102 * t96 * t178 + (t159 + (qJD(2) * t91 + t89) * t193) * t136) * t133 + (t91 * t136 * t103 * t171 + ((-qJD(2) * t102 + t91 * t169) * t136 + (t91 * t178 + (-(-t107 * t176 - t134 * t90) * t117 - ((t107 * t134 - 0.1e1) * t95 + (-t107 + t134) * qJD(2)) * t116) * t135 * t136) * t103) * t96 - ((-t90 + t176) * t116 + (t154 * t107 - t155) * t117) * t182 * t193) * t135, 0, 0, 0, 0; (-t109 * t114 + t115 * t189) * t170 + (t115 * t158 - t109 * t123 * t148 + t196 * t190 + (t112 * t124 * t148 - t114 * t94 - t115 * t93 + t196 * t188) * t110) * t98, t133 * t144 * t172 + (t144 * t178 + (t146 * t170 + ((qJD(5) * t109 + t158) * t123 + (-t123 * t93 + (-qJD(5) * t112 + t94) * t124) * t110) * t98) * t136) * t135, 0, 0, -0.2e1 * t194 + 0.2e1 * (t110 * t93 * t98 + (-t110 * t194 - t98 * t191) * t112) * t112, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:52:06
	% EndTime: 2019-10-10 09:52:07
	% DurationCPUTime: 1.06s
	% Computational Cost: add. (1820->94), mult. (2734->208), div. (498->12), fcn. (3199->9), ass. (0->96)
	t154 = sin(qJ(2));
	t148 = 0.1e1 / t154 ^ 2;
	t156 = cos(qJ(2));
	t152 = t156 ^ 2;
	t202 = t148 * t152;
	t220 = t156 * t202;
	t155 = sin(qJ(1));
	t177 = 0.1e1 + t202;
	t219 = t155 * t177;
	t146 = qJD(5) + qJD(6);
	t173 = qJD(1) * t154 + t146;
	t157 = cos(qJ(1));
	t189 = qJD(2) * t157;
	t218 = t155 * t173 - t156 * t189;
	t190 = qJD(2) * t156;
	t217 = t155 * t190 + t157 * t173;
	t196 = t155 * t156;
	t139 = atan2(-t196, t154);
	t137 = cos(t139);
	t136 = sin(t139);
	t182 = t136 * t196;
	t125 = t137 * t154 - t182;
	t122 = 0.1e1 / t125;
	t145 = pkin(10) + qJ(5) + qJ(6);
	t143 = sin(t145);
	t144 = cos(t145);
	t197 = t155 * t144;
	t198 = t154 * t157;
	t133 = t143 * t198 + t197;
	t129 = 0.1e1 / t133;
	t147 = 0.1e1 / t154;
	t123 = 0.1e1 / t125 ^ 2;
	t130 = 0.1e1 / t133 ^ 2;
	t150 = t155 ^ 2;
	t142 = t150 * t202 + 0.1e1;
	t140 = 0.1e1 / t142;
	t216 = t140 - 0.1e1;
	t174 = t146 * t154 + qJD(1);
	t168 = t174 * t157;
	t113 = t143 * t168 + t144 * t218;
	t132 = t143 * t155 - t144 * t198;
	t128 = t132 ^ 2;
	t118 = t128 * t130 + 0.1e1;
	t205 = t130 * t132;
	t114 = -t143 * t218 + t144 * t168;
	t213 = t114 * t129 * t130;
	t215 = (t113 * t205 - t128 * t213) / t118 ^ 2;
	t193 = qJD(1) * t157;
	t180 = t156 * t193;
	t191 = qJD(2) * t155;
	t115 = ((t154 * t191 - t180) * t147 + t191 * t202) * t140;
	t203 = t137 * t156;
	t109 = (-t115 * t155 + qJD(2)) * t203 + (-t180 + (-t115 + t191) * t154) * t136;
	t214 = t109 * t122 * t123;
	t212 = t115 * t136;
	t211 = t115 * t156;
	t210 = t123 * t156;
	t209 = t123 * t157;
	t166 = qJD(2) * (-t156 - t220) * t147;
	t200 = t152 * t155;
	t171 = t193 * t200;
	t208 = (t148 * t171 + t150 * t166) / t142 ^ 2;
	t127 = t140 * t219;
	t207 = t127 * t155;
	t206 = t129 * t144;
	t204 = t132 * t143;
	t153 = t157 ^ 2;
	t201 = t152 * t153;
	t199 = t154 * t155;
	t195 = qJD(1) * t155;
	t194 = qJD(1) * t156;
	t192 = qJD(2) * t154;
	t121 = t123 * t201 + 0.1e1;
	t188 = 0.2e1 * (-t201 * t214 + (-t153 * t154 * t190 - t171) * t123) / t121 ^ 2;
	t187 = 0.2e1 * t215;
	t186 = 0.2e1 * t214;
	t185 = -0.2e1 * t208;
	t184 = t156 * t209;
	t183 = t156 * t208;
	t181 = t147 * t200;
	t176 = t156 * t188;
	t175 = 0.2e1 * t132 * t213;
	t172 = t140 * t181;
	t170 = t177 * t157;
	t169 = t174 * t155;
	t167 = t130 * t204 + t206;
	t165 = t167 * t157;
	t135 = -t143 * t199 + t144 * t157;
	t134 = t143 * t157 + t154 * t197;
	t119 = 0.1e1 / t121;
	t116 = 0.1e1 / t118;
	t112 = (t136 * t156 * t216 + t137 * t172) * t157;
	t111 = t136 * t199 + t203 + (-t136 * t154 - t137 * t196) * t127;
	t110 = t185 * t219 + (qJD(1) * t170 + 0.2e1 * t155 * t166) * t140;
	t106 = -0.2e1 * t215 + 0.2e1 * (t113 * t116 * t130 + (-t116 * t213 - t130 * t215) * t132) * t132;
	t1 = [0.2e1 * t147 * t157 * t183 + (t147 * t155 * t194 + qJD(2) * t170) * t140, t110, 0, 0, 0, 0; (t122 * t176 + (t122 * t192 + (qJD(1) * t112 + t109) * t210) * t119) * t155 + (t123 * t176 * t112 + (-((-t115 * t172 - t192 * t216 - 0.2e1 * t183) * t136 + (t181 * t185 - t211 + (t211 + (-0.2e1 * t156 - t220) * t191) * t140) * t137) * t184 + (t123 * t192 + t156 * t186) * t112 + (-t122 + ((t150 - t153) * t152 * t147 * t140 * t137 + t216 * t182) * t123) * t194) * t119) * t157, (t111 * t210 + t122 * t154) * t157 * t188 + ((t122 * t195 + (qJD(2) * t111 + t109) * t209) * t154 + (-t122 * t189 - (-t110 * t137 * t155 + t136 * t191 + t207 * t212 - t212 + (-qJD(2) * t136 - t137 * t193) * t127) * t184 + (t123 * t195 + t157 * t186) * t111 - ((-t110 + t193) * t136 + ((-0.1e1 + t207) * qJD(2) + (-t127 + t155) * t115) * t137) * t123 * t198) * t156) * t119, 0, 0, 0, 0; (-t129 * t134 + t135 * t205) * t187 + (t135 * t175 - t129 * t143 * t169 + t217 * t206 + (t132 * t144 * t169 - t135 * t113 - t134 * t114 + t204 * t217) * t130) * t116, t156 * t165 * t187 + (t165 * t192 + (t167 * t195 + ((t129 * t146 + t175) * t143 + (-t113 * t143 + (-t132 * t146 + t114) * t144) * t130) * t157) * t156) * t116, 0, 0, t106, t106;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end