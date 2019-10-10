% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR11
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
%   Wie in S6RRPRPR11_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:22
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR11_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR11_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:22:42
	% EndTime: 2019-10-10 10:22:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:22:42
	% EndTime: 2019-10-10 10:22:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:22:43
	% EndTime: 2019-10-10 10:22:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:22:43
	% EndTime: 2019-10-10 10:22:43
	% DurationCPUTime: 0.79s
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
	% StartTime: 2019-10-10 10:22:43
	% EndTime: 2019-10-10 10:22:44
	% DurationCPUTime: 0.95s
	% Computational Cost: add. (813->90), mult. (2519->205), div. (480->12), fcn. (2968->9), ass. (0->93)
	t131 = cos(qJ(2));
	t129 = sin(qJ(1));
	t123 = t129 ^ 2;
	t128 = sin(qJ(2));
	t121 = 0.1e1 / t128 ^ 2;
	t125 = t131 ^ 2;
	t179 = t121 * t125;
	t118 = t123 * t179 + 0.1e1;
	t116 = 0.1e1 / t118;
	t120 = 0.1e1 / t128;
	t169 = qJD(2) * t129;
	t157 = t121 * t169;
	t132 = cos(qJ(1));
	t171 = qJD(1) * t132;
	t158 = t131 * t171;
	t90 = ((t128 * t169 - t158) * t120 + t125 * t157) * t116;
	t193 = t131 * t90;
	t148 = -t90 + t169;
	t146 = qJD(1) * t128 + qJD(4);
	t167 = qJD(2) * t132;
	t192 = t146 * t129 - t131 * t167;
	t176 = t129 * t131;
	t115 = atan2(-t176, t128);
	t114 = cos(t115);
	t113 = sin(t115);
	t161 = t113 * t176;
	t99 = t114 * t128 - t161;
	t96 = 0.1e1 / t99;
	t127 = sin(qJ(4));
	t175 = t132 * t127;
	t130 = cos(qJ(4));
	t177 = t129 * t130;
	t110 = t128 * t175 + t177;
	t106 = 0.1e1 / t110;
	t107 = 0.1e1 / t110 ^ 2;
	t97 = 0.1e1 / t99 ^ 2;
	t191 = t116 - 0.1e1;
	t149 = -t129 * t90 + qJD(2);
	t181 = t114 * t131;
	t85 = t149 * t181 + (t148 * t128 - t158) * t113;
	t190 = t85 * t96 * t97;
	t126 = t132 ^ 2;
	t95 = t126 * t125 * t97 + 0.1e1;
	t93 = 0.1e1 / t95;
	t189 = t93 * t97;
	t188 = t96 * t93;
	t174 = t132 * t130;
	t178 = t129 * t127;
	t109 = -t128 * t174 + t178;
	t105 = t109 ^ 2;
	t104 = t105 * t107 + 0.1e1;
	t184 = t107 * t109;
	t147 = qJD(4) * t128 + qJD(1);
	t143 = t147 * t132;
	t92 = -t192 * t127 + t130 * t143;
	t186 = t106 * t107 * t92;
	t91 = t127 * t143 + t192 * t130;
	t187 = 0.1e1 / t104 ^ 2 * (-t105 * t186 + t91 * t184);
	t185 = t132 * t97;
	t183 = t109 * t127;
	t182 = t113 * t129;
	t180 = t120 * t125;
	t173 = qJD(1) * t129;
	t172 = qJD(1) * t131;
	t170 = qJD(2) * t128;
	t168 = qJD(2) * t131;
	t159 = t129 * t171;
	t162 = t97 * t170;
	t166 = 0.2e1 * (-t126 * t131 * t162 + (-t126 * t190 - t97 * t159) * t125) / t95 ^ 2;
	t165 = 0.2e1 * t190;
	t164 = 0.2e1 * t187;
	t124 = t131 * t125;
	t141 = qJD(2) * (-t121 * t124 - t131) * t120;
	t163 = 0.2e1 * (t123 * t141 + t159 * t179) / t118 ^ 2;
	t160 = t116 * t180;
	t155 = t96 * t166;
	t154 = t97 * t166;
	t153 = 0.1e1 + t179;
	t152 = 0.2e1 * t109 * t186;
	t151 = t129 * t163;
	t150 = t131 * t163;
	t145 = t129 * t160;
	t144 = t153 * t132;
	t142 = t106 * t130 + t107 * t183;
	t140 = t142 * t132;
	t112 = -t128 * t178 + t174;
	t111 = t128 * t177 + t175;
	t103 = t153 * t129 * t116;
	t101 = 0.1e1 / t104;
	t89 = (t191 * t131 * t113 + t114 * t145) * t132;
	t88 = t128 * t182 + t181 + (-t113 * t128 - t114 * t176) * t103;
	t86 = -t153 * t151 + (qJD(1) * t144 + 0.2e1 * t129 * t141) * t116;
	t1 = [t132 * t120 * t150 + (t120 * t129 * t172 + qJD(2) * t144) * t116, t86, 0, 0, 0, 0; (t170 * t188 + (t155 + (qJD(1) * t89 + t85) * t189) * t131) * t129 + (t89 * t154 * t131 + (t89 * t162 + (t89 * t165 + ((t90 * t145 + t191 * t170 + t150) * t113 + (t151 * t180 + t193 + (t124 * t157 - (t90 - 0.2e1 * t169) * t131) * t116) * t114) * t185) * t131 + (-t96 + (-(-t123 + t126) * t114 * t160 + t191 * t161) * t97) * t172) * t93) * t132, (t173 * t188 + (t155 + (qJD(2) * t88 + t85) * t189) * t132) * t128 + (t88 * t132 * t154 + (-t96 * t167 + (t132 * t165 + t97 * t173) * t88 + (-t103 * t182 * t193 + (-(-t103 * t171 - t129 * t86) * t131 - (t148 * t103 - t149) * t128) * t114 + (-(-qJD(2) * t103 + t148) * t131 - (-t86 + t171) * t128) * t113) * t185) * t93) * t131, 0, 0, 0, 0; (-t106 * t111 + t112 * t184) * t164 + (t112 * t152 + (-t111 * t92 - t112 * t91 + t147 * t109 * t177 - (-t129 * t168 - t146 * t132) * t183) * t107 + (t146 * t174 + (-t147 * t127 + t130 * t168) * t129) * t106) * t101, t131 * t140 * t164 + (t140 * t170 + (t142 * t173 + ((qJD(4) * t106 + t152) * t127 + (-t127 * t91 + (-qJD(4) * t109 + t92) * t130) * t107) * t132) * t131) * t101, 0, -0.2e1 * t187 + 0.2e1 * (t91 * t107 * t101 + (-t101 * t186 - t107 * t187) * t109) * t109, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:22:43
	% EndTime: 2019-10-10 10:22:44
	% DurationCPUTime: 0.96s
	% Computational Cost: add. (1161->91), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->91)
	t136 = sin(qJ(1));
	t138 = cos(qJ(1));
	t135 = sin(qJ(2));
	t154 = qJD(1) * t135 + qJD(4);
	t137 = cos(qJ(2));
	t174 = qJD(2) * t137;
	t198 = t154 * t136 - t138 * t174;
	t197 = t136 * t174 + t154 * t138;
	t181 = t136 * t137;
	t120 = atan2(-t181, t135);
	t119 = cos(t120);
	t118 = sin(t120);
	t168 = t118 * t181;
	t107 = t119 * t135 - t168;
	t104 = 0.1e1 / t107;
	t127 = qJ(4) + pkin(10);
	t125 = sin(t127);
	t126 = cos(t127);
	t182 = t136 * t126;
	t183 = t135 * t138;
	t115 = t125 * t183 + t182;
	t111 = 0.1e1 / t115;
	t128 = 0.1e1 / t135;
	t105 = 0.1e1 / t107 ^ 2;
	t112 = 0.1e1 / t115 ^ 2;
	t129 = 0.1e1 / t135 ^ 2;
	t131 = t136 ^ 2;
	t133 = t137 ^ 2;
	t186 = t129 * t133;
	t124 = t131 * t186 + 0.1e1;
	t122 = 0.1e1 / t124;
	t196 = t122 - 0.1e1;
	t114 = t125 * t136 - t126 * t183;
	t110 = t114 ^ 2;
	t103 = t110 * t112 + 0.1e1;
	t190 = t112 * t114;
	t155 = qJD(4) * t135 + qJD(1);
	t149 = t155 * t138;
	t96 = -t198 * t125 + t126 * t149;
	t192 = t111 * t112 * t96;
	t95 = t125 * t149 + t198 * t126;
	t195 = 0.1e1 / t103 ^ 2 * (-t110 * t192 + t95 * t190);
	t134 = t138 ^ 2;
	t185 = t133 * t134;
	t102 = t105 * t185 + 0.1e1;
	t98 = 0.1e1 / t102;
	t194 = t105 * t98;
	t175 = qJD(2) * t136;
	t165 = t129 * t175;
	t177 = qJD(1) * t138;
	t166 = t137 * t177;
	t97 = ((t135 * t175 - t166) * t128 + t133 * t165) * t122;
	t156 = -t97 + t175;
	t157 = -t136 * t97 + qJD(2);
	t188 = t119 * t137;
	t91 = t157 * t188 + (t156 * t135 - t166) * t118;
	t193 = t104 * t105 * t91;
	t191 = t111 * t126;
	t189 = t114 * t125;
	t187 = t128 * t133;
	t184 = t135 * t136;
	t179 = qJD(1) * t136;
	t178 = qJD(1) * t137;
	t176 = qJD(2) * t135;
	t152 = t133 * t136 * t177;
	t173 = 0.2e1 * (-t185 * t193 + (-t134 * t135 * t174 - t152) * t105) / t102 ^ 2;
	t172 = 0.2e1 * t195;
	t171 = 0.2e1 * t193;
	t132 = t137 * t133;
	t147 = qJD(2) * (-t129 * t132 - t137) * t128;
	t170 = 0.2e1 * (t129 * t152 + t131 * t147) / t124 ^ 2;
	t169 = t98 * t176;
	t167 = t122 * t187;
	t162 = 0.1e1 + t186;
	t161 = t104 * t173;
	t160 = 0.2e1 * t114 * t192;
	t159 = t136 * t170;
	t158 = t137 * t170;
	t153 = t136 * t167;
	t151 = t162 * t138;
	t150 = t155 * t136;
	t148 = t112 * t189 + t191;
	t146 = t148 * t138;
	t117 = -t125 * t184 + t126 * t138;
	t116 = t125 * t138 + t135 * t182;
	t109 = t162 * t136 * t122;
	t100 = 0.1e1 / t103;
	t94 = (t196 * t137 * t118 + t119 * t153) * t138;
	t93 = t118 * t184 + t188 + (-t118 * t135 - t119 * t181) * t109;
	t92 = -t162 * t159 + (qJD(1) * t151 + 0.2e1 * t136 * t147) * t122;
	t1 = [t128 * t138 * t158 + (t128 * t136 * t178 + qJD(2) * t151) * t122, t92, 0, 0, 0, 0; (t104 * t169 + (t161 + (qJD(1) * t94 + t91) * t194) * t137) * t136 + (t94 * t137 * t98 * t171 + (t94 * t169 + (t94 * t173 + ((t97 * t153 + t196 * t176 + t158) * t118 + (t159 * t187 + t137 * t97 + (t132 * t165 - (t97 - 0.2e1 * t175) * t137) * t122) * t119) * t98 * t138) * t137) * t105 + (-t104 + ((t131 - t134) * t119 * t167 + t196 * t168) * t105) * t98 * t178) * t138, (t104 * t98 * t179 + (t161 + (qJD(2) * t93 + t91) * t194) * t138) * t135 + (t93 * t138 * t105 * t173 + ((-qJD(2) * t104 + t93 * t171) * t138 + (t93 * t179 + (-(-t109 * t177 - t136 * t92) * t119 - ((t109 * t136 - 0.1e1) * t97 + (-t109 + t136) * qJD(2)) * t118) * t137 * t138) * t105) * t98 - ((-t92 + t177) * t118 + (t156 * t109 - t157) * t119) * t183 * t194) * t137, 0, 0, 0, 0; (-t111 * t116 + t117 * t190) * t172 + (t117 * t160 - t111 * t125 * t150 + t197 * t191 + (t114 * t126 * t150 - t116 * t96 - t117 * t95 + t197 * t189) * t112) * t100, t137 * t146 * t172 + (t146 * t176 + (t148 * t179 + ((qJD(4) * t111 + t160) * t125 + (-t125 * t95 + (-qJD(4) * t114 + t96) * t126) * t112) * t138) * t137) * t100, 0, -0.2e1 * t195 + 0.2e1 * (t100 * t112 * t95 + (-t100 * t192 - t112 * t195) * t114) * t114, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:22:43
	% EndTime: 2019-10-10 10:22:44
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
	t146 = qJD(4) + qJD(6);
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
	t145 = qJ(4) + pkin(10) + qJ(6);
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
	t1 = [0.2e1 * t147 * t157 * t183 + (t147 * t155 * t194 + qJD(2) * t170) * t140, t110, 0, 0, 0, 0; (t122 * t176 + (t122 * t192 + (qJD(1) * t112 + t109) * t210) * t119) * t155 + (t123 * t176 * t112 + (-((-t115 * t172 - t192 * t216 - 0.2e1 * t183) * t136 + (t181 * t185 - t211 + (t211 + (-0.2e1 * t156 - t220) * t191) * t140) * t137) * t184 + (t123 * t192 + t156 * t186) * t112 + (-t122 + ((t150 - t153) * t152 * t147 * t140 * t137 + t216 * t182) * t123) * t194) * t119) * t157, (t111 * t210 + t122 * t154) * t157 * t188 + ((t122 * t195 + (qJD(2) * t111 + t109) * t209) * t154 + (-t122 * t189 - (-t110 * t137 * t155 + t136 * t191 + t207 * t212 - t212 + (-qJD(2) * t136 - t137 * t193) * t127) * t184 + (t123 * t195 + t157 * t186) * t111 - ((-t110 + t193) * t136 + ((-0.1e1 + t207) * qJD(2) + (-t127 + t155) * t115) * t137) * t123 * t198) * t156) * t119, 0, 0, 0, 0; (-t129 * t134 + t135 * t205) * t187 + (t135 * t175 - t129 * t143 * t169 + t217 * t206 + (t132 * t144 * t169 - t135 * t113 - t134 * t114 + t204 * t217) * t130) * t116, t156 * t165 * t187 + (t165 * t192 + (t167 * t195 + ((t129 * t146 + t175) * t143 + (-t113 * t143 + (-t132 * t146 + t114) * t144) * t130) * t157) * t156) * t116, 0, t106, 0, t106;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end