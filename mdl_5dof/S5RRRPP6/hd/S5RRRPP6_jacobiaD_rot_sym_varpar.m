% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRRPP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRPP6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:03:04
	% EndTime: 2019-12-31 21:03:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:03:05
	% EndTime: 2019-12-31 21:03:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:03:05
	% EndTime: 2019-12-31 21:03:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:03:05
	% EndTime: 2019-12-31 21:03:06
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (1002->94), mult. (2519->211), div. (480->12), fcn. (2968->9), ass. (0->92)
	t126 = sin(qJ(1));
	t119 = t126 ^ 2;
	t125 = sin(qJ(2));
	t118 = t125 ^ 2;
	t128 = cos(qJ(2));
	t121 = 0.1e1 / t128 ^ 2;
	t173 = t118 * t121;
	t113 = t119 * t173 + 0.1e1;
	t117 = t125 * t118;
	t120 = 0.1e1 / t128;
	t172 = t120 * t125;
	t137 = qJD(2) * (t117 * t120 * t121 + t172);
	t129 = cos(qJ(1));
	t163 = qJD(1) * t129;
	t151 = t126 * t163;
	t181 = 0.1e1 / t113 ^ 2 * (t119 * t137 + t151 * t173);
	t193 = -0.2e1 * t181;
	t111 = 0.1e1 / t113;
	t146 = 0.1e1 + t173;
	t190 = t126 * t146;
	t98 = t111 * t190;
	t192 = t126 * t98 - 0.1e1;
	t124 = sin(qJ(3));
	t127 = cos(qJ(3));
	t165 = t129 * t127;
	t107 = t126 * t124 + t128 * t165;
	t102 = 0.1e1 / t107 ^ 2;
	t166 = t129 * t124;
	t168 = t126 * t127;
	t106 = t128 * t166 - t168;
	t175 = t106 * t127;
	t101 = 0.1e1 / t107;
	t177 = t101 * t124;
	t139 = t102 * t175 - t177;
	t100 = t106 ^ 2;
	t99 = t100 * t102 + 0.1e1;
	t96 = 0.1e1 / t99;
	t191 = t139 * t96;
	t169 = t126 * t125;
	t110 = atan2(-t169, -t128);
	t109 = cos(t110);
	t108 = sin(t110);
	t154 = t108 * t169;
	t94 = -t109 * t128 - t154;
	t91 = 0.1e1 / t94;
	t92 = 0.1e1 / t94 ^ 2;
	t189 = t111 - 0.1e1;
	t123 = t129 ^ 2;
	t161 = qJD(2) * t128;
	t155 = t92 * t161;
	t150 = t125 * t163;
	t162 = qJD(2) * t126;
	t174 = t109 * t125;
	t149 = t121 * t162;
	t85 = (-(-t126 * t161 - t150) * t120 + t118 * t149) * t111;
	t80 = (-t126 * t85 + qJD(2)) * t174 + (-t150 + (t85 - t162) * t128) * t108;
	t187 = t80 * t91 * t92;
	t90 = t123 * t118 * t92 + 0.1e1;
	t188 = (t123 * t125 * t155 + (-t123 * t187 - t92 * t151) * t118) / t90 ^ 2;
	t176 = t102 * t106;
	t143 = -qJD(1) * t128 + qJD(3);
	t144 = qJD(3) * t128 - qJD(1);
	t160 = qJD(2) * t129;
	t148 = t125 * t160;
	t87 = -t144 * t166 + (t143 * t126 - t148) * t127;
	t183 = t101 * t102 * t87;
	t167 = t126 * t128;
	t138 = t124 * t167 + t165;
	t86 = t138 * qJD(1) - t107 * qJD(3) + t124 * t148;
	t186 = (-t100 * t183 - t86 * t176) / t99 ^ 2;
	t88 = 0.1e1 / t90;
	t185 = t88 * t91;
	t184 = t88 * t92;
	t179 = t129 * t92;
	t178 = qJD(2) * t98;
	t171 = t125 * t129;
	t164 = qJD(1) * t126;
	t159 = 0.2e1 * t187;
	t158 = -0.2e1 * t186;
	t157 = t91 * t188;
	t156 = t106 * t183;
	t153 = t111 * t118 * t120;
	t147 = 0.2e1 * t92 * t188;
	t145 = t120 * t193;
	t142 = t126 * t153;
	t141 = t146 * t129;
	t140 = t143 * t129;
	t105 = -t127 * t167 + t166;
	t84 = (t189 * t125 * t108 - t109 * t142) * t129;
	t83 = -t192 * t174 + (-t126 + t98) * t128 * t108;
	t81 = t190 * t193 + (qJD(1) * t141 + 0.2e1 * t126 * t137) * t111;
	t1 = [t145 * t171 + (qJD(2) * t141 - t164 * t172) * t111, t81, 0, 0, 0; (-t161 * t185 + (0.2e1 * t157 + (qJD(1) * t84 + t80) * t184) * t125) * t126 + (t84 * t147 * t125 + (-t84 * t155 + (t84 * t159 + ((0.2e1 * t125 * t181 - t85 * t142 - t189 * t161) * t108 + (t118 * t126 * t145 + t125 * t85 + (t117 * t149 - (t85 - 0.2e1 * t162) * t125) * t111) * t109) * t179) * t125 + (-t91 + (-(t119 - t123) * t109 * t153 + t189 * t154) * t92) * t125 * qJD(1)) * t88) * t129, (-t164 * t185 + (-0.2e1 * t157 + (-qJD(2) * t83 - t80) * t184) * t129) * t128 + (t83 * t129 * t147 + (-t91 * t160 - ((-t126 * t81 - t163 * t98) * t109 + (t192 * t85 + t162 - t178) * t108) * t92 * t171 + (t129 * t159 + t92 * t164) * t83 - ((t81 - t163) * t108 + (t85 * t98 + qJD(2) + (-t85 - t178) * t126) * t109) * t128 * t179) * t88) * t125, 0, 0, 0; 0.2e1 * (t101 * t138 + t105 * t176) * t186 + (0.2e1 * t105 * t156 - t144 * t101 * t168 + (t125 * t162 + t140) * t177 + (t138 * t87 + t105 * t86 - t140 * t175 - (qJD(2) * t125 * t127 + t144 * t124) * t106 * t126) * t102) * t96, t128 * t160 * t191 + (-t164 * t191 + (t139 * t158 + ((-qJD(3) * t101 - 0.2e1 * t156) * t127 + (-t127 * t86 + (-qJD(3) * t106 + t87) * t124) * t102) * t96) * t129) * t125, t158 + 0.2e1 * (-t86 * t102 * t96 + (-t102 * t186 - t96 * t183) * t106) * t106, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:03:05
	% EndTime: 2019-12-31 21:03:06
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (1350->93), mult. (2519->209), div. (480->12), fcn. (2968->9), ass. (0->94)
	t135 = sin(qJ(1));
	t129 = t135 ^ 2;
	t134 = sin(qJ(2));
	t128 = t134 ^ 2;
	t136 = cos(qJ(2));
	t131 = 0.1e1 / t136 ^ 2;
	t184 = t128 * t131;
	t121 = t129 * t184 + 0.1e1;
	t127 = t134 * t128;
	t130 = 0.1e1 / t136;
	t181 = t130 * t134;
	t146 = qJD(2) * (t127 * t130 * t131 + t181);
	t137 = cos(qJ(1));
	t173 = qJD(1) * t137;
	t182 = t128 * t135;
	t150 = t173 * t182;
	t190 = (t129 * t146 + t131 * t150) / t121 ^ 2;
	t200 = -0.2e1 * t190;
	t126 = qJ(3) + pkin(8);
	t125 = cos(t126);
	t175 = t136 * t137;
	t124 = sin(t126);
	t179 = t135 * t124;
	t114 = t125 * t175 + t179;
	t109 = 0.1e1 / t114 ^ 2;
	t178 = t135 * t125;
	t113 = t124 * t175 - t178;
	t187 = t113 * t125;
	t108 = 0.1e1 / t114;
	t189 = t108 * t124;
	t148 = t109 * t187 - t189;
	t107 = t113 ^ 2;
	t100 = t107 * t109 + 0.1e1;
	t98 = 0.1e1 / t100;
	t199 = t148 * t98;
	t157 = 0.1e1 + t184;
	t198 = t135 * t157;
	t177 = t135 * t134;
	t118 = atan2(-t177, -t136);
	t116 = cos(t118);
	t115 = sin(t118);
	t163 = t115 * t177;
	t104 = -t116 * t136 - t163;
	t101 = 0.1e1 / t104;
	t102 = 0.1e1 / t104 ^ 2;
	t119 = 0.1e1 / t121;
	t197 = t119 - 0.1e1;
	t133 = t137 ^ 2;
	t171 = qJD(2) * t136;
	t183 = t128 * t133;
	t172 = qJD(2) * t135;
	t159 = t131 * t172;
	t160 = t134 * t173;
	t94 = (-(-t135 * t171 - t160) * t130 + t128 * t159) * t119;
	t154 = t94 - t172;
	t155 = -t135 * t94 + qJD(2);
	t186 = t116 * t134;
	t88 = t155 * t186 + (t154 * t136 - t160) * t115;
	t193 = t101 * t102 * t88;
	t97 = t102 * t183 + 0.1e1;
	t196 = (-t183 * t193 + (t133 * t134 * t171 - t150) * t102) / t97 ^ 2;
	t188 = t109 * t113;
	t152 = -qJD(1) * t136 + qJD(3);
	t153 = qJD(3) * t136 - qJD(1);
	t170 = qJD(2) * t137;
	t158 = t134 * t170;
	t185 = t124 * t137;
	t93 = -t153 * t185 + (t152 * t135 - t158) * t125;
	t192 = t108 * t109 * t93;
	t176 = t135 * t136;
	t147 = t124 * t176 + t125 * t137;
	t92 = t147 * qJD(1) - qJD(3) * t114 + t124 * t158;
	t195 = (-t107 * t192 - t92 * t188) / t100 ^ 2;
	t95 = 0.1e1 / t97;
	t194 = t102 * t95;
	t180 = t134 * t137;
	t174 = qJD(1) * t135;
	t169 = 0.2e1 * t196;
	t168 = -0.2e1 * t195;
	t167 = 0.2e1 * t193;
	t166 = t101 * t196;
	t165 = t113 * t192;
	t164 = t95 * t171;
	t162 = t119 * t128 * t130;
	t156 = t130 * t200;
	t151 = t135 * t162;
	t149 = t157 * t137;
	t145 = t134 * t172 + t152 * t137;
	t112 = -t125 * t176 + t185;
	t106 = t119 * t198;
	t91 = (t197 * t134 * t115 - t116 * t151) * t137;
	t90 = -t115 * t176 + t186 + (t115 * t136 - t116 * t177) * t106;
	t89 = t198 * t200 + (qJD(1) * t149 + 0.2e1 * t135 * t146) * t119;
	t1 = [t156 * t180 + (qJD(2) * t149 - t174 * t181) * t119, t89, 0, 0, 0; (-t101 * t164 + (0.2e1 * t166 + (qJD(1) * t91 + t88) * t194) * t134) * t135 + ((-t91 * t164 + (t91 * t169 + ((0.2e1 * t134 * t190 - t94 * t151 - t197 * t171) * t115 + (t156 * t182 + t134 * t94 + (t127 * t159 - (t94 - 0.2e1 * t172) * t134) * t119) * t116) * t95 * t137) * t134) * t102 + (t91 * t167 + (-t101 + ((-t129 + t133) * t116 * t162 + t197 * t163) * t102) * qJD(1)) * t134 * t95) * t137, (-t101 * t95 * t174 + (-0.2e1 * t166 + (-qJD(2) * t90 - t88) * t194) * t137) * t136 + (t90 * t137 * t102 * t169 + ((-qJD(2) * t101 + t90 * t167) * t137 + (t90 * t174 + (-(-t106 * t173 - t135 * t89) * t116 - ((t106 * t135 - 0.1e1) * t94 + (-t106 + t135) * qJD(2)) * t115) * t180) * t102) * t95 - ((t89 - t173) * t115 + (t154 * t106 + t155) * t116) * t175 * t194) * t134, 0, 0, 0; 0.2e1 * (t108 * t147 + t112 * t188) * t195 + (0.2e1 * t112 * t165 - t153 * t108 * t178 + t145 * t189 + (-t153 * t113 * t179 + t112 * t92 - t145 * t187 + t147 * t93) * t109) * t98, t136 * t170 * t199 + (-t174 * t199 + (t148 * t168 + ((-qJD(3) * t108 - 0.2e1 * t165) * t125 + (-t125 * t92 + (-qJD(3) * t113 + t93) * t124) * t109) * t98) * t137) * t134, t168 + 0.2e1 * (-t109 * t92 * t98 + (-t109 * t195 - t98 * t192) * t113) * t113, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:03:05
	% EndTime: 2019-12-31 21:03:06
	% DurationCPUTime: 1.05s
	% Computational Cost: add. (4821->123), mult. (6168->268), div. (1114->15), fcn. (7752->9), ass. (0->114)
	t158 = qJ(3) + pkin(8);
	t157 = cos(t158);
	t240 = 0.2e1 * t157;
	t156 = sin(t158);
	t165 = cos(qJ(2));
	t166 = cos(qJ(1));
	t216 = t166 * t157;
	t236 = sin(qJ(1));
	t142 = t236 * t156 + t165 * t216;
	t136 = 0.1e1 / t142 ^ 2;
	t164 = sin(qJ(2));
	t159 = t164 ^ 2;
	t163 = t166 ^ 2;
	t221 = t159 * t163;
	t201 = t136 * t221;
	t132 = 0.1e1 + t201;
	t190 = qJD(1) * t236;
	t214 = qJD(2) * t165;
	t174 = t159 * t166 * t190 - t163 * t164 * t214;
	t213 = qJD(2) * t166;
	t194 = t164 * t213;
	t177 = t165 * t190 + t194;
	t189 = t236 * qJD(3);
	t217 = t166 * t156;
	t121 = (-qJD(3) * t165 + qJD(1)) * t217 + (t189 - t177) * t157;
	t135 = 0.1e1 / t142;
	t231 = t121 * t135 * t136;
	t184 = t221 * t231;
	t239 = (-t174 * t136 - t184) / t132 ^ 2;
	t219 = t164 * t166;
	t196 = t236 * t165;
	t138 = t156 * t196 + t216;
	t182 = t156 * t189;
	t210 = qJD(3) * t166;
	t192 = t157 * t210;
	t120 = t138 * qJD(1) + t156 * t194 - t165 * t192 - t182;
	t141 = -t236 * t157 + t165 * t217;
	t153 = 0.1e1 / t156;
	t154 = 0.1e1 / t156 ^ 2;
	t160 = 0.1e1 / t164;
	t161 = 0.1e1 / t164 ^ 2;
	t195 = t161 * t214;
	t212 = qJD(3) * t157;
	t224 = t153 * t160;
	t238 = (t154 * t160 * t212 + t153 * t195) * t141 + t120 * t224;
	t220 = t164 * t156;
	t128 = atan2(-t138, t220);
	t125 = cos(t128);
	t124 = sin(t128);
	t230 = t124 * t138;
	t119 = t125 * t220 - t230;
	t116 = 0.1e1 / t119;
	t117 = 0.1e1 / t119 ^ 2;
	t237 = 0.2e1 * t141;
	t133 = t138 ^ 2;
	t222 = t154 * t161;
	t129 = t133 * t222 + 0.1e1;
	t126 = 0.1e1 / t129;
	t211 = qJD(3) * t164;
	t178 = t156 * t214 + t157 * t211;
	t199 = t138 * t222;
	t197 = t236 * t164;
	t183 = qJD(2) * t197;
	t215 = qJD(1) * t166;
	t122 = (t189 * t165 - t190) * t157 + (t215 * t165 - t183 - t210) * t156;
	t202 = t122 * t224;
	t108 = (t178 * t199 - t202) * t126;
	t175 = -t108 * t138 + t178;
	t104 = (-t108 * t220 - t122) * t124 + t175 * t125;
	t118 = t116 * t117;
	t235 = t104 * t118;
	t155 = t153 * t154;
	t162 = t160 / t159;
	t193 = t161 * t212;
	t234 = (t122 * t199 + (-t154 * t162 * t214 - t155 * t193) * t133) / t129 ^ 2;
	t233 = t117 * t141;
	t232 = t120 * t117;
	t229 = t124 * t141;
	t228 = t124 * t164;
	t227 = t125 * t138;
	t226 = t125 * t141;
	t225 = t125 * t165;
	t223 = t154 * t157;
	t218 = t165 * t166;
	t134 = t141 ^ 2;
	t114 = t117 * t134 + 0.1e1;
	t209 = 0.2e1 * (-t134 * t235 - t141 * t232) / t114 ^ 2;
	t208 = -0.2e1 * t234;
	t207 = 0.2e1 * t239;
	t206 = t118 * t237;
	t205 = t160 * t234;
	t204 = t117 * t229;
	t200 = t138 * t224;
	t198 = t153 * t161 * t165;
	t180 = t138 * t198 + t236;
	t115 = t180 * t126;
	t191 = t236 - t115;
	t188 = t116 * t209;
	t187 = t117 * t209;
	t186 = t219 * t237;
	t185 = t153 * t205;
	t140 = t157 * t196 - t217;
	t181 = t138 * t223 - t140 * t153;
	t179 = t136 * t140 * t166 - t236 * t135;
	t130 = 0.1e1 / t132;
	t123 = t142 * qJD(1) - t157 * t183 - t165 * t182 - t192;
	t112 = 0.1e1 / t114;
	t111 = t181 * t160 * t126;
	t107 = (-t124 + (t125 * t200 + t124) * t126) * t141;
	t106 = -t115 * t227 + (t191 * t228 + t225) * t156;
	t105 = t125 * t157 * t164 - t124 * t140 + (-t124 * t220 - t227) * t111;
	t103 = t180 * t208 + (t122 * t198 + t215 + (-t154 * t165 * t193 + (-0.2e1 * t162 * t165 ^ 2 - t160) * t153 * qJD(2)) * t138) * t126;
	t101 = -0.2e1 * t181 * t205 + (-t181 * t195 + (t122 * t223 - t123 * t153 + (t140 * t223 + (-0.2e1 * t155 * t157 ^ 2 - t153) * t138) * qJD(3)) * t160) * t126;
	t1 = [t238 * t126 + t185 * t237, t103, t101, 0, 0; t138 * t188 + (-t122 * t116 + (t104 * t138 + t107 * t120) * t117) * t112 + (t107 * t187 + (0.2e1 * t107 * t235 + (t120 * t126 - t120 - (-t108 * t126 * t200 + t208) * t141) * t117 * t124 + (-(-0.2e1 * t138 * t185 - t108) * t233 + (-(t108 + t202) * t141 + t238 * t138) * t117 * t126) * t125) * t112) * t141, t106 * t141 * t187 + (-(-t103 * t227 + (t108 * t230 - t122 * t125) * t115) * t233 + (t104 * t206 + t232) * t106 + (-t116 * t219 - (-t115 * t228 + t124 * t197 + t225) * t233) * t212) * t112 + (t188 * t219 + ((-t116 * t213 - (t191 * qJD(2) - t108) * t204) * t165 + (t116 * t190 + (t166 * t104 - (-t103 + t215) * t229 - (t191 * t108 - qJD(2)) * t226) * t117) * t164) * t112) * t156, (t105 * t233 - t116 * t142) * t209 + (t105 * t232 + t121 * t116 + (t105 * t206 - t117 * t142) * t104 - (t157 * t214 - t156 * t211 - t101 * t138 - t111 * t122 + (-t111 * t220 - t140) * t108) * t117 * t226 - (-t123 + (-t101 * t156 - t108 * t157) * t164 - t175 * t111) * t204) * t112, 0, 0; t179 * t164 * t207 + (-t179 * t214 + ((qJD(1) * t135 + 0.2e1 * t140 * t231) * t166 + (-t236 * t121 - t123 * t166 + t140 * t190) * t136) * t164) * t130, (t135 * t218 + t157 * t201) * t207 + (t184 * t240 + t177 * t135 + (qJD(3) * t156 * t221 + t121 * t218 + t174 * t240) * t136) * t130, t136 * t186 * t239 + (t186 * t231 + (t120 * t219 + (t164 * t190 - t165 * t213) * t141) * t136) * t130, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end