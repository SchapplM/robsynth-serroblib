% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPR11
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
%   Wie in S5RRRPR11_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRPR11_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR11_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:25
	% EndTime: 2019-12-29 20:16:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:14
	% EndTime: 2019-12-29 20:16:14
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:19
	% EndTime: 2019-12-29 20:16:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:19
	% EndTime: 2019-12-29 20:16:20
	% DurationCPUTime: 1.51s
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
	t188 = (t123 * t125 * t155 + (-t123 * t187 - t151 * t92) * t118) / t90 ^ 2;
	t176 = t102 * t106;
	t143 = -qJD(1) * t128 + qJD(3);
	t144 = qJD(3) * t128 - qJD(1);
	t160 = qJD(2) * t129;
	t148 = t125 * t160;
	t87 = -t144 * t166 + (t126 * t143 - t148) * t127;
	t183 = t101 * t102 * t87;
	t167 = t126 * t128;
	t138 = t124 * t167 + t165;
	t86 = qJD(1) * t138 - t107 * qJD(3) + t124 * t148;
	t186 = (-t100 * t183 - t176 * t86) / t99 ^ 2;
	t88 = 0.1e1 / t90;
	t185 = t88 * t92;
	t184 = t91 * t88;
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
	t84 = (t108 * t125 * t189 - t109 * t142) * t129;
	t83 = -t192 * t174 + (-t126 + t98) * t128 * t108;
	t81 = t190 * t193 + (qJD(1) * t141 + 0.2e1 * t126 * t137) * t111;
	t1 = [t145 * t171 + (qJD(2) * t141 - t164 * t172) * t111, t81, 0, 0, 0; (-t161 * t184 + (0.2e1 * t157 + (qJD(1) * t84 + t80) * t185) * t125) * t126 + (t84 * t147 * t125 + (-t84 * t155 + (t84 * t159 + ((0.2e1 * t125 * t181 - t142 * t85 - t161 * t189) * t108 + (t118 * t126 * t145 + t125 * t85 + (t117 * t149 - (t85 - 0.2e1 * t162) * t125) * t111) * t109) * t179) * t125 + (-t91 + (-(t119 - t123) * t109 * t153 + t189 * t154) * t92) * t125 * qJD(1)) * t88) * t129, (-t164 * t184 + (-0.2e1 * t157 + (-qJD(2) * t83 - t80) * t185) * t129) * t128 + (t83 * t129 * t147 + (-t91 * t160 - ((-t126 * t81 - t163 * t98) * t109 + (t192 * t85 + t162 - t178) * t108) * t92 * t171 + (t129 * t159 + t164 * t92) * t83 - ((t81 - t163) * t108 + (t85 * t98 + qJD(2) + (-t85 - t178) * t126) * t109) * t128 * t179) * t88) * t125, 0, 0, 0; 0.2e1 * (t101 * t138 + t105 * t176) * t186 + (0.2e1 * t105 * t156 - t144 * t101 * t168 + (t125 * t162 + t140) * t177 + (t138 * t87 + t105 * t86 - t140 * t175 - (qJD(2) * t125 * t127 + t124 * t144) * t106 * t126) * t102) * t96, t128 * t160 * t191 + (-t164 * t191 + (t139 * t158 + ((-qJD(3) * t101 - 0.2e1 * t156) * t127 + (-t127 * t86 + (-qJD(3) * t106 + t87) * t124) * t102) * t96) * t129) * t125, t158 + 0.2e1 * (-t86 * t102 * t96 + (-t102 * t186 - t183 * t96) * t106) * t106, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:26
	% EndTime: 2019-12-29 20:16:28
	% DurationCPUTime: 1.99s
	% Computational Cost: add. (1824->122), mult. (6168->269), div. (1114->15), fcn. (7752->9), ass. (0->112)
	t156 = sin(qJ(3));
	t159 = cos(qJ(2));
	t158 = cos(qJ(3));
	t160 = cos(qJ(1));
	t212 = t160 * t158;
	t231 = sin(qJ(1));
	t137 = t231 * t156 + t159 * t212;
	t131 = 0.1e1 / t137 ^ 2;
	t157 = sin(qJ(2));
	t151 = t157 ^ 2;
	t155 = t160 ^ 2;
	t216 = t151 * t155;
	t194 = t131 * t216;
	t126 = 0.1e1 + t194;
	t185 = qJD(1) * t231;
	t209 = qJD(2) * t160;
	t189 = t157 * t209;
	t170 = t159 * t185 + t189;
	t184 = t231 * qJD(3);
	t213 = t160 * t156;
	t116 = (-qJD(3) * t159 + qJD(1)) * t213 + (t184 - t170) * t158;
	t130 = 0.1e1 / t137;
	t226 = t116 * t130 * t131;
	t179 = t216 * t226;
	t190 = qJD(2) * t155 * t157;
	t234 = (-t179 + (-t151 * t160 * t185 + t159 * t190) * t131) / t126 ^ 2;
	t214 = t157 * t160;
	t192 = t231 * t159;
	t133 = t156 * t192 + t212;
	t176 = t156 * t184;
	t206 = qJD(3) * t160;
	t187 = t158 * t206;
	t115 = t133 * qJD(1) + t156 * t189 - t159 * t187 - t176;
	t136 = -t231 * t158 + t159 * t213;
	t148 = 0.1e1 / t156;
	t149 = 0.1e1 / t156 ^ 2;
	t152 = 0.1e1 / t157;
	t153 = 0.1e1 / t157 ^ 2;
	t210 = qJD(2) * t159;
	t191 = t153 * t210;
	t207 = qJD(3) * t158;
	t219 = t148 * t152;
	t233 = (t149 * t152 * t207 + t148 * t191) * t136 + t115 * t219;
	t215 = t157 * t156;
	t125 = atan2(-t133, t215);
	t120 = cos(t125);
	t119 = sin(t125);
	t225 = t119 * t133;
	t114 = t120 * t215 - t225;
	t111 = 0.1e1 / t114;
	t112 = 0.1e1 / t114 ^ 2;
	t232 = 0.2e1 * t136;
	t128 = t133 ^ 2;
	t218 = t149 * t153;
	t127 = t128 * t218 + 0.1e1;
	t123 = 0.1e1 / t127;
	t171 = t156 * t210 + t157 * t207;
	t196 = t133 * t218;
	t193 = t231 * t157;
	t177 = qJD(2) * t193;
	t178 = t158 * t185;
	t211 = qJD(1) * t160;
	t117 = t158 * t184 * t159 - t178 + (t211 * t159 - t177 - t206) * t156;
	t198 = t117 * t219;
	t103 = (t171 * t196 - t198) * t123;
	t168 = -t103 * t133 + t171;
	t99 = (-t103 * t215 - t117) * t119 + t168 * t120;
	t230 = t111 * t112 * t99;
	t150 = t148 * t149;
	t154 = t152 / t151;
	t188 = t153 * t207;
	t229 = (t117 * t196 + (-t149 * t154 * t210 - t150 * t188) * t128) / t127 ^ 2;
	t228 = t112 * t136;
	t227 = t115 * t112;
	t224 = t119 * t136;
	t223 = t119 * t157;
	t222 = t120 * t133;
	t221 = t120 * t136;
	t220 = t120 * t159;
	t217 = t149 * t158;
	t208 = qJD(3) * t156;
	t129 = t136 ^ 2;
	t109 = t129 * t112 + 0.1e1;
	t205 = 0.2e1 / t109 ^ 2 * (-t129 * t230 - t136 * t227);
	t204 = 0.2e1 * t230;
	t203 = 0.2e1 * t234;
	t202 = -0.2e1 * t229;
	t201 = t152 * t229;
	t200 = t112 * t224;
	t197 = t133 * t219;
	t195 = t148 * t153 * t159;
	t173 = t133 * t195 + t231;
	t110 = t173 * t123;
	t186 = t231 - t110;
	t183 = t111 * t205;
	t182 = t112 * t205;
	t181 = t214 * t232;
	t180 = t148 * t201;
	t175 = t136 * t204 + t227;
	t135 = t158 * t192 - t213;
	t174 = t133 * t217 - t135 * t148;
	t172 = t131 * t135 * t160 - t231 * t130;
	t121 = 0.1e1 / t126;
	t118 = t137 * qJD(1) - t158 * t177 - t159 * t176 - t187;
	t107 = 0.1e1 / t109;
	t106 = t174 * t152 * t123;
	t102 = (-t119 + (t120 * t197 + t119) * t123) * t136;
	t101 = -t110 * t222 + (t186 * t223 + t220) * t156;
	t100 = t120 * t157 * t158 - t119 * t135 + (-t119 * t215 - t222) * t106;
	t98 = t173 * t202 + (t117 * t195 + t211 + (-t149 * t159 * t188 + (-0.2e1 * t154 * t159 ^ 2 - t152) * t148 * qJD(2)) * t133) * t123;
	t96 = -0.2e1 * t174 * t201 + (-t174 * t191 + (t117 * t217 - t118 * t148 + (t135 * t217 + (-0.2e1 * t150 * t158 ^ 2 - t148) * t133) * qJD(3)) * t152) * t123;
	t1 = [t233 * t123 + t180 * t232, t98, t96, 0, 0; t133 * t183 + (-t117 * t111 + (t102 * t115 + t133 * t99) * t112) * t107 + (t102 * t182 + (t102 * t204 + (t115 * t123 - t115 - (-t103 * t123 * t197 + t202) * t136) * t112 * t119 + (-(-0.2e1 * t133 * t180 - t103) * t228 + (-(t103 + t198) * t136 + t233 * t133) * t112 * t123) * t120) * t107) * t136, t101 * t136 * t182 + (-(-t98 * t222 + (t103 * t225 - t117 * t120) * t110) * t228 + t175 * t101 + (-t111 * t214 - (-t110 * t223 + t119 * t193 + t220) * t228) * t207) * t107 + (t183 * t214 + ((-t111 * t209 - (t186 * qJD(2) - t103) * t200) * t159 + (t111 * t185 + (t160 * t99 - (-t98 + t211) * t224 - (t186 * t103 - qJD(2)) * t221) * t112) * t157) * t107) * t156, (t100 * t228 - t111 * t137) * t205 + (t116 * t111 + t175 * t100 - (-t118 + (-t103 * t158 - t156 * t96) * t157 - t168 * t106) * t200 + (-t137 * t99 - (t158 * t210 - t157 * t208 - t106 * t117 - t133 * t96 + (-t106 * t215 - t135) * t103) * t221) * t112) * t107, 0, 0; t172 * t157 * t203 + (-t172 * t210 + ((qJD(1) * t130 + 0.2e1 * t135 * t226) * t160 + (-t231 * t116 - t118 * t160 + t135 * t185) * t131) * t157) * t121, (t130 * t159 * t160 + t158 * t194) * t203 + (0.2e1 * t158 * t179 + t170 * t130 + ((t116 * t160 - 0.2e1 * t158 * t190) * t159 + (t155 * t208 + 0.2e1 * t160 * t178) * t151) * t131) * t121, t131 * t181 * t234 + (t181 * t226 + (t115 * t214 + (t157 * t185 - t159 * t209) * t136) * t131) * t121, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:26
	% EndTime: 2019-12-29 20:16:28
	% DurationCPUTime: 2.14s
	% Computational Cost: add. (1265->118), mult. (4276->253), div. (493->12), fcn. (5073->11), ass. (0->110)
	t180 = sin(qJ(2));
	t172 = t180 ^ 2;
	t184 = cos(qJ(2));
	t175 = 0.1e1 / t184 ^ 2;
	t237 = t172 * t175;
	t181 = sin(qJ(1));
	t173 = t181 ^ 2;
	t167 = t173 * t237 + 0.1e1;
	t174 = 0.1e1 / t184;
	t234 = t174 * t180;
	t256 = t180 * t237;
	t193 = qJD(2) * (t174 * t256 + t234);
	t185 = cos(qJ(1));
	t226 = qJD(1) * t185;
	t235 = t172 * t181;
	t202 = t226 * t235;
	t240 = (t173 * t193 + t175 * t202) / t167 ^ 2;
	t257 = -0.2e1 * t240;
	t209 = 0.1e1 + t237;
	t255 = t181 * t209;
	t183 = cos(qJ(3));
	t254 = (-qJD(3) + qJD(5)) * t183;
	t179 = sin(qJ(3));
	t222 = qJD(3) * t179;
	t253 = -qJD(5) * t179 + t222;
	t229 = t181 * t184;
	t194 = t179 * t229 + t183 * t185;
	t223 = qJD(2) * t185;
	t210 = t180 * t223;
	t228 = t184 * t185;
	t212 = t183 * t228;
	t132 = t194 * qJD(1) - qJD(3) * t212 + t179 * t210 - t181 * t222;
	t204 = -qJD(1) * t184 + qJD(3);
	t205 = qJD(3) * t184 - qJD(1);
	t233 = t179 * t185;
	t133 = -t205 * t233 + (t204 * t181 - t210) * t183;
	t178 = sin(qJ(5));
	t182 = cos(qJ(5));
	t230 = t181 * t183;
	t160 = t179 * t228 - t230;
	t161 = t181 * t179 + t212;
	t198 = t160 * t182 - t161 * t178;
	t125 = t198 * qJD(5) - t132 * t178 + t133 * t182;
	t148 = t160 * t178 + t161 * t182;
	t140 = 0.1e1 / t148;
	t196 = t178 * t179 + t182 * t183;
	t197 = t178 * t183 - t179 * t182;
	t141 = 0.1e1 / t148 ^ 2;
	t243 = t141 * t198;
	t252 = t197 * t140 + t196 * t243;
	t231 = t181 * t180;
	t166 = atan2(t231, t184);
	t163 = cos(t166);
	t162 = sin(t166);
	t214 = t162 * t231;
	t152 = t163 * t184 + t214;
	t149 = 0.1e1 / t152;
	t150 = 0.1e1 / t152 ^ 2;
	t251 = 0.2e1 * t180;
	t164 = 0.1e1 / t167;
	t250 = t164 - 0.1e1;
	t124 = t148 * qJD(5) + t132 * t182 + t133 * t178;
	t139 = t198 ^ 2;
	t130 = t139 * t141 + 0.1e1;
	t142 = t140 * t141;
	t246 = t125 * t142;
	t249 = (-t124 * t243 - t139 * t246) / t130 ^ 2;
	t177 = t185 ^ 2;
	t236 = t172 * t177;
	t138 = t150 * t236 + 0.1e1;
	t224 = qJD(2) * t184;
	t211 = t180 * t226;
	t225 = qJD(2) * t181;
	t131 = ((t181 * t224 + t211) * t174 + t225 * t237) * t164;
	t238 = t163 * t180;
	t122 = (t131 * t181 - qJD(2)) * t238 + (t211 + (-t131 + t225) * t184) * t162;
	t247 = t122 * t149 * t150;
	t248 = (-t236 * t247 + (t177 * t180 * t224 - t202) * t150) / t138 ^ 2;
	t245 = t131 * t162;
	t244 = t131 * t180;
	t232 = t180 * t185;
	t156 = t196 * t232;
	t242 = t141 * t156;
	t241 = t150 * t180;
	t154 = t164 * t255;
	t239 = t154 * t181;
	t227 = qJD(1) * t181;
	t218 = 0.2e1 * t249;
	t217 = -0.2e1 * t247;
	t216 = -0.2e1 * t142 * t198;
	t215 = t150 * t232;
	t213 = t164 * t172 * t174;
	t208 = -0.2e1 * t180 * t248;
	t207 = t125 * t216;
	t206 = t174 * t257;
	t203 = t181 * t213;
	t201 = t209 * t185;
	t159 = -t183 * t229 + t233;
	t199 = -t159 * t178 - t182 * t194;
	t144 = t159 * t182 - t178 * t194;
	t195 = t204 * t185;
	t155 = t197 * t232;
	t136 = 0.1e1 / t138;
	t135 = t183 * t195 + (qJD(2) * t180 * t183 + t205 * t179) * t181;
	t134 = -t205 * t230 + (t180 * t225 + t195) * t179;
	t128 = 0.1e1 / t130;
	t127 = (-t250 * t180 * t162 + t163 * t203) * t185;
	t126 = t162 * t229 - t238 + (-t162 * t184 + t163 * t231) * t154;
	t123 = t255 * t257 + (qJD(1) * t201 + 0.2e1 * t181 * t193) * t164;
	t1 = [t206 * t232 + (qJD(2) * t201 - t227 * t234) * t164, t123, 0, 0, 0; (t149 * t208 + (t149 * t224 + (-qJD(1) * t127 - t122) * t241) * t136) * t181 + (t150 * t208 * t127 + (((-t131 * t203 - t250 * t224 + t240 * t251) * t162 + (t206 * t235 + t244 + (-t244 + (t251 + t256) * t225) * t164) * t163) * t215 + (t150 * t224 + t180 * t217) * t127 + (t149 + ((-t173 + t177) * t163 * t213 + t250 * t214) * t150) * t180 * qJD(1)) * t136) * t185, 0.2e1 * (-t126 * t241 + t149 * t184) * t185 * t248 + ((t149 * t227 + (qJD(2) * t126 + t122) * t185 * t150) * t184 + (t149 * t223 + (t123 * t163 * t181 - t162 * t225 - t239 * t245 + t245 + (qJD(2) * t162 + t163 * t226) * t154) * t215 + (-t150 * t227 + t185 * t217) * t126 + ((-t123 + t226) * t162 + ((-0.1e1 + t239) * qJD(2) + (-t154 + t181) * t131) * t163) * t150 * t228) * t180) * t136, 0, 0, 0; (t140 * t199 - t144 * t243) * t218 + ((t144 * qJD(5) - t134 * t182 + t135 * t178) * t140 + t144 * t207 + (t199 * t125 + (t199 * qJD(5) + t134 * t178 + t135 * t182) * t198 - t144 * t124) * t141) * t128, (t140 * t155 + t198 * t242) * t218 + (t124 * t242 + (t141 * t155 - t156 * t216) * t125 - t252 * t184 * t223 + (t252 * t227 + ((-t254 * t140 + t253 * t243) * t182 + (t253 * t140 + t254 * t243) * t178) * t185) * t180) * t128, (t140 * t148 + t198 * t243) * t218 + (-t125 * t140 - t198 * t207 + (0.2e1 * t198 * t124 + t148 * t125) * t141) * t128, 0, -0.2e1 * t249 - 0.2e1 * (t124 * t128 * t141 - (-t128 * t246 - t141 * t249) * t198) * t198;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end