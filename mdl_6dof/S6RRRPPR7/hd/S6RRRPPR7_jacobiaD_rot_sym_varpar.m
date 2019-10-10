% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR7
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
%   Wie in S6RRRPPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:27
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:45
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:45
	% EndTime: 2019-10-10 11:27:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:46
	% EndTime: 2019-10-10 11:27:47
	% DurationCPUTime: 1.02s
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
	t139 = t175 * t102 - t177;
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
	t1 = [t145 * t171 + (qJD(2) * t141 - t164 * t172) * t111, t81, 0, 0, 0, 0; (-t161 * t184 + (0.2e1 * t157 + (qJD(1) * t84 + t80) * t185) * t125) * t126 + (t84 * t147 * t125 + (-t84 * t155 + (t84 * t159 + ((0.2e1 * t125 * t181 - t142 * t85 - t161 * t189) * t108 + (t118 * t126 * t145 + t125 * t85 + (t117 * t149 - (t85 - 0.2e1 * t162) * t125) * t111) * t109) * t179) * t125 + (-t91 + (-(t119 - t123) * t109 * t153 + t189 * t154) * t92) * t125 * qJD(1)) * t88) * t129, (-t164 * t184 + (-0.2e1 * t157 + (-qJD(2) * t83 - t80) * t185) * t129) * t128 + (t83 * t129 * t147 + (-t91 * t160 - ((-t126 * t81 - t163 * t98) * t109 + (t192 * t85 + t162 - t178) * t108) * t92 * t171 + (t129 * t159 + t164 * t92) * t83 - ((t81 - t163) * t108 + (t85 * t98 + qJD(2) + (-t85 - t178) * t126) * t109) * t128 * t179) * t88) * t125, 0, 0, 0, 0; 0.2e1 * (t101 * t138 + t105 * t176) * t186 + (0.2e1 * t105 * t156 - t144 * t101 * t168 + (t162 * t125 + t140) * t177 + (t138 * t87 + t105 * t86 - t140 * t175 - (qJD(2) * t125 * t127 + t124 * t144) * t106 * t126) * t102) * t96, t128 * t160 * t191 + (-t164 * t191 + (t139 * t158 + ((-qJD(3) * t101 - 0.2e1 * t156) * t127 + (-t127 * t86 + (-qJD(3) * t106 + t87) * t124) * t102) * t96) * t129) * t125, t158 + 0.2e1 * (-t86 * t102 * t96 + (-t102 * t186 - t183 * t96) * t106) * t106, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:46
	% EndTime: 2019-10-10 11:27:47
	% DurationCPUTime: 1.35s
	% Computational Cost: add. (1824->123), mult. (6168->270), div. (1114->15), fcn. (7752->9), ass. (0->112)
	t155 = sin(qJ(3));
	t158 = cos(qJ(2));
	t157 = cos(qJ(3));
	t159 = cos(qJ(1));
	t211 = t159 * t157;
	t230 = sin(qJ(1));
	t136 = t230 * t155 + t158 * t211;
	t130 = 0.1e1 / t136 ^ 2;
	t156 = sin(qJ(2));
	t150 = t156 ^ 2;
	t154 = t159 ^ 2;
	t215 = t150 * t154;
	t193 = t130 * t215;
	t125 = 0.1e1 + t193;
	t184 = qJD(1) * t230;
	t208 = qJD(2) * t159;
	t188 = t156 * t208;
	t169 = t158 * t184 + t188;
	t183 = t230 * qJD(3);
	t212 = t159 * t155;
	t115 = (-qJD(3) * t158 + qJD(1)) * t212 + (t183 - t169) * t157;
	t129 = 0.1e1 / t136;
	t225 = t115 * t129 * t130;
	t177 = t215 * t225;
	t189 = qJD(2) * t154 * t156;
	t233 = (-t177 + (-t150 * t159 * t184 + t158 * t189) * t130) / t125 ^ 2;
	t213 = t156 * t159;
	t191 = t230 * t158;
	t132 = t155 * t191 + t211;
	t174 = t155 * t183;
	t205 = qJD(3) * t159;
	t186 = t157 * t205;
	t114 = t132 * qJD(1) + t155 * t188 - t158 * t186 - t174;
	t135 = -t230 * t157 + t158 * t212;
	t147 = 0.1e1 / t155;
	t148 = 0.1e1 / t155 ^ 2;
	t151 = 0.1e1 / t156;
	t152 = 0.1e1 / t156 ^ 2;
	t209 = qJD(2) * t158;
	t190 = t152 * t209;
	t206 = qJD(3) * t157;
	t218 = t147 * t151;
	t232 = (t148 * t151 * t206 + t147 * t190) * t135 + t114 * t218;
	t214 = t156 * t155;
	t124 = atan2(-t132, t214);
	t119 = cos(t124);
	t118 = sin(t124);
	t224 = t118 * t132;
	t113 = t119 * t214 - t224;
	t110 = 0.1e1 / t113;
	t111 = 0.1e1 / t113 ^ 2;
	t231 = 0.2e1 * t135;
	t127 = t132 ^ 2;
	t217 = t148 * t152;
	t126 = t127 * t217 + 0.1e1;
	t122 = 0.1e1 / t126;
	t170 = t155 * t209 + t156 * t206;
	t195 = t132 * t217;
	t192 = t230 * t156;
	t175 = qJD(2) * t192;
	t176 = t157 * t184;
	t210 = qJD(1) * t159;
	t116 = t157 * t183 * t158 - t176 + (t210 * t158 - t175 - t205) * t155;
	t197 = t116 * t218;
	t102 = (t170 * t195 - t197) * t122;
	t167 = -t102 * t132 + t170;
	t98 = (-t102 * t214 - t116) * t118 + t167 * t119;
	t229 = t110 * t111 * t98;
	t149 = t147 * t148;
	t153 = t151 / t150;
	t187 = t152 * t206;
	t228 = (t116 * t195 + (-t148 * t153 * t209 - t149 * t187) * t127) / t126 ^ 2;
	t227 = t111 * t135;
	t226 = t114 * t111;
	t223 = t118 * t135;
	t222 = t118 * t156;
	t221 = t119 * t132;
	t220 = t119 * t135;
	t219 = t119 * t158;
	t216 = t148 * t157;
	t207 = qJD(3) * t155;
	t128 = t135 ^ 2;
	t108 = t128 * t111 + 0.1e1;
	t204 = 0.2e1 / t108 ^ 2 * (-t128 * t229 - t135 * t226);
	t203 = 0.2e1 * t229;
	t202 = 0.2e1 * t233;
	t201 = -0.2e1 * t228;
	t200 = t151 * t228;
	t199 = t111 * t223;
	t196 = t132 * t218;
	t194 = t147 * t152 * t158;
	t172 = t132 * t194 + t230;
	t109 = t172 * t122;
	t185 = t230 - t109;
	t182 = t110 * t204;
	t181 = t111 * t204;
	t180 = t135 * t203;
	t179 = t213 * t231;
	t178 = t147 * t200;
	t134 = t157 * t191 - t212;
	t173 = t132 * t216 - t134 * t147;
	t171 = t130 * t134 * t159 - t230 * t129;
	t120 = 0.1e1 / t125;
	t117 = t136 * qJD(1) - t157 * t175 - t158 * t174 - t186;
	t106 = 0.1e1 / t108;
	t105 = t173 * t151 * t122;
	t101 = (-t118 + (t119 * t196 + t118) * t122) * t135;
	t100 = -t109 * t221 + (t185 * t222 + t219) * t155;
	t99 = t119 * t156 * t157 - t118 * t134 + (-t118 * t214 - t221) * t105;
	t97 = t172 * t201 + (t116 * t194 + t210 + (-t148 * t158 * t187 + (-0.2e1 * t153 * t158 ^ 2 - t151) * t147 * qJD(2)) * t132) * t122;
	t95 = -0.2e1 * t173 * t200 + (-t173 * t190 + (t116 * t216 - t117 * t147 + (t134 * t216 + (-0.2e1 * t149 * t157 ^ 2 - t147) * t132) * qJD(3)) * t151) * t122;
	t1 = [t232 * t122 + t178 * t231, t97, t95, 0, 0, 0; t132 * t182 + (-t116 * t110 + (t101 * t114 + t132 * t98) * t111) * t106 + (t101 * t181 + (t101 * t203 + (t114 * t122 - t114 - (-t102 * t122 * t196 + t201) * t135) * t111 * t118 + (-(-0.2e1 * t132 * t178 - t102) * t227 + (-(t102 + t197) * t135 + t232 * t132) * t111 * t122) * t119) * t106) * t135, t100 * t135 * t181 + (-(-t97 * t221 + (t102 * t224 - t116 * t119) * t109) * t227 + (t180 + t226) * t100 + (-t110 * t213 - (-t109 * t222 + t118 * t192 + t219) * t227) * t206) * t106 + (t182 * t213 + ((-t110 * t208 - (t185 * qJD(2) - t102) * t199) * t158 + (t110 * t184 + (t159 * t98 - (-t97 + t210) * t223 - (t185 * t102 - qJD(2)) * t220) * t111) * t156) * t106) * t155, (-t110 * t136 + t99 * t227) * t204 + (t99 * t180 + t115 * t110 - (-t117 + (-t102 * t157 - t155 * t95) * t156 - t167 * t105) * t199 + (t99 * t114 - t136 * t98 - (t157 * t209 - t156 * t207 - t105 * t116 - t132 * t95 + (-t105 * t214 - t134) * t102) * t220) * t111) * t106, 0, 0, 0; t171 * t156 * t202 + (-t171 * t209 + ((qJD(1) * t129 + 0.2e1 * t134 * t225) * t159 + (-t230 * t115 - t117 * t159 + t134 * t184) * t130) * t156) * t120, (t129 * t158 * t159 + t157 * t193) * t202 + (0.2e1 * t157 * t177 + t169 * t129 + ((t115 * t159 - 0.2e1 * t157 * t189) * t158 + (t154 * t207 + 0.2e1 * t159 * t176) * t150) * t130) * t120, t130 * t179 * t233 + (t179 * t225 + (t114 * t213 + (t156 * t184 - t158 * t208) * t135) * t130) * t120, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:46
	% EndTime: 2019-10-10 11:27:47
	% DurationCPUTime: 1.26s
	% Computational Cost: add. (975->110), mult. (3480->241), div. (475->12), fcn. (4141->11), ass. (0->108)
	t169 = sin(qJ(2));
	t160 = t169 ^ 2;
	t172 = cos(qJ(2));
	t163 = 0.1e1 / t172 ^ 2;
	t225 = t160 * t163;
	t170 = sin(qJ(1));
	t161 = t170 ^ 2;
	t155 = t161 * t225 + 0.1e1;
	t162 = 0.1e1 / t172;
	t222 = t162 * t169;
	t241 = t169 * t225;
	t181 = qJD(2) * (t162 * t241 + t222);
	t173 = cos(qJ(1));
	t213 = qJD(1) * t173;
	t223 = t160 * t170;
	t188 = t213 * t223;
	t228 = (t161 * t181 + t163 * t188) / t155 ^ 2;
	t242 = -0.2e1 * t228;
	t196 = 0.1e1 + t225;
	t240 = t170 * t196;
	t168 = sin(qJ(3));
	t171 = cos(qJ(3));
	t209 = qJD(3) * t173;
	t197 = t171 * t209;
	t214 = qJD(1) * t170;
	t239 = -t168 * t214 + t197;
	t238 = t168 * t209 + t171 * t214;
	t215 = t172 * t173;
	t217 = t170 * t171;
	t148 = t168 * t215 - t217;
	t219 = t170 * t168;
	t149 = t171 * t215 + t219;
	t166 = sin(pkin(10));
	t167 = cos(pkin(10));
	t136 = t148 * t166 + t149 * t167;
	t128 = 0.1e1 / t136;
	t218 = t170 * t169;
	t154 = atan2(t218, t172);
	t151 = cos(t154);
	t150 = sin(t154);
	t204 = t150 * t218;
	t140 = t151 * t172 + t204;
	t137 = 0.1e1 / t140;
	t129 = 0.1e1 / t136 ^ 2;
	t138 = 0.1e1 / t140 ^ 2;
	t237 = 0.2e1 * t169;
	t152 = 0.1e1 / t155;
	t236 = t152 - 0.1e1;
	t165 = t173 ^ 2;
	t224 = t160 * t165;
	t126 = t138 * t224 + 0.1e1;
	t211 = qJD(2) * t172;
	t202 = t169 * t213;
	t212 = qJD(2) * t170;
	t118 = ((t170 * t211 + t202) * t162 + t212 * t225) * t152;
	t226 = t151 * t169;
	t109 = (t118 * t170 - qJD(2)) * t226 + (t202 + (-t118 + t212) * t172) * t150;
	t234 = t109 * t137 * t138;
	t235 = (-t224 * t234 + (t165 * t169 * t211 - t188) * t138) / t126 ^ 2;
	t233 = t118 * t150;
	t232 = t118 * t169;
	t192 = -t148 * t167 + t149 * t166;
	t231 = t129 * t192;
	t184 = -t166 * t168 - t167 * t171;
	t220 = t169 * t173;
	t144 = t184 * t220;
	t230 = t129 * t144;
	t229 = t138 * t169;
	t142 = t152 * t240;
	t227 = t142 * t170;
	t221 = t168 * t173;
	t216 = t170 * t172;
	t210 = qJD(2) * t173;
	t182 = t168 * t216 + t171 * t173;
	t199 = t169 * t210;
	t120 = t182 * qJD(1) - qJD(3) * t219 + t168 * t199 - t172 * t197;
	t190 = -qJD(1) * t172 + qJD(3);
	t191 = qJD(3) * t172 - qJD(1);
	t121 = -t191 * t221 + (t190 * t170 - t199) * t171;
	t112 = t120 * t167 + t121 * t166;
	t127 = t192 ^ 2;
	t116 = t127 * t129 + 0.1e1;
	t130 = t128 * t129;
	t186 = t120 * t166 - t121 * t167;
	t208 = 0.2e1 * (t127 * t130 * t186 + t112 * t231) / t116 ^ 2;
	t207 = -0.2e1 * t234;
	t206 = 0.2e1 * t130 * t192;
	t205 = t138 * t220;
	t203 = t152 * t160 * t162;
	t195 = -0.2e1 * t169 * t235;
	t194 = t186 * t206;
	t193 = t162 * t242;
	t189 = t170 * t203;
	t187 = t196 * t173;
	t185 = -t166 * t171 + t167 * t168;
	t183 = t190 * t173;
	t147 = -t171 * t216 + t221;
	t143 = t185 * t220;
	t132 = t147 * t167 - t166 * t182;
	t131 = t147 * t166 + t167 * t182;
	t124 = 0.1e1 / t126;
	t123 = t171 * t183 + (qJD(2) * t169 * t171 + t191 * t168) * t170;
	t122 = -t191 * t217 + (t169 * t212 + t183) * t168;
	t117 = (-t236 * t169 * t150 + t151 * t189) * t173;
	t114 = 0.1e1 / t116;
	t111 = t150 * t216 - t226 + (-t150 * t172 + t151 * t218) * t142;
	t110 = t240 * t242 + (qJD(1) * t187 + 0.2e1 * t170 * t181) * t152;
	t1 = [t193 * t220 + (qJD(2) * t187 - t214 * t222) * t152, t110, 0, 0, 0, 0; (t137 * t195 + (t137 * t211 + (-qJD(1) * t117 - t109) * t229) * t124) * t170 + (t138 * t195 * t117 + (((-t118 * t189 - t236 * t211 + t228 * t237) * t150 + (t193 * t223 + t232 + (-t232 + (t237 + t241) * t212) * t152) * t151) * t205 + (t138 * t211 + t169 * t207) * t117 + (t137 + ((-t161 + t165) * t151 * t203 + t236 * t204) * t138) * t169 * qJD(1)) * t124) * t173, 0.2e1 * (-t111 * t229 + t137 * t172) * t173 * t235 + ((t137 * t214 + (qJD(2) * t111 + t109) * t173 * t138) * t172 + (t137 * t210 + (t110 * t151 * t170 - t150 * t212 - t227 * t233 + t233 + (qJD(2) * t150 + t151 * t213) * t142) * t205 + (-t138 * t214 + t173 * t207) * t111 + ((-t110 + t213) * t150 + ((-0.1e1 + t227) * qJD(2) + (-t142 + t170) * t118) * t151) * t138 * t215) * t169) * t124, 0, 0, 0, 0; (-t128 * t131 + t132 * t231) * t208 + ((-t122 * t167 + t123 * t166) * t128 - t132 * t194 + (t131 * t186 - (t122 * t166 + t123 * t167) * t192 - t132 * t112) * t129) * t114, (-t128 * t143 + t192 * t230) * t208 + (-t112 * t230 - (-t129 * t143 + t144 * t206) * t186 + (t185 * t128 - t184 * t231) * t172 * t210 + ((t239 * t128 - t238 * t231) * t167 + (t238 * t128 + t239 * t231) * t166) * t169) * t114, (t128 * t136 + t192 * t231) * t208 + (t186 * t128 - t192 * t194 + (-0.2e1 * t192 * t112 - t136 * t186) * t129) * t114, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:27:46
	% EndTime: 2019-10-10 11:27:47
	% DurationCPUTime: 1.46s
	% Computational Cost: add. (1687->119), mult. (4276->252), div. (493->12), fcn. (5073->11), ass. (0->112)
	t190 = sin(qJ(2));
	t183 = t190 ^ 2;
	t193 = cos(qJ(2));
	t186 = 0.1e1 / t193 ^ 2;
	t246 = t183 * t186;
	t191 = sin(qJ(1));
	t184 = t191 ^ 2;
	t175 = t184 * t246 + 0.1e1;
	t185 = 0.1e1 / t193;
	t243 = t185 * t190;
	t266 = t190 * t246;
	t202 = qJD(2) * (t185 * t266 + t243);
	t194 = cos(qJ(1));
	t235 = qJD(1) * t194;
	t244 = t183 * t191;
	t211 = t235 * t244;
	t249 = (t184 * t202 + t186 * t211) / t175 ^ 2;
	t267 = -0.2e1 * t249;
	t218 = 0.1e1 + t246;
	t265 = t191 * t218;
	t192 = cos(qJ(3));
	t264 = (-qJD(3) + qJD(6)) * t192;
	t189 = sin(qJ(3));
	t231 = qJD(3) * t189;
	t263 = -qJD(6) * t189 + t231;
	t237 = t194 * t192;
	t239 = t191 * t193;
	t203 = t189 * t239 + t237;
	t232 = qJD(2) * t194;
	t219 = t190 * t232;
	t221 = t193 * t237;
	t140 = t203 * qJD(1) - qJD(3) * t221 + t189 * t219 - t191 * t231;
	t213 = -qJD(1) * t193 + qJD(3);
	t214 = qJD(3) * t193 - qJD(1);
	t238 = t194 * t189;
	t141 = -t214 * t238 + (t213 * t191 - t219) * t192;
	t181 = pkin(10) + qJ(6);
	t179 = sin(t181);
	t180 = cos(t181);
	t240 = t191 * t192;
	t168 = t193 * t238 - t240;
	t169 = t191 * t189 + t221;
	t207 = t168 * t180 - t169 * t179;
	t133 = t207 * qJD(6) - t140 * t179 + t141 * t180;
	t153 = t168 * t179 + t169 * t180;
	t145 = 0.1e1 / t153;
	t205 = t179 * t189 + t180 * t192;
	t206 = t179 * t192 - t180 * t189;
	t146 = 0.1e1 / t153 ^ 2;
	t253 = t146 * t207;
	t262 = t206 * t145 + t205 * t253;
	t241 = t191 * t190;
	t174 = atan2(t241, t193);
	t171 = cos(t174);
	t170 = sin(t174);
	t223 = t170 * t241;
	t160 = t171 * t193 + t223;
	t157 = 0.1e1 / t160;
	t158 = 0.1e1 / t160 ^ 2;
	t261 = 0.2e1 * t190;
	t172 = 0.1e1 / t175;
	t260 = t172 - 0.1e1;
	t132 = t153 * qJD(6) + t140 * t180 + t141 * t179;
	t144 = t207 ^ 2;
	t137 = t144 * t146 + 0.1e1;
	t147 = t145 * t146;
	t256 = t133 * t147;
	t259 = (-t132 * t253 - t144 * t256) / t137 ^ 2;
	t188 = t194 ^ 2;
	t245 = t183 * t188;
	t156 = t158 * t245 + 0.1e1;
	t233 = qJD(2) * t193;
	t220 = t190 * t235;
	t234 = qJD(2) * t191;
	t139 = ((t191 * t233 + t220) * t185 + t234 * t246) * t172;
	t247 = t171 * t190;
	t130 = (t139 * t191 - qJD(2)) * t247 + (t220 + (-t139 + t234) * t193) * t170;
	t257 = t130 * t157 * t158;
	t258 = (-t245 * t257 + (t188 * t190 * t233 - t211) * t158) / t156 ^ 2;
	t255 = t139 * t170;
	t254 = t139 * t190;
	t242 = t190 * t194;
	t164 = t205 * t242;
	t252 = t146 * t164;
	t251 = t158 * t190;
	t250 = t158 * t194;
	t162 = t172 * t265;
	t248 = t162 * t191;
	t236 = qJD(1) * t191;
	t227 = 0.2e1 * t259;
	t226 = -0.2e1 * t257;
	t225 = -0.2e1 * t147 * t207;
	t224 = t158 * t242;
	t222 = t172 * t183 * t185;
	t217 = -0.2e1 * t190 * t258;
	t216 = t133 * t225;
	t215 = t185 * t267;
	t212 = t191 * t222;
	t210 = t218 * t194;
	t167 = -t192 * t239 + t238;
	t208 = -t167 * t179 - t180 * t203;
	t149 = t167 * t180 - t179 * t203;
	t204 = t213 * t194;
	t163 = t206 * t242;
	t154 = 0.1e1 / t156;
	t143 = t192 * t204 + (qJD(2) * t190 * t192 + t214 * t189) * t191;
	t142 = -t214 * t240 + (t190 * t234 + t204) * t189;
	t138 = (-t260 * t190 * t170 + t171 * t212) * t194;
	t135 = 0.1e1 / t137;
	t134 = t170 * t239 - t247 + (-t170 * t193 + t171 * t241) * t162;
	t131 = t265 * t267 + (qJD(1) * t210 + 0.2e1 * t191 * t202) * t172;
	t1 = [t215 * t242 + (qJD(2) * t210 - t236 * t243) * t172, t131, 0, 0, 0, 0; (t157 * t217 + (t157 * t233 + (-qJD(1) * t138 - t130) * t251) * t154) * t191 + (t158 * t217 * t138 + (((-t139 * t212 - t260 * t233 + t249 * t261) * t170 + (t215 * t244 + t254 + (-t254 + (t261 + t266) * t234) * t172) * t171) * t224 + (t158 * t233 + t190 * t226) * t138 + (t157 + ((-t184 + t188) * t171 * t222 + t260 * t223) * t158) * t190 * qJD(1)) * t154) * t194, 0.2e1 * (-t134 * t251 + t157 * t193) * t194 * t258 + ((t157 * t236 + (qJD(2) * t134 + t130) * t250) * t193 + (t157 * t232 + (t131 * t171 * t191 - t170 * t234 - t248 * t255 + t255 + (qJD(2) * t170 + t171 * t235) * t162) * t224 + (-t158 * t236 + t194 * t226) * t134 + ((-t131 + t235) * t170 + ((-0.1e1 + t248) * qJD(2) + (-t162 + t191) * t139) * t171) * t193 * t250) * t190) * t154, 0, 0, 0, 0; (t145 * t208 - t149 * t253) * t227 + ((t149 * qJD(6) - t142 * t180 + t143 * t179) * t145 + t149 * t216 + (t208 * t133 + (t208 * qJD(6) + t142 * t179 + t143 * t180) * t207 - t149 * t132) * t146) * t135, (t145 * t163 + t207 * t252) * t227 + (t132 * t252 + (t163 * t146 - t164 * t225) * t133 - t262 * t193 * t232 + (t262 * t236 + ((-t264 * t145 + t263 * t253) * t180 + (t263 * t145 + t264 * t253) * t179) * t194) * t190) * t135, (t145 * t153 + t207 * t253) * t227 + (-t133 * t145 - t207 * t216 + (0.2e1 * t207 * t132 + t153 * t133) * t146) * t135, 0, 0, -0.2e1 * t259 - 0.2e1 * (t132 * t146 * t135 - (-t135 * t256 - t146 * t259) * t207) * t207;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end