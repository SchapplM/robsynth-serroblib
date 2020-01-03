% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRP9
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
%   Wie in S5RRRRP9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRRP9_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP9_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:07:14
	% EndTime: 2019-12-31 22:07:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:07:14
	% EndTime: 2019-12-31 22:07:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:07:14
	% EndTime: 2019-12-31 22:07:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:07:14
	% EndTime: 2019-12-31 22:07:15
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
	t86 = t138 * qJD(1) - t107 * qJD(3) + t124 * t148;
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
	% StartTime: 2019-12-31 22:07:14
	% EndTime: 2019-12-31 22:07:15
	% DurationCPUTime: 0.76s
	% Computational Cost: add. (1570->96), mult. (2734->204), div. (498->12), fcn. (3199->9), ass. (0->95)
	t152 = sin(qJ(2));
	t145 = t152 ^ 2;
	t154 = cos(qJ(2));
	t148 = 0.1e1 / t154 ^ 2;
	t199 = t145 * t148;
	t153 = sin(qJ(1));
	t217 = 0.2e1 * t153;
	t216 = t152 * t199;
	t151 = qJ(3) + qJ(4);
	t142 = cos(t151);
	t155 = cos(qJ(1));
	t191 = t154 * t155;
	t141 = sin(t151);
	t195 = t153 * t141;
	t131 = t142 * t191 + t195;
	t193 = t153 * t152;
	t135 = atan2(-t193, -t154);
	t134 = cos(t135);
	t133 = sin(t135);
	t180 = t133 * t193;
	t121 = -t134 * t154 - t180;
	t118 = 0.1e1 / t121;
	t125 = 0.1e1 / t131;
	t147 = 0.1e1 / t154;
	t119 = 0.1e1 / t121 ^ 2;
	t126 = 0.1e1 / t131 ^ 2;
	t215 = -0.2e1 * t152;
	t146 = t153 ^ 2;
	t139 = t146 * t199 + 0.1e1;
	t137 = 0.1e1 / t139;
	t214 = t137 - 0.1e1;
	t143 = qJD(3) + qJD(4);
	t192 = t153 * t154;
	t165 = t141 * t192 + t142 * t155;
	t186 = qJD(2) * t155;
	t176 = t152 * t186;
	t109 = t165 * qJD(1) - t131 * t143 + t141 * t176;
	t194 = t153 * t142;
	t130 = t141 * t191 - t194;
	t124 = t130 ^ 2;
	t117 = t124 * t126 + 0.1e1;
	t204 = t126 * t130;
	t170 = -qJD(1) * t154 + t143;
	t171 = t143 * t154 - qJD(1);
	t201 = t141 * t155;
	t110 = -t171 * t201 + (t170 * t153 - t176) * t142;
	t211 = t110 * t125 * t126;
	t213 = (-t109 * t204 - t124 * t211) / t117 ^ 2;
	t189 = qJD(1) * t155;
	t177 = t152 * t189;
	t187 = qJD(2) * t154;
	t188 = qJD(2) * t153;
	t111 = (-(-t153 * t187 - t177) * t147 + t188 * t199) * t137;
	t202 = t134 * t152;
	t105 = (-t111 * t153 + qJD(2)) * t202 + (-t177 + (t111 - t188) * t154) * t133;
	t212 = t105 * t118 * t119;
	t210 = t111 * t133;
	t209 = t111 * t152;
	t208 = t119 * t152;
	t197 = t147 * t152;
	t164 = qJD(2) * (t147 * t216 + t197);
	t168 = t145 * t153 * t189;
	t207 = (t146 * t164 + t148 * t168) / t139 ^ 2;
	t175 = 0.1e1 + t199;
	t123 = t175 * t153 * t137;
	t206 = t123 * t153;
	t205 = t125 * t141;
	t203 = t130 * t142;
	t200 = t145 * t147;
	t150 = t155 ^ 2;
	t198 = t145 * t150;
	t196 = t152 * t155;
	t190 = qJD(1) * t153;
	t114 = t119 * t198 + 0.1e1;
	t185 = 0.2e1 * (-t198 * t212 + (t150 * t152 * t187 - t168) * t119) / t114 ^ 2;
	t184 = -0.2e1 * t213;
	t183 = 0.2e1 * t212;
	t182 = t130 * t211;
	t181 = t119 * t196;
	t179 = t137 * t200;
	t174 = t152 * t185;
	t173 = t207 * t215;
	t172 = t207 * t217;
	t169 = t153 * t179;
	t167 = t175 * t155;
	t166 = t126 * t203 - t205;
	t163 = t152 * t188 + t170 * t155;
	t129 = -t142 * t192 + t201;
	t115 = 0.1e1 / t117;
	t112 = 0.1e1 / t114;
	t108 = (t214 * t152 * t133 - t134 * t169) * t155;
	t107 = -t133 * t192 + t202 + (t133 * t154 - t134 * t193) * t123;
	t106 = -t175 * t172 + (qJD(1) * t167 + t164 * t217) * t137;
	t102 = t184 + 0.2e1 * (-t109 * t115 * t126 + (-t115 * t211 - t126 * t213) * t130) * t130;
	t1 = [t147 * t155 * t173 + (qJD(2) * t167 - t190 * t197) * t137, t106, 0, 0, 0; (t118 * t174 + (-t118 * t187 + (qJD(1) * t108 + t105) * t208) * t112) * t153 + (t119 * t174 * t108 + (-((t111 * t169 + t214 * t187 + t173) * t133 + (t172 * t200 - t209 + (t209 + (t215 - t216) * t188) * t137) * t134) * t181 + (-t119 * t187 + t152 * t183) * t108 + (-t118 + ((-t146 + t150) * t134 * t179 + t214 * t180) * t119) * t152 * qJD(1)) * t112) * t155, (t107 * t208 - t118 * t154) * t155 * t185 + ((-t118 * t190 + (-qJD(2) * t107 - t105) * t155 * t119) * t154 + (-t118 * t186 - (-t106 * t134 * t153 + t133 * t188 + t206 * t210 - t210 + (-qJD(2) * t133 - t134 * t189) * t123) * t181 + (t119 * t190 + t155 * t183) * t107 - ((t106 - t189) * t133 + ((0.1e1 - t206) * qJD(2) + (t123 - t153) * t111) * t134) * t119 * t191) * t152) * t112, 0, 0, 0; 0.2e1 * (t125 * t165 + t129 * t204) * t213 + (0.2e1 * t129 * t182 - t171 * t125 * t194 + t163 * t205 + (-t171 * t130 * t195 + t129 * t109 + t110 * t165 - t163 * t203) * t126) * t115, t166 * t184 * t196 + (t166 * t154 * t186 + (-t166 * t190 + ((-t125 * t143 - 0.2e1 * t182) * t142 + (-t109 * t142 + (-t130 * t143 + t110) * t141) * t126) * t155) * t152) * t115, t102, t102, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 22:07:14
	% EndTime: 2019-12-31 22:07:15
	% DurationCPUTime: 1.18s
	% Computational Cost: add. (6985->123), mult. (8378->264), div. (1558->15), fcn. (10537->9), ass. (0->116)
	t191 = qJ(3) + qJ(4);
	t184 = cos(t191);
	t268 = 0.2e1 * t184;
	t183 = sin(t191);
	t193 = cos(qJ(2));
	t194 = cos(qJ(1));
	t241 = t194 * t184;
	t224 = t193 * t241;
	t262 = sin(qJ(1));
	t169 = t262 * t183 + t224;
	t163 = 0.1e1 / t169 ^ 2;
	t192 = sin(qJ(2));
	t186 = t192 ^ 2;
	t190 = t194 ^ 2;
	t246 = t186 * t190;
	t225 = t163 * t246;
	t159 = 0.1e1 + t225;
	t217 = qJD(1) * t262;
	t239 = qJD(2) * t193;
	t202 = t186 * t194 * t217 - t190 * t192 * t239;
	t185 = qJD(3) + qJD(4);
	t238 = qJD(2) * t194;
	t219 = t192 * t238;
	t205 = t193 * t217 + t219;
	t223 = t262 * t185;
	t242 = t194 * t183;
	t148 = (-t185 * t193 + qJD(1)) * t242 + (t223 - t205) * t184;
	t162 = 0.1e1 / t169;
	t257 = t148 * t162 * t163;
	t212 = t246 * t257;
	t267 = (-t202 * t163 - t212) / t159 ^ 2;
	t221 = t262 * t193;
	t165 = t183 * t221 + t241;
	t266 = t165 * t185;
	t244 = t192 * t194;
	t147 = t165 * qJD(1) - t185 * t224 + (t219 - t223) * t183;
	t168 = -t262 * t184 + t193 * t242;
	t180 = 0.1e1 / t183;
	t181 = 0.1e1 / t183 ^ 2;
	t187 = 0.1e1 / t192;
	t188 = 0.1e1 / t192 ^ 2;
	t220 = t188 * t239;
	t248 = t184 * t185;
	t250 = t180 * t187;
	t265 = (t181 * t187 * t248 + t180 * t220) * t168 + t147 * t250;
	t245 = t192 * t183;
	t155 = atan2(-t165, t245);
	t152 = cos(t155);
	t151 = sin(t155);
	t256 = t151 * t165;
	t146 = t152 * t245 - t256;
	t143 = 0.1e1 / t146;
	t144 = 0.1e1 / t146 ^ 2;
	t264 = -0.2e1 * t165;
	t263 = 0.2e1 * t168;
	t160 = t165 ^ 2;
	t249 = t181 * t188;
	t156 = t160 * t249 + 0.1e1;
	t153 = 0.1e1 / t156;
	t247 = t184 * t192;
	t206 = t183 * t239 + t185 * t247;
	t228 = t165 * t249;
	t222 = t262 * t192;
	t210 = qJD(2) * t222;
	t240 = qJD(1) * t194;
	t149 = -t183 * t210 - t185 * t242 - t184 * t217 + (t183 * t240 + t184 * t223) * t193;
	t230 = t149 * t250;
	t135 = (t206 * t228 - t230) * t153;
	t203 = -t135 * t165 + t206;
	t130 = (-t135 * t245 - t149) * t151 + t203 * t152;
	t145 = t143 * t144;
	t261 = t130 * t145;
	t182 = t180 * t181;
	t189 = t187 / t186;
	t226 = t188 * t248;
	t260 = (t149 * t228 + (-t181 * t189 * t239 - t182 * t226) * t160) / t156 ^ 2;
	t259 = t144 * t168;
	t258 = t147 * t144;
	t255 = t151 * t168;
	t254 = t151 * t192;
	t253 = t152 * t165;
	t252 = t152 * t168;
	t251 = t152 * t193;
	t243 = t193 * t194;
	t161 = t168 ^ 2;
	t141 = t161 * t144 + 0.1e1;
	t237 = 0.2e1 * (-t161 * t261 - t168 * t258) / t141 ^ 2;
	t236 = -0.2e1 * t260;
	t235 = 0.2e1 * t267;
	t234 = t145 * t263;
	t233 = t187 * t260;
	t232 = t144 * t255;
	t229 = t165 * t250;
	t227 = t180 * t188 * t193;
	t208 = t165 * t227 + t262;
	t142 = t208 * t153;
	t218 = t262 - t142;
	t216 = t143 * t237;
	t215 = t144 * t237;
	t214 = t244 * t263;
	t213 = t180 * t233;
	t167 = t184 * t221 - t242;
	t209 = t165 * t181 * t184 - t167 * t180;
	t207 = t163 * t167 * t194 - t262 * t162;
	t157 = 0.1e1 / t159;
	t150 = t169 * qJD(1) - t184 * t210 - t266;
	t139 = 0.1e1 / t141;
	t138 = t209 * t187 * t153;
	t134 = (-t151 + (t152 * t229 + t151) * t153) * t168;
	t133 = -t142 * t253 + (t218 * t254 + t251) * t183;
	t132 = t152 * t247 - t151 * t167 + (-t151 * t245 - t253) * t138;
	t131 = t163 * t214 * t267 + (t214 * t257 + (t147 * t244 + (t192 * t217 - t193 * t238) * t168) * t163) * t157;
	t129 = t208 * t236 + (t149 * t227 + t240 + (-t181 * t193 * t226 + (-0.2e1 * t189 * t193 ^ 2 - t187) * t180 * qJD(2)) * t165) * t153;
	t127 = -0.2e1 * t209 * t233 + (-t209 * t220 + ((-t150 - t266) * t180 + (t182 * t248 * t264 + (t167 * t185 + t149) * t181) * t184) * t187) * t153;
	t126 = (t132 * t259 - t143 * t169) * t237 + (t132 * t258 + t148 * t143 + (t132 * t234 - t169 * t144) * t130 - (t184 * t239 - t185 * t245 - t127 * t165 - t138 * t149 + (-t138 * t245 - t167) * t135) * t144 * t252 - (-t150 + (-t127 * t183 - t135 * t184) * t192 - t203 * t138) * t232) * t139;
	t1 = [t265 * t153 + t213 * t263, t129, t127, t127, 0; t165 * t216 + (-t149 * t143 + (t130 * t165 + t134 * t147) * t144) * t139 + (t134 * t215 + (0.2e1 * t134 * t261 + (t147 * t153 - t147 - (-t135 * t153 * t229 + t236) * t168) * t144 * t151 + (-(t213 * t264 - t135) * t259 + (-(t135 + t230) * t168 + t265 * t165) * t144 * t153) * t152) * t139) * t168, t133 * t168 * t215 + (-(-t129 * t253 + (t135 * t256 - t149 * t152) * t142) * t259 + (-t143 * t244 - (-t142 * t254 + t151 * t222 + t251) * t259) * t248 + (t130 * t234 + t258) * t133) * t139 + (t216 * t244 + ((-t143 * t238 - (t218 * qJD(2) - t135) * t232) * t193 + (t143 * t217 + (t194 * t130 - (-t129 + t240) * t255 - (t218 * t135 - qJD(2)) * t252) * t144) * t192) * t139) * t183, t126, t126, 0; t207 * t192 * t235 + (-t207 * t239 + ((qJD(1) * t162 + 0.2e1 * t167 * t257) * t194 + (-t262 * t148 - t150 * t194 + t167 * t217) * t163) * t192) * t157, (t162 * t243 + t184 * t225) * t235 + (t212 * t268 + t205 * t162 + (t183 * t185 * t246 + t148 * t243 + t202 * t268) * t163) * t157, t131, t131, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end