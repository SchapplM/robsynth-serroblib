% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRP3
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
%   Wie in S6RRRPRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:38
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRP3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:40
	% EndTime: 2019-10-10 11:38:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:40
	% EndTime: 2019-10-10 11:38:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:40
	% EndTime: 2019-10-10 11:38:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:40
	% EndTime: 2019-10-10 11:38:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:41
	% EndTime: 2019-10-10 11:38:42
	% DurationCPUTime: 0.99s
	% Computational Cost: add. (3413->87), mult. (3414->191), div. (723->12), fcn. (4027->9), ass. (0->91)
	t151 = sin(qJ(1));
	t206 = 0.2e1 * t151;
	t148 = qJ(2) + qJ(3);
	t144 = cos(t148);
	t150 = cos(pkin(10));
	t180 = t151 * t150;
	t149 = sin(pkin(10));
	t152 = cos(qJ(1));
	t184 = t149 * t152;
	t130 = t144 * t184 - t180;
	t124 = t130 ^ 2;
	t181 = t151 * t149;
	t183 = t150 * t152;
	t131 = t144 * t183 + t181;
	t126 = 0.1e1 / t131 ^ 2;
	t121 = t124 * t126 + 0.1e1;
	t128 = -t144 * t181 - t183;
	t143 = sin(t148);
	t145 = qJD(2) + qJD(3);
	t185 = t145 * t152;
	t170 = t143 * t185;
	t122 = t128 * qJD(1) - t149 * t170;
	t196 = t126 * t130;
	t129 = -t144 * t180 + t184;
	t123 = t129 * qJD(1) - t150 * t170;
	t125 = 0.1e1 / t131;
	t198 = t123 * t125 * t126;
	t205 = (t122 * t196 - t124 * t198) / t121 ^ 2;
	t146 = t151 ^ 2;
	t139 = t143 ^ 2;
	t141 = 0.1e1 / t144 ^ 2;
	t192 = t139 * t141;
	t137 = t146 * t192 + 0.1e1;
	t135 = 0.1e1 / t137;
	t140 = 0.1e1 / t144;
	t178 = qJD(1) * t152;
	t169 = t143 * t178;
	t186 = t145 * t151;
	t172 = t141 * t186;
	t109 = (-(-t144 * t186 - t169) * t140 + t139 * t172) * t135;
	t204 = t109 - t186;
	t182 = t151 * t143;
	t134 = atan2(-t182, -t144);
	t133 = cos(t134);
	t132 = sin(t134);
	t173 = t132 * t182;
	t117 = -t133 * t144 - t173;
	t114 = 0.1e1 / t117;
	t115 = 0.1e1 / t117 ^ 2;
	t203 = t135 - 0.1e1;
	t194 = t133 * t143;
	t104 = (-t109 * t151 + t145) * t194 + (t204 * t144 - t169) * t132;
	t202 = t104 * t114 * t115;
	t138 = t143 * t139;
	t189 = t140 * t143;
	t160 = t145 * (t138 * t140 * t141 + t189);
	t190 = t139 * t151;
	t163 = t178 * t190;
	t201 = (t141 * t163 + t146 * t160) / t137 ^ 2;
	t200 = t115 * t143;
	t199 = t115 * t152;
	t197 = t125 * t149;
	t195 = t132 * t151;
	t193 = t139 * t140;
	t147 = t152 ^ 2;
	t191 = t139 * t147;
	t188 = t143 * t152;
	t187 = t144 * t145;
	t179 = qJD(1) * t151;
	t112 = t115 * t191 + 0.1e1;
	t177 = 0.2e1 * (-t191 * t202 + (t143 * t147 * t187 - t163) * t115) / t112 ^ 2;
	t176 = 0.2e1 * t202;
	t175 = t115 * t188;
	t174 = t130 * t198;
	t171 = t145 * t182;
	t168 = 0.1e1 + t192;
	t167 = t143 * t177;
	t166 = -0.2e1 * t143 * t201;
	t165 = t201 * t206;
	t164 = t133 * t135 * t193;
	t162 = t168 * t152;
	t161 = t150 * t196 - t197;
	t119 = 0.1e1 / t121;
	t118 = t168 * t151 * t135;
	t110 = 0.1e1 / t112;
	t108 = (t203 * t143 * t132 - t151 * t164) * t152;
	t106 = -t144 * t195 + t194 + (t132 * t144 - t133 * t182) * t118;
	t105 = -t168 * t165 + (qJD(1) * t162 + t160 * t206) * t135;
	t102 = -0.2e1 * t161 * t188 * t205 + (t161 * t144 * t185 + (-0.2e1 * t174 * t183 + t179 * t197 + (t123 * t184 + (t122 * t152 - t130 * t179) * t150) * t126) * t143) * t119;
	t101 = (t106 * t200 - t114 * t144) * t152 * t177 + ((-t114 * t179 + (-t106 * t145 - t104) * t199) * t144 + (-t114 * t185 - (-t105 * t133 * t151 - t204 * t132 + (t109 * t195 - t132 * t145 - t133 * t178) * t118) * t175 + (t115 * t179 + t152 * t176) * t106 - ((t105 - t178) * t132 + ((-t118 * t151 + 0.1e1) * t145 + (t118 - t151) * t109) * t133) * t144 * t199) * t143) * t110;
	t1 = [t140 * t152 * t166 + (t145 * t162 - t179 * t189) * t135, t105, t105, 0, 0, 0; (t114 * t167 + (-t114 * t187 + (qJD(1) * t108 + t104) * t200) * t110) * t151 + (t115 * t167 * t108 + (-((t166 - t187 + (t109 * t140 * t190 + t187) * t135) * t132 + (t165 * t193 - t109 * t143 + (-t138 * t172 + (t109 - 0.2e1 * t186) * t143) * t135) * t133) * t175 + (-t115 * t187 + t143 * t176) * t108 + (-t114 + ((-t146 + t147) * t164 + t203 * t173) * t115) * t143 * qJD(1)) * t110) * t152, t101, t101, 0, 0, 0; 0.2e1 * (-t125 * t128 + t129 * t196) * t205 + ((-t130 * qJD(1) + t149 * t171) * t125 + 0.2e1 * t129 * t174 + (-t128 * t123 - (-t131 * qJD(1) + t150 * t171) * t130 - t129 * t122) * t126) * t119, t102, t102, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:41
	% EndTime: 2019-10-10 11:38:42
	% DurationCPUTime: 1.14s
	% Computational Cost: add. (4131->97), mult. (3810->206), div. (753->12), fcn. (4455->9), ass. (0->95)
	t170 = sin(qJ(1));
	t230 = 0.2e1 * t170;
	t167 = t170 ^ 2;
	t169 = qJ(2) + qJ(3);
	t163 = sin(t169);
	t157 = t163 ^ 2;
	t164 = cos(t169);
	t159 = 0.1e1 / t164 ^ 2;
	t216 = t157 * t159;
	t153 = t167 * t216 + 0.1e1;
	t150 = 0.1e1 / t153;
	t158 = 0.1e1 / t164;
	t171 = cos(qJ(1));
	t202 = qJD(1) * t171;
	t192 = t163 * t202;
	t166 = qJD(2) + qJD(3);
	t208 = t166 * t170;
	t195 = t159 * t208;
	t124 = (-(-t164 * t208 - t192) * t158 + t157 * t195) * t150;
	t229 = t124 - t208;
	t165 = pkin(10) + qJ(5);
	t162 = cos(t165);
	t161 = sin(t165);
	t206 = t170 * t161;
	t209 = t164 * t171;
	t146 = t162 * t209 + t206;
	t204 = t170 * t163;
	t149 = atan2(-t204, -t164);
	t148 = cos(t149);
	t147 = sin(t149);
	t196 = t147 * t204;
	t134 = -t148 * t164 - t196;
	t131 = 0.1e1 / t134;
	t140 = 0.1e1 / t146;
	t132 = 0.1e1 / t134 ^ 2;
	t141 = 0.1e1 / t146 ^ 2;
	t228 = t150 - 0.1e1;
	t218 = t148 * t163;
	t119 = (-t124 * t170 + t166) * t218 + (t229 * t164 - t192) * t147;
	t227 = t119 * t131 * t132;
	t181 = t162 * t171 + t164 * t206;
	t207 = t166 * t171;
	t193 = t163 * t207;
	t125 = t181 * qJD(1) - t146 * qJD(5) + t161 * t193;
	t205 = t170 * t162;
	t145 = t161 * t209 - t205;
	t139 = t145 ^ 2;
	t137 = t139 * t141 + 0.1e1;
	t221 = t141 * t145;
	t186 = -qJD(1) * t164 + qJD(5);
	t187 = qJD(5) * t164 - qJD(1);
	t212 = t161 * t171;
	t126 = -t187 * t212 + (t186 * t170 - t193) * t162;
	t225 = t126 * t140 * t141;
	t226 = (-t125 * t221 - t139 * t225) / t137 ^ 2;
	t156 = t163 * t157;
	t213 = t158 * t163;
	t180 = t166 * (t156 * t158 * t159 + t213);
	t214 = t157 * t170;
	t184 = t202 * t214;
	t224 = (t159 * t184 + t167 * t180) / t153 ^ 2;
	t223 = t132 * t163;
	t222 = t140 * t161;
	t220 = t145 * t162;
	t219 = t147 * t170;
	t217 = t157 * t158;
	t168 = t171 ^ 2;
	t215 = t157 * t168;
	t211 = t163 * t171;
	t210 = t164 * t166;
	t203 = qJD(1) * t170;
	t129 = t132 * t215 + 0.1e1;
	t201 = 0.2e1 * (-t215 * t227 + (t163 * t168 * t210 - t184) * t132) / t129 ^ 2;
	t200 = 0.2e1 * t227;
	t199 = -0.2e1 * t226;
	t198 = t132 * t211;
	t197 = t145 * t225;
	t191 = 0.1e1 + t216;
	t190 = t163 * t201;
	t189 = -0.2e1 * t163 * t224;
	t188 = t224 * t230;
	t185 = t148 * t150 * t217;
	t183 = t191 * t171;
	t182 = t141 * t220 - t222;
	t179 = t166 * t204 + t186 * t171;
	t144 = -t164 * t205 + t212;
	t138 = t191 * t170 * t150;
	t135 = 0.1e1 / t137;
	t127 = 0.1e1 / t129;
	t123 = (t228 * t163 * t147 - t170 * t185) * t171;
	t122 = -t164 * t219 + t218 + (t147 * t164 - t148 * t204) * t138;
	t120 = -t191 * t188 + (qJD(1) * t183 + t180 * t230) * t150;
	t117 = t182 * t199 * t211 + (t182 * t164 * t207 + (-t182 * t203 + ((-qJD(5) * t140 - 0.2e1 * t197) * t162 + (-t125 * t162 + (-qJD(5) * t145 + t126) * t161) * t141) * t171) * t163) * t135;
	t116 = (t122 * t223 - t131 * t164) * t171 * t201 + ((-t131 * t203 + (-t122 * t166 - t119) * t171 * t132) * t164 + (-t131 * t207 - (-t120 * t148 * t170 - t229 * t147 + (t124 * t219 - t147 * t166 - t148 * t202) * t138) * t198 + (t132 * t203 + t171 * t200) * t122 - ((t120 - t202) * t147 + ((-t138 * t170 + 0.1e1) * t166 + (t138 - t170) * t124) * t148) * t132 * t209) * t163) * t127;
	t1 = [t171 * t158 * t189 + (t166 * t183 - t203 * t213) * t150, t120, t120, 0, 0, 0; (t131 * t190 + (-t131 * t210 + (qJD(1) * t123 + t119) * t223) * t127) * t170 + (t132 * t190 * t123 + (-((t189 - t210 + (t124 * t158 * t214 + t210) * t150) * t147 + (t188 * t217 - t124 * t163 + (-t156 * t195 + (t124 - 0.2e1 * t208) * t163) * t150) * t148) * t198 + (-t132 * t210 + t163 * t200) * t123 + (-t131 + ((-t167 + t168) * t185 + t228 * t196) * t132) * t163 * qJD(1)) * t127) * t171, t116, t116, 0, 0, 0; 0.2e1 * (t140 * t181 + t144 * t221) * t226 + (0.2e1 * t144 * t197 - t187 * t140 * t205 + t179 * t222 + (-t187 * t145 * t206 + t144 * t125 + t126 * t181 - t179 * t220) * t141) * t135, t117, t117, 0, t199 + 0.2e1 * (-t125 * t141 * t135 + (-t135 * t225 - t141 * t226) * t145) * t145, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:38:41
	% EndTime: 2019-10-10 11:38:42
	% DurationCPUTime: 1.71s
	% Computational Cost: add. (9575->126), mult. (8382->274), div. (1515->15), fcn. (10508->9), ass. (0->119)
	t201 = qJ(2) + qJ(3);
	t197 = cos(t201);
	t198 = pkin(10) + qJ(5);
	t194 = sin(t198);
	t274 = sin(qJ(1));
	t231 = t274 * t194;
	t195 = cos(t198);
	t202 = cos(qJ(1));
	t252 = t202 * t195;
	t176 = t197 * t252 + t231;
	t170 = 0.1e1 / t176 ^ 2;
	t196 = sin(t201);
	t190 = t196 ^ 2;
	t200 = t202 ^ 2;
	t259 = t190 * t200;
	t239 = t170 * t259;
	t166 = 0.1e1 + t239;
	t226 = qJD(1) * t274;
	t199 = qJD(2) + qJD(3);
	t255 = t199 * t202;
	t233 = t196 * t255;
	t212 = t197 * t226 + t233;
	t225 = t274 * qJD(5);
	t253 = t202 * t194;
	t155 = (-qJD(5) * t197 + qJD(1)) * t253 + (t225 - t212) * t195;
	t169 = 0.1e1 / t176;
	t269 = t155 * t169 * t170;
	t220 = t259 * t269;
	t234 = t196 * t199 * t200;
	t277 = (-t220 + (-t190 * t202 * t226 + t197 * t234) * t170) / t166 ^ 2;
	t257 = t196 * t202;
	t172 = t197 * t231 + t252;
	t217 = t194 * t225;
	t248 = qJD(5) * t202;
	t228 = t195 * t248;
	t154 = t172 * qJD(1) + t194 * t233 - t197 * t228 - t217;
	t230 = t274 * t195;
	t175 = t197 * t253 - t230;
	t187 = 0.1e1 / t194;
	t188 = 0.1e1 / t194 ^ 2;
	t191 = 0.1e1 / t196;
	t192 = 0.1e1 / t196 ^ 2;
	t256 = t197 * t199;
	t235 = t192 * t256;
	t250 = qJD(5) * t195;
	t262 = t187 * t191;
	t276 = (t188 * t191 * t250 + t187 * t235) * t175 + t154 * t262;
	t258 = t196 * t194;
	t162 = atan2(-t172, t258);
	t159 = cos(t162);
	t158 = sin(t162);
	t268 = t158 * t172;
	t153 = t159 * t258 - t268;
	t150 = 0.1e1 / t153;
	t151 = 0.1e1 / t153 ^ 2;
	t275 = 0.2e1 * t175;
	t167 = t172 ^ 2;
	t261 = t188 * t192;
	t163 = t167 * t261 + 0.1e1;
	t160 = 0.1e1 / t163;
	t249 = qJD(5) * t196;
	t213 = t194 * t256 + t195 * t249;
	t237 = t172 * t261;
	t218 = t195 * t226;
	t232 = t196 * t274;
	t219 = t199 * t232;
	t251 = qJD(1) * t202;
	t156 = t195 * t225 * t197 - t218 + (t251 * t197 - t219 - t248) * t194;
	t240 = t156 * t262;
	t142 = (t213 * t237 - t240) * t160;
	t210 = -t142 * t172 + t213;
	t138 = (-t142 * t258 - t156) * t158 + t210 * t159;
	t152 = t150 * t151;
	t273 = t138 * t152;
	t189 = t187 * t188;
	t193 = t191 / t190;
	t229 = t192 * t250;
	t272 = (t156 * t237 + (-t188 * t193 * t256 - t189 * t229) * t167) / t163 ^ 2;
	t271 = t151 * t175;
	t270 = t154 * t151;
	t267 = t158 * t175;
	t266 = t158 * t196;
	t265 = t159 * t172;
	t264 = t159 * t175;
	t263 = t159 * t197;
	t260 = t188 * t195;
	t254 = t202 * t150;
	t168 = t175 ^ 2;
	t148 = t151 * t168 + 0.1e1;
	t247 = 0.2e1 * (-t168 * t273 - t175 * t270) / t148 ^ 2;
	t246 = -0.2e1 * t272;
	t245 = 0.2e1 * t277;
	t244 = t152 * t275;
	t243 = t191 * t272;
	t242 = t151 * t267;
	t238 = t172 * t262;
	t236 = t187 * t192 * t197;
	t215 = t172 * t236 + t274;
	t149 = t215 * t160;
	t227 = t274 - t149;
	t224 = t150 * t247;
	t223 = t151 * t247;
	t222 = t257 * t275;
	t221 = t187 * t243;
	t174 = t197 * t230 - t253;
	t216 = t172 * t260 - t174 * t187;
	t214 = t170 * t174 * t202 - t274 * t169;
	t164 = 0.1e1 / t166;
	t157 = t176 * qJD(1) - t195 * t219 - t197 * t217 - t228;
	t146 = 0.1e1 / t148;
	t145 = t216 * t191 * t160;
	t141 = (-t158 + (t159 * t238 + t158) * t160) * t175;
	t140 = -t149 * t265 + (t227 * t266 + t263) * t194;
	t139 = t159 * t195 * t196 - t158 * t174 + (-t158 * t258 - t265) * t145;
	t137 = t215 * t246 + (t156 * t236 + t251 + (-t188 * t197 * t229 + (-0.2e1 * t193 * t197 ^ 2 - t191) * t199 * t187) * t172) * t160;
	t135 = (t169 * t197 * t202 + t195 * t239) * t245 + (0.2e1 * t195 * t220 + t212 * t169 + ((t155 * t202 - 0.2e1 * t195 * t234) * t197 + (qJD(5) * t194 * t200 + 0.2e1 * t202 * t218) * t190) * t170) * t164;
	t134 = -0.2e1 * t216 * t243 + (-t216 * t235 + (t156 * t260 - t157 * t187 + (t174 * t260 + (-0.2e1 * t189 * t195 ^ 2 - t187) * t172) * qJD(5)) * t191) * t160;
	t133 = t140 * t175 * t223 + (-(-t137 * t265 + (t142 * t268 - t156 * t159) * t149) * t271 + (t138 * t244 + t270) * t140 + (-t196 * t254 - (-t149 * t266 + t158 * t232 + t263) * t271) * t250) * t146 + (t224 * t257 + ((-t199 * t254 - (t227 * t199 - t142) * t242) * t197 + (t150 * t226 + (t202 * t138 - (-t137 + t251) * t267 - (t227 * t142 - t199) * t264) * t151) * t196) * t146) * t194;
	t1 = [t276 * t160 + t221 * t275, t137, t137, 0, t134, 0; t172 * t224 + (-t156 * t150 + (t138 * t172 + t141 * t154) * t151) * t146 + (t141 * t223 + (0.2e1 * t141 * t273 + (t154 * t160 - t154 - (-t142 * t160 * t238 + t246) * t175) * t151 * t158 + (-(-0.2e1 * t172 * t221 - t142) * t271 + (-(t142 + t240) * t175 + t276 * t172) * t151 * t160) * t159) * t146) * t175, t133, t133, 0, (t139 * t271 - t150 * t176) * t247 + (t139 * t270 + t155 * t150 + (t139 * t244 - t151 * t176) * t138 - (-t194 * t249 + t195 * t256 - t134 * t172 - t145 * t156 + (-t145 * t258 - t174) * t142) * t151 * t264 - (-t157 + (-t134 * t194 - t142 * t195) * t196 - t210 * t145) * t242) * t146, 0; t214 * t196 * t245 + (-t214 * t256 + ((qJD(1) * t169 + 0.2e1 * t174 * t269) * t202 + (-t274 * t155 - t157 * t202 + t174 * t226) * t170) * t196) * t164, t135, t135, 0, t170 * t222 * t277 + (t222 * t269 + (t154 * t257 + (t196 * t226 - t197 * t255) * t175) * t170) * t164, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end