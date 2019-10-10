% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRP2
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
%   Wie in S6RRRPRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:36
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRP2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:36:55
	% EndTime: 2019-10-10 11:36:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:36:55
	% EndTime: 2019-10-10 11:36:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:36:55
	% EndTime: 2019-10-10 11:36:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:36:55
	% EndTime: 2019-10-10 11:36:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:36:55
	% EndTime: 2019-10-10 11:36:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:36:55
	% EndTime: 2019-10-10 11:36:56
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (5587->97), mult. (3810->208), div. (753->12), fcn. (4455->9), ass. (0->95)
	t179 = sin(qJ(1));
	t241 = 0.2e1 * t179;
	t176 = t179 ^ 2;
	t174 = qJ(2) + qJ(3) + pkin(10);
	t171 = sin(t174);
	t167 = t171 ^ 2;
	t172 = cos(t174);
	t169 = 0.1e1 / t172 ^ 2;
	t226 = t167 * t169;
	t163 = t176 * t226 + 0.1e1;
	t161 = 0.1e1 / t163;
	t168 = 0.1e1 / t172;
	t181 = cos(qJ(1));
	t212 = qJD(1) * t181;
	t202 = t171 * t212;
	t175 = qJD(2) + qJD(3);
	t220 = t175 * t179;
	t205 = t169 * t220;
	t135 = (-(-t172 * t220 - t202) * t168 + t167 * t205) * t161;
	t240 = t135 - t220;
	t180 = cos(qJ(5));
	t214 = t181 * t180;
	t178 = sin(qJ(5));
	t217 = t179 * t178;
	t160 = t172 * t214 + t217;
	t218 = t179 * t171;
	t156 = atan2(-t218, -t172);
	t151 = cos(t156);
	t150 = sin(t156);
	t207 = t150 * t218;
	t145 = -t151 * t172 - t207;
	t142 = 0.1e1 / t145;
	t153 = 0.1e1 / t160;
	t143 = 0.1e1 / t145 ^ 2;
	t154 = 0.1e1 / t160 ^ 2;
	t239 = t161 - 0.1e1;
	t231 = t151 * t171;
	t130 = (-t135 * t179 + t175) * t231 + (t240 * t172 - t202) * t150;
	t238 = t130 * t142 * t143;
	t190 = t172 * t217 + t214;
	t219 = t175 * t181;
	t203 = t171 * t219;
	t140 = t190 * qJD(1) - t160 * qJD(5) + t178 * t203;
	t215 = t181 * t178;
	t216 = t179 * t180;
	t159 = t172 * t215 - t216;
	t152 = t159 ^ 2;
	t149 = t152 * t154 + 0.1e1;
	t229 = t154 * t159;
	t196 = -qJD(1) * t172 + qJD(5);
	t197 = qJD(5) * t172 - qJD(1);
	t141 = -t197 * t215 + (t196 * t179 - t203) * t180;
	t235 = t141 * t153 * t154;
	t237 = (-t140 * t229 - t152 * t235) / t149 ^ 2;
	t166 = t171 * t167;
	t223 = t168 * t171;
	t189 = t175 * (t166 * t168 * t169 + t223);
	t224 = t167 * t179;
	t194 = t212 * t224;
	t236 = (t169 * t194 + t176 * t189) / t163 ^ 2;
	t234 = t143 * t171;
	t233 = t143 * t181;
	t232 = t150 * t179;
	t230 = t153 * t178;
	t228 = t159 * t180;
	t227 = t167 * t168;
	t177 = t181 ^ 2;
	t225 = t167 * t177;
	t222 = t171 * t181;
	t221 = t172 * t175;
	t213 = qJD(1) * t179;
	t138 = t143 * t225 + 0.1e1;
	t211 = 0.2e1 * (-t225 * t238 + (t171 * t177 * t221 - t194) * t143) / t138 ^ 2;
	t210 = 0.2e1 * t238;
	t209 = -0.2e1 * t237;
	t208 = t143 * t222;
	t206 = t159 * t235;
	t201 = 0.1e1 + t226;
	t200 = t171 * t211;
	t199 = -0.2e1 * t171 * t236;
	t198 = t236 * t241;
	t195 = t151 * t161 * t227;
	t193 = t201 * t181;
	t192 = t196 * t181;
	t191 = t154 * t228 - t230;
	t158 = -t172 * t216 + t215;
	t147 = 0.1e1 / t149;
	t146 = t201 * t179 * t161;
	t136 = 0.1e1 / t138;
	t134 = (t239 * t171 * t150 - t179 * t195) * t181;
	t132 = -t172 * t232 + t231 + (t150 * t172 - t151 * t218) * t146;
	t131 = -t201 * t198 + (qJD(1) * t193 + t189 * t241) * t161;
	t128 = t191 * t209 * t222 + (t191 * t172 * t219 + (-t191 * t213 + ((-qJD(5) * t153 - 0.2e1 * t206) * t180 + (-t140 * t180 + (-qJD(5) * t159 + t141) * t178) * t154) * t181) * t171) * t147;
	t127 = (t132 * t234 - t142 * t172) * t181 * t211 + ((-t142 * t213 + (-t132 * t175 - t130) * t233) * t172 + (-t142 * t219 - (-t131 * t151 * t179 - t240 * t150 + (t135 * t232 - t150 * t175 - t151 * t212) * t146) * t208 + (t143 * t213 + t181 * t210) * t132 - ((t131 - t212) * t150 + ((-t146 * t179 + 0.1e1) * t175 + (t146 - t179) * t135) * t151) * t172 * t233) * t171) * t136;
	t1 = [t181 * t168 * t199 + (t175 * t193 - t213 * t223) * t161, t131, t131, 0, 0, 0; (t142 * t200 + (-t142 * t221 + (qJD(1) * t134 + t130) * t234) * t136) * t179 + (t143 * t200 * t134 + (-((t199 - t221 + (t135 * t168 * t224 + t221) * t161) * t150 + (t198 * t227 - t135 * t171 + (-t166 * t205 + (t135 - 0.2e1 * t220) * t171) * t161) * t151) * t208 + (-t143 * t221 + t171 * t210) * t134 + (-t142 + ((-t176 + t177) * t195 + t239 * t207) * t143) * t171 * qJD(1)) * t136) * t181, t127, t127, 0, 0, 0; 0.2e1 * (t153 * t190 + t158 * t229) * t237 + (0.2e1 * t158 * t206 - t197 * t153 * t216 + (t175 * t218 + t192) * t230 + (t158 * t140 + t190 * t141 - t192 * t228 - (t171 * t175 * t180 + t197 * t178) * t159 * t179) * t154) * t147, t128, t128, 0, t209 + 0.2e1 * (-t140 * t154 * t147 + (-t147 * t235 - t154 * t237) * t159) * t159, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:36:55
	% EndTime: 2019-10-10 11:36:57
	% DurationCPUTime: 1.65s
	% Computational Cost: add. (8326->125), mult. (8382->273), div. (1515->15), fcn. (10508->9), ass. (0->119)
	t202 = qJ(2) + qJ(3) + pkin(10);
	t199 = cos(t202);
	t208 = sin(qJ(5));
	t283 = sin(qJ(1));
	t239 = t283 * t208;
	t209 = cos(qJ(5));
	t210 = cos(qJ(1));
	t260 = t210 * t209;
	t185 = t199 * t260 + t239;
	t179 = 0.1e1 / t185 ^ 2;
	t198 = sin(t202);
	t194 = t198 ^ 2;
	t207 = t210 ^ 2;
	t271 = t194 * t207;
	t247 = t179 * t271;
	t174 = 0.1e1 + t247;
	t234 = qJD(1) * t283;
	t203 = qJD(2) + qJD(3);
	t264 = t203 * t210;
	t242 = t198 * t264;
	t220 = t199 * t234 + t242;
	t233 = t283 * qJD(5);
	t261 = t210 * t208;
	t164 = (-qJD(5) * t199 + qJD(1)) * t261 + (t233 - t220) * t209;
	t178 = 0.1e1 / t185;
	t278 = t164 * t178 * t179;
	t228 = t271 * t278;
	t243 = t198 * t203 * t207;
	t286 = (-t228 + (-t194 * t210 * t234 + t199 * t243) * t179) / t174 ^ 2;
	t267 = t198 * t210;
	t181 = t199 * t239 + t260;
	t225 = t208 * t233;
	t256 = qJD(5) * t210;
	t236 = t209 * t256;
	t163 = t181 * qJD(1) - t199 * t236 + t208 * t242 - t225;
	t238 = t283 * t209;
	t184 = t199 * t261 - t238;
	t195 = 0.1e1 / t198;
	t196 = 0.1e1 / t198 ^ 2;
	t205 = 0.1e1 / t208 ^ 2;
	t257 = qJD(5) * t209;
	t237 = t205 * t257;
	t204 = 0.1e1 / t208;
	t265 = t203 * t204;
	t241 = t199 * t265;
	t270 = t195 * t204;
	t285 = (t195 * t237 + t196 * t241) * t184 + t163 * t270;
	t268 = t198 * t208;
	t173 = atan2(-t181, t268);
	t168 = cos(t173);
	t167 = sin(t173);
	t277 = t167 * t181;
	t162 = t168 * t268 - t277;
	t159 = 0.1e1 / t162;
	t160 = 0.1e1 / t162 ^ 2;
	t284 = 0.2e1 * t184;
	t176 = t181 ^ 2;
	t269 = t196 * t205;
	t175 = t176 * t269 + 0.1e1;
	t171 = 0.1e1 / t175;
	t266 = t199 * t203;
	t221 = t198 * t257 + t208 * t266;
	t245 = t181 * t269;
	t226 = t209 * t234;
	t240 = t283 * t198;
	t227 = t203 * t240;
	t259 = qJD(1) * t210;
	t165 = t209 * t233 * t199 - t226 + (t259 * t199 - t227 - t256) * t208;
	t248 = t165 * t270;
	t151 = (t221 * t245 - t248) * t171;
	t218 = -t151 * t181 + t221;
	t147 = (-t151 * t268 - t165) * t167 + t218 * t168;
	t161 = t159 * t160;
	t282 = t147 * t161;
	t197 = t195 / t194;
	t206 = t204 * t205;
	t281 = (t165 * t245 + (-t196 * t206 * t257 - t197 * t205 * t266) * t176) / t175 ^ 2;
	t280 = t160 * t184;
	t279 = t163 * t160;
	t276 = t167 * t184;
	t275 = t167 * t198;
	t274 = t168 * t181;
	t273 = t168 * t184;
	t272 = t168 * t199;
	t263 = t205 * t209;
	t262 = t210 * t159;
	t258 = qJD(5) * t208;
	t177 = t184 ^ 2;
	t157 = t160 * t177 + 0.1e1;
	t255 = 0.2e1 * (-t177 * t282 - t184 * t279) / t157 ^ 2;
	t254 = 0.2e1 * t286;
	t253 = -0.2e1 * t281;
	t252 = t161 * t284;
	t251 = t195 * t281;
	t250 = t160 * t276;
	t246 = t181 * t270;
	t244 = t196 * t199 * t204;
	t223 = t181 * t244 + t283;
	t158 = t223 * t171;
	t235 = t283 - t158;
	t232 = t159 * t255;
	t231 = t160 * t255;
	t230 = t267 * t284;
	t229 = t204 * t251;
	t183 = t199 * t238 - t261;
	t224 = t181 * t263 - t183 * t204;
	t222 = t179 * t183 * t210 - t283 * t178;
	t169 = 0.1e1 / t174;
	t166 = t185 * qJD(1) - t199 * t225 - t209 * t227 - t236;
	t155 = 0.1e1 / t157;
	t154 = t224 * t195 * t171;
	t150 = (-t167 + (t168 * t246 + t167) * t171) * t184;
	t149 = -t158 * t274 + (t235 * t275 + t272) * t208;
	t148 = t168 * t198 * t209 - t167 * t183 + (-t167 * t268 - t274) * t154;
	t146 = t223 * t253 + (t165 * t244 + t259 + (-t195 * t265 + (-t196 * t237 - 0.2e1 * t197 * t241) * t199) * t181) * t171;
	t144 = -0.2e1 * t224 * t251 + (-t224 * t196 * t266 + (t165 * t263 - t166 * t204 + (t183 * t263 + (-0.2e1 * t206 * t209 ^ 2 - t204) * t181) * qJD(5)) * t195) * t171;
	t143 = (t178 * t199 * t210 + t209 * t247) * t254 + (0.2e1 * t209 * t228 + t220 * t178 + ((t164 * t210 - 0.2e1 * t209 * t243) * t199 + (t207 * t258 + 0.2e1 * t210 * t226) * t194) * t179) * t169;
	t142 = t149 * t184 * t231 + (-(-t146 * t274 + (t151 * t277 - t165 * t168) * t158) * t280 + (t147 * t252 + t279) * t149 + (-t198 * t262 - (-t158 * t275 + t167 * t240 + t272) * t280) * t257) * t155 + (t232 * t267 + ((-t203 * t262 - (t235 * t203 - t151) * t250) * t199 + (t159 * t234 + (t210 * t147 - (-t146 + t259) * t276 - (t235 * t151 - t203) * t273) * t160) * t198) * t155) * t208;
	t1 = [t285 * t171 + t229 * t284, t146, t146, 0, t144, 0; t181 * t232 + (-t165 * t159 + (t147 * t181 + t150 * t163) * t160) * t155 + (t150 * t231 + (0.2e1 * t150 * t282 + (t163 * t171 - t163 - (-t151 * t171 * t246 + t253) * t184) * t160 * t167 + (-(-0.2e1 * t181 * t229 - t151) * t280 + (-(t151 + t248) * t184 + t285 * t181) * t160 * t171) * t168) * t155) * t184, t142, t142, 0, (t148 * t280 - t159 * t185) * t255 + (t148 * t279 + t164 * t159 + (t148 * t252 - t160 * t185) * t147 - (-t198 * t258 + t209 * t266 - t144 * t181 - t154 * t165 + (-t154 * t268 - t183) * t151) * t160 * t273 - (-t166 + (-t144 * t208 - t151 * t209) * t198 - t218 * t154) * t250) * t155, 0; t222 * t198 * t254 + (-t222 * t266 + ((qJD(1) * t178 + 0.2e1 * t183 * t278) * t210 + (-t283 * t164 - t166 * t210 + t183 * t234) * t179) * t198) * t169, t143, t143, 0, t179 * t230 * t286 + (t230 * t278 + (t163 * t267 + (t198 * t234 - t199 * t264) * t184) * t179) * t169, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end