% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR10V2
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
%   Wie in S6RRRRRR10V2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:38
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR10V2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR10V2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:03
	% DurationCPUTime: 1.13s
	% Computational Cost: add. (3645->97), mult. (3810->208), div. (753->12), fcn. (4455->9), ass. (0->95)
	t172 = sin(qJ(1));
	t234 = 0.2e1 * t172;
	t168 = t172 ^ 2;
	t170 = qJ(2) + qJ(3);
	t165 = sin(t170);
	t161 = t165 ^ 2;
	t166 = cos(t170);
	t163 = 0.1e1 / t166 ^ 2;
	t219 = t161 * t163;
	t156 = t168 * t219 + 0.1e1;
	t154 = 0.1e1 / t156;
	t162 = 0.1e1 / t166;
	t174 = cos(qJ(1));
	t205 = qJD(1) * t174;
	t195 = t165 * t205;
	t167 = qJD(2) + qJD(3);
	t213 = t167 * t172;
	t198 = t163 * t213;
	t128 = (-(-t166 * t213 - t195) * t162 + t161 * t198) * t154;
	t233 = t128 - t213;
	t173 = cos(qJ(4));
	t207 = t173 * t174;
	t171 = sin(qJ(4));
	t209 = t172 * t171;
	t150 = t166 * t207 + t209;
	t210 = t172 * t165;
	t153 = atan2(-t210, -t166);
	t152 = cos(t153);
	t151 = sin(t153);
	t199 = t151 * t210;
	t138 = -t152 * t166 - t199;
	t135 = 0.1e1 / t138;
	t144 = 0.1e1 / t150;
	t136 = 0.1e1 / t138 ^ 2;
	t145 = 0.1e1 / t150 ^ 2;
	t232 = t154 - 0.1e1;
	t221 = t152 * t165;
	t123 = (-t128 * t172 + t167) * t221 + (t233 * t166 - t195) * t151;
	t231 = t123 * t135 * t136;
	t183 = t166 * t209 + t207;
	t212 = t167 * t174;
	t196 = t165 * t212;
	t132 = t183 * qJD(1) - t150 * qJD(4) + t171 * t196;
	t208 = t172 * t173;
	t211 = t171 * t174;
	t149 = t166 * t211 - t208;
	t143 = t149 ^ 2;
	t142 = t143 * t145 + 0.1e1;
	t224 = t145 * t149;
	t189 = -qJD(1) * t166 + qJD(4);
	t190 = qJD(4) * t166 - qJD(1);
	t133 = -t190 * t211 + (t189 * t172 - t196) * t173;
	t229 = t133 * t144 * t145;
	t230 = (-t132 * t224 - t143 * t229) / t142 ^ 2;
	t160 = t165 * t161;
	t216 = t162 * t165;
	t182 = t167 * (t160 * t162 * t163 + t216);
	t217 = t161 * t172;
	t187 = t205 * t217;
	t228 = (t163 * t187 + t168 * t182) / t156 ^ 2;
	t227 = t136 * t165;
	t226 = t136 * t174;
	t225 = t144 * t171;
	t223 = t149 * t173;
	t222 = t151 * t172;
	t220 = t161 * t162;
	t169 = t174 ^ 2;
	t218 = t161 * t169;
	t215 = t165 * t174;
	t214 = t166 * t167;
	t206 = qJD(1) * t172;
	t131 = t136 * t218 + 0.1e1;
	t204 = 0.2e1 * (-t218 * t231 + (t165 * t169 * t214 - t187) * t136) / t131 ^ 2;
	t203 = 0.2e1 * t231;
	t202 = -0.2e1 * t230;
	t201 = t136 * t215;
	t200 = t149 * t229;
	t194 = 0.1e1 + t219;
	t193 = t165 * t204;
	t192 = -0.2e1 * t165 * t228;
	t191 = t228 * t234;
	t188 = t152 * t154 * t220;
	t186 = t194 * t174;
	t185 = t189 * t174;
	t184 = t145 * t223 - t225;
	t148 = -t166 * t208 + t211;
	t140 = 0.1e1 / t142;
	t139 = t194 * t172 * t154;
	t129 = 0.1e1 / t131;
	t127 = (t232 * t165 * t151 - t172 * t188) * t174;
	t125 = -t166 * t222 + t221 + (t151 * t166 - t152 * t210) * t139;
	t124 = -t194 * t191 + (qJD(1) * t186 + t182 * t234) * t154;
	t121 = t184 * t202 * t215 + (t184 * t166 * t212 + (-t184 * t206 + ((-qJD(4) * t144 - 0.2e1 * t200) * t173 + (-t132 * t173 + (-qJD(4) * t149 + t133) * t171) * t145) * t174) * t165) * t140;
	t120 = (t125 * t227 - t135 * t166) * t174 * t204 + ((-t135 * t206 + (-t125 * t167 - t123) * t226) * t166 + (-t135 * t212 - (-t124 * t152 * t172 - t233 * t151 + (t128 * t222 - t151 * t167 - t152 * t205) * t139) * t201 + (t136 * t206 + t174 * t203) * t125 - ((t124 - t205) * t151 + ((-t139 * t172 + 0.1e1) * t167 + (t139 - t172) * t128) * t152) * t166 * t226) * t165) * t129;
	t1 = [t162 * t174 * t192 + (t167 * t186 - t206 * t216) * t154, t124, t124, 0, 0, 0; (t135 * t193 + (-t135 * t214 + (qJD(1) * t127 + t123) * t227) * t129) * t172 + (t136 * t193 * t127 + (-((t192 - t214 + (t128 * t162 * t217 + t214) * t154) * t151 + (t191 * t220 - t128 * t165 + (-t160 * t198 + (t128 - 0.2e1 * t213) * t165) * t154) * t152) * t201 + (-t136 * t214 + t165 * t203) * t127 + (-t135 + ((-t168 + t169) * t188 + t232 * t199) * t136) * t165 * qJD(1)) * t129) * t174, t120, t120, 0, 0, 0; 0.2e1 * (t144 * t183 + t148 * t224) * t230 + (0.2e1 * t148 * t200 - t190 * t144 * t208 + (t167 * t210 + t185) * t225 + (t148 * t132 + t183 * t133 - t185 * t223 - (t165 * t167 * t173 + t190 * t171) * t149 * t172) * t145) * t140, t121, t121, t202 + 0.2e1 * (-t132 * t140 * t145 + (-t140 * t229 - t145 * t230) * t149) * t149, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:02
	% EndTime: 2019-10-10 13:38:04
	% DurationCPUTime: 2.29s
	% Computational Cost: add. (6566->153), mult. (9976->323), div. (1522->14), fcn. (12462->11), ass. (0->137)
	t240 = qJ(2) + qJ(3);
	t235 = cos(t240);
	t242 = sin(qJ(4));
	t331 = sin(qJ(1));
	t276 = t331 * t242;
	t244 = cos(qJ(4));
	t245 = cos(qJ(1));
	t302 = t245 * t244;
	t219 = t235 * t302 + t276;
	t241 = sin(qJ(5));
	t243 = cos(qJ(5));
	t234 = sin(t240);
	t311 = t234 * t245;
	t203 = t219 * t241 - t243 * t311;
	t335 = 0.2e1 * t203;
	t236 = qJD(2) + qJD(3);
	t334 = -qJD(5) * t244 + t236;
	t204 = t219 * t243 + t241 * t311;
	t198 = 0.1e1 / t204;
	t199 = 0.1e1 / t204 ^ 2;
	t323 = t199 * t203;
	t261 = t241 * t198 - t243 * t323;
	t215 = t235 * t276 + t302;
	t270 = t331 * qJD(4);
	t264 = t242 * t270;
	t298 = qJD(4) * t245;
	t273 = t244 * t298;
	t307 = t236 * t245;
	t281 = t234 * t307;
	t193 = qJD(1) * t215 - t235 * t273 + t242 * t281 - t264;
	t275 = t331 * t244;
	t303 = t245 * t242;
	t218 = t235 * t303 - t275;
	t231 = 0.1e1 / t234;
	t232 = 0.1e1 / t234 ^ 2;
	t238 = 0.1e1 / t242 ^ 2;
	t299 = qJD(4) * t244;
	t274 = t238 * t299;
	t237 = 0.1e1 / t242;
	t309 = t236 * t237;
	t280 = t235 * t309;
	t315 = t231 * t237;
	t333 = (t231 * t274 + t232 * t280) * t218 + t193 * t315;
	t313 = t234 * t242;
	t209 = atan2(-t215, t313);
	t206 = cos(t209);
	t205 = sin(t209);
	t321 = t205 * t215;
	t189 = t206 * t313 - t321;
	t186 = 0.1e1 / t189;
	t187 = 0.1e1 / t189 ^ 2;
	t332 = 0.2e1 * t218;
	t213 = t215 ^ 2;
	t314 = t232 * t238;
	t210 = t213 * t314 + 0.1e1;
	t207 = 0.1e1 / t210;
	t310 = t235 * t236;
	t258 = t234 * t299 + t242 * t310;
	t283 = t215 * t314;
	t271 = qJD(1) * t331;
	t301 = qJD(1) * t245;
	t195 = -t236 * t234 * t276 - t242 * t298 - t244 * t271 + (t242 * t301 + t244 * t270) * t235;
	t286 = t195 * t315;
	t177 = (t258 * t283 - t286) * t207;
	t256 = -t177 * t215 + t258;
	t172 = (-t177 * t313 - t195) * t205 + t256 * t206;
	t188 = t186 * t187;
	t330 = t172 * t188;
	t255 = qJD(5) * t219 + t234 * t271 - t235 * t307;
	t194 = (-qJD(4) * t235 + qJD(1)) * t303 + (-t235 * t271 + t270 - t281) * t244;
	t263 = qJD(5) * t311 + t194;
	t179 = t241 * t263 + t243 * t255;
	t197 = t203 ^ 2;
	t192 = t197 * t199 + 0.1e1;
	t180 = -t241 * t255 + t243 * t263;
	t200 = t198 * t199;
	t327 = t180 * t200;
	t329 = (t179 * t323 - t197 * t327) / t192 ^ 2;
	t233 = t231 * t232;
	t239 = t237 * t238;
	t328 = (t195 * t283 + (-t232 * t239 * t299 - t233 * t238 * t310) * t213) / t210 ^ 2;
	t326 = t187 * t218;
	t190 = 0.1e1 / t192;
	t325 = t190 * t199;
	t324 = t193 * t187;
	t312 = t234 * t244;
	t212 = (t235 * t241 - t243 * t312) * t245;
	t322 = t203 * t212;
	t320 = t205 * t218;
	t319 = t205 * t234;
	t318 = t206 * t215;
	t317 = t206 * t218;
	t316 = t206 * t235;
	t308 = t236 * t244;
	t306 = t238 * t244;
	t304 = t245 * t186;
	t300 = qJD(4) * t242;
	t214 = t218 ^ 2;
	t184 = t187 * t214 + 0.1e1;
	t296 = 0.2e1 * (-t214 * t330 - t218 * t324) / t184 ^ 2;
	t295 = -0.2e1 * t329;
	t294 = 0.2e1 * t329;
	t293 = -0.2e1 * t328;
	t292 = t188 * t332;
	t291 = t199 * t329;
	t290 = t231 * t328;
	t289 = t179 * t325;
	t288 = t187 * t320;
	t285 = t203 * t327;
	t284 = t215 * t315;
	t282 = t232 * t235 * t237;
	t278 = t331 * t234;
	t277 = t331 * t235;
	t259 = t215 * t282 + t331;
	t185 = t259 * t207;
	t272 = t331 - t185;
	t269 = t186 * t296;
	t268 = t187 * t296;
	t266 = t237 * t290;
	t265 = t234 * t275;
	t196 = qJD(1) * t219 - t235 * t264 - t236 * t265 - t273;
	t262 = -qJD(5) * t278 - t196;
	t217 = t235 * t275 - t303;
	t260 = t215 * t306 - t217 * t237;
	t254 = -qJD(5) * t217 + t234 * t301 + t236 * t277;
	t211 = (-t235 * t243 - t241 * t312) * t245;
	t202 = -t217 * t243 - t241 * t278;
	t182 = 0.1e1 / t184;
	t181 = t260 * t231 * t207;
	t176 = (-t205 + (t206 * t284 + t205) * t207) * t218;
	t175 = -t185 * t318 + (t272 * t319 + t316) * t242;
	t173 = t206 * t312 - t205 * t217 + (-t205 * t313 - t318) * t181;
	t171 = t259 * t293 + (t195 * t282 + t301 + (-t231 * t309 + (-t232 * t274 - 0.2e1 * t233 * t280) * t235) * t215) * t207;
	t169 = -0.2e1 * t260 * t290 + (-t260 * t232 * t310 + (t195 * t306 - t196 * t237 + (t217 * t306 + (-0.2e1 * t239 * t244 ^ 2 - t237) * t215) * qJD(4)) * t231) * t207;
	t168 = (-t198 * t211 + t199 * t322) * t294 + (-t212 * t179 * t199 + (-t199 * t211 + 0.2e1 * t200 * t322) * t180 + ((t241 * t265 + t243 * t277) * t198 - (-t241 * t277 + t243 * t265) * t323) * qJD(1) + (((t241 * t300 + t243 * t334) * t198 - (-t241 * t334 + t243 * t300) * t323) * t234 + t261 * t235 * (qJD(5) - t308)) * t245) * t190;
	t167 = t175 * t218 * t268 + (-(-t171 * t318 + (t177 * t321 - t195 * t206) * t185) * t326 + (t172 * t292 + t324) * t175 + (-t234 * t304 - (-t185 * t319 + t205 * t278 + t316) * t326) * t299) * t182 + (t269 * t311 + ((-t236 * t304 - (t236 * t272 - t177) * t288) * t235 + (t186 * t271 + (t245 * t172 - (-t171 + t301) * t320 - (t177 * t272 - t236) * t317) * t187) * t234) * t182) * t242;
	t1 = [t207 * t333 + t266 * t332, t171, t171, t169, 0, 0; t215 * t269 + (-t195 * t186 + (t172 * t215 + t176 * t193) * t187) * t182 + (t176 * t268 + (0.2e1 * t176 * t330 + (t193 * t207 - t193 - (-t177 * t207 * t284 + t293) * t218) * t187 * t205 + (-(-0.2e1 * t215 * t266 - t177) * t326 + (-(t177 + t286) * t218 + t333 * t215) * t187 * t207) * t206) * t182) * t218, t167, t167, (t173 * t326 - t186 * t219) * t296 + (t173 * t324 + t194 * t186 + (t173 * t292 - t187 * t219) * t172 - (-t234 * t300 + t235 * t308 - t169 * t215 - t181 * t195 + (-t181 * t313 - t217) * t177) * t187 * t317 - (-t196 + (-t169 * t242 - t177 * t244) * t234 - t256 * t181) * t288) * t182, 0, 0; (t291 * t335 - t289) * t202 + (-t180 * t325 + t198 * t295) * (-t217 * t241 + t243 * t278) + ((t241 * t262 + t243 * t254) * t198 - (-t241 * t254 + t243 * t262) * t323 + 0.2e1 * t202 * t285) * t190, t168, t168, t261 * t218 * t294 + (t261 * t193 + ((-qJD(5) * t198 - 0.2e1 * t285) * t243 + (t179 * t243 + (-qJD(5) * t203 + t180) * t241) * t199) * t218) * t190, t295 + (t289 + (-t190 * t327 - t291) * t203) * t335, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:38:03
	% EndTime: 2019-10-10 13:38:06
	% DurationCPUTime: 3.69s
	% Computational Cost: add. (17109->206), mult. (24753->390), div. (1249->12), fcn. (30668->13), ass. (0->164)
	t351 = qJ(2) + qJ(3);
	t348 = sin(t351);
	t353 = sin(qJ(5));
	t357 = cos(qJ(5));
	t350 = qJD(2) + qJD(3);
	t358 = cos(qJ(4));
	t389 = qJD(5) * t358 - t350;
	t354 = sin(qJ(4));
	t419 = qJD(4) * t354;
	t460 = t348 * (t389 * t353 + t357 * t419);
	t349 = cos(t351);
	t355 = sin(qJ(1));
	t430 = t350 * t355;
	t359 = cos(qJ(1));
	t432 = t348 * t359;
	t459 = qJD(1) * t432 + t349 * t430;
	t421 = t359 * t354;
	t425 = t355 * t358;
	t333 = t349 * t425 - t421;
	t433 = t348 * t357;
	t315 = t333 * t353 - t355 * t433;
	t427 = t353 * t358;
	t431 = t349 * t357;
	t329 = t348 * t427 + t431;
	t308 = atan2(-t315, t329);
	t299 = sin(t308);
	t300 = cos(t308);
	t277 = -t299 * t315 + t300 * t329;
	t275 = 0.1e1 / t277 ^ 2;
	t422 = t358 * t359;
	t426 = t355 * t354;
	t335 = t349 * t422 + t426;
	t320 = t335 * t353 - t357 * t432;
	t314 = t320 ^ 2;
	t273 = t275 * t314 + 0.1e1;
	t387 = -qJD(1) * t349 + qJD(4);
	t388 = -qJD(4) * t349 + qJD(1);
	t429 = t350 * t359;
	t304 = t388 * t421 + (-t348 * t429 + t387 * t355) * t358;
	t321 = t335 * t357 + t353 * t432;
	t420 = qJD(1) * t355;
	t455 = t348 * t420 - t349 * t429;
	t278 = t321 * qJD(5) + t304 * t353 + t357 * t455;
	t444 = t278 * t275;
	t313 = t315 ^ 2;
	t327 = 0.1e1 / t329 ^ 2;
	t307 = t313 * t327 + 0.1e1;
	t301 = 0.1e1 / t307;
	t418 = qJD(4) * t358;
	t394 = t359 * t418;
	t395 = t355 * t419;
	t404 = t348 * t430;
	t306 = t335 * qJD(1) - t349 * t395 - t358 * t404 - t394;
	t317 = t348 * t353 * t355 + t333 * t357;
	t280 = t317 * qJD(5) + t306 * t353 - t459 * t357;
	t390 = t350 * t358 - qJD(5);
	t382 = t349 * t390;
	t367 = -t353 * t419 + t389 * t357;
	t456 = t367 * t348;
	t289 = t353 * t382 + t456;
	t326 = 0.1e1 / t329;
	t437 = t315 * t327;
	t379 = -t280 * t326 + t289 * t437;
	t264 = t379 * t301;
	t383 = -t299 * t329 - t300 * t315;
	t258 = t383 * t264 - t280 * t299 + t289 * t300;
	t274 = 0.1e1 / t277;
	t276 = t274 * t275;
	t449 = t258 * t276;
	t414 = 0.2e1 * (-t314 * t449 + t320 * t444) / t273 ^ 2;
	t458 = t289 * t327;
	t375 = t349 * t426 + t422;
	t434 = t348 * t354;
	t405 = t315 * t434;
	t374 = -t326 * t375 + t327 * t405;
	t457 = t353 * t374;
	t330 = -t349 * t353 + t358 * t433;
	t324 = t330 * t359;
	t400 = t350 * t421;
	t454 = (t354 * t420 - t394) * t348 + qJD(6) * t324 - t349 * t400;
	t417 = qJD(5) * t348;
	t281 = (-qJD(5) * t333 + t459) * t353 + (t355 * t417 + t306) * t357;
	t334 = t349 * t421 - t425;
	t352 = sin(qJ(6));
	t356 = cos(qJ(6));
	t298 = t321 * t356 + t334 * t352;
	t292 = 0.1e1 / t298;
	t293 = 0.1e1 / t298 ^ 2;
	t452 = -0.2e1 * t315;
	t451 = 0.2e1 * t320;
	t279 = (t359 * t417 + t304) * t357 + (-qJD(5) * t335 - t455) * t353;
	t303 = t375 * qJD(1) + t348 * t400 - t349 * t394 - t395;
	t266 = t298 * qJD(6) + t279 * t352 + t303 * t356;
	t297 = t321 * t352 - t334 * t356;
	t291 = t297 ^ 2;
	t285 = t291 * t293 + 0.1e1;
	t442 = t293 * t297;
	t415 = qJD(6) * t297;
	t267 = t279 * t356 - t303 * t352 - t415;
	t446 = t267 * t292 * t293;
	t448 = (t266 * t442 - t291 * t446) / t285 ^ 2;
	t443 = t326 * t458;
	t447 = (t280 * t437 - t313 * t443) / t307 ^ 2;
	t445 = t275 * t320;
	t441 = t297 * t352;
	t440 = t299 * t320;
	t439 = t300 * t320;
	t438 = t315 * t326;
	t436 = t334 * t353;
	t435 = t334 * t357;
	t428 = t352 * t292;
	t424 = t356 * t292;
	t423 = t356 * t297;
	t416 = qJD(5) * t357;
	t413 = -0.2e1 * t448;
	t412 = 0.2e1 * t448;
	t411 = -0.2e1 * t447;
	t410 = t276 * t451;
	t409 = t326 * t447;
	t408 = t275 * t440;
	t407 = t275 * t439;
	t406 = t297 * t446;
	t403 = t348 * t421;
	t393 = t258 * t410;
	t392 = t443 * t452;
	t391 = 0.2e1 * t406;
	t384 = qJD(6) * t435 + t304;
	t296 = -t317 * t356 - t352 * t375;
	t295 = -t317 * t352 + t356 * t375;
	t381 = t390 * t353;
	t380 = -qJD(6) * t403 + t330 * t420 + (-t390 * t431 + t460) * t359;
	t378 = t293 * t423 - t428;
	t377 = -t317 * t326 + t330 * t437;
	t322 = t329 * t355;
	t331 = t349 * t427 - t433;
	t376 = t322 * t326 + t331 * t437;
	t373 = -t349 * t350 * t354 - t348 * t418;
	t370 = -t299 + (t300 * t438 + t299) * t301;
	t369 = qJD(1) * t329;
	t368 = qJD(5) * t436 + qJD(6) * t335 + t303 * t357;
	t323 = t329 * t359;
	t312 = -t324 * t356 - t352 * t403;
	t311 = -t324 * t352 + t356 * t403;
	t310 = t335 * t352 - t356 * t435;
	t309 = -t335 * t356 - t352 * t435;
	t305 = t388 * t425 + (t387 * t359 + t404) * t354;
	t290 = t357 * t382 - t460;
	t288 = -t348 * t381 + t367 * t349;
	t287 = t289 * t355 + t359 * t369;
	t283 = 0.1e1 / t285;
	t271 = 0.1e1 / t273;
	t270 = t301 * t457;
	t269 = t376 * t301;
	t268 = t377 * t301;
	t263 = t370 * t320;
	t261 = (t299 * t375 - t300 * t434) * t353 - t383 * t270;
	t260 = t383 * t269 + t299 * t322 + t300 * t331;
	t259 = t383 * t268 - t299 * t317 + t300 * t330;
	t256 = t376 * t411 + (t331 * t392 + t287 * t326 + (t280 * t331 + t288 * t315 - t289 * t322) * t327) * t301;
	t255 = t377 * t411 + (t330 * t392 - t281 * t326 + (t280 * t330 + t289 * t317 + t290 * t315) * t327) * t301;
	t254 = 0.2e1 * t447 * t457 + (-t374 * t416 + (0.2e1 * t405 * t443 - t305 * t326 + (-t280 * t434 - t289 * t375 + t373 * t315) * t327) * t353) * t301;
	t253 = (-t292 * t311 + t312 * t442) * t412 + (t312 * t391 + t380 * t428 - t454 * t424 + (-t312 * t266 - t311 * t267 - t380 * t423 - t441 * t454) * t293) * t283;
	t252 = (t260 * t445 + t274 * t323) * t414 + (t260 * t393 + (t323 * t258 - t260 * t278 - (-t256 * t315 - t269 * t280 + t288 + (-t269 * t329 + t322) * t264) * t439 - (-t256 * t329 - t269 * t289 + t287 + (t269 * t315 - t331) * t264) * t440) * t275 + (t355 * t369 + (-t349 * t381 - t456) * t359) * t274) * t271;
	t1 = [t409 * t451 + (-t278 * t326 + t320 * t458) * t301, t256, t256, t254, t255, 0; t315 * t274 * t414 + (-t280 * t274 + (t258 * t315 - t263 * t278) * t275) * t271 + (t263 * t275 * t414 + (0.2e1 * t263 * t449 - (-t264 * t301 * t438 + t411) * t408 - (t409 * t452 - t264 + (t264 - t379) * t301) * t407 - t370 * t444) * t271) * t320, t252, t252, (t261 * t445 + t274 * t436) * t414 + (-t261 * t444 + (t303 * t353 - t334 * t416) * t274 + (t261 * t410 + t275 * t436) * t258 - (t375 * t416 - t254 * t329 + t270 * t289 - t305 * t353 + (-t270 * t315 + t353 * t434) * t264) * t408 - (-t416 * t434 - t254 * t315 - (-t264 * t329 - t280) * t270 + (t264 * t375 + t373) * t353) * t407) * t271, (t259 * t445 - t274 * t321) * t414 + (t259 * t393 + t279 * t274 + (-t321 * t258 - t259 * t278 - (-t255 * t315 - t268 * t280 + t290 + (-t268 * t329 - t317) * t264) * t439 - (-t255 * t329 - t268 * t289 - t281 + (t268 * t315 - t330) * t264) * t440) * t275) * t271, 0; (-t292 * t295 + t296 * t442) * t412 + ((t296 * qJD(6) - t281 * t352 - t305 * t356) * t292 + t296 * t391 + (-t295 * t267 - (-t295 * qJD(6) - t281 * t356 + t305 * t352) * t297 - t296 * t266) * t293) * t283, t253, t253, (-t292 * t309 + t310 * t442) * t412 + (t310 * t391 - t384 * t424 + t368 * t428 + (-t310 * t266 - t309 * t267 - t368 * t423 - t384 * t441) * t293) * t283, t378 * t320 * t413 + (t378 * t278 + ((-qJD(6) * t292 - 0.2e1 * t406) * t356 + (t266 * t356 + (t267 - t415) * t352) * t293) * t320) * t283, t413 + 0.2e1 * (t266 * t293 * t283 + (-t283 * t446 - t293 * t448) * t297) * t297;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end