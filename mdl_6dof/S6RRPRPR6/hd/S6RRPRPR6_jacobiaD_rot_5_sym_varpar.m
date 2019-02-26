% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR6_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:37
% EndTime: 2019-02-26 21:40:39
% DurationCPUTime: 1.48s
% Computational Cost: add. (6793->127), mult. (19653->266), div. (691->12), fcn. (25511->13), ass. (0->115)
t241 = sin(pkin(11));
t242 = cos(pkin(11));
t245 = sin(qJ(2));
t248 = cos(qJ(2));
t231 = t241 * t245 - t248 * t242;
t243 = cos(pkin(6));
t228 = t231 * t243;
t223 = qJD(2) * t228;
t264 = t241 * t248 + t242 * t245;
t230 = t264 * qJD(2);
t246 = sin(qJ(1));
t249 = cos(qJ(1));
t306 = sin(pkin(6));
t273 = t249 * t306;
t229 = t264 * t243;
t289 = t246 * t229;
t311 = -t246 * t230 - qJD(1) * t289 + (-qJD(1) * t231 - t223) * t249 - qJD(4) * t273;
t209 = t229 * t249 - t246 * t231;
t244 = sin(qJ(4));
t247 = cos(qJ(4));
t198 = t209 * t244 + t247 * t273;
t274 = t248 * t306;
t276 = t245 * t306;
t256 = -t241 * t274 - t242 * t276;
t217 = -t243 * t247 - t244 * t256;
t192 = atan2(-t198, t217);
t187 = sin(t192);
t188 = cos(t192);
t175 = -t187 * t198 + t188 * t217;
t173 = 0.1e1 / t175 ^ 2;
t265 = -t231 * t249 - t289;
t275 = t246 * t306;
t259 = -t244 * t265 + t247 * t275;
t196 = t259 ^ 2;
t171 = t173 * t196 + 0.1e1;
t184 = -t209 * qJD(1) + t246 * t223 - t249 * t230;
t204 = t244 * t275 + t247 * t265;
t272 = qJD(1) * t306;
t269 = t249 * t272;
t176 = t204 * qJD(4) + t184 * t244 - t247 * t269;
t300 = t173 * t259;
t195 = t198 ^ 2;
t215 = 0.1e1 / t217 ^ 2;
t191 = t195 * t215 + 0.1e1;
t189 = 0.1e1 / t191;
t267 = t246 * t272;
t286 = qJD(4) * t247;
t178 = t209 * t286 + t244 * t311 - t247 * t267;
t218 = t243 * t244 - t247 * t256;
t226 = -t241 * t276 + t242 * t274;
t222 = t226 * qJD(2);
t193 = t218 * qJD(4) + t222 * t244;
t214 = 0.1e1 / t217;
t294 = t198 * t215;
t263 = -t178 * t214 + t193 * t294;
t164 = t263 * t189;
t266 = -t187 * t217 - t188 * t198;
t160 = t266 * t164 - t187 * t178 + t188 * t193;
t172 = 0.1e1 / t175;
t174 = t172 * t173;
t304 = t160 * t174;
t285 = 0.2e1 * (-t176 * t300 - t196 * t304) / t171 ^ 2;
t310 = t193 * t215;
t208 = -t228 * t249 - t246 * t264;
t261 = -t208 * t214 + t226 * t294;
t309 = t244 * t261;
t179 = (-qJD(4) * t209 + t267) * t244 + t311 * t247;
t211 = t246 * t228 - t249 * t264;
t205 = 0.1e1 / t211;
t206 = 0.1e1 / t211 ^ 2;
t308 = -0.2e1 * t198;
t307 = -0.2e1 * t259;
t177 = t259 * qJD(4) + t184 * t247 + t244 * t269;
t197 = t204 ^ 2;
t182 = t197 * t206 + 0.1e1;
t224 = t243 * t230;
t260 = t231 * qJD(2);
t183 = t208 * qJD(1) - t246 * t224 - t249 * t260;
t207 = t205 * t206;
t293 = t204 * t206;
t303 = (t183 * t197 * t207 + t177 * t293) / t182 ^ 2;
t296 = t214 * t310;
t302 = (t178 * t294 - t195 * t296) / t191 ^ 2;
t301 = t173 * t176;
t299 = t183 * t206;
t298 = t187 * t259;
t297 = t188 * t259;
t295 = t198 * t214;
t292 = t206 * t208;
t291 = t211 * t244;
t284 = 0.2e1 * t303;
t283 = -0.2e1 * t302;
t282 = t174 * t307;
t281 = -0.2e1 * t204 * t207;
t280 = t205 * t303;
t279 = t214 * t302;
t278 = t173 * t298;
t277 = t173 * t297;
t271 = t296 * t308;
t200 = t209 * t247 - t244 * t273;
t262 = -t200 * t214 + t218 * t294;
t258 = -t187 + (t188 * t295 + t187) * t189;
t221 = t256 * qJD(2);
t194 = -t217 * qJD(4) + t222 * t247;
t185 = t211 * qJD(1) - t249 * t224 + t246 * t260;
t180 = 0.1e1 / t182;
t169 = 0.1e1 / t171;
t168 = t189 * t309;
t167 = t262 * t189;
t163 = t258 * t259;
t162 = (-t187 * t208 + t188 * t226) * t244 + t266 * t168;
t161 = t266 * t167 - t187 * t200 + t188 * t218;
t159 = t262 * t283 + (t218 * t271 - t179 * t214 + (t178 * t218 + t193 * t200 + t194 * t198) * t215) * t189;
t157 = t283 * t309 + (t261 * t286 + (t226 * t271 - t185 * t214 + (t178 * t226 + t193 * t208 + t198 * t221) * t215) * t244) * t189;
t1 = [t279 * t307 + (-t176 * t214 - t259 * t310) * t189, t157, 0, t159, 0, 0; t198 * t172 * t285 + (-t178 * t172 + (t160 * t198 + t163 * t176) * t173) * t169 - (-t163 * t173 * t285 + (-0.2e1 * t163 * t304 + (-t164 * t189 * t295 + t283) * t278 + (t279 * t308 - t164 + (t164 - t263) * t189) * t277 - t258 * t301) * t169) * t259 (-t162 * t300 - t172 * t291) * t285 + (-t162 * t301 + (-t183 * t244 + t211 * t286) * t172 + (t162 * t282 - t173 * t291) * t160 + (t226 * t286 - t157 * t198 - t168 * t178 + t221 * t244 + (-t168 * t217 - t208 * t244) * t164) * t277 + (-t208 * t286 - t157 * t217 - t168 * t193 - t185 * t244 + (t168 * t198 - t226 * t244) * t164) * t278) * t169, 0 (-t161 * t300 - t172 * t204) * t285 + (t161 * t160 * t282 + t177 * t172 + (-t204 * t160 - t161 * t176 + (-t159 * t198 - t167 * t178 + t194 + (-t167 * t217 - t200) * t164) * t297 + (-t159 * t217 - t167 * t193 - t179 + (t167 * t198 - t218) * t164) * t298) * t173) * t169, 0, 0; t204 * t284 * t292 - 0.2e1 * t200 * t280 + (t208 * t183 * t281 - t177 * t292 + t179 * t205 - t185 * t293 + t200 * t299) * t180 (t205 * t211 * t247 + t265 * t293) * t284 + (qJD(4) * t205 * t291 + (-t177 * t265 - t184 * t204) * t206 + (t265 * t281 + (-t206 * t211 + t205) * t247) * t183) * t180, 0, 0.2e1 * t259 * t280 + (t176 * t205 - t259 * t299) * t180, 0, 0;];
JaD_rot  = t1;
