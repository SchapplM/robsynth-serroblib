% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR9_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobiaD_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:28
% EndTime: 2019-02-26 22:20:29
% DurationCPUTime: 0.96s
% Computational Cost: add. (3107->118), mult. (9489->259), div. (449->12), fcn. (12082->15), ass. (0->115)
t241 = sin(pkin(13));
t244 = cos(pkin(13));
t246 = sin(qJ(3));
t249 = cos(qJ(3));
t231 = t246 * t241 - t249 * t244;
t242 = sin(pkin(7));
t217 = t231 * t242;
t310 = cos(pkin(7));
t276 = t249 * t310;
t277 = t246 * t310;
t219 = -t241 * t277 + t244 * t276;
t245 = cos(pkin(6));
t250 = cos(qJ(2));
t251 = cos(qJ(1));
t288 = t251 * t250;
t247 = sin(qJ(2));
t248 = sin(qJ(1));
t292 = t248 * t247;
t264 = t245 * t292 - t288;
t289 = t251 * t247;
t291 = t248 * t250;
t266 = t245 * t291 + t289;
t267 = t249 * t241 + t246 * t244;
t243 = sin(pkin(6));
t295 = t243 * t248;
t190 = -t217 * t295 - t219 * t266 + t264 * t267;
t183 = t190 ^ 2;
t218 = t267 * t242;
t258 = -t241 * t276 - t244 * t277;
t260 = t218 * t295 + t231 * t264 + t258 * t266;
t185 = 0.1e1 / t260 ^ 2;
t312 = t183 * t185;
t297 = t242 * t250;
t224 = -t243 * t297 + t245 * t310;
t221 = 0.1e1 / t224;
t265 = t245 * t289 + t291;
t296 = t243 * t247;
t222 = 0.1e1 / t224 ^ 2;
t227 = -t245 * t288 + t292;
t278 = t243 * t310;
t262 = -t227 * t242 + t251 * t278;
t299 = t262 * t222;
t311 = t242 * (t221 * t265 + t296 * t299);
t202 = atan2(t262, t224);
t197 = sin(t202);
t198 = cos(t202);
t182 = t197 * t262 + t198 * t224;
t179 = 0.1e1 / t182;
t184 = 0.1e1 / t260;
t180 = 0.1e1 / t182 ^ 2;
t211 = -t242 * t266 - t248 * t278;
t208 = t211 ^ 2;
t177 = t208 * t180 + 0.1e1;
t203 = t227 * qJD(1) + t264 * qJD(2);
t270 = qJD(1) * t278;
t195 = t203 * t242 - t251 * t270;
t303 = t195 * t180;
t207 = t262 ^ 2;
t201 = t207 * t222 + 0.1e1;
t199 = 0.1e1 / t201;
t205 = t266 * qJD(1) + t265 * qJD(2);
t196 = -t205 * t242 - t248 * t270;
t285 = qJD(2) * t243;
t298 = t242 * t247;
t273 = t285 * t298;
t269 = t222 * t273;
t259 = t196 * t221 - t262 * t269;
t170 = t259 * t199;
t268 = -t197 * t224 + t198 * t262;
t165 = t268 * t170 + t197 * t196 + t198 * t273;
t308 = t165 * t179 * t180;
t309 = (-t208 * t308 + t211 * t303) / t177 ^ 2;
t204 = t265 * qJD(1) + t266 * qJD(2);
t213 = qJD(3) * t217;
t215 = t219 * qJD(3);
t226 = t267 * qJD(3);
t286 = qJD(1) * t251;
t169 = -t203 * t258 + t204 * t231 - t266 * t215 + t264 * t226 + (-t213 * t248 + t218 * t286) * t243;
t186 = t184 * t185;
t307 = t169 * t186;
t223 = t221 * t222;
t306 = (-t207 * t223 * t273 + t196 * t299) / t201 ^ 2;
t305 = t180 * t211;
t304 = t185 * t190;
t302 = t197 * t211;
t301 = t198 * t211;
t300 = t262 * t221;
t294 = t243 * t251;
t287 = qJD(1) * t248;
t283 = -0.2e1 * t309;
t174 = 0.1e1 + t312;
t214 = t242 * t226;
t216 = t258 * qJD(3);
t225 = t231 * qJD(3);
t168 = t203 * t219 + t204 * t267 - t266 * t216 - t264 * t225 + (-t214 * t248 - t217 * t286) * t243;
t279 = t168 * t304;
t282 = 0.2e1 * (-t183 * t307 + t279) / t174 ^ 2;
t281 = -0.2e1 * t308;
t280 = 0.2e1 * t306;
t275 = -0.2e1 * t190 * t307;
t274 = -0.2e1 * t221 * t306;
t261 = t197 + (t198 * t300 - t197) * t199;
t240 = t242 ^ 2;
t206 = t264 * qJD(1) + t227 * qJD(2);
t194 = t231 * t266 - t258 * t264;
t193 = -t219 * t264 - t266 * t267;
t188 = t218 * t294 - t227 * t258 + t231 * t265;
t187 = t217 * t294 - t227 * t219 - t265 * t267;
t175 = 0.1e1 / t177;
t172 = 0.1e1 / t174;
t171 = t199 * t311;
t167 = t261 * t211;
t166 = (-t197 * t265 + t198 * t296) * t242 - t268 * t171;
t163 = t280 * t311 + (t206 * t221 * t242 + (-t196 * t222 * t298 + (t222 * t265 * t240 * t247 + (0.2e1 * t223 * t240 * t243 * t247 ^ 2 - t222 * t297) * t262) * qJD(2)) * t243) * t199;
t1 = [t211 * t274 + (t195 * t221 - t211 * t269) * t199, t163, 0, 0, 0, 0; t262 * t179 * t283 + (t196 * t179 + (-t165 * t262 + t167 * t195) * t180) * t175 + ((t167 * t281 + t261 * t303) * t175 + (t167 * t283 + ((-t170 * t199 * t300 + t280) * t302 + (t262 * t274 + t170 + (-t170 + t259) * t199) * t301) * t175) * t180) * t211, 0.2e1 * (t179 * t242 * t264 - t166 * t305) * t309 + ((t268 * t163 - (-t182 * t170 + t196 * t198) * t171) * t305 + (t211 * t281 + t303) * t166 + (-t204 * t179 + (t264 * t165 + (-t170 * t265 + t250 * t285) * t301 + (t206 + (qJD(2) * t171 - t170) * t296) * t302) * t180) * t242) * t175, 0, 0, 0, 0; (-t184 * t187 - t188 * t304) * t282 + ((-t205 * t219 + t206 * t267 - t227 * t216 + t265 * t225 + (t214 * t251 - t217 * t287) * t243) * t184 + t188 * t275 + (-t187 * t169 + (-t205 * t258 - t206 * t231 + t227 * t215 + t265 * t226 + (-t213 * t251 - t218 * t287) * t243) * t190 + t188 * t168) * t185) * t172 (-t184 * t193 - t194 * t304) * t282 + ((t203 * t267 - t204 * t219 - t216 * t264 + t225 * t266) * t184 + t194 * t275 + (-t193 * t169 + (-t203 * t231 - t204 * t258 + t215 * t264 + t226 * t266) * t190 + t194 * t168) * t185) * t172 (-t184 * t260 - t312) * t282 + (0.2e1 * t279 + (-0.2e1 * t183 * t186 - t185 * t260 + t184) * t169) * t172, 0, 0, 0;];
JaD_rot  = t1;
