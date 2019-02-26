% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR14_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:45
% EndTime: 2019-02-26 22:23:46
% DurationCPUTime: 1.32s
% Computational Cost: add. (4522->150), mult. (13478->299), div. (726->12), fcn. (17045->13), ass. (0->128)
t258 = sin(qJ(1));
t254 = cos(pkin(6));
t279 = qJD(2) * t254 + qJD(1);
t257 = sin(qJ(2));
t304 = t258 * t257;
t287 = t254 * t304;
t299 = qJD(2) * t257;
t261 = cos(qJ(2));
t262 = cos(qJ(1));
t301 = t261 * t262;
t253 = sin(pkin(6));
t308 = t253 * t262;
t333 = -qJD(1) * t287 - qJD(3) * t308 - t258 * t299 + t279 * t301;
t243 = -t287 + t301;
t256 = sin(qJ(3));
t260 = cos(qJ(3));
t310 = t253 * t260;
t231 = t243 * t256 - t258 * t310;
t303 = t258 * t261;
t305 = t257 * t262;
t242 = t254 * t303 + t305;
t255 = sin(qJ(5));
t259 = cos(qJ(5));
t275 = t231 * t259 - t242 * t255;
t332 = t275 * qJD(5);
t241 = t254 * t305 + t303;
t306 = t256 * t262;
t227 = t241 * t260 - t253 * t306;
t239 = t254 * t256 + t257 * t310;
t214 = atan2(-t227, t239);
t209 = sin(t214);
t210 = cos(t214);
t192 = -t209 * t227 + t210 * t239;
t190 = 0.1e1 / t192 ^ 2;
t311 = t253 * t256;
t232 = t243 * t260 + t258 * t311;
t223 = t232 ^ 2;
t188 = t190 * t223 + 0.1e1;
t218 = -t241 * qJD(1) - t242 * qJD(2);
t296 = qJD(3) * t260;
t297 = qJD(3) * t256;
t197 = t218 * t260 - t243 * t297 + (qJD(1) * t306 + t258 * t296) * t253;
t321 = t197 * t190;
t222 = t227 ^ 2;
t236 = 0.1e1 / t239 ^ 2;
t213 = t222 * t236 + 0.1e1;
t211 = 0.1e1 / t213;
t280 = t333 * t260;
t300 = qJD(1) * t253;
t286 = t258 * t300;
t199 = -t241 * t297 + t256 * t286 + t280;
t238 = t254 * t260 - t257 * t311;
t298 = qJD(2) * t261;
t284 = t253 * t298;
t225 = t238 * qJD(3) + t260 * t284;
t235 = 0.1e1 / t239;
t315 = t227 * t236;
t274 = -t199 * t235 + t225 * t315;
t180 = t274 * t211;
t277 = -t209 * t239 - t210 * t227;
t175 = t277 * t180 - t199 * t209 + t210 * t225;
t189 = 0.1e1 / t192;
t191 = t189 * t190;
t326 = t175 * t191;
t295 = 0.2e1 * (-t223 * t326 + t232 * t321) / t188 ^ 2;
t331 = t225 * t236;
t288 = t254 * t301;
t240 = t288 - t304;
t309 = t253 * t261;
t271 = -t235 * t240 + t309 * t315;
t330 = t260 * t271;
t313 = t242 * t259;
t208 = t231 * t255 + t313;
t202 = 0.1e1 / t208;
t203 = 0.1e1 / t208 ^ 2;
t329 = -0.2e1 * t227;
t328 = 0.2e1 * t232;
t285 = t260 * t300;
t196 = t232 * qJD(3) + t218 * t256 - t262 * t285;
t217 = -qJD(1) * t288 - t262 * t298 + t279 * t304;
t184 = t208 * qJD(5) - t196 * t259 - t217 * t255;
t201 = t275 ^ 2;
t195 = t201 * t203 + 0.1e1;
t320 = t203 * t275;
t185 = t196 * t255 - t217 * t259 + t332;
t323 = t185 * t202 * t203;
t325 = (-t184 * t320 - t201 * t323) / t195 ^ 2;
t317 = t235 * t331;
t324 = (t199 * t315 - t222 * t317) / t213 ^ 2;
t322 = t190 * t232;
t319 = t209 * t232;
t318 = t210 * t232;
t316 = t227 * t235;
t314 = t242 * t256;
t312 = t242 * t260;
t307 = t255 * t275;
t302 = t259 * t202;
t294 = 0.2e1 * t325;
t293 = -0.2e1 * t324;
t292 = t191 * t328;
t291 = t235 * t324;
t290 = t190 * t319;
t289 = t190 * t318;
t282 = -0.2e1 * t275 * t323;
t281 = t317 * t329;
t278 = -qJD(5) * t314 + t218;
t226 = t241 * t256 + t260 * t308;
t276 = -t226 * t259 - t240 * t255;
t206 = -t226 * t255 + t240 * t259;
t273 = -t203 * t307 + t302;
t272 = t226 * t235 + t238 * t315;
t270 = -t209 + (t210 * t316 + t209) * t211;
t198 = t241 * t296 + t256 * t333 - t258 * t285;
t269 = qJD(5) * t243 - t217 * t256 + t242 * t296;
t224 = -t239 * qJD(3) - t256 * t284;
t219 = -t242 * qJD(1) - t241 * qJD(2);
t216 = t243 * t259 - t255 * t314;
t215 = t243 * t255 + t256 * t313;
t193 = 0.1e1 / t195;
t186 = 0.1e1 / t188;
t183 = t211 * t330;
t182 = t272 * t211;
t179 = t270 * t232;
t177 = (-t209 * t240 + t210 * t309) * t260 + t277 * t183;
t176 = t277 * t182 + t209 * t226 + t210 * t238;
t174 = t272 * t293 + (t238 * t281 + t198 * t235 + (t199 * t238 + t224 * t227 - t225 * t226) * t236) * t211;
t172 = t293 * t330 + (-t271 * t297 + (t281 * t309 - t219 * t235 + (t225 * t240 + (t199 * t261 - t227 * t299) * t253) * t236) * t260) * t211;
t1 = [t291 * t328 + (-t197 * t235 + t232 * t331) * t211, t172, t174, 0, 0, 0; t227 * t189 * t295 + (((qJD(3) * t241 - t286) * t256 - t280) * t189 + (t175 * t227 - t179 * t197) * t190) * t186 + (t179 * t190 * t295 + (0.2e1 * t179 * t326 - (-t180 * t211 * t316 + t293) * t290 - (t291 * t329 - t180 + (t180 - t274) * t211) * t289 - t270 * t321) * t186) * t232 (t177 * t322 + t189 * t312) * t295 + (-t177 * t321 + (t217 * t260 + t242 * t297) * t189 + (t177 * t292 + t190 * t312) * t175 - (-t172 * t227 - t183 * t199 + (-t260 * t299 - t261 * t297) * t253 + (-t183 * t239 - t240 * t260) * t180) * t289 - (t240 * t297 - t172 * t239 - t183 * t225 - t219 * t260 + (t183 * t227 - t260 * t309) * t180) * t290) * t186 (t176 * t322 + t189 * t231) * t295 + (t176 * t175 * t292 - t196 * t189 + (t231 * t175 - t176 * t197 - (-t174 * t227 - t182 * t199 + t224 + (-t182 * t239 + t226) * t180) * t318 - (-t174 * t239 - t182 * t225 + t198 + (t182 * t227 - t238) * t180) * t319) * t190) * t186, 0, 0, 0; (t202 * t276 - t206 * t320) * t294 + ((t206 * qJD(5) + t198 * t259 + t219 * t255) * t202 + t206 * t282 + (t276 * t185 + (t276 * qJD(5) - t198 * t255 + t219 * t259) * t275 - t206 * t184) * t203) * t193 (-t202 * t215 - t216 * t320) * t294 + (t216 * t282 + t278 * t202 * t255 + t269 * t302 + (t259 * t275 * t278 - t216 * t184 - t215 * t185 - t269 * t307) * t203) * t193, t273 * t232 * t294 + (-t273 * t197 + ((qJD(5) * t202 + t282) * t255 + (-t184 * t255 + (t185 + t332) * t259) * t203) * t232) * t193, 0, -0.2e1 * t325 - 0.2e1 * (t184 * t203 * t193 - (-t193 * t323 - t203 * t325) * t275) * t275, 0;];
JaD_rot  = t1;
