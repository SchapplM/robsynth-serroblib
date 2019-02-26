% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:29:35
% EndTime: 2019-02-26 22:29:36
% DurationCPUTime: 1.43s
% Computational Cost: add. (4552->150), mult. (13562->301), div. (730->12), fcn. (17144->13), ass. (0->131)
t274 = cos(pkin(6));
t277 = sin(qJ(2));
t281 = cos(qJ(1));
t354 = sin(qJ(1));
t320 = t354 * t277;
t309 = t274 * t320;
t353 = sin(pkin(6));
t313 = qJD(3) * t353;
t318 = qJD(2) * t354;
t280 = cos(qJ(2));
t333 = t281 * t280;
t360 = -qJD(1) * t309 - t277 * t318 + (qJD(2) * t274 + qJD(1)) * t333 - t281 * t313;
t276 = sin(qJ(3));
t279 = cos(qJ(3));
t319 = t354 * t280;
t334 = t281 * t277;
t296 = -t274 * t334 - t319;
t307 = t354 * t353;
t300 = qJD(1) * t307;
t332 = qJD(3) * t279;
t221 = t360 * t276 - t279 * t300 - t296 * t332;
t315 = t281 * t353;
t294 = -t276 * t296 + t279 * t315;
t244 = t294 ^ 2;
t317 = t277 * t353;
t293 = -t274 * t279 + t276 * t317;
t256 = 0.1e1 / t293 ^ 2;
t236 = t244 * t256 + 0.1e1;
t339 = t294 * t256;
t255 = 0.1e1 / t293;
t259 = -t274 * t276 - t279 * t317;
t314 = qJD(2) * t353;
t304 = t280 * t314;
t246 = t259 * qJD(3) - t276 * t304;
t357 = t246 * t256;
t341 = t255 * t357;
t350 = (t221 * t339 + t244 * t341) / t236 ^ 2;
t359 = -0.2e1 * t350;
t295 = -t274 * t319 - t334;
t241 = t296 * qJD(1) + t295 * qJD(2);
t263 = -t309 + t333;
t290 = -t263 * t276 + t279 * t307;
t306 = qJD(1) * t315;
t220 = t290 * qJD(3) + t241 * t279 + t276 * t306;
t321 = t274 * t333;
t240 = -qJD(1) * t321 - qJD(2) * t333 + (t354 * qJD(1) + t274 * t318) * t277;
t275 = sin(qJ(4));
t278 = cos(qJ(4));
t254 = t263 * t279 + t276 * t307;
t230 = -t254 * t275 - t278 * t295;
t331 = qJD(4) * t230;
t208 = t220 * t278 - t240 * t275 + t331;
t224 = t230 ^ 2;
t231 = t254 * t278 - t275 * t295;
t226 = 0.1e1 / t231 ^ 2;
t347 = t224 * t226;
t218 = 0.1e1 + t347;
t207 = -t231 * qJD(4) - t220 * t275 - t240 * t278;
t345 = t226 * t230;
t323 = t207 * t345;
t225 = 0.1e1 / t231;
t227 = t225 * t226;
t346 = t224 * t227;
t358 = (-t208 * t346 + t323) / t218 ^ 2;
t260 = -t320 + t321;
t316 = t280 * t353;
t308 = t294 * t316;
t291 = -t255 * t260 + t256 * t308;
t356 = t276 * t291;
t222 = (qJD(3) * t296 + t300) * t276 + t360 * t279;
t237 = atan2(t294, -t293);
t232 = sin(t237);
t233 = cos(t237);
t215 = t232 * t294 - t233 * t293;
t212 = 0.1e1 / t215;
t213 = 0.1e1 / t215 ^ 2;
t355 = 0.2e1 * t290;
t245 = t290 ^ 2;
t211 = t245 * t213 + 0.1e1;
t219 = t254 * qJD(3) + t241 * t276 - t279 * t306;
t348 = t219 * t213;
t234 = 0.1e1 / t236;
t299 = -t221 * t255 - t246 * t339;
t203 = t299 * t234;
t301 = t232 * t293 + t233 * t294;
t198 = t301 * t203 + t232 * t221 + t233 * t246;
t214 = t212 * t213;
t351 = t198 * t214;
t352 = (-t245 * t351 - t290 * t348) / t211 ^ 2;
t349 = t213 * t290;
t344 = t230 * t278;
t343 = t232 * t290;
t342 = t233 * t290;
t340 = t294 * t255;
t338 = t294 * t259;
t337 = t295 * t276;
t336 = t295 * t279;
t335 = t275 * t225;
t330 = -0.2e1 * t352;
t329 = 0.2e1 * t352;
t328 = 0.2e1 * t358;
t327 = 0.2e1 * t350;
t326 = t214 * t355;
t325 = t213 * t343;
t324 = t213 * t342;
t322 = t230 * t227 * t208;
t312 = t255 * t359;
t311 = 0.2e1 * t322;
t249 = -t276 * t315 - t279 * t296;
t305 = t277 * t314;
t302 = -qJD(4) * t336 + t241;
t229 = -t249 * t278 + t260 * t275;
t228 = t249 * t275 + t260 * t278;
t298 = t226 * t344 + t335;
t297 = t249 * t255 + t256 * t338;
t292 = t232 + (-t233 * t340 - t232) * t234;
t288 = -qJD(3) * t337 + qJD(4) * t263 + t240 * t279;
t247 = t293 * qJD(3) - t279 * t304;
t242 = t295 * qJD(1) + t296 * qJD(2);
t239 = t263 * t275 + t278 * t336;
t238 = t263 * t278 - t275 * t336;
t216 = 0.1e1 / t218;
t209 = 0.1e1 / t211;
t206 = t234 * t356;
t205 = t297 * t234;
t202 = t292 * t290;
t200 = (t232 * t260 - t233 * t316) * t276 + t301 * t206;
t199 = -t301 * t205 + t232 * t249 + t233 * t259;
t197 = t297 * t327 + (-0.2e1 * t338 * t341 - t222 * t255 + (-t221 * t259 - t246 * t249 - t247 * t294) * t256) * t234;
t195 = t356 * t359 + (t291 * t332 + (0.2e1 * t308 * t341 - t242 * t255 + (t221 * t316 - t246 * t260 - t294 * t305) * t256) * t276) * t234;
t1 = [t290 * t312 + (-t219 * t255 + t290 * t357) * t234, t195, t197, 0, 0, 0; t294 * t212 * t330 + (t221 * t212 + (-t198 * t294 - t202 * t219) * t213) * t209 - (-t202 * t213 * t330 + (0.2e1 * t202 * t351 - (t203 * t234 * t340 + t327) * t325 - (-t294 * t312 + t203 + (-t203 + t299) * t234) * t324 + t292 * t348) * t209) * t290 (t200 * t349 + t212 * t337) * t329 + (t200 * t348 + (-t240 * t276 - t295 * t332) * t212 + (t200 * t326 + t213 * t337) * t198 - (t276 * t305 - t279 * t280 * t313 + t195 * t294 + t206 * t221 + (t206 * t293 + t260 * t276) * t203) * t324 - (t260 * t332 + t195 * t293 - t206 * t246 + t242 * t276 + (-t206 * t294 + t276 * t316) * t203) * t325) * t209 (t199 * t349 + t212 * t254) * t329 + (t199 * t198 * t326 - t220 * t212 + (t254 * t198 + t199 * t219 - (t197 * t294 - t205 * t221 + t247 + (-t205 * t293 + t249) * t203) * t342 - (t197 * t293 + t205 * t246 + t222 + (t205 * t294 - t259) * t203) * t343) * t213) * t209, 0, 0, 0; (-t225 * t228 + t229 * t345) * t328 + ((-t229 * qJD(4) + t222 * t275 + t242 * t278) * t225 + t229 * t311 + (-t228 * t208 - (t228 * qJD(4) - t222 * t278 + t242 * t275) * t230 - t229 * t207) * t226) * t216 (-t225 * t238 + t239 * t345) * t328 + (t239 * t311 + t302 * t225 * t278 - t288 * t335 + (-t302 * t230 * t275 - t239 * t207 - t238 * t208 - t288 * t344) * t226) * t216, t298 * t355 * t358 + (t298 * t219 - ((qJD(4) * t225 - 0.2e1 * t322) * t278 + (t207 * t278 + (-t208 - t331) * t275) * t226) * t290) * t216 (t225 * t231 + t347) * t328 + (-0.2e1 * t323 + (t226 * t231 - t225 + 0.2e1 * t346) * t208) * t216, 0, 0;];
JaD_rot  = t1;
