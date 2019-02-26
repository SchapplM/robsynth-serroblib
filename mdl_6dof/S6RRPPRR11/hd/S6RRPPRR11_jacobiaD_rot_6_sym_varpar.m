% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR11_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:20
% EndTime: 2019-02-26 21:34:21
% DurationCPUTime: 1.86s
% Computational Cost: add. (8428->148), mult. (13478->295), div. (726->12), fcn. (17045->13), ass. (0->123)
t266 = sin(qJ(2));
t267 = sin(qJ(1));
t269 = cos(qJ(2));
t270 = cos(qJ(1));
t348 = cos(pkin(6));
t300 = t270 * t348;
t285 = -t266 * t300 - t267 * t269;
t301 = t267 * t348;
t287 = t270 * t266 + t269 * t301;
t264 = sin(pkin(6));
t326 = t264 * t270;
t361 = t287 * qJD(1) - t285 * qJD(2) + qJD(5) * t326;
t252 = t267 * t266 - t269 * t300;
t263 = pkin(11) + qJ(5);
t261 = sin(t263);
t262 = cos(t263);
t289 = t252 * t262 + t261 * t326;
t237 = t289 ^ 2;
t327 = t264 * t269;
t286 = -t348 * t261 - t262 * t327;
t247 = 0.1e1 / t286 ^ 2;
t226 = t237 * t247 + 0.1e1;
t220 = 0.1e1 / t226;
t250 = -t261 * t327 + t348 * t262;
t322 = qJD(2) * t266;
t305 = t264 * t322;
t235 = t250 * qJD(5) - t262 * t305;
t246 = 0.1e1 / t286;
t334 = t289 * t247;
t323 = qJD(1) * t267;
t349 = t252 * qJD(5) + t264 * t323;
t359 = -t261 * t349 + t361 * t262;
t293 = -t235 * t334 - t246 * t359;
t193 = t293 * t220;
t227 = atan2(t289, -t286);
t214 = sin(t227);
t215 = cos(t227);
t295 = t214 * t286 + t215 * t289;
t188 = t295 * t193 + t214 * t359 + t215 * t235;
t205 = t214 * t289 - t215 * t286;
t203 = 0.1e1 / t205 ^ 2;
t360 = t188 * t203;
t210 = t361 * t261 + t262 * t349;
t328 = t264 * t267;
t239 = t287 * t261 + t262 * t328;
t296 = t266 * t301;
t325 = t270 * t269;
t254 = -t296 + t325;
t265 = sin(qJ(6));
t268 = cos(qJ(6));
t222 = t239 * t265 - t254 * t268;
t358 = 0.2e1 * t222;
t202 = 0.1e1 / t205;
t357 = t202 * t360;
t281 = (t348 * qJD(1) + qJD(2)) * t325 - qJD(2) * t296 - t266 * t323;
t306 = qJD(1) * t326;
t212 = t239 * qJD(5) + t261 * t306 - t281 * t262;
t350 = -t261 * t328 + t287 * t262;
t299 = -0.2e1 * t350 * t357;
t356 = -t203 * t212 + t299;
t354 = t235 * t247;
t329 = t264 * t266;
t309 = t289 * t329;
t288 = t246 * t285 + t247 * t309;
t353 = t262 * t288;
t352 = -t254 * t261 * qJD(6) - t281;
t332 = t254 * t262;
t351 = qJD(5) * t332 - qJD(6) * t287;
t223 = t239 * t268 + t254 * t265;
t217 = 0.1e1 / t223;
t218 = 0.1e1 / t223 ^ 2;
t236 = t350 ^ 2;
t199 = t236 * t203 + 0.1e1;
t341 = t203 * t350;
t347 = (-t212 * t341 - t236 * t357) / t199 ^ 2;
t213 = qJD(5) * t350 + t281 * t261 + t262 * t306;
t232 = t285 * qJD(1) - t287 * qJD(2);
t200 = t223 * qJD(6) + t213 * t265 - t232 * t268;
t216 = t222 ^ 2;
t208 = t216 * t218 + 0.1e1;
t337 = t218 * t222;
t319 = qJD(6) * t222;
t201 = t213 * t268 + t232 * t265 - t319;
t343 = t201 * t217 * t218;
t346 = (t200 * t337 - t216 * t343) / t208 ^ 2;
t336 = t246 * t354;
t344 = (t237 * t336 + t334 * t359) / t226 ^ 2;
t206 = 0.1e1 / t208;
t340 = t206 * t218;
t339 = t214 * t350;
t338 = t215 * t350;
t335 = t289 * t246;
t333 = t289 * t250;
t331 = t261 * t265;
t330 = t261 * t268;
t321 = qJD(2) * t269;
t320 = qJD(5) * t261;
t317 = 0.2e1 * t347;
t316 = -0.2e1 * t346;
t315 = -0.2e1 * t344;
t314 = 0.2e1 * t344;
t312 = t218 * t346;
t311 = t200 * t340;
t310 = t222 * t343;
t298 = t246 * t314;
t297 = 0.2e1 * t310;
t290 = -t252 * t261 + t262 * t326;
t225 = t265 * t285 + t268 * t290;
t224 = t265 * t290 - t268 * t285;
t292 = -t265 * t217 + t268 * t337;
t291 = t246 * t290 + t247 * t333;
t284 = -t214 + (t215 * t335 + t214) * t220;
t234 = t286 * qJD(5) + t261 * t305;
t233 = -qJD(1) * t296 - t267 * t322 + (qJD(2) * t348 + qJD(1)) * t325;
t229 = t254 * t330 - t287 * t265;
t197 = 0.1e1 / t199;
t196 = t220 * t353;
t194 = t291 * t220;
t190 = (-t214 * t285 - t215 * t329) * t262 + t295 * t196;
t189 = -t295 * t194 + t214 * t290 + t215 * t250;
t187 = t291 * t314 + (-0.2e1 * t333 * t336 + t210 * t246 + (-t234 * t289 - t235 * t290 - t250 * t359) * t247) * t220;
t185 = t315 * t353 + (-t288 * t320 + (0.2e1 * t309 * t336 - t233 * t246 + (t235 * t285 + (t266 * t359 + t289 * t321) * t264) * t247) * t262) * t220;
t1 = [t350 * t298 + (t212 * t246 - t350 * t354) * t220, t185, 0, 0, t187, 0; -0.2e1 * t289 * t202 * t347 + (t359 * t202 - t289 * t360 + (t284 * t212 - ((-t193 * t220 * t335 + t315) * t214 + (-t289 * t298 - t193 + (t193 - t293) * t220) * t215) * t350) * t341) * t197 - (t356 * t197 - t341 * t317) * t284 * t350 (-t190 * t341 + t202 * t332) * t317 + ((-t232 * t262 + t254 * t320) * t202 + t356 * t190 + (t332 * t188 + (t185 * t289 + t196 * t359 + (-t262 * t321 + t266 * t320) * t264 + (t196 * t286 - t262 * t285) * t193) * t338 + (t285 * t320 + t185 * t286 - t196 * t235 + t233 * t262 + (-t196 * t289 + t262 * t329) * t193) * t339) * t203) * t197, 0, 0 (-t189 * t341 - t202 * t239) * t317 + (t189 * t299 + t213 * t202 + (-t239 * t188 - t189 * t212 + (t187 * t289 - t194 * t359 + t234 + (-t194 * t286 + t290) * t193) * t338 + (t187 * t286 + t194 * t235 - t210 + (t194 * t289 - t250) * t193) * t339) * t203) * t197, 0; 0.2e1 * (-t217 * t224 + t225 * t337) * t346 + ((t225 * qJD(6) - t210 * t265 + t233 * t268) * t217 + t225 * t297 + (-t224 * t201 - (-t224 * qJD(6) - t210 * t268 - t233 * t265) * t222 - t225 * t200) * t218) * t206 (t312 * t358 - t311) * t229 + (-t201 * t340 + t217 * t316) * (t254 * t331 + t287 * t268) + (t229 * t297 + (t331 * t217 - t330 * t337) * t232 + (-t217 * t352 - t337 * t351) * t268 + (t217 * t351 - t337 * t352) * t265) * t206, 0, 0, -t292 * t350 * t316 + (t292 * t212 - ((-qJD(6) * t217 - 0.2e1 * t310) * t268 + (t200 * t268 + (t201 - t319) * t265) * t218) * t350) * t206, t316 + (t311 + (-t206 * t343 - t312) * t222) * t358;];
JaD_rot  = t1;
