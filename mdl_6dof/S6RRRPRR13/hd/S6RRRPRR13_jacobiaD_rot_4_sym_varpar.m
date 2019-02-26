% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR13
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
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR13_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobiaD_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:04
% EndTime: 2019-02-26 22:23:07
% DurationCPUTime: 2.00s
% Computational Cost: add. (7376->164), mult. (23291->306), div. (681->12), fcn. (29052->15), ass. (0->155)
t301 = cos(pkin(7));
t303 = sin(qJ(3));
t394 = sin(qJ(2));
t395 = sin(qJ(1));
t297 = t395 * t394;
t302 = cos(pkin(6));
t341 = t302 * t297;
t304 = cos(qJ(2));
t305 = cos(qJ(1));
t375 = t305 * t304;
t324 = t341 - t375;
t359 = t395 * t304;
t361 = t305 * t394;
t325 = t302 * t359 + t361;
t299 = sin(pkin(6));
t363 = t299 * t395;
t393 = sin(pkin(7));
t339 = t393 * t363;
t396 = cos(qJ(3));
t263 = -t324 * t396 + (-t325 * t301 + t339) * t303;
t350 = t301 * t363;
t281 = t325 * t393 + t350;
t298 = sin(pkin(13));
t300 = cos(pkin(13));
t243 = t263 * t298 - t281 * t300;
t286 = t302 * t361 + t359;
t272 = t286 * qJD(1) + t325 * qJD(2);
t355 = t303 * t393;
t379 = t299 * t305;
t293 = t355 * t379;
t403 = -t302 * t375 + t297;
t314 = t403 * qJD(1) + t324 * qJD(2);
t313 = t314 * t303;
t319 = t325 * t396;
t331 = t396 * t339;
t401 = -t301 * t319 + t303 * t324 + t331;
t231 = qJD(1) * t293 + t401 * qJD(3) - t272 * t396 + t301 * t313;
t366 = t301 * t379;
t264 = qJD(1) * t366 - t314 * t393;
t226 = t231 * t300 + t264 * t298;
t244 = t263 * t300 + t281 * t298;
t238 = 0.1e1 / t244;
t239 = 0.1e1 / t244 ^ 2;
t387 = t226 * t238 * t239;
t352 = 0.2e1 * t243 * t387;
t347 = t396 * t393;
t333 = t347 * t379;
t320 = -t286 * t303 - t333;
t332 = t403 * t396;
t328 = t301 * t332;
t257 = t328 - t320;
t255 = t257 ^ 2;
t358 = t394 * t303;
t360 = t396 * t304;
t326 = t301 * t360 - t358;
t340 = t302 * t347;
t278 = -t326 * t299 - t340;
t276 = 0.1e1 / t278 ^ 2;
t249 = t255 * t276 + 0.1e1;
t247 = 0.1e1 / t249;
t337 = t403 * t303;
t364 = t286 * t396;
t317 = -t301 * t337 + t364;
t273 = t325 * qJD(1) + t286 * qJD(2);
t274 = -qJD(1) * t341 - qJD(2) * t297 + (qJD(2) * t302 + qJD(1)) * t375;
t362 = t301 * t396;
t321 = -qJD(1) * t331 - qJD(3) * t293 + t273 * t362 + t274 * t303;
t232 = t317 * qJD(3) + t321;
t349 = t394 * t396;
t376 = t303 * t304;
t322 = t301 * t349 + t376;
t327 = t301 * t376 + t349;
t348 = t302 * t355;
t253 = qJD(3) * t348 + (t322 * qJD(2) + t327 * qJD(3)) * t299;
t275 = 0.1e1 / t278;
t382 = t257 * t276;
t336 = -t232 * t275 + t253 * t382;
t214 = t336 * t247;
t250 = atan2(-t257, t278);
t245 = sin(t250);
t246 = cos(t250);
t342 = -t245 * t278 - t246 * t257;
t209 = t342 * t214 - t232 * t245 + t246 * t253;
t224 = -t245 * t257 + t246 * t278;
t222 = 0.1e1 / t224 ^ 2;
t405 = t209 * t222;
t404 = t253 * t276;
t256 = t401 ^ 2;
t220 = t222 * t256 + 0.1e1;
t218 = 0.1e1 / t220;
t221 = 0.1e1 / t224;
t312 = t314 * t396;
t230 = -qJD(1) * t333 + t263 * qJD(3) - t272 * t303 - t301 * t312;
t386 = t230 * t222;
t391 = t221 * t405;
t392 = (-t256 * t391 - t386 * t401) / t220 ^ 2;
t402 = -t218 * t405 - 0.2e1 * t221 * t392;
t397 = -0.2e1 * t401;
t353 = t391 * t397;
t374 = 0.2e1 * t392;
t388 = t222 * t401;
t400 = -t374 * t388 + (t353 - t386) * t218;
t399 = (qJD(1) * t339 - t286 * qJD(3) - t273 * t301) * t303 - qJD(3) * t333 + t274 * t396;
t398 = -0.2e1 * t257;
t384 = t275 * t404;
t390 = (t232 * t382 - t255 * t384) / t249 ^ 2;
t389 = t218 * t221;
t385 = t239 * t243;
t383 = t257 * t275;
t380 = t298 * t238;
t378 = t300 * t243;
t377 = t301 * t303;
t225 = t231 * t298 - t264 * t300;
t237 = t243 ^ 2;
t229 = t237 * t239 + 0.1e1;
t373 = 0.2e1 * (t225 * t385 - t237 * t387) / t229 ^ 2;
t372 = -0.2e1 * t390;
t370 = t275 * t390;
t368 = t218 * t388;
t357 = t298 * t393;
t356 = t300 * t393;
t351 = t384 * t398;
t338 = t301 * t403;
t259 = -t293 + t317;
t279 = t327 * t299 + t348;
t335 = -t259 * t275 + t279 * t382;
t269 = t286 * t362 - t337;
t285 = t322 * t299;
t334 = -t269 * t275 + t285 * t382;
t330 = -t245 + (t246 * t383 + t245) * t247;
t329 = t396 * t338;
t323 = -t301 * t358 + t360;
t318 = t303 * t338 - t364;
t270 = -t325 * t303 - t324 * t362;
t271 = t324 * t377 - t319;
t280 = -t393 * t403 + t366;
t266 = (t326 * qJD(2) + t323 * qJD(3)) * t299;
t265 = -qJD(1) * t350 - t273 * t393;
t261 = t293 + t318;
t254 = qJD(3) * t340 + (t323 * qJD(2) + t326 * qJD(3)) * t299;
t252 = t271 * t300 - t324 * t357;
t251 = t271 * t298 + t324 * t356;
t242 = t261 * t300 + t280 * t298;
t241 = t261 * t298 - t280 * t300;
t236 = t274 * t362 - t273 * t303 + (-t286 * t377 - t332) * qJD(3);
t235 = -t270 * qJD(3) + t272 * t377 + t312;
t234 = qJD(3) * t329 - t399;
t233 = -qJD(3) * t328 + t399;
t227 = 0.1e1 / t229;
t217 = t334 * t247;
t216 = t335 * t247;
t210 = t342 * t216 - t245 * t259 + t246 * t279;
t208 = t334 * t372 + (t285 * t351 - t236 * t275 + (t232 * t285 + t253 * t269 + t257 * t266) * t276) * t247;
t207 = t335 * t372 + (t279 * t351 - t233 * t275 + (t232 * t279 + t253 * t259 + t254 * t257) * t276) * t247;
t1 = [t370 * t397 + (-t230 * t275 - t401 * t404) * t247, t208, t207, 0, 0, 0; (t318 * qJD(3) - t321) * t389 + (t330 * t230 - ((-t214 * t247 * t383 + t372) * t245 + (t370 * t398 - t214 + (t214 - t336) * t247) * t246) * t401) * t368 + t402 * (-t329 + t320) - t400 * t330 * t401 (t271 * qJD(3) - t272 * t362 + t313) * t389 + ((-t208 * t257 - t217 * t232 + t266 + (-t217 * t278 - t269) * t214) * t246 + (-t208 * t278 - t217 * t253 - t236 + (t217 * t257 - t285) * t214) * t245) * t368 + t402 * t270 + t400 * (t342 * t217 - t245 * t269 + t246 * t285) (-t210 * t388 - t221 * t263) * t374 + (t210 * t353 + t231 * t221 + (-t263 * t209 - t210 * t230 - (-(-t207 * t257 - t216 * t232 + t254 + (-t216 * t278 - t259) * t214) * t246 - (-t207 * t278 - t216 * t253 - t233 + (t216 * t257 - t279) * t214) * t245) * t401) * t222) * t218, 0, 0, 0; (-t238 * t241 + t242 * t385) * t373 + ((t234 * t298 - t265 * t300) * t238 + t242 * t352 + (-t241 * t226 - (t234 * t300 + t265 * t298) * t243 - t242 * t225) * t239) * t227 (-t238 * t251 + t252 * t385) * t373 + ((t235 * t298 + t272 * t356) * t238 + t252 * t352 + (-t251 * t226 - (t235 * t300 - t272 * t357) * t243 - t252 * t225) * t239) * t227 -(-t239 * t378 + t380) * t401 * t373 + (t401 * t300 * t352 - t230 * t380 + (t230 * t378 - (t225 * t300 + t226 * t298) * t401) * t239) * t227, 0, 0, 0;];
JaD_rot  = t1;
