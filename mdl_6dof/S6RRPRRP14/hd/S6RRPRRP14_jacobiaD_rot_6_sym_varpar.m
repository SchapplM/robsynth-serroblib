% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP14_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:53:25
% EndTime: 2019-02-26 21:53:28
% DurationCPUTime: 2.65s
% Computational Cost: add. (9833->178), mult. (28394->345), div. (951->12), fcn. (36086->13), ass. (0->149)
t305 = cos(pkin(6));
t309 = sin(qJ(1));
t312 = cos(qJ(2));
t377 = t309 * t312;
t308 = sin(qJ(2));
t313 = cos(qJ(1));
t378 = t308 * t313;
t298 = t305 * t377 + t378;
t333 = t305 * t378 + t377;
t276 = t298 * qJD(1) + t333 * qJD(2);
t376 = t312 * t313;
t355 = t305 * t376;
t380 = t308 * t309;
t297 = -t355 + t380;
t307 = sin(qJ(4));
t311 = cos(qJ(4));
t304 = sin(pkin(6));
t382 = t304 * t313;
t287 = t297 * t311 + t307 * t382;
t374 = qJD(1) * t304;
t354 = t309 * t374;
t250 = t287 * qJD(4) + t276 * t307 + t311 * t354;
t289 = -t297 * t307 + t311 * t382;
t306 = sin(qJ(5));
t310 = cos(qJ(5));
t264 = t289 * t306 + t333 * t310;
t346 = qJD(2) * t305 + qJD(1);
t356 = t305 * t380;
t373 = qJD(2) * t308;
t323 = -qJD(1) * t356 - t309 * t373 + t346 * t376;
t236 = t264 * qJD(5) + t250 * t310 + t323 * t306;
t328 = t333 * t306;
t265 = t289 * t310 - t328;
t418 = t265 * qJD(5) - t250 * t306 + t323 * t310;
t261 = t264 ^ 2;
t383 = t304 * t312;
t334 = -t305 * t311 + t307 * t383;
t379 = t308 * t310;
t335 = t304 * t379 + t306 * t334;
t278 = 0.1e1 / t335 ^ 2;
t246 = t261 * t278 + 0.1e1;
t242 = 0.1e1 / t246;
t295 = -t305 * t307 - t311 * t383;
t352 = t304 * t373;
t281 = t295 * qJD(4) + t307 * t352;
t381 = t306 * t308;
t284 = t304 * t381 - t310 * t334;
t372 = qJD(2) * t312;
t351 = t304 * t372;
t254 = t284 * qJD(5) + t281 * t306 - t310 * t351;
t277 = 0.1e1 / t335;
t386 = t264 * t278;
t339 = -t254 * t386 - t277 * t418;
t222 = t339 * t242;
t247 = atan2(t264, -t335);
t240 = sin(t247);
t241 = cos(t247);
t340 = t240 * t335 + t241 * t264;
t217 = t340 * t222 + t240 * t418 + t241 * t254;
t234 = t240 * t264 - t241 * t335;
t232 = 0.1e1 / t234 ^ 2;
t417 = t217 * t232;
t274 = -qJD(1) * t355 - t313 * t372 + t346 * t380;
t384 = t304 * t309;
t285 = t298 * t311 - t307 * t384;
t353 = t313 * t374;
t253 = t285 * qJD(4) - t274 * t307 + t311 * t353;
t286 = t298 * t307 + t311 * t384;
t332 = t356 - t376;
t327 = t332 * t306;
t263 = t286 * t310 - t327;
t275 = -t333 * qJD(1) - t298 * qJD(2);
t237 = t263 * qJD(5) + t253 * t306 - t275 * t310;
t326 = t332 * t310;
t262 = t286 * t306 + t326;
t403 = 0.2e1 * t262;
t231 = 0.1e1 / t234;
t412 = t231 * t417;
t349 = t403 * t412;
t416 = -t232 * t237 + t349;
t259 = 0.1e1 / t263 ^ 2;
t280 = t285 ^ 2;
t389 = t259 * t280;
t248 = 0.1e1 + t389;
t252 = -t286 * qJD(4) - t274 * t311 - t307 * t353;
t238 = -t262 * qJD(5) + t253 * t310 + t275 * t306;
t258 = 0.1e1 / t263;
t395 = t238 * t258 * t259;
t360 = t280 * t395;
t388 = t259 * t285;
t400 = (t252 * t388 - t360) / t248 ^ 2;
t413 = 0.2e1 * t400;
t411 = -0.2e1 * t262;
t410 = t254 * t278;
t336 = t277 * t287 - t295 * t386;
t409 = t306 * t336;
t244 = 0.1e1 / t248;
t391 = t244 * t259;
t407 = t238 * t391 + t258 * t413;
t257 = t262 ^ 2;
t230 = t232 * t257 + 0.1e1;
t228 = 0.1e1 / t230;
t396 = t232 * t262;
t401 = (t237 * t396 - t257 * t412) / t230 ^ 2;
t406 = -t228 * t417 - 0.2e1 * t231 * t401;
t345 = t388 * t400;
t359 = t285 * t395;
t405 = 0.2e1 * t244 * t359 - t252 * t391 + 0.2e1 * t345;
t368 = 0.2e1 * t401;
t404 = t416 * t228 + t368 * t396;
t249 = t289 * qJD(4) + t276 * t311 - t307 * t354;
t402 = 0.2e1 * t264;
t390 = t277 * t410;
t399 = (t261 * t390 + t386 * t418) / t246 ^ 2;
t398 = t228 * t231;
t394 = t240 * t262;
t393 = t241 * t262;
t392 = t244 * t258;
t387 = t264 * t277;
t385 = t285 * t306;
t371 = qJD(4) * t311;
t370 = qJD(5) * t307;
t369 = qJD(5) * t310;
t367 = -0.2e1 * t399;
t364 = t277 * t399;
t362 = t228 * t396;
t357 = t244 * t388;
t347 = t390 * t402;
t338 = -t265 * t277 - t284 * t386;
t269 = t297 * t310 + t307 * t328;
t291 = (t307 * t381 - t310 * t312) * t304;
t337 = t269 * t277 - t291 * t386;
t325 = -t240 + (t241 * t387 + t240) * t242;
t324 = qJD(4) * t332;
t321 = -t332 * t370 - t274;
t320 = -t298 * qJD(5) + t275 * t307 - t311 * t324;
t282 = t334 * qJD(4) + t311 * t352;
t256 = ((qJD(2) + t370) * t379 + (t308 * t371 + (qJD(2) * t307 + qJD(5)) * t312) * t306) * t304;
t255 = t335 * qJD(5) + t281 * t310 + t306 * t351;
t239 = (t333 * t370 + t276) * t310 + (-t297 * qJD(5) + t323 * t307 + t333 * t371) * t306;
t227 = t242 * t409;
t226 = t337 * t242;
t225 = t338 * t242;
t220 = (-t240 * t287 + t241 * t295) * t306 + t340 * t227;
t218 = t340 * t225 + t240 * t265 + t241 * t284;
t216 = t337 * t367 + (-t291 * t347 + t239 * t277 + (t254 * t269 - t256 * t264 - t291 * t418) * t278) * t242;
t214 = t338 * t367 + (-t284 * t347 + t236 * t277 + (-t254 * t265 - t255 * t264 - t284 * t418) * t278) * t242;
t213 = t367 * t409 + (t336 * t369 + (-t295 * t347 + t249 * t277 + (t254 * t287 - t264 * t282 - t295 * t418) * t278) * t306) * t242;
t1 = [-t364 * t403 + (t237 * t277 + t262 * t410) * t242, t216, 0, t213, t214, 0; t418 * t398 - (t325 * t237 + ((-t222 * t242 * t387 + t367) * t240 + (-t364 * t402 - t222 + (t222 - t339) * t242) * t241) * t262) * t362 + t406 * t264 + t404 * t325 * t262 (t320 * t306 + t321 * t310) * t398 - ((t216 * t264 + t226 * t418 + t256 + (t226 * t335 - t269) * t222) * t241 + (t216 * t335 - t226 * t254 - t239 + (-t226 * t264 - t291) * t222) * t240) * t362 + t406 * (t298 * t310 - t307 * t327) + t404 * (t340 * t226 - t240 * t269 + t241 * t291) 0 (t220 * t396 - t231 * t385) * t368 + ((t252 * t306 + t285 * t369) * t231 + t416 * t220 + (-t385 * t217 - (t295 * t369 + t213 * t264 + t227 * t418 + t282 * t306 + (t227 * t335 - t287 * t306) * t222) * t393 - (-t287 * t369 + t213 * t335 - t227 * t254 - t249 * t306 + (-t227 * t264 - t295 * t306) * t222) * t394) * t232) * t228 (t218 * t396 - t231 * t263) * t368 + (t218 * t349 + t238 * t231 + (-t263 * t217 - t218 * t237 - (t214 * t264 + t225 * t418 + t255 + (t225 * t335 + t265) * t222) * t393 - (t214 * t335 - t225 * t254 - t236 + (-t225 * t264 - t284) * t222) * t394) * t232) * t228, 0; t236 * t357 - t249 * t392 + t405 * t265 + t407 * t287 -(-t321 * t306 + t320 * t310) * t357 + (t275 * t311 + t307 * t324) * t392 + t407 * t311 * t332 + t405 * (-t298 * t306 - t307 * t326) 0 (t258 * t286 + t310 * t389) * t413 + (0.2e1 * t310 * t360 - t253 * t258 + (qJD(5) * t280 * t306 - 0.2e1 * t252 * t285 * t310 + t238 * t286) * t259) * t244, t345 * t411 + (t359 * t411 + (t237 * t285 + t252 * t262) * t259) * t244, 0;];
JaD_rot  = t1;
