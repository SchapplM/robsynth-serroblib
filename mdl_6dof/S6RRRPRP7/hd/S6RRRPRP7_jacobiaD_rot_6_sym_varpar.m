% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRP7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:12:39
% EndTime: 2019-02-26 22:12:43
% DurationCPUTime: 2.92s
% Computational Cost: add. (15414->179), mult. (28394->346), div. (951->12), fcn. (36086->13), ass. (0->147)
t314 = cos(pkin(6));
t317 = sin(qJ(1));
t316 = sin(qJ(2));
t375 = qJD(2) * t316;
t353 = t317 * t375;
t376 = qJD(1) * t317;
t356 = t316 * t376;
t319 = cos(qJ(2));
t320 = cos(qJ(1));
t377 = t319 * t320;
t282 = -t314 * t356 - t353 + (qJD(2) * t314 + qJD(1)) * t377;
t379 = t317 * t319;
t380 = t316 * t320;
t303 = t314 * t380 + t379;
t312 = qJ(3) + pkin(11);
t310 = sin(t312);
t311 = cos(t312);
t313 = sin(pkin(6));
t383 = t313 * t320;
t293 = t303 * t310 + t311 * t383;
t358 = t313 * t376;
t259 = t293 * qJD(3) - t282 * t311 - t310 * t358;
t294 = -t303 * t311 + t310 * t383;
t315 = sin(qJ(5));
t318 = cos(qJ(5));
t381 = t316 * t317;
t349 = -t314 * t377 + t381;
t271 = t294 * t315 + t349 * t318;
t337 = t314 * t379 + t380;
t327 = t337 * qJD(1) + t303 * qJD(2);
t244 = t271 * qJD(5) - t259 * t318 + t327 * t315;
t343 = t349 * t315;
t272 = t294 * t318 - t343;
t419 = t272 * qJD(5) + t259 * t315 + t327 * t318;
t263 = t271 ^ 2;
t385 = t313 * t316;
t299 = t310 * t314 + t311 * t385;
t378 = t318 * t319;
t289 = t299 * t315 + t313 * t378;
t286 = 0.1e1 / t289 ^ 2;
t252 = t263 * t286 + 0.1e1;
t248 = 0.1e1 / t252;
t298 = -t310 * t385 + t311 * t314;
t354 = qJD(2) * t313 * t319;
t284 = t298 * qJD(3) + t311 * t354;
t382 = t315 * t319;
t290 = t299 * t318 - t313 * t382;
t355 = t313 * t375;
t260 = t290 * qJD(5) + t284 * t315 - t318 * t355;
t285 = 0.1e1 / t289;
t387 = t271 * t286;
t341 = -t260 * t387 + t285 * t419;
t228 = t341 * t248;
t253 = atan2(t271, t289);
t246 = sin(t253);
t247 = cos(t253);
t342 = -t246 * t289 + t247 * t271;
t223 = t342 * t228 + t246 * t419 + t247 * t260;
t240 = t246 * t271 + t247 * t289;
t238 = 0.1e1 / t240 ^ 2;
t418 = t223 * t238;
t281 = -t303 * qJD(1) - t337 * qJD(2);
t304 = -t314 * t381 + t377;
t384 = t313 * t317;
t295 = -t304 * t310 + t311 * t384;
t357 = qJD(1) * t383;
t256 = t295 * qJD(3) + t281 * t311 + t310 * t357;
t296 = t304 * t311 + t310 * t384;
t274 = t296 * t318 + t337 * t315;
t280 = t314 * t353 + t356 + (-qJD(1) * t314 - qJD(2)) * t377;
t241 = t274 * qJD(5) + t256 * t315 + t280 * t318;
t273 = t296 * t315 - t337 * t318;
t403 = 0.2e1 * t273;
t237 = 0.1e1 / t240;
t413 = t237 * t418;
t351 = t403 * t413;
t417 = -t238 * t241 + t351;
t266 = 0.1e1 / t274 ^ 2;
t288 = t295 ^ 2;
t390 = t266 * t288;
t254 = 0.1e1 + t390;
t255 = -t296 * qJD(3) - t281 * t310 + t311 * t357;
t242 = -t273 * qJD(5) + t256 * t318 - t280 * t315;
t265 = 0.1e1 / t274;
t396 = t242 * t265 * t266;
t362 = t288 * t396;
t389 = t266 * t295;
t401 = (t255 * t389 - t362) / t254 ^ 2;
t414 = 0.2e1 * t401;
t412 = -0.2e1 * t273;
t411 = t260 * t286;
t338 = t285 * t293 - t298 * t387;
t410 = t315 * t338;
t250 = 0.1e1 / t254;
t392 = t250 * t266;
t408 = t242 * t392 + t265 * t414;
t264 = t273 ^ 2;
t236 = t238 * t264 + 0.1e1;
t234 = 0.1e1 / t236;
t397 = t238 * t273;
t402 = (t241 * t397 - t264 * t413) / t236 ^ 2;
t407 = -t234 * t418 - 0.2e1 * t237 * t402;
t348 = t389 * t401;
t361 = t295 * t396;
t406 = 0.2e1 * t250 * t361 - t255 * t392 + 0.2e1 * t348;
t371 = 0.2e1 * t402;
t405 = t417 * t234 + t371 * t397;
t257 = t294 * qJD(3) - t282 * t310 + t311 * t358;
t404 = 0.2e1 * t271;
t391 = t285 * t411;
t400 = (-t263 * t391 + t387 * t419) / t252 ^ 2;
t399 = t234 * t237;
t395 = t246 * t273;
t394 = t247 * t273;
t393 = t250 * t265;
t388 = t271 * t285;
t386 = t295 * t315;
t374 = qJD(3) * t310;
t373 = qJD(5) * t311;
t372 = qJD(5) * t318;
t370 = -0.2e1 * t400;
t366 = t285 * t400;
t364 = t234 * t397;
t359 = t250 * t389;
t350 = t391 * t404;
t340 = t272 * t285 - t290 * t387;
t275 = -t303 * t318 - t311 * t343;
t297 = (t311 * t382 - t316 * t318) * t313;
t339 = -t275 * t285 - t297 * t387;
t333 = t311 * t337;
t332 = -t246 + (-t247 * t388 + t246) * t248;
t331 = qJD(3) * t337;
t330 = -qJD(5) * t333 - t281;
t329 = t304 * qJD(5) + t280 * t311 + t310 * t331;
t283 = -t299 * qJD(3) - t310 * t354;
t262 = ((-qJD(2) + t373) * t378 + (-t319 * t374 + (-qJD(2) * t311 + qJD(5)) * t316) * t315) * t313;
t261 = -t289 * qJD(5) + t284 * t318 + t315 * t355;
t245 = (-t349 * t373 - t282) * t318 + (t303 * qJD(5) - t311 * t327 + t349 * t374) * t315;
t233 = t248 * t410;
t232 = t339 * t248;
t231 = t340 * t248;
t226 = (t246 * t293 + t247 * t298) * t315 + t342 * t233;
t224 = t342 * t231 + t246 * t272 + t247 * t290;
t222 = t339 * t370 + (t297 * t350 - t245 * t285 + (t260 * t275 - t262 * t271 - t297 * t419) * t286) * t248;
t220 = t340 * t370 + (t290 * t350 - t244 * t285 + (-t260 * t272 - t261 * t271 - t290 * t419) * t286) * t248;
t219 = t370 * t410 + (t338 * t372 + (t298 * t350 - t257 * t285 + (-t260 * t293 - t271 * t283 - t298 * t419) * t286) * t315) * t248;
t1 = [t366 * t403 + (-t241 * t285 + t273 * t411) * t248, t222, t219, 0, t220, 0; t419 * t399 - (t332 * t241 + ((t228 * t248 * t388 + t370) * t246 + (t366 * t404 - t228 + (t228 - t341) * t248) * t247) * t273) * t364 + t407 * t271 + t405 * t332 * t273 (t329 * t315 + t330 * t318) * t399 - ((t222 * t271 + t232 * t419 + t262 + (-t232 * t289 - t275) * t228) * t247 + (-t222 * t289 - t232 * t260 - t245 + (-t232 * t271 - t297) * t228) * t246) * t364 + t407 * (-t304 * t318 - t315 * t333) + t405 * (t342 * t232 - t246 * t275 + t247 * t297) (t226 * t397 - t237 * t386) * t371 + ((t255 * t315 + t295 * t372) * t237 + t417 * t226 + (-t386 * t223 - (t298 * t372 + t219 * t271 + t233 * t419 + t283 * t315 + (-t233 * t289 + t293 * t315) * t228) * t394 - (t293 * t372 - t219 * t289 - t233 * t260 - t257 * t315 + (-t233 * t271 - t298 * t315) * t228) * t395) * t238) * t234, 0 (t224 * t397 - t237 * t274) * t371 + (t224 * t351 + t242 * t237 + (-t274 * t223 - t224 * t241 - (t220 * t271 + t231 * t419 + t261 + (-t231 * t289 + t272) * t228) * t394 - (-t220 * t289 - t231 * t260 - t244 + (-t231 * t271 - t290) * t228) * t395) * t238) * t234, 0; t244 * t359 - t257 * t393 + t406 * t272 - t408 * t293 -(-t330 * t315 + t329 * t318) * t359 + (-t280 * t310 + t311 * t331) * t393 - t408 * t310 * t337 + t406 * (t304 * t315 - t318 * t333) (t265 * t296 + t318 * t390) * t414 + (0.2e1 * t318 * t362 - t256 * t265 + (qJD(5) * t288 * t315 - 0.2e1 * t255 * t295 * t318 + t242 * t296) * t266) * t250, 0, t348 * t412 + (t361 * t412 + (t241 * t295 + t255 * t273) * t266) * t250, 0;];
JaD_rot  = t1;
