% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRRR6
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRR6_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobiaD_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:52
% EndTime: 2019-02-26 20:21:54
% DurationCPUTime: 1.76s
% Computational Cost: add. (8706->151), mult. (26595->311), div. (532->12), fcn. (34071->17), ass. (0->132)
t327 = sin(pkin(8));
t328 = sin(pkin(7));
t330 = cos(pkin(14));
t331 = cos(pkin(8));
t326 = sin(pkin(14));
t336 = sin(qJ(2));
t333 = cos(pkin(6));
t339 = cos(qJ(2));
t375 = t333 * t339;
t353 = -t326 * t336 + t330 * t375;
t329 = sin(pkin(6));
t332 = cos(pkin(7));
t381 = t329 * t332;
t376 = t333 * t336;
t321 = t326 * t339 + t330 * t376;
t335 = sin(qJ(3));
t338 = cos(qJ(3));
t385 = t328 * t329;
t347 = -t330 * t385 + t353 * t332;
t402 = -t321 * t335 + t347 * t338;
t294 = t402 * t327 - (-t353 * t328 - t330 * t381) * t331;
t290 = t294 ^ 2;
t371 = t338 * t339;
t374 = t335 * t336;
t350 = t332 * t371 - t374;
t382 = t328 * t339;
t383 = t328 * t333;
t305 = -(t350 * t329 + t338 * t383) * t327 + (-t329 * t382 + t333 * t332) * t331;
t303 = 0.1e1 / t305 ^ 2;
t278 = t290 * t303 + 0.1e1;
t306 = -t321 * t338 - t347 * t335;
t317 = t353 * qJD(2);
t318 = t321 * qJD(2);
t370 = qJD(3) * t327;
t377 = t332 * t338;
t384 = t328 * t331;
t280 = (-t317 * t335 - t318 * t377) * t327 - t318 * t384 + t306 * t370;
t372 = t336 * t338;
t373 = t335 * t339;
t346 = -(-t332 * t372 - t373) * t327 + t336 * t384;
t349 = -t332 * t373 - t372;
t364 = qJD(3) * t383;
t296 = t335 * t327 * t364 + (t346 * qJD(2) - t349 * t370) * t329;
t302 = 0.1e1 / t305;
t394 = t296 * t302 * t303;
t395 = t294 * t303;
t403 = (t280 * t395 - t290 * t394) / t278 ^ 2;
t351 = t326 * t376 - t330 * t339;
t352 = t326 * t375 + t330 * t336;
t354 = t326 * t385 - t332 * t352;
t308 = t354 * t335 - t338 * t351;
t279 = atan2(t294, t305);
t274 = sin(t279);
t275 = cos(t279);
t262 = t274 * t294 + t275 * t305;
t259 = 0.1e1 / t262;
t334 = sin(qJ(4));
t337 = cos(qJ(4));
t390 = t351 * t335;
t307 = t354 * t338 + t390;
t315 = t326 * t381 + t328 * t352;
t360 = t307 * t331 + t315 * t327;
t273 = t308 * t337 + t360 * t334;
t269 = 0.1e1 / t273;
t260 = 0.1e1 / t262 ^ 2;
t270 = 0.1e1 / t273 ^ 2;
t276 = 0.1e1 / t278;
t252 = (t280 * t302 - t296 * t395) * t276;
t361 = -t274 * t305 + t275 * t294;
t248 = t361 * t252 + t274 * t280 + t275 * t296;
t401 = t248 * t259 * t260;
t319 = t352 * qJD(2);
t320 = t351 * qJD(2);
t378 = t332 * t335;
t289 = t320 * t378 + qJD(3) * t390 + (t354 * qJD(3) - t319) * t338;
t288 = -t308 * qJD(3) + t319 * t335 + t320 * t377;
t387 = t327 * t328;
t359 = -t288 * t331 + t320 * t387;
t263 = t273 * qJD(4) + t289 * t334 + t359 * t337;
t379 = t331 * t337;
t386 = t327 * t337;
t393 = t308 * t334;
t272 = -t307 * t379 - t315 * t386 + t393;
t268 = t272 ^ 2;
t267 = t268 * t270 + 0.1e1;
t397 = t270 * t272;
t264 = t289 * t337 - t359 * t334 + (t360 * t337 - t393) * qJD(4);
t398 = t264 * t269 * t270;
t400 = (t263 * t397 - t268 * t398) / t267 ^ 2;
t295 = -t307 * t327 + t315 * t331;
t399 = t260 * t295;
t281 = -t288 * t327 - t320 * t384;
t396 = t281 * t260;
t355 = t338 * t352 - t351 * t378;
t392 = t355 * t334;
t380 = t331 * t334;
t291 = t295 ^ 2;
t257 = t260 * t291 + 0.1e1;
t369 = 0.2e1 * (-t291 * t401 + t295 * t396) / t257 ^ 2;
t368 = 0.2e1 * t400;
t367 = t294 * t394;
t366 = t328 * t386;
t363 = 0.2e1 * t272 * t398;
t362 = 0.2e1 * t295 * t401;
t297 = (-t321 * t377 - t353 * t335) * t327 - t321 * t384;
t311 = t346 * t329;
t358 = -t297 * t302 + t311 * t395;
t314 = t349 * t329 - t335 * t383;
t357 = t302 * t306 + t314 * t395;
t285 = t307 * t334 + t308 * t379;
t286 = t307 * t337 - t308 * t380;
t309 = t335 * t352 + t351 * t377;
t356 = t309 * t331 - t351 * t387;
t348 = t332 * t374 - t371;
t284 = t356 * t334 - t337 * t355;
t301 = -t338 * t364 + (t348 * qJD(2) - t350 * qJD(3)) * t329;
t299 = (-t348 * t370 + (t350 * t327 + t331 * t382) * qJD(2)) * t329;
t298 = -t309 * t327 - t351 * t384;
t293 = t309 * qJD(3) + t319 * t378 + t320 * t338;
t292 = t355 * qJD(3) + t319 * t377 - t320 * t335;
t287 = -t402 * qJD(3) - t317 * t338 + t318 * t378;
t283 = -t356 * t337 - t392;
t282 = (-t317 * t377 + t318 * t335 + (t321 * t378 - t353 * t338) * qJD(3)) * t327 - t317 * t384;
t265 = 0.1e1 / t267;
t255 = 0.1e1 / t257;
t254 = t357 * t327 * t276;
t253 = t358 * t276;
t250 = (t274 * t306 - t275 * t314) * t327 + t361 * t254;
t249 = -t361 * t253 + t274 * t297 + t275 * t311;
t247 = (-0.2e1 * t357 * t403 + (-0.2e1 * t314 * t367 + t287 * t302 + (t280 * t314 + t294 * t301 - t296 * t306) * t303) * t276) * t327;
t246 = 0.2e1 * t358 * t403 + (0.2e1 * t311 * t367 + t282 * t302 + (-t280 * t311 - t294 * t299 - t296 * t297) * t303) * t276;
t1 = [0, t246, t247, 0, 0, 0; 0 (t249 * t399 - t259 * t298) * t369 + ((-t292 * t327 - t319 * t384) * t259 + t249 * t362 + (-t298 * t248 - t249 * t281 + (-(t246 * t294 - t253 * t280 + t299 + (t253 * t305 + t297) * t252) * t275 - (-t246 * t305 + t253 * t296 + t282 + (t253 * t294 - t311) * t252) * t274) * t295) * t260) * t255 (-t259 * t308 * t327 + t250 * t399) * t369 + (-(t361 * t247 + (-t262 * t252 - t274 * t296 + t275 * t280) * t254) * t399 + (t362 - t396) * t250 + (t289 * t259 + (-t308 * t248 - (t274 * t287 - t275 * t301 + (t274 * t314 + t275 * t306) * t252) * t295) * t260) * t327) * t255, 0, 0, 0; 0 (-t269 * t283 + t284 * t397) * t368 + ((-t292 * t379 + t293 * t334 + t319 * t366) * t269 + t284 * t363 + (-t283 * t264 - (-t319 * t334 * t387 + t292 * t380 + t293 * t337) * t272 - t284 * t263) * t270 + (t284 * t269 - (t309 * t379 - t351 * t366 + t392) * t397) * qJD(4)) * t265 (-t269 * t285 + t286 * t397) * t368 + ((t286 * qJD(4) + t288 * t334 + t289 * t379) * t269 + t286 * t363 + (-t285 * t264 - (-t285 * qJD(4) + t288 * t337 - t289 * t380) * t272 - t286 * t263) * t270) * t265, -0.2e1 * t400 + 0.2e1 * (t263 * t265 * t270 + (-t265 * t398 - t270 * t400) * t272) * t272, 0, 0;];
JaD_rot  = t1;
