% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR14_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:55:46
% EndTime: 2019-02-26 22:55:48
% DurationCPUTime: 1.90s
% Computational Cost: add. (7994->144), mult. (24079->300), div. (442->12), fcn. (30676->17), ass. (0->138)
t349 = sin(pkin(8));
t353 = cos(pkin(8));
t357 = cos(qJ(2));
t358 = cos(qJ(1));
t430 = cos(pkin(6));
t431 = sin(qJ(2));
t386 = t430 * t431;
t432 = sin(qJ(1));
t342 = t432 * t357 + t358 * t386;
t391 = t357 * t430;
t372 = t358 * t431 + t432 * t391;
t336 = t372 * qJD(1) + t342 * qJD(2);
t350 = sin(pkin(7));
t354 = cos(pkin(7));
t351 = sin(pkin(6));
t393 = t351 * t432;
t388 = qJD(1) * t393;
t374 = t336 * t350 + t354 * t388;
t347 = t432 * t431;
t379 = t432 * t386;
t401 = t358 * t357;
t337 = -qJD(1) * t379 - qJD(2) * t347 + (qJD(2) * t430 + qJD(1)) * t401;
t348 = sin(pkin(14));
t352 = cos(pkin(14));
t373 = -t336 * t354 + t350 * t388;
t433 = -t337 * t348 + t373 * t352;
t285 = t433 * t349 - t374 * t353;
t383 = -t358 * t391 + t347;
t405 = t351 * t358;
t370 = t350 * t405 + t383 * t354;
t324 = t342 * t348 + t370 * t352;
t434 = t383 * t350 - t354 * t405;
t437 = -t324 * t349 - t434 * t353;
t402 = t354 * t352;
t378 = -t431 * t348 + t357 * t402;
t406 = t350 * t357;
t321 = -(t430 * t350 * t352 + t378 * t351) * t349 + (-t351 * t406 + t430 * t354) * t353;
t319 = 0.1e1 / t321 ^ 2;
t407 = t350 * t353;
t333 = (-(-t348 * t357 - t431 * t402) * t349 + t431 * t407) * t351;
t330 = qJD(2) * t333;
t417 = t319 * t330;
t299 = atan2(t437, t321);
t294 = sin(t299);
t295 = cos(t299);
t278 = t294 * t437 + t295 * t321;
t275 = 0.1e1 / t278;
t367 = t350 * t393 - t372 * t354;
t371 = t379 - t401;
t327 = t367 * t348 - t352 * t371;
t355 = sin(qJ(4));
t356 = cos(qJ(4));
t326 = t348 * t371 + t367 * t352;
t368 = -t372 * t350 - t354 * t393;
t366 = t368 * t349;
t365 = t326 * t353 - t366;
t293 = t327 * t356 + t365 * t355;
t287 = 0.1e1 / t293;
t318 = 0.1e1 / t321;
t276 = 0.1e1 / t278 ^ 2;
t288 = 0.1e1 / t293 ^ 2;
t311 = t326 * t349 + t368 * t353;
t308 = t311 ^ 2;
t273 = t308 * t276 + 0.1e1;
t335 = t342 * qJD(1) + t372 * qJD(2);
t334 = t383 * qJD(1) + t371 * qJD(2);
t392 = qJD(1) * t405;
t376 = t334 * t354 + t350 * t392;
t303 = t335 * t348 + t376 * t352;
t377 = -t334 * t350 + t354 * t392;
t284 = t303 * t349 - t377 * t353;
t423 = t284 * t276;
t307 = t437 ^ 2;
t298 = t307 * t319 + 0.1e1;
t296 = 0.1e1 / t298;
t382 = t285 * t318 - t417 * t437;
t269 = t382 * t296;
t385 = -t294 * t321 + t295 * t437;
t265 = t385 * t269 + t294 * t285 + t295 * t330;
t428 = t265 * t275 * t276;
t429 = (-t308 * t428 + t311 * t423) / t273 ^ 2;
t315 = (-t342 * t402 + t383 * t348) * t349 - t342 * t407;
t418 = t437 * t333;
t381 = -t315 * t318 + t319 * t418;
t270 = t381 * t296;
t266 = -t385 * t270 + t294 * t315 + t295 * t333;
t427 = t266 * t311;
t304 = -t335 * t352 + t376 * t348;
t369 = t303 * t353 + t377 * t349;
t279 = t293 * qJD(4) + t304 * t355 - t369 * t356;
t403 = t353 * t356;
t414 = t327 * t355;
t292 = -t326 * t403 + t356 * t366 + t414;
t286 = t292 ^ 2;
t283 = t286 * t288 + 0.1e1;
t422 = t288 * t292;
t280 = t304 * t356 + t369 * t355 + (t365 * t356 - t414) * qJD(4);
t424 = t280 * t287 * t288;
t426 = (t279 * t422 - t286 * t424) / t283 ^ 2;
t416 = t318 * t417;
t425 = (t285 * t319 * t437 - t307 * t416) / t298 ^ 2;
t421 = t294 * t311;
t420 = t295 * t311;
t419 = t437 * t318;
t325 = -t342 * t352 + t370 * t348;
t415 = t325 * t355;
t410 = t348 * t354;
t332 = -t372 * t352 + t371 * t410;
t413 = t332 * t355;
t409 = t349 * t355;
t408 = t349 * t356;
t404 = t353 * t355;
t400 = -0.2e1 * t429;
t399 = -0.2e1 * t428;
t398 = 0.2e1 * t426;
t397 = 0.2e1 * t425;
t396 = t350 * t408;
t390 = -0.2e1 * t318 * t425;
t389 = 0.2e1 * t292 * t424;
t384 = t324 * t353 - t349 * t434;
t331 = t372 * t348 + t371 * t402;
t380 = -t349 * t350 * t371 + t331 * t353;
t375 = t294 + (t295 * t419 - t294) * t296;
t291 = t325 * t356 + t384 * t355;
t302 = t332 * t356 + t380 * t355;
t329 = (t378 * t349 + t353 * t406) * t351 * qJD(2);
t316 = -t331 * t349 - t371 * t407;
t314 = t334 * t352 + t335 * t410;
t313 = -t334 * t348 + t335 * t402;
t306 = -t337 * t352 - t373 * t348;
t301 = -t380 * t356 + t413;
t300 = (t336 * t348 - t337 * t402) * t349 - t337 * t407;
t290 = -t384 * t356 + t415;
t281 = 0.1e1 / t283;
t271 = 0.1e1 / t273;
t268 = t375 * t311;
t264 = t381 * t397 + (0.2e1 * t416 * t418 + t300 * t318 + (-t285 * t333 - t315 * t330 - t329 * t437) * t319) * t296;
t1 = [t311 * t390 + (t284 * t318 - t311 * t417) * t296, t264, 0, 0, 0, 0; t437 * t275 * t400 + (t285 * t275 + (-t265 * t437 + t268 * t284) * t276) * t271 + ((t268 * t399 + t375 * t423) * t271 + (t268 * t400 + ((-t269 * t296 * t419 + t397) * t421 + (t437 * t390 + t269 + (-t269 + t382) * t296) * t420) * t271) * t276) * t311, 0.2e1 * (-t275 * t316 - t276 * t427) * t429 + ((-t313 * t349 - t335 * t407) * t275 + t399 * t427 + (-t316 * t265 + t266 * t284 + (t264 * t437 - t270 * t285 + t329 + (t270 * t321 + t315) * t269) * t420 + (-t264 * t321 + t270 * t330 + t300 + (t270 * t437 - t333) * t269) * t421) * t276) * t271, 0, 0, 0, 0; (-t287 * t290 + t291 * t422) * t398 + ((t306 * t355 + t374 * t408 + t403 * t433) * t287 + t291 * t389 + (-t290 * t280 - (t306 * t356 - t374 * t409 - t404 * t433) * t292 - t291 * t279) * t288 + (t291 * t287 - (t324 * t403 - t408 * t434 - t415) * t422) * qJD(4)) * t281 (-t287 * t301 + t302 * t422) * t398 + ((-t313 * t403 + t314 * t355 + t335 * t396) * t287 + t302 * t389 + (-t301 * t280 - (-t335 * t350 * t409 + t313 * t404 + t314 * t356) * t292 - t302 * t279) * t288 + (t302 * t287 - (t331 * t403 - t371 * t396 - t413) * t422) * qJD(4)) * t281, 0, -0.2e1 * t426 + 0.2e1 * (t279 * t288 * t281 + (-t281 * t424 - t288 * t426) * t292) * t292, 0, 0;];
JaD_rot  = t1;
