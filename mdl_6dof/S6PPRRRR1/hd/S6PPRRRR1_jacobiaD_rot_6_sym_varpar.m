% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRRR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:43
% EndTime: 2019-02-26 19:42:46
% DurationCPUTime: 2.26s
% Computational Cost: add. (18272->124), mult. (40377->248), div. (822->12), fcn. (53065->17), ass. (0->126)
t415 = sin(pkin(13));
t416 = sin(pkin(12));
t381 = t416 * t415;
t419 = cos(pkin(13));
t420 = cos(pkin(12));
t387 = t420 * t419;
t422 = cos(pkin(6));
t348 = -t422 * t381 + t387;
t355 = sin(qJ(3));
t357 = cos(qJ(3));
t383 = t416 * t419;
t386 = t420 * t415;
t373 = t422 * t383 + t386;
t418 = sin(pkin(6));
t382 = t416 * t418;
t417 = sin(pkin(7));
t421 = cos(pkin(7));
t424 = t373 * t421 - t417 * t382;
t337 = t348 * t355 + t424 * t357;
t347 = t422 * t386 + t383;
t372 = -t422 * t387 + t381;
t385 = t418 * t417;
t367 = -t372 * t421 - t420 * t385;
t336 = t347 * t357 + t367 * t355;
t353 = qJ(4) + qJ(5);
t350 = sin(t353);
t335 = -t347 * t355 + t367 * t357;
t352 = qJD(4) + qJD(5);
t388 = t421 * t418;
t366 = t372 * t417 - t420 * t388;
t365 = t335 * qJD(3) + t366 * t352;
t351 = cos(t353);
t397 = t351 * t352;
t298 = t336 * t397 + t365 * t350;
t320 = t336 * t350 - t366 * t351;
t318 = t320 ^ 2;
t371 = t419 * t388 + t422 * t417;
t384 = t418 * t415;
t344 = t371 * t355 + t357 * t384;
t346 = -t419 * t385 + t422 * t421;
t333 = t344 * t350 - t346 * t351;
t330 = 0.1e1 / t333 ^ 2;
t312 = t318 * t330 + 0.1e1;
t310 = 0.1e1 / t312;
t343 = -t355 * t384 + t371 * t357;
t390 = t343 * qJD(3) + t346 * t352;
t316 = t344 * t397 + t390 * t350;
t329 = 0.1e1 / t333;
t404 = t320 * t330;
t282 = (-t298 * t329 + t316 * t404) * t310;
t313 = atan2(-t320, t333);
t308 = sin(t313);
t309 = cos(t313);
t380 = -t308 * t333 - t309 * t320;
t278 = t380 * t282 - t308 * t298 + t309 * t316;
t292 = -t308 * t320 + t309 * t333;
t289 = 0.1e1 / t292;
t290 = 0.1e1 / t292 ^ 2;
t427 = t278 * t289 * t290;
t338 = t348 * t357 - t355 * t424;
t368 = t373 * t417 + t421 * t382;
t323 = t338 * t350 - t368 * t351;
t426 = 0.2e1 * t323 * t427;
t376 = -t329 * t335 + t343 * t404;
t425 = t350 * t376;
t405 = t316 * t329 * t330;
t423 = -0.2e1 * (t298 * t404 - t318 * t405) / t312 ^ 2;
t324 = t338 * t351 + t368 * t350;
t356 = cos(qJ(6));
t354 = sin(qJ(6));
t402 = t337 * t354;
t307 = t324 * t356 + t402;
t303 = 0.1e1 / t307;
t304 = 0.1e1 / t307 ^ 2;
t327 = t337 * qJD(3);
t364 = t368 * t352 - t327;
t398 = t350 * t352;
t301 = -t338 * t398 + t364 * t351;
t328 = t338 * qJD(3);
t293 = t307 * qJD(6) + t301 * t354 - t328 * t356;
t401 = t337 * t356;
t306 = t324 * t354 - t401;
t302 = t306 ^ 2;
t297 = t302 * t304 + 0.1e1;
t409 = t304 * t306;
t396 = qJD(6) * t306;
t294 = t301 * t356 + t328 * t354 - t396;
t412 = t294 * t303 * t304;
t414 = (t293 * t409 - t302 * t412) / t297 ^ 2;
t413 = t290 * t323;
t300 = t338 * t397 + t364 * t350;
t411 = t300 * t290;
t410 = t303 * t354;
t408 = t306 * t356;
t407 = t308 * t323;
t406 = t309 * t323;
t403 = t337 * t350;
t319 = t323 ^ 2;
t288 = t319 * t290 + 0.1e1;
t395 = 0.2e1 * (-t319 * t427 + t323 * t411) / t288 ^ 2;
t394 = -0.2e1 * t414;
t392 = t306 * t412;
t391 = -0.2e1 * t320 * t405;
t389 = qJD(6) * t337 * t351 - t327;
t378 = t304 * t408 - t410;
t322 = t336 * t351 + t366 * t350;
t334 = t344 * t351 + t346 * t350;
t377 = -t322 * t329 + t334 * t404;
t374 = qJD(6) * t338 - t328 * t351 + t337 * t398;
t342 = t344 * qJD(3);
t326 = t336 * qJD(3);
t317 = -t344 * t398 + t390 * t351;
t315 = t338 * t354 - t351 * t401;
t314 = -t338 * t356 - t351 * t402;
t299 = -t336 * t398 + t365 * t351;
t295 = 0.1e1 / t297;
t286 = 0.1e1 / t288;
t284 = t310 * t425;
t283 = t377 * t310;
t280 = (-t308 * t335 + t309 * t343) * t350 + t380 * t284;
t279 = t380 * t283 - t308 * t322 + t309 * t334;
t276 = t377 * t423 + (t334 * t391 - t299 * t329 + (t298 * t334 + t316 * t322 + t317 * t320) * t330) * t310;
t275 = t423 * t425 + (t376 * t397 + (t343 * t391 + t326 * t329 + (t298 * t343 + t316 * t335 - t320 * t342) * t330) * t350) * t310;
t274 = t378 * t323 * t394 + (t378 * t300 + ((-qJD(6) * t303 - 0.2e1 * t392) * t356 + (t293 * t356 + (t294 - t396) * t354) * t304) * t323) * t295;
t273 = (t279 * t413 - t289 * t324) * t395 + (t279 * t426 + t301 * t289 + (-t324 * t278 - t279 * t300 - (-t276 * t320 - t283 * t298 + t317 + (-t283 * t333 - t322) * t282) * t406 - (-t276 * t333 - t283 * t316 - t299 + (t283 * t320 - t334) * t282) * t407) * t290) * t286;
t1 = [0, 0, t275, t276, t276, 0; 0, 0 (t280 * t413 + t289 * t403) * t395 + ((-t328 * t350 - t337 * t397) * t289 + (-t411 + t426) * t280 + (t403 * t278 - (t343 * t397 - t275 * t320 - t284 * t298 - t342 * t350 + (-t284 * t333 - t335 * t350) * t282) * t406 - (-t335 * t397 - t275 * t333 - t284 * t316 + t326 * t350 + (t284 * t320 - t343 * t350) * t282) * t407) * t290) * t286, t273, t273, 0; 0, 0, 0.2e1 * (-t303 * t314 + t315 * t409) * t414 + (0.2e1 * t315 * t392 - t389 * t303 * t356 + t374 * t410 + (-t389 * t306 * t354 - t315 * t293 - t314 * t294 - t374 * t408) * t304) * t295, t274, t274, t394 + 0.2e1 * (t293 * t304 * t295 + (-t295 * t412 - t304 * t414) * t306) * t306;];
JaD_rot  = t1;
