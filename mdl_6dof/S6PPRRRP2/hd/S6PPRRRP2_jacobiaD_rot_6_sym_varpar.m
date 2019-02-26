% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRP2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRRP2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:05
% EndTime: 2019-02-26 19:42:09
% DurationCPUTime: 3.17s
% Computational Cost: add. (21168->168), mult. (62476->326), div. (784->12), fcn. (81881->17), ass. (0->139)
t416 = sin(pkin(12));
t417 = sin(pkin(11));
t383 = t417 * t416;
t420 = cos(pkin(12));
t421 = cos(pkin(11));
t390 = t421 * t420;
t423 = cos(pkin(6));
t373 = -t423 * t390 + t383;
t418 = sin(pkin(7));
t419 = sin(pkin(6));
t388 = t419 * t418;
t422 = cos(pkin(7));
t433 = -t373 * t422 - t421 * t388;
t354 = sin(qJ(3));
t385 = t417 * t420;
t389 = t421 * t416;
t374 = t423 * t389 + t385;
t424 = cos(qJ(3));
t334 = t354 * t433 + t374 * t424;
t391 = t422 * t419;
t343 = t373 * t418 - t421 * t391;
t353 = sin(qJ(4));
t356 = cos(qJ(4));
t320 = -t334 * t353 + t343 * t356;
t428 = -t374 * t354 + t433 * t424;
t328 = t428 * qJD(3);
t298 = qJD(4) * t320 + t328 * t356;
t321 = t334 * t356 + t343 * t353;
t355 = cos(qJ(5));
t352 = sin(qJ(5));
t430 = t352 * t428;
t308 = t321 * t355 - t430;
t363 = qJD(3) * t334;
t280 = qJD(5) * t308 + t298 * t352 - t355 * t363;
t364 = t428 * t355;
t306 = t321 * t352 + t364;
t301 = t306 ^ 2;
t387 = t419 * t416;
t426 = t420 * t391 + t418 * t423;
t342 = t354 * t426 + t424 * t387;
t347 = -t420 * t388 + t423 * t422;
t338 = t342 * t356 + t347 * t353;
t369 = -t354 * t387 + t424 * t426;
t324 = t338 * t352 + t355 * t369;
t318 = 0.1e1 / t324 ^ 2;
t292 = t301 * t318 + 0.1e1;
t289 = 0.1e1 / t292;
t337 = -t342 * t353 + t347 * t356;
t339 = t369 * qJD(3);
t315 = qJD(4) * t337 + t339 * t356;
t325 = t338 * t355 - t352 * t369;
t340 = t342 * qJD(3);
t294 = qJD(5) * t325 + t315 * t352 - t340 * t355;
t317 = 0.1e1 / t324;
t406 = t306 * t318;
t267 = (-t280 * t317 + t294 * t406) * t289;
t293 = atan2(-t306, t324);
t285 = sin(t293);
t286 = cos(t293);
t382 = -t285 * t324 - t286 * t306;
t263 = t382 * t267 - t280 * t285 + t286 * t294;
t279 = -t285 * t306 + t286 * t324;
t276 = 0.1e1 / t279;
t277 = 0.1e1 / t279 ^ 2;
t432 = t263 * t276 * t277;
t348 = -t423 * t383 + t390;
t375 = t423 * t385 + t389;
t384 = t417 * t419;
t427 = t375 * t422 - t418 * t384;
t336 = t348 * t424 - t354 * t427;
t344 = t375 * t418 + t422 * t384;
t323 = t336 * t356 + t344 * t353;
t335 = t348 * t354 + t424 * t427;
t309 = t323 * t352 - t335 * t355;
t431 = -0.2e1 * t309;
t392 = 0.2e1 * t309 * t432;
t379 = -t317 * t320 + t337 * t406;
t429 = t352 * t379;
t409 = t294 * t317 * t318;
t425 = -0.2e1 * (t280 * t406 - t301 * t409) / t292 ^ 2;
t310 = t323 * t355 + t335 * t352;
t303 = 0.1e1 / t310;
t304 = 0.1e1 / t310 ^ 2;
t322 = -t336 * t353 + t344 * t356;
t316 = t322 ^ 2;
t408 = t304 * t316;
t291 = 0.1e1 + t408;
t329 = t335 * qJD(3);
t299 = -qJD(4) * t323 + t329 * t353;
t300 = qJD(4) * t322 - t329 * t356;
t330 = t336 * qJD(3);
t283 = -t309 * qJD(5) + t300 * t355 + t330 * t352;
t412 = t283 * t303 * t304;
t395 = t316 * t412;
t407 = t304 * t322;
t415 = (t299 * t407 - t395) / t291 ^ 2;
t414 = t277 * t309;
t282 = t310 * qJD(5) + t300 * t352 - t330 * t355;
t413 = t282 * t277;
t411 = t285 * t309;
t410 = t286 * t309;
t405 = t322 * t352;
t404 = t352 * t356;
t403 = t355 * t356;
t402 = qJD(4) * t353;
t401 = qJD(5) * t352;
t400 = qJD(5) * t355;
t399 = qJD(5) * t356;
t302 = t309 ^ 2;
t275 = t277 * t302 + 0.1e1;
t398 = 0.2e1 * (-t302 * t432 + t309 * t413) / t275 ^ 2;
t397 = 0.2e1 * t415;
t394 = t322 * t412;
t393 = -0.2e1 * t306 * t409;
t381 = -t308 * t317 + t325 * t406;
t365 = t356 * t428;
t311 = -t334 * t355 + t352 * t365;
t326 = -t342 * t355 + t369 * t404;
t380 = -t311 * t317 + t326 * t406;
t314 = -qJD(4) * t338 - t339 * t353;
t313 = -t335 * t403 + t336 * t352;
t312 = -t335 * t404 - t336 * t355;
t297 = -qJD(4) * t321 - t328 * t353;
t296 = (t369 * t399 - t339) * t355 + (qJD(5) * t342 - t340 * t356 - t369 * t402) * t352;
t295 = -qJD(5) * t324 + t315 * t355 + t340 * t352;
t287 = 0.1e1 / t291;
t284 = -t328 * t355 + t334 * t401 - t363 * t404 + t365 * t400 - t402 * t430;
t281 = -qJD(5) * t364 + t298 * t355 - t321 * t401 + t352 * t363;
t273 = 0.1e1 / t275;
t272 = t289 * t429;
t271 = t380 * t289;
t269 = t381 * t289;
t266 = (-t285 * t320 + t286 * t337) * t352 + t382 * t272;
t265 = t382 * t271 - t285 * t311 + t286 * t326;
t264 = t382 * t269 - t285 * t308 + t286 * t325;
t262 = t380 * t425 + (t326 * t393 - t284 * t317 + (t280 * t326 + t294 * t311 + t296 * t306) * t318) * t289;
t260 = t381 * t425 + (t325 * t393 - t281 * t317 + (t280 * t325 + t294 * t308 + t295 * t306) * t318) * t289;
t259 = t425 * t429 + (t379 * t400 + (t337 * t393 - t297 * t317 + (t280 * t337 + t294 * t320 + t306 * t314) * t318) * t352) * t289;
t1 = [0, 0, t262, t259, t260, 0; 0, 0 (t265 * t414 - t276 * t312) * t398 + (t265 * t392 + (-t312 * t263 - t265 * t282 - (-t262 * t306 - t271 * t280 + t296 + (-t271 * t324 - t311) * t267) * t410 - (-t262 * t324 - t271 * t294 - t284 + (t271 * t306 - t326) * t267) * t411) * t277 + ((-t335 * t399 + t329) * t355 + (qJD(5) * t336 - t330 * t356 + t335 * t402) * t352) * t276) * t273 (t266 * t414 - t276 * t405) * t398 + ((t299 * t352 + t322 * t400) * t276 + (-t413 + t392) * t266 + (-t405 * t263 - (t337 * t400 - t259 * t306 - t272 * t280 + t314 * t352 + (-t272 * t324 - t320 * t352) * t267) * t410 - (-t320 * t400 - t259 * t324 - t272 * t294 - t297 * t352 + (t272 * t306 - t337 * t352) * t267) * t411) * t277) * t273 (t264 * t414 - t276 * t310) * t398 + (t264 * t392 + t283 * t276 + (-t310 * t263 - t264 * t282 - (-t260 * t306 - t269 * t280 + t295 + (-t269 * t324 - t308) * t267) * t410 - (-t260 * t324 - t269 * t294 - t281 + (t269 * t306 - t325) * t267) * t411) * t277) * t273, 0; 0, 0 (-t303 * t335 * t353 + t313 * t407) * t397 + (0.2e1 * t313 * t394 + (t335 * qJD(4) * t356 + t330 * t353) * t303 + (-(-t329 * t352 - t330 * t403 + t336 * t400) * t322 - t313 * t299 + (-t353 * t283 - (t352 * t399 + t355 * t402) * t322) * t335) * t304) * t287 (t303 * t323 + t355 * t408) * t397 + (0.2e1 * t355 * t395 - t300 * t303 + (-0.2e1 * t299 * t322 * t355 + t283 * t323 + t316 * t401) * t304) * t287, t407 * t415 * t431 + (t394 * t431 + (t282 * t322 + t299 * t309) * t304) * t287, 0;];
JaD_rot  = t1;
