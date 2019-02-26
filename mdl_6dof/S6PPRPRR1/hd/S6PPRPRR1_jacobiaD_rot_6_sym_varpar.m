% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRPRR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:49
% EndTime: 2019-02-26 19:39:51
% DurationCPUTime: 2.10s
% Computational Cost: add. (14279->140), mult. (41349->282), div. (559->12), fcn. (54415->19), ass. (0->131)
t357 = sin(pkin(11));
t360 = cos(pkin(12));
t361 = cos(pkin(6));
t356 = sin(pkin(12));
t421 = cos(pkin(11));
t393 = t421 * t356;
t345 = t357 * t360 + t361 * t393;
t355 = sin(pkin(13));
t359 = cos(pkin(13));
t364 = sin(qJ(3));
t367 = cos(qJ(3));
t350 = t364 * t355 - t367 * t359;
t420 = sin(pkin(7));
t394 = t367 * t420;
t396 = t364 * t420;
t375 = -t355 * t394 - t359 * t396;
t422 = cos(pkin(7));
t395 = t367 * t422;
t397 = t364 * t422;
t376 = -t355 * t395 - t359 * t397;
t392 = t421 * t360;
t381 = t357 * t356 - t361 * t392;
t358 = sin(pkin(6));
t398 = t358 * t421;
t321 = -t345 * t350 + t375 * t398 + t376 * t381;
t363 = sin(qJ(5));
t366 = cos(qJ(5));
t374 = t381 * t420 - t422 * t398;
t311 = t321 * t366 + t374 * t363;
t340 = -t355 * t396 + t359 * t394;
t336 = t340 * qJD(3);
t342 = -t355 * t397 + t359 * t395;
t338 = t342 * qJD(3);
t351 = -t367 * t355 - t364 * t359;
t349 = t351 * qJD(3);
t315 = -t336 * t398 - t381 * t338 + t345 * t349;
t287 = t311 * qJD(5) + t315 * t363;
t309 = t321 * t363 - t374 * t366;
t307 = t309 ^ 2;
t333 = -t361 * t375 + (-t350 * t356 - t360 * t376) * t358;
t344 = -t358 * t360 * t420 + t361 * t422;
t328 = t333 * t363 - t344 * t366;
t326 = 0.1e1 / t328 ^ 2;
t301 = t307 * t326 + 0.1e1;
t299 = 0.1e1 / t301;
t329 = t333 * t366 + t344 * t363;
t331 = t361 * t336 + (t338 * t360 + t349 * t356) * t358;
t305 = t329 * qJD(5) + t331 * t363;
t325 = 0.1e1 / t328;
t409 = t309 * t326;
t271 = (-t287 * t325 + t305 * t409) * t299;
t302 = atan2(-t309, t328);
t297 = sin(t302);
t298 = cos(t302);
t385 = -t297 * t328 - t298 * t309;
t267 = t385 * t271 - t297 * t287 + t298 * t305;
t281 = -t297 * t309 + t298 * t328;
t278 = 0.1e1 / t281;
t279 = 0.1e1 / t281 ^ 2;
t426 = t267 * t278 * t279;
t405 = t357 * t361;
t346 = -t360 * t405 - t393;
t406 = t357 * t358;
t378 = -t346 * t420 + t422 * t406;
t347 = -t356 * t405 + t392;
t379 = -t346 * t376 - t347 * t350 - t375 * t406;
t312 = t363 * t379 - t378 * t366;
t425 = 0.2e1 * t312 * t426;
t320 = -t340 * t398 - t381 * t342 + t345 * t351;
t332 = t361 * t340 + (t342 * t360 + t351 * t356) * t358;
t382 = -t320 * t325 + t332 * t409;
t424 = t363 * t382;
t410 = t305 * t325 * t326;
t423 = -0.2e1 * (t287 * t409 - t307 * t410) / t301 ^ 2;
t313 = t378 * t363 + t366 * t379;
t323 = t340 * t406 + t346 * t342 + t347 * t351;
t362 = sin(qJ(6));
t365 = cos(qJ(6));
t296 = t313 * t365 - t323 * t362;
t292 = 0.1e1 / t296;
t293 = 0.1e1 / t296 ^ 2;
t380 = t336 * t406 + t346 * t338 + t347 * t349;
t290 = -t312 * qJD(5) + t366 * t380;
t337 = t375 * qJD(3);
t339 = t376 * qJD(3);
t348 = t350 * qJD(3);
t316 = t337 * t406 + t346 * t339 + t347 * t348;
t282 = t296 * qJD(6) + t290 * t362 + t316 * t365;
t295 = t313 * t362 + t323 * t365;
t291 = t295 ^ 2;
t286 = t291 * t293 + 0.1e1;
t414 = t293 * t295;
t403 = qJD(6) * t295;
t283 = t290 * t365 - t316 * t362 - t403;
t417 = t283 * t292 * t293;
t419 = (t282 * t414 - t291 * t417) / t286 ^ 2;
t418 = t279 * t312;
t289 = t313 * qJD(5) + t363 * t380;
t416 = t289 * t279;
t415 = t292 * t362;
t413 = t295 * t365;
t412 = t297 * t312;
t411 = t298 * t312;
t408 = t323 * t363;
t407 = t323 * t366;
t404 = qJD(5) * t366;
t308 = t312 ^ 2;
t277 = t308 * t279 + 0.1e1;
t402 = 0.2e1 * (-t308 * t426 + t312 * t416) / t277 ^ 2;
t401 = -0.2e1 * t419;
t399 = t295 * t417;
t391 = -0.2e1 * t309 * t410;
t386 = qJD(6) * t407 - t380;
t384 = t293 * t413 - t415;
t383 = -t311 * t325 + t329 * t409;
t377 = -qJD(5) * t408 + qJD(6) * t379 + t316 * t366;
t330 = t361 * t337 + (t339 * t360 + t348 * t356) * t358;
t314 = -t337 * t398 - t381 * t339 + t345 * t348;
t306 = -t328 * qJD(5) + t331 * t366;
t304 = t362 * t379 + t365 * t407;
t303 = t362 * t407 - t365 * t379;
t288 = -t309 * qJD(5) + t315 * t366;
t284 = 0.1e1 / t286;
t275 = 0.1e1 / t277;
t273 = t299 * t424;
t272 = t383 * t299;
t269 = (-t297 * t320 + t298 * t332) * t363 + t385 * t273;
t268 = t385 * t272 - t297 * t311 + t298 * t329;
t265 = t383 * t423 + (t329 * t391 - t288 * t325 + (t287 * t329 + t305 * t311 + t306 * t309) * t326) * t299;
t264 = t423 * t424 + (t382 * t404 + (t332 * t391 - t314 * t325 + (t287 * t332 + t305 * t320 + t309 * t330) * t326) * t363) * t299;
t1 = [0, 0, t264, 0, t265, 0; 0, 0 (t269 * t418 - t278 * t408) * t402 + ((t316 * t363 + t323 * t404) * t278 + (-t416 + t425) * t269 + (-t408 * t267 - (t332 * t404 - t264 * t309 - t273 * t287 + t330 * t363 + (-t273 * t328 - t320 * t363) * t271) * t411 - (-t320 * t404 - t264 * t328 - t273 * t305 - t314 * t363 + (t273 * t309 - t332 * t363) * t271) * t412) * t279) * t275, 0 (t268 * t418 - t278 * t313) * t402 + (t268 * t425 + t290 * t278 + (-t313 * t267 - t268 * t289 - (-t265 * t309 - t272 * t287 + t306 + (-t272 * t328 - t311) * t271) * t411 - (-t265 * t328 - t272 * t305 - t288 + (t272 * t309 - t329) * t271) * t412) * t279) * t275, 0; 0, 0, 0.2e1 * (-t292 * t303 + t304 * t414) * t419 + (0.2e1 * t304 * t399 + t386 * t292 * t365 + t377 * t415 + (t386 * t295 * t362 - t304 * t282 - t303 * t283 - t377 * t413) * t293) * t284, 0, t384 * t312 * t401 + (t384 * t289 + ((-qJD(6) * t292 - 0.2e1 * t399) * t365 + (t282 * t365 + (t283 - t403) * t362) * t293) * t312) * t284, t401 + 0.2e1 * (t282 * t293 * t284 + (-t284 * t417 - t293 * t419) * t295) * t295;];
JaD_rot  = t1;
