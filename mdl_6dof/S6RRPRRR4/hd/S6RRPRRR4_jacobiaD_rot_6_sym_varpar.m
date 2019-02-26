% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:47
% EndTime: 2019-02-26 21:55:49
% DurationCPUTime: 2.22s
% Computational Cost: add. (17152->158), mult. (32251->311), div. (989->12), fcn. (41704->15), ass. (0->139)
t359 = qJ(4) + qJ(5);
t356 = sin(t359);
t362 = cos(pkin(6));
t360 = sin(pkin(12));
t361 = cos(pkin(12));
t364 = sin(qJ(2));
t366 = cos(qJ(2));
t388 = t366 * t360 + t364 * t361;
t344 = t388 * t362;
t348 = t364 * t360 - t366 * t361;
t367 = cos(qJ(1));
t436 = sin(qJ(1));
t381 = -t367 * t344 + t436 * t348;
t357 = cos(t359);
t435 = sin(pkin(6));
t398 = t367 * t435;
t393 = t357 * t398;
t316 = -t356 * t381 + t393;
t399 = t366 * t435;
t400 = t364 * t435;
t375 = -t360 * t399 - t361 * t400;
t332 = -t356 * t375 - t362 * t357;
t304 = atan2(-t316, t332);
t299 = sin(t304);
t300 = cos(t304);
t282 = -t299 * t316 + t300 * t332;
t280 = 0.1e1 / t282 ^ 2;
t401 = t436 * t344;
t382 = -t367 * t348 - t401;
t392 = t436 * t435;
t321 = t356 * t382 - t357 * t392;
t315 = t321 ^ 2;
t278 = t315 * t280 + 0.1e1;
t358 = qJD(4) + qJD(5);
t343 = t348 * t362;
t338 = qJD(2) * t343;
t345 = t388 * qJD(2);
t374 = t381 * qJD(1) + t436 * t338 - t367 * t345;
t379 = t358 * t392 + t374;
t391 = qJD(1) * t398;
t417 = t357 * t358;
t286 = t356 * t379 - t357 * t391 + t382 * t417;
t428 = t286 * t280;
t314 = t316 ^ 2;
t330 = 0.1e1 / t332 ^ 2;
t303 = t314 * t330 + 0.1e1;
t301 = 0.1e1 / t303;
t311 = -t436 * t345 - qJD(1) * t401 + (-qJD(1) * t348 - t338) * t367;
t353 = t356 * t398;
t387 = qJD(1) * t392;
t288 = t311 * t356 - t358 * t353 - t357 * t387 - t381 * t417;
t341 = -t360 * t400 + t361 * t399;
t394 = t341 * qJD(2) + t358 * t362;
t312 = t394 * t356 - t375 * t417;
t329 = 0.1e1 / t332;
t422 = t316 * t330;
t386 = -t288 * t329 + t312 * t422;
t270 = t386 * t301;
t389 = -t299 * t332 - t300 * t316;
t265 = t270 * t389 - t299 * t288 + t300 * t312;
t279 = 0.1e1 / t282;
t281 = t279 * t280;
t433 = t265 * t281;
t410 = 0.2e1 * (-t315 * t433 + t321 * t428) / t278 ^ 2;
t440 = t312 * t330;
t324 = -t367 * t343 - t388 * t436;
t383 = -t324 * t329 + t341 * t422;
t439 = t356 * t383;
t289 = (t358 * t381 + t387) * t356 + t311 * t357 - t358 * t393;
t322 = t356 * t392 + t357 * t382;
t365 = cos(qJ(6));
t327 = t436 * t343 - t367 * t388;
t363 = sin(qJ(6));
t420 = t327 * t363;
t298 = t322 * t365 - t420;
t292 = 0.1e1 / t298;
t293 = 0.1e1 / t298 ^ 2;
t438 = -0.2e1 * t316;
t437 = 0.2e1 * t321;
t287 = t379 * t357 + (-t358 * t382 + t391) * t356;
t339 = t362 * t345;
t380 = t348 * qJD(2);
t308 = t324 * qJD(1) - t436 * t339 - t367 * t380;
t274 = qJD(6) * t298 + t287 * t363 - t308 * t365;
t419 = t327 * t365;
t297 = t322 * t363 + t419;
t291 = t297 ^ 2;
t285 = t291 * t293 + 0.1e1;
t427 = t293 * t297;
t411 = qJD(6) * t297;
t275 = t287 * t365 + t308 * t363 - t411;
t430 = t275 * t292 * t293;
t432 = (t274 * t427 - t291 * t430) / t285 ^ 2;
t424 = t329 * t440;
t431 = (t288 * t422 - t314 * t424) / t303 ^ 2;
t429 = t280 * t321;
t426 = t299 * t321;
t425 = t300 * t321;
t423 = t316 * t329;
t421 = t327 * t356;
t418 = t356 * t358;
t416 = t363 * t292;
t414 = t365 * t297;
t409 = -0.2e1 * t432;
t408 = 0.2e1 * t432;
t407 = -0.2e1 * t431;
t406 = t281 * t437;
t405 = t329 * t431;
t404 = t280 * t426;
t403 = t280 * t425;
t402 = t297 * t430;
t397 = 0.2e1 * t402;
t396 = t424 * t438;
t318 = -t357 * t381 - t353;
t390 = qJD(6) * t327 * t357 - t374;
t296 = -t318 * t365 + t324 * t363;
t295 = -t318 * t363 - t324 * t365;
t385 = t293 * t414 - t416;
t333 = t362 * t356 - t357 * t375;
t384 = -t318 * t329 + t333 * t422;
t378 = -t299 + (t300 * t423 + t299) * t301;
t376 = qJD(6) * t382 - t308 * t357 - t327 * t418;
t336 = t375 * qJD(2);
t313 = t394 * t357 + t375 * t418;
t310 = t327 * qJD(1) - t367 * t339 + t436 * t380;
t306 = t357 * t419 + t363 * t382;
t305 = t357 * t420 - t365 * t382;
t283 = 0.1e1 / t285;
t276 = 0.1e1 / t278;
t273 = t301 * t439;
t272 = t384 * t301;
t269 = t378 * t321;
t267 = (-t299 * t324 + t300 * t341) * t356 + t389 * t273;
t266 = t272 * t389 - t299 * t318 + t300 * t333;
t263 = t384 * t407 + (t333 * t396 - t289 * t329 + (t288 * t333 + t312 * t318 + t313 * t316) * t330) * t301;
t262 = t407 * t439 + (t383 * t417 + (t341 * t396 - t310 * t329 + (t288 * t341 + t312 * t324 + t316 * t336) * t330) * t356) * t301;
t261 = t385 * t321 * t409 + (t385 * t286 + ((-qJD(6) * t292 - 0.2e1 * t402) * t365 + (t274 * t365 + (t275 - t411) * t363) * t293) * t321) * t283;
t260 = (t266 * t429 - t279 * t322) * t410 + (t266 * t265 * t406 + t287 * t279 + (-t322 * t265 - t266 * t286 - (-t263 * t316 - t272 * t288 + t313 + (-t272 * t332 - t318) * t270) * t425 - (-t263 * t332 - t272 * t312 - t289 + (t272 * t316 - t333) * t270) * t426) * t280) * t276;
t1 = [t405 * t437 + (-t286 * t329 + t321 * t440) * t301, t262, 0, t263, t263, 0; t316 * t279 * t410 + (-t288 * t279 + (t265 * t316 - t269 * t286) * t280) * t276 + (t269 * t280 * t410 + (0.2e1 * t269 * t433 - (-t270 * t301 * t423 + t407) * t404 - (t405 * t438 - t270 + (t270 - t386) * t301) * t403 - t378 * t428) * t276) * t321 (t267 * t429 - t279 * t421) * t410 + (-t267 * t428 + (-t308 * t356 + t327 * t417) * t279 + (t267 * t406 - t280 * t421) * t265 - (t341 * t417 - t262 * t316 - t273 * t288 + t336 * t356 + (-t273 * t332 - t324 * t356) * t270) * t403 - (-t324 * t417 - t262 * t332 - t273 * t312 - t310 * t356 + (t273 * t316 - t341 * t356) * t270) * t404) * t276, 0, t260, t260, 0; (-t292 * t295 + t296 * t427) * t408 + ((qJD(6) * t296 - t289 * t363 - t310 * t365) * t292 + t296 * t397 + (-t295 * t275 - (-qJD(6) * t295 - t289 * t365 + t310 * t363) * t297 - t296 * t274) * t293) * t283 (-t292 * t305 + t306 * t427) * t408 + (t306 * t397 + t390 * t292 * t365 + t376 * t416 + (t390 * t297 * t363 - t306 * t274 - t305 * t275 - t376 * t414) * t293) * t283, 0, t261, t261, t409 + 0.2e1 * (t274 * t293 * t283 + (-t283 * t430 - t293 * t432) * t297) * t297;];
JaD_rot  = t1;
