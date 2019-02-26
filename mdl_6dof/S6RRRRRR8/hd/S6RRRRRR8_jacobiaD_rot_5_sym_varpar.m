% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR8_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:51:33
% EndTime: 2019-02-26 22:51:35
% DurationCPUTime: 1.97s
% Computational Cost: add. (8994->176), mult. (26325->322), div. (723->12), fcn. (32770->15), ass. (0->162)
t366 = cos(pkin(7));
t367 = sin(qJ(3));
t466 = cos(pkin(6));
t467 = sin(qJ(2));
t416 = t466 * t467;
t468 = sin(qJ(1));
t400 = t468 * t416;
t368 = cos(qJ(2));
t369 = cos(qJ(1));
t447 = t369 * t368;
t384 = t400 - t447;
t429 = t368 * t466;
t385 = t369 * t467 + t468 * t429;
t365 = sin(pkin(6));
t465 = sin(pkin(7));
t430 = t365 * t465;
t410 = t468 * t430;
t469 = cos(qJ(3));
t325 = -t384 * t469 + (-t385 * t366 + t410) * t367;
t450 = t365 * t366;
t420 = t468 * t450;
t344 = t385 * t465 + t420;
t364 = qJ(4) + qJ(5);
t361 = sin(t364);
t362 = cos(t364);
t305 = t325 * t361 - t344 * t362;
t477 = 0.2e1 * t305;
t419 = t469 * t467;
t448 = t367 * t368;
t392 = t366 * t419 + t448;
t394 = t366 * t448 + t419;
t412 = t466 * t465;
t409 = t367 * t412;
t315 = qJD(3) * t409 + (t392 * qJD(2) + t394 * qJD(3)) * t365;
t432 = t467 * t367;
t433 = t469 * t368;
t393 = t366 * t433 - t432;
t399 = t469 * t412;
t341 = -t393 * t365 - t399;
t339 = 0.1e1 / t341 ^ 2;
t476 = t315 * t339;
t360 = t468 * t467;
t408 = -t369 * t429 + t360;
t349 = t468 * t368 + t369 * t416;
t418 = t369 * t430;
t402 = t469 * t418;
t386 = -t349 * t367 - t402;
t390 = t408 * t469;
t388 = t366 * t390;
t319 = t388 - t386;
t317 = t319 ^ 2;
t311 = t317 * t339 + 0.1e1;
t309 = 0.1e1 / t311;
t397 = t408 * t367;
t435 = t349 * t469;
t380 = -t366 * t397 + t435;
t335 = t385 * qJD(1) + t349 * qJD(2);
t336 = -qJD(1) * t400 - qJD(2) * t360 + (qJD(2) * t466 + qJD(1)) * t447;
t356 = t367 * t418;
t396 = t469 * t410;
t434 = t366 * t469;
t389 = -qJD(1) * t396 - qJD(3) * t356 + t335 * t434 + t336 * t367;
t294 = t380 * qJD(3) + t389;
t338 = 0.1e1 / t341;
t452 = t319 * t339;
t406 = -t294 * t338 + t315 * t452;
t276 = t406 * t309;
t312 = atan2(-t319, t341);
t307 = sin(t312);
t308 = cos(t312);
t411 = -t307 * t341 - t308 * t319;
t271 = t411 * t276 - t307 * t294 + t308 * t315;
t288 = -t307 * t319 + t308 * t341;
t285 = 0.1e1 / t288;
t286 = 0.1e1 / t288 ^ 2;
t383 = t385 * t469;
t474 = -t366 * t383 + t367 * t384 + t396;
t318 = t474 ^ 2;
t282 = t286 * t318 + 0.1e1;
t280 = 0.1e1 / t282;
t459 = t280 * t286;
t334 = t349 * qJD(1) + t385 * qJD(2);
t378 = qJD(1) * t408 + t384 * qJD(2);
t376 = t378 * t469;
t292 = -qJD(1) * t402 + t325 * qJD(3) - t334 * t367 - t366 * t376;
t457 = t286 * t474;
t463 = t271 * t285 * t286;
t464 = (-t292 * t457 - t318 * t463) / t282 ^ 2;
t475 = -t271 * t459 - 0.2e1 * t285 * t464;
t470 = -0.2e1 * t474;
t427 = t463 * t470;
t446 = 0.2e1 * t464;
t473 = t280 * t427 - t292 * t459 - t446 * t457;
t472 = (qJD(1) * t410 - t349 * qJD(3) - t335 * t366) * t367 - qJD(3) * t402 + t336 * t469;
t306 = t325 * t362 + t344 * t361;
t300 = 0.1e1 / t306;
t301 = 0.1e1 / t306 ^ 2;
t471 = -0.2e1 * t319;
t363 = qJD(4) + qJD(5);
t436 = t369 * t450;
t421 = -qJD(1) * t436 + t325 * t363 + t378 * t465;
t377 = t378 * t367;
t293 = qJD(1) * t356 + qJD(3) * t474 - t334 * t469 + t366 * t377;
t424 = t344 * t363 + t293;
t283 = t424 * t361 + t421 * t362;
t299 = t305 ^ 2;
t291 = t299 * t301 + 0.1e1;
t455 = t301 * t305;
t284 = -t421 * t361 + t424 * t362;
t458 = t284 * t300 * t301;
t462 = (t283 * t455 - t299 * t458) / t291 ^ 2;
t454 = t338 * t476;
t461 = (t294 * t452 - t317 * t454) / t311 ^ 2;
t460 = t280 * t285;
t289 = 0.1e1 / t291;
t456 = t289 * t301;
t453 = t319 * t338;
t449 = t366 * t367;
t445 = -0.2e1 * t462;
t444 = -0.2e1 * t461;
t442 = t301 * t462;
t441 = t338 * t461;
t439 = t280 * t457;
t438 = t283 * t456;
t437 = t305 * t458;
t431 = t384 * t465;
t426 = 0.2e1 * t437;
t425 = t454 * t471;
t343 = -t408 * t465 + t436;
t398 = t366 * t408;
t387 = t469 * t398;
t423 = qJD(3) * t387 + t343 * t363 - t472;
t381 = t367 * t398 - t435;
t323 = t356 + t381;
t422 = qJD(1) * t420 + t323 * t363 + t335 * t465;
t332 = -t385 * t367 - t384 * t434;
t407 = -t332 * qJD(3) + t334 * t449 - t363 * t431 + t376;
t405 = -t300 * t361 + t362 * t455;
t321 = -t356 + t380;
t342 = t394 * t365 + t409;
t404 = -t321 * t338 + t342 * t452;
t331 = t349 * t434 - t397;
t348 = t392 * t365;
t403 = -t331 * t338 + t348 * t452;
t333 = t384 * t449 - t383;
t401 = t333 * t363 + t465 * t334;
t395 = -t307 + (t308 * t453 + t307) * t309;
t391 = -t366 * t432 + t433;
t328 = (t393 * qJD(2) + t391 * qJD(3)) * t365;
t316 = qJD(3) * t399 + (t391 * qJD(2) + t393 * qJD(3)) * t365;
t314 = t333 * t362 - t361 * t431;
t304 = t323 * t362 + t343 * t361;
t303 = t323 * t361 - t343 * t362;
t298 = t336 * t434 - t335 * t367 + (-t349 * t449 - t390) * qJD(3);
t295 = -qJD(3) * t388 + t472;
t279 = t403 * t309;
t278 = t404 * t309;
t272 = t411 * t278 - t307 * t321 + t308 * t342;
t270 = t403 * t444 + (t348 * t425 - t298 * t338 + (t294 * t348 + t315 * t331 + t319 * t328) * t339) * t309;
t269 = t404 * t444 + (t342 * t425 - t295 * t338 + (t294 * t342 + t315 * t321 + t316 * t319) * t339) * t309;
t267 = t445 + (t438 + (-t289 * t458 - t442) * t305) * t477;
t1 = [t441 * t470 + (-t292 * t338 - t474 * t476) * t309, t270, t269, 0, 0, 0; (t381 * qJD(3) - t389) * t460 + (t395 * t292 - ((-t276 * t309 * t453 + t444) * t307 + (t441 * t471 - t276 + (t276 - t406) * t309) * t308) * t474) * t439 + t475 * (-t387 + t386) - t473 * t395 * t474 (t333 * qJD(3) - t334 * t434 + t377) * t460 + ((-t270 * t319 - t279 * t294 + t328 + (-t279 * t341 - t331) * t276) * t308 + (-t270 * t341 - t279 * t315 - t298 + (t279 * t319 - t348) * t276) * t307) * t439 + t475 * t332 + t473 * (t411 * t279 - t307 * t331 + t308 * t348) (-t272 * t457 - t285 * t325) * t446 + (t272 * t427 + t293 * t285 + (-t325 * t271 - t272 * t292 - (-(-t269 * t319 - t278 * t294 + t316 + (-t278 * t341 - t321) * t276) * t308 - (-t269 * t341 - t278 * t315 - t295 + (t278 * t319 - t342) * t276) * t307) * t474) * t286) * t280, 0, 0, 0; 0.2e1 * (-t300 * t303 + t304 * t455) * t462 + ((t423 * t361 + t422 * t362) * t300 + t304 * t426 + (-t303 * t284 - (-t422 * t361 + t423 * t362) * t305 - t304 * t283) * t301) * t289 (t442 * t477 - t438) * t314 + (-t284 * t456 + t300 * t445) * (t333 * t361 + t362 * t431) + ((t407 * t361 + t401 * t362) * t300 - (-t401 * t361 + t407 * t362) * t455 + t314 * t426) * t289, -t405 * t474 * t445 + (t405 * t292 - ((-t300 * t363 - 0.2e1 * t437) * t362 + (t283 * t362 + (-t305 * t363 + t284) * t361) * t301) * t474) * t289, t267, t267, 0;];
JaD_rot  = t1;
