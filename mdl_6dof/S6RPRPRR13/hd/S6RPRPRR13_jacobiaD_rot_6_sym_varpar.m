% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR13_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:44
% EndTime: 2019-02-26 20:55:47
% DurationCPUTime: 3.46s
% Computational Cost: add. (13786->169), mult. (41237->335), div. (726->12), fcn. (53008->17), ass. (0->145)
t388 = sin(qJ(3));
t478 = sin(pkin(12));
t482 = cos(pkin(6));
t436 = t482 * t478;
t480 = cos(pkin(12));
t483 = sin(qJ(1));
t485 = cos(qJ(1));
t411 = t485 * t436 + t483 * t480;
t484 = cos(qJ(3));
t406 = t411 * t484;
t385 = sin(pkin(6));
t479 = sin(pkin(7));
t441 = t385 * t479;
t430 = t485 * t441;
t481 = cos(pkin(7));
t438 = t482 * t480;
t493 = -t485 * t438 + t483 * t478;
t494 = t493 * t481;
t364 = (t494 + t430) * t388 - t406;
t413 = -t483 * t436 + t485 * t480;
t379 = t413 * qJD(1);
t428 = t483 * t441;
t382 = t484 * t428;
t412 = t483 * t438 + t485 * t478;
t378 = t412 * qJD(1);
t443 = t378 * t481;
t339 = qJD(1) * t382 + t364 * qJD(3) - t379 * t388 - t484 * t443;
t442 = t385 * t481;
t429 = t483 * t442;
t370 = qJD(1) * t429 + t378 * t479;
t387 = sin(qJ(5));
t390 = cos(qJ(5));
t374 = -t485 * t442 + t479 * t493;
t498 = -t411 * t388 - t484 * t430;
t399 = t484 * t494 - t498;
t504 = t374 * t387 - t399 * t390;
t318 = qJD(5) * t504 + t339 * t387 - t370 * t390;
t351 = t374 * t390 + t399 * t387;
t319 = t351 * qJD(5) + t339 * t390 + t370 * t387;
t340 = t399 * qJD(3) - (qJD(1) * t428 - t443) * t388 - t379 * t484;
t375 = t412 * t479 + t429;
t405 = t412 * t481;
t403 = t413 * t388 + t484 * t405 - t382;
t353 = t375 * t390 + t403 * t387;
t365 = t413 * t484 + (-t405 + t428) * t388;
t386 = sin(qJ(6));
t389 = cos(qJ(6));
t327 = t353 * t386 - t365 * t389;
t501 = 0.2e1 * t327;
t380 = -t480 * t441 + t482 * t481;
t435 = t481 * t480;
t437 = t482 * t479;
t488 = (-t478 * t388 + t484 * t435) * t385 + t484 * t437;
t432 = -t380 * t387 - t390 * t488;
t334 = atan2(-t504, -t432);
t329 = sin(t334);
t330 = cos(t334);
t312 = -t329 * t504 - t330 * t432;
t310 = 0.1e1 / t312 ^ 2;
t401 = t403 * t390;
t352 = t375 * t387 - t401;
t346 = t352 ^ 2;
t306 = t310 * t346 + 0.1e1;
t368 = t374 * qJD(1);
t404 = qJD(1) * t494;
t397 = t498 * qJD(1) + t365 * qJD(3) - t484 * t404;
t316 = t353 * qJD(5) - t368 * t387 - t397 * t390;
t472 = t310 * t352;
t345 = t504 ^ 2;
t358 = 0.1e1 / t432 ^ 2;
t333 = t345 * t358 + 0.1e1;
t331 = 0.1e1 / t333;
t361 = t380 * t390 - t387 * t488;
t372 = t388 * t437 + (t388 * t435 + t478 * t484) * t385;
t367 = t372 * qJD(3);
t344 = t361 * qJD(5) - t367 * t390;
t357 = 0.1e1 / t432;
t464 = t504 * t358;
t427 = t319 * t357 + t344 * t464;
t300 = t427 * t331;
t434 = t329 * t432 - t330 * t504;
t295 = t434 * t300 - t329 * t319 + t330 * t344;
t309 = 0.1e1 / t312;
t311 = t309 * t310;
t476 = t295 * t311;
t457 = 0.2e1 * (t316 * t472 - t346 * t476) / t306 ^ 2;
t496 = t344 * t358;
t424 = -t357 * t364 + t372 * t464;
t495 = t390 * t424;
t492 = t365 * t387 * qJD(6) + t397;
t463 = t365 * t390;
t490 = qJD(5) * t463 - qJD(6) * t403;
t328 = t353 * t389 + t365 * t386;
t322 = 0.1e1 / t328;
t323 = 0.1e1 / t328 ^ 2;
t486 = 0.2e1 * t352;
t466 = t357 * t496;
t475 = (t319 * t464 + t345 * t466) / t333 ^ 2;
t459 = qJD(5) * t387;
t317 = qJD(5) * t401 - t368 * t390 - t375 * t459 + t397 * t387;
t337 = t388 * t404 + (t388 * t430 - t406) * qJD(1) - t403 * qJD(3);
t458 = qJD(6) * t327;
t308 = t317 * t389 + t337 * t386 - t458;
t474 = t308 * t322 * t323;
t473 = t310 * t316;
t321 = t327 ^ 2;
t315 = t321 * t323 + 0.1e1;
t313 = 0.1e1 / t315;
t471 = t313 * t323;
t307 = t328 * qJD(6) + t317 * t386 - t337 * t389;
t469 = t323 * t327;
t470 = 0.1e1 / t315 ^ 2 * (t307 * t469 - t321 * t474);
t468 = t329 * t352;
t467 = t330 * t352;
t465 = t504 * t357;
t462 = t386 * t387;
t461 = t387 * t389;
t456 = -0.2e1 * t475;
t455 = t311 * t486;
t454 = -0.2e1 * t470;
t453 = t323 * t470;
t452 = t357 * t475;
t451 = t307 * t471;
t450 = t327 * t474;
t449 = t310 * t468;
t448 = t310 * t467;
t447 = t504 * t466;
t440 = 0.2e1 * t450;
t326 = -t351 * t389 + t364 * t386;
t325 = -t351 * t386 - t364 * t389;
t426 = -t322 * t386 + t389 * t469;
t425 = t351 * t357 + t361 * t464;
t419 = -t329 + (-t330 * t465 + t329) * t331;
t366 = t488 * qJD(3);
t343 = t432 * qJD(5) + t367 * t387;
t342 = t365 * t461 - t403 * t386;
t304 = 0.1e1 / t306;
t303 = t331 * t495;
t301 = t425 * t331;
t299 = t419 * t352;
t297 = (-t329 * t364 - t330 * t372) * t390 - t434 * t303;
t296 = t434 * t301 - t329 * t351 + t330 * t361;
t293 = t425 * t456 + (0.2e1 * t361 * t447 - t318 * t357 + (t319 * t361 + t343 * t504 + t344 * t351) * t358) * t331;
t292 = 0.2e1 * t475 * t495 + (t424 * t459 + (-0.2e1 * t372 * t447 + t340 * t357 + (-t319 * t372 + t344 * t364 - t366 * t504) * t358) * t390) * t331;
t1 = [-t452 * t486 + (t316 * t357 + t352 * t496) * t331, 0, t292, 0, t293, 0; t504 * t309 * t457 + (-t319 * t309 + (t295 * t504 - t299 * t316) * t310) * t304 + (t299 * t310 * t457 + (0.2e1 * t299 * t476 - (t300 * t331 * t465 + t456) * t449 - (0.2e1 * t504 * t452 - t300 + (t300 - t427) * t331) * t448 - t419 * t473) * t304) * t352, 0 (t297 * t472 + t309 * t463) * t457 + (-t297 * t473 + (-t337 * t390 + t365 * t459) * t309 + (t297 * t455 + t310 * t463) * t295 - (t372 * t459 - t292 * t504 + t303 * t319 - t366 * t390 + (-t303 * t432 - t364 * t390) * t300) * t448 - (t364 * t459 + t292 * t432 + t303 * t344 - t340 * t390 + (-t303 * t504 + t372 * t390) * t300) * t449) * t304, 0 (t296 * t472 - t309 * t353) * t457 + (t296 * t295 * t455 + t317 * t309 + (-t353 * t295 - t296 * t316 - (-t293 * t504 - t301 * t319 + t343 + (t301 * t432 - t351) * t300) * t467 - (t293 * t432 - t301 * t344 + t318 + (t301 * t504 - t361) * t300) * t468) * t310) * t304, 0; 0.2e1 * (-t322 * t325 + t326 * t469) * t470 + ((t326 * qJD(6) + t318 * t386 - t340 * t389) * t322 + t326 * t440 + (-t325 * t308 - (-t325 * qJD(6) + t318 * t389 + t340 * t386) * t327 - t326 * t307) * t323) * t313, 0 (t453 * t501 - t451) * t342 + (-t308 * t471 + t322 * t454) * (t365 * t462 + t403 * t389) + (t342 * t440 + (t462 * t322 - t461 * t469) * t337 + (t492 * t322 - t490 * t469) * t389 + (t490 * t322 + t492 * t469) * t386) * t313, 0, t426 * t352 * t454 + (t426 * t316 + ((-qJD(6) * t322 - 0.2e1 * t450) * t389 + (t307 * t389 + (t308 - t458) * t386) * t323) * t352) * t313, t454 + (t451 + (-t313 * t474 - t453) * t327) * t501;];
JaD_rot  = t1;
