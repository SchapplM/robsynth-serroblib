% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP12_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:14:14
% EndTime: 2019-02-26 21:14:17
% DurationCPUTime: 3.55s
% Computational Cost: add. (13786->165), mult. (41237->320), div. (726->12), fcn. (53008->17), ass. (0->146)
t482 = sin(pkin(12));
t487 = cos(pkin(6));
t442 = t487 * t482;
t485 = cos(pkin(12));
t488 = sin(qJ(1));
t490 = cos(qJ(1));
t383 = -t488 * t442 + t490 * t485;
t380 = t383 * qJD(1);
t392 = sin(qJ(3));
t394 = cos(qJ(3));
t382 = t490 * t442 + t488 * t485;
t483 = sin(pkin(7));
t484 = sin(pkin(6));
t440 = t484 * t483;
t427 = t490 * t440;
t486 = cos(pkin(7));
t443 = t487 * t485;
t497 = -t490 * t443 + t488 * t482;
t498 = t497 * t486;
t404 = t498 + t427;
t493 = t382 * t392 + t404 * t394;
t413 = t488 * t443 + t490 * t482;
t379 = t413 * qJD(1);
t424 = t488 * t440;
t503 = qJD(1) * t424 - t379 * t486;
t340 = qJD(3) * t493 - t380 * t394 - t503 * t392;
t391 = sin(qJ(4));
t441 = t486 * t484;
t425 = t488 * t441;
t414 = qJD(1) * t425 + t379 * t483;
t489 = cos(qJ(4));
t463 = t382 * t394;
t363 = t404 * t392 - t463;
t405 = -t490 * t441 + t483 * t497;
t512 = -t363 * t391 - t405 * t489;
t320 = qJD(4) * t512 + t340 * t489 - t414 * t391;
t351 = t363 * t489 - t405 * t391;
t514 = t351 * qJD(4) + t340 * t391 + t414 * t489;
t345 = t512 ^ 2;
t410 = t485 * t441 + t487 * t483;
t439 = t484 * t482;
t374 = t410 * t392 + t394 * t439;
t381 = -t485 * t440 + t487 * t486;
t430 = -t374 * t391 + t381 * t489;
t356 = 0.1e1 / t430 ^ 2;
t333 = t345 * t356 + 0.1e1;
t331 = 0.1e1 / t333;
t359 = t374 * t489 + t381 * t391;
t373 = -t392 * t439 + t410 * t394;
t368 = t373 * qJD(3);
t343 = t359 * qJD(4) + t368 * t391;
t355 = 0.1e1 / t430;
t468 = t512 * t356;
t435 = t343 * t468 - t355 * t514;
t300 = t435 * t331;
t334 = atan2(-t512, -t430);
t329 = sin(t334);
t330 = cos(t334);
t438 = t329 * t430 - t330 * t512;
t295 = t438 * t300 + t329 * t514 + t330 * t343;
t312 = -t329 * t512 - t330 * t430;
t310 = 0.1e1 / t312 ^ 2;
t511 = t295 * t310;
t508 = -t413 * t486 + t424;
t365 = t383 * t394 + t508 * t392;
t403 = t413 * t483 + t425;
t353 = t365 * t489 + t403 * t391;
t495 = t508 * t394;
t364 = t383 * t392 - t495;
t390 = sin(qJ(5));
t393 = cos(qJ(5));
t327 = t353 * t390 - t364 * t393;
t507 = 0.2e1 * t327;
t309 = 0.1e1 / t312;
t506 = t309 * t511;
t401 = t403 * t489;
t352 = t365 * t391 - t401;
t491 = 0.2e1 * t352;
t447 = t491 * t506;
t378 = t382 * qJD(1);
t462 = qJD(3) * t392;
t496 = t404 * qJD(1);
t336 = qJD(3) * t495 - t378 * t394 - t383 * t462 + t392 * t496;
t402 = t405 * qJD(1);
t316 = t353 * qJD(4) + t336 * t391 + t489 * t402;
t474 = t316 * t310;
t502 = -t474 + t447;
t501 = (t392 * t498 - t463) * qJD(3) - t380 * t392 + t503 * t394 + t427 * t462;
t346 = t352 ^ 2;
t306 = t346 * t310 + 0.1e1;
t459 = 0.2e1 * (-t346 * t506 + t352 * t474) / t306 ^ 2;
t500 = t343 * t356;
t432 = -t355 * t493 + t373 * t468;
t499 = t391 * t432;
t449 = qJD(4) * t489;
t328 = t353 * t393 + t364 * t390;
t322 = 0.1e1 / t328;
t323 = 0.1e1 / t328 ^ 2;
t492 = -0.2e1 * t512;
t461 = qJD(4) * t391;
t317 = qJD(4) * t401 + t336 * t489 - t365 * t461 - t391 * t402;
t335 = t365 * qJD(3) - t378 * t392 - t394 * t496;
t307 = t328 * qJD(5) + t317 * t390 - t335 * t393;
t321 = t327 ^ 2;
t315 = t321 * t323 + 0.1e1;
t473 = t323 * t327;
t460 = qJD(5) * t327;
t308 = t317 * t393 + t335 * t390 - t460;
t477 = t308 * t322 * t323;
t480 = (t307 * t473 - t321 * t477) / t315 ^ 2;
t470 = t355 * t500;
t478 = (t345 * t470 - t468 * t514) / t333 ^ 2;
t476 = t310 * t352;
t313 = 0.1e1 / t315;
t475 = t313 * t323;
t472 = t329 * t352;
t471 = t330 * t352;
t469 = t512 * t355;
t466 = t364 * t391;
t458 = -0.2e1 * t480;
t457 = -0.2e1 * t478;
t455 = t323 * t480;
t454 = t355 * t478;
t453 = t307 * t475;
t452 = t327 * t477;
t450 = t364 * t489;
t446 = 0.2e1 * t452;
t445 = t470 * t492;
t326 = t351 * t393 - t390 * t493;
t325 = t351 * t390 + t393 * t493;
t436 = qJD(5) * t450 + t336;
t434 = -t390 * t322 + t393 * t473;
t433 = -t351 * t355 + t359 * t468;
t423 = -t329 + (-t330 * t469 + t329) * t331;
t417 = qJD(5) * t365 - t489 * t335 + t364 * t461;
t369 = t374 * qJD(3);
t344 = t430 * qJD(4) + t368 * t489;
t342 = t365 * t390 - t393 * t450;
t304 = 0.1e1 / t306;
t303 = t331 * t499;
t301 = t433 * t331;
t297 = (t329 * t493 + t330 * t373) * t391 + t438 * t303;
t296 = t438 * t301 + t329 * t351 + t330 * t359;
t293 = t433 * t457 + (-t359 * t445 - t320 * t355 + (-t343 * t351 + t344 * t512 - t359 * t514) * t356) * t331;
t292 = t457 * t499 + ((-t373 * t445 + t501 * t355 + (-t343 * t493 - t369 * t512 - t373 * t514) * t356) * t391 + t432 * t449) * t331;
t1 = [-t454 * t491 + (t316 * t355 + t352 * t500) * t331, 0, t292, t293, 0, 0; t512 * t309 * t459 + (t514 * t309 + t512 * t511 - (t423 * t316 + ((t300 * t331 * t469 + t457) * t329 + (-t454 * t492 - t300 + (t300 - t435) * t331) * t330) * t352) * t476) * t304 + (t502 * t304 + t476 * t459) * t423 * t352, 0 (t297 * t476 + t309 * t466) * t459 + ((-t335 * t391 - t364 * t449) * t309 + t502 * t297 + (t466 * t295 - (t373 * t449 - t292 * t512 + t303 * t514 - t369 * t391 + (t303 * t430 + t391 * t493) * t300) * t471 - (t493 * t449 + t292 * t430 - t303 * t343 - t501 * t391 + (t303 * t512 - t373 * t391) * t300) * t472) * t310) * t304 (t296 * t476 - t309 * t353) * t459 + (t296 * t447 + t317 * t309 + (-t353 * t295 - t296 * t316 - (-t293 * t512 + t301 * t514 + t344 + (t301 * t430 + t351) * t300) * t471 - (t293 * t430 - t301 * t343 + t320 + (t301 * t512 - t359) * t300) * t472) * t310) * t304, 0, 0; 0.2e1 * (-t322 * t325 + t326 * t473) * t480 + ((t326 * qJD(5) + t320 * t390 - t393 * t501) * t322 + t326 * t446 + (-t325 * t308 - (-t325 * qJD(5) + t320 * t393 + t390 * t501) * t327 - t326 * t307) * t323) * t313, 0 (t455 * t507 - t453) * t342 + (-t308 * t475 + t322 * t458) * (-t365 * t393 - t390 * t450) + ((t417 * t390 - t436 * t393) * t322 - (t436 * t390 + t417 * t393) * t473 + t342 * t446) * t313, t434 * t352 * t458 + (t434 * t316 + ((-qJD(5) * t322 - 0.2e1 * t452) * t393 + (t307 * t393 + (t308 - t460) * t390) * t323) * t352) * t313, t458 + (t453 + (-t313 * t477 - t455) * t327) * t507, 0;];
JaD_rot  = t1;
