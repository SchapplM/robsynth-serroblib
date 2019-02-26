% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP11
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
% Datum: 2019-02-26 21:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP11_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:33
% EndTime: 2019-02-26 21:13:36
% DurationCPUTime: 3.46s
% Computational Cost: add. (13786->165), mult. (41237->320), div. (726->12), fcn. (53008->17), ass. (0->146)
t493 = sin(pkin(12));
t498 = cos(pkin(6));
t453 = t498 * t493;
t496 = cos(pkin(12));
t499 = sin(qJ(1));
t501 = cos(qJ(1));
t394 = -t499 * t453 + t501 * t496;
t391 = t394 * qJD(1);
t403 = sin(qJ(3));
t405 = cos(qJ(3));
t393 = t501 * t453 + t499 * t496;
t494 = sin(pkin(7));
t495 = sin(pkin(6));
t451 = t495 * t494;
t438 = t501 * t451;
t497 = cos(pkin(7));
t454 = t498 * t496;
t508 = -t501 * t454 + t499 * t493;
t509 = t508 * t497;
t415 = t509 + t438;
t504 = t393 * t403 + t415 * t405;
t424 = t499 * t454 + t501 * t493;
t390 = t424 * qJD(1);
t435 = t499 * t451;
t514 = qJD(1) * t435 - t390 * t497;
t351 = qJD(3) * t504 - t391 * t405 - t514 * t403;
t402 = sin(qJ(4));
t452 = t497 * t495;
t436 = t499 * t452;
t425 = qJD(1) * t436 + t390 * t494;
t500 = cos(qJ(4));
t474 = t393 * t405;
t374 = t415 * t403 - t474;
t416 = -t501 * t452 + t494 * t508;
t523 = -t374 * t402 - t416 * t500;
t331 = qJD(4) * t523 + t351 * t500 - t425 * t402;
t362 = t374 * t500 - t416 * t402;
t525 = t362 * qJD(4) + t351 * t402 + t425 * t500;
t356 = t523 ^ 2;
t421 = t496 * t452 + t498 * t494;
t450 = t495 * t493;
t385 = t421 * t403 + t405 * t450;
t392 = -t496 * t451 + t498 * t497;
t441 = -t385 * t402 + t392 * t500;
t367 = 0.1e1 / t441 ^ 2;
t344 = t356 * t367 + 0.1e1;
t342 = 0.1e1 / t344;
t370 = t385 * t500 + t392 * t402;
t384 = -t403 * t450 + t421 * t405;
t379 = t384 * qJD(3);
t354 = t370 * qJD(4) + t379 * t402;
t366 = 0.1e1 / t441;
t479 = t523 * t367;
t446 = t354 * t479 - t366 * t525;
t311 = t446 * t342;
t345 = atan2(-t523, -t441);
t340 = sin(t345);
t341 = cos(t345);
t449 = t340 * t441 - t341 * t523;
t306 = t449 * t311 + t340 * t525 + t341 * t354;
t323 = -t340 * t523 - t341 * t441;
t321 = 0.1e1 / t323 ^ 2;
t522 = t306 * t321;
t519 = -t424 * t497 + t435;
t376 = t394 * t405 + t519 * t403;
t414 = t424 * t494 + t436;
t364 = t376 * t500 + t414 * t402;
t506 = t519 * t405;
t375 = t394 * t403 - t506;
t401 = sin(qJ(5));
t404 = cos(qJ(5));
t338 = t364 * t401 - t375 * t404;
t518 = 0.2e1 * t338;
t320 = 0.1e1 / t323;
t517 = t320 * t522;
t412 = t414 * t500;
t363 = t376 * t402 - t412;
t502 = 0.2e1 * t363;
t456 = t502 * t517;
t389 = t393 * qJD(1);
t473 = qJD(3) * t403;
t507 = t415 * qJD(1);
t347 = qJD(3) * t506 - t389 * t405 - t394 * t473 + t403 * t507;
t413 = t416 * qJD(1);
t327 = t364 * qJD(4) + t347 * t402 + t500 * t413;
t485 = t327 * t321;
t513 = -t485 + t456;
t512 = (t403 * t509 - t474) * qJD(3) - t391 * t403 + t514 * t405 + t438 * t473;
t357 = t363 ^ 2;
t317 = t357 * t321 + 0.1e1;
t470 = 0.2e1 * (-t357 * t517 + t363 * t485) / t317 ^ 2;
t511 = t354 * t367;
t443 = -t366 * t504 + t384 * t479;
t510 = t402 * t443;
t460 = qJD(4) * t500;
t339 = t364 * t404 + t375 * t401;
t333 = 0.1e1 / t339;
t334 = 0.1e1 / t339 ^ 2;
t503 = -0.2e1 * t523;
t472 = qJD(4) * t402;
t328 = qJD(4) * t412 + t347 * t500 - t376 * t472 - t402 * t413;
t346 = t376 * qJD(3) - t389 * t403 - t405 * t507;
t318 = t339 * qJD(5) + t328 * t401 - t346 * t404;
t332 = t338 ^ 2;
t326 = t332 * t334 + 0.1e1;
t484 = t334 * t338;
t471 = qJD(5) * t338;
t319 = t328 * t404 + t346 * t401 - t471;
t488 = t319 * t333 * t334;
t491 = (t318 * t484 - t332 * t488) / t326 ^ 2;
t481 = t366 * t511;
t489 = (t356 * t481 - t479 * t525) / t344 ^ 2;
t487 = t321 * t363;
t324 = 0.1e1 / t326;
t486 = t324 * t334;
t483 = t340 * t363;
t482 = t341 * t363;
t480 = t523 * t366;
t477 = t375 * t402;
t469 = -0.2e1 * t491;
t468 = -0.2e1 * t489;
t466 = t334 * t491;
t465 = t366 * t489;
t464 = t318 * t486;
t463 = t338 * t488;
t461 = t375 * t500;
t458 = 0.2e1 * t463;
t457 = t481 * t503;
t337 = t362 * t404 - t401 * t504;
t336 = t362 * t401 + t404 * t504;
t447 = qJD(5) * t461 + t347;
t445 = -t401 * t333 + t404 * t484;
t444 = -t362 * t366 + t370 * t479;
t434 = -t340 + (-t341 * t480 + t340) * t342;
t428 = qJD(5) * t376 - t500 * t346 + t375 * t472;
t380 = t385 * qJD(3);
t355 = t441 * qJD(4) + t379 * t500;
t353 = t376 * t401 - t404 * t461;
t315 = 0.1e1 / t317;
t314 = t342 * t510;
t312 = t444 * t342;
t308 = (t340 * t504 + t341 * t384) * t402 + t449 * t314;
t307 = t449 * t312 + t340 * t362 + t341 * t370;
t304 = t444 * t468 + (-t370 * t457 - t331 * t366 + (-t354 * t362 + t355 * t523 - t370 * t525) * t367) * t342;
t303 = t468 * t510 + ((-t384 * t457 + t512 * t366 + (-t354 * t504 - t380 * t523 - t384 * t525) * t367) * t402 + t443 * t460) * t342;
t1 = [-t465 * t502 + (t327 * t366 + t363 * t511) * t342, 0, t303, t304, 0, 0; t523 * t320 * t470 + (t525 * t320 + t523 * t522 - (t434 * t327 + ((t311 * t342 * t480 + t468) * t340 + (-t465 * t503 - t311 + (t311 - t446) * t342) * t341) * t363) * t487) * t315 + (t513 * t315 + t487 * t470) * t434 * t363, 0 (t308 * t487 + t320 * t477) * t470 + ((-t346 * t402 - t375 * t460) * t320 + t513 * t308 + (t477 * t306 - (t384 * t460 - t303 * t523 + t314 * t525 - t380 * t402 + (t314 * t441 + t402 * t504) * t311) * t482 - (t504 * t460 + t303 * t441 - t314 * t354 - t512 * t402 + (t314 * t523 - t384 * t402) * t311) * t483) * t321) * t315 (t307 * t487 - t320 * t364) * t470 + (t307 * t456 + t328 * t320 + (-t364 * t306 - t307 * t327 - (-t304 * t523 + t312 * t525 + t355 + (t312 * t441 + t362) * t311) * t482 - (t304 * t441 - t312 * t354 + t331 + (t312 * t523 - t370) * t311) * t483) * t321) * t315, 0, 0; 0.2e1 * (-t333 * t336 + t337 * t484) * t491 + ((t337 * qJD(5) + t331 * t401 - t404 * t512) * t333 + t337 * t458 + (-t336 * t319 - (-t336 * qJD(5) + t331 * t404 + t401 * t512) * t338 - t337 * t318) * t334) * t324, 0 (t466 * t518 - t464) * t353 + (-t319 * t486 + t333 * t469) * (-t376 * t404 - t401 * t461) + ((t428 * t401 - t447 * t404) * t333 - (t447 * t401 + t428 * t404) * t484 + t353 * t458) * t324, t445 * t363 * t469 + (t445 * t327 + ((-qJD(5) * t333 - 0.2e1 * t463) * t404 + (t318 * t404 + (t319 - t471) * t401) * t334) * t363) * t324, t469 + (t464 + (-t324 * t488 - t466) * t338) * t518, 0;];
JaD_rot  = t1;
