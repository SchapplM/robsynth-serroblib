% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR11_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:32
% EndTime: 2019-02-26 20:54:36
% DurationCPUTime: 3.58s
% Computational Cost: add. (17692->166), mult. (41237->320), div. (726->12), fcn. (53008->17), ass. (0->146)
t495 = sin(pkin(12));
t500 = cos(pkin(6));
t455 = t500 * t495;
t498 = cos(pkin(12));
t501 = sin(qJ(1));
t502 = cos(qJ(1));
t394 = -t501 * t455 + t502 * t498;
t391 = t394 * qJD(1);
t403 = sin(qJ(3));
t405 = cos(qJ(3));
t393 = t502 * t455 + t501 * t498;
t496 = sin(pkin(7));
t497 = sin(pkin(6));
t453 = t497 * t496;
t441 = t502 * t453;
t499 = cos(pkin(7));
t456 = t500 * t498;
t509 = -t502 * t456 + t501 * t495;
t510 = t509 * t499;
t415 = t510 + t441;
t505 = t393 * t403 + t415 * t405;
t424 = t501 * t456 + t502 * t495;
t390 = t424 * qJD(1);
t438 = t501 * t453;
t515 = qJD(1) * t438 - t390 * t499;
t353 = t505 * qJD(3) - t391 * t405 - t515 * t403;
t472 = pkin(13) + qJ(5);
t401 = sin(t472);
t454 = t499 * t497;
t439 = t501 * t454;
t425 = qJD(1) * t439 + t390 * t496;
t461 = cos(t472);
t476 = t393 * t405;
t374 = t415 * t403 - t476;
t416 = -t502 * t454 + t496 * t509;
t524 = -t374 * t401 - t416 * t461;
t331 = qJD(5) * t524 + t353 * t461 - t425 * t401;
t362 = t374 * t461 - t416 * t401;
t526 = t362 * qJD(5) + t353 * t401 + t425 * t461;
t356 = t524 ^ 2;
t421 = t498 * t454 + t500 * t496;
t452 = t497 * t495;
t385 = t421 * t403 + t405 * t452;
t392 = -t498 * t453 + t500 * t499;
t431 = -t385 * t401 + t392 * t461;
t366 = 0.1e1 / t431 ^ 2;
t340 = t356 * t366 + 0.1e1;
t334 = 0.1e1 / t340;
t369 = t385 * t461 + t392 * t401;
t384 = -t403 * t452 + t421 * t405;
t379 = t384 * qJD(3);
t354 = t369 * qJD(5) + t379 * t401;
t365 = 0.1e1 / t431;
t481 = t524 * t366;
t447 = t354 * t481 - t365 * t526;
t311 = t447 * t334;
t341 = atan2(-t524, -t431);
t332 = sin(t341);
t333 = cos(t341);
t449 = t332 * t431 - t333 * t524;
t306 = t449 * t311 + t332 * t526 + t333 * t354;
t323 = -t332 * t524 - t333 * t431;
t321 = 0.1e1 / t323 ^ 2;
t523 = t306 * t321;
t520 = -t424 * t499 + t438;
t376 = t394 * t405 + t520 * t403;
t414 = t424 * t496 + t439;
t364 = t376 * t461 + t414 * t401;
t507 = t520 * t405;
t375 = t394 * t403 - t507;
t402 = sin(qJ(6));
t404 = cos(qJ(6));
t344 = t364 * t402 - t375 * t404;
t519 = 0.2e1 * t344;
t320 = 0.1e1 / t323;
t518 = t320 * t523;
t389 = t393 * qJD(1);
t475 = qJD(3) * t403;
t508 = t415 * qJD(1);
t349 = t507 * qJD(3) - t389 * t405 - t394 * t475 + t508 * t403;
t413 = qJD(1) * t416;
t327 = t364 * qJD(5) + t349 * t401 + t461 * t413;
t412 = t414 * t461;
t363 = t376 * t401 - t412;
t503 = 0.2e1 * t363;
t460 = t503 * t518;
t514 = -t321 * t327 + t460;
t513 = (t403 * t510 - t476) * qJD(3) - t391 * t403 + t515 * t405 + t441 * t475;
t357 = t363 ^ 2;
t317 = t321 * t357 + 0.1e1;
t489 = t321 * t363;
t471 = 0.2e1 * (t327 * t489 - t357 * t518) / t317 ^ 2;
t512 = t354 * t366;
t444 = -t365 * t505 + t384 * t481;
t511 = t401 * t444;
t450 = qJD(5) * t461;
t345 = t364 * t404 + t375 * t402;
t337 = 0.1e1 / t345;
t338 = 0.1e1 / t345 ^ 2;
t504 = -0.2e1 * t524;
t483 = t365 * t512;
t492 = (t356 * t483 - t481 * t526) / t340 ^ 2;
t474 = qJD(5) * t401;
t328 = qJD(5) * t412 + t349 * t461 - t376 * t474 - t401 * t413;
t348 = t376 * qJD(3) - t389 * t403 - t508 * t405;
t473 = qJD(6) * t344;
t319 = t328 * t404 + t348 * t402 - t473;
t491 = t319 * t337 * t338;
t336 = t344 ^ 2;
t326 = t336 * t338 + 0.1e1;
t324 = 0.1e1 / t326;
t488 = t324 * t338;
t318 = t345 * qJD(6) + t328 * t402 - t348 * t404;
t484 = t338 * t344;
t487 = 0.1e1 / t326 ^ 2 * (t318 * t484 - t336 * t491);
t486 = t332 * t363;
t485 = t333 * t363;
t482 = t524 * t365;
t479 = t375 * t401;
t470 = -0.2e1 * t492;
t468 = -0.2e1 * t487;
t467 = t338 * t487;
t466 = t365 * t492;
t465 = t318 * t488;
t464 = t344 * t491;
t459 = 0.2e1 * t464;
t458 = t483 * t504;
t451 = t375 * t461;
t343 = t362 * t404 - t402 * t505;
t342 = t362 * t402 + t404 * t505;
t446 = -t337 * t402 + t404 * t484;
t445 = -t362 * t365 + t369 * t481;
t437 = -t332 + (-t333 * t482 + t332) * t334;
t436 = qJD(6) * t451 + t349;
t426 = qJD(6) * t376 - t461 * t348 + t375 * t474;
t380 = t385 * qJD(3);
t355 = t431 * qJD(5) + t379 * t461;
t347 = t376 * t402 - t404 * t451;
t315 = 0.1e1 / t317;
t314 = t334 * t511;
t312 = t445 * t334;
t308 = (t332 * t505 + t333 * t384) * t401 + t449 * t314;
t307 = t449 * t312 + t332 * t362 + t333 * t369;
t304 = t445 * t470 + (-t369 * t458 - t331 * t365 + (-t354 * t362 + t355 * t524 - t369 * t526) * t366) * t334;
t303 = t470 * t511 + ((-t384 * t458 + t513 * t365 + (-t354 * t505 - t380 * t524 - t384 * t526) * t366) * t401 + t444 * t450) * t334;
t1 = [-t466 * t503 + (t327 * t365 + t363 * t512) * t334, 0, t303, 0, t304, 0; t524 * t320 * t471 + (t526 * t320 + t524 * t523 - (t437 * t327 + ((t311 * t334 * t482 + t470) * t332 + (-t466 * t504 - t311 + (t311 - t447) * t334) * t333) * t363) * t489) * t315 + (t514 * t315 + t489 * t471) * t437 * t363, 0 (t308 * t489 + t320 * t479) * t471 + ((-t348 * t401 - t375 * t450) * t320 + t514 * t308 + (t479 * t306 - (t384 * t450 - t303 * t524 + t314 * t526 - t380 * t401 + (t314 * t431 + t401 * t505) * t311) * t485 - (t505 * t450 + t303 * t431 - t314 * t354 - t513 * t401 + (t314 * t524 - t384 * t401) * t311) * t486) * t321) * t315, 0 (t307 * t489 - t320 * t364) * t471 + (t307 * t460 + t328 * t320 + (-t364 * t306 - t307 * t327 - (-t304 * t524 + t312 * t526 + t355 + (t312 * t431 + t362) * t311) * t485 - (t304 * t431 - t312 * t354 + t331 + (t312 * t524 - t369) * t311) * t486) * t321) * t315, 0; 0.2e1 * (-t337 * t342 + t343 * t484) * t487 + ((t343 * qJD(6) + t331 * t402 - t404 * t513) * t337 + t343 * t459 + (-t342 * t319 - (-t342 * qJD(6) + t331 * t404 + t402 * t513) * t344 - t343 * t318) * t338) * t324, 0 (t467 * t519 - t465) * t347 + (-t319 * t488 + t337 * t468) * (-t376 * t404 - t402 * t451) + ((t426 * t402 - t436 * t404) * t337 - (t436 * t402 + t426 * t404) * t484 + t347 * t459) * t324, 0, t446 * t363 * t468 + (t446 * t327 + ((-qJD(6) * t337 - 0.2e1 * t464) * t404 + (t318 * t404 + (t319 - t473) * t402) * t338) * t363) * t324, t468 + (t465 + (-t324 * t491 - t467) * t344) * t519;];
JaD_rot  = t1;
