% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JaD_rot = S6RPRRRP12_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:14:19
% EndTime: 2019-02-26 21:14:25
% DurationCPUTime: 5.90s
% Computational Cost: add. (28471->207), mult. (84790->400), div. (951->12), fcn. (109142->17), ass. (0->168)
t536 = sin(pkin(12));
t540 = cos(pkin(6));
t489 = t540 * t536;
t538 = cos(pkin(12));
t541 = sin(qJ(1));
t543 = cos(qJ(1));
t429 = -t541 * t489 + t543 * t538;
t427 = t429 * qJD(1);
t441 = sin(qJ(3));
t438 = sin(pkin(6));
t537 = sin(pkin(7));
t501 = t438 * t537;
t483 = t541 * t501;
t491 = t540 * t538;
t462 = t541 * t491 + t543 * t536;
t426 = t462 * qJD(1);
t539 = cos(pkin(7));
t503 = t426 * t539;
t542 = cos(qJ(3));
t461 = t543 * t489 + t541 * t538;
t485 = t543 * t501;
t470 = t542 * t485;
t549 = -t543 * t491 + t541 * t536;
t551 = t549 * t539;
t557 = -t461 * t441 - t542 * t551 - t470;
t382 = qJD(3) * t557 + (qJD(1) * t483 - t503) * t441 + t427 * t542;
t552 = t461 * t542;
t410 = -t552 + (t551 + t485) * t441;
t502 = t438 * t539;
t420 = -t543 * t502 + t537 * t549;
t440 = sin(qJ(4));
t443 = cos(qJ(4));
t396 = t410 * t440 + t420 * t443;
t484 = t541 * t502;
t416 = qJD(1) * t484 + t426 * t537;
t352 = t396 * qJD(4) + t382 * t443 + t416 * t440;
t399 = t410 * t443 - t420 * t440;
t439 = sin(qJ(5));
t442 = cos(qJ(5));
t374 = t399 * t439 - t442 * t557;
t434 = t542 * t483;
t475 = t441 * t485;
t383 = qJD(1) * t434 - t427 * t441 - t542 * t503 + (t441 * t551 + t475 - t552) * qJD(3);
t347 = t374 * qJD(5) + t352 * t442 - t383 * t439;
t375 = t399 * t442 + t439 * t557;
t571 = t375 * qJD(5) - t352 * t439 - t383 * t442;
t351 = t399 * qJD(4) - t382 * t440 + t416 * t443;
t488 = t539 * t538;
t490 = t540 * t537;
t559 = (-t536 * t441 + t542 * t488) * t438 + t542 * t490;
t456 = t462 * t539;
t411 = t429 * t542 + (-t456 + t483) * t441;
t421 = t462 * t537 + t484;
t401 = t411 * t443 + t421 * t440;
t547 = -t429 * t441 - t542 * t456 + t434;
t376 = t401 * t439 + t442 * t547;
t556 = -0.2e1 * t376;
t425 = t461 * qJD(1);
t455 = qJD(1) * t551;
t381 = qJD(1) * t475 + t547 * qJD(3) - t425 * t542 + t441 * t455;
t400 = -t411 * t440 + t421 * t443;
t414 = t420 * qJD(1);
t350 = t400 * qJD(4) + t381 * t443 - t414 * t440;
t380 = -qJD(1) * t470 + t411 * qJD(3) - t425 * t441 - t542 * t455;
t345 = -t376 * qJD(5) + t350 * t442 + t380 * t439;
t377 = t401 * t442 - t439 * t547;
t369 = 0.1e1 / t377 ^ 2;
t555 = t345 * t369;
t418 = t441 * t490 + (t441 * t488 + t542 * t536) * t438;
t428 = -t538 * t501 + t540 * t539;
t406 = -t418 * t440 + t428 * t443;
t412 = t559 * qJD(3);
t389 = t406 * qJD(4) + t412 * t443;
t407 = t418 * t443 + t428 * t440;
t394 = t407 * t442 - t439 * t559;
t413 = t418 * qJD(3);
t363 = t394 * qJD(5) + t389 * t439 - t413 * t442;
t393 = t407 * t439 + t442 * t559;
t391 = 0.1e1 / t393 ^ 2;
t554 = t363 * t391;
t390 = 0.1e1 / t393;
t520 = t374 * t391;
t477 = -t390 * t396 - t406 * t520;
t553 = t439 * t477;
t359 = atan2(t374, t393);
t354 = sin(t359);
t355 = cos(t359);
t343 = t354 * t374 + t355 * t393;
t340 = 0.1e1 / t343;
t368 = 0.1e1 / t377;
t341 = 0.1e1 / t343 ^ 2;
t545 = 0.2e1 * t374;
t544 = 0.2e1 * t376;
t367 = t376 ^ 2;
t339 = t367 * t341 + 0.1e1;
t344 = t377 * qJD(5) + t350 * t439 - t380 * t442;
t529 = t341 * t376;
t366 = t374 ^ 2;
t358 = t366 * t391 + 0.1e1;
t356 = 0.1e1 / t358;
t480 = -t363 * t520 + t390 * t571;
t331 = t480 * t356;
t487 = -t354 * t393 + t355 * t374;
t326 = t487 * t331 + t354 * t571 + t355 * t363;
t342 = t340 * t341;
t534 = t326 * t342;
t535 = (t344 * t529 - t367 * t534) / t339 ^ 2;
t349 = -t401 * qJD(4) - t381 * t440 - t414 * t443;
t395 = t400 ^ 2;
t519 = t395 * t369;
t362 = 0.1e1 + t519;
t528 = t368 * t555;
t505 = t395 * t528;
t522 = t369 * t400;
t532 = (t349 * t522 - t505) / t362 ^ 2;
t523 = t390 * t554;
t531 = (-t366 * t523 + t520 * t571) / t358 ^ 2;
t530 = t341 * t344;
t527 = t354 * t376;
t526 = t355 * t376;
t521 = t374 * t390;
t518 = t400 * t439;
t517 = t559 * t443;
t515 = qJD(4) * t440;
t514 = qJD(5) * t442;
t513 = 0.2e1 * t535;
t512 = 0.2e1 * t532;
t511 = -0.2e1 * t531;
t510 = t342 * t544;
t509 = t390 * t531;
t508 = t341 * t527;
t507 = t341 * t526;
t506 = t400 * t528;
t499 = -0.2e1 * t340 * t535;
t498 = t341 * t513;
t497 = t326 * t510;
t496 = 0.2e1 * t506;
t495 = t523 * t545;
t494 = t522 * t532;
t482 = t440 * t547;
t481 = t443 * t547;
t479 = t375 * t390 - t394 * t520;
t452 = t443 * t557;
t385 = t410 * t442 + t439 * t452;
t402 = -t418 * t442 + t439 * t517;
t478 = -t385 * t390 - t402 * t520;
t476 = qJD(4) * t547;
t469 = -t354 + (-t355 * t521 + t354) * t356;
t467 = qJD(5) * t481 - t381;
t459 = t411 * qJD(5) - t380 * t443 - t440 * t476;
t388 = -t407 * qJD(4) - t412 * t440;
t387 = t411 * t439 + t442 * t481;
t386 = -t411 * t442 + t439 * t481;
t365 = (qJD(5) * t517 - t412) * t442 + (qJD(5) * t418 - t413 * t443 - t515 * t559) * t439;
t364 = -t393 * qJD(5) + t389 * t442 + t413 * t439;
t360 = 0.1e1 / t362;
t348 = (qJD(5) * t452 - t382) * t442 + (-qJD(5) * t410 + t383 * t443 - t515 * t557) * t439;
t337 = 0.1e1 / t339;
t336 = t356 * t553;
t335 = t478 * t356;
t334 = t479 * t356;
t330 = t469 * t376;
t329 = (-t354 * t396 + t355 * t406) * t439 + t487 * t336;
t327 = t487 * t334 + t354 * t375 + t355 * t394;
t325 = t478 * t511 + (t402 * t495 - t348 * t390 + (t363 * t385 - t365 * t374 - t402 * t571) * t391) * t356;
t323 = t479 * t511 + (t394 * t495 - t347 * t390 + (-t363 * t375 - t364 * t374 - t394 * t571) * t391) * t356;
t322 = t511 * t553 + (t477 * t514 + (t406 * t495 - t351 * t390 + (t363 * t396 - t374 * t388 - t406 * t571) * t391) * t439) * t356;
t1 = [t509 * t544 + (-t344 * t390 + t376 * t554) * t356, 0, t325, t322, t323, 0; t374 * t499 + (t571 * t340 + (-t374 * t326 - t330 * t344) * t341) * t337 + (t330 * t498 + (0.2e1 * t330 * t534 - (t331 * t356 * t521 + t511) * t508 - (t509 * t545 - t331 + (t331 - t480) * t356) * t507 - t469 * t530) * t337) * t376, 0, t386 * t499 + ((t459 * t439 + t467 * t442) * t340 - t386 * t341 * t326 - ((t325 * t374 + t335 * t571 + t365 + (-t335 * t393 - t385) * t331) * t355 + (-t325 * t393 - t335 * t363 - t348 + (-t335 * t374 - t402) * t331) * t354) * t529) * t337 + (t376 * t498 + (-t530 + t497) * t337) * (t487 * t335 - t354 * t385 + t355 * t402) (t329 * t529 - t340 * t518) * t513 + (-t329 * t530 + (t349 * t439 + t400 * t514) * t340 + (t329 * t510 - t341 * t518) * t326 - (t406 * t514 + t322 * t374 + t336 * t571 + t388 * t439 + (-t336 * t393 - t396 * t439) * t331) * t507 - (-t396 * t514 - t322 * t393 - t336 * t363 - t351 * t439 + (-t336 * t374 - t406 * t439) * t331) * t508) * t337 (t327 * t529 - t340 * t377) * t513 + (t327 * t497 + t345 * t340 + (-t377 * t326 - t327 * t344 - (t323 * t374 + t334 * t571 + t364 + (-t334 * t393 + t375) * t331) * t526 - (-t323 * t393 - t334 * t363 - t347 + (-t334 * t374 - t394) * t331) * t527) * t341) * t337, 0; (t368 * t396 + t375 * t522) * t512 + (-t351 * t368 + t375 * t496 + (t396 * t345 + t347 * t400 - t375 * t349) * t369) * t360, 0, t368 * t482 * t512 + 0.2e1 * t387 * t494 + ((t380 * t440 - t443 * t476) * t368 + (-t349 * t369 + t496) * t387 + t482 * t555 - (-t467 * t439 + t459 * t442) * t522) * t360 (t368 * t401 + t442 * t519) * t512 + (0.2e1 * t442 * t505 - t350 * t368 + (qJD(5) * t395 * t439 - 0.2e1 * t349 * t400 * t442 + t345 * t401) * t369) * t360, t494 * t556 + (t506 * t556 + (t344 * t400 + t349 * t376) * t369) * t360, 0;];
JaD_rot  = t1;
