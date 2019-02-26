% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR10_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:51
% EndTime: 2019-02-26 21:19:55
% DurationCPUTime: 4.45s
% Computational Cost: add. (24403->171), mult. (55387->326), div. (989->12), fcn. (71220->17), ass. (0->154)
t526 = sin(pkin(13));
t531 = cos(pkin(6));
t484 = t531 * t526;
t529 = cos(pkin(13));
t532 = sin(qJ(1));
t533 = cos(qJ(1));
t424 = -t532 * t484 + t533 * t529;
t421 = t424 * qJD(1);
t434 = sin(qJ(3));
t436 = cos(qJ(3));
t423 = t533 * t484 + t532 * t529;
t527 = sin(pkin(7));
t528 = sin(pkin(6));
t482 = t528 * t527;
t472 = t533 * t482;
t530 = cos(pkin(7));
t485 = t531 * t529;
t540 = -t533 * t485 + t532 * t526;
t541 = t540 * t530;
t447 = t541 + t472;
t536 = t423 * t434 + t447 * t436;
t456 = t532 * t485 + t533 * t526;
t420 = t456 * qJD(1);
t469 = t532 * t482;
t546 = qJD(1) * t469 - t420 * t530;
t383 = t536 * qJD(3) - t421 * t436 - t434 * t546;
t507 = t423 * t436;
t404 = t447 * t434 - t507;
t505 = qJ(4) + qJ(5);
t431 = sin(t505);
t432 = qJD(4) + qJD(5);
t483 = t530 * t528;
t448 = -t533 * t483 + t527 * t540;
t470 = t532 * t483;
t457 = qJD(1) * t470 + t420 * t527;
t493 = cos(t505);
t486 = t432 * t493;
t559 = t404 * t486 + (-t448 * t432 + t383) * t431 + t457 * t493;
t557 = -t383 * t493 + (t404 * t432 + t457) * t431;
t445 = t448 * t493;
t391 = t404 * t431 + t445;
t392 = t404 * t493 - t448 * t431;
t386 = t391 ^ 2;
t453 = t529 * t483 + t531 * t527;
t481 = t528 * t526;
t415 = t453 * t434 + t436 * t481;
t422 = -t529 * t482 + t531 * t530;
t466 = -t415 * t431 + t422 * t493;
t396 = 0.1e1 / t466 ^ 2;
t374 = t386 * t396 + 0.1e1;
t368 = 0.1e1 / t374;
t414 = -t434 * t481 + t453 * t436;
t409 = t414 * qJD(3);
t384 = t415 * t486 + (t422 * t432 + t409) * t431;
t395 = 0.1e1 / t466;
t512 = t391 * t396;
t478 = -t384 * t512 - t395 * t559;
t341 = t478 * t368;
t375 = atan2(t391, -t466);
t362 = sin(t375);
t363 = cos(t375);
t480 = t362 * t466 + t363 * t391;
t336 = t480 * t341 + t362 * t559 + t363 * t384;
t353 = t362 * t391 - t363 * t466;
t351 = 0.1e1 / t353 ^ 2;
t553 = t336 * t351;
t550 = -t456 * t530 + t469;
t406 = t424 * t436 + t550 * t434;
t446 = t456 * t527 + t470;
t394 = t406 * t493 + t446 * t431;
t538 = t550 * t436;
t405 = t424 * t434 - t538;
t433 = sin(qJ(6));
t435 = cos(qJ(6));
t372 = t394 * t433 - t405 * t435;
t549 = 0.2e1 * t372;
t350 = 0.1e1 / t353;
t548 = t350 * t553;
t443 = t446 * t493;
t393 = t406 * t431 - t443;
t534 = 0.2e1 * t393;
t491 = t534 * t548;
t419 = t423 * qJD(1);
t504 = qJD(3) * t434;
t539 = t447 * qJD(1);
t379 = t538 * qJD(3) - t419 * t436 - t424 * t504 + t539 * t434;
t444 = qJD(1) * t448;
t357 = t406 * t486 + t493 * t444 + (t446 * t432 + t379) * t431;
t518 = t357 * t351;
t545 = -t518 + t491;
t544 = (t434 * t541 - t507) * qJD(3) - t421 * t434 + t546 * t436 + t472 * t504;
t543 = t384 * t396;
t475 = -t395 * t536 - t414 * t512;
t542 = t431 * t475;
t373 = t394 * t435 + t405 * t433;
t365 = 0.1e1 / t373;
t366 = 0.1e1 / t373 ^ 2;
t535 = 0.2e1 * t391;
t387 = t393 ^ 2;
t347 = t351 * t387 + 0.1e1;
t525 = (-t387 * t548 + t393 * t518) / t347 ^ 2;
t506 = t431 * t432;
t358 = t379 * t493 - t406 * t506 - t431 * t444 + t432 * t443;
t378 = t406 * qJD(3) - t419 * t434 - t539 * t436;
t348 = t373 * qJD(6) + t358 * t433 - t378 * t435;
t364 = t372 ^ 2;
t356 = t364 * t366 + 0.1e1;
t515 = t366 * t372;
t503 = qJD(6) * t372;
t349 = t358 * t435 + t378 * t433 - t503;
t521 = t349 * t365 * t366;
t524 = (t348 * t515 - t364 * t521) / t356 ^ 2;
t514 = t395 * t543;
t522 = (t386 * t514 + t512 * t559) / t374 ^ 2;
t520 = t351 * t393;
t354 = 0.1e1 / t356;
t519 = t354 * t366;
t517 = t362 * t393;
t516 = t363 * t393;
t513 = t391 * t395;
t510 = t405 * t431;
t502 = 0.2e1 * t525;
t501 = -0.2e1 * t524;
t500 = -0.2e1 * t522;
t498 = t366 * t524;
t497 = t395 * t522;
t496 = t348 * t519;
t495 = t372 * t521;
t490 = 0.2e1 * t495;
t489 = t514 * t535;
t487 = t405 * t493;
t371 = t392 * t435 - t433 * t536;
t370 = t392 * t433 + t435 * t536;
t477 = -t433 * t365 + t435 * t515;
t399 = t415 * t493 + t422 * t431;
t476 = -t392 * t395 - t399 * t512;
t468 = qJD(6) * t487 + t379;
t467 = -t362 + (t363 * t513 + t362) * t368;
t458 = qJD(6) * t406 - t493 * t378 + t405 * t506;
t410 = t415 * qJD(3);
t385 = t409 * t493 + t466 * t432;
t377 = t406 * t433 - t435 * t487;
t361 = -t448 * t486 - t557;
t360 = t432 * t445 + t557;
t345 = 0.1e1 / t347;
t344 = t368 * t542;
t343 = t476 * t368;
t338 = (t362 * t536 + t363 * t414) * t431 + t480 * t344;
t337 = t480 * t343 + t362 * t392 + t363 * t399;
t334 = t476 * t500 + (-t399 * t489 + t360 * t395 + (-t384 * t392 - t385 * t391 - t399 * t559) * t396) * t368;
t333 = t500 * t542 + ((-t414 * t489 + t544 * t395 + (-t384 * t536 + t391 * t410 - t414 * t559) * t396) * t431 + t475 * t486) * t368;
t332 = t477 * t393 * t501 + (t477 * t357 + ((-qJD(6) * t365 - 0.2e1 * t495) * t435 + (t348 * t435 + (t349 - t503) * t433) * t366) * t393) * t354;
t331 = (t337 * t520 - t350 * t394) * t502 + (t337 * t491 + t358 * t350 + (-t394 * t336 - t337 * t357 - (t334 * t391 + t343 * t559 + t385 + (t343 * t466 + t392) * t341) * t516 - (t334 * t466 - t343 * t384 - t360 + (-t343 * t391 - t399) * t341) * t517) * t351) * t345;
t1 = [-t497 * t534 + (t357 * t395 + t393 * t543) * t368, 0, t333, t334, t334, 0; -0.2e1 * t391 * t350 * t525 + (t559 * t350 - t391 * t553 - (t467 * t357 + ((-t341 * t368 * t513 + t500) * t362 + (-t497 * t535 - t341 + (t341 - t478) * t368) * t363) * t393) * t520) * t345 + (t345 * t545 + t520 * t502) * t467 * t393, 0 (t338 * t520 + t350 * t510) * t502 + ((-t378 * t431 - t405 * t486) * t350 + t545 * t338 + (t510 * t336 - (t414 * t486 + t333 * t391 + t344 * t559 - t410 * t431 + (t344 * t466 + t431 * t536) * t341) * t516 - (t536 * t486 + t333 * t466 - t344 * t384 - t544 * t431 + (-t344 * t391 - t414 * t431) * t341) * t517) * t351) * t345, t331, t331, 0; 0.2e1 * (-t365 * t370 + t371 * t515) * t524 + ((t371 * qJD(6) + t361 * t433 - t435 * t544) * t365 + t371 * t490 + (-t370 * t349 - (-t370 * qJD(6) + t361 * t435 + t433 * t544) * t372 - t371 * t348) * t366) * t354, 0 (t498 * t549 - t496) * t377 + (-t349 * t519 + t365 * t501) * (-t406 * t435 - t433 * t487) + ((t458 * t433 - t468 * t435) * t365 - (t468 * t433 + t458 * t435) * t515 + t377 * t490) * t354, t332, t332, t501 + (t496 + (-t354 * t521 - t498) * t372) * t549;];
JaD_rot  = t1;
