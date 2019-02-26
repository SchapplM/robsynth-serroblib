% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRP6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:03
% EndTime: 2019-02-26 20:18:08
% DurationCPUTime: 4.78s
% Computational Cost: add. (30342->261), mult. (90603->485), div. (1041->12), fcn. (116353->17), ass. (0->184)
t532 = cos(pkin(12));
t533 = cos(pkin(6));
t486 = t533 * t532;
t531 = sin(pkin(12));
t534 = sin(qJ(2));
t536 = cos(qJ(2));
t468 = -t536 * t486 + t531 * t534;
t432 = t468 * qJD(2);
t469 = t486 * t534 + t531 * t536;
t433 = t469 * qJD(2);
t445 = cos(pkin(7));
t448 = sin(qJ(3));
t518 = t445 * t448;
t535 = cos(qJ(3));
t463 = t468 * t535;
t443 = sin(pkin(7));
t444 = sin(pkin(6));
t498 = t444 * t532;
t489 = t443 * t498;
t539 = -t445 * t463 - t448 * t469 - t489 * t535;
t388 = qJD(3) * t539 - t432 * t535 - t433 * t518;
t413 = t448 * (-t445 * t468 - t489) + t469 * t535;
t425 = t443 * t468 - t445 * t498;
t447 = sin(qJ(4));
t450 = cos(qJ(4));
t399 = -t413 * t447 + t425 * t450;
t520 = t443 * t447;
t364 = qJD(4) * t399 + t388 * t450 + t433 * t520;
t400 = t413 * t450 + t425 * t447;
t446 = sin(qJ(5));
t449 = cos(qJ(5));
t376 = t400 * t449 - t446 * t539;
t503 = t445 * t535;
t538 = t413 * qJD(3) - t432 * t448 + t433 * t503;
t345 = t376 * qJD(5) + t364 * t446 - t449 * t538;
t492 = t536 * t535;
t501 = t534 * t448;
t472 = -t445 * t501 + t492;
t473 = t445 * t492 - t501;
t499 = t443 * t533;
t481 = t535 * t499;
t409 = qJD(3) * t481 + (qJD(2) * t472 + qJD(3) * t473) * t444;
t491 = t534 * t535;
t502 = t536 * t448;
t474 = t445 * t502 + t491;
t490 = t448 * t499;
t424 = t444 * t474 + t490;
t505 = t443 * t536;
t436 = -t444 * t505 + t445 * t533;
t416 = -t424 * t447 + t436 * t450;
t504 = t443 * t534;
t493 = t444 * t504;
t482 = qJD(2) * t493;
t380 = qJD(4) * t416 + t409 * t450 + t447 * t482;
t417 = t424 * t450 + t436 * t447;
t423 = -t444 * t473 - t481;
t398 = t417 * t449 + t423 * t446;
t471 = t445 * t491 + t502;
t408 = qJD(3) * t490 + (qJD(2) * t471 + qJD(3) * t474) * t444;
t351 = qJD(5) * t398 + t380 * t446 - t408 * t449;
t374 = t400 * t446 + t449 * t539;
t369 = t374 ^ 2;
t397 = t417 * t446 - t423 * t449;
t392 = 0.1e1 / t397 ^ 2;
t359 = t369 * t392 + 0.1e1;
t355 = 0.1e1 / t359;
t391 = 0.1e1 / t397;
t522 = t374 * t392;
t331 = (-t345 * t391 + t351 * t522) * t355;
t361 = atan2(-t374, t397);
t353 = sin(t361);
t354 = cos(t361);
t484 = -t353 * t397 - t354 * t374;
t326 = t331 * t484 - t353 * t345 + t354 * t351;
t344 = -t353 * t374 + t354 * t397;
t341 = 0.1e1 / t344;
t342 = 0.1e1 / t344 ^ 2;
t542 = t326 * t341 * t342;
t485 = t533 * t531;
t437 = -t485 * t534 + t532 * t536;
t470 = t485 * t536 + t532 * t534;
t497 = t444 * t531;
t488 = t443 * t497;
t415 = t437 * t535 + (-t445 * t470 + t488) * t448;
t426 = t443 * t470 + t445 * t497;
t402 = t415 * t450 + t426 * t447;
t465 = t470 * t535;
t414 = t437 * t448 + t445 * t465 - t488 * t535;
t377 = t402 * t446 - t414 * t449;
t541 = -0.2e1 * t377;
t496 = 0.2e1 * t377 * t542;
t477 = -t391 * t399 + t416 * t522;
t540 = t446 * t477;
t527 = t351 * t391 * t392;
t537 = -0.2e1 * (t345 * t522 - t369 * t527) / t359 ^ 2;
t378 = t402 * t449 + t414 * t446;
t371 = 0.1e1 / t378;
t372 = 0.1e1 / t378 ^ 2;
t401 = -t415 * t447 + t426 * t450;
t396 = t401 ^ 2;
t524 = t372 * t396;
t360 = 0.1e1 + t524;
t434 = t470 * qJD(2);
t435 = t437 * qJD(2);
t390 = -qJD(3) * t414 - t434 * t535 - t435 * t518;
t519 = t443 * t450;
t365 = -qJD(4) * t402 - t390 * t447 + t435 * t519;
t366 = qJD(4) * t401 + t390 * t450 + t435 * t520;
t389 = qJD(3) * t415 - t434 * t448 + t435 * t503;
t348 = -qJD(5) * t377 + t366 * t449 + t389 * t446;
t528 = t348 * t371 * t372;
t507 = t396 * t528;
t523 = t372 * t401;
t530 = (t365 * t523 - t507) / t360 ^ 2;
t529 = t342 * t377;
t526 = t353 * t377;
t525 = t354 * t377;
t521 = t401 * t446;
t517 = t446 * t450;
t516 = t449 * t450;
t515 = qJD(4) * t447;
t514 = qJD(4) * t450;
t513 = qJD(5) * t446;
t512 = qJD(5) * t449;
t511 = qJD(5) * t450;
t370 = t377 ^ 2;
t340 = t342 * t370 + 0.1e1;
t347 = qJD(5) * t378 + t366 * t446 - t389 * t449;
t510 = 0.2e1 * (t347 * t529 - t370 * t542) / t340 ^ 2;
t509 = 0.2e1 * t530;
t506 = t401 * t528;
t500 = t446 * t515;
t495 = 0.2e1 * t506;
t494 = -0.2e1 * t374 * t527;
t421 = -t437 * t518 - t465;
t406 = t421 * t450 + t437 * t520;
t420 = t437 * t503 - t448 * t470;
t386 = t406 * t449 + t420 * t446;
t385 = t406 * t446 - t420 * t449;
t480 = -t376 * t391 + t398 * t522;
t459 = t450 * t539;
t381 = -t413 * t449 + t446 * t459;
t403 = -t423 * t517 - t424 * t449;
t479 = -t381 * t391 + t403 * t522;
t466 = t445 * t469;
t419 = -t448 * t466 - t463;
t467 = t443 * t469;
t404 = t419 * t450 + t447 * t467;
t458 = -t448 * t468 + t466 * t535;
t384 = t404 * t446 - t449 * t458;
t431 = t472 * t444;
t422 = t431 * t450 + t447 * t493;
t430 = t471 * t444;
t407 = t422 * t446 - t430 * t449;
t478 = -t384 * t391 + t407 * t522;
t405 = -t421 * t447 + t437 * t519;
t395 = -qJD(3) * t420 + t434 * t518 - t435 * t535;
t394 = qJD(3) * t421 - t434 * t503 - t435 * t448;
t383 = -t414 * t516 + t415 * t446;
t382 = -t414 * t517 - t415 * t449;
t379 = -qJD(4) * t417 - t409 * t447 + t450 * t482;
t368 = qJD(4) * t405 + t395 * t450 - t434 * t520;
t367 = -t431 * t500 + (t422 * t449 + t430 * t446) * qJD(5) + (t446 * t504 * t514 + (-t449 * t472 - t471 * t517) * qJD(3) + ((t447 * t505 - t450 * t474) * t446 - t473 * t449) * qJD(2)) * t444;
t363 = -qJD(4) * t400 - t388 * t447 + t433 * t519;
t362 = (-t423 * t511 - t409) * t449 + (qJD(5) * t424 - t408 * t450 + t423 * t515) * t446;
t357 = 0.1e1 / t360;
t352 = -qJD(5) * t397 + t380 * t449 + t408 * t446;
t350 = -t388 * t449 + t413 * t513 + t459 * t512 - t500 * t539 - t517 * t538;
t349 = ((-qJD(3) * t458 + t432 * t518 - t433 * t535) * t450 - t419 * t515 - t432 * t520 + t467 * t514) * t446 + t404 * t512 - (qJD(3) * t419 - t432 * t503 - t433 * t448) * t449 + t458 * t513;
t346 = -qJD(5) * t374 + t364 * t449 + t446 * t538;
t338 = 0.1e1 / t340;
t337 = t355 * t540;
t336 = t478 * t355;
t335 = t479 * t355;
t334 = t480 * t355;
t330 = (-t353 * t399 + t354 * t416) * t446 + t484 * t337;
t329 = t336 * t484 - t353 * t384 + t354 * t407;
t328 = t335 * t484 - t353 * t381 + t354 * t403;
t327 = t334 * t484 - t353 * t376 + t354 * t398;
t325 = t478 * t537 + (t407 * t494 - t349 * t391 + (t345 * t407 + t351 * t384 + t367 * t374) * t392) * t355;
t323 = t479 * t537 + (t403 * t494 - t350 * t391 + (t345 * t403 + t351 * t381 + t362 * t374) * t392) * t355;
t322 = t480 * t537 + (t398 * t494 - t346 * t391 + (t345 * t398 + t351 * t376 + t352 * t374) * t392) * t355;
t321 = t537 * t540 + (t477 * t512 + (t416 * t494 - t363 * t391 + (t345 * t416 + t351 * t399 + t374 * t379) * t392) * t446) * t355;
t1 = [0, t325, t323, t321, t322, 0; 0 (t329 * t529 - t341 * t385) * t510 + ((qJD(5) * t386 + t368 * t446 - t394 * t449) * t341 + t329 * t496 + (-t385 * t326 - t329 * t347 - (-t325 * t374 - t336 * t345 + t367 + (-t336 * t397 - t384) * t331) * t525 - (-t325 * t397 - t336 * t351 - t349 + (t336 * t374 - t407) * t331) * t526) * t342) * t338 (t328 * t529 - t341 * t382) * t510 + (t328 * t496 + (-t382 * t326 - t328 * t347 - (-t323 * t374 - t335 * t345 + t362 + (-t335 * t397 - t381) * t331) * t525 - (-t323 * t397 - t335 * t351 - t350 + (t335 * t374 - t403) * t331) * t526) * t342 + ((-t414 * t511 - t390) * t449 + (qJD(5) * t415 - t389 * t450 + t414 * t515) * t446) * t341) * t338 (t330 * t529 - t341 * t521) * t510 + ((t365 * t446 + t401 * t512) * t341 + t330 * t496 + (-t330 * t347 - t521 * t326 - (t416 * t512 - t321 * t374 - t337 * t345 + t379 * t446 + (-t337 * t397 - t399 * t446) * t331) * t525 - (-t399 * t512 - t321 * t397 - t337 * t351 - t363 * t446 + (t337 * t374 - t416 * t446) * t331) * t526) * t342) * t338 (t327 * t529 - t341 * t378) * t510 + (t327 * t496 + t348 * t341 + (-t378 * t326 - t327 * t347 - (-t322 * t374 - t334 * t345 + t352 + (-t334 * t397 - t376) * t331) * t525 - (-t322 * t397 - t334 * t351 - t346 + (t334 * t374 - t398) * t331) * t526) * t342) * t338, 0; 0 (-t371 * t405 + t386 * t523) * t509 + ((-qJD(4) * t406 - t395 * t447 - t434 * t519) * t371 + t386 * t495 + (-t405 * t348 - (-qJD(5) * t385 + t368 * t449 + t394 * t446) * t401 - t386 * t365) * t372) * t357 (-t371 * t414 * t447 + t383 * t523) * t509 + (t383 * t495 + (t389 * t447 + t414 * t514) * t371 + (-(-t389 * t516 + t390 * t446 + t415 * t512) * t401 - t383 * t365 + (-t447 * t348 - (t446 * t511 + t449 * t515) * t401) * t414) * t372) * t357 (t371 * t402 + t449 * t524) * t509 + (0.2e1 * t449 * t507 - t366 * t371 + (-0.2e1 * t365 * t401 * t449 + t348 * t402 + t396 * t513) * t372) * t357, t523 * t530 * t541 + (t506 * t541 + (t347 * t401 + t365 * t377) * t372) * t357, 0;];
JaD_rot  = t1;
