% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR9
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
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR9_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:25
% EndTime: 2019-02-26 20:53:29
% DurationCPUTime: 3.98s
% Computational Cost: add. (20977->187), mult. (61154->360), div. (726->12), fcn. (79375->19), ass. (0->159)
t435 = sin(qJ(3));
t527 = sin(pkin(7));
t528 = cos(pkin(13));
t476 = t528 * t527;
t429 = sin(pkin(13));
t488 = t429 * t527;
t532 = cos(qJ(3));
t410 = -t435 * t488 + t532 * t476;
t406 = t410 * qJD(3);
t529 = cos(pkin(7));
t477 = t529 * t528;
t489 = t429 * t529;
t412 = -t435 * t489 + t532 * t477;
t408 = t412 * qJD(3);
t432 = cos(pkin(6));
t437 = cos(qJ(1));
t526 = sin(pkin(12));
t486 = t437 * t526;
t431 = cos(pkin(12));
t530 = sin(qJ(1));
t492 = t530 * t431;
t456 = t432 * t492 + t486;
t415 = t456 * qJD(1);
t427 = t530 * t526;
t473 = t432 * t427;
t507 = qJD(1) * t437;
t416 = -qJD(1) * t473 + t431 * t507;
t423 = -t532 * t429 - t435 * t528;
t421 = t423 * qJD(3);
t422 = t435 * t429 - t532 * t528;
t430 = sin(pkin(6));
t450 = -t435 * t476 - t532 * t488;
t451 = -t435 * t477 - t532 * t489;
t457 = -t432 * t486 - t492;
t508 = t437 * t431;
t458 = -t432 * t508 + t427;
t368 = t458 * t408 - t415 * t451 + t416 * t422 + t457 * t421 - (-qJD(1) * t450 * t530 - t437 * t406) * t430;
t509 = t430 * t437;
t391 = -t422 * t457 - t450 * t509 - t451 * t458;
t487 = t430 * t529;
t404 = t437 * t487 - t458 * t527;
t434 = sin(qJ(5));
t531 = cos(qJ(5));
t379 = t391 * t531 + t404 * t434;
t474 = t530 * t487;
t453 = qJD(1) * t474 + t415 * t527;
t548 = t379 * qJD(5) + t368 * t434 + t453 * t531;
t544 = -t391 * t434 + t404 * t531;
t347 = qJD(5) * t544 + t368 * t531 - t453 * t434;
t373 = t544 ^ 2;
t398 = -t432 * t450 + (-t526 * t422 - t431 * t451) * t430;
t417 = -t430 * t431 * t527 + t432 * t529;
t464 = -t398 * t434 + t417 * t531;
t384 = 0.1e1 / t464 ^ 2;
t360 = t373 * t384 + 0.1e1;
t358 = 0.1e1 / t360;
t387 = t398 * t531 + t417 * t434;
t396 = t432 * t406 + (t408 * t431 + t526 * t421) * t430;
t371 = t387 * qJD(5) + t396 * t434;
t383 = 0.1e1 / t464;
t512 = t544 * t384;
t471 = t371 * t512 - t383 * t548;
t327 = t471 * t358;
t361 = atan2(-t544, -t464);
t356 = sin(t361);
t357 = cos(t361);
t475 = t356 * t464 - t357 * t544;
t322 = t475 * t327 + t356 * t548 + t357 * t371;
t339 = -t356 * t544 - t357 * t464;
t337 = 0.1e1 / t339 ^ 2;
t547 = t322 * t337;
t447 = t456 * t527 + t474;
t419 = -t473 + t508;
t493 = t430 * t530;
t448 = -t419 * t422 - t450 * t493 + t451 * t456;
t381 = t447 * t434 + t448 * t531;
t482 = t410 * t493;
t393 = -t456 * t412 + t419 * t423 + t482;
t433 = sin(qJ(6));
t436 = cos(qJ(6));
t354 = t381 * t433 + t393 * t436;
t543 = 0.2e1 * t354;
t336 = 0.1e1 / t339;
t542 = t336 * t547;
t414 = t457 * qJD(1);
t454 = qJD(1) * t458;
t491 = t430 * t507;
t444 = t406 * t493 - t456 * t408 - t414 * t422 + t419 * t421 - t450 * t491 - t451 * t454;
t446 = qJD(1) * t404;
t343 = t381 * qJD(5) + t434 * t444 - t531 * t446;
t445 = t447 * t531;
t380 = t434 * t448 - t445;
t533 = 0.2e1 * t380;
t485 = t533 * t542;
t539 = -t337 * t343 + t485;
t374 = t380 ^ 2;
t335 = t337 * t374 + 0.1e1;
t520 = t337 * t380;
t504 = 0.2e1 * (t343 * t520 - t374 * t542) / t335 ^ 2;
t538 = t371 * t384;
t397 = t432 * t410 + (t412 * t431 + t526 * t423) * t430;
t535 = -t410 * t509 - t412 * t458 - t423 * t457;
t468 = t383 * t535 + t397 * t512;
t537 = t434 * t468;
t490 = qJD(5) * t531;
t407 = t450 * qJD(3);
t409 = t451 * qJD(3);
t420 = t422 * qJD(3);
t536 = qJD(1) * t482 - t407 * t509 - t409 * t458 - t415 * t412 + t416 * t423 - t420 * t457;
t355 = t381 * t436 - t393 * t433;
t349 = 0.1e1 / t355;
t350 = 0.1e1 / t355 ^ 2;
t534 = -0.2e1 * t544;
t514 = t383 * t538;
t523 = (t373 * t514 - t512 * t548) / t360 ^ 2;
t506 = qJD(5) * t434;
t344 = qJD(5) * t445 + t434 * t446 + t444 * t531 - t448 * t506;
t363 = t407 * t493 - t456 * t409 + t410 * t491 + t412 * t454 + t414 * t423 + t419 * t420;
t505 = qJD(6) * t354;
t332 = t344 * t436 - t363 * t433 - t505;
t522 = t332 * t349 * t350;
t348 = t354 ^ 2;
t342 = t348 * t350 + 0.1e1;
t340 = 0.1e1 / t342;
t519 = t340 * t350;
t331 = t355 * qJD(6) + t344 * t433 + t363 * t436;
t517 = t350 * t354;
t518 = 0.1e1 / t342 ^ 2 * (t331 * t517 - t348 * t522);
t516 = t356 * t380;
t515 = t357 * t380;
t513 = t544 * t383;
t510 = t393 * t434;
t503 = -0.2e1 * t523;
t501 = -0.2e1 * t518;
t500 = t350 * t518;
t499 = t383 * t523;
t498 = t331 * t519;
t497 = t354 * t522;
t494 = t393 * t531;
t484 = 0.2e1 * t497;
t483 = t514 * t534;
t353 = t379 * t436 + t433 * t535;
t352 = t379 * t433 - t436 * t535;
t472 = qJD(6) * t494 - t444;
t470 = -t349 * t433 + t436 * t517;
t469 = -t379 * t383 + t387 * t512;
t460 = -t356 + (-t357 * t513 + t356) * t358;
t455 = qJD(6) * t448 + t531 * t363 - t393 * t506;
t395 = t432 * t407 + (t409 * t431 + t526 * t420) * t430;
t372 = t464 * qJD(5) + t396 * t531;
t370 = t433 * t448 + t436 * t494;
t333 = 0.1e1 / t335;
t330 = t358 * t537;
t329 = t469 * t358;
t324 = (-t356 * t535 + t357 * t397) * t434 + t475 * t330;
t323 = t475 * t329 + t356 * t379 + t357 * t387;
t320 = t469 * t503 + (-t387 * t483 - t347 * t383 + (-t371 * t379 + t372 * t544 - t387 * t548) * t384) * t358;
t319 = t503 * t537 + ((-t397 * t483 + t536 * t383 + (t371 * t535 + t395 * t544 - t397 * t548) * t384) * t434 + t468 * t490) * t358;
t1 = [-t499 * t533 + (t343 * t383 + t380 * t538) * t358, 0, t319, 0, t320, 0; t544 * t336 * t504 + (t548 * t336 + t544 * t547 - (t460 * t343 + ((t327 * t358 * t513 + t503) * t356 + (-t499 * t534 - t327 + (t327 - t471) * t358) * t357) * t380) * t520) * t333 + (t539 * t333 + t520 * t504) * t460 * t380, 0 (t324 * t520 - t336 * t510) * t504 + ((t363 * t434 + t393 * t490) * t336 + t539 * t324 + (-t510 * t322 - (t397 * t490 - t319 * t544 + t330 * t548 + t395 * t434 + (t330 * t464 - t434 * t535) * t327) * t515 - (-t535 * t490 + t319 * t464 - t330 * t371 - t536 * t434 + (t330 * t544 - t397 * t434) * t327) * t516) * t337) * t333, 0 (t323 * t520 - t336 * t381) * t504 + (t323 * t485 + t344 * t336 + (-t381 * t322 - t323 * t343 - (-t320 * t544 + t329 * t548 + t372 + (t329 * t464 + t379) * t327) * t515 - (t320 * t464 - t329 * t371 + t347 + (t329 * t544 - t387) * t327) * t516) * t337) * t333, 0; 0.2e1 * (-t349 * t352 + t353 * t517) * t518 + ((t353 * qJD(6) + t347 * t433 - t436 * t536) * t349 + t353 * t484 + (-t352 * t332 - (-t352 * qJD(6) + t347 * t436 + t433 * t536) * t354 - t353 * t331) * t350) * t340, 0 (t500 * t543 - t498) * t370 + (-t332 * t519 + t349 * t501) * (t433 * t494 - t436 * t448) + ((t455 * t433 + t472 * t436) * t349 - (-t472 * t433 + t455 * t436) * t517 + t370 * t484) * t340, 0, t470 * t380 * t501 + (t470 * t343 + ((-qJD(6) * t349 - 0.2e1 * t497) * t436 + (t331 * t436 + (t332 - t505) * t433) * t350) * t380) * t340, t501 + (t498 + (-t340 * t522 - t500) * t354) * t543;];
JaD_rot  = t1;
