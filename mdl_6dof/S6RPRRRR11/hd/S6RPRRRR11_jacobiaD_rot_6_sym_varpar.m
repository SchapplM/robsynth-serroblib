% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR11
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

function JaD_rot = S6RPRRRR11_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:20:27
% EndTime: 2019-02-26 21:20:30
% DurationCPUTime: 3.73s
% Computational Cost: add. (15019->167), mult. (43174->321), div. (744->12), fcn. (55416->17), ass. (0->152)
t511 = sin(pkin(13));
t516 = cos(pkin(6));
t468 = t516 * t511;
t514 = cos(pkin(13));
t517 = sin(qJ(1));
t519 = cos(qJ(1));
t407 = -t517 * t468 + t519 * t514;
t404 = t407 * qJD(1);
t419 = sin(qJ(3));
t420 = cos(qJ(3));
t406 = t519 * t468 + t517 * t514;
t512 = sin(pkin(7));
t513 = sin(pkin(6));
t466 = t513 * t512;
t453 = t519 * t466;
t515 = cos(pkin(7));
t469 = t516 * t514;
t526 = -t519 * t469 + t517 * t511;
t527 = t526 * t515;
t430 = t527 + t453;
t522 = t406 * t419 + t430 * t420;
t439 = t517 * t469 + t519 * t511;
t403 = t439 * qJD(1);
t450 = t517 * t466;
t531 = qJD(1) * t450 - t403 * t515;
t366 = qJD(3) * t522 - t404 * t420 - t531 * t419;
t418 = sin(qJ(4));
t467 = t515 * t513;
t451 = t517 * t467;
t440 = qJD(1) * t451 + t403 * t512;
t518 = cos(qJ(4));
t492 = t406 * t420;
t387 = t430 * t419 - t492;
t431 = -t519 * t467 + t512 * t526;
t540 = -t387 * t418 - t431 * t518;
t343 = -qJD(4) * t540 - t366 * t518 + t440 * t418;
t375 = t387 * t518 - t431 * t418;
t542 = t375 * qJD(4) + t366 * t418 + t440 * t518;
t369 = t540 ^ 2;
t436 = t514 * t467 + t516 * t512;
t465 = t513 * t511;
t398 = t436 * t419 + t420 * t465;
t405 = -t514 * t466 + t516 * t515;
t456 = -t398 * t418 + t405 * t518;
t380 = 0.1e1 / t456 ^ 2;
t357 = t369 * t380 + 0.1e1;
t355 = 0.1e1 / t357;
t383 = t398 * t518 + t405 * t418;
t397 = -t419 * t465 + t436 * t420;
t392 = t397 * qJD(3);
t367 = t383 * qJD(4) + t392 * t418;
t379 = 0.1e1 / t456;
t497 = t540 * t380;
t461 = t367 * t497 - t379 * t542;
t324 = t461 * t355;
t358 = atan2(-t540, -t456);
t353 = sin(t358);
t354 = cos(t358);
t464 = t353 * t456 - t354 * t540;
t319 = t464 * t324 + t353 * t542 + t354 * t367;
t336 = -t353 * t540 - t354 * t456;
t334 = 0.1e1 / t336 ^ 2;
t539 = t334 * t319;
t536 = -t439 * t515 + t450;
t389 = t407 * t420 + t536 * t419;
t429 = t439 * t512 + t451;
t377 = t389 * t518 + t429 * t418;
t524 = t536 * t420;
t388 = t407 * t419 - t524;
t417 = qJ(5) + qJ(6);
t414 = sin(t417);
t415 = cos(t417);
t351 = t377 * t414 - t388 * t415;
t535 = 0.2e1 * t351;
t333 = 0.1e1 / t336;
t534 = t333 * t539;
t427 = t429 * t518;
t376 = t389 * t418 - t427;
t520 = 0.2e1 * t376;
t477 = t520 * t534;
t402 = t406 * qJD(1);
t491 = qJD(3) * t419;
t525 = t430 * qJD(1);
t362 = qJD(3) * t524 - t402 * t420 - t407 * t491 + t419 * t525;
t428 = qJD(1) * t431;
t340 = t377 * qJD(4) + t362 * t418 + t518 * t428;
t503 = t340 * t334;
t530 = -t503 + t477;
t363 = (t419 * t527 - t492) * qJD(3) - t404 * t419 + t531 * t420 + t453 * t491;
t370 = t376 ^ 2;
t332 = t334 * t370 + 0.1e1;
t489 = 0.2e1 * (-t370 * t534 + t376 * t503) / t332 ^ 2;
t529 = t367 * t380;
t458 = -t379 * t522 + t397 * t497;
t528 = t418 * t458;
t479 = qJD(4) * t518;
t352 = t377 * t415 + t388 * t414;
t346 = 0.1e1 / t352;
t347 = 0.1e1 / t352 ^ 2;
t521 = -0.2e1 * t540;
t361 = t389 * qJD(3) - t402 * t419 - t420 * t525;
t416 = qJD(5) + qJD(6);
t472 = t377 * t416 - t361;
t490 = qJD(4) * t418;
t341 = qJD(4) * t427 + t362 * t518 - t389 * t490 - t418 * t428;
t474 = t388 * t416 + t341;
t328 = t474 * t414 + t472 * t415;
t345 = t351 ^ 2;
t339 = t345 * t347 + 0.1e1;
t502 = t347 * t351;
t329 = -t472 * t414 + t474 * t415;
t506 = t329 * t346 * t347;
t509 = (t328 * t502 - t345 * t506) / t339 ^ 2;
t499 = t379 * t529;
t507 = (t369 * t499 - t497 * t542) / t357 ^ 2;
t505 = t334 * t376;
t337 = 0.1e1 / t339;
t504 = t337 * t347;
t501 = t353 * t376;
t500 = t354 * t376;
t498 = t540 * t379;
t495 = t388 * t418;
t488 = -0.2e1 * t509;
t487 = -0.2e1 * t507;
t485 = t347 * t509;
t484 = t379 * t507;
t483 = t328 * t504;
t482 = t351 * t506;
t480 = t388 * t518;
t476 = 0.2e1 * t482;
t475 = t499 * t521;
t473 = -t416 * t522 - t343;
t471 = t375 * t416 - t363;
t463 = t416 * t480 + t362;
t460 = -t346 * t414 + t415 * t502;
t459 = -t375 * t379 + t383 * t497;
t449 = -t353 + (-t354 * t498 + t353) * t355;
t443 = -t518 * t361 + t388 * t490 + t389 * t416;
t393 = t398 * qJD(3);
t368 = t456 * qJD(4) + t392 * t518;
t360 = t389 * t414 - t415 * t480;
t350 = t375 * t415 - t414 * t522;
t349 = t375 * t414 + t415 * t522;
t330 = 0.1e1 / t332;
t327 = t355 * t528;
t325 = t459 * t355;
t321 = (t353 * t522 + t354 * t397) * t418 + t464 * t327;
t320 = t464 * t325 + t353 * t375 + t354 * t383;
t317 = t459 * t487 + (-t383 * t475 + t343 * t379 + (-t367 * t375 + t368 * t540 - t383 * t542) * t380) * t355;
t316 = t487 * t528 + ((-t397 * t475 + t363 * t379 + (-t367 * t522 - t393 * t540 - t397 * t542) * t380) * t418 + t458 * t479) * t355;
t315 = t488 + (t483 + (-t337 * t506 - t485) * t351) * t535;
t1 = [-t484 * t520 + (t340 * t379 + t376 * t529) * t355, 0, t316, t317, 0, 0; t540 * t333 * t489 + (t542 * t333 + t540 * t539 - (t449 * t340 + ((t324 * t355 * t498 + t487) * t353 + (-t484 * t521 - t324 + (t324 - t461) * t355) * t354) * t376) * t505) * t330 + (t530 * t330 + t505 * t489) * t449 * t376, 0 (t321 * t505 + t333 * t495) * t489 + ((-t361 * t418 - t388 * t479) * t333 + t530 * t321 + (t495 * t319 - (t397 * t479 - t316 * t540 + t327 * t542 - t393 * t418 + (t327 * t456 + t418 * t522) * t324) * t500 - (t522 * t479 + t316 * t456 - t327 * t367 - t363 * t418 + (t327 * t540 - t397 * t418) * t324) * t501) * t334) * t330 (t320 * t505 - t333 * t377) * t489 + (t320 * t477 + t341 * t333 + (-t377 * t319 - t320 * t340 - (-t317 * t540 + t325 * t542 + t368 + (t325 * t456 + t375) * t324) * t500 - (t317 * t456 - t325 * t367 - t343 + (t325 * t540 - t383) * t324) * t501) * t334) * t330, 0, 0; 0.2e1 * (-t346 * t349 + t350 * t502) * t509 + ((t473 * t414 + t471 * t415) * t346 + t350 * t476 + (-t349 * t329 - (-t471 * t414 + t473 * t415) * t351 - t350 * t328) * t347) * t337, 0 (t485 * t535 - t483) * t360 + (-t329 * t504 + t346 * t488) * (-t389 * t415 - t414 * t480) + ((t443 * t414 - t463 * t415) * t346 - (t463 * t414 + t443 * t415) * t502 + t360 * t476) * t337, t460 * t376 * t488 + (t460 * t340 + ((-t346 * t416 - 0.2e1 * t482) * t415 + (t328 * t415 + (-t351 * t416 + t329) * t414) * t347) * t376) * t337, t315, t315;];
JaD_rot  = t1;
