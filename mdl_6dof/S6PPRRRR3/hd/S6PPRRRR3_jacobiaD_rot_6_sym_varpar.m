% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRR3
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRRR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobiaD_rot_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:43:57
% EndTime: 2019-02-26 19:44:02
% DurationCPUTime: 5.34s
% Computational Cost: add. (40633->214), mult. (120610->400), div. (816->12), fcn. (158008->21), ass. (0->169)
t573 = sin(pkin(14));
t574 = sin(pkin(13));
t533 = t574 * t573;
t577 = cos(pkin(14));
t578 = cos(pkin(13));
t539 = t578 * t577;
t580 = cos(pkin(6));
t513 = -t539 * t580 + t533;
t575 = sin(pkin(7));
t576 = sin(pkin(6));
t537 = t576 * t575;
t579 = cos(pkin(7));
t589 = t513 * t579 + t578 * t537;
t535 = t574 * t577;
t538 = t578 * t573;
t515 = t535 * t580 + t538;
t534 = t574 * t576;
t588 = t515 * t579 - t575 * t534;
t540 = t579 * t576;
t587 = t577 * t540 + t580 * t575;
t482 = sin(qJ(4));
t486 = cos(qJ(4));
t478 = sin(pkin(8));
t479 = cos(pkin(8));
t483 = sin(qJ(3));
t514 = t538 * t580 + t535;
t581 = cos(qJ(3));
t498 = t483 * t514 + t581 * t589;
t502 = t513 * t575 - t540 * t578;
t493 = t478 * t502 - t479 * t498;
t585 = t483 * t589 - t514 * t581;
t442 = t482 * t585 + t486 * t493;
t459 = t498 * qJD(3);
t460 = t585 * qJD(3);
t555 = t479 * t482;
t421 = qJD(4) * t442 - t459 * t486 + t460 * t555;
t443 = t482 * t493 - t486 * t585;
t481 = sin(qJ(5));
t485 = cos(qJ(5));
t494 = t478 * t498 + t479 * t502;
t428 = t443 * t485 + t481 * t494;
t556 = t478 * t485;
t395 = qJD(5) * t428 + t421 * t481 + t460 * t556;
t426 = t443 * t481 - t485 * t494;
t424 = t426 ^ 2;
t536 = t576 * t573;
t470 = t483 * t587 + t581 * t536;
t469 = -t483 * t536 + t581 * t587;
t473 = -t537 * t577 + t579 * t580;
t531 = t469 * t479 + t473 * t478;
t455 = t470 * t486 + t482 * t531;
t465 = -t469 * t478 + t473 * t479;
t446 = t455 * t481 - t465 * t485;
t440 = 0.1e1 / t446 ^ 2;
t411 = t424 * t440 + 0.1e1;
t409 = 0.1e1 / t411;
t454 = -t470 * t482 + t486 * t531;
t467 = t469 * qJD(3);
t468 = t470 * qJD(3);
t437 = qJD(4) * t454 + t467 * t486 - t468 * t555;
t447 = t455 * t485 + t465 * t481;
t413 = qJD(5) * t447 + t437 * t481 - t468 * t556;
t439 = 0.1e1 / t446;
t563 = t426 * t440;
t378 = (-t395 * t439 + t413 * t563) * t409;
t412 = atan2(-t426, t446);
t407 = sin(t412);
t408 = cos(t412);
t532 = -t407 * t446 - t408 * t426;
t373 = t378 * t532 - t407 * t395 + t408 * t413;
t391 = -t407 * t426 + t408 * t446;
t388 = 0.1e1 / t391;
t389 = 0.1e1 / t391 ^ 2;
t586 = t373 * t388 * t389;
t516 = -t533 * t580 + t539;
t464 = -t483 * t588 + t516 * t581;
t499 = t516 * t483 + t581 * t588;
t503 = t515 * t575 + t534 * t579;
t500 = t503 * t478;
t495 = -t479 * t499 + t500;
t445 = t464 * t486 + t482 * t495;
t496 = t478 * t499 + t479 * t503;
t429 = t445 * t481 - t485 * t496;
t545 = 0.2e1 * t429 * t586;
t526 = -t439 * t442 + t454 * t563;
t584 = t481 * t526;
t564 = t413 * t439 * t440;
t582 = -0.2e1 * (t395 * t563 - t424 * t564) / t411 ^ 2;
t430 = t445 * t485 + t481 * t496;
t497 = t499 * t486;
t559 = t464 * t482;
t444 = t479 * t497 - t486 * t500 + t559;
t480 = sin(qJ(6));
t484 = cos(qJ(6));
t406 = t430 * t484 + t444 * t480;
t402 = 0.1e1 / t406;
t403 = 0.1e1 / t406 ^ 2;
t461 = t499 * qJD(3);
t462 = t464 * qJD(3);
t423 = -t462 * t555 - t461 * t486 + (t486 * t495 - t559) * qJD(4);
t557 = t478 * t481;
t398 = -qJD(5) * t429 + t423 * t485 + t462 * t557;
t554 = t479 * t486;
t422 = qJD(4) * t445 - t461 * t482 + t462 * t554;
t405 = t430 * t480 - t444 * t484;
t551 = qJD(6) * t405;
t387 = t398 * t484 + t422 * t480 - t551;
t572 = t387 * t402 * t403;
t571 = t389 * t429;
t386 = qJD(6) * t406 + t398 * t480 - t422 * t484;
t401 = t405 ^ 2;
t394 = t401 * t403 + 0.1e1;
t568 = t403 * t405;
t570 = 0.1e1 / t394 ^ 2 * (t386 * t568 - t401 * t572);
t569 = t402 * t480;
t567 = t405 * t484;
t566 = t407 * t429;
t565 = t408 * t429;
t562 = t444 * t481;
t561 = t444 * t485;
t553 = qJD(5) * t481;
t552 = qJD(5) * t485;
t425 = t429 ^ 2;
t385 = t389 * t425 + 0.1e1;
t397 = qJD(5) * t430 + t423 * t481 - t462 * t556;
t550 = 0.2e1 * (t397 * t571 - t425 * t586) / t385 ^ 2;
t548 = -0.2e1 * t570;
t547 = 0.2e1 * t570;
t546 = t405 * t572;
t544 = 0.2e1 * t546;
t543 = -0.2e1 * t426 * t564;
t542 = qJD(6) * t561 + t423;
t450 = -t464 * t555 - t497;
t435 = t450 * t485 + t464 * t557;
t449 = t464 * t554 - t482 * t499;
t418 = t435 * t484 + t449 * t480;
t417 = t435 * t480 - t449 * t484;
t529 = t403 * t567 - t569;
t528 = -t428 * t439 + t447 * t563;
t448 = -t486 * t498 + t555 * t585;
t433 = t448 * t481 + t556 * t585;
t456 = t469 * t486 - t470 * t555;
t451 = t456 * t481 - t470 * t556;
t527 = -t433 * t439 + t451 * t563;
t525 = -t450 * t481 + t464 * t556;
t519 = qJD(6) * t445 - t422 * t485 + t444 * t553;
t436 = -t455 * qJD(4) - t467 * t482 - t468 * t554;
t432 = -qJD(4) * t449 + t461 * t555 - t462 * t486;
t431 = qJD(4) * t450 - t461 * t554 - t462 * t482;
t420 = -t443 * qJD(4) + t459 * t482 + t460 * t554;
t419 = (-t467 * t555 - t468 * t486 + (-t469 * t482 - t470 * t554) * qJD(4)) * t481 - t467 * t556 + (t456 * t485 + t470 * t557) * qJD(5);
t416 = t445 * t480 - t484 * t561;
t415 = -t445 * t484 - t480 * t561;
t414 = -qJD(5) * t446 + t437 * t485 + t468 * t557;
t400 = qJD(5) * t525 + t432 * t485 - t461 * t557;
t399 = (t459 * t555 + t460 * t486 + (t482 * t498 + t554 * t585) * qJD(4)) * t481 + t448 * t552 + t459 * t556 - t585 * t478 * t553;
t396 = -qJD(5) * t426 + t421 * t485 - t460 * t557;
t392 = 0.1e1 / t394;
t383 = 0.1e1 / t385;
t382 = t409 * t584;
t381 = t527 * t409;
t380 = t528 * t409;
t376 = (-t407 * t442 + t408 * t454) * t481 + t532 * t382;
t375 = t381 * t532 - t407 * t433 + t408 * t451;
t374 = t380 * t532 - t407 * t428 + t408 * t447;
t372 = t527 * t582 + (t451 * t543 - t399 * t439 + (t395 * t451 + t413 * t433 + t419 * t426) * t440) * t409;
t370 = t528 * t582 + (t447 * t543 - t396 * t439 + (t395 * t447 + t413 * t428 + t414 * t426) * t440) * t409;
t369 = t582 * t584 + (t526 * t552 + (t454 * t543 - t420 * t439 + (t395 * t454 + t413 * t442 + t426 * t436) * t440) * t481) * t409;
t1 = [0, 0, t372, t369, t370, 0; 0, 0 (t375 * t571 + t388 * t525) * t550 + ((qJD(5) * t435 + t432 * t481 + t461 * t556) * t388 + t375 * t545 + (t525 * t373 - t375 * t397 - (-t372 * t426 - t381 * t395 + t419 + (-t381 * t446 - t433) * t378) * t565 - (-t372 * t446 - t381 * t413 - t399 + (t381 * t426 - t451) * t378) * t566) * t389) * t383 (t376 * t571 + t388 * t562) * t550 + ((-t422 * t481 - t444 * t552) * t388 + t376 * t545 + (-t376 * t397 + t562 * t373 - (t454 * t552 - t369 * t426 - t382 * t395 + t436 * t481 + (-t382 * t446 - t442 * t481) * t378) * t565 - (-t442 * t552 - t369 * t446 - t382 * t413 - t420 * t481 + (t382 * t426 - t454 * t481) * t378) * t566) * t389) * t383 (t374 * t571 - t388 * t430) * t550 + (t374 * t545 + t398 * t388 + (-t430 * t373 - t374 * t397 - (-t370 * t426 - t380 * t395 + t414 + (-t380 * t446 - t428) * t378) * t565 - (-t370 * t446 - t380 * t413 - t396 + (t380 * t426 - t447) * t378) * t566) * t389) * t383, 0; 0, 0 (-t402 * t417 + t418 * t568) * t547 + ((qJD(6) * t418 + t400 * t480 - t431 * t484) * t402 + t418 * t544 + (-t417 * t387 - (-qJD(6) * t417 + t400 * t484 + t431 * t480) * t405 - t418 * t386) * t403) * t392 (-t402 * t415 + t416 * t568) * t547 + (t416 * t544 - t542 * t402 * t484 + t519 * t569 + (-t405 * t480 * t542 - t416 * t386 - t415 * t387 - t519 * t567) * t403) * t392, t529 * t429 * t548 + (t529 * t397 + ((-qJD(6) * t402 - 0.2e1 * t546) * t484 + (t386 * t484 + (t387 - t551) * t480) * t403) * t429) * t392, t548 + 0.2e1 * (t386 * t403 * t392 + (-t392 * t572 - t403 * t570) * t405) * t405;];
JaD_rot  = t1;
