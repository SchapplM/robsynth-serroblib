% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR12_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobiaD_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:11
% EndTime: 2019-02-26 21:21:16
% DurationCPUTime: 4.82s
% Computational Cost: add. (23019->195), mult. (70693->359), div. (705->12), fcn. (90087->19), ass. (0->173)
t461 = sin(qJ(3));
t579 = sin(pkin(14));
t585 = cos(pkin(6));
t533 = t585 * t579;
t583 = cos(pkin(14));
t586 = sin(qJ(1));
t589 = cos(qJ(1));
t502 = t589 * t533 + t586 * t583;
t581 = sin(pkin(7));
t582 = sin(pkin(6));
t531 = t582 * t581;
t520 = t589 * t531;
t588 = cos(qJ(3));
t511 = t588 * t520;
t584 = cos(pkin(7));
t535 = t585 * t583;
t601 = -t589 * t535 + t586 * t579;
t602 = t601 * t584;
t435 = t502 * t461 + t588 * t602 + t511;
t504 = -t586 * t533 + t589 * t583;
t451 = t504 * qJD(1);
t516 = t586 * t531;
t503 = t586 * t535 + t589 * t579;
t450 = t503 * qJD(1);
t547 = t450 * t584;
t418 = -t451 * t588 + (-qJD(1) * t516 + t547) * t461 + t435 * qJD(3);
t458 = cos(pkin(8));
t460 = sin(qJ(4));
t532 = t584 * t582;
t517 = t586 * t532;
t505 = qJD(1) * t517 + t450 * t581;
t580 = sin(pkin(8));
t491 = t505 * t580;
t587 = cos(qJ(4));
t436 = (t602 + t520) * t461 - t502 * t588;
t510 = t588 * t516;
t598 = qJD(1) * t510 + t436 * qJD(3) - t451 * t461 - t588 * t547;
t619 = -t418 * t587 + (t598 * t458 + t491) * t460;
t593 = -t589 * t532 + t581 * t601;
t484 = t593 * t580;
t410 = t436 * t587 + (t435 * t458 - t484) * t460;
t551 = t458 * t587;
t618 = -t410 * qJD(4) - t418 * t460 - t598 * t551;
t613 = t436 * t460;
t610 = t583 * t532 + t585 * t581;
t490 = t503 * t584;
t437 = t504 * t588 + (-t490 + t516) * t461;
t475 = t504 * t461 + t588 * t490 - t510;
t483 = t503 * t581 + t517;
t481 = t483 * t580;
t412 = t437 * t587 + (-t475 * t458 + t481) * t460;
t429 = t483 * t458 + t475 * t580;
t459 = sin(qJ(5));
t462 = cos(qJ(5));
t389 = t412 * t459 - t429 * t462;
t609 = 0.2e1 * t389;
t476 = t435 * t587;
t406 = t458 * t476 - t587 * t484 - t613;
t404 = t406 ^ 2;
t530 = t582 * t579;
t446 = t610 * t461 + t588 * t530;
t452 = -t583 * t531 + t585 * t584;
t539 = t587 * t580;
t445 = -t461 * t530 + t610 * t588;
t550 = t587 * t445;
t498 = -t446 * t460 + t452 * t539 + t458 * t550;
t423 = 0.1e1 / t498 ^ 2;
t395 = t404 * t423 + 0.1e1;
t393 = 0.1e1 / t395;
t378 = -t587 * t491 + t618;
t426 = t446 * t587 + (t445 * t458 + t580 * t452) * t460;
t442 = t445 * qJD(3);
t443 = t446 * qJD(3);
t399 = t426 * qJD(4) + t442 * t460 + t443 * t551;
t422 = 0.1e1 / t498;
t566 = t406 * t423;
t527 = t378 * t422 + t399 * t566;
t360 = t527 * t393;
t396 = atan2(-t406, -t498);
t391 = sin(t396);
t392 = cos(t396);
t529 = t391 * t498 - t392 * t406;
t355 = t529 * t360 - t378 * t391 + t392 * t399;
t372 = -t391 * t406 - t392 * t498;
t370 = 0.1e1 / t372 ^ 2;
t605 = t355 * t370;
t604 = t399 * t423;
t449 = t502 * qJD(1);
t489 = qJD(1) * t602;
t471 = -qJD(1) * t511 + qJD(3) * t437 - t449 * t461 - t588 * t489;
t473 = t475 * t587;
t600 = -qJD(4) * t473 - t471 * t460;
t599 = -t437 * t551 + t475 * t460;
t479 = t587 * t481;
t565 = t437 * t460;
t411 = t458 * t473 - t479 + t565;
t405 = t411 ^ 2;
t366 = t370 * t405 + 0.1e1;
t364 = 0.1e1 / t366;
t369 = 0.1e1 / t372;
t415 = -t449 * t588 + (qJD(1) * t520 + t489) * t461 - t475 * qJD(3);
t469 = t471 * t587;
t482 = qJD(1) * t593;
t480 = t580 * t482;
t376 = t412 * qJD(4) + t415 * t460 + t458 * t469 + t587 * t480;
t570 = t376 * t370;
t577 = t369 * t605;
t578 = (-t405 * t577 + t411 * t570) / t366 ^ 2;
t597 = -t364 * t605 - 0.2e1 * t369 * t578;
t590 = 0.2e1 * t411;
t544 = t577 * t590;
t561 = 0.2e1 * t578;
t572 = t370 * t411;
t596 = t561 * t572 + (t544 - t570) * t364;
t390 = t412 * t462 + t429 * t459;
t384 = 0.1e1 / t390;
t385 = 0.1e1 / t390 ^ 2;
t591 = -0.2e1 * t406;
t549 = qJD(4) * t565;
t377 = qJD(4) * t479 + t415 * t587 + t600 * t458 - t460 * t480 - t549;
t401 = -t458 * t482 + t471 * t580;
t367 = t390 * qJD(5) + t377 * t459 - t401 * t462;
t383 = t389 ^ 2;
t375 = t383 * t385 + 0.1e1;
t569 = t385 * t389;
t562 = qJD(5) * t389;
t368 = t377 * t462 + t401 * t459 - t562;
t573 = t368 * t384 * t385;
t576 = (t367 * t569 - t383 * t573) / t375 ^ 2;
t568 = t422 * t604;
t575 = (t378 * t566 + t404 * t568) / t395 ^ 2;
t574 = t364 * t369;
t373 = 0.1e1 / t375;
t571 = t373 * t385;
t567 = t406 * t422;
t564 = t458 * t460;
t560 = -0.2e1 * t576;
t559 = -0.2e1 * t575;
t557 = t385 * t576;
t556 = t422 * t575;
t554 = t364 * t572;
t553 = t367 * t571;
t552 = t389 * t573;
t548 = t437 * t580;
t546 = t580 * t415;
t543 = 0.2e1 * t552;
t542 = t568 * t591;
t428 = -t435 * t580 - t458 * t593;
t388 = t410 * t462 + t428 * t459;
t387 = t410 * t459 - t428 * t462;
t526 = -t459 * t384 + t462 * t569;
t525 = -t410 * t422 + t426 * t566;
t419 = -t435 * t460 - t436 * t551;
t430 = t445 * t460 + t446 * t551;
t524 = t419 * t422 + t430 * t566;
t421 = -t437 * t564 - t473;
t398 = t421 * t462 + t459 * t548;
t514 = -t421 * t459 + t462 * t548;
t513 = -t391 + (-t392 * t567 + t391) * t393;
t499 = t435 * t551 - t539 * t593 - t613;
t403 = t442 * t551 - t443 * t460 + (-t446 * t564 + t550) * qJD(4);
t402 = -t458 * t505 + t580 * t598;
t400 = t498 * qJD(4) + t442 * t587 - t443 * t564;
t382 = -t418 * t551 + (t436 * t564 - t476) * qJD(4) + t598 * t460;
t381 = t599 * qJD(4) - t415 * t564 - t469;
t380 = t499 * qJD(4) - t619;
t379 = -t406 * qJD(4) + t619;
t363 = t524 * t393;
t362 = t525 * t393;
t356 = t529 * t362 + t391 * t410 + t392 * t426;
t354 = t524 * t559 + (-t430 * t542 + t382 * t422 + (t378 * t430 + t399 * t419 + t403 * t406) * t423) * t393;
t352 = t525 * t559 + (-t426 * t542 + t379 * t422 + (t378 * t426 - t399 * t410 + t400 * t406) * t423) * t393;
t1 = [-t556 * t590 + (t376 * t422 + t411 * t604) * t393, 0, t354, t352, 0, 0; (t505 * t539 - t618) * t574 - (t513 * t376 + ((t360 * t393 * t567 + t559) * t391 + (-t556 * t591 - t360 + (t360 - t527) * t393) * t392) * t411) * t554 - t597 * t499 + t596 * t513 * t411, 0 (t415 * t551 - t458 * t549 + t600) * t574 - ((-t354 * t406 - t363 * t378 + t403 + (t363 * t498 - t419) * t360) * t392 + (t354 * t498 - t363 * t399 - t382 + (t363 * t406 - t430) * t360) * t391) * t554 - t597 * t599 + t596 * (t529 * t363 - t391 * t419 + t392 * t430) (t356 * t572 - t369 * t412) * t561 + (t356 * t544 + t377 * t369 + (-t412 * t355 - t356 * t376 + (-(-t352 * t406 - t362 * t378 + t400 + (t362 * t498 + t410) * t360) * t392 - (t352 * t498 - t362 * t399 - t379 + (t362 * t406 - t426) * t360) * t391) * t411) * t370) * t364, 0, 0; 0.2e1 * (-t384 * t387 + t388 * t569) * t576 + ((t388 * qJD(5) + t380 * t459 - t402 * t462) * t384 + t388 * t543 + (-t387 * t368 - (-t387 * qJD(5) + t380 * t462 + t402 * t459) * t389 - t388 * t367) * t385) * t373, 0 (t557 * t609 - t553) * t398 - (-t368 * t571 + t384 * t560) * t514 + ((t398 * qJD(5) + t381 * t459 - t462 * t546) * t384 - (t514 * qJD(5) + t381 * t462 + t459 * t546) * t569 + t398 * t543) * t373, t526 * t411 * t560 + (t526 * t376 + ((-qJD(5) * t384 - 0.2e1 * t552) * t462 + (t367 * t462 + (t368 - t562) * t459) * t385) * t411) * t373, t560 + (t553 + (-t373 * t573 - t557) * t389) * t609, 0;];
JaD_rot  = t1;
