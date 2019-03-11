% Calculate vector of inverse dynamics joint torques for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:30:08
% EndTime: 2019-03-08 20:30:17
% DurationCPUTime: 6.73s
% Computational Cost: add. (3419->501), mult. (7838->700), div. (0->0), fcn. (6538->16), ass. (0->220)
t498 = cos(pkin(12));
t504 = sin(qJ(2));
t497 = sin(pkin(6));
t584 = qJD(1) * t497;
t563 = t504 * t584;
t462 = t498 * t563;
t495 = sin(pkin(12));
t508 = cos(qJ(2));
t562 = t508 * t584;
t413 = t495 * t562 + t462;
t503 = sin(qJ(4));
t507 = cos(qJ(4));
t544 = pkin(4) * t503 - pkin(9) * t507;
t456 = t544 * qJD(4);
t635 = t413 - t456;
t502 = sin(qJ(5));
t506 = cos(qJ(5));
t570 = t506 * qJD(4);
t582 = qJD(2) * t503;
t445 = t502 * t582 - t570;
t505 = cos(qJ(6));
t579 = qJD(4) * t502;
t447 = t506 * t582 + t579;
t501 = sin(qJ(6));
t614 = t447 * t501;
t396 = t505 * t445 + t614;
t581 = qJD(2) * t507;
t477 = -qJD(5) + t581;
t472 = -qJD(6) + t477;
t634 = t396 * t472;
t534 = t445 * t501 - t505 * t447;
t633 = t472 * t534;
t583 = qJD(2) * t497;
t556 = qJD(1) * t583;
t568 = qJDD(1) * t497;
t632 = t504 * t568 + t508 * t556;
t574 = qJD(5) * t503;
t631 = -qJD(2) * t574 + qJDD(4);
t458 = qJD(2) * pkin(2) + t562;
t410 = t495 * t458 + t462;
t407 = qJD(2) * pkin(8) + t410;
t500 = cos(pkin(6));
t476 = qJD(1) * t500 + qJD(3);
t376 = t507 * t407 + t503 * t476;
t373 = qJD(4) * pkin(9) + t376;
t461 = t495 * t563;
t409 = t458 * t498 - t461;
t532 = -pkin(4) * t507 - pkin(9) * t503 - pkin(3);
t379 = qJD(2) * t532 - t409;
t341 = t373 * t506 + t379 * t502;
t335 = -pkin(10) * t445 + t341;
t572 = qJD(6) * t501;
t333 = t335 * t572;
t625 = -t503 * t407 + t476 * t507;
t372 = -qJD(4) * pkin(4) - t625;
t350 = pkin(5) * t445 + t372;
t499 = cos(pkin(11));
t533 = t495 * t508 + t498 * t504;
t422 = t533 * t500;
t597 = t508 * t498;
t443 = t495 * t504 - t597;
t496 = sin(pkin(11));
t536 = t422 * t499 - t443 * t496;
t607 = t497 * t503;
t363 = -t499 * t607 + t507 * t536;
t535 = -t422 * t496 - t443 * t499;
t365 = t496 * t607 + t507 * t535;
t526 = t443 * t500;
t381 = -t496 * t533 - t499 * t526;
t384 = t496 * t526 - t499 * t533;
t421 = t533 * t497;
t405 = t421 * t507 + t500 * t503;
t606 = t497 * t504;
t420 = t495 * t606 - t497 * t597;
t494 = qJ(5) + qJ(6);
t490 = sin(t494);
t491 = cos(t494);
t630 = t350 * t396 - g(1) * (-t365 * t491 + t384 * t490) - g(2) * (-t363 * t491 + t381 * t490) - g(3) * (-t405 * t491 - t420 * t490) + t333;
t569 = qJD(2) * qJD(4);
t555 = t507 * t569;
t567 = qJDD(2) * t503;
t391 = qJD(5) * t570 + (t555 + t567) * t506 + t631 * t502;
t488 = t507 * qJDD(2);
t442 = t503 * t569 + qJDD(5) - t488;
t473 = t508 * t568;
t417 = qJDD(2) * pkin(2) - t504 * t556 + t473;
t370 = t495 * t417 + t498 * t632;
t368 = qJDD(2) * pkin(8) + t370;
t474 = t500 * qJDD(1) + qJDD(3);
t613 = t474 * t503;
t331 = qJDD(4) * pkin(9) + qJD(4) * t625 + t368 * t507 + t613;
t369 = t417 * t498 - t495 * t632;
t347 = qJD(2) * t456 + qJDD(2) * t532 - t369;
t344 = t506 * t347;
t514 = -qJD(5) * t341 - t502 * t331 + t344;
t320 = pkin(5) * t442 - pkin(10) * t391 + t514;
t392 = t502 * (qJD(4) * (qJD(5) + t581) + t567) - t631 * t506;
t573 = qJD(5) * t506;
t566 = t506 * t331 + t502 * t347 + t379 * t573;
t575 = qJD(5) * t502;
t524 = t373 * t575 - t566;
t321 = -pkin(10) * t392 - t524;
t552 = t505 * t320 - t501 * t321;
t629 = t350 * t534 - g(1) * (-t365 * t490 - t384 * t491) - g(2) * (-t363 * t490 - t381 * t491) - g(3) * (-t405 * t490 + t420 * t491) + t552;
t439 = qJDD(6) + t442;
t628 = t439 * MDP(24) + (-t396 ^ 2 + t534 ^ 2) * MDP(21) - t396 * MDP(20) * t534;
t449 = t501 * t506 + t502 * t505;
t424 = t449 * t503;
t416 = t498 * t562 - t461;
t578 = qJD(4) * t503;
t600 = t502 * t507;
t480 = pkin(2) * t495 + pkin(8);
t610 = t480 * t502;
t627 = -t416 * t600 + t635 * t506 - t578 * t610;
t620 = pkin(2) * t498;
t440 = t532 - t620;
t598 = t506 * t507;
t626 = -t416 * t598 + t440 * t573 - t635 * t502;
t520 = -t502 * t574 + t507 * t570;
t624 = qJD(5) + qJD(6);
t521 = g(1) * t384 + g(2) * t381 - g(3) * t420;
t622 = qJD(5) * (t477 * t480 + t373) - t521;
t550 = t391 * t501 + t505 * t392;
t330 = -qJD(6) * t534 + t550;
t621 = pkin(9) + pkin(10);
t340 = -t373 * t502 + t506 * t379;
t334 = -pkin(10) * t447 + t340;
t328 = -pkin(5) * t477 + t334;
t619 = t328 * t505;
t618 = t335 * t505;
t617 = t391 * t502;
t616 = t445 * t477;
t615 = t447 * t477;
t611 = t477 * t506;
t609 = t490 * t507;
t608 = t491 * t507;
t605 = t497 * t507;
t604 = t500 * t504;
t603 = t500 * t508;
t602 = t501 * t320;
t601 = t502 * t503;
t599 = t503 * t506;
t596 = qJDD(1) - g(3);
t452 = t480 * t598;
t531 = pkin(5) * t503 - pkin(10) * t598;
t595 = -t531 * qJD(4) - (-t452 + (pkin(10) * t503 - t440) * t502) * qJD(5) + t627;
t577 = qJD(4) * t507;
t559 = t502 * t577;
t519 = t503 * t573 + t559;
t594 = (-t503 * t570 - t507 * t575) * t480 - t519 * pkin(10) + t626;
t352 = -t572 * t601 + (t599 * t624 + t559) * t505 + t520 * t501;
t593 = t352 * t472 - t424 * t439;
t453 = t544 * qJD(2);
t592 = t502 * t453 + t506 * t625;
t448 = t501 * t502 - t505 * t506;
t525 = t448 * t507;
t591 = qJD(2) * t525 - t448 * t624;
t590 = (-t581 + t624) * t449;
t587 = t502 * t440 + t452;
t492 = t503 ^ 2;
t586 = -t507 ^ 2 + t492;
t580 = qJD(4) * t445;
t576 = qJD(5) * t477;
t571 = qJD(6) * t505;
t565 = t505 * t391 - t501 * t392 - t445 * t571;
t564 = qJD(5) * t621;
t561 = t502 * t581;
t560 = t477 * t579;
t549 = pkin(5) * t519 - t416 * t503 + t480 * t577;
t548 = qJD(6) * t328 + t321;
t545 = -t376 + (-t561 + t575) * pkin(5);
t436 = t506 * t453;
t470 = t621 * t506;
t543 = qJD(2) * t531 + qJD(6) * t470 - t502 * t625 + t506 * t564 + t436;
t469 = t621 * t502;
t542 = pkin(10) * t561 - qJD(6) * t469 - t502 * t564 - t592;
t323 = t328 * t501 + t618;
t351 = -qJD(4) * t525 - t424 * t624;
t425 = t448 * t503;
t541 = t351 * t472 + t425 * t439;
t357 = -t405 * t502 + t420 * t506;
t358 = t405 * t506 + t420 * t502;
t540 = t357 * t505 - t358 * t501;
t539 = t357 * t501 + t358 * t505;
t427 = t506 * t440;
t378 = -pkin(10) * t599 + t427 + (-pkin(5) - t610) * t507;
t390 = -pkin(10) * t601 + t587;
t538 = t378 * t501 + t390 * t505;
t404 = t421 * t503 - t500 * t507;
t530 = t503 * t368 + t407 * t577 - t474 * t507 + t476 * t578;
t528 = -t442 * t502 + t477 * t573;
t527 = -t442 * t506 - t477 * t575;
t329 = -t447 * t572 + t565;
t523 = g(1) * (t496 * t605 - t503 * t535) + g(2) * (-t499 * t605 - t503 * t536) - g(3) * t404;
t522 = -g(1) * t535 - g(2) * t536 - g(3) * t421;
t332 = -qJDD(4) * pkin(4) + t530;
t517 = -pkin(9) * t442 - t372 * t477;
t406 = -qJD(2) * pkin(3) - t409;
t481 = -pkin(3) - t620;
t516 = -qJDD(4) * t480 + (qJD(2) * t481 + t406 + t416) * qJD(4);
t513 = pkin(9) * t576 - t332 - t523;
t509 = qJD(4) ^ 2;
t512 = -qJD(2) * t413 + t480 * t509 - t369 + t521 + (-pkin(3) + t481) * qJDD(2);
t511 = -g(1) * (-t496 * t603 - t499 * t504) - g(2) * (-t496 * t504 + t499 * t603) - g(3) * t497 * t508;
t510 = qJD(2) ^ 2;
t484 = -pkin(5) * t506 - pkin(4);
t465 = qJDD(4) * t507 - t503 * t509;
t464 = qJDD(4) * t503 + t507 * t509;
t430 = (pkin(5) * t502 + t480) * t503;
t428 = t447 * t578;
t415 = t443 * t583;
t414 = qJD(2) * t421;
t389 = t534 * t578;
t356 = -qJD(4) * t404 - t415 * t507;
t355 = qJD(4) * t405 - t415 * t503;
t326 = qJD(5) * t357 + t356 * t506 + t414 * t502;
t325 = -qJD(5) * t358 - t356 * t502 + t414 * t506;
t324 = pkin(5) * t392 + t332;
t322 = -t335 * t501 + t619;
t1 = [t596 * MDP(1) + (-t369 * t420 + t370 * t421 - t409 * t414 - t410 * t415 + t474 * t500 - g(3)) * MDP(5) + (-qJD(4) * t355 - qJDD(4) * t404 - t420 * t488) * MDP(11) + (-qJD(4) * t356 - qJDD(4) * t405 + t420 * t567) * MDP(12) + (-t325 * t477 + t355 * t445 + t357 * t442 + t392 * t404) * MDP(18) + (t326 * t477 + t355 * t447 - t358 * t442 + t391 * t404) * MDP(19) + (-(-qJD(6) * t539 + t325 * t505 - t326 * t501) * t472 + t540 * t439 + t355 * t396 + t404 * t330) * MDP(25) + ((qJD(6) * t540 + t325 * t501 + t326 * t505) * t472 - t539 * t439 - t355 * t534 + t404 * t329) * MDP(26) + ((-t414 * t507 + t420 * t578) * MDP(11) + (t414 * t503 + t420 * t577) * MDP(12)) * qJD(2) + ((qJDD(2) * t508 - t504 * t510) * MDP(3) + (-qJDD(2) * t504 - t508 * t510) * MDP(4)) * t497; qJDD(2) * MDP(2) + (t473 + t511) * MDP(3) + (-g(1) * (t496 * t604 - t499 * t508) - g(2) * (-t496 * t508 - t499 * t604) - t596 * t606) * MDP(4) + (t409 * t413 - t410 * t416 + (t369 * t498 + t370 * t495 + t511) * pkin(2)) * MDP(5) + (qJDD(2) * t492 + 0.2e1 * t503 * t555) * MDP(6) + 0.2e1 * (t488 * t503 - t569 * t586) * MDP(7) + t464 * MDP(8) + t465 * MDP(9) + (t503 * t516 - t507 * t512) * MDP(11) + (t503 * t512 + t507 * t516) * MDP(12) + (t391 * t599 + t447 * t520) * MDP(13) + ((-t445 * t506 - t447 * t502) * t577 + (-t617 - t392 * t506 + (t445 * t502 - t447 * t506) * qJD(5)) * t503) * MDP(14) + (-t391 * t507 + t442 * t599 - t477 * t520 + t428) * MDP(15) + ((t392 + t560) * t507 + (t528 - t580) * t503) * MDP(16) + (-t442 * t507 - t477 * t578) * MDP(17) + (t427 * t442 + t627 * t477 + (t440 * t576 + t522) * t502 + (t480 * t580 - t344 + (qJD(4) * t372 + qJD(5) * t379 - t442 * t480 + t331) * t502 + t622 * t506) * t507 + (qJD(4) * t340 + t332 * t502 + t372 * t573 + t392 * t480 - t416 * t445) * t503) * MDP(18) + (-t587 * t442 + t626 * t477 + t522 * t506 + ((t372 * t506 + t447 * t480) * qJD(4) - t622 * t502 + t566) * t507 + (-t372 * t575 + t332 * t506 + t480 * t391 - t416 * t447 + (-t480 * t611 - t341) * qJD(4)) * t503) * MDP(19) + (-t329 * t425 - t351 * t534) * MDP(20) + (-t329 * t424 + t330 * t425 - t351 * t396 + t352 * t534) * MDP(21) + (-t329 * t507 - t389 - t541) * MDP(22) + (t330 * t507 - t396 * t578 + t593) * MDP(23) + (-t439 * t507 - t472 * t578) * MDP(24) + ((t378 * t505 - t390 * t501) * t439 - t552 * t507 + t322 * t578 + t430 * t330 + t324 * t424 + t350 * t352 - g(1) * (t384 * t608 + t490 * t535) - g(2) * (t381 * t608 + t490 * t536) - g(3) * (-t420 * t608 + t421 * t490) + (t501 * t594 + t505 * t595) * t472 + t549 * t396 + (t323 * t507 + t472 * t538) * qJD(6)) * MDP(25) + (-t538 * t439 + (t548 * t505 - t333 + t602) * t507 - t323 * t578 + t430 * t329 - t324 * t425 + t350 * t351 - g(1) * (-t384 * t609 + t491 * t535) - g(2) * (-t381 * t609 + t491 * t536) - g(3) * (t420 * t609 + t421 * t491) + ((qJD(6) * t378 + t594) * t505 + (-qJD(6) * t390 - t595) * t501) * t472 - t549 * t534) * MDP(26); (-g(3) * t500 + (-g(1) * t496 + g(2) * t499) * t497 + t474) * MDP(5) + t465 * MDP(11) - t464 * MDP(12) + t428 * MDP(19) + t593 * MDP(25) + (-t389 + t541) * MDP(26) + ((-t392 + t560) * MDP(18) + (t477 * t570 - t391) * MDP(19) - t330 * MDP(25) - t329 * MDP(26)) * t507 + ((t528 + t580) * MDP(18) + t527 * MDP(19) + qJD(4) * t396 * MDP(25)) * t503; MDP(8) * t567 + MDP(9) * t488 + qJDD(4) * MDP(10) + (qJD(4) * t376 - t523 - t530) * MDP(11) + (g(1) * t365 + g(2) * t363 + g(3) * t405 - t613 + (-qJD(2) * t406 - t368) * t507) * MDP(12) + (-t447 * t611 + t617) * MDP(13) + ((t391 + t616) * t506 + (-t392 + t615) * t502) * MDP(14) + ((-t447 * t503 + t477 * t598) * qJD(2) - t528) * MDP(15) + ((t445 * t503 - t477 * t600) * qJD(2) - t527) * MDP(16) + (-pkin(4) * t392 - t376 * t445 + t436 * t477 + (-t477 * t625 + t517) * t502 + t513 * t506) * MDP(18) + (-pkin(4) * t391 - t376 * t447 - t477 * t592 - t502 * t513 + t506 * t517) * MDP(19) + (t329 * t449 - t534 * t591) * MDP(20) + (-t329 * t448 - t330 * t449 - t396 * t591 + t534 * t590) * MDP(21) + (t439 * t449 - t472 * t591) * MDP(22) + (-t439 * t448 + t472 * t590) * MDP(23) + ((-t469 * t505 - t470 * t501) * t439 + t484 * t330 + t324 * t448 + (t501 * t542 + t505 * t543) * t472 + t545 * t396 + t590 * t350 - t523 * t491) * MDP(25) + (-(-t469 * t501 + t470 * t505) * t439 + t484 * t329 + t324 * t449 + (-t501 * t543 + t505 * t542) * t472 - t545 * t534 + t591 * t350 + t523 * t490) * MDP(26) + (-MDP(11) * t406 + t477 * MDP(17) - t340 * MDP(18) + MDP(19) * t341 + MDP(22) * t534 + t396 * MDP(23) + t472 * MDP(24) - t322 * MDP(25) + t323 * MDP(26)) * t582 + (-MDP(6) * t503 * t507 + MDP(7) * t586) * t510; t447 * t445 * MDP(13) + (-t445 ^ 2 + t447 ^ 2) * MDP(14) + (t391 - t616) * MDP(15) + (-t392 - t615) * MDP(16) + t442 * MDP(17) + (-t341 * t477 - t372 * t447 - g(1) * (-t365 * t502 - t384 * t506) - g(2) * (-t363 * t502 - t381 * t506) - g(3) * t357 + t514) * MDP(18) + (-t340 * t477 + t372 * t445 - g(1) * (-t365 * t506 + t384 * t502) - g(2) * (-t363 * t506 + t381 * t502) + g(3) * t358 + t524) * MDP(19) + (t329 - t634) * MDP(22) + (-t330 + t633) * MDP(23) + ((-t334 * t501 - t618) * t472 - t323 * qJD(6) + (-t396 * t447 + t439 * t505 + t472 * t572) * pkin(5) + t629) * MDP(25) + ((t335 * t472 - t320) * t501 + (-t334 * t472 - t548) * t505 + (-t439 * t501 + t447 * t534 + t472 * t571) * pkin(5) + t630) * MDP(26) + t628; (t565 - t634) * MDP(22) + (-t550 + t633) * MDP(23) + (-t323 * t472 + t629) * MDP(25) + (-t505 * t321 - t322 * t472 - t602 + t630) * MDP(26) + (-MDP(22) * t614 + MDP(23) * t534 - MDP(25) * t323 - MDP(26) * t619) * qJD(6) + t628;];
tau  = t1;
