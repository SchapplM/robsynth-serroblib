% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:13:48
% EndTime: 2019-03-09 05:13:59
% DurationCPUTime: 8.35s
% Computational Cost: add. (6510->507), mult. (15679->607), div. (0->0), fcn. (12502->14), ass. (0->225)
t496 = qJD(3) + qJD(4);
t501 = sin(qJ(4));
t498 = sin(pkin(10));
t499 = cos(pkin(10));
t624 = sin(qJ(3));
t626 = cos(qJ(3));
t452 = t498 * t626 + t499 * t624;
t524 = t452 * qJDD(1);
t572 = t624 * t498;
t574 = t626 * t499;
t535 = t572 - t574;
t446 = t535 * qJD(1);
t636 = t446 * qJD(3);
t512 = t524 - t636;
t568 = qJD(3) * t624;
t569 = qJD(3) * t626;
t575 = -qJDD(1) * t572 + (-t498 * t569 - t499 * t568) * qJD(1);
t532 = -qJDD(1) * t574 - t575;
t625 = cos(qJ(4));
t510 = t501 * t512 + t625 * t532;
t447 = t452 * qJD(1);
t540 = -t501 * t446 + t447 * t625;
t378 = qJD(4) * t540 + t510;
t418 = t625 * t446 + t447 * t501;
t491 = qJDD(3) + qJDD(4);
t500 = sin(qJ(6));
t503 = cos(qJ(6));
t584 = qJD(6) * t503;
t576 = t500 * t378 + t418 * t584 + t503 * t491;
t585 = qJD(6) * t500;
t351 = -t496 * t585 + t576;
t350 = t351 * t503;
t373 = t503 * t378;
t409 = t418 * t500 + t496 * t503;
t352 = t409 * qJD(6) + t491 * t500 - t373;
t567 = qJD(4) * t625;
t588 = qJD(4) * t501;
t509 = t446 * t588 - t447 * t567 - t510;
t377 = t446 * t567 + t447 * t588 + t501 * t532 - t625 * t512;
t609 = t418 * t496;
t523 = -t377 + t609;
t611 = t540 * t496;
t635 = qJD(6) + t540;
t651 = t635 * t503;
t652 = t635 * t500;
t407 = -t503 * t418 + t496 * t500;
t653 = t407 * t635;
t656 = t523 * MDP(17) + t491 * MDP(19) + (t509 + t611) * MDP(18) + (-t409 * t652 + t350) * MDP(26) + (-t409 * t651 + (-t351 + t653) * t500 - t503 * t352) * MDP(27);
t618 = qJDD(1) * pkin(1);
t502 = sin(qJ(1));
t504 = cos(qJ(1));
t637 = -g(1) * t502 + g(2) * t504;
t542 = -qJDD(2) + t618 - t637;
t374 = -qJDD(6) + t377;
t368 = t503 * t374;
t650 = -t635 * t652 - t368;
t616 = t374 * t500;
t533 = -t635 * t651 + t616;
t622 = pkin(5) * t418;
t619 = qJ(5) * t418;
t620 = pkin(7) + qJ(2);
t464 = t620 * t498;
t453 = qJD(1) * t464;
t465 = t620 * t499;
t454 = qJD(1) * t465;
t406 = -t446 * pkin(8) - t453 * t624 + t454 * t626;
t403 = t625 * t406;
t639 = -t626 * t453 - t624 * t454;
t405 = -t447 * pkin(8) + t639;
t404 = qJD(3) * pkin(3) + t405;
t381 = t501 * t404 + t403;
t379 = -qJ(5) * t496 - t381;
t355 = -t379 - t622;
t654 = t355 * t635;
t495 = pkin(10) + qJ(3);
t488 = cos(t495);
t489 = qJ(4) + t495;
t478 = sin(t489);
t479 = cos(t489);
t592 = t479 * pkin(4) + t478 * qJ(5);
t649 = pkin(3) * t488 + t592;
t648 = t499 * MDP(4) - t498 * MDP(5);
t402 = t501 * t406;
t380 = -t625 * t404 + t402;
t581 = qJD(5) + t380;
t645 = t540 ^ 2;
t644 = pkin(4) * t540;
t643 = pkin(5) * t540;
t642 = t355 * t540;
t627 = pkin(4) + pkin(9);
t640 = t540 * t627;
t383 = t501 * t405 + t403;
t554 = pkin(3) * t588 - t383;
t384 = t405 * t625 - t402;
t594 = -pkin(3) * t567 - qJD(5) + t384;
t484 = t491 * qJ(5);
t638 = -t496 * qJD(5) - t484;
t582 = t643 + t581;
t634 = qJ(2) * qJDD(1);
t481 = -t499 * pkin(2) - pkin(1);
t457 = qJD(1) * t481 + qJD(2);
t428 = pkin(3) * t446 + t457;
t518 = -qJ(5) * t540 + t428;
t382 = pkin(4) * t418 + t518;
t474 = g(3) * t479;
t578 = qJD(1) * qJD(2);
t631 = qJDD(1) * t620 + t578;
t432 = t631 * t498;
t433 = t631 * t499;
t550 = -t626 * t432 - t433 * t624;
t371 = qJDD(3) * pkin(3) - pkin(8) * t512 + t453 * t568 - t454 * t569 + t550;
t536 = -t432 * t624 + t433 * t626;
t376 = -t532 * pkin(8) + qJD(3) * t639 + t536;
t559 = -t625 * t371 + t501 * t376 + t404 * t588 + t406 * t567;
t605 = t478 * t504;
t606 = t478 * t502;
t538 = -g(1) * t605 - g(2) * t606 + t474 + t559;
t519 = t382 * t540 + qJDD(5) + t538;
t483 = -pkin(3) * t625 - pkin(4);
t476 = -pkin(9) + t483;
t633 = t635 * (t554 + t622) - t476 * t374;
t632 = -t428 * t540 - t538;
t353 = -t496 * t627 + t582;
t358 = t418 * t627 + t518;
t341 = t353 * t503 - t358 * t500;
t342 = t353 * t500 + t358 * t503;
t630 = MDP(15) * t540 + t428 * MDP(21) - t382 * MDP(24) + MDP(30) * t635 + t341 * MDP(31) - t342 * MDP(32);
t628 = qJD(3) ^ 2;
t623 = pkin(4) * t491;
t473 = g(3) * t478;
t520 = t625 * t535;
t426 = t452 * t501 + t520;
t529 = t501 * t535;
t427 = t452 * t625 - t529;
t430 = pkin(3) * t535 + t481;
t513 = -t427 * qJ(5) + t430;
t362 = t426 * t627 + t513;
t617 = t362 * t374;
t615 = t381 * t496;
t614 = t407 * t418;
t613 = t409 * t418;
t610 = t418 ^ 2;
t608 = t426 * t500;
t604 = t479 * t502;
t603 = t479 * t504;
t600 = t500 * t502;
t599 = t500 * t504;
t598 = t502 * t503;
t597 = t503 * t504;
t558 = -t501 * t371 - t625 * t376 - t404 * t567 + t406 * t588;
t338 = t558 + t638;
t336 = -pkin(5) * t378 - t338;
t596 = t336 * t500 + t355 * t584;
t595 = -t624 * t464 + t626 * t465;
t593 = -t594 + t643;
t591 = t498 ^ 2 + t499 ^ 2;
t587 = qJD(6) * t358;
t586 = qJD(6) * t496;
t583 = t447 * qJD(3);
t571 = qJD(2) * t626;
t570 = qJD(2) * t624;
t565 = t591 * qJD(1) ^ 2;
t410 = qJDD(2) - t575 * pkin(3) + (-pkin(1) + (-pkin(3) * t626 - pkin(2)) * t499) * qJDD(1);
t511 = t377 * qJ(5) - qJD(5) * t540 + t410;
t337 = t378 * t627 + t511;
t560 = qJD(6) * t353 + t337;
t557 = 0.2e1 * t591;
t487 = sin(t495);
t553 = -pkin(3) * t487 - pkin(4) * t478;
t552 = g(1) * t504 + g(2) * t502;
t549 = -t626 * t464 - t465 * t624;
t548 = pkin(3) * t447 + t619;
t547 = -t587 + t474;
t516 = -t464 * t569 - t465 * t568 - t498 * t570 + t499 * t571;
t526 = qJD(3) * t452;
t394 = -pkin(8) * t526 + t516;
t525 = qJD(3) * t535;
t395 = pkin(8) * t525 + t464 * t568 - t465 * t569 - t498 * t571 - t499 * t570;
t414 = -t452 * pkin(8) + t549;
t415 = -pkin(8) * t535 + t595;
t346 = -t625 * t394 - t501 * t395 - t414 * t567 + t415 * t588;
t388 = t501 * t414 + t415 * t625;
t546 = -t346 * t496 + t388 * t491;
t347 = qJD(4) * t388 + t501 * t394 - t395 * t625;
t387 = -t414 * t625 + t501 * t415;
t545 = t347 * t496 + t387 * t491;
t544 = (t381 - t622) * t635 - t627 * t374;
t543 = qJDD(5) + t559;
t541 = -t481 + t649;
t391 = -qJD(4) * t529 + t452 * t567 - t501 * t525 + t526 * t625;
t539 = t391 * t500 + t426 * t584;
t537 = g(1) * t603 + g(2) * t604 + t473 + t558;
t363 = t427 * pkin(5) + t387;
t534 = t336 * t426 + t355 * t391 + t363 * t374;
t528 = -t479 * t552 - t473;
t527 = -t537 - t638;
t521 = pkin(3) * t526;
t517 = t557 * t578 - t552;
t515 = (-qJD(6) * t476 + t548 + t640) * t635 + t528;
t514 = (qJD(6) * t627 + t619 + t640) * t635 + t528;
t390 = t452 * t588 + t496 * t520 + t501 * t526;
t348 = t391 * pkin(4) + t390 * qJ(5) - t427 * qJD(5) + t521;
t340 = t378 * pkin(4) + t511;
t492 = -pkin(8) - t620;
t480 = pkin(3) * t501 + qJ(5);
t459 = qJ(5) * t603;
t458 = qJ(5) * t604;
t456 = qJDD(1) * t481 + qJDD(2);
t442 = -t478 * t600 + t597;
t441 = t478 * t598 + t599;
t440 = t478 * t599 + t598;
t439 = t478 * t597 - t600;
t389 = t619 + t644;
t386 = t426 * pkin(4) + t513;
t385 = t548 + t644;
t375 = -pkin(4) * t496 + t581;
t364 = -t426 * pkin(5) + t388;
t345 = t391 * pkin(9) + t348;
t344 = -t390 * pkin(5) + t347;
t343 = -pkin(5) * t391 - t346;
t339 = t543 - t623;
t335 = -pkin(5) * t377 - t491 * t627 + t543;
t334 = t336 * t503;
t332 = t503 * t335;
t1 = [(-t338 * t388 + t339 * t387 + t340 * t386 + t379 * t346 + t375 * t347 + t382 * t348 + (g(1) * t492 - g(2) * t541) * t504 + (g(1) * t541 + g(2) * t492) * t502) * MDP(25) + (t452 * qJDD(3) - t535 * t628) * MDP(10) + (-qJDD(3) * t535 - t452 * t628) * MDP(11) + ((-t407 * t500 + t409 * t503) * t391 + (t350 - t352 * t500 + (-t407 * t503 - t409 * t500) * qJD(6)) * t426) * MDP(27) + (t351 * t608 + t409 * t539) * MDP(26) - t637 * MDP(2) + qJDD(1) * MDP(1) + t648 * (t542 + t618) + (-qJD(3) * t516 - qJDD(3) * t595 + t456 * t452 - t457 * t525 + t481 * t512 + t487 * t637) * MDP(14) + (t549 * qJDD(3) + t456 * t535 + t481 * t532 - t637 * t488 + (-t595 * qJD(3) + (-qJD(2) + t457) * t452) * qJD(3)) * MDP(13) + (pkin(1) * t542 + (t591 * t634 + t517) * qJ(2)) * MDP(7) + (t557 * t634 + t517) * MDP(6) + t552 * MDP(3) + (g(1) * t604 - g(2) * t603 + t430 * t378 + t428 * t391 + t410 * t426 + t418 * t521 - t545) * MDP(20) + (-g(1) * t606 + g(2) * t605 - t430 * t377 - t428 * t390 + t410 * t427 + t521 * t540 - t546) * MDP(21) + (t446 * t525 - t447 * t526 - t452 * t532 - t512 * t535) * MDP(9) + (t338 * t426 + t339 * t427 + t346 * t418 + t347 * t540 - t375 * t390 - t377 * t387 - t378 * t388 + t379 * t391 - t552) * MDP(22) + (-t377 * t427 - t390 * t540) * MDP(15) + (t377 * t426 - t378 * t427 + t390 * t418 - t391 * t540) * MDP(16) + (t452 * t524 + (-t446 * t452 - t447 * t535) * qJD(3)) * MDP(8) + (-g(1) * t442 - g(2) * t440 + t332 * t427 - t341 * t390 + t343 * t407 + t364 * t352 + (-t337 * t427 - t345 * t635 + t617) * t500 + (t344 * t635 - t534) * t503 + ((-t362 * t503 - t363 * t500) * t635 - t342 * t427 + t355 * t608) * qJD(6)) * MDP(31) + (g(1) * t441 - g(2) * t439 + t342 * t390 + t343 * t409 + t364 * t351 + (-(qJD(6) * t363 + t345) * t635 + t617 - t560 * t427 + t355 * qJD(6) * t426) * t503 + (-(-qJD(6) * t362 + t344) * t635 - (t335 - t587) * t427 + t534) * t500) * MDP(32) + (t351 * t427 - t374 * t608 - t390 * t409 + t539 * t635) * MDP(28) + (-t426 * t368 - t352 * t427 + t390 * t407 + (t391 * t503 - t426 * t585) * t635) * MDP(29) + (-t374 * t427 - t390 * t635) * MDP(30) + (-t391 * t496 - t426 * t491) * MDP(18) + (-t390 * t496 + t427 * t491) * MDP(17) + (-t340 * t426 - t348 * t418 - t378 * t386 - t382 * t391 + t479 * t637 + t545) * MDP(23) + (-t340 * t427 - t348 * t540 + t377 * t386 + t382 * t390 - t478 * t637 + t546) * MDP(24); -MDP(6) * t565 + (-qJ(2) * t565 - t542) * MDP(7) + (t532 + t583) * MDP(13) + (t524 - 0.2e1 * t636) * MDP(14) + (-t610 - t645) * MDP(22) + (-t375 * t540 - t379 * t418 + t340 + t637) * MDP(25) + (t533 + t614) * MDP(31) + (t613 - t650) * MDP(32) + (-MDP(21) + MDP(24)) * (t377 + t609) + (MDP(20) - MDP(23)) * (-t509 + t611) - t648 * qJDD(1); t524 * MDP(10) + (t383 * t496 + (-t418 * t447 + t491 * t625 - t496 * t588) * pkin(3) + t632) * MDP(20) + (t385 * t418 + t554 * t496 + (-pkin(4) + t483) * t491 + t519) * MDP(23) + (-t338 * t480 + t339 * t483 - t382 * t385 - g(1) * (t504 * t553 + t459) - g(2) * (t502 * t553 + t458) - g(3) * t649 + t594 * t379 + t554 * t375) * MDP(25) + (t613 + t650) * MDP(28) + (-g(3) * t488 - t457 * t447 + t487 * t552 + t550) * MDP(13) + (t480 * t351 + t334 + t593 * t409 + t515 * t503 + (-t633 - t654) * t500) * MDP(32) + qJDD(3) * MDP(12) + (t533 - t614) * MDP(29) + (t385 * t540 + t480 * t491 - t496 * t594 + t527) * MDP(24) + (t384 * t496 + (-t447 * t540 - t491 * t501 - t496 * t567) * pkin(3) + t537) * MDP(21) + (g(3) * t487 + t457 * t446 + t488 * t552 - t536) * MDP(14) - (-MDP(22) * t375 - t630) * t418 + (-t610 + t645) * MDP(16) + (-t532 + t583) * MDP(11) + (t480 * t352 + t593 * t407 + (t633 + t642) * t503 + t515 * t500 + t596) * MDP(31) + t447 * t446 * MDP(8) + (-t446 ^ 2 + t447 ^ 2) * MDP(9) + (-t377 * t483 - t378 * t480 + t418 * t594 + (-t379 + t554) * t540) * MDP(22) + t656; t645 * MDP(16) + (t615 + t632) * MDP(20) + (-t380 * t496 + t537) * MDP(21) + (pkin(4) * t377 - qJ(5) * t378 + (-t379 - t381) * t540) * MDP(22) + (t519 - t615 - 0.2e1 * t623) * MDP(23) + (t389 * t540 + t496 * t581 + t484 + t527) * MDP(24) + (-t338 * qJ(5) - t339 * pkin(4) - t382 * t389 - t375 * t381 - g(1) * (-pkin(4) * t605 + t459) - g(2) * (-pkin(4) * t606 + t458) - g(3) * t592 - t581 * t379) * MDP(25) + t650 * MDP(28) + t533 * MDP(29) + (qJ(5) * t352 + t582 * t407 + (-t544 + t642) * t503 + t514 * t500 + t596) * MDP(31) + (qJ(5) * t351 + t334 + t582 * t409 + (t544 - t654) * t500 + t514 * t503) * MDP(32) + ((t375 - t581) * MDP(22) + t389 * MDP(23) + t409 * MDP(28) - t407 * MDP(29) - MDP(16) * t418 + t630) * t418 + t656; t523 * MDP(22) + (-t418 * t540 + t491) * MDP(23) + (-t496 ^ 2 - t645) * MDP(24) + (t379 * t496 + t519 - t623) * MDP(25) + (-t407 * t496 - t368) * MDP(31) + (-t409 * t496 + t616) * MDP(32) + (-MDP(31) * t652 - MDP(32) * t651) * t635; t409 * t407 * MDP(26) + (-t407 ^ 2 + t409 ^ 2) * MDP(27) + (t576 + t653) * MDP(28) + (t409 * t635 + t373) * MDP(29) - t374 * MDP(30) + (-g(1) * t439 - g(2) * t441 + t342 * t635 - t355 * t409 + t332) * MDP(31) + (g(1) * t440 - g(2) * t442 + t341 * t635 + t355 * t407) * MDP(32) + (-MDP(29) * t586 + MDP(31) * t547 - MDP(32) * t560) * t503 + (-MDP(28) * t586 + (-qJD(6) * t418 - t491) * MDP(29) - t560 * MDP(31) + (-t335 - t547) * MDP(32)) * t500;];
tau  = t1;
