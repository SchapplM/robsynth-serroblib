% Calculate vector of inverse dynamics joint torques for
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RRPPPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:08:32
% EndTime: 2019-03-09 08:08:42
% DurationCPUTime: 8.75s
% Computational Cost: add. (5207->557), mult. (12211->694), div. (0->0), fcn. (9048->12), ass. (0->240)
t535 = sin(pkin(9));
t539 = sin(qJ(2));
t542 = cos(qJ(2));
t648 = cos(pkin(9));
t495 = t535 * t542 + t539 * t648;
t534 = sin(pkin(10));
t536 = cos(pkin(10));
t538 = sin(qJ(6));
t541 = cos(qJ(6));
t667 = -t534 * t541 + t536 * t538;
t422 = t667 * t495;
t479 = t495 * qJD(1);
t455 = t536 * qJD(2) - t479 * t534;
t569 = qJD(2) * t534 + t536 * t479;
t639 = t569 * t538;
t395 = t455 * t541 + t639;
t599 = t648 * t542;
t513 = qJD(1) * t599;
t623 = qJD(1) * t539;
t476 = t535 * t623 - t513;
t616 = qJD(6) - t476;
t678 = t395 * t616;
t397 = -t455 * t538 + t541 * t569;
t677 = t397 * t616;
t676 = t455 * t476;
t619 = qJD(6) * t541;
t620 = qJD(6) * t538;
t626 = t667 * t476 + t534 * t619 - t536 * t620;
t531 = qJ(2) + pkin(9);
t528 = cos(t531);
t519 = g(3) * t528;
t527 = sin(t531);
t540 = sin(qJ(1));
t543 = cos(qJ(1));
t589 = g(1) * t543 + g(2) * t540;
t672 = t589 * t527;
t675 = t672 - t519;
t666 = MDP(15) + MDP(18);
t674 = t569 ^ 2;
t529 = t542 * pkin(2);
t524 = t529 + pkin(1);
t478 = t495 * qJD(2);
t614 = qJDD(1) * t539;
t585 = -qJDD(1) * t599 + t535 * t614;
t438 = qJD(1) * t478 + t585;
t615 = qJD(1) * qJD(2);
t605 = t539 * t615;
t439 = qJD(2) * t513 + qJDD(1) * t495 - t535 * t605;
t612 = pkin(2) * t605 + qJDD(3);
t554 = -qJDD(1) * t524 + t612;
t366 = pkin(3) * t438 - qJ(4) * t439 - qJD(4) * t479 + t554;
t650 = qJ(3) + pkin(7);
t601 = qJD(2) * t650;
t474 = -qJD(3) * t539 - t542 * t601;
t503 = t650 * t539;
t431 = qJDD(2) * pkin(2) + qJD(1) * t474 - qJDD(1) * t503;
t473 = qJD(3) * t542 - t539 * t601;
t504 = t650 * t542;
t441 = qJD(1) * t473 + qJDD(1) * t504;
t389 = t535 * t431 + t648 * t441;
t383 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t389;
t348 = t366 * t536 - t534 * t383;
t584 = qJDD(5) - t348;
t346 = -pkin(4) * t438 + t584;
t502 = -qJD(1) * t524 + qJD(3);
t406 = pkin(3) * t476 - qJ(4) * t479 + t502;
t498 = qJD(1) * t503;
t649 = qJD(2) * pkin(2);
t491 = -t498 + t649;
t499 = qJD(1) * t504;
t600 = t648 * t499;
t437 = t535 * t491 + t600;
t427 = qJD(2) * qJ(4) + t437;
t376 = t534 * t406 + t536 * t427;
t365 = t476 * qJ(5) + t376;
t671 = -t365 * t476 + t346;
t418 = pkin(2) * t623 + pkin(3) * t479 + qJ(4) * t476;
t484 = t535 * t499;
t443 = -t498 * t648 - t484;
t385 = t534 * t418 + t536 * t443;
t371 = t479 * qJ(5) + t385;
t622 = qJD(4) * t536;
t670 = -t371 + t622;
t669 = -t528 * pkin(3) - t527 * qJ(4);
t657 = g(1) * t540;
t604 = -g(2) * t543 + t657;
t413 = qJDD(2) * t534 + t439 * t536;
t621 = qJD(5) * t569;
t668 = qJ(5) * t413 + t621;
t432 = -qJDD(6) + t438;
t494 = t534 * t538 + t536 * t541;
t627 = t616 * t494;
t664 = t432 * t667 - t616 * t627;
t646 = qJ(5) * t536;
t661 = pkin(4) + pkin(5);
t663 = t534 * t661 - t646;
t412 = -t536 * qJDD(2) + t439 * t534;
t596 = -t541 * t412 + t538 * t413;
t354 = t397 * qJD(6) + t596;
t388 = t648 * t431 - t535 * t441;
t565 = qJDD(2) * pkin(3) - qJDD(4) + t388;
t548 = t565 + t668;
t345 = -t412 * t661 + t548;
t662 = t345 + t675;
t475 = t476 ^ 2;
t660 = pkin(2) * t539;
t659 = pkin(4) * t412;
t658 = pkin(8) * t534;
t654 = g(3) * t527;
t653 = g(3) * t542;
t652 = t536 * pkin(4);
t517 = pkin(2) * t535 + qJ(4);
t651 = -pkin(8) + t517;
t645 = qJDD(1) * pkin(1);
t643 = t395 * t479;
t642 = t397 * t479;
t640 = t438 * t534;
t638 = t517 * t536;
t637 = t527 * t543;
t636 = t528 * t543;
t635 = t534 * qJ(5);
t634 = t534 * t540;
t631 = t536 * t543;
t630 = t650 * t543;
t629 = t540 * t536;
t628 = t543 * t534;
t349 = t534 * t366 + t536 * t383;
t558 = -t535 * t539 + t599;
t481 = t558 * qJD(2);
t611 = t539 * t649;
t398 = pkin(3) * t478 - qJ(4) * t481 - qJD(4) * t495 + t611;
t420 = t473 * t648 + t535 * t474;
t368 = t534 * t398 + t536 * t420;
t435 = -pkin(3) * t558 - qJ(4) * t495 - t524;
t450 = -t535 * t503 + t504 * t648;
t392 = t534 * t435 + t536 * t450;
t436 = t648 * t491 - t484;
t625 = (g(1) * t631 + g(2) * t629) * t527;
t532 = t539 ^ 2;
t624 = -t542 ^ 2 + t532;
t566 = qJD(2) * pkin(3) - qJD(4) + t436;
t618 = -qJD(4) - t566;
t552 = qJ(5) * t569 + t566;
t374 = -pkin(4) * t455 - t552;
t617 = qJD(4) - t374;
t613 = qJDD(1) * t542;
t610 = t438 * t638;
t609 = t538 * t412 + t541 * t413 - t455 * t619;
t381 = -qJ(5) * t558 + t392;
t511 = t543 * t524;
t608 = pkin(3) * t636 + qJ(4) * t637 + t511;
t607 = t529 - t669;
t347 = -t548 + t659;
t603 = -t347 - t519;
t602 = t565 - t519;
t340 = -pkin(8) * t413 - t438 * t661 + t584;
t344 = t438 * qJ(5) + t476 * qJD(5) + t349;
t343 = pkin(8) * t412 + t344;
t598 = t541 * t340 - t538 * t343;
t408 = t534 * t420;
t367 = t398 * t536 - t408;
t375 = t406 * t536 - t534 * t427;
t433 = t534 * t443;
t384 = t418 * t536 - t433;
t445 = t534 * t450;
t391 = t435 * t536 - t445;
t419 = t473 * t535 - t648 * t474;
t442 = -t498 * t535 + t600;
t449 = t648 * t503 + t504 * t535;
t595 = t616 ^ 2;
t594 = qJD(5) * t534 - t476 * t663 + t442;
t355 = t478 * qJ(5) - qJD(5) * t558 + t368;
t593 = -g(2) * t637 + t527 * t657;
t523 = -pkin(2) * t648 - pkin(3);
t592 = -pkin(3) * t527 - t660;
t464 = t528 * t634 + t631;
t466 = t528 * t628 - t629;
t591 = -g(1) * t464 + g(2) * t466;
t465 = t528 * t629 - t628;
t467 = t528 * t631 + t634;
t590 = g(1) * t465 - g(2) * t467;
t588 = t494 * t432 - t616 * t626;
t587 = qJD(5) - t375;
t586 = pkin(4) * t534 - t646;
t583 = t538 * t340 + t541 * t343;
t582 = -t344 * t534 + t346 * t536;
t581 = t347 * t495 + t374 * t481;
t580 = -t348 * t536 - t349 * t534;
t352 = -pkin(8) * t569 - t476 * t661 + t587;
t357 = -pkin(8) * t455 + t365;
t341 = t352 * t541 - t357 * t538;
t342 = t352 * t538 + t357 * t541;
t362 = t445 + (-pkin(8) * t495 - t435) * t536 + t661 * t558;
t369 = t495 * t658 + t381;
t579 = t362 * t541 - t369 * t538;
t578 = t362 * t538 + t369 * t541;
t577 = -t375 * t534 + t376 * t536;
t576 = -t481 * t566 - t495 * t565;
t575 = t438 * t536 + t455 * t479;
t574 = t479 * t569 + t640;
t571 = t464 * t541 - t465 * t538;
t570 = t464 * t538 + t465 * t541;
t568 = qJD(4) * t569 + t413 * t517;
t562 = -0.2e1 * pkin(1) * t615 - pkin(7) * qJDD(2);
t489 = t651 * t536;
t560 = qJD(4) * t534 - qJD(6) * t489 - t433 - (pkin(8) * t476 - t418) * t536 + t661 * t479;
t488 = t651 * t534;
t559 = qJD(6) * t488 + t476 * t658 + t670;
t423 = t494 * t495;
t557 = t569 * t620 - t609;
t556 = t523 - t635;
t555 = -qJD(5) * t495 * t536 + t419;
t544 = qJD(2) ^ 2;
t551 = -pkin(7) * t544 + t604 + 0.2e1 * t645;
t545 = qJD(1) ^ 2;
t550 = pkin(1) * t545 - pkin(7) * qJDD(1) + t589;
t549 = -t412 * t638 + t455 * t622 - t528 * t589 - t654;
t547 = (-g(1) * (-t524 + t669) - g(2) * t650) * t540;
t546 = -t565 - t675;
t507 = qJ(4) * t636;
t505 = t540 * t528 * qJ(4);
t487 = t556 - t652;
t459 = t536 * t661 - t556;
t415 = t466 * t538 + t467 * t541;
t414 = t466 * t541 - t467 * t538;
t393 = t495 * t586 + t449;
t390 = -t476 * t586 + t442;
t387 = -t495 * t663 - t449;
t382 = pkin(4) * t558 - t391;
t380 = qJD(6) * t423 + t481 * t667;
t379 = -qJD(6) * t422 + t481 * t494;
t372 = -pkin(4) * t479 - t384;
t370 = t481 * t586 + t555;
t364 = -pkin(4) * t476 + t587;
t360 = t455 * t661 + t552;
t359 = -t481 * t663 - t555;
t358 = -pkin(4) * t478 - t367;
t351 = t481 * t658 + t355;
t350 = t408 + (-pkin(8) * t481 - t398) * t536 - t661 * t478;
t1 = [0.2e1 * (t539 * t613 - t615 * t624) * MDP(5) + t604 * MDP(2) + (qJDD(1) * t532 + 0.2e1 * t542 * t605) * MDP(4) + t589 * MDP(3) + (t539 * t562 + t542 * t551) * MDP(9) + (-t539 * t551 + t542 * t562) * MDP(10) + (-t354 * t423 - t379 * t395 - t380 * t397 + t422 * t557) * MDP(22) + (t379 * t397 - t423 * t557) * MDP(21) + (-g(1) * t630 - g(2) * t608 + t348 * t391 + t349 * t392 + t375 * t367 + t376 * t368 - t419 * t566 - t449 * t565 + t547) * MDP(16) + (-t367 * t569 + t368 * t455 - t391 * t413 - t392 * t412 + t580 * t495 + (-t375 * t536 - t376 * t534) * t481 + t593) * MDP(15) + (t355 * t455 + t358 * t569 - t381 * t412 + t382 * t413 + t582 * t495 + (t364 * t536 - t365 * t534) * t481 + t593) * MDP(18) + (-t348 * t558 + t367 * t476 + t375 * t478 + t391 * t438 + t412 * t449 - t419 * t455 + t534 * t576 + t590) * MDP(13) + (t346 * t558 - t358 * t476 - t364 * t478 - t370 * t455 - t382 * t438 + t393 * t412 + t534 * t581 + t590) * MDP(17) + (-(t350 * t538 + t351 * t541) * t616 + t578 * t432 - t583 * t558 + t342 * t478 + t359 * t397 - t387 * t557 + t345 * t423 + t360 * t379 + g(1) * t571 - g(2) * t414 + (-t341 * t558 - t579 * t616) * qJD(6)) * MDP(27) + (t379 * t616 - t397 * t478 - t423 * t432 - t557 * t558) * MDP(23) + ((t350 * t541 - t351 * t538) * t616 - t579 * t432 + t598 * t558 - t341 * t478 + t359 * t395 + t387 * t354 + t345 * t422 + t360 * t380 + g(1) * t570 - g(2) * t415 + (-t342 * t558 - t578 * t616) * qJD(6)) * MDP(26) + (-t354 * t558 - t380 * t616 + t395 * t478 + t422 * t432) * MDP(24) + (-t432 * t558 - t478 * t616) * MDP(25) + (-t388 * t495 + t389 * t558 + t419 * t479 - t420 * t476 - t436 * t481 - t437 * t478 - t438 * t450 + t439 * t449 - t589) * MDP(11) + (t349 * t558 - t368 * t476 - t376 * t478 - t392 * t438 + t413 * t449 + t419 * t569 + t536 * t576 + t591) * MDP(14) + (-t344 * t558 + t355 * t476 + t365 * t478 - t370 * t569 + t381 * t438 - t393 * t413 - t536 * t581 - t591) * MDP(19) + (t389 * t450 + t437 * t420 - t388 * t449 - t436 * t419 - t554 * t524 + t502 * t611 - g(1) * (-t524 * t540 + t630) - g(2) * (t540 * t650 + t511)) * MDP(12) + (t344 * t381 + t365 * t355 + t347 * t393 + t374 * t370 + t346 * t382 + t364 * t358 - g(1) * (-pkin(4) * t465 - qJ(5) * t464 + t630) - g(2) * (pkin(4) * t467 + qJ(5) * t466 + t608) + t547) * MDP(20) + (qJDD(2) * t539 + t542 * t544) * MDP(6) + (qJDD(2) * t542 - t539 * t544) * MDP(7) + qJDD(1) * MDP(1); t616 * t479 * MDP(25) + MDP(6) * t614 + MDP(7) * t613 + qJDD(2) * MDP(8) + (t539 * t550 - t653) * MDP(9) + (g(3) * t539 + t542 * t550) * MDP(10) + ((t437 - t442) * t479 + (-t436 + t443) * t476 + (-t438 * t535 - t439 * t648) * pkin(2)) * MDP(11) + (t436 * t442 - t437 * t443 + (t648 * t388 - t653 + t389 * t535 + (-qJD(1) * t502 + t589) * t539) * pkin(2)) * MDP(12) + (-t517 * t640 - t375 * t479 + t412 * t523 + t442 * t455 + t602 * t536 + (t534 * t618 - t384) * t476 + t625) * MDP(13) + (-t610 + t376 * t479 + t413 * t523 - t442 * t569 + (t536 * t618 + t385) * t476 + (-t672 - t602) * t534) * MDP(14) + (t384 * t569 - t385 * t455 + (-t375 * t476 + t349) * t536 + (-t376 * t476 - t348 + t568) * t534 + t549) * MDP(15) + (-t565 * t523 - t376 * t385 - t375 * t384 + t566 * t442 - g(1) * (t543 * t592 + t507) - g(2) * (t540 * t592 + t505) - g(3) * t607 + (-t348 * t534 + t349 * t536) * t517 + t577 * qJD(4)) * MDP(16) + (t364 * t479 + t372 * t476 + t390 * t455 + t412 * t487 + t603 * t536 + (qJD(5) * t455 - t438 * t517 - t476 * t617) * t534 + t625) * MDP(17) + (-t371 * t455 - t372 * t569 + (t364 * t476 + t344) * t536 + (t568 + t671) * t534 + t549) * MDP(18) + (t610 - t365 * t479 + t390 * t569 - t413 * t487 + (t536 * t617 - t371) * t476 + (t603 + t621 + t672) * t534) * MDP(19) + (t344 * t638 + t347 * t487 - t374 * t390 - t364 * t372 - g(1) * (-t543 * t660 + t507) - g(2) * (-t540 * t660 + t505) - g(3) * (t528 * t652 + t607) + t670 * t365 + (-qJ(5) * t519 + qJD(4) * t364 - qJD(5) * t374 + t346 * t517) * t534 + (pkin(3) + t635 + t652) * t672) * MDP(20) + (-t397 * t627 + t557 * t667) * MDP(21) + (t354 * t667 + t395 * t627 - t397 * t626 + t494 * t557) * MDP(22) + (t642 + t664) * MDP(23) + (t588 - t643) * MDP(24) + (-(t488 * t541 - t489 * t538) * t432 + t459 * t354 + t341 * t479 - (t538 * t559 - t541 * t560) * t616 + t594 * t395 + t626 * t360 + t662 * t494) * MDP(26) + ((t488 * t538 + t489 * t541) * t432 - t459 * t557 - t342 * t479 - (t538 * t560 + t541 * t559) * t616 + t594 * t397 - t627 * t360 - t662 * t667) * MDP(27) + (-MDP(4) * t539 * t542 + MDP(5) * t624) * t545; (-t479 ^ 2 - t475) * MDP(11) + (-pkin(2) * t613 + t436 * t479 - t604 + t612 - t645) * MDP(12) + t575 * MDP(13) + (-t475 * t536 - t574) * MDP(14) + (t479 * t566 - t580 - t604) * MDP(16) + (-t475 * t534 + t575) * MDP(17) + t574 * MDP(19) + (-t374 * t479 - t582 - t604) * MDP(20) + (t588 + t643) * MDP(26) + (t642 - t664) * MDP(27) + (t437 * MDP(12) + t577 * MDP(16) + (t364 * t534 + t365 * t536) * MDP(20) + (-MDP(13) * t534 + MDP(19) * t536) * t476 + t666 * (t455 * t536 + t534 * t569)) * t476 + t666 * (-t534 * t412 - t413 * t536); (t375 * t569 - t376 * t455 + t546) * MDP(16) + (-t364 * t569 - t365 * t455 + t546 + t659 - t668) * MDP(20) + (-t354 - t677) * MDP(26) + (t557 + t678) * MDP(27) + t666 * (-t455 ^ 2 - t674) + (MDP(13) + MDP(17)) * (t476 * t569 + t412) + (MDP(14) - MDP(19)) * (t413 + t676); (-qJD(2) * t479 - t455 * t569 - t585) * MDP(17) + (t413 - t676) * MDP(18) + (-t475 - t674) * MDP(19) + (-g(1) * t466 - g(2) * t464 + t374 * t569 - t534 * t654 + t671) * MDP(20) + (-t395 * t569 - t541 * t432 - t538 * t595) * MDP(26) + (-t397 * t569 + t538 * t432 - t541 * t595) * MDP(27); t397 * t395 * MDP(21) + (-t395 ^ 2 + t397 ^ 2) * MDP(22) + (t609 + t678) * MDP(23) + (-t596 + t677) * MDP(24) - t432 * MDP(25) + (-g(1) * t414 - g(2) * t571 + t342 * t616 - t360 * t397 + t598) * MDP(26) + (g(1) * t415 + g(2) * t570 + t341 * t616 + t360 * t395 - t583) * MDP(27) + (MDP(26) * t667 + MDP(27) * t494) * t654 + (-MDP(23) * t639 - MDP(24) * t397 - MDP(26) * t342 - MDP(27) * t341) * qJD(6);];
tau  = t1;
