% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP13_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP13_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP13_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:01:46
% EndTime: 2019-03-09 13:02:01
% DurationCPUTime: 8.72s
% Computational Cost: add. (6169->563), mult. (15633->752), div. (0->0), fcn. (11400->8), ass. (0->239)
t537 = cos(qJ(2));
t530 = sin(pkin(6));
t534 = sin(qJ(2));
t628 = qJD(1) * t534;
t608 = t530 * t628;
t531 = cos(pkin(6));
t629 = qJD(1) * t531;
t613 = pkin(1) * t629;
t472 = pkin(8) * t608 - t537 * t613;
t615 = qJD(3) + t472;
t533 = sin(qJ(4));
t536 = cos(qJ(4));
t577 = pkin(4) * t536 + pkin(10) * t533;
t678 = (-pkin(3) - t577) * t608 - qJD(4) * t577 - t615;
t532 = sin(qJ(5));
t630 = qJD(1) * t530;
t535 = cos(qJ(5));
t645 = t534 * t535;
t461 = (t532 * t537 + t533 * t645) * t630;
t623 = qJD(4) * t533;
t574 = -t535 * t623 - t461;
t519 = qJD(2) + t629;
t614 = qJD(1) * qJD(2);
t595 = t530 * t614;
t580 = t534 * t595;
t612 = pkin(1) * qJD(2) * t531;
t584 = qJD(1) * t612;
t633 = -pkin(8) * t580 + t537 * t584;
t440 = -t519 * qJD(3) - t633;
t408 = -pkin(3) * t580 - t440;
t538 = -pkin(2) - pkin(9);
t572 = qJD(4) + t608;
t670 = qJD(4) * t572;
t677 = -t538 * t670 + t408;
t528 = t534 ^ 2;
t676 = MDP(5) * (-t537 ^ 2 + t528);
t607 = t537 * t630;
t466 = t519 * t536 - t533 * t607;
t562 = t535 * t572;
t411 = t466 * t532 - t562;
t675 = t411 * t572;
t674 = t537 * t538;
t648 = t530 * t534;
t521 = pkin(8) * t648;
t609 = -pkin(1) * t537 - pkin(2);
t433 = pkin(3) * t648 + t521 + (-pkin(9) + t609) * t531;
t664 = qJ(3) * t534;
t594 = -pkin(1) - t664;
t450 = (t594 + t674) * t530;
t637 = t533 * t433 + t536 * t450;
t387 = pkin(10) * t648 + t637;
t667 = pkin(1) * t534;
t524 = t531 * t667;
t647 = t530 * t537;
t672 = pkin(8) * t647 + t524;
t468 = -t531 * qJ(3) - t672;
t449 = pkin(3) * t647 - t468;
t479 = t531 * t533 + t536 * t647;
t480 = t531 * t536 - t533 * t647;
t394 = pkin(4) * t479 - pkin(10) * t480 + t449;
t639 = t535 * t387 + t532 * t394;
t673 = t678 * t535;
t492 = pkin(4) * t533 - pkin(10) * t536 + qJ(3);
t632 = t535 * t533 * t538 + t532 * t492;
t517 = pkin(2) * t608;
t573 = pkin(9) * t534 - qJ(3) * t537;
t452 = t573 * t630 + t517;
t473 = pkin(8) * t607 + t534 * t613;
t454 = pkin(3) * t607 + t473;
t636 = t536 * t452 + t533 * t454;
t390 = pkin(10) * t607 + t636;
t622 = qJD(4) * t536;
t600 = t535 * t622;
t619 = qJD(5) * t535;
t671 = t535 * t390 - t492 * t619 + t678 * t532 - t538 * t600;
t616 = pkin(3) * t608 + t615;
t413 = t535 * t466 + t532 * t572;
t507 = t537 * t595;
t621 = qJD(4) * t537;
t599 = t536 * t621;
t626 = qJD(2) * t534;
t604 = t533 * t626;
t550 = -t599 + t604;
t546 = t550 * t530;
t602 = t519 * t623;
t543 = -qJD(1) * t546 + t602;
t374 = qJD(5) * t413 - t535 * t507 - t532 * t543;
t669 = t413 ^ 2;
t668 = pkin(3) + pkin(8);
t666 = pkin(5) * t532;
t665 = -qJ(6) - pkin(10);
t663 = qJ(6) * t532;
t404 = t519 * t538 + t616;
t498 = pkin(2) * t580;
t624 = qJD(3) * t534;
t544 = (qJD(2) * t573 - t624) * t530;
t407 = qJD(1) * t544 + t498;
t429 = qJD(1) * t450;
t467 = pkin(8) * t507 + t534 * t584;
t434 = pkin(3) * t507 + t467;
t582 = -t404 * t623 - t533 * t407 - t429 * t622 + t536 * t434;
t356 = -pkin(4) * t507 - t582;
t662 = t356 * t532;
t661 = t356 * t535;
t620 = qJD(5) * t532;
t373 = -qJD(5) * t562 + t466 * t620 - t532 * t507 + t535 * t543;
t660 = t373 * t532;
t659 = t374 * t535;
t650 = t519 * t533;
t464 = t536 * t607 + t650;
t463 = qJD(5) + t464;
t658 = t411 * t463;
t657 = t411 * t532;
t656 = t411 * t535;
t655 = t413 * t463;
t654 = t413 * t532;
t653 = t413 * t535;
t603 = t530 * t621;
t579 = qJD(1) * t603;
t418 = t519 * t622 - t533 * t579 - t536 * t580;
t652 = t418 * t532;
t651 = t418 * t535;
t508 = t519 * qJ(3);
t527 = t530 ^ 2;
t649 = t527 * qJD(1) ^ 2;
t646 = t532 * t538;
t644 = t535 * t536;
t378 = t533 * t404 + t536 * t429;
t369 = pkin(10) * t572 + t378;
t423 = t508 + t454;
t381 = pkin(4) * t464 - pkin(10) * t466 + t423;
t352 = -t369 * t532 + t535 * t381;
t348 = -qJ(6) * t413 + t352;
t347 = pkin(5) * t463 + t348;
t643 = t347 - t348;
t377 = t536 * t404 - t533 * t429;
t403 = pkin(4) * t466 + pkin(10) * t464;
t642 = t535 * t377 + t532 * t403;
t581 = t536 * t608;
t593 = pkin(5) - t646;
t617 = qJD(6) * t535;
t641 = pkin(5) * t581 + t390 * t532 - t632 * qJD(5) + (qJD(4) * t593 - t617) * t536 - t673 + (t536 * t620 - t574) * qJ(6);
t460 = t532 * t533 * t608 - t535 * t607;
t618 = qJD(5) * t536;
t597 = t535 * t618;
t640 = qJ(6) * t460 - qJ(6) * t597 + (-qJD(6) * t536 + (qJ(6) * qJD(4) - qJD(5) * t538) * t533) * t532 - t671;
t592 = qJD(5) * t665;
t635 = -t464 * t663 + t532 * t592 + t617 - t642;
t398 = t535 * t403;
t634 = -pkin(5) * t466 - t398 + (-qJ(6) * t464 + t592) * t535 + (-qJD(6) + t377) * t532;
t627 = qJD(2) * t533;
t625 = qJD(2) * t537;
t611 = t537 * t649;
t606 = t530 * t626;
t605 = t530 * t625;
t598 = t463 * t620;
t596 = t527 * t614;
t591 = -t387 * t532 + t535 * t394;
t590 = t433 * t536 - t533 * t450;
t589 = -t533 * t452 + t536 * t454;
t588 = t463 * t535;
t586 = t464 - t627;
t583 = t534 * t611;
t578 = -0.2e1 * pkin(1) * t596;
t576 = t473 * t519 - t467;
t575 = t532 * t623 + t460;
t353 = t369 * t535 + t381 * t532;
t349 = -qJ(6) * t411 + t353;
t571 = t347 * t535 + t349 * t532;
t570 = t347 * t532 - t349 * t535;
t474 = t672 * qJD(2);
t569 = t467 * t531 + t474 * t519;
t516 = t537 * t612;
t568 = -pkin(8) * t606 + t516;
t512 = pkin(2) * t606;
t425 = t512 + t544;
t455 = (t647 * t668 + t524) * qJD(2);
t567 = -t533 * t425 - t433 * t623 - t450 * t622 + t455 * t536;
t566 = -t519 * t607 + t507;
t565 = -t480 * t532 + t530 * t645;
t447 = t480 * t535 + t532 * t648;
t564 = -t463 * t619 - t652;
t563 = -t598 + t651;
t561 = t466 * t572;
t560 = t572 * t413;
t559 = t572 * t530;
t558 = t572 * t536;
t368 = -pkin(4) * t572 - t377;
t554 = -pkin(10) * t418 + t368 * t463;
t553 = -t404 * t622 - t536 * t407 + t429 * t623 - t533 * t434;
t552 = t536 * t425 + t433 * t622 - t450 * t623 + t533 * t455;
t355 = pkin(10) * t507 - t553;
t367 = pkin(10) * t602 + t418 * pkin(4) + (-pkin(3) * t626 - pkin(10) * t550) * t630 - t440;
t344 = t535 * t355 + t532 * t367 - t369 * t620 + t381 * t619;
t362 = pkin(10) * t605 + t552;
t525 = t531 * qJD(3);
t432 = -t606 * t668 + t516 + t525;
t444 = -qJD(4) * t479 + t530 * t604;
t445 = t531 * t622 - t533 * t603 - t536 * t606;
t372 = pkin(4) * t445 - pkin(10) * t444 + t432;
t551 = t535 * t362 + t532 * t372 - t387 * t620 + t394 * t619;
t386 = -pkin(4) * t648 - t590;
t469 = (-pkin(2) * t537 + t594) * t530;
t549 = qJD(1) * t559;
t548 = (-qJ(3) * t625 - t624) * t530;
t547 = pkin(4) * t607 + t589;
t363 = -pkin(4) * t605 - t567;
t345 = -qJD(5) * t353 - t355 * t532 + t535 * t367;
t545 = -qJD(5) * t639 - t362 * t532 + t535 * t372;
t346 = pkin(5) * t374 + t356;
t542 = -t534 * t549 - t670;
t540 = t536 * t507 + t533 * t542;
t505 = t665 * t535;
t504 = t665 * t532;
t484 = t535 * t492;
t471 = -qJ(3) * t607 + t517;
t470 = t531 * t609 + t521;
t462 = -t525 - t568;
t459 = t519 * t532 - t535 * t581;
t458 = t519 * t535 + t532 * t581;
t457 = qJD(1) * t469;
t456 = t512 + t548;
t451 = -t508 - t473;
t448 = -pkin(2) * t519 + t615;
t441 = -t536 * t663 + t632;
t437 = qJD(1) * t548 + t498;
t430 = t457 * t608;
t426 = -qJ(6) * t644 + t533 * t593 + t484;
t410 = t411 ^ 2;
t385 = qJD(5) * t565 + t444 * t535 + t532 * t605;
t384 = qJD(5) * t447 + t444 * t532 - t535 * t605;
t364 = t411 * pkin(5) + qJD(6) + t368;
t357 = qJ(6) * t565 + t639;
t350 = pkin(5) * t479 - qJ(6) * t447 + t591;
t343 = -qJ(6) * t384 + qJD(6) * t565 + t551;
t342 = pkin(5) * t445 - qJ(6) * t385 - qJD(6) * t447 + t545;
t341 = -qJ(6) * t374 - qJD(6) * t411 + t344;
t340 = pkin(5) * t418 + qJ(6) * t373 - qJD(6) * t413 + t345;
t1 = [(t527 * t628 + t559) * MDP(19) * t625 + (t345 * t479 + t352 * t445 - t356 * t565 + t363 * t411 + t368 * t384 + t386 * t374 + t418 * t591 + t463 * t545) * MDP(27) + (-t373 * t565 - t374 * t447 - t384 * t413 - t385 * t411) * MDP(23) + (-t340 * t447 + t341 * t565 - t342 * t413 - t343 * t411 - t347 * t385 - t349 * t384 + t350 * t373 - t357 * t374) * MDP(29) + (t341 * t357 + t349 * t343 + t340 * t350 + t347 * t342 + t346 * (-pkin(5) * t565 + t386) + t364 * (pkin(5) * t384 + t363)) * MDP(30) + (-t374 * t479 - t384 * t463 - t411 * t445 + t418 * t565) * MDP(25) + (t466 * t444 - t480 * t543) * MDP(15) + (-t480 * t418 - t444 * t464 - t466 * t445 + t479 * t543) * MDP(16) + (-t440 * t537 + t467 * t534 + (t448 * t537 + t451 * t534) * qJD(2) + (-t462 * t537 + t474 * t534 + (t468 * t534 + t470 * t537) * qJD(2)) * qJD(1)) * t530 * MDP(11) + (t444 * qJD(4) + ((qJD(1) * t480 + t466) * t625 + (-t602 + (t444 + t546) * qJD(1)) * t534) * t530) * MDP(17) + (t567 * t572 + t582 * t648 + t432 * t464 + t449 * t418 + t408 * t479 + t423 * t445 + (qJD(1) * t590 + t377) * t605) * MDP(20) + (-t552 * t572 + t553 * t648 + t432 * t466 + t449 * (-t536 * t579 - t602) + t408 * t480 + t423 * t444 + (-t378 * t537 + (t449 * t533 * t534 - t537 * t637) * qJD(1)) * t530 * qJD(2)) * MDP(21) + (-t344 * t479 - t353 * t445 + t356 * t447 + t363 * t413 + t368 * t385 - t386 * t373 - t418 * t639 - t463 * t551) * MDP(28) + (-t519 * t568 - t531 * t633 + t537 * t578) * MDP(10) + ((-t457 * t626 + t437 * t537 + (t456 * t537 - t469 * t626) * qJD(1)) * t530 + t569) * MDP(12) + (-t445 * qJD(4) + (-t464 * t625 - t418 * t534 + (-t445 * t534 - t479 * t625) * qJD(1)) * t530) * MDP(18) + (-t440 * t531 - t462 * t519 + (-t457 * t625 - t437 * t534 + (-t456 * t534 - t469 * t625) * qJD(1)) * t530) * MDP(13) + (t534 * t578 - t569) * MDP(9) + (-t373 * t447 + t385 * t413) * MDP(22) + (t437 * t469 + t440 * t468 + t448 * t474 + t451 * t462 + t456 * t457 + t467 * t470) * MDP(14) + (t418 * t479 + t445 * t463) * MDP(26) + (-t373 * t479 + t385 * t463 + t413 * t445 + t418 * t447) * MDP(24) + (MDP(6) * t605 - MDP(7) * t606) * (t519 + t629) + 0.2e1 * (t534 * t537 * MDP(4) - t676) * t596; ((-t418 - t561) * t536 + ((t464 + t650) * qJD(4) + (t534 * t586 + t599) * t630) * t533) * MDP(16) + (-t536 ^ 2 * t579 + ((-qJD(4) * t519 + t580) * t536 - t561) * t533) * MDP(15) + (t636 * t572 + t378 * t607 + t616 * t466 + (-qJ(3) * t579 + t677) * t536 + ((-t423 - t508) * qJD(4) + (-t423 * t534 + (t664 - t674) * qJD(2)) * t630) * t533) * MDP(21) + (qJ(3) * t418 - t589 * t572 + (qJD(2) * t536 * t538 - t377) * t607 + t677 * t533 + t616 * t464 + t423 * t558) * MDP(20) - t537 * MDP(19) * t549 + (-t368 * t460 + t547 * t411 + t484 * t418 + ((-qJD(5) * t492 + t390) * t532 - t673) * t463 + (-t368 * t532 * qJD(4) + t345 + (qJD(4) * t411 + t564) * t538) * t533 + (t352 * t608 + t368 * t619 + t662 - t538 * t374 + (-t463 * t646 + t352) * qJD(4)) * t536) * MDP(27) + (-t374 * t533 + t575 * t463 + (t564 - t675) * t536) * MDP(25) + (-t632 * t418 + t547 * t413 - t368 * t461 + t671 * t463 + (t538 * t598 - t344 + (-t368 * t535 + t413 * t538) * qJD(4)) * t533 + (-t353 * t572 - t368 * t620 + t538 * t373 + t661) * t536) * MDP(28) + ((-qJ(3) * qJD(2) - t451 - t473) * t534 + (-pkin(2) * qJD(2) - t448 + t615) * t537) * MDP(11) * t630 + (t649 * t667 + t576) * MDP(9) + (t340 * t426 + t341 * t441 + t346 * (-t538 + t666) * t536 + (t538 * t623 + (-t575 + t597) * pkin(5) + t547) * t364 + t640 * t349 + t641 * t347) * MDP(30) + (t411 * t461 + t413 * t460 + (t654 + t656) * t623 + (t660 - t659 + (-t653 + t657) * qJD(5)) * t536) * MDP(23) + (-t373 * t644 + (-t532 * t618 + t574) * t413) * MDP(22) + (t347 * t461 + t349 * t460 + t373 * t426 - t374 * t441 - t641 * t413 - t640 * t411 + t571 * t623 + (qJD(5) * t570 - t340 * t535 - t341 * t532) * t536) * MDP(29) + (pkin(1) * t611 - t472 * t519 - t633) * MDP(10) + (t615 * t519 + (t457 * t537 + t471 * t534) * t630 - t440) * MDP(13) + (-pkin(2) * t467 - qJ(3) * t440 - t448 * t473 - t451 * t615 - t457 * t471) * MDP(14) + (-t466 * t607 + t540) * MDP(17) + (-t471 * t607 + t430 - t576) * MDP(12) + (t536 * t542 + t586 * t607) * MDP(18) + t649 * t676 + (-t373 * t533 + t574 * t463 + (t560 + t563) * t536) * MDP(24) + (t418 * t533 + t463 * t558) * MDP(26) + (-qJD(2) + t519) * MDP(7) * t608 + t566 * MDP(6) - MDP(4) * t583; t566 * MDP(11) + MDP(12) * t583 + (-t519 ^ 2 - t528 * t649) * MDP(13) + (t451 * t519 + t430 + t467) * MDP(14) + (-t519 * t464 + t540) * MDP(20) + (-t533 * t507 - t519 * t466 + (-t572 * t608 - t670) * t536) * MDP(21) + (-t374 * t536 + (-t532 * t622 - t458) * t463 + (t564 + t675) * t533) * MDP(27) + (t373 * t536 + (t459 - t600) * t463 + (t560 - t563) * t533) * MDP(28) + (t411 * t459 + t413 * t458 + (t654 - t656) * t622 + (-t660 - t659 + (t653 + t657) * qJD(5)) * t533) * MDP(29) + (-t347 * t458 - t349 * t459 + (-qJD(4) * t570 - t346) * t536 + (-qJD(5) * t571 - t340 * t532 + t341 * t535 + t364 * t572) * t533) * MDP(30); -t464 ^ 2 * MDP(16) + ((t464 - t650) * qJD(4) + (-t599 + (t464 + t627) * t534) * t630) * MDP(17) - t418 * MDP(18) + MDP(19) * t507 + (t378 * t572 + t582) * MDP(20) + (t377 * t572 + t423 * t464 + t553) * MDP(21) + (t413 * t588 - t660) * MDP(22) + ((-t373 - t658) * t535 + (-t374 - t655) * t532) * MDP(23) + (t463 * t588 + t652) * MDP(24) + (-t463 ^ 2 * t532 + t651) * MDP(25) + (-pkin(4) * t374 - t661 - t378 * t411 + (-pkin(10) * t619 - t398) * t463 + (t377 * t463 + t554) * t532) * MDP(27) + (pkin(4) * t373 + t662 - t378 * t413 + (pkin(10) * t620 + t642) * t463 + t554 * t535) * MDP(28) + (t373 * t504 + t374 * t505 - t634 * t413 - t635 * t411 + (-t347 * t463 + t341) * t535 + (-t349 * t463 - t340) * t532) * MDP(29) + (-t341 * t505 + t340 * t504 + t346 * (-pkin(5) * t535 - pkin(4)) + (t463 * t666 - t378) * t364 + t635 * t349 + t634 * t347) * MDP(30) + (MDP(15) * t464 + t466 * MDP(16) + MDP(18) * t572 - t423 * MDP(20) - t413 * MDP(24) + t411 * MDP(25) - t463 * MDP(26) - t352 * MDP(27) + t353 * MDP(28)) * t466; t413 * t411 * MDP(22) + (-t410 + t669) * MDP(23) + (-t373 + t658) * MDP(24) + (-t374 + t655) * MDP(25) + t418 * MDP(26) + (t353 * t463 - t368 * t413 + t345) * MDP(27) + (t352 * t463 + t368 * t411 - t344) * MDP(28) + (pkin(5) * t373 - t411 * t643) * MDP(29) + (t643 * t349 + (-t364 * t413 + t340) * pkin(5)) * MDP(30); (-t410 - t669) * MDP(29) + (t347 * t413 + t349 * t411 + t346) * MDP(30);];
tauc  = t1;
