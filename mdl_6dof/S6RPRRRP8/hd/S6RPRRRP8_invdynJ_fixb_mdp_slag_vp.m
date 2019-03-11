% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRRRP8_invdynJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:24:47
% EndTime: 2019-03-09 06:24:58
% DurationCPUTime: 7.92s
% Computational Cost: add. (7078->553), mult. (13442->676), div. (0->0), fcn. (8989->10), ass. (0->240)
t538 = cos(qJ(3));
t686 = cos(qJ(4));
t610 = t686 * qJD(4);
t564 = t686 * qJD(3) + t610;
t535 = sin(qJ(3));
t685 = sin(qJ(4));
t614 = t685 * t535;
t599 = qJD(1) * t614;
t616 = t686 * t535;
t626 = qJD(3) + qJD(4);
t408 = qJDD(1) * t616 + (qJD(1) * t564 + t685 * qJDD(1)) * t538 - t626 * t599;
t407 = qJDD(5) + t408;
t469 = t538 * t685 + t616;
t463 = t469 * qJD(1);
t698 = qJD(5) + t463;
t704 = t698 ^ 2;
t615 = t686 * t538;
t462 = -qJD(1) * t615 + t599;
t534 = sin(qJ(5));
t537 = cos(qJ(5));
t605 = t537 * t626;
t438 = -t462 * t534 - t605;
t703 = t438 * t698;
t569 = t537 * t462 - t534 * t626;
t702 = t569 * t698;
t635 = qJD(5) * t537;
t661 = t463 * t537;
t701 = t635 + t661;
t636 = qJD(5) * t534;
t700 = t463 * t534 + t636;
t539 = cos(qJ(1));
t527 = g(2) * t539;
t536 = sin(qJ(1));
t681 = g(1) * t536;
t608 = -t527 + t681;
t632 = qJD(1) * qJD(3);
t609 = t538 * t632;
t631 = qJDD(1) * t535;
t699 = t609 + t631;
t669 = t407 * t537;
t697 = pkin(9) * (t636 * t698 - t669);
t696 = t463 * t626;
t541 = -pkin(1) - pkin(7);
t491 = qJDD(1) * t541 + qJDD(2);
t470 = t538 * t491;
t493 = qJD(1) * t541 + qJD(2);
t630 = qJDD(1) * t538;
t638 = qJD(3) * t535;
t424 = -t493 * t638 + qJDD(3) * pkin(3) + t470 + (t535 * t632 - t630) * pkin(8);
t637 = qJD(3) * t538;
t427 = -pkin(8) * t699 + t491 * t535 + t493 * t637;
t639 = qJD(1) * t538;
t453 = -pkin(8) * t639 + t538 * t493;
t449 = qJD(3) * pkin(3) + t453;
t640 = qJD(1) * t535;
t452 = -pkin(8) * t640 + t493 * t535;
t612 = qJD(4) * t685;
t555 = t424 * t685 + t427 * t686 + t449 * t610 - t452 * t612;
t625 = qJDD(3) + qJDD(4);
t367 = pkin(9) * t625 + t555;
t528 = qJDD(1) * qJ(2);
t529 = qJD(1) * qJD(2);
t450 = pkin(3) * t699 + t528 + t529;
t468 = t614 - t615;
t545 = -t468 * qJDD(1) - t696;
t375 = t408 * pkin(4) - pkin(9) * t545 + t450;
t448 = t686 * t452;
t417 = t685 * t449 + t448;
t410 = pkin(9) * t626 + t417;
t480 = pkin(3) * t640 + qJD(1) * qJ(2);
t419 = pkin(4) * t463 + pkin(9) * t462 + t480;
t572 = t537 * t367 + t534 * t375 - t410 * t636 + t419 * t635;
t675 = qJ(6) * t407;
t353 = qJD(6) * t698 + t572 + t675;
t601 = -t534 * t367 + t537 * t375 - t410 * t635 - t419 * t636;
t684 = pkin(5) * t407;
t355 = qJDD(6) - t601 - t684;
t589 = t353 * t537 + t355 * t534;
t525 = t535 * pkin(3);
t506 = qJ(2) + t525;
t431 = pkin(4) * t469 + pkin(9) * t468 + t506;
t678 = pkin(8) - t541;
t477 = t678 * t535;
t478 = t678 * t538;
t437 = -t477 * t686 - t478 * t685;
t646 = t534 * t431 + t537 * t437;
t695 = pkin(5) * t700 - qJ(6) * t701 - t534 * qJD(6);
t533 = qJ(3) + qJ(4);
t522 = cos(t533);
t507 = t522 * pkin(9);
t694 = t507 - t525;
t543 = qJD(1) ^ 2;
t566 = -qJ(2) * t543 - t608;
t421 = t453 * t685 + t448;
t580 = pkin(3) * t612 - t421;
t692 = -qJD(3) * t685 - t612;
t691 = MDP(26) + MDP(28);
t627 = MDP(27) - MDP(30);
t594 = g(1) * t539 + g(2) * t536;
t622 = 0.2e1 * t529;
t690 = 0.2e1 * t528 + t622 - t594;
t388 = -qJD(5) * t569 + t534 * t545 - t537 * t625;
t688 = t569 ^ 2;
t683 = pkin(5) * t462;
t682 = pkin(9) * t407;
t521 = sin(t533);
t508 = g(3) * t521;
t509 = g(3) * t522;
t680 = g(3) * t534;
t679 = g(3) * t537;
t677 = pkin(9) * qJD(5);
t676 = pkin(1) * qJDD(1);
t447 = t685 * t452;
t416 = t449 * t686 - t447;
t409 = -pkin(4) * t626 - t416;
t381 = t438 * pkin(5) + qJ(6) * t569 + t409;
t674 = t381 * t463;
t386 = t410 * t537 + t419 * t534;
t673 = t386 * t698;
t387 = -qJD(5) * t605 - t462 * t636 - t534 * t625 - t537 * t545;
t672 = t387 * t534;
t671 = t388 * t537;
t513 = pkin(3) * t685 + pkin(9);
t670 = t407 * t513;
t668 = t409 * t463;
t667 = t438 * t534;
t666 = t438 * t537;
t665 = t569 * t438;
t664 = t569 * t534;
t663 = t569 * t537;
t660 = t468 * t534;
t659 = t468 * t537;
t658 = t521 * t536;
t657 = t522 * t536;
t656 = t522 * t539;
t655 = t534 * t536;
t654 = t534 * t539;
t653 = t536 * t537;
t652 = t539 * t537;
t650 = -t580 - t695;
t430 = -pkin(4) * t462 + pkin(9) * t463;
t649 = t537 * t416 + t534 * t430;
t422 = t453 * t686 - t447;
t425 = pkin(3) * t639 + t430;
t648 = t537 * t422 + t534 * t425;
t434 = -t535 * t564 + t538 * t692;
t647 = t434 * t626 - t468 * t625;
t645 = -t417 + t695;
t489 = g(2) * t656;
t644 = t537 * t489 + t521 * t679;
t643 = t539 * pkin(1) + t536 * qJ(2);
t532 = t538 ^ 2;
t642 = t535 ^ 2 - t532;
t542 = qJD(3) ^ 2;
t641 = -t542 - t543;
t494 = pkin(3) * t637 + qJD(2);
t385 = -t410 * t534 + t419 * t537;
t633 = qJD(6) - t385;
t629 = qJDD(3) * t535;
t624 = g(1) * t657;
t488 = t521 * t527;
t623 = t686 * pkin(3);
t621 = t522 * t655;
t620 = g(1) * t621 - t534 * t489 - t521 * t680;
t619 = g(1) * t658 - t488 + t509;
t618 = t534 * t686;
t617 = t537 * t686;
t606 = t537 * t698;
t604 = qJDD(2) - t676;
t603 = pkin(3) * t610;
t600 = -t686 * t424 + t685 * t427 + t449 * t612 + t452 * t610;
t368 = -pkin(4) * t625 + t600;
t356 = t388 * pkin(5) + t387 * qJ(6) + qJD(6) * t569 + t368;
t598 = -t356 - t624;
t597 = -t368 - t624;
t455 = t521 * t655 - t652;
t457 = t521 * t654 + t653;
t596 = g(1) * t457 + g(2) * t455;
t456 = t521 * t653 + t654;
t458 = t521 * t652 - t655;
t595 = -g(1) * t458 - g(2) * t456;
t592 = -t537 * pkin(5) - t534 * qJ(6);
t591 = pkin(5) * t534 - qJ(6) * t537;
t376 = -pkin(5) * t698 + t633;
t377 = qJ(6) * t698 + t386;
t588 = t376 * t537 - t377 * t534;
t587 = t376 * t534 + t377 * t537;
t586 = -t671 - t672;
t584 = t668 - t670;
t583 = -t663 + t667;
t579 = pkin(4) - t592;
t578 = -t376 * t462 + t381 * t636 + t644;
t577 = t385 * t462 + t409 * t636 + t644;
t575 = -t434 * t534 + t468 * t635;
t574 = t434 * t537 + t468 * t636;
t573 = 0.2e1 * qJ(2) * t632 + qJDD(3) * t541;
t435 = t535 * t692 + t564 * t538;
t396 = pkin(4) * t435 - pkin(9) * t434 + t494;
t436 = -t477 * t685 + t478 * t686;
t466 = t678 * t638;
t467 = qJD(3) * t478;
t398 = -qJD(4) * t436 + t466 * t685 - t467 * t686;
t571 = t534 * t396 + t537 * t398 + t431 * t635 - t437 * t636;
t570 = pkin(4) * t521 - t694;
t567 = t579 * t522;
t565 = t368 * t534 - t386 * t462 + t409 * t635 + t620;
t563 = -t356 * t534 + t377 * t462 - t381 * t661 - t620;
t561 = -t513 * t636 + t537 * t603;
t560 = g(1) * t455 - g(2) * t457 + t522 * t680 + t601;
t559 = t480 * t462 + t489 + t508 - t600 - t624;
t558 = qJD(5) * t588 + t589;
t557 = qJD(5) * t583 + t586;
t556 = -t541 * t542 + t690;
t554 = -t381 * t569 + qJDD(6) - t560;
t553 = -g(1) * (t522 * pkin(5) * t653 + pkin(4) * t657 + pkin(9) * t658 + qJ(6) * t621) + t579 * t508;
t552 = t376 * t701 - t377 * t700 + t589 - t619;
t551 = -g(1) * t456 + g(2) * t458 - t522 * t679 + t572;
t550 = MDP(29) * t583 + MDP(31) * t588;
t549 = ((-t387 - t703) * t537 + (-t388 + t702) * t534) * MDP(22) + (-t569 * t606 - t672) * MDP(21) + (-t438 * t462 - t534 * t704 + t669) * MDP(24) + (t407 * t534 - t462 * t569 + t606 * t698) * MDP(23) + (t545 + t696) * MDP(16) + (-t462 * t626 - t408) * MDP(17) + (t462 ^ 2 - t463 ^ 2) * MDP(15) + t625 * MDP(18) + (-MDP(14) * t463 + MDP(25) * t698) * t462;
t399 = qJD(4) * t437 - t466 * t686 - t467 * t685;
t547 = t480 * t463 - t555 + t619;
t540 = -pkin(8) - pkin(7);
t524 = t539 * qJ(2);
t520 = qJDD(3) * t538;
t514 = -t623 - pkin(4);
t465 = -t623 - t579;
t454 = t462 * qJ(6);
t402 = -pkin(5) * t569 + qJ(6) * t438;
t400 = -t468 * t591 + t436;
t391 = -pkin(5) * t469 - t431 * t537 + t437 * t534;
t390 = qJ(6) * t469 + t646;
t383 = t416 * t534 - t430 * t537 + t683;
t382 = -t454 + t649;
t380 = t422 * t534 - t425 * t537 + t683;
t379 = -t454 + t648;
t370 = -t387 + t703;
t360 = t591 * t434 + (qJD(5) * t592 + qJD(6) * t537) * t468 + t399;
t359 = -pkin(5) * t435 + qJD(5) * t646 - t396 * t537 + t398 * t534;
t358 = qJ(6) * t435 + qJD(6) * t469 + t571;
t1 = [t608 * MDP(2) + (qJDD(1) * t532 - 0.2e1 * t535 * t609) * MDP(7) + t690 * MDP(5) + (t468 * t408 - t434 * t463 + t462 * t435 - t469 * t545) * MDP(15) + (-t462 * t434 - t468 * t545) * MDP(14) + (-g(1) * t656 - g(2) * t657 - t398 * t626 + t480 * t434 - t437 * t625 - t450 * t468 - t494 * t462 + t506 * t545) * MDP(20) + (t353 * t390 + t377 * t358 + t356 * t400 + t381 * t360 + t355 * t391 + t376 * t359 - g(1) * (pkin(5) * t458 + qJ(6) * t457 + t524) - g(2) * (pkin(5) * t456 + qJ(6) * t455 + t643) + (-g(1) * t570 + g(2) * t540) * t539 + (-g(1) * (-pkin(1) + t540) - g(2) * t570) * t536) * MDP(31) + (-t604 * pkin(1) - g(1) * (-pkin(1) * t536 + t524) - g(2) * t643 + (t622 + t528) * qJ(2)) * MDP(6) + t647 * MDP(16) + 0.2e1 * (-t535 * t630 + t632 * t642) * MDP(8) + (qJDD(2) - t608 - 0.2e1 * t676) * MDP(4) + ((t664 - t666) * t434 + (-t672 + t671 + (-t663 - t667) * qJD(5)) * t468) * MDP(22) + (-t399 * t626 + t506 * t408 + t480 * t435 - t436 * t625 + t450 * t469 + t494 * t463 - t521 * t594) * MDP(19) + (-t538 * t542 - t629) * MDP(10) + t594 * MDP(3) + (-t355 * t469 - t356 * t660 - t359 * t698 + t360 * t438 - t376 * t435 - t381 * t575 + t388 * t400 - t391 * t407 + t595) * MDP(28) + (-t388 * t469 + t407 * t660 - t435 * t438 + t575 * t698) * MDP(24) + (t601 * t469 + t385 * t435 + t399 * t438 + t436 * t388 + ((-qJD(5) * t437 + t396) * t698 + t431 * t407 - t409 * qJD(5) * t468) * t537 + ((-qJD(5) * t431 - t398) * t698 - t437 * t407 - t368 * t468 + t409 * t434) * t534 + t595) * MDP(26) + (t407 * t469 + t435 * t698) * MDP(25) + (-t387 * t469 - t407 * t659 - t435 * t569 + t574 * t698) * MDP(23) + (-t368 * t659 - t386 * t435 - t436 * t387 - t399 * t569 - t407 * t646 + t409 * t574 - t469 * t572 - t571 * t698 + t596) * MDP(27) + (t353 * t469 + t356 * t659 + t358 * t698 + t360 * t569 + t377 * t435 - t381 * t574 + t387 * t400 + t390 * t407 - t596) * MDP(30) + (t387 * t659 - t569 * t574) * MDP(21) + (-t358 * t438 - t359 * t569 - t387 * t391 - t388 * t390 + t594 * t522 + t588 * t434 + (qJD(5) * t587 + t353 * t534 - t355 * t537) * t468) * MDP(29) + (-t435 * t626 - t469 * t625) * MDP(17) + qJDD(1) * MDP(1) + (-t535 * t542 + t520) * MDP(9) + (t535 * t556 + t538 * t573) * MDP(12) + (-t535 * t573 + t538 * t556) * MDP(13); qJDD(1) * MDP(4) - t543 * MDP(5) + (t604 + t566) * MDP(6) + (t535 * t641 + t520) * MDP(12) + (t538 * t641 - t629) * MDP(13) + t647 * MDP(19) + (t356 * t468 - t381 * t434 - t608) * MDP(31) + (-t626 * MDP(20) + (-t664 - t666) * MDP(29) + t587 * MDP(31)) * t435 + (-t463 * MDP(19) + t462 * MDP(20) + t550) * qJD(1) + (t691 * (-qJD(1) * t537 - t435 * t534) + t627 * (qJD(1) * t534 - t435 * t537)) * t698 + (-t625 * MDP(20) + t586 * MDP(29) + t589 * MDP(31) + (-t534 * t691 - t537 * t627) * t407 + ((t534 * t627 - t537 * t691) * t698 + t550) * qJD(5)) * t469 + t691 * (t468 * t388 - t434 * t438) - t627 * (t387 * t468 - t434 * t569); (t465 * t388 + t598 * t537 + (-t670 + t674) * t534 - t650 * t438 + (-t513 * t635 - t534 * t603 + t380) * t698 + t578) * MDP(28) + (t421 * t626 + (-t463 * t639 - t612 * t626 + t625 * t686) * pkin(3) + t559) * MDP(19) + (t465 * t387 + (-qJD(5) * t381 + t670) * t537 - t650 * t569 + (-t379 + t561) * t698 + t563) * MDP(30) + (t514 * t388 + t597 * t537 + t584 * t534 + t580 * t438 + ((-qJD(5) * t513 - t425) * t537 + (-t603 + t422) * t534) * t698 + t577) * MDP(26) + (-t514 * t387 + t584 * t537 - t580 * t569 + (-t561 + t648) * t698 + t565) * MDP(27) + (t379 * t438 + t380 * t569 + (-t438 * t617 - t569 * t618) * qJD(4) * pkin(3) + t557 * t513 + t552) * MDP(29) + t549 - MDP(10) * t631 + MDP(9) * t630 + (t422 * t626 + (t462 * t639 - t610 * t626 - t625 * t685) * pkin(3) + t547) * MDP(20) + (t356 * t465 - t377 * t379 - t376 * t380 - g(3) * t694 - t650 * t381 + (-t538 * t681 + (t376 * t618 + t377 * t617) * qJD(4)) * pkin(3) + t558 * t513 - (-pkin(3) * t538 - pkin(9) * t521 - t567) * t527 + t553) * MDP(31) + (g(3) * t535 + t538 * t566 + t470) * MDP(12) + qJDD(3) * MDP(11) + (g(3) * t538 + (-t491 - t566) * t535) * MDP(13) + (t538 * t535 * MDP(7) - MDP(8) * t642) * t543; (-t381 * t635 - t382 * t698 - t387 * t579 + t569 * t645 + t563 - t697) * MDP(30) + (-pkin(4) * t388 - t417 * t438 + (t416 * t698 + t668 - t682) * t534 + ((-t430 - t677) * t698 + t597) * t537 + t577) * MDP(26) + (pkin(4) * t387 + t409 * t661 + t417 * t569 + t649 * t698 + t565 + t697) * MDP(27) + (t383 * t698 - t388 * t579 + (t674 - t682) * t534 + t645 * t438 + (-t677 * t698 + t598) * t537 + t578) * MDP(28) + (pkin(9) * t557 + t382 * t438 + t383 * t569 + t552) * MDP(29) + t549 + (-t356 * t579 - t377 * t382 - t376 * t383 - g(3) * t507 + t645 * t381 + t567 * t527 + (t558 + t488) * pkin(9) + t553) * MDP(31) + (t416 * t626 + t547) * MDP(20) + (t417 * t626 + t559) * MDP(19); -MDP(21) * t665 + (-t438 ^ 2 + t688) * MDP(22) + t370 * MDP(23) + (-t388 - t702) * MDP(24) + t407 * MDP(25) + (t409 * t569 + t560 + t673) * MDP(26) + (t385 * t698 + t409 * t438 - t551) * MDP(27) + (-t402 * t438 - t554 + t673 + 0.2e1 * t684) * MDP(28) + (pkin(5) * t387 - qJ(6) * t388 - (t377 - t386) * t569 + (t376 - t633) * t438) * MDP(29) + (0.2e1 * t675 - t381 * t438 - t402 * t569 + (0.2e1 * qJD(6) - t385) * t698 + t551) * MDP(30) + (t353 * qJ(6) - t355 * pkin(5) - t381 * t402 - t376 * t386 - g(1) * (-pkin(5) * t455 + qJ(6) * t456) - g(2) * (pkin(5) * t457 - qJ(6) * t458) + t591 * t509 + t633 * t377) * MDP(31); t370 * MDP(29) + (-t688 - t704) * MDP(30) + (-t377 * t698 + t554 - t684) * MDP(31) + (-t665 - t407) * MDP(28);];
tau  = t1;
