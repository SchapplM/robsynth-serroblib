% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:09:25
% EndTime: 2019-03-09 21:09:38
% DurationCPUTime: 8.17s
% Computational Cost: add. (7742->593), mult. (18551->705), div. (0->0), fcn. (12498->6), ass. (0->240)
t580 = cos(qJ(3));
t579 = sin(qJ(2));
t581 = cos(qJ(2));
t686 = t580 * t581;
t607 = pkin(3) * t579 - pkin(9) * t686;
t706 = -pkin(9) - pkin(8);
t640 = qJD(3) * t706;
t613 = pkin(2) * t579 - pkin(8) * t581;
t520 = t613 * qJD(1);
t578 = sin(qJ(3));
t663 = qJD(1) * t579;
t633 = t578 * t663;
t666 = pkin(7) * t633 + t580 * t520;
t737 = qJD(1) * t607 - t580 * t640 + t666;
t500 = t578 * t520;
t688 = t579 * t580;
t736 = -t500 - (-pkin(9) * t578 * t581 - pkin(7) * t688) * qJD(1) + t578 * t640;
t577 = sin(qJ(4));
t654 = qJD(3) * t579;
t636 = t578 * t654;
t656 = qJD(2) * t581;
t594 = t580 * t656 - t636;
t646 = qJD(2) * qJD(3);
t587 = qJD(1) * t594 + t580 * t646;
t657 = qJD(2) * t580;
t513 = -t633 + t657;
t659 = qJD(2) * t578;
t514 = t580 * t663 + t659;
t705 = cos(qJ(4));
t604 = t577 * t513 + t514 * t705;
t628 = qJD(1) * t654;
t647 = qJD(1) * qJD(2);
t629 = t581 * t647;
t641 = t580 * t628 + (t629 + t646) * t578;
t412 = qJD(4) * t604 + t577 * t587 + t705 * t641;
t462 = t604 ^ 2;
t632 = qJD(4) * t705;
t652 = qJD(4) * t577;
t411 = -t513 * t632 + t514 * t652 + t577 * t641 - t705 * t587;
t464 = -t705 * t513 + t514 * t577;
t662 = qJD(1) * t581;
t559 = -qJD(3) + t662;
t542 = -qJD(4) + t559;
t696 = t464 * t542;
t590 = t411 + t696;
t627 = MDP(22) * t663;
t695 = t604 * t542;
t708 = t464 ^ 2;
t734 = t464 * t604;
t735 = -t590 * MDP(20) + qJD(2) * t627 + MDP(18) * t734 + (-t412 - t695) * MDP(21) + (t462 - t708) * MDP(19);
t638 = t578 * t662;
t639 = t705 * t580;
t691 = t577 * t578;
t709 = qJD(3) + qJD(4);
t710 = t705 * qJD(3) + t632;
t672 = -t577 * t638 - t580 * t710 + t639 * t662 + t691 * t709;
t699 = qJ(5) * t464;
t732 = qJ(6) * t464;
t701 = qJD(2) * pkin(2);
t534 = pkin(7) * t663 - t701;
t605 = pkin(3) * t513 - t534;
t593 = qJ(5) * t604 + t605;
t707 = pkin(4) + pkin(5);
t385 = -t464 * t707 + qJD(6) + t593;
t731 = t385 * t464;
t403 = pkin(4) * t464 - t593;
t730 = t403 * t464;
t729 = t464 * t605;
t536 = t706 * t578;
t537 = t706 * t580;
t728 = t536 * t632 + t537 * t652 - t577 * t737 + t736 * t705;
t516 = t577 * t580 + t578 * t705;
t477 = t709 * t516;
t671 = -t516 * t662 + t477;
t571 = pkin(7) * t662;
t655 = qJD(3) * t578;
t616 = -t571 + (-t638 + t655) * pkin(3);
t630 = t579 * t647;
t712 = qJ(5) * t630 - t542 * qJD(5);
t653 = qJD(3) * t580;
t635 = t579 * t653;
t637 = t578 * t656;
t727 = t635 + t637;
t725 = -0.2e1 * t647;
t724 = pkin(4) * t604;
t723 = MDP(4) * t579;
t575 = t579 ^ 2;
t722 = MDP(5) * (-t581 ^ 2 + t575);
t721 = qJ(6) * t604;
t720 = t403 * t604;
t719 = t605 * t604;
t718 = t604 * t707;
t717 = t412 * qJ(6) + t464 * qJD(6);
t676 = qJ(5) * t663 - t728;
t484 = t577 * t536 - t705 * t537;
t716 = qJD(4) * t484 + t736 * t577 + t705 * t737;
t535 = qJD(2) * pkin(8) + t571;
t528 = -pkin(2) * t581 - pkin(8) * t579 - pkin(1);
t505 = t528 * qJD(1);
t690 = t578 * t505;
t473 = t580 * t535 + t690;
t447 = pkin(9) * t513 + t473;
t550 = pkin(3) * t632 + qJD(5);
t563 = pkin(3) * t577 + qJ(5);
t715 = -t550 * t542 + t563 * t630;
t714 = 0.2e1 * t712;
t711 = t672 * qJ(5) - qJD(5) * t516 + t616;
t650 = t575 * qJD(1);
t472 = t580 * t505 - t535 * t578;
t446 = -pkin(9) * t514 + t472;
t437 = -pkin(3) * t559 + t446;
t692 = t577 * t447;
t391 = t705 * t437 - t692;
t648 = qJD(5) - t391;
t557 = pkin(4) * t630;
t595 = t607 * qJD(2);
t523 = t613 * qJD(2);
t506 = qJD(1) * t523;
t620 = pkin(7) * t630;
t669 = t580 * t506 + t578 * t620;
t402 = qJD(1) * t595 - qJD(3) * t447 + t669;
t670 = t505 * t653 + t578 * t506;
t588 = -t535 * t655 - t580 * t620 + t670;
t413 = -pkin(9) * t641 + t588;
t618 = -t705 * t402 + t577 * t413 + t437 * t652 + t447 * t632;
t370 = -t557 + t618;
t367 = -pkin(5) * t630 + t411 * qJ(6) - qJD(6) * t604 + t370;
t589 = -t385 * t604 + t367;
t704 = pkin(7) * t578;
t703 = pkin(8) * t559;
t700 = qJ(5) * t412;
t444 = t705 * t447;
t392 = t577 * t437 + t444;
t698 = t392 * t542;
t697 = t412 * t563;
t694 = t514 * t559;
t693 = t534 * t578;
t689 = t578 * t579;
t583 = qJD(2) ^ 2;
t687 = t579 * t583;
t685 = t581 * t559;
t684 = t581 * t583;
t584 = qJD(1) ^ 2;
t683 = t581 * t584;
t398 = t705 * t446 - t692;
t383 = t398 + t721;
t682 = -t383 + t550;
t681 = -t671 * t707 - t711;
t643 = t579 * t707;
t680 = -qJ(6) * t672 - qJD(1) * t643 + t516 * qJD(6) - t716;
t515 = -t639 + t691;
t679 = -qJ(6) * t671 - qJD(6) * t515 + t676;
t678 = t398 - t550;
t677 = pkin(4) * t671 + t711;
t675 = pkin(4) * t663 + t716;
t512 = t580 * t528;
t471 = -pkin(9) * t688 + t512 + (-pkin(3) - t704) * t581;
t561 = pkin(7) * t686;
t665 = t578 * t528 + t561;
t478 = -pkin(9) * t689 + t665;
t673 = t577 * t471 + t705 * t478;
t668 = t578 * t523 + t528 * t653;
t658 = qJD(2) * t579;
t667 = t580 * t523 + t658 * t704;
t461 = t641 * pkin(3) + pkin(7) * t629;
t524 = pkin(3) * t689 + t579 * pkin(7);
t483 = -t536 * t705 - t577 * t537;
t661 = qJD(2) * t483;
t660 = qJD(2) * t484;
t380 = t391 + t721;
t649 = qJD(5) - t380;
t644 = pkin(3) * t652;
t642 = t577 * t689;
t486 = pkin(3) * t727 + pkin(7) * t656;
t569 = -pkin(3) * t580 - pkin(2);
t634 = t542 * t652;
t626 = MDP(15) * t658;
t625 = pkin(1) * t725;
t624 = t559 + t662;
t623 = -t513 + t657;
t622 = -t514 + t659;
t621 = qJD(3) + t662;
t619 = -t577 * t402 - t705 * t413 - t437 * t632 + t447 * t652;
t568 = -pkin(3) * t705 - pkin(4);
t617 = t705 * t656;
t397 = t577 * t446 + t444;
t382 = t397 + t732;
t615 = -t382 + t644;
t614 = -t397 + t644;
t420 = -qJ(5) * t581 + t673;
t494 = t579 * t639 - t642;
t611 = qJ(5) * t494 - t524;
t381 = t392 + t732;
t609 = t471 * t705 - t577 * t478;
t608 = -pkin(3) * t514 - t699;
t373 = t412 * pkin(4) + t411 * qJ(5) - qJD(5) * t604 + t461;
t606 = qJ(5) * t516 - t569;
t421 = t581 * pkin(4) - t609;
t369 = -t619 + t712;
t602 = t621 * t659;
t601 = -t391 * t542 + t619;
t600 = -t397 * t542 - t618;
t599 = -t398 * t542 + t619;
t598 = t618 + t720;
t424 = t595 + (-t561 + (pkin(9) * t579 - t528) * t578) * qJD(3) + t667;
t427 = -t727 * pkin(9) + (-t579 * t657 - t581 * t655) * pkin(7) + t668;
t597 = -t424 * t705 + t577 * t427 + t471 * t652 + t478 * t632;
t596 = t577 * t424 + t705 * t427 + t471 * t632 - t478 * t652;
t368 = -pkin(5) * t412 - t373;
t366 = t369 + t717;
t428 = t477 * t579 + t577 * t637 - t580 * t617;
t592 = -qJ(5) * t428 + qJD(5) * t494 - t486;
t375 = qJ(5) * t658 - qJD(5) * t581 + t596;
t562 = -pkin(5) + t568;
t527 = t542 * qJ(5);
t509 = pkin(3) * t634;
t493 = t516 * t579;
t460 = pkin(4) * t515 - t606;
t450 = qJ(6) * t515 + t484;
t449 = -t516 * qJ(6) + t483;
t445 = -t515 * t707 + t606;
t438 = pkin(4) * t493 - t611;
t429 = t578 * t617 - t577 * t636 - qJD(4) * t642 + (t577 * t656 + t579 * t710) * t580;
t423 = t699 + t724;
t422 = -t493 * t707 + t611;
t415 = -t608 + t724;
t408 = qJ(6) * t493 + t420;
t401 = t581 * pkin(5) - t494 * qJ(6) + t421;
t396 = -t699 - t718;
t388 = -t527 + t392;
t387 = pkin(4) * t542 + t648;
t386 = t608 - t718;
t379 = t381 - t527;
t378 = t542 * t707 + t649;
t377 = pkin(4) * t429 - t592;
t376 = -pkin(4) * t658 + t597;
t374 = -t429 * t707 + t592;
t372 = qJ(6) * t429 + qJD(6) * t493 + t375;
t371 = t428 * qJ(6) - qJD(2) * t643 - t494 * qJD(6) + t597;
t1 = [(t366 * t408 + t367 * t401 + t368 * t422 + t371 * t378 + t372 * t379 + t374 * t385) * MDP(32) - MDP(7) * t687 + (pkin(7) * t687 + t581 * t625) * MDP(10) + (t514 * t594 + t587 * t688) * MDP(11) + 0.2e1 * t629 * t723 + (t411 * t581 + t428 * t542 + (qJD(1) * t494 + t604) * t658) * MDP(20) + (-t369 * t581 - t373 * t494 - t375 * t542 - t377 * t604 + t403 * t428 + t411 * t438 + (qJD(1) * t420 + t388) * t658) * MDP(27) + (-t366 * t581 + t368 * t494 - t372 * t542 + t374 * t604 - t385 * t428 - t411 * t422 + (qJD(1) * t408 + t379) * t658) * MDP(30) + (-t411 * t494 - t428 * t604) * MDP(18) + (-t369 * t493 + t370 * t494 - t375 * t464 + t376 * t604 - t387 * t428 - t388 * t429 - t411 * t421 - t412 * t420) * MDP(26) + (t411 * t493 - t412 * t494 + t428 * t464 - t429 * t604) * MDP(19) + (t366 * t493 - t367 * t494 - t371 * t604 + t372 * t464 + t378 * t428 + t379 * t429 + t401 * t411 + t408 * t412) * MDP(31) + (-pkin(7) * t684 + t579 * t625) * MDP(9) + (t559 * t635 + t641 * t581 + (t513 * t579 + (-t650 + t685) * t578) * qJD(2)) * MDP(14) - t624 * t626 + ((t580 * t513 - t514 * t578) * t656 + ((-t513 + t633) * t655 + (-t514 * qJD(3) - t602 - t641) * t580) * t579) * MDP(12) + (t412 * t581 + t429 * t542 + (-qJD(1) * t493 - t464) * t658) * MDP(21) + (t370 * t581 + t373 * t493 + t376 * t542 + t377 * t464 + t403 * t429 + t412 * t438 + (-qJD(1) * t421 - t387) * t658) * MDP(25) + (t367 * t581 - t368 * t493 + t371 * t542 - t374 * t464 - t385 * t429 - t412 * t422 + (-qJD(1) * t401 - t378) * t658) * MDP(29) + (t596 * t542 - t619 * t581 + t486 * t604 - t524 * t411 + t461 * t494 + t605 * t428 + (-qJD(1) * t673 - t392) * t658) * MDP(24) + (t597 * t542 + t618 * t581 + t486 * t464 + t524 * t412 + t461 * t493 - t605 * t429 + (qJD(1) * t609 + t391) * t658) * MDP(23) + (-t542 - t662) * MDP(22) * t658 + (t624 * t636 + (t514 * t579 + (t650 + (-t559 - t621) * t581) * t580) * qJD(2)) * MDP(13) + (t369 * t420 + t370 * t421 + t373 * t438 + t375 * t388 + t376 * t387 + t377 * t403) * MDP(28) + t722 * t725 + (-(-t528 * t655 + t667) * t559 + (pkin(7) * t641 + t534 * t653 + (t512 * qJD(1) + t472) * qJD(2)) * t579 + ((-pkin(7) * t513 + t693) * qJD(2) + (t690 + (pkin(7) * t559 + t535) * t580) * qJD(3) - t669) * t581) * MDP(16) + MDP(6) * t684 + (t668 * t559 + t670 * t581 + (-t534 * t579 - t535 * t581 + (-t685 - t650) * pkin(7)) * t655 + ((pkin(7) * t514 + t534 * t580) * t581 + (-t665 * qJD(1) - t473 + (-t559 + t621) * pkin(7) * t580) * t579) * qJD(2)) * MDP(17); t559 * MDP(15) * t663 + t542 * t627 - t683 * t723 + (-t369 * t515 + t370 * t516 - t387 * t672 - t388 * t671 - t411 * t483 - t412 * t484 + t464 * t676 + t604 * t675) * MDP(26) + (-t373 * t516 + t411 * t460 + t676 * t542 - t677 * t604 + t672 * t403 + (-t388 + t660) * t663) * MDP(27) + (t366 * t515 - t367 * t516 + t378 * t672 + t379 * t671 + t411 * t449 + t412 * t450 - t464 * t679 + t604 * t680) * MDP(31) + (t368 * t516 - t411 * t445 + t679 * t542 + t681 * t604 - t672 * t385 + (qJD(2) * t450 - t379) * t663) * MDP(30) + (t411 * t515 - t412 * t516 + t464 * t672 - t604 * t671) * MDP(19) + (t672 * t542 + (qJD(2) * t516 - t604) * t663) * MDP(20) + (-t411 * t516 - t604 * t672) * MDP(18) + (t559 * t655 + (-t685 * t578 + t579 * t623) * qJD(1)) * MDP(14) + (-t559 * t653 + (t579 * t622 + t580 * t685) * qJD(1)) * MDP(13) + (t373 * t515 + t412 * t460 + t675 * t542 + t677 * t464 + t671 * t403 + (t387 - t661) * t663) * MDP(25) + (t369 * t484 + t370 * t483 + t373 * t460 + t387 * t675 - t388 * t676 + t403 * t677) * MDP(28) + (t366 * t450 + t367 * t449 + t368 * t445 - t378 * t680 - t379 * t679 + t385 * t681) * MDP(32) + (-t368 * t515 - t412 * t445 - t680 * t542 - t681 * t464 - t671 * t385 + (-qJD(2) * t449 + t378) * t663) * MDP(29) + (t671 * t542 + (-qJD(2) * t515 + t464) * t663) * MDP(21) + (-t569 * t411 + t461 * t516 + t728 * t542 + t672 * t605 + t616 * t604 + (t392 - t660) * t663) * MDP(24) + (t569 * t412 + t461 * t515 + t716 * t542 - t671 * t605 + t616 * t464 + (-t391 - t661) * t663) * MDP(23) + (MDP(9) * t579 * t584 + MDP(10) * t683) * pkin(1) + (-t578 ^ 2 * t628 + (t602 - t694) * t580) * MDP(11) + ((-t641 + t694) * t578 + ((t513 + t657) * qJD(3) + (t581 * t623 - t636) * qJD(1)) * t580) * MDP(12) + t584 * t722 + (-t500 * t559 + (-t578 * t703 + (t534 - t701) * t580) * qJD(3) + ((-t534 - t701) * t686 + (pkin(2) * t655 - pkin(8) * t657 + t473) * t579 + (t559 * t688 + t581 * t622) * pkin(7)) * qJD(1)) * MDP(17) + (-pkin(2) * t641 + t666 * t559 + (t580 * t703 + t693) * qJD(3) + ((-pkin(8) * t659 - t472) * t579 + (-pkin(7) * t623 - t693) * t581) * qJD(1)) * MDP(16); qJD(1) * t626 + (t369 * t563 + t370 * t568 + t387 * t614 - t388 * t678 - t403 * t415) * MDP(28) + (t366 * t563 + t367 * t562 + t378 * t615 + t379 * t682 - t385 * t386) * MDP(32) + (-t415 * t464 - t568 * t630 + t509 + t557 + t600 - t720) * MDP(25) + (t719 + (-t464 * t514 + t630 * t705 + t634) * pkin(3) + t600) * MDP(23) + (-t382 * t542 + t386 * t464 - t562 * t630 + t509 - t589) * MDP(29) + (-t411 * t568 - t697 + (t388 + t614) * t604 + (t387 + t678) * t464) * MDP(26) + (t415 * t604 - t599 + t712 + t715 - t730) * MDP(27) + (t383 * t542 - t386 * t604 + t366 + t715 + t731) * MDP(30) + (t411 * t562 + t697 + (-t379 - t615) * t604 + (-t378 + t682) * t464) * MDP(31) + (-t729 + (-t514 * t604 + t542 * t632 - t577 * t630) * pkin(3) + t599) * MDP(24) + (-t514 * t534 + t669 + (-qJD(3) - t559) * t473) * MDP(16) + (-t513 ^ 2 + t514 ^ 2) * MDP(12) + (t513 * t559 + t587) * MDP(13) + (-t472 * t559 - t513 * t534 - t588) * MDP(17) + (-t641 - t694) * MDP(14) - t514 * t513 * MDP(11) + t735; (-t618 - t698 + t719) * MDP(23) + (t601 - t729) * MDP(24) + (-t423 * t464 + 0.2e1 * t557 - t598 - t698) * MDP(25) + (pkin(4) * t411 - t700 + (t388 - t392) * t604 + (t387 - t648) * t464) * MDP(26) + (t423 * t604 - t601 + t714 - t730) * MDP(27) + (-pkin(4) * t370 + qJ(5) * t369 - t387 * t392 + t388 * t648 - t403 * t423) * MDP(28) + (-t381 * t542 + t396 * t464 + t630 * t707 - t589) * MDP(29) + (t380 * t542 - t396 * t604 - t619 + t714 + t717 + t731) * MDP(30) + (t700 - t411 * t707 + (-t379 + t381) * t604 + (-t378 + t649) * t464) * MDP(31) + (qJ(5) * t366 - t367 * t707 - t378 * t381 + t379 * t649 - t385 * t396) * MDP(32) + t735; (t388 * t542 - t557 + t598) * MDP(28) + (t379 * t542 + t589) * MDP(32) + (MDP(25) + MDP(29)) * (-t630 + t734) + (MDP(27) + MDP(30)) * (-t542 ^ 2 - t462) + (-MDP(26) + MDP(31)) * t590; (-t412 + t695) * MDP(29) + (-t411 + t696) * MDP(30) + (-t462 - t708) * MDP(31) + (t378 * t604 - t379 * t464 + t368) * MDP(32);];
tauc  = t1;
