% Calculate vector of inverse dynamics joint torques for
% S6RRPPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:44:13
% EndTime: 2019-03-09 08:44:24
% DurationCPUTime: 9.97s
% Computational Cost: add. (5587->601), mult. (11484->712), div. (0->0), fcn. (7364->10), ass. (0->257)
t740 = MDP(24) + MDP(26);
t715 = cos(qJ(5));
t643 = qJD(5) * t715;
t568 = sin(qJ(2));
t672 = qJD(1) * t568;
t738 = t715 * t672 + t643;
t567 = sin(qJ(5));
t664 = qJD(5) * t567;
t737 = -t567 * t672 - t664;
t563 = sin(pkin(9));
t707 = pkin(2) + qJ(4);
t705 = -pkin(8) - t707;
t501 = t705 * t563;
t564 = cos(pkin(9));
t502 = t705 * t564;
t447 = t501 * t715 + t567 * t502;
t570 = cos(qJ(2));
t658 = qJD(1) * qJD(2);
t640 = t570 * t658;
t657 = qJDD(1) * t568;
t597 = t640 + t657;
t491 = qJDD(5) + t597;
t557 = pkin(9) + qJ(5);
t547 = cos(t557);
t569 = sin(qJ(1));
t571 = cos(qJ(1));
t625 = g(1) * t571 + g(2) * t569;
t611 = t625 * t570;
t709 = g(3) * t568;
t731 = t611 + t709;
t736 = t447 * t491 + t547 * t731;
t528 = qJD(5) + t672;
t602 = -t563 * t567 + t715 * t564;
t677 = -t563 * t738 + t564 * t737;
t622 = t602 * t491 + t528 * t677;
t671 = qJD(1) * t570;
t645 = t564 * t671;
t662 = t563 * qJD(2);
t488 = -t645 - t662;
t646 = t563 * t671;
t668 = qJD(2) * t564;
t489 = -t646 + t668;
t437 = -t715 * t488 + t489 * t567;
t735 = t437 ^ 2;
t605 = -t567 * t488 - t489 * t715;
t717 = t605 ^ 2;
t716 = pkin(3) + pkin(7);
t734 = t437 * t528;
t733 = t528 * t605;
t470 = t602 * t570;
t676 = t563 * t737 + t564 * t738;
t539 = pkin(7) * t672;
t730 = qJD(3) + t539;
t729 = MDP(25) - MDP(28);
t543 = pkin(2) * t672;
t703 = qJ(3) * t570;
t620 = qJ(4) * t568 - t703;
t474 = qJD(1) * t620 + t543;
t540 = pkin(7) * t671;
t499 = pkin(3) * t671 + t540;
t430 = -t474 * t563 + t564 * t499;
t689 = t563 * t568;
t614 = pkin(4) * t570 - pkin(8) * t689;
t412 = qJD(1) * t614 + t430;
t431 = t564 * t474 + t563 * t499;
t421 = pkin(8) * t564 * t672 + t431;
t603 = -t563 * t715 - t567 * t564;
t604 = -t567 * t501 + t502 * t715;
t728 = -qJD(4) * t603 - qJD(5) * t604 + t567 * t412 + t715 * t421;
t727 = -qJD(4) * t602 - qJD(5) * t447 - t412 * t715 + t567 * t421;
t549 = t568 * qJ(3);
t553 = t570 * pkin(2);
t674 = t553 + t549;
t702 = qJ(4) * t570;
t623 = t674 + t702;
t487 = -pkin(1) - t623;
t511 = t716 * t568;
t495 = t564 * t511;
t423 = pkin(4) * t568 + t495 + (pkin(8) * t570 - t487) * t563;
t443 = t564 * t487 + t563 * t511;
t688 = t564 * t570;
t429 = -pkin(8) * t688 + t443;
t726 = t567 * t423 + t715 * t429;
t535 = pkin(4) * t564 + pkin(3);
t661 = t535 * t672 + t730;
t560 = qJD(2) * qJ(3);
t724 = qJD(4) + t560;
t656 = qJDD(1) * t570;
t515 = t563 * t656;
t635 = qJDD(2) * t564 - t515;
t641 = t568 * t658;
t587 = t563 * t641 + t635;
t675 = t563 * qJDD(2) + t564 * t656;
t594 = t564 * t641 - t675;
t393 = -t488 * t643 + t489 * t664 - t567 * t594 - t715 * t587;
t394 = -qJD(5) * t605 + t567 * t587 - t715 * t594;
t723 = t394 * pkin(5) + t393 * qJ(6) + qJD(6) * t605;
t722 = -t393 * t602 - t605 * t677;
t598 = t491 * t603 - t528 * t676;
t478 = t499 + t724;
t666 = qJD(2) * t570;
t721 = t568 * (-t478 + t724) + t707 * t666;
t637 = -pkin(1) - t549;
t593 = -t570 * t707 + t637;
t460 = t593 * qJD(1);
t660 = pkin(3) * t672 + t730;
t469 = -qJD(2) * t707 + t660;
t416 = -t460 * t563 + t564 * t469;
t417 = t564 * t460 + t563 * t469;
t619 = t416 * t563 - t417 * t564;
t556 = g(3) * t570;
t686 = t568 * t571;
t687 = t568 * t569;
t653 = g(1) * t686 + g(2) * t687 - t556;
t527 = pkin(2) * t641;
t665 = qJD(3) * t568;
t578 = qJD(2) * t620 - qJD(4) * t570 - t665;
t409 = qJD(1) * t578 + qJDD(1) * t593 + t527;
t526 = pkin(7) * t640;
t536 = pkin(7) * t657;
t638 = qJDD(3) + t526 + t536;
t432 = pkin(3) * t597 - qJD(2) * qJD(4) - qJDD(2) * t707 + t638;
t389 = -t563 * t409 + t564 * t432;
t698 = t389 * t564;
t720 = t619 * t672 + t653 - t698;
t667 = qJD(2) * t568;
t542 = pkin(2) * t667;
t448 = t542 + t578;
t500 = t716 * t666;
t419 = -t448 * t563 + t564 * t500;
t403 = qJD(2) * t614 + t419;
t420 = t564 * t448 + t563 * t500;
t648 = t564 * t667;
t408 = pkin(8) * t648 + t420;
t719 = -qJD(5) * t726 + t403 * t715 - t567 * t408;
t380 = pkin(4) * t597 - pkin(8) * t587 + t389;
t390 = t564 * t409 + t563 * t432;
t385 = pkin(8) * t594 + t390;
t398 = pkin(4) * t672 - pkin(8) * t489 + t416;
t402 = pkin(8) * t488 + t417;
t596 = t567 * t380 + t715 * t385 + t398 * t643 - t402 * t664;
t701 = qJ(6) * t491;
t369 = qJD(6) * t528 + t596 + t701;
t634 = -t715 * t380 + t567 * t385 + t398 * t664 + t402 * t643;
t714 = pkin(5) * t491;
t370 = qJDD(6) + t634 - t714;
t378 = t398 * t715 - t567 * t402;
t659 = qJD(6) - t378;
t374 = -t528 * pkin(5) + t659;
t379 = t567 * t398 + t402 * t715;
t375 = t528 * qJ(6) + t379;
t718 = -t369 * t603 - t370 * t602 - t374 * t677 + t375 * t676 - t653;
t713 = g(1) * t569;
t710 = g(2) * t571;
t548 = t563 * pkin(4);
t704 = pkin(7) * qJDD(2);
t700 = qJDD(2) * pkin(2);
t699 = t379 * t528;
t697 = t390 * t563;
t694 = t437 * t605;
t561 = t568 ^ 2;
t573 = qJD(1) ^ 2;
t691 = t561 * t573;
t685 = t568 * t573;
t684 = t569 * t570;
t683 = t570 * t571;
t533 = qJ(3) + t548;
t682 = pkin(5) * t676 - qJ(6) * t677 - qJD(6) * t602 + t661;
t680 = -qJ(6) * t671 - t728;
t679 = pkin(5) * t671 - t727;
t512 = t716 * t570;
t562 = t570 ^ 2;
t673 = t561 - t562;
t558 = qJDD(2) * qJ(3);
t655 = pkin(4) * t689;
t654 = t570 * t685;
t537 = pkin(7) * t656;
t559 = qJD(2) * qJD(3);
t652 = t537 + t558 + t559;
t479 = pkin(4) * t688 + t512;
t651 = t716 * qJD(2);
t644 = g(3) * t674;
t642 = t675 * pkin(4);
t639 = t707 * t657;
t636 = -qJD(2) * pkin(2) + qJD(3);
t633 = t571 * pkin(1) + pkin(2) * t683 + t569 * pkin(7) + qJ(3) * t686;
t632 = -t536 + t653;
t630 = qJD(1) * t651;
t572 = qJD(2) ^ 2;
t629 = pkin(7) * t572 + t710;
t546 = sin(t557);
t465 = t546 * t569 - t547 * t686;
t467 = t546 * t571 + t547 * t687;
t628 = g(1) * t467 + g(2) * t465;
t466 = t546 * t686 + t547 * t569;
t468 = -t546 * t687 + t547 * t571;
t627 = -g(1) * t468 - g(2) * t466;
t522 = qJ(3) * t684;
t525 = qJ(3) * t683;
t626 = -g(1) * t525 - g(2) * t522;
t624 = -t710 + t713;
t504 = t539 + t636;
t510 = -t540 - t560;
t618 = t504 * t570 + t510 * t568;
t617 = pkin(3) * t656 + qJDD(4) + t652;
t616 = qJD(2) * (-pkin(7) - t535);
t615 = t564 * t635;
t613 = t637 - t553;
t610 = -0.2e1 * pkin(1) * t658 - t704;
t608 = t423 * t715 - t567 * t429;
t486 = t613 * qJD(1);
t601 = t486 * t672 + qJDD(3) - t632;
t463 = t568 * t616;
t600 = -qJ(3) * t666 - t665;
t444 = -pkin(4) * t488 + t478;
t599 = qJD(1) * t616;
t595 = t567 * t403 + t715 * t408 + t423 * t643 - t429 * t664;
t592 = 0.2e1 * qJDD(1) * pkin(1) - t629;
t591 = pkin(5) * t546 - qJ(6) * t547 + t548;
t435 = -t568 * t630 + t617;
t590 = -t435 * t570 + t625;
t471 = t603 * t570;
t505 = -pkin(1) - t674;
t588 = t704 + (-qJD(1) * t505 - t486) * qJD(2);
t584 = g(1) * t465 - g(2) * t467 + t547 * t556 - t634;
t583 = t435 - t731;
t434 = qJD(1) * t600 + qJDD(1) * t613 + t527;
t476 = t542 + t600;
t582 = qJD(1) * t476 + qJDD(1) * t505 + t434 + t629;
t580 = -t611 + t617;
t461 = pkin(7) * t641 - t652;
t473 = t638 - t700;
t579 = qJD(2) * t618 - t461 * t570 + t473 * t568;
t388 = pkin(5) * t437 + qJ(6) * t605 + t444;
t577 = -t388 * t605 + qJDD(6) - t584;
t576 = t491 * t604 - t546 * t731;
t575 = -g(1) * t466 + g(2) * t468 + t546 * t556 + t596;
t406 = t568 * t599 + t617 + t642;
t565 = -pkin(8) - qJ(4);
t554 = t571 * pkin(7);
t531 = g(1) * t684;
t498 = t568 * t651;
t496 = -qJ(3) * t671 + t543;
t442 = -t487 * t563 + t495;
t433 = -pkin(5) * t603 - qJ(6) * t602 + t533;
t426 = -t567 * t568 * t662 - qJD(5) * t471 + t648 * t715;
t425 = -qJD(5) * t470 - t603 * t667;
t407 = pkin(5) * t470 - qJ(6) * t471 + t479;
t397 = -pkin(5) * t605 + qJ(6) * t437;
t392 = -t568 * pkin(5) - t608;
t391 = qJ(6) * t568 + t726;
t384 = -pkin(5) * t426 - qJ(6) * t425 - qJD(6) * t471 + t463;
t381 = -t393 + t734;
t373 = -pkin(5) * t666 - t719;
t372 = qJ(6) * t666 + qJD(6) * t568 + t595;
t371 = t406 + t723;
t1 = [qJDD(1) * MDP(1) + (pkin(7) * t579 - g(1) * t554 - g(2) * t633 + t434 * t505 + t486 * t476 - t613 * t713) * MDP(14) + (t610 * t570 + (-t592 - t713) * t568) * MDP(10) + (t588 * t570 + (-t582 + t713) * t568) * MDP(13) + (t420 * t488 - t443 * t675 - t419 * t489 - t442 * t635 + t531 + ((-t442 * t563 + t443 * t564) * qJD(1) - t619) * t667 + (t389 * t563 - t390 * t564 - t710) * t570) * MDP(17) + (t369 * t391 + t375 * t372 + t371 * t407 + t388 * t384 + t370 * t392 + t374 * t373 - g(1) * (pkin(5) * t468 + qJ(6) * t467 + t535 * t571 + t554) - g(2) * (pkin(5) * t466 + qJ(6) * t465 - t565 * t683 + t571 * t655 + t633) + (-g(1) * (t565 * t570 + t613 - t655) - g(2) * t535) * t569) * MDP(29) + (qJDD(2) * t568 + t570 * t572) * MDP(6) + (qJDD(2) * t570 - t568 * t572) * MDP(7) + (-t393 * t568 + t425 * t528 + t471 * t491 - t605 * t666) * MDP(21) + (t369 * t568 - t371 * t471 + t372 * t528 + t375 * t666 + t384 * t605 - t388 * t425 + t391 * t491 + t393 * t407 - t628) * MDP(28) + (-t393 * t471 - t425 * t605) * MDP(19) + (-t379 * t666 - t479 * t393 + t406 * t471 + t444 * t425 - t463 * t605 - t491 * t726 - t528 * t595 - t568 * t596 + t628) * MDP(25) + (t378 * t666 + t479 * t394 + t406 * t470 - t444 * t426 + t463 * t437 + t608 * t491 + t528 * t719 - t634 * t568 + t627) * MDP(24) + (-g(2) * t683 - t369 * t470 + t370 * t471 - t372 * t437 - t373 * t605 + t374 * t425 + t375 * t426 - t391 * t394 - t392 * t393 + t531) * MDP(27) + (t393 * t470 - t394 * t471 - t425 * t437 - t426 * t605) * MDP(20) + (-t394 * t568 + t426 * t528 - t437 * t666 - t470 * t491) * MDP(22) + (-t370 * t568 + t371 * t470 - t373 * t528 - t374 * t666 + t384 * t437 - t388 * t426 - t392 * t491 + t394 * t407 + t627) * MDP(26) + (t390 * t443 + t417 * t420 + t389 * t442 + t416 * t419 + t435 * t512 - t478 * t498 - g(1) * (pkin(3) * t571 + t554) - g(2) * (qJ(4) * t683 + t633) + (-g(1) * (t613 - t702) - g(2) * pkin(3)) * t569) * MDP(18) + (t498 * t488 + t512 * t675 + (qJD(1) * t442 + t416) * t666 - t590 * t564 + (-t478 * t668 + t442 * qJDD(1) + t389 + t624 * t563 + (-t512 * t668 + t419) * qJD(1)) * t568) * MDP(15) + (-t498 * t489 + t512 * t635 + (-qJD(1) * t443 - t417) * t666 + t590 * t563 + (t478 * t662 - t443 * qJDD(1) - t390 + t624 * t564 + (t512 * t662 - t420) * qJD(1)) * t568) * MDP(16) + 0.2e1 * (t568 * t656 - t658 * t673) * MDP(5) + (t491 * t568 + t528 * t666) * MDP(23) + (qJDD(1) * t561 + 0.2e1 * t568 * t640) * MDP(4) + t624 * MDP(2) + ((t561 + t562) * qJDD(1) * pkin(7) + t579 - t625) * MDP(11) + t625 * MDP(3) + (t568 * t588 + t570 * t582 - t531) * MDP(12) + (t568 * t610 + t570 * t592 + t531) * MDP(9); (t379 * t671 - t533 * t393 + t406 * t602 + t677 * t444 + t528 * t728 - t605 * t661 - t736) * MDP(25) + (-t371 * t602 - t375 * t671 - t677 * t388 + t393 * t433 + t680 * t528 + t682 * t605 + t736) * MDP(28) + (t563 * t639 - qJ(3) * t515 + t660 * t489 + (t583 + t558) * t564 + (t417 * t570 + t431 * t568 + t563 * t721) * qJD(1)) * MDP(16) + (-t564 * t639 + qJ(3) * t675 - t660 * t488 + t583 * t563 + (-t416 * t570 - t430 * t568 - t564 * t721) * qJD(1)) * MDP(15) + t722 * MDP(19) + qJDD(2) * MDP(8) + (t709 - t537 + (pkin(1) * t573 + t625) * t570) * MDP(10) + ((-pkin(2) * t568 + t703) * qJDD(1) + ((-t510 - t560) * t568 + (-t504 + t636) * t570) * qJD(1)) * MDP(11) + (-t496 * t671 + t601 - 0.2e1 * t700) * MDP(12) + (-t461 * qJ(3) - t510 * qJD(3) - t473 * pkin(2) - t486 * t496 - g(1) * (-pkin(2) * t686 + t525) - g(2) * (-pkin(2) * t687 + t522) - t644 - t618 * qJD(1) * pkin(7)) * MDP(14) + (pkin(1) * t685 + t632) * MDP(9) + MDP(7) * t656 + MDP(6) * t657 + (t605 * t671 + t622) * MDP(21) + (t393 * t604 - t394 * t447 - t437 * t680 - t605 * t679 - t718) * MDP(27) + (t369 * t447 + t371 * t433 - t370 * t604 - t644 + t682 * t388 + t680 * t375 + t679 * t374 + (g(3) * t565 - t591 * t625) * t570 + (-g(3) * t591 + t625 * (pkin(2) - t565)) * t568 + t626) * MDP(29) + (-t393 * t603 - t394 * t602 - t437 * t677 + t605 * t676) * MDP(20) + (-t371 * t603 + t374 * t671 + t388 * t676 + t394 * t433 + t437 * t682 - t528 * t679 + t576) * MDP(26) + (-t378 * t671 + t533 * t394 - t406 * t603 + t661 * t437 + t676 * t444 + t528 * t727 + t576) * MDP(24) + (t707 * t615 - t431 * t488 + (qJD(4) * t564 + t430) * t489 + (-qJD(4) * t488 + t675 * t707 - t390) * t563 + t720) * MDP(17) + (t435 * qJ(3) - t417 * t431 - t416 * t430 - g(3) * t623 + t660 * t478 + (-t416 * t564 - t417 * t563) * qJD(4) + t626 + (t568 * t625 - t697 - t698) * t707) * MDP(18) + (t437 * t671 + t598) * MDP(22) - MDP(4) * t654 + (t537 + 0.2e1 * t558 + 0.2e1 * t559 + (qJD(1) * t496 - g(3)) * t568 + (qJD(1) * t486 - t625) * t570) * MDP(13) - t528 * MDP(23) * t671 + t673 * MDP(5) * t573; MDP(11) * t657 + (qJDD(2) + t654) * MDP(12) + (-t572 - t691) * MDP(13) + (qJD(2) * t510 + t526 + t601 - t700) * MDP(14) + (t564 * t657 - t563 * t691 + (t488 + t645) * qJD(2)) * MDP(15) + (-t563 * t657 - t564 * t691 + (-t489 - t646) * qJD(2)) * MDP(16) + (-t563 * t675 - t615 + (t488 * t564 + t489 * t563) * t672) * MDP(17) + (-qJD(2) * t478 + t697 - t720) * MDP(18) + (t394 * t603 - t437 * t676 - t722) * MDP(27) + (-qJD(2) * t388 + t718) * MDP(29) + t729 * (qJD(2) * t605 + t598) + t740 * (-qJD(2) * t437 + t622); ((t489 - t668) * t672 + t675) * MDP(15) + ((t488 + t662) * t672 + t635) * MDP(16) + (-t488 ^ 2 - t489 ^ 2) * MDP(17) + (t416 * t489 - t417 * t488 + (-g(3) - t630) * t568 + t580) * MDP(18) + (-t717 - t735) * MDP(27) + (t642 + t375 * t437 + t374 * t605 + (-g(3) + t599) * t568 + t580 + t723) * MDP(29) - t729 * (t393 + t734) + t740 * (t394 - t733); -MDP(19) * t694 + (t717 - t735) * MDP(20) + t381 * MDP(21) + (-t394 - t733) * MDP(22) + t491 * MDP(23) + (t444 * t605 + t584 + t699) * MDP(24) + (t378 * t528 + t437 * t444 - t575) * MDP(25) + (-t397 * t437 - t577 + t699 + 0.2e1 * t714) * MDP(26) + (pkin(5) * t393 - qJ(6) * t394 - (t375 - t379) * t605 + (t374 - t659) * t437) * MDP(27) + (0.2e1 * t701 - t388 * t437 - t397 * t605 + (0.2e1 * qJD(6) - t378) * t528 + t575) * MDP(28) + (t369 * qJ(6) - t370 * pkin(5) - t388 * t397 - t374 * t379 - g(1) * (-pkin(5) * t465 + qJ(6) * t466) - g(2) * (pkin(5) * t467 - qJ(6) * t468) - (-pkin(5) * t547 - qJ(6) * t546) * t556 + t659 * t375) * MDP(29); (-t491 - t694) * MDP(26) + t381 * MDP(27) + (-t528 ^ 2 - t717) * MDP(28) + (-t375 * t528 + t577 - t714) * MDP(29);];
tau  = t1;
