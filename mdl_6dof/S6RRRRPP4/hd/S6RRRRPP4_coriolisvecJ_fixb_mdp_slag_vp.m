% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRRPP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:02:52
% EndTime: 2019-03-09 21:03:07
% DurationCPUTime: 10.09s
% Computational Cost: add. (12011->572), mult. (29303->738), div. (0->0), fcn. (20438->8), ass. (0->255)
t616 = sin(qJ(3));
t617 = sin(qJ(2));
t686 = qJD(1) * t617;
t669 = t616 * t686;
t619 = cos(qJ(3));
t677 = t619 * qJD(2);
t559 = t669 - t677;
t618 = cos(qJ(4));
t684 = qJD(2) * t616;
t561 = t619 * t686 + t684;
t615 = sin(qJ(4));
t721 = t561 * t615;
t505 = t618 * t559 + t721;
t613 = sin(pkin(10));
t614 = cos(pkin(10));
t638 = t559 * t615 - t618 * t561;
t459 = t614 * t505 - t613 * t638;
t578 = -qJD(2) * pkin(2) + pkin(7) * t686;
t523 = pkin(3) * t559 + t578;
t474 = pkin(4) * t505 + qJD(5) + t523;
t737 = -t505 * t613 - t614 * t638;
t412 = pkin(5) * t459 - qJ(6) * t737 + t474;
t752 = t412 * MDP(29);
t562 = t615 * t616 - t618 * t619;
t730 = qJD(3) + qJD(4);
t514 = t730 * t562;
t620 = cos(qJ(2));
t630 = t562 * t620;
t527 = qJD(1) * t630;
t698 = -t514 + t527;
t563 = t615 * t619 + t616 * t618;
t515 = t730 * t563;
t685 = qJD(1) * t620;
t526 = t563 * t685;
t697 = t515 - t526;
t751 = pkin(5) * t737 + qJ(6) * t459;
t573 = -pkin(2) * t620 - pkin(8) * t617 - pkin(1);
t552 = t573 * qJD(1);
t606 = pkin(7) * t685;
t579 = qJD(2) * pkin(8) + t606;
t511 = t619 * t552 - t579 * t616;
t484 = -pkin(9) * t561 + t511;
t596 = -qJD(3) + t685;
t477 = -pkin(3) * t596 + t484;
t716 = t616 * t552;
t512 = t579 * t619 + t716;
t485 = -pkin(9) * t559 + t512;
t481 = t615 * t485;
t436 = t618 * t477 - t481;
t743 = qJ(5) * t638;
t425 = t436 + t743;
t584 = -qJD(4) + t596;
t422 = -pkin(4) * t584 + t425;
t483 = t618 * t485;
t437 = t477 * t615 + t483;
t744 = qJ(5) * t505;
t426 = t437 - t744;
t423 = t614 * t426;
t394 = t613 * t422 + t423;
t392 = -qJ(6) * t584 + t394;
t749 = t392 * t737;
t748 = t394 * t737;
t681 = qJD(3) * t616;
t645 = -t606 + (-t616 * t685 + t681) * pkin(3);
t747 = pkin(4) * t697 + t645;
t675 = qJD(1) * qJD(2);
t662 = t620 * t675;
t666 = t617 * t681;
t674 = qJD(2) * qJD(3);
t520 = -qJD(1) * t666 + (t662 + t674) * t619;
t680 = qJD(3) * t619;
t664 = t617 * t680;
t682 = qJD(2) * t620;
t668 = t616 * t682;
t627 = t664 + t668;
t521 = qJD(1) * t627 + t616 * t674;
t678 = qJD(4) * t618;
t671 = t618 * t520 - t615 * t521 - t559 * t678;
t679 = qJD(4) * t615;
t448 = -t561 * t679 + t671;
t644 = pkin(2) * t617 - pkin(8) * t620;
t570 = t644 * qJD(2);
t553 = qJD(1) * t570;
t663 = t617 * t675;
t648 = pkin(7) * t663;
t696 = -t619 * t553 - t616 * t648;
t626 = -qJD(3) * t512 - t696;
t444 = pkin(3) * t663 - pkin(9) * t520 + t626;
t633 = t552 * t680 + t616 * t553 - t579 * t681;
t624 = -t619 * t648 + t633;
t455 = -pkin(9) * t521 + t624;
t656 = t618 * t444 - t615 * t455;
t625 = -qJD(4) * t437 + t656;
t386 = pkin(4) * t663 - qJ(5) * t448 + qJD(5) * t638 + t625;
t652 = t520 * t615 + t618 * t521;
t449 = -qJD(4) * t638 + t652;
t647 = -t615 * t444 - t618 * t455 - t477 * t678 + t485 * t679;
t389 = -qJ(5) * t449 - qJD(5) * t505 - t647;
t377 = t614 * t386 - t613 * t389;
t376 = -pkin(5) * t663 - t377;
t729 = t412 * t737 + t376;
t746 = pkin(4) * t638;
t742 = t459 * t584;
t741 = t505 * t584;
t740 = t523 * t638;
t739 = t584 * t638;
t703 = t613 * t698 + t697 * t614;
t702 = -t697 * t613 + t614 * t698;
t738 = t505 * t523 + t647;
t659 = MDP(22) * t686;
t736 = qJD(2) * t659 + (-t505 ^ 2 + t638 ^ 2) * MDP(19) - t505 * MDP(18) * t638;
t735 = t737 ^ 2;
t734 = -0.2e1 * t675;
t733 = MDP(4) * t617;
t611 = t617 ^ 2;
t732 = MDP(5) * (-t620 ^ 2 + t611);
t537 = t563 * t617;
t711 = t619 * t620;
t634 = pkin(3) * t617 - pkin(9) * t711;
t567 = t644 * qJD(1);
t690 = pkin(7) * t669 + t619 * t567;
t492 = qJD(1) * t634 + t690;
t490 = t618 * t492;
t546 = t616 * t567;
t713 = t617 * t619;
t714 = t616 * t620;
t509 = t546 + (-pkin(7) * t713 - pkin(9) * t714) * qJD(1);
t438 = pkin(4) * t686 + qJ(5) * t527 - t509 * t615 + t490;
t728 = pkin(8) + pkin(9);
t670 = qJD(3) * t728;
t568 = t616 * t670;
t569 = t619 * t670;
t580 = t728 * t616;
t581 = t728 * t619;
t628 = -t618 * t568 - t615 * t569 - t580 * t678 - t581 * t679;
t439 = -qJ(5) * t515 - qJD(5) * t562 + t628;
t700 = t615 * t492 + t618 * t509;
t441 = -qJ(5) * t526 + t700;
t549 = t618 * t569;
t691 = -t615 * t580 + t618 * t581;
t623 = qJ(5) * t514 - qJD(4) * t691 - qJD(5) * t563 + t568 * t615 - t549;
t705 = (t438 - t623) * t614 + (t439 - t441) * t613;
t558 = t619 * t573;
t727 = pkin(7) * t616;
t510 = -pkin(9) * t713 + t558 + (-pkin(3) - t727) * t620;
t598 = pkin(7) * t711;
t689 = t616 * t573 + t598;
t715 = t616 * t617;
t516 = -pkin(9) * t715 + t689;
t699 = t615 * t510 + t618 * t516;
t701 = t618 * t484 - t481;
t429 = t701 + t743;
t654 = -t484 * t615 - t483;
t632 = t654 + t744;
t673 = pkin(3) * t613 * t615;
t694 = -qJD(4) * t673 - t613 * t632 + (pkin(3) * t678 - t429) * t614;
t731 = t620 * t677 - t666;
t399 = t425 * t613 + t423;
t726 = t399 * t737;
t725 = t426 * t613;
t724 = t520 * t616;
t723 = t559 * t596;
t722 = t561 * t596;
t720 = t578 * t616;
t719 = t578 * t619;
t718 = t596 * t619;
t717 = t614 * t615;
t621 = qJD(2) ^ 2;
t712 = t617 * t621;
t710 = t620 * t621;
t622 = qJD(1) ^ 2;
t709 = t620 * t622;
t378 = t613 * t386 + t614 * t389;
t472 = -qJD(2) * t630 - t537 * t730;
t538 = t562 * t617;
t683 = qJD(2) * t617;
t692 = t619 * t570 + t683 * t727;
t466 = t634 * qJD(2) + (-t598 + (pkin(9) * t617 - t573) * t616) * qJD(3) + t692;
t665 = t620 * t681;
t693 = t616 * t570 + t573 * t680;
t471 = -t627 * pkin(9) + (-t617 * t677 - t665) * pkin(7) + t693;
t655 = t618 * t466 - t471 * t615;
t396 = pkin(4) * t683 - qJ(5) * t472 - qJD(4) * t699 + qJD(5) * t538 + t655;
t473 = -t679 * t715 + (t713 * t730 + t668) * t618 + t731 * t615;
t629 = t615 * t466 + t618 * t471 + t510 * t678 - t516 * t679;
t401 = -qJ(5) * t473 - qJD(5) * t537 + t629;
t382 = t613 * t396 + t614 * t401;
t708 = -qJD(6) - t694;
t410 = t613 * t438 + t614 * t441;
t405 = qJ(6) * t686 + t410;
t408 = t614 * t439 + t613 * t623;
t707 = t405 - t408;
t706 = pkin(5) * t686 + t705;
t503 = -t562 * t613 + t563 * t614;
t704 = pkin(5) * t703 - qJ(6) * t702 - qJD(6) * t503 + t747;
t653 = t618 * t510 - t516 * t615;
t450 = -pkin(4) * t620 + qJ(5) * t538 + t653;
t456 = -qJ(5) * t537 + t699;
t420 = t613 * t450 + t614 * t456;
t695 = t429 * t613 - t614 * t632 - (t613 * t618 + t717) * qJD(4) * pkin(3);
t603 = pkin(3) * t618 + pkin(4);
t541 = pkin(3) * t717 + t613 * t603;
t571 = pkin(3) * t715 + t617 * pkin(7);
t400 = t425 * t614 - t725;
t676 = qJD(6) - t400;
t672 = qJ(6) * t663 + t378;
t524 = pkin(3) * t627 + pkin(7) * t682;
t604 = -pkin(3) * t619 - pkin(2);
t660 = MDP(15) * t686;
t501 = pkin(3) * t521 + pkin(7) * t662;
t658 = t695 * t737;
t657 = pkin(1) * t734;
t415 = t448 * t613 + t614 * t449;
t651 = -t618 * t580 - t581 * t615;
t650 = t559 + t677;
t649 = -t561 + t684;
t572 = t584 * qJD(6);
t375 = -t572 + t672;
t646 = pkin(4) * t537 + t571;
t643 = pkin(3) * t561 - t746;
t416 = t448 * t614 - t449 * t613;
t491 = -qJ(5) * t562 + t691;
t631 = -qJ(5) * t563 + t651;
t453 = t491 * t613 - t614 * t631;
t454 = t614 * t491 + t613 * t631;
t642 = -t408 * t459 - t454 * t415 + t416 * t453;
t381 = t396 * t614 - t401 * t613;
t393 = t422 * t614 - t725;
t419 = t450 * t614 - t456 * t613;
t637 = qJD(1) * t611 - t596 * t620;
t636 = pkin(4) * t562 + t604;
t635 = pkin(4) * t473 + t524;
t432 = pkin(4) * t449 + t501;
t540 = t603 * t614 - t673;
t383 = pkin(5) * t415 - qJ(6) * t416 - qJD(6) * t737 + t432;
t600 = -pkin(4) * t614 - pkin(5);
t599 = pkin(4) * t613 + qJ(6);
t531 = -pkin(5) - t540;
t530 = qJ(6) + t541;
t502 = t614 * t562 + t563 * t613;
t487 = -t537 * t613 - t538 * t614;
t486 = t614 * t537 - t538 * t613;
t452 = pkin(5) * t502 - qJ(6) * t503 + t636;
t433 = pkin(5) * t486 - qJ(6) * t487 + t646;
t431 = t472 * t614 - t473 * t613;
t430 = t472 * t613 + t614 * t473;
t418 = -t746 + t751;
t417 = pkin(5) * t620 - t419;
t414 = -qJ(6) * t620 + t420;
t413 = t643 + t751;
t391 = pkin(5) * t584 + qJD(6) - t393;
t390 = pkin(5) * t430 - qJ(6) * t431 - qJD(6) * t487 + t635;
t380 = -pkin(5) * t683 - t381;
t379 = qJ(6) * t683 - qJD(6) * t620 + t382;
t1 = [MDP(6) * t710 + (t520 * t713 + t561 * t731) * MDP(11) + ((-qJD(1) * t699 - t437) * MDP(24) + (-qJD(1) * t538 - t638) * MDP(20) + (-qJD(1) * t537 - t505) * MDP(21) + (qJD(1) * t653 + t436) * MDP(23) + (qJD(1) * t414 + t392) * MDP(29) + (-qJD(1) * t417 - t391) * MDP(27) + (-t596 - t685) * MDP(15) + (-t584 - t685) * MDP(22)) * t683 + (t571 * t448 + t523 * t472 - t501 * t538 - t524 * t638 + t629 * t584 - t647 * t620) * MDP(24) + (-t448 * t537 + t449 * t538 - t472 * t505 + t473 * t638) * MDP(19) + (-t448 * t538 - t472 * t638) * MDP(18) + 0.2e1 * t662 * t733 + (-t375 * t620 - t379 * t584 - t383 * t487 - t390 * t737 - t412 * t431 - t416 * t433) * MDP(29) + (-t655 * t584 - t656 * t620 + t524 * t505 + t571 * t449 + t501 * t537 + t523 * t473 + (t437 * t620 + t584 * t699) * qJD(4)) * MDP(23) + (-t375 * t486 + t376 * t487 - t379 * t459 + t380 * t737 + t391 * t431 - t392 * t430 - t414 * t415 + t416 * t417) * MDP(28) + (-t377 * t487 - t378 * t486 - t381 * t737 - t382 * t459 - t393 * t431 - t394 * t430 - t415 * t420 - t416 * t419) * MDP(25) + t732 * t734 + ((-t559 * t619 - t561 * t616) * t682 + (-t724 - t521 * t619 + (t559 * t616 - t561 * t619) * qJD(3)) * t617) * MDP(12) + ((-pkin(7) * t665 + t693) * t596 + t633 * t620 + (pkin(7) * t520 - t578 * t681) * t617 + ((pkin(7) * t561 + t719) * t620 + (-pkin(7) * t718 - qJD(1) * t689 - t512) * t617) * qJD(2)) * MDP(17) + (-(-t573 * t681 + t692) * t596 + (t578 * t680 + pkin(7) * t521 + (t558 * qJD(1) + t511) * qJD(2)) * t617 + ((pkin(7) * t559 + t720) * qJD(2) + (t716 + (pkin(7) * t596 + t579) * t619) * qJD(3) + t696) * t620) * MDP(16) + (pkin(7) * t712 + t620 * t657) * MDP(10) - MDP(7) * t712 + (-pkin(7) * t710 + t617 * t657) * MDP(9) + (t375 * t414 + t376 * t417 + t379 * t392 + t380 * t391 + t383 * t433 + t390 * t412) * MDP(30) + (t596 * t666 - t520 * t620 + (t561 * t617 + t619 * t637) * qJD(2)) * MDP(13) + (t596 * t664 + t521 * t620 + (-t559 * t617 - t616 * t637) * qJD(2)) * MDP(14) + (t377 * t419 + t378 * t420 + t393 * t381 + t394 * t382 + t432 * t646 + t474 * t635) * MDP(26) + (t449 * t620 + t473 * t584) * MDP(21) + (-t448 * t620 - t472 * t584) * MDP(20) + (t376 * t620 + t380 * t584 + t383 * t486 + t390 * t459 + t412 * t430 + t415 * t433) * MDP(27); t622 * t732 + (-t561 * t718 + t724) * MDP(11) + ((t520 + t723) * t619 + (-t521 + t722) * t616) * MDP(12) + (-t596 * t680 + (t596 * t711 + t617 * t649) * qJD(1)) * MDP(13) + (t596 * t681 + (-t596 * t714 + t617 * t650) * qJD(1)) * MDP(14) + (-pkin(2) * t521 + t690 * t596 + (pkin(8) * t718 + t720) * qJD(3) + ((-pkin(8) * t684 - t511) * t617 + (-pkin(7) * t650 - t720) * t620) * qJD(1)) * MDP(16) + (-pkin(2) * t520 - t546 * t596 + (-pkin(8) * t596 * t616 + t719) * qJD(3) + (-t578 * t711 + (-pkin(8) * t677 + t512) * t617 + (t596 * t713 + t620 * t649) * pkin(7)) * qJD(1)) * MDP(17) + (t448 * t563 - t638 * t698) * MDP(18) + (-t448 * t562 - t449 * t563 - t505 * t698 + t638 * t697) * MDP(19) + (-t698 * t584 + (qJD(2) * t563 + t638) * t686) * MDP(20) + (t697 * t584 + (-qJD(2) * t562 + t505) * t686) * MDP(21) + (t604 * t449 + t501 * t562 + (t581 * t678 + t490 + t549 + (-qJD(4) * t580 - t509 - t568) * t615) * t584 + t697 * t523 + t645 * t505 + (qJD(2) * t651 - t436) * t686) * MDP(23) + (t604 * t448 + t501 * t563 + (t628 - t700) * t584 + t698 * t523 - t645 * t638 + (-qJD(2) * t691 + t437) * t686) * MDP(24) + (-t377 * t503 - t378 * t502 - t393 * t702 - t394 * t703 + t410 * t459 + t705 * t737 + t642) * MDP(25) + (t378 * t454 - t377 * t453 + t432 * t636 + t747 * t474 + (t408 - t410) * t394 - t705 * t393) * MDP(26) + (t383 * t502 + t415 * t452 + t706 * t584 + t704 * t459 + t703 * t412 + (-qJD(2) * t453 + t391) * t686) * MDP(27) + (-t375 * t502 + t376 * t503 + t391 * t702 - t392 * t703 + t405 * t459 + t706 * t737 + t642) * MDP(28) + (-t383 * t503 - t416 * t452 + t707 * t584 - t704 * t737 - t702 * t412 + (qJD(2) * t454 - t392) * t686) * MDP(29) + (t375 * t454 + t376 * t453 + t383 * t452 + t391 * t706 - t392 * t707 + t412 * t704) * MDP(30) - t709 * t733 + t596 * t660 + t584 * t659 + (MDP(9) * t617 * t622 + MDP(10) * t709) * pkin(1); t561 * t559 * MDP(11) + (-t559 ^ 2 + t561 ^ 2) * MDP(12) + (t520 - t723) * MDP(13) + (-t521 - t722) * MDP(14) + qJD(2) * t660 + (-t512 * t596 - t561 * t578 + t626) * MDP(16) + (-t511 * t596 + t559 * t578 - t624) * MDP(17) + (t448 - t741) * MDP(20) + (-t449 + t739) * MDP(21) + (t654 * t584 + t740 + (-t505 * t561 + t584 * t679 + t618 * t663) * pkin(3) + t625) * MDP(23) + (-t701 * t584 + (t561 * t638 + t584 * t678 - t615 * t663) * pkin(3) + t738) * MDP(24) + (-t415 * t541 - t416 * t540 - t658 + t748) * MDP(25) + (t377 * t540 + t378 * t541 + t393 * t695 + t394 * t694 - t474 * t643) * MDP(26) + (-t531 * t663 - t584 * t695 - t729) * MDP(27) + (-t415 * t530 + t416 * t531 - t658 + t749) * MDP(28) + (t413 * t737 + t530 * t663 + t584 * t708 + t375) * MDP(29) + (t375 * t530 + t376 * t531 - t391 * t695 - t392 * t708 - t412 * t413) * MDP(30) + t736 + ((-t393 - t694) * MDP(25) - t413 * MDP(27) + (t391 + t708) * MDP(28) - t752) * t459; (t671 - t741) * MDP(20) + (-t652 + t739) * MDP(21) + (-t437 * t584 + t656 + t740) * MDP(23) + (-t436 * t584 + t738) * MDP(24) + (-t726 + t748) * MDP(25) + (t393 * t399 - t394 * t400) * MDP(26) + (-t399 * t584 - t600 * t663 - t729) * MDP(27) + (-t415 * t599 + t416 * t600 - t726 + t749) * MDP(28) + (t400 * t584 + t418 * t737 + t599 * t663 - 0.2e1 * t572 + t672) * MDP(29) + (t375 * t599 + t376 * t600 - t391 * t399 + t392 * t676 - t412 * t418) * MDP(30) + (-MDP(20) * t721 + MDP(21) * t638 - MDP(23) * t437) * qJD(4) + ((-t415 * t613 - t416 * t614) * MDP(25) + (t377 * t614 + t378 * t613 + t474 * t638) * MDP(26)) * pkin(4) + t736 + ((-t393 + t400) * MDP(25) - t418 * MDP(27) + (t391 - t676) * MDP(28) - t752) * t459; (t393 * t737 + t394 * t459 + t432) * MDP(26) + (-t584 * t737 + t415) * MDP(27) + (-t416 - t742) * MDP(29) + (-t391 * t737 + t392 * t459 + t383) * MDP(30) + (MDP(25) + MDP(28)) * (-t459 ^ 2 - t735); (t459 * t737 - t663) * MDP(27) + (t416 - t742) * MDP(28) + (-t584 ^ 2 - t735) * MDP(29) + (t392 * t584 + t729) * MDP(30);];
tauc  = t1;
