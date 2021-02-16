% Calculate vector of inverse dynamics joint torques for
% S6PRRPRP1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:39:08
% EndTime: 2021-01-16 02:39:31
% DurationCPUTime: 11.05s
% Computational Cost: add. (5707->574), mult. (13336->736), div. (0->0), fcn. (10499->14), ass. (0->256)
t562 = sin(qJ(3));
t563 = sin(qJ(2));
t557 = sin(pkin(6));
t656 = qJD(1) * t557;
t636 = t563 * t656;
t704 = qJD(3) * pkin(3);
t726 = -t562 * t704 + t636;
t565 = cos(qJ(3));
t560 = qJ(4) + pkin(8);
t626 = qJD(3) * t560;
t502 = qJD(4) * t565 - t562 * t626;
t555 = sin(pkin(11));
t585 = -qJD(4) * t562 - t565 * t626;
t702 = cos(pkin(11));
t440 = t502 * t702 + t555 * t585;
t620 = t702 * t565;
t682 = t555 * t562;
t593 = t620 - t682;
t566 = cos(qJ(2));
t635 = t566 * t656;
t478 = t593 * t635;
t733 = t440 - t478;
t621 = t702 * t562;
t518 = t555 * t565 + t621;
t509 = t518 * qJD(3);
t512 = t593 * qJD(3);
t734 = pkin(4) * t509 - pkin(9) * t512 - t726;
t510 = t518 * qJD(2);
t561 = sin(qJ(5));
t564 = cos(qJ(5));
t484 = qJD(3) * t561 + t510 * t564;
t691 = t484 * t561;
t558 = cos(pkin(6));
t678 = t557 * t563;
t513 = t558 * t565 - t562 * t678;
t703 = cos(pkin(10));
t624 = t703 * t563;
t556 = sin(pkin(10));
t679 = t556 * t566;
t505 = t558 * t624 + t679;
t552 = qJ(3) + pkin(11);
t547 = sin(t552);
t548 = cos(t552);
t625 = t557 * t703;
t461 = t505 * t547 + t548 * t625;
t623 = t703 * t566;
t680 = t556 * t563;
t503 = t558 * t680 - t623;
t681 = t556 * t557;
t463 = t503 * t547 + t548 * t681;
t496 = t547 * t678 - t558 * t548;
t587 = -g(1) * t463 + g(2) * t461 + g(3) * t496;
t732 = g(3) * t557;
t535 = qJD(2) * t620;
t654 = qJD(2) * t562;
t507 = t555 * t654 - t535;
t501 = qJD(5) + t507;
t731 = t501 * t691;
t645 = t558 * qJDD(1);
t534 = t565 * t645;
t650 = qJD(1) * qJD(2);
t494 = qJDD(2) * pkin(8) + (qJDD(1) * t563 + t566 * t650) * t557;
t655 = qJD(1) * t558;
t578 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t655 + t494;
t616 = t560 * qJD(2) + t636;
t603 = t616 * qJD(3);
t405 = qJDD(3) * pkin(3) - t562 * t578 - t565 * t603 + t534;
t406 = (-t603 + t645) * t562 + t578 * t565;
t382 = t555 * t405 + t702 * t406;
t380 = qJDD(3) * pkin(9) + t382;
t475 = -t562 * t616 + t565 * t655;
t471 = t475 + t704;
t476 = t562 * t655 + t565 * t616;
t622 = t702 * t476;
t415 = t555 * t471 + t622;
t410 = qJD(3) * pkin(9) + t415;
t546 = t565 * pkin(3) + pkin(2);
t500 = -qJD(2) * t546 + qJD(4) - t635;
t426 = pkin(4) * t507 - pkin(9) * t510 + t500;
t392 = t410 * t564 + t426 * t561;
t647 = qJDD(2) * t562;
t608 = -qJDD(2) * t620 + t555 * t647;
t459 = qJD(2) * t509 + t608;
t676 = t557 * t566;
t627 = qJDD(1) * t676;
t630 = t563 * t650;
t609 = t557 * t630 - t627;
t699 = qJDD(2) * pkin(2);
t493 = t609 - t699;
t646 = qJDD(2) * t565;
t649 = qJD(2) * qJD(3);
t629 = t562 * t649;
t722 = pkin(3) * t629 + qJDD(4);
t580 = qJDD(2) * t518 - t555 * t629;
t727 = qJD(3) * t535 + t580;
t395 = -pkin(3) * t646 + t459 * pkin(4) - pkin(9) * t727 + t493 + t722;
t394 = t564 * t395;
t576 = -qJD(5) * t392 - t380 * t561 + t394;
t648 = qJD(3) * qJD(5);
t652 = qJD(5) * t561;
t592 = t561 * qJDD(3) - t510 * t652 + (t727 + t648) * t564;
t701 = qJ(6) * t592;
t455 = qJDD(5) + t459;
t716 = pkin(5) * t455;
t370 = -qJD(6) * t484 + t576 - t701 + t716;
t482 = -t564 * qJD(3) + t510 * t561;
t384 = -qJ(6) * t482 + t392;
t698 = t384 * t501;
t730 = t370 + t698;
t381 = t702 * t405 - t555 * t406;
t379 = -qJDD(3) * pkin(4) - t381;
t570 = -t564 * qJDD(3) + t561 * t727;
t412 = t484 * qJD(5) + t570;
t374 = t412 * pkin(5) + qJDD(6) + t379;
t729 = -t374 + t587;
t718 = pkin(3) * t555;
t542 = pkin(9) + t718;
t670 = qJ(6) + t542;
t618 = qJD(5) * t670;
t467 = t555 * t476;
t419 = t475 * t702 - t467;
t644 = pkin(3) * t654;
t441 = pkin(4) * t510 + pkin(9) * t507 + t644;
t662 = t564 * t419 + t561 * t441;
t688 = t507 * t561;
t665 = qJ(6) * t688 - qJD(6) * t564 + t561 * t618 + t662;
t728 = t478 * t561 + t734 * t564;
t661 = t502 * t555 - t518 * t635 - t702 * t585;
t432 = t564 * t441;
t687 = t507 * t564;
t658 = -pkin(5) * t510 - qJ(6) * t687 - t564 * t618 - t432 + (-qJD(6) + t419) * t561;
t504 = -t558 * t623 + t680;
t506 = t558 * t679 + t624;
t611 = g(1) * t506 + g(2) * t504;
t583 = -g(3) * t676 + t611;
t457 = -pkin(4) * t593 - pkin(9) * t518 - t546;
t651 = qJD(5) * t564;
t725 = t457 * t651 + t734 * t561 + t733 * t564;
t724 = MDP(21) + MDP(23);
t723 = MDP(22) + MDP(24);
t460 = t503 * t548 - t547 * t681;
t462 = t505 * t548 - t547 * t625;
t497 = t547 * t558 + t548 * t678;
t671 = t564 * t566;
t640 = t557 * t671;
t721 = -g(3) * (-t497 * t561 - t640) - g(2) * (-t462 * t561 + t504 * t564) - g(1) * (t460 * t561 + t506 * t564);
t567 = qJD(3) ^ 2;
t720 = -pkin(8) * t567 + t557 * (-g(3) * t566 + t630) - t493 + t611 + t699;
t719 = t484 ^ 2;
t717 = pkin(3) * t562;
t707 = t482 * pkin(5);
t706 = t564 * pkin(5);
t705 = qJD(2) * pkin(2);
t700 = qJ(6) * t412;
t697 = t592 * t561;
t696 = t482 * t501;
t695 = t482 * t507;
t694 = t482 * t510;
t693 = t484 * t501;
t692 = t484 * t510;
t689 = t503 * t561;
t686 = t518 * t561;
t685 = t518 * t564;
t684 = t548 * t561;
t683 = t548 * t564;
t677 = t557 * t565;
t674 = t561 * t455;
t673 = t561 * t563;
t672 = t561 * t566;
t448 = t564 * t455;
t527 = t560 * t565;
t480 = t527 * t702 - t560 * t682;
t472 = t564 * t480;
t669 = qJDD(1) - g(3);
t605 = -qJ(6) * t512 - qJD(6) * t518;
t668 = pkin(5) * t509 - t440 * t561 + t605 * t564 + (-t472 + (qJ(6) * t518 - t457) * t561) * qJD(5) + t728;
t632 = t518 * t651;
t667 = -qJ(6) * t632 + (-qJD(5) * t480 + t605) * t561 + t725;
t391 = -t410 * t561 + t564 * t426;
t383 = -qJ(6) * t484 + t391;
t377 = pkin(5) * t501 + t383;
t666 = -t383 + t377;
t595 = t512 * t561 + t632;
t664 = pkin(5) * t595 + t661;
t663 = -t561 * t412 - t482 * t651;
t660 = t561 * t457 + t472;
t659 = -t503 * t560 - t506 * t546;
t553 = t562 ^ 2;
t657 = -t565 ^ 2 + t553;
t653 = qJD(2) * t563;
t642 = t557 * t672;
t638 = t587 * t561;
t637 = t702 * pkin(3);
t634 = t557 * t653;
t633 = qJD(2) * t676;
t631 = g(3) * (t546 * t676 + t560 * t678);
t628 = t565 * t649;
t619 = qJD(6) + t707;
t417 = t475 * t555 + t622;
t479 = t527 * t555 + t560 * t621;
t617 = t501 * t564;
t615 = t564 * t380 + t561 * t395 - t410 * t652 + t426 * t651;
t614 = t562 * t633;
t613 = t556 * pkin(3) * t677 + t503 * t717;
t543 = -t637 - pkin(4);
t612 = g(1) * t503 - g(2) * t505;
t414 = t471 * t702 - t467;
t545 = pkin(4) + t706;
t559 = -qJ(6) - pkin(9);
t607 = t545 * t548 - t547 * t559;
t371 = -qJD(6) * t482 + t615 - t700;
t606 = -t377 * t501 + t371;
t604 = t513 * pkin(3);
t568 = qJD(2) ^ 2;
t602 = qJDD(2) * t566 - t563 * t568;
t599 = t448 + (-t652 - t688) * t501;
t598 = -g(1) * t556 + g(2) * t703;
t514 = t558 * t562 + t563 * t677;
t447 = t555 * t513 + t514 * t702;
t427 = -t447 * t561 - t640;
t596 = -t447 * t564 + t642;
t594 = t512 * t564 - t518 * t652;
t409 = -qJD(3) * pkin(4) - t414;
t591 = t409 * t501 - t542 * t455;
t590 = -g(1) * (-t503 * t564 + t506 * t684) - g(2) * (t504 * t684 + t505 * t564) - (-t548 * t672 + t563 * t564) * t732;
t589 = -g(1) * (-t506 * t683 - t689) - g(2) * (-t504 * t683 + t505 * t561) - (t548 * t671 + t673) * t732;
t588 = -g(1) * t460 + g(2) * t462 + g(3) * t497;
t586 = t602 * t557;
t584 = t514 * qJD(3);
t524 = -t635 - t705;
t582 = -qJD(2) * t524 - t494 - t612;
t581 = (-t505 * t562 - t565 * t625) * pkin(3);
t577 = -g(1) * (t460 * t564 - t506 * t561) - g(2) * (-t462 * t564 - t504 * t561) - g(3) * (-t497 * t564 + t642) - t615;
t575 = -pkin(8) * qJDD(3) + (t524 + t635 - t705) * qJD(3);
t458 = -qJDD(2) * t546 + t609 + t722;
t574 = -t584 - t614;
t571 = t576 + t721;
t569 = t510 * t651 + t561 * t648 + t570;
t525 = t543 - t706;
t516 = t670 * t564;
t515 = t670 * t561;
t489 = t504 * t546;
t481 = t482 ^ 2;
t474 = qJD(3) * t513 + t565 * t633;
t450 = t564 * t457;
t446 = -t513 * t702 + t514 * t555;
t438 = pkin(5) * t686 + t479;
t418 = t474 * t702 + t555 * t574;
t416 = t474 * t555 - t574 * t702;
t399 = -pkin(5) * t688 + t417;
t398 = -qJ(6) * t686 + t660;
t397 = t409 + t619;
t396 = -pkin(5) * t593 - qJ(6) * t685 - t480 * t561 + t450;
t390 = qJD(5) * t596 - t418 * t561 + t564 * t634;
t389 = qJD(5) * t427 + t418 * t564 + t561 * t634;
t1 = [t669 * MDP(1) + MDP(3) * t586 + (-qJDD(2) * t563 - t566 * t568) * t557 * MDP(4) + (t513 * qJDD(3) + t565 * t586 + (-t584 - 0.2e1 * t614) * qJD(3)) * MDP(10) + (-qJD(3) * t474 - qJDD(3) * t514 + (-t562 * t602 - t566 * t628) * t557) * MDP(11) + (-qJD(3) * t416 - qJDD(3) * t446 + (-t459 * t566 + t507 * t653) * t557) * MDP(12) + (-t418 * qJD(3) - t447 * qJDD(3) + (t510 * t653 - t566 * t727) * t557) * MDP(13) + (t416 * t510 - t418 * t507 + t446 * t727 - t447 * t459) * MDP(14) + (-t381 * t446 + t382 * t447 - t414 * t416 + t415 * t418 - g(3) + (-t458 * t566 + t500 * t653) * t557) * MDP(15) + (-t389 * t482 - t390 * t484 + t412 * t596 - t427 * t592) * MDP(25) + (t370 * t427 - t371 * t596 + t374 * t446 + t377 * t390 + t384 * t389 + t397 * t416 - g(3)) * MDP(26) + t724 * (t390 * t501 + t412 * t446 + t416 * t482 + t427 * t455) + t723 * (-t389 * t501 + t416 * t484 + t446 * t592 + t455 * t596); (-t396 * t592 - t398 * t412 + (-t377 * t564 - t384 * t561) * t512 - t668 * t484 - t667 * t482 + (-t370 * t564 - t371 * t561 + (t377 * t561 - t384 * t564) * qJD(5)) * t518) * MDP(25) + (t484 * t594 + t592 * t685) * MDP(16) + (-t660 * t455 + t615 * t593 - t392 * t509 + t479 * t592 + t379 * t685 + (t480 * t652 - t725) * t501 + t661 * t484 + t594 * t409 + t590) * MDP(22) + (t371 * t593 + t374 * t685 - t384 * t509 + t397 * t594 - t398 * t455 + t438 * t592 + t484 * t664 - t501 * t667 + t590) * MDP(24) + (t448 * t518 + t484 * t509 + t501 * t594 - t592 * t593) * MDP(18) + (-t507 * t636 - qJDD(3) * t479 - t458 * t593 - t459 * t546 + t500 * t509 + (t507 * t717 - t661) * qJD(3)) * MDP(12) + (-t455 * t593 + t501 * t509) * MDP(20) + (-t370 * t593 + t374 * t686 + t377 * t509 + t396 * t455 + t397 * t595 + t412 * t438 + t482 * t664 + t501 * t668 + t589) * MDP(23) + (t412 * t593 - t482 * t509 - t501 * t595 - t518 * t674) * MDP(19) + (t450 * t455 - (-t410 * t651 + t394) * t593 + t391 * t509 + t479 * t412 + t409 * t632 + (-t480 * t651 + t728) * t501 + t661 * t482 + ((-qJD(5) * t457 - t440) * t501 - t480 * t455 - (-qJD(5) * t426 - t380) * t593 + t379 * t518 + t409 * t512) * t561 + t589) * MDP(21) + (t382 * t480 - t381 * t479 - t458 * t546 - g(1) * t659 - g(2) * (t505 * t560 - t489) - t631 - t726 * t500 + t733 * t415 - t661 * t414) * MDP(15) + (-t733 * qJD(3) - t480 * qJDD(3) + t458 * t518 + t500 * t512 - t546 * t727) * MDP(13) + (-g(3) * t678 - t381 * t518 + t382 * t593 - t414 * t512 - t415 * t509 - t480 * t459 + t479 * t727 - t733 * t507 + t612) * MDP(14) + (t583 + t627) * MDP(3) + ((-MDP(13) + MDP(25)) * t547 + t548 * MDP(12)) * t583 + (t575 * t562 + t565 * t720) * MDP(10) + (-t562 * t720 + t575 * t565) * MDP(11) + ((-t482 * t564 - t691) * t512 + (-t697 - t412 * t564 + (t482 * t561 - t484 * t564) * qJD(5)) * t518) * MDP(17) + (qJDD(3) * t562 + t565 * t567) * MDP(7) + (qJDD(3) * t565 - t562 * t567) * MDP(8) + qJDD(2) * MDP(2) + (-t669 * t678 - t612) * MDP(4) + 0.2e1 * (t562 * t646 - t649 * t657) * MDP(6) + (qJDD(2) * t553 + 0.2e1 * t562 * t628) * MDP(5) - (MDP(13) * t726 - MDP(14) * t661) * t510 + (t371 * t398 + t370 * t396 + t374 * t438 - g(1) * (-pkin(5) * t689 - t506 * t607 + t659) - g(2) * (-t489 + (pkin(5) * t561 + t560) * t505 - t607 * t504) - t631 - (pkin(5) * t673 + t566 * t607) * t732 + t664 * t397 + t667 * t384 + t668 * t377) * MDP(26); MDP(7) * t647 + MDP(8) * t646 + qJDD(3) * MDP(9) + (-g(3) * t513 + t562 * t582 + t598 * t677 + t534) * MDP(10) + (g(3) * t514 + (-t557 * t598 - t645) * t562 + t582 * t565) * MDP(11) + (t417 * qJD(3) - t500 * t510 + (qJDD(3) * t702 - t507 * t654) * pkin(3) + t587 + t381) * MDP(12) + (qJD(3) * t419 + t500 * t507 + (-qJDD(3) * t555 - t510 * t654) * pkin(3) + t588 - t382) * MDP(13) + (-t459 * t718 - t727 * t637 - (-t415 + t417) * t510 + (-t414 + t419) * t507) * MDP(14) + (-g(1) * t613 - g(2) * t581 - g(3) * t604 + t381 * t637 + t382 * t718 + t414 * t417 - t415 * t419 - t500 * t644) * MDP(15) + (t484 * t617 + t697) * MDP(16) + ((t592 - t695) * t564 - t731 + t663) * MDP(17) + (t501 * t617 + t674 - t692) * MDP(18) + (t599 + t694) * MDP(19) - t501 * t510 * MDP(20) + (-t391 * t510 + t543 * t412 - t417 * t482 - t432 * t501 + (t419 * t501 + t591) * t561 + (-qJD(5) * t501 * t542 - t379 + t587) * t564) * MDP(21) + (t379 * t561 + t392 * t510 + t543 * t592 - t417 * t484 + (t542 * t652 + t662) * t501 + t591 * t564 - t638) * MDP(22) + (-t377 * t510 - t399 * t482 + t412 * t525 - t455 * t515 + t658 * t501 + (t397 * t507 + (t397 + t707) * qJD(5)) * t561 + t729 * t564) * MDP(23) + (t397 * t687 + t374 * t561 + t384 * t510 - t399 * t484 + t592 * t525 - t455 * t516 + t665 * t501 + (pkin(5) * t691 + t397 * t564) * qJD(5) - t638) * MDP(24) + (-t412 * t516 + t665 * t482 - t658 * t484 + t515 * t592 - t561 * t730 + t606 * t564 - t588) * MDP(25) + (t371 * t516 - t370 * t515 + t374 * t525 - g(1) * (t460 * t559 + t463 * t545 + t613) - g(2) * (-t461 * t545 - t462 * t559 + t581) - g(3) * (-t496 * t545 - t497 * t559 + t604) + (pkin(5) * t652 - t399) * t397 - t665 * t384 + t658 * t377) * MDP(26) + (-MDP(5) * t562 * t565 + MDP(6) * t657) * t568; (0.2e1 * qJD(3) * t510 + t608) * MDP(12) + ((t535 - t507) * qJD(3) + t580) * MDP(13) + (-t507 ^ 2 - t510 ^ 2) * MDP(14) + (t414 * t510 + t415 * t507 + t458 - t583) * MDP(15) + ((-t592 - t695) * t564 + t731 + t663) * MDP(25) + (-t397 * t510 + t606 * t561 + t564 * t730 - t583) * MDP(26) + t723 * (-t501 ^ 2 * t564 - t674 - t692) + t724 * (t599 - t694); t484 * t482 * MDP(16) + (-t481 + t719) * MDP(17) + (t592 + t696) * MDP(18) + (-t569 + t693) * MDP(19) + t455 * MDP(20) + (t392 * t501 - t409 * t484 + t571) * MDP(21) + (t391 * t501 + t409 * t482 + t577) * MDP(22) + (0.2e1 * t716 - t701 + t698 + (-t397 - t619) * t484 + t571) * MDP(23) + (-pkin(5) * t719 + t700 + t383 * t501 + (qJD(6) + t397) * t482 + t577) * MDP(24) + (-pkin(5) * t592 - t482 * t666) * MDP(25) + (t666 * t384 + (-t397 * t484 + t370 + t721) * pkin(5)) * MDP(26); (t569 + t693) * MDP(23) + (t592 - t696) * MDP(24) + (-t481 - t719) * MDP(25) + (t377 * t484 + t384 * t482 - t729) * MDP(26);];
tau = t1;
