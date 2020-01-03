% Calculate vector of inverse dynamics joint torques for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR9_invdynJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:30:10
% EndTime: 2019-12-31 22:30:31
% DurationCPUTime: 12.26s
% Computational Cost: add. (5820->555), mult. (13069->740), div. (0->0), fcn. (9569->14), ass. (0->231)
t581 = sin(qJ(3));
t582 = sin(qJ(2));
t665 = qJD(1) * t582;
t644 = t581 * t665;
t586 = cos(qJ(3));
t651 = t586 * qJD(2);
t521 = t644 - t651;
t661 = qJD(2) * t581;
t523 = t586 * t665 + t661;
t580 = sin(qJ(4));
t585 = cos(qJ(4));
t455 = t585 * t521 + t523 * t580;
t584 = cos(qJ(5));
t579 = sin(qJ(5));
t611 = t521 * t580 - t585 * t523;
t696 = t611 * t579;
t407 = -t584 * t455 + t696;
t587 = cos(qJ(2));
t570 = t587 * qJDD(1);
t650 = qJD(1) * qJD(2);
t707 = -t582 * t650 + t570;
t518 = qJDD(3) - t707;
t513 = qJDD(4) + t518;
t498 = qJDD(5) + t513;
t613 = t455 * t579 + t584 * t611;
t730 = t498 * MDP(29) + (-t407 ^ 2 + t613 ^ 2) * MDP(26) + t407 * MDP(25) * t613;
t635 = t587 * t650;
t649 = qJDD(1) * t582;
t657 = qJD(3) * t582;
t717 = -qJD(1) * t657 + qJDD(2);
t443 = qJD(3) * t651 + (t635 + t649) * t586 + t717 * t581;
t664 = qJD(1) * t587;
t444 = (qJD(2) * (qJD(3) + t664) + t649) * t581 - t717 * t586;
t654 = qJD(4) * t585;
t655 = qJD(4) * t580;
t389 = t585 * t443 - t580 * t444 - t521 * t654 - t523 * t655;
t592 = qJD(4) * t611 - t443 * t580 - t585 * t444;
t652 = qJD(5) * t584;
t646 = t584 * t389 - t455 * t652 + t579 * t592;
t653 = qJD(5) * t579;
t366 = t611 * t653 + t646;
t555 = -qJD(3) + t664;
t547 = -qJD(4) + t555;
t631 = t389 * t579 - t584 * t592;
t593 = qJD(5) * t613 - t631;
t539 = -qJD(5) + t547;
t721 = t539 * t613;
t722 = t407 * t539;
t729 = t513 * MDP(22) + (-t455 ^ 2 + t611 ^ 2) * MDP(19) + (-t455 * t547 + t389) * MDP(20) + (t547 * t611 + t592) * MDP(21) - t455 * MDP(18) * t611 + (t593 + t721) * MDP(28) + (t366 + t722) * MDP(27) + t730;
t534 = -pkin(2) * t587 - pkin(7) * t582 - pkin(1);
t514 = t534 * qJD(1);
t566 = pkin(6) * t664;
t541 = qJD(2) * pkin(7) + t566;
t465 = t586 * t514 - t541 * t581;
t429 = -pkin(8) * t523 + t465;
t421 = -pkin(3) * t555 + t429;
t466 = t514 * t581 + t541 * t586;
t430 = -pkin(8) * t521 + t466;
t426 = t585 * t430;
t392 = t421 * t580 + t426;
t724 = pkin(9) * t455;
t380 = t392 - t724;
t376 = t380 * t653;
t540 = -qJD(2) * pkin(2) + pkin(6) * t665;
t474 = pkin(3) * t521 + t540;
t418 = pkin(4) * t455 + t474;
t578 = qJ(3) + qJ(4);
t574 = qJ(5) + t578;
t560 = sin(t574);
t561 = cos(t574);
t588 = cos(qJ(1));
t583 = sin(qJ(1));
t687 = t583 * t587;
t477 = t560 * t588 - t561 * t687;
t683 = t587 * t588;
t479 = t560 * t583 + t561 * t683;
t700 = g(3) * t582;
t716 = g(1) * t479 - g(2) * t477 - t418 * t407 + t561 * t700 + t376;
t476 = t560 * t687 + t561 * t588;
t478 = -t560 * t683 + t561 * t583;
t620 = pkin(2) * t582 - pkin(7) * t587;
t532 = t620 * qJD(2);
t471 = qJD(1) * t532 + qJDD(1) * t534;
t460 = t586 * t471;
t496 = t707 * pkin(6) + qJDD(2) * pkin(7);
t381 = pkin(3) * t518 - pkin(8) * t443 - qJD(3) * t466 - t496 * t581 + t460;
t656 = qJD(3) * t586;
t658 = qJD(3) * t581;
t603 = t581 * t471 + t586 * t496 + t514 * t656 - t541 * t658;
t385 = -pkin(8) * t444 + t603;
t632 = t585 * t381 - t580 * t385;
t594 = -t392 * qJD(4) + t632;
t364 = pkin(4) * t513 - pkin(9) * t389 + t594;
t622 = -t580 * t381 - t585 * t385 - t421 * t654 + t430 * t655;
t365 = pkin(9) * t592 - t622;
t633 = t584 * t364 - t579 * t365;
t715 = -g(1) * t478 + g(2) * t476 + t418 * t613 + t560 * t700 + t633;
t684 = t586 * t587;
t610 = pkin(3) * t582 - pkin(8) * t684;
t703 = pkin(7) + pkin(8);
t645 = qJD(3) * t703;
t529 = t620 * qJD(1);
t670 = pkin(6) * t644 + t586 * t529;
t726 = qJD(1) * t610 + t586 * t645 + t670;
t507 = t581 * t529;
t688 = t582 * t586;
t689 = t581 * t587;
t719 = -t507 - (-pkin(6) * t688 - pkin(8) * t689) * qJD(1) - t581 * t645;
t564 = pkin(6) * t649;
t497 = -qJDD(2) * pkin(2) + pkin(6) * t635 + t564;
t619 = g(1) * t588 + g(2) * t583;
t699 = g(3) * t587;
t599 = t582 * t619 - t699;
t725 = qJD(3) * pkin(7) * t555 - t497 + t599;
t723 = pkin(9) * t611;
t524 = t580 * t581 - t585 * t586;
t604 = t524 * t587;
t705 = qJD(3) + qJD(4);
t676 = qJD(1) * t604 - t705 * t524;
t525 = t580 * t586 + t581 * t585;
t675 = (-t664 + t705) * t525;
t659 = qJD(2) * t587;
t642 = t581 * t659;
t718 = t582 * t656 + t642;
t572 = sin(t578);
t573 = cos(t578);
t487 = t572 * t588 - t573 * t687;
t489 = t572 * t583 + t573 * t683;
t714 = g(1) * t489 - g(2) * t487 + t455 * t474 + t573 * t700 + t622;
t486 = t572 * t687 + t573 * t588;
t488 = -t572 * t683 + t573 * t583;
t713 = -g(1) * t488 + g(2) * t486 + t474 * t611 + t572 * t700 + t594;
t492 = t525 * t582;
t710 = t726 * t585;
t520 = t586 * t534;
t701 = pkin(6) * t581;
t464 = -pkin(8) * t688 + t520 + (-pkin(3) - t701) * t587;
t557 = pkin(6) * t684;
t669 = t581 * t534 + t557;
t690 = t581 * t582;
t472 = -pkin(8) * t690 + t669;
t677 = t580 * t464 + t585 * t472;
t621 = -t566 + (-t581 * t664 + t658) * pkin(3);
t542 = t703 * t581;
t543 = t703 * t586;
t671 = -t580 * t542 + t585 * t543;
t708 = -t542 * t654 - t543 * t655 - t726 * t580 + t719 * t585;
t706 = -t581 * t657 + t587 * t651;
t702 = pkin(3) * t580;
t424 = t580 * t430;
t391 = t585 * t421 - t424;
t379 = t391 + t723;
t375 = -pkin(4) * t547 + t379;
t698 = t375 * t584;
t697 = t443 * t581;
t695 = t521 * t555;
t694 = t523 * t555;
t693 = t523 * t586;
t692 = t579 * t498;
t691 = t580 * t584;
t686 = t584 * t380;
t685 = t584 * t498;
t461 = t584 * t524 + t525 * t579;
t682 = -qJD(5) * t461 - t675 * t579 + t676 * t584;
t462 = -t524 * t579 + t525 * t584;
t681 = qJD(5) * t462 + t676 * t579 + t675 * t584;
t680 = t585 * t429 - t424;
t679 = t675 * pkin(4) + t621;
t673 = t581 * t532 + t534 * t656;
t660 = qJD(2) * t582;
t672 = t586 * t532 + t660 * t701;
t533 = pkin(3) * t690 + t582 * pkin(6);
t576 = t582 ^ 2;
t668 = -t587 ^ 2 + t576;
t663 = qJD(2) * t521;
t662 = qJD(2) * t523;
t475 = t718 * pkin(3) + pkin(6) * t659;
t563 = -pkin(3) * t586 - pkin(2);
t643 = t555 * t651;
t640 = t555 * t658;
t639 = t555 * t656;
t411 = t610 * qJD(2) + (-t557 + (pkin(8) * t582 - t534) * t581) * qJD(3) + t672;
t414 = -t718 * pkin(8) + (-t582 * t651 - t587 * t658) * pkin(6) + t673;
t630 = t585 * t411 - t414 * t580;
t629 = -t429 * t580 - t426;
t627 = t585 * t464 - t472 * t580;
t626 = -t585 * t542 - t543 * t580;
t625 = -qJD(3) * t514 - t496;
t624 = qJD(5) * t375 + t365;
t618 = g(1) * t583 - g(2) * t588;
t617 = t541 * t656 - t460;
t439 = -pkin(9) * t524 + t671;
t616 = pkin(4) * t665 + t676 * pkin(9) + t671 * qJD(4) + qJD(5) * t439 + t719 * t580 + t710;
t438 = -pkin(9) * t525 + t626;
t615 = -t675 * pkin(9) + qJD(5) * t438 + t708;
t614 = -pkin(7) * t518 + qJD(3) * t540;
t369 = t579 * t375 + t686;
t493 = t524 * t582;
t431 = t584 * t492 - t493 * t579;
t432 = -t492 * t579 - t493 * t584;
t608 = -0.2e1 * pkin(1) * t650 - pkin(6) * qJDD(2);
t607 = t518 * t581 - t639;
t606 = t518 * t586 + t640;
t602 = t580 * t411 + t585 * t414 + t464 * t654 - t472 * t655;
t590 = qJD(1) ^ 2;
t600 = pkin(1) * t590 + t619;
t420 = pkin(3) * t444 + t497;
t589 = qJD(2) ^ 2;
t596 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t589 + t618;
t562 = pkin(3) * t585 + pkin(4);
t505 = t581 * t583 + t586 * t683;
t504 = -t581 * t683 + t583 * t586;
t503 = t581 * t588 - t583 * t684;
t502 = t581 * t687 + t586 * t588;
t483 = pkin(4) * t524 + t563;
t467 = pkin(4) * t492 + t533;
t423 = pkin(3) * t523 - pkin(4) * t611;
t416 = -t655 * t690 + (t705 * t688 + t642) * t585 + t706 * t580;
t415 = -qJD(2) * t604 - t705 * t492;
t401 = pkin(4) * t416 + t475;
t400 = -pkin(9) * t492 + t677;
t399 = -pkin(4) * t587 + pkin(9) * t493 + t627;
t383 = t680 + t723;
t382 = t629 + t724;
t374 = -pkin(4) * t592 + t420;
t373 = qJD(5) * t432 + t415 * t579 + t584 * t416;
t372 = -qJD(5) * t431 + t415 * t584 - t416 * t579;
t371 = -pkin(9) * t416 + t602;
t370 = pkin(4) * t660 - pkin(9) * t415 - qJD(4) * t677 + t630;
t368 = -t380 * t579 + t698;
t1 = [(-t389 * t493 - t415 * t611) * MDP(18) + (-g(1) * t486 - g(2) * t488 + t533 * t389 - t392 * t660 + t474 * t415 - t420 * t493 - t475 * t611 - t513 * t677 + t547 * t602 - t587 * t622) * MDP(24) + (-t389 * t587 - t415 * t547 - t493 * t513 - t611 * t660) * MDP(20) + (t366 * t432 - t372 * t613) * MDP(25) + (-t366 * t587 - t372 * t539 + t432 * t498 - t613 * t660) * MDP(27) + (-t369 * t660 - g(1) * t476 - g(2) * t478 + t467 * t366 + t418 * t372 + t374 * t432 - t376 * t587 - t401 * t613 + ((-qJD(5) * t400 + t370) * t539 - t399 * t498 + t364 * t587) * t579 + ((qJD(5) * t399 + t371) * t539 - t400 * t498 + t624 * t587) * t584) * MDP(31) + (t443 * t688 + t706 * t523) * MDP(11) + (qJDD(2) * t582 + t587 * t589) * MDP(6) + (qJDD(2) * t587 - t582 * t589) * MDP(7) + (-t389 * t492 - t415 * t455 + t416 * t611 - t493 * t592) * MDP(19) + (t416 * t547 - t455 * t660 - t492 * t513 - t587 * t592) * MDP(21) + (-t630 * t547 + t627 * t513 - t632 * t587 + t391 * t660 + t475 * t455 - t533 * t592 + t420 * t492 + t474 * t416 - g(1) * t487 - g(2) * t489 + (t392 * t587 + t547 * t677) * qJD(4)) * MDP(23) + (-t366 * t431 + t372 * t407 + t373 * t613 + t432 * t593) * MDP(26) + (t373 * t539 + t407 * t660 - t431 * t498 - t587 * t593) * MDP(28) + (-(t370 * t584 - t371 * t579) * t539 + (t399 * t584 - t400 * t579) * t498 - t633 * t587 + t368 * t660 - t401 * t407 - t467 * t593 + t374 * t431 + t418 * t373 - g(1) * t477 - g(2) * t479 + (-(-t399 * t579 - t400 * t584) * t539 + t369 * t587) * qJD(5)) * MDP(30) + qJDD(1) * MDP(1) + 0.2e1 * (t570 * t582 - t650 * t668) * MDP(5) + (-(-t534 * t658 + t672) * t555 + t520 * t518 - g(1) * t503 - g(2) * t505 + ((t639 + t663) * pkin(6) + (-pkin(6) * t518 + qJD(2) * t540 - t625) * t581 + t617) * t587 + (pkin(6) * t444 + qJD(2) * t465 + t497 * t581 + t540 * t656) * t582) * MDP(16) + (t673 * t555 - t669 * t518 - g(1) * t502 - g(2) * t504 + (t540 * t651 + (-t640 + t662) * pkin(6) + t603) * t587 + (-t540 * t658 - t466 * qJD(2) + t497 * t586 + (t443 - t643) * pkin(6)) * t582) * MDP(17) + (-t513 * t587 - t547 * t660) * MDP(22) + (-t498 * t587 - t539 * t660) * MDP(29) + ((-t443 - t643) * t587 + (t606 + t662) * t582) * MDP(13) + ((t555 * t661 + t444) * t587 + (-t607 - t663) * t582) * MDP(14) + (-t518 * t587 - t555 * t660) * MDP(15) + (qJDD(1) * t576 + 0.2e1 * t582 * t635) * MDP(4) + ((-t521 * t586 - t523 * t581) * t659 + (-t697 - t444 * t586 + (t521 * t581 - t693) * qJD(3)) * t582) * MDP(12) + (t582 * t608 + t587 * t596) * MDP(9) + (-t582 * t596 + t587 * t608) * MDP(10) + t618 * MDP(2) + t619 * MDP(3); (t366 * t462 - t613 * t682) * MDP(25) + (t389 * t525 - t611 * t676) * MDP(18) + (-pkin(2) * t444 + t670 * t555 + t614 * t581 + (-t465 * t582 + (-pkin(6) * t521 - t540 * t581) * t587) * qJD(1) + t725 * t586) * MDP(16) + (-pkin(2) * t443 - t507 * t555 + t614 * t586 + (-t540 * t684 + t466 * t582 + (-t523 * t587 + t555 * t688) * pkin(6)) * qJD(1) - t725 * t581) * MDP(17) + (-MDP(4) * t582 * t587 + MDP(5) * t668) * t590 + (t462 * t498 - t539 * t682) * MDP(27) + (t513 * t525 - t676 * t547) * MDP(20) + (-t389 * t524 - t455 * t676 + t525 * t592 + t611 * t675) * MDP(19) + (t626 * t513 - t563 * t592 + t420 * t524 + (t543 * t654 + (-qJD(4) * t542 + t719) * t580 + t710) * t547 + t675 * t474 + t621 * t455 + t599 * t573) * MDP(23) + (-(t438 * t579 + t439 * t584) * t498 + t483 * t366 + t374 * t462 + (-t579 * t616 + t584 * t615) * t539 + t682 * t418 - t679 * t613 - t599 * t560) * MDP(31) + (t563 * t389 + t420 * t525 + t676 * t474 - t671 * t513 + t708 * t547 - t572 * t599 - t611 * t621) * MDP(24) + (-t513 * t524 + t547 * t675) * MDP(21) + (t555 * MDP(15) + MDP(20) * t611 + t455 * MDP(21) + t547 * MDP(22) - t391 * MDP(23) + t392 * MDP(24) + MDP(27) * t613 - MDP(28) * t407 + t539 * MDP(29) - t368 * MDP(30) + t369 * MDP(31)) * t665 + (-t366 * t461 + t407 * t682 + t462 * t593 + t613 * t681) * MDP(26) + ((t438 * t584 - t439 * t579) * t498 - t483 * t593 + t374 * t461 + (t579 * t615 + t584 * t616) * t539 + t681 * t418 - t679 * t407 + t599 * t561) * MDP(30) + qJDD(2) * MDP(8) + (-t461 * t498 + t539 * t681) * MDP(28) + ((-t523 * t582 + t555 * t684) * qJD(1) + t607) * MDP(13) + MDP(7) * t570 + MDP(6) * t649 + (t582 * t600 - t564 - t699) * MDP(9) + (t700 + (-pkin(6) * qJDD(1) + t600) * t587) * MDP(10) + (-t555 * t693 + t697) * MDP(11) + ((t443 + t695) * t586 + (-t444 + t694) * t581) * MDP(12) + ((t521 * t582 - t555 * t689) * qJD(1) + t606) * MDP(14); (t423 * t613 + (-t562 * t498 - t364 + (-t382 + (-qJD(4) - qJD(5)) * t702) * t539) * t579 + (-t498 * t702 + (pkin(3) * t654 + qJD(5) * t562 - t383) * t539 - t624) * t584 + t716) * MDP(31) + (t629 * t547 + (-t455 * t523 + t513 * t585 + t547 * t655) * pkin(3) + t713) * MDP(23) + (-t680 * t547 + (-t580 * t513 + t523 * t611 + t547 * t654) * pkin(3) + t714) * MDP(24) + (-t444 - t694) * MDP(14) + (-t521 ^ 2 + t523 ^ 2) * MDP(12) + t518 * MDP(15) + t523 * t521 * MDP(11) + (-g(1) * t504 + g(2) * t502 - t466 * t555 - t523 * t540 + (t625 + t700) * t581 - t617) * MDP(16) + (t443 - t695) * MDP(13) + (t562 * t685 + (t382 * t584 - t383 * t579) * t539 + t423 * t407 + (-t580 * t692 - (-t579 * t585 - t691) * t539 * qJD(4)) * pkin(3) + (-(-pkin(3) * t691 - t562 * t579) * t539 - t369) * qJD(5) + t715) * MDP(30) + (g(1) * t505 - g(2) * t503 + g(3) * t688 - t465 * t555 + t521 * t540 - t603) * MDP(17) + t729; (-t392 * t547 + t713) * MDP(23) + (-t391 * t547 + t714) * MDP(24) + ((-t379 * t579 - t686) * t539 - t369 * qJD(5) + (-t407 * t611 + t539 * t653 + t685) * pkin(4) + t715) * MDP(30) + ((t380 * t539 - t364) * t579 + (-t379 * t539 - t624) * t584 + (t539 * t652 - t611 * t613 - t692) * pkin(4) + t716) * MDP(31) + t729; (t646 + t722) * MDP(27) + (-t631 + t721) * MDP(28) + (-t369 * t539 + t715) * MDP(30) + (-t579 * t364 - t584 * t365 - t368 * t539 + t716) * MDP(31) + (MDP(27) * t696 + MDP(28) * t613 - MDP(30) * t369 - MDP(31) * t698) * qJD(5) + t730;];
tau = t1;
