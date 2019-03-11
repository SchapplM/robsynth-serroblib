% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% MDP [34x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(34,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [34 1]), ...
  'S6RPRRRR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [34x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:18:10
% EndTime: 2019-03-09 07:18:19
% DurationCPUTime: 6.98s
% Computational Cost: add. (7108->470), mult. (14370->603), div. (0->0), fcn. (10589->14), ass. (0->219)
t557 = sin(qJ(4));
t562 = cos(qJ(4));
t563 = cos(qJ(3));
t646 = qJD(1) * t563;
t558 = sin(qJ(3));
t647 = qJD(1) * t558;
t482 = t557 * t647 - t562 * t646;
t483 = -t557 * t646 - t562 * t647;
t556 = sin(qJ(5));
t561 = cos(qJ(5));
t441 = t482 * t556 + t561 * t483;
t560 = cos(qJ(6));
t638 = qJD(6) * t560;
t711 = -t441 * t560 + t638;
t594 = t561 * t482 - t483 * t556;
t490 = t557 * t563 + t558 * t562;
t548 = qJD(3) + qJD(4);
t452 = t548 * t490;
t633 = qJDD(1) * t563;
t634 = qJDD(1) * t558;
t424 = -t452 * qJD(1) - t557 * t634 + t562 * t633;
t643 = qJD(4) * t557;
t626 = t558 * t643;
t635 = qJD(1) * qJD(3);
t697 = t558 * t635 - t633;
t590 = -qJD(1) * t626 - t697 * t557;
t604 = t548 * t563;
t425 = (qJD(1) * t604 + t634) * t562 + t590;
t640 = qJD(5) * t561;
t641 = qJD(5) * t556;
t385 = t561 * t424 - t556 * t425 + t482 * t641 + t483 * t640;
t547 = qJDD(3) + qJDD(4);
t539 = qJDD(5) + t547;
t540 = qJD(5) + t548;
t555 = sin(qJ(6));
t628 = t560 * t385 + t555 * t539 + t540 * t638;
t639 = qJD(6) * t555;
t373 = t594 * t639 + t628;
t371 = t373 * t555;
t372 = t373 * t560;
t428 = t540 * t555 - t560 * t594;
t509 = t560 * t539;
t374 = t428 * qJD(6) + t385 * t555 - t509;
t579 = t594 * qJD(5) - t424 * t556 - t561 * t425;
t384 = qJDD(6) - t579;
t381 = t555 * t384;
t382 = t560 * t384;
t426 = -t560 * t540 - t555 * t594;
t695 = qJD(6) - t441;
t709 = t695 * t555;
t710 = t539 * MDP(25) - t441 ^ 2 * MDP(22) + (-t441 * t540 + t385) * MDP(23) + t579 * MDP(24) + (t711 * t695 + t381) * MDP(30) + (t441 * MDP(21) + MDP(22) * t594 - t540 * MDP(24) + t428 * MDP(30) + MDP(32) * t695) * t594 + (t711 * t428 + t371) * MDP(28) + (-t555 * t374 - t711 * t426 - t428 * t709 + t372) * MDP(29) + (-t426 * t594 - t695 * t709 + t382) * MDP(31);
t564 = cos(qJ(1));
t546 = g(2) * t564;
t559 = sin(qJ(1));
t677 = g(1) * t559;
t699 = t546 - t677;
t645 = qJD(3) * t558;
t453 = -t557 * t645 + t562 * t604 - t626;
t491 = -t557 * t558 + t562 * t563;
t593 = t561 * t490 + t491 * t556;
t401 = t593 * qJD(5) + t561 * t452 + t453 * t556;
t689 = -t490 * t556 + t561 * t491;
t708 = -t401 * t540 + t689 * t539;
t705 = pkin(5) * t594;
t476 = t482 * pkin(9);
t565 = -pkin(1) - pkin(7);
t512 = t565 * qJD(1) + qJD(2);
t471 = -pkin(8) * t647 + t512 * t558;
t458 = t557 * t471;
t472 = -pkin(8) * t646 + t563 * t512;
t461 = qJD(3) * pkin(3) + t472;
t616 = t562 * t461 - t458;
t415 = t476 + t616;
t413 = pkin(4) * t548 + t415;
t459 = t562 * t471;
t595 = -t461 * t557 - t459;
t678 = pkin(9) * t483;
t416 = -t595 + t678;
t670 = t416 * t556;
t391 = t413 * t561 - t670;
t389 = -pkin(5) * t540 - t391;
t673 = t389 * t441;
t510 = t565 * qJDD(1) + qJDD(2);
t494 = t563 * t510;
t434 = qJDD(3) * pkin(3) + t697 * pkin(8) - t512 * t645 + t494;
t644 = qJD(3) * t563;
t624 = t563 * t635;
t698 = t624 + t634;
t437 = -t698 * pkin(8) + t510 * t558 + t512 * t644;
t581 = t595 * qJD(4) + t562 * t434 - t557 * t437;
t378 = pkin(4) * t547 - pkin(9) * t424 + t581;
t685 = (qJD(4) * t461 + t437) * t562 + t557 * t434 - t471 * t643;
t379 = -pkin(9) * t425 + t685;
t669 = t416 * t561;
t392 = t413 * t556 + t669;
t684 = t392 * qJD(5) - t561 * t378 + t556 * t379;
t359 = -pkin(5) * t539 + t684;
t554 = qJ(3) + qJ(4);
t544 = qJ(5) + t554;
t530 = cos(t544);
t631 = t530 * t677;
t702 = t359 + t631;
t529 = sin(t544);
t701 = g(3) * t529 + t530 * t546;
t700 = qJDD(2) + t699;
t549 = qJDD(1) * qJ(2);
t550 = qJD(1) * qJD(2);
t600 = g(1) * t564 + g(2) * t559;
t583 = -t600 + 0.2e1 * t550;
t696 = 0.2e1 * t549 + t583;
t390 = pkin(10) * t540 + t392;
t501 = pkin(3) * t647 + qJD(1) * qJ(2);
t454 = -pkin(4) * t483 + t501;
t397 = -pkin(5) * t441 + pkin(10) * t594 + t454;
t367 = -t390 * t555 + t397 * t560;
t694 = t367 * t594 + t389 * t639 + t701 * t560;
t368 = t390 * t560 + t397 * t555;
t693 = -t368 * t594 + t389 * t638 + t702 * t555;
t522 = g(3) * t530;
t686 = (qJD(5) * t413 + t379) * t561 + t556 * t378 - t416 * t641;
t572 = -t454 * t441 - t699 * t529 + t522 - t686;
t570 = t454 * t594 - t631 - t684 + t701;
t690 = (t695 * pkin(10) - t705) * t695;
t676 = pkin(8) - t565;
t496 = t676 * t558;
t497 = t676 * t563;
t652 = -t562 * t496 - t557 * t497;
t675 = pkin(1) * qJDD(1);
t687 = t675 - t700;
t577 = t689 * qJD(5) - t452 * t556 + t561 * t453;
t683 = -t539 * t593 - t540 * t577;
t488 = t676 * t645;
t489 = qJD(3) * t497;
t642 = qJD(4) * t562;
t586 = t557 * t488 - t562 * t489 + t496 * t643 - t497 * t642;
t404 = -pkin(9) * t453 + t586;
t580 = -t652 * qJD(4) + t562 * t488 + t489 * t557;
t405 = pkin(9) * t452 + t580;
t614 = t496 * t557 - t562 * t497;
t432 = -pkin(9) * t491 + t614;
t433 = -pkin(9) * t490 + t652;
t597 = t432 * t561 - t433 * t556;
t365 = t597 * qJD(5) + t404 * t561 + t405 * t556;
t407 = t432 * t556 + t433 * t561;
t525 = t558 * pkin(3) + qJ(2);
t463 = pkin(4) * t490 + t525;
t408 = pkin(5) * t593 - pkin(10) * t689 + t463;
t358 = pkin(10) * t539 + t686;
t609 = qJD(6) * t397 + t358;
t682 = t359 * t689 - t407 * t384 - t389 * t401 - (qJD(6) * t408 + t365) * t695 - t593 * t609;
t679 = pkin(4) * t482;
t567 = qJD(1) ^ 2;
t674 = qJ(2) * t567;
t672 = t389 * t689;
t671 = t408 * t384;
t665 = t555 * t559;
t664 = t555 * t564;
t663 = t556 * t557;
t662 = t557 * t561;
t661 = t559 * t560;
t659 = t560 * t564;
t615 = -t472 * t557 - t459;
t417 = t615 - t678;
t653 = t562 * t472 - t458;
t418 = t476 + t653;
t533 = pkin(3) * t562 + pkin(4);
t656 = t417 * t556 + t418 * t561 - t533 * t640 - (-t557 * t641 + (t561 * t562 - t663) * qJD(4)) * pkin(3);
t655 = t417 * t561 - t418 * t556 + t533 * t641 + (t557 * t640 + (t556 * t562 + t662) * qJD(4)) * pkin(3);
t654 = -t452 * t548 + t491 * t547;
t651 = pkin(3) * t662 + t556 * t533;
t553 = t563 ^ 2;
t650 = t558 ^ 2 - t553;
t566 = qJD(3) ^ 2;
t649 = -t566 - t567;
t513 = pkin(3) * t644 + qJD(2);
t632 = qJDD(3) * t558;
t462 = t698 * pkin(3) + t549 + t550;
t410 = pkin(4) * t425 + t462;
t361 = -pkin(5) * t579 - pkin(10) * t385 + t410;
t608 = qJD(6) * t390 - t361;
t403 = -pkin(10) * t441 - t679 - t705;
t474 = pkin(10) + t651;
t537 = pkin(3) * t646;
t606 = qJD(6) * t474 + t403 + t537;
t531 = pkin(4) * t556 + pkin(10);
t605 = qJD(6) * t531 + t403;
t393 = t415 * t556 + t669;
t602 = pkin(4) * t641 - t393;
t394 = t415 * t561 - t670;
t601 = -pkin(4) * t640 + t394;
t436 = pkin(4) * t453 + t513;
t599 = -t384 * t474 - t673;
t598 = -t384 * t531 - t673;
t596 = -t453 * t548 - t490 * t547;
t591 = -pkin(3) * t663 + t533 * t561;
t588 = -t401 * t560 - t639 * t689;
t587 = 0.2e1 * qJ(2) * t635 + qJDD(3) * t565;
t585 = -t674 + t699;
t584 = -pkin(10) * t384 + t391 * t695 - t673;
t578 = -t560 * t702 + t694;
t576 = -t565 * t566 + t696;
t575 = -t555 * t701 + t693;
t542 = sin(t554);
t543 = cos(t554);
t571 = g(3) * t543 - t501 * t483 - t699 * t542 - t685;
t569 = g(3) * t542 + t501 * t482 + t699 * t543 + t581;
t568 = t482 * t483 * MDP(14) + (-t483 * t548 + t424) * MDP(16) + (-t482 * t548 + (-t548 * t646 - t634) * t562 - t590) * MDP(17) + (t482 ^ 2 - t483 ^ 2) * MDP(15) + t547 * MDP(18) + t710;
t541 = qJDD(3) * t563;
t532 = -pkin(4) * t561 - pkin(5);
t473 = -pkin(5) - t591;
t470 = t529 * t659 - t665;
t469 = t529 * t664 + t661;
t468 = t529 * t661 + t664;
t467 = -t529 * t665 + t659;
t456 = t537 - t679;
t369 = pkin(5) * t577 + pkin(10) * t401 + t436;
t366 = t407 * qJD(5) + t404 * t556 - t405 * t561;
t360 = t560 * t361;
t1 = [t683 * MDP(24) + (-t424 * t490 - t425 * t491 - t452 * t483 + t453 * t482) * MDP(15) + (t424 * t491 + t452 * t482) * MDP(14) + (t525 * t424 - t501 * t452 + t462 * t491 - t513 * t482 - t600 * t543 - t652 * t547 - t586 * t548) * MDP(20) + t708 * MDP(23) + t596 * MDP(17) + (-t558 * t566 + t541) * MDP(9) + (-t366 * t540 + t410 * t593 - t436 * t441 + t454 * t577 - t463 * t579 - t600 * t529 + t539 * t597) * MDP(26) + (t372 * t689 + t588 * t428) * MDP(28) + (-t563 * t566 - t632) * MDP(10) + (t385 * t689 + t401 * t594) * MDP(21) + (-t365 * t540 + t385 * t463 - t401 * t454 - t407 * t539 + t410 * t689 - t436 * t594 - t600 * t530) * MDP(27) + (-(-t426 * t560 - t428 * t555) * t401 + (-t371 - t374 * t560 + (t426 * t555 - t428 * t560) * qJD(6)) * t689) * MDP(29) + (-t385 * t593 - t401 * t441 + t577 * t594 + t579 * t689) * MDP(22) + t696 * MDP(5) + (-0.2e1 * t675 + t700) * MDP(4) + (t576 * t558 + t587 * t563) * MDP(12) + (-t587 * t558 + t576 * t563) * MDP(13) + (-t689 * t381 - t374 * t593 - t577 * t426 + (t401 * t555 - t638 * t689) * t695) * MDP(31) + (g(1) * t469 - g(2) * t467 + t366 * t428 - t368 * t577 - t597 * t373 + (-(-qJD(6) * t407 + t369) * t695 - t671 + t608 * t593 - qJD(6) * t672) * t555 + t682 * t560) * MDP(34) + (-g(1) * t470 - g(2) * t468 + t360 * t593 + t366 * t426 + t367 * t577 - t597 * t374 + (t369 * t695 + t671 + (-t390 * t593 - t407 * t695 + t672) * qJD(6)) * t560 + t682 * t555) * MDP(33) + (t384 * t593 + t577 * t695) * MDP(32) + (t373 * t593 + t382 * t689 + t428 * t577 + t588 * t695) * MDP(30) - t699 * MDP(2) + (t525 * t425 + t501 * t453 + t462 * t490 - t513 * t483 - t600 * t542 + t614 * t547 + t580 * t548) * MDP(19) + 0.2e1 * (-t558 * t633 + t650 * t635) * MDP(8) + t654 * MDP(16) + t600 * MDP(3) + (qJDD(1) * t553 - 0.2e1 * t558 * t624) * MDP(7) + qJDD(1) * MDP(1) + (t687 * pkin(1) + (t583 + t549) * qJ(2)) * MDP(6); qJDD(1) * MDP(4) - t567 * MDP(5) + (-t674 - t687) * MDP(6) + (t649 * t558 + t541) * MDP(12) + (t649 * t563 - t632) * MDP(13) + (qJD(1) * t483 + t654) * MDP(19) + (qJD(1) * t482 + t596) * MDP(20) + (qJD(1) * t441 + t708) * MDP(26) + (qJD(1) * t594 + t683) * MDP(27) + (-t374 * t689 - t381 * t593 + t401 * t426) * MDP(33) + (-t373 * t689 - t382 * t593 + t401 * t428) * MDP(34) + ((-qJD(1) * t560 - t555 * t577 - t593 * t638) * MDP(33) + (qJD(1) * t555 - t560 * t577 + t593 * t639) * MDP(34)) * t695; (t653 * t548 + (t482 * t646 - t557 * t547 - t548 * t642) * pkin(3) + t571) * MDP(20) + (-t615 * t548 + (t483 * t646 + t547 * t562 - t548 * t643) * pkin(3) + t569) * MDP(19) + (g(3) * t563 + (-t510 - t585) * t558) * MDP(13) + (g(3) * t558 + t585 * t563 + t494) * MDP(12) + t568 + MDP(9) * t633 + (t473 * t374 + t599 * t555 + t655 * t426 + (t656 * t555 - t606 * t560) * t695 + t578) * MDP(33) + (t473 * t373 + t599 * t560 + t655 * t428 + (t606 * t555 + t656 * t560) * t695 + t575) * MDP(34) - MDP(10) * t634 + qJDD(3) * MDP(11) + (t441 * t456 + t591 * t539 - t655 * t540 + t570) * MDP(26) + (t456 * t594 - t651 * t539 + t656 * t540 + t572) * MDP(27) + (t563 * t558 * MDP(7) - t650 * MDP(8)) * t567; (t394 * t540 + (-t482 * t594 - t539 * t556 - t540 * t640) * pkin(4) + t572) * MDP(27) + (t393 * t540 + (-t441 * t482 + t539 * t561 - t540 * t641) * pkin(4) + t570) * MDP(26) + (t532 * t374 + t598 * t555 + t602 * t426 + (t601 * t555 - t605 * t560) * t695 + t578) * MDP(33) + t568 + (-t595 * t548 + t569) * MDP(19) + (t532 * t373 + t598 * t560 + t602 * t428 + (t605 * t555 + t601 * t560) * t695 + t575) * MDP(34) + (t616 * t548 + t571) * MDP(20); (t392 * t540 + t570) * MDP(26) + (t391 * t540 + t572) * MDP(27) + (-pkin(5) * t374 - t392 * t426 + t584 * t555 + (-t702 - t690) * t560 + t694) * MDP(33) + (-pkin(5) * t373 - t392 * t428 + t584 * t560 + (-t701 + t690) * t555 + t693) * MDP(34) + t710; t428 * t426 * MDP(28) + (-t426 ^ 2 + t428 ^ 2) * MDP(29) + (t426 * t695 + t628) * MDP(30) + (t428 * t695 + t509) * MDP(31) + t384 * MDP(32) + (-g(1) * t467 - g(2) * t469 + t368 * t695 - t389 * t428 + t360) * MDP(33) + (g(1) * t468 - g(2) * t470 + t367 * t695 + t389 * t426) * MDP(34) + ((-t358 + t522) * MDP(34) + (MDP(31) * t594 - MDP(33) * t390 - MDP(34) * t397) * qJD(6)) * t560 + (qJD(6) * t594 * MDP(30) + (-qJD(6) * t540 - t385) * MDP(31) + (-t609 + t522) * MDP(33) + t608 * MDP(34)) * t555;];
tau  = t1;
