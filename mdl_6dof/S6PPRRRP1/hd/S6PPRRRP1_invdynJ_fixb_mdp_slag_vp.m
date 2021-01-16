% Calculate vector of inverse dynamics joint torques for
% S6PPRRRP1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:51
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PPRRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:50:32
% EndTime: 2021-01-16 00:50:47
% DurationCPUTime: 9.01s
% Computational Cost: add. (4579->524), mult. (11411->711), div. (0->0), fcn. (10404->14), ass. (0->238)
t542 = cos(pkin(6));
t516 = qJDD(1) * t542 + qJDD(2);
t519 = qJD(1) * t542 + qJD(2);
t537 = sin(pkin(7));
t538 = sin(pkin(6));
t546 = sin(qJ(3));
t539 = cos(pkin(12));
t541 = cos(pkin(7));
t644 = t539 * t541;
t535 = sin(pkin(12));
t549 = cos(qJ(3));
t649 = t535 * t549;
t574 = t546 * t644 + t649;
t604 = t549 * t644;
t650 = t535 * t546;
t575 = t604 - t650;
t623 = qJD(3) * t549;
t404 = qJDD(3) * pkin(9) + (t516 * t546 + t519 * t623) * t537 + (qJD(1) * qJD(3) * t575 + qJDD(1) * t574) * t538;
t560 = t574 * t538;
t647 = t537 * t546;
t442 = qJD(1) * t560 + t519 * t647;
t440 = qJD(3) * pkin(9) + t442;
t610 = qJDD(1) * t538;
t591 = t539 * t610;
t469 = t516 * t541 - t537 * t591;
t628 = qJD(1) * t538;
t600 = t539 * t628;
t471 = t519 * t541 - t537 * t600;
t545 = sin(qJ(4));
t548 = cos(qJ(4));
t618 = qJD(4) * t548;
t619 = qJD(4) * t545;
t581 = t404 * t545 + t440 * t618 - t548 * t469 + t471 * t619;
t367 = -qJDD(4) * pkin(4) + t581;
t544 = sin(qJ(5));
t547 = cos(qJ(5));
t616 = qJD(5) * t545;
t592 = qJD(3) * t616;
t577 = (-qJDD(4) + t592) * t547;
t609 = qJDD(3) * t545;
t624 = qJD(3) * t548;
t444 = ((qJD(5) + t624) * qJD(4) + t609) * t544 + t577;
t362 = t444 * pkin(5) + qJDD(6) + t367;
t643 = t542 * t541;
t645 = t538 * t537;
t478 = t539 * t643 - t645;
t536 = sin(pkin(11));
t540 = cos(pkin(11));
t652 = t535 * t541;
t457 = -t478 * t540 + t536 * t652;
t651 = t535 * t542;
t480 = t536 * t539 + t540 * t651;
t430 = t457 * t546 - t480 * t549;
t648 = t537 * t542;
t476 = t538 * t541 + t539 * t648;
t653 = t535 * t537;
t453 = t476 * t540 - t536 * t653;
t393 = t430 * t545 - t453 * t548;
t454 = t478 * t536 + t540 * t652;
t481 = -t536 * t651 + t539 * t540;
t424 = t454 * t546 - t481 * t549;
t452 = t476 * t536 + t540 * t653;
t395 = t424 * t545 + t452 * t548;
t479 = t538 * t644 + t648;
t458 = t479 * t546 + t538 * t649;
t477 = t539 * t645 - t643;
t661 = t477 * t548;
t703 = t458 * t545 + t661;
t704 = g(1) * t395 + g(2) * t393 - g(3) * t703;
t714 = -t362 - t704;
t550 = qJD(4) ^ 2;
t605 = t538 * t650;
t461 = -t479 * t549 + t605;
t425 = t454 * t549 + t481 * t546;
t431 = t457 * t549 + t480 * t546;
t699 = g(1) * t425 + g(2) * t431;
t563 = g(3) * t461 + t699;
t584 = t541 * t600;
t625 = qJD(3) * t546;
t599 = t537 * t625;
t601 = t535 * t628;
t553 = -(t516 * t537 + t541 * t591) * t549 + qJDD(1) * t605 + t519 * t599 + t584 * t625 + t601 * t623;
t705 = qJD(3) * t442 - t553;
t713 = -0.2e1 * qJDD(3) * pkin(3) + pkin(9) * t550 - t563 - t705;
t525 = pkin(5) * t547 + pkin(4);
t685 = qJ(6) + pkin(10);
t712 = -t525 * t548 - t685 * t545;
t392 = t424 * t548 - t452 * t545;
t397 = t430 * t548 + t453 * t545;
t583 = pkin(4) * t545 - pkin(10) * t548;
t506 = t583 * qJD(4);
t708 = -t442 + t506;
t612 = qJD(3) * qJD(4);
t593 = t548 * t612;
t569 = -t593 - t609;
t706 = qJD(4) * qJD(5) - t569;
t662 = t477 * t545;
t433 = t458 * t548 - t662;
t441 = -t546 * t601 + (t519 * t537 + t584) * t549;
t510 = -pkin(4) * t548 - pkin(10) * t545 - pkin(3);
t615 = qJD(5) * t547;
t639 = t547 * t548;
t702 = -t441 * t639 + t510 * t615 + t544 * t708;
t701 = -t440 * t545 + t471 * t548;
t451 = t542 * t647 + t560;
t423 = t451 * t548 - t662;
t646 = t537 * t549;
t603 = t542 * t646;
t450 = -t538 * t604 - t603 + t605;
t447 = t450 * t547;
t390 = -t423 * t544 + t447;
t641 = t544 * t548;
t693 = pkin(9) * t544;
t700 = t441 * t641 + t547 * t708 + t619 * t693;
t697 = -g(2) * (t397 * t544 + t431 * t547) - g(1) * (t392 * t544 + t425 * t547) - g(3) * (-t433 * t544 + t447);
t620 = qJD(4) * t544;
t626 = qJD(3) * t545;
t497 = t547 * t626 + t620;
t696 = t497 ^ 2;
t608 = t548 * qJDD(3);
t568 = -t545 * t612 + t608;
t492 = qJDD(5) - t568;
t695 = pkin(5) * t492;
t613 = t547 * qJD(4);
t495 = t544 * t626 - t613;
t694 = pkin(5) * t495;
t684 = qJD(3) * pkin(3);
t585 = t544 * qJDD(4) + t547 * t706;
t443 = t544 * t592 - t585;
t683 = qJ(6) * t443;
t682 = qJ(6) * t444;
t681 = qJ(6) * t545;
t679 = t362 * t544;
t678 = t367 * t544;
t418 = t440 * t548 + t471 * t545;
t413 = qJD(4) * pkin(10) + t418;
t435 = qJD(3) * t510 - t441;
t376 = t413 * t547 + t435 * t544;
t374 = -qJ(6) * t495 + t376;
t520 = -qJD(5) + t624;
t677 = t374 * t520;
t675 = t424 * t544;
t674 = t430 * t544;
t672 = t441 * t495;
t671 = t441 * t497;
t670 = t443 * t544;
t669 = t450 * t544;
t668 = t451 * t544;
t664 = t469 * t545;
t658 = t495 * t520;
t657 = t497 * t520;
t656 = t497 * t544;
t655 = t497 * t547;
t640 = t545 * t547;
t375 = -t413 * t544 + t435 * t547;
t373 = -qJ(6) * t497 + t375;
t368 = -pkin(5) * t520 + t373;
t638 = -t373 + t368;
t521 = pkin(9) * t639;
t579 = pkin(5) * t545 - qJ(6) * t639;
t614 = qJD(6) * t547;
t637 = -t545 * t614 + t579 * qJD(4) + (-t521 + (-t510 + t681) * t544) * qJD(5) + t700;
t505 = t583 * qJD(3);
t636 = -t505 * t544 - t547 * t701;
t635 = (-pkin(9) * qJD(4) - qJ(6) * qJD(5)) * t640 + (-qJD(6) * t545 + (-pkin(9) * qJD(5) - qJ(6) * qJD(4)) * t548) * t544 + t702;
t590 = qJD(5) * t685;
t597 = t544 * t624;
t634 = qJ(6) * t597 - t544 * t590 + t614 + t636;
t486 = t547 * t505;
t633 = -t579 * qJD(3) - t547 * t590 - t486 + (-qJD(6) + t701) * t544;
t630 = t510 * t544 + t521;
t533 = t545 ^ 2;
t629 = -t548 ^ 2 + t533;
t622 = qJD(4) * t495;
t621 = qJD(4) * t497;
t617 = qJD(5) * t544;
t607 = MDP(18) + MDP(20);
t606 = MDP(19) + MDP(21);
t602 = t704 * t544;
t598 = t537 * t623;
t596 = t520 * t613;
t595 = t520 * t617;
t594 = t520 * t615;
t589 = -qJD(6) - t694;
t366 = qJDD(4) * pkin(10) + qJD(4) * t701 + t404 * t548 + t664;
t383 = qJD(3) * t506 + qJDD(3) * t510 + t553;
t587 = t366 * t547 + t383 * t544 - t413 * t617 + t435 * t615;
t391 = t423 * t547 + t669;
t422 = t451 * t545 + t661;
t502 = pkin(3) * t644 + pkin(9) * t535;
t580 = pkin(3) * t645 - t502 * t542;
t501 = -pkin(3) * t535 + pkin(9) * t644;
t578 = pkin(9) * t645 - t501 * t542;
t412 = -qJD(4) * pkin(4) - t701;
t483 = t541 * t545 + t548 * t647;
t464 = -t483 * t544 - t547 * t646;
t576 = -t483 * t547 + t544 * t646;
t482 = -t541 * t548 + t545 * t647;
t573 = t492 * t544 - t594;
t572 = t492 * t547 + t595;
t571 = -MDP(11) * t548 + MDP(12) * t545 - MDP(4);
t570 = -g(3) * t542 + (-g(1) * t536 + g(2) * t540) * t538;
t567 = -g(1) * (-t424 * t547 + t425 * t641) - g(2) * (-t430 * t547 + t431 * t641) - g(3) * (t451 * t547 + t461 * t641);
t566 = -g(1) * (-t425 * t639 - t675) - g(2) * (-t431 * t639 - t674) - g(3) * (-t461 * t639 + t668);
t565 = -g(1) * t392 - g(2) * t397 + g(3) * t433;
t559 = -pkin(10) * t492 - t412 * t520;
t439 = -t441 - t684;
t557 = -pkin(9) * qJDD(4) + (t439 + t441 - t684) * qJD(4);
t555 = -g(1) * (t392 * t547 - t425 * t544) - g(2) * (t397 * t547 - t431 * t544) - g(3) * (-t433 * t547 - t669) - t587;
t382 = t547 * t383;
t554 = -qJD(5) * t376 - t366 * t544 + t382;
t552 = t554 + t697;
t551 = qJD(3) ^ 2;
t512 = t685 * t547;
t511 = t685 * t544;
t507 = (pkin(5) * t544 + pkin(9)) * t545;
t500 = pkin(3) * t652 - pkin(9) * t539;
t499 = pkin(3) * t539 + pkin(9) * t652;
t494 = t547 * t510;
t491 = t495 ^ 2;
t470 = pkin(9) * t618 + (t544 * t618 + t545 * t615) * pkin(5);
t466 = -t544 * t681 + t630;
t463 = qJD(4) * t483 + t545 * t598;
t462 = -qJD(4) * t482 + t548 * t598;
t449 = -qJ(6) * t640 + t494 + (-pkin(5) - t693) * t548;
t446 = t451 * qJD(3);
t445 = (t538 * t575 + t603) * qJD(3);
t411 = qJD(5) * t576 - t462 * t544 + t547 * t599;
t410 = qJD(5) * t464 + t462 * t547 + t544 * t599;
t406 = pkin(5) * t597 + t418;
t400 = t412 - t589;
t389 = qJD(4) * t423 + t445 * t545;
t388 = -qJD(4) * t422 + t445 * t548;
t361 = -qJD(5) * t391 - t388 * t544 + t446 * t547;
t360 = qJD(5) * t390 + t388 * t547 + t446 * t544;
t359 = -qJD(6) * t495 + t587 - t682;
t358 = -qJD(6) * t497 + t554 + t683 + t695;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t516 * t542 - g(3) + (t535 ^ 2 + t539 ^ 2) * t538 ^ 2 * qJDD(1)) * MDP(2) + (-qJD(4) * t389 - qJDD(4) * t422) * MDP(11) + (-qJD(4) * t388 - qJDD(4) * t423) * MDP(12) + (-t360 * t495 - t361 * t497 + t390 * t443 - t391 * t444) * MDP(22) + (t358 * t390 + t359 * t391 + t360 * t374 + t361 * t368 + t362 * t422 + t389 * t400 - g(3)) * MDP(23) + t607 * (-t361 * t520 + t389 * t495 + t390 * t492 + t422 * t444) + t606 * (t360 * t520 + t389 * t497 - t391 * t492 - t422 * t443) + (-MDP(5) * t451 + t450 * t571) * qJDD(3) + (-t446 * MDP(4) - t445 * MDP(5) + (-t446 * t548 + t450 * t619) * MDP(11) + (t446 * t545 + t450 * t618) * MDP(12)) * qJD(3); (t570 + t516) * MDP(2) + (-qJD(4) * t463 - qJDD(4) * t482) * MDP(11) + (-qJD(4) * t462 - qJDD(4) * t483) * MDP(12) + (-t410 * t495 - t411 * t497 + t443 * t464 + t444 * t576) * MDP(22) + (t358 * t464 - t359 * t576 + t362 * t482 + t368 * t411 + t374 * t410 + t400 * t463 + t570) * MDP(23) + t607 * (-t411 * t520 + t444 * t482 + t463 * t495 + t464 * t492) + t606 * (t410 * t520 - t443 * t482 + t463 * t497 + t492 * t576) + ((-MDP(5) * qJDD(3) + t551 * t571) * t546 + (MDP(11) * t568 + MDP(12) * t569 + MDP(4) * qJDD(3) - MDP(5) * t551) * t549) * t537; qJDD(3) * MDP(3) + (g(3) * t450 + t699 + t705) * MDP(4) + (-t516 * t647 - g(1) * t424 - g(2) * t430 + g(3) * t451 - t574 * t610 + (-t519 * t646 - t575 * t628 + t441) * qJD(3)) * MDP(5) + (qJDD(3) * t533 + 0.2e1 * t545 * t593) * MDP(6) + 0.2e1 * (t545 * t608 - t612 * t629) * MDP(7) + (qJDD(4) * t545 + t548 * t550) * MDP(8) + (qJDD(4) * t548 - t545 * t550) * MDP(9) + (t557 * t545 - t548 * t713) * MDP(11) + (t545 * t713 + t557 * t548) * MDP(12) + (-t443 * t640 + (-t544 * t616 + t548 * t613) * t497) * MDP(13) + ((-t495 * t547 - t656) * t618 + (t670 - t444 * t547 + (t495 * t544 - t655) * qJD(5)) * t545) * MDP(14) + ((t443 - t596) * t548 + (t572 + t621) * t545) * MDP(15) + ((t520 * t620 + t444) * t548 + (-t573 - t622) * t545) * MDP(16) + (-t492 * t548 - t520 * t619) * MDP(17) + (t494 * t492 + (t510 * t617 - t700) * t520 + (t413 * t615 - t382 + (t594 + t622) * pkin(9) + (-pkin(9) * t492 + qJD(4) * t412 + qJD(5) * t435 + t366) * t544) * t548 + (pkin(9) * t444 + qJD(4) * t375 + t412 * t615 - t672 + t678) * t545 + t566) * MDP(18) + (-t630 * t492 + t702 * t520 + (t412 * t613 + (-t595 + t621) * pkin(9) + t587) * t548 + (-t412 * t617 - t376 * qJD(4) + t367 * t547 - t671 + (-t443 - t596) * pkin(9)) * t545 + t567) * MDP(19) + (t444 * t507 + t449 * t492 + t470 * t495 + (t400 * t620 - t358) * t548 - t637 * t520 + (qJD(4) * t368 + t400 * t615 - t672 + t679) * t545 + t566) * MDP(20) + (-t443 * t507 - t466 * t492 + t470 * t497 + (t400 * t613 + t359) * t548 + t635 * t520 + (-qJD(4) * t374 + t362 * t547 - t400 * t617 - t671) * t545 + t567) * MDP(21) + (t443 * t449 - t444 * t466 - t637 * t497 - t635 * t495 + (-t368 * t547 - t374 * t544) * t618 + (-t358 * t547 - t359 * t544 + (t368 * t544 - t374 * t547) * qJD(5) + t563) * t545) * MDP(22) + (t359 * t466 + t358 * t449 + t362 * t507 - g(1) * (-pkin(5) * t675 - (t540 * t499 - t536 * t578) * t546 + (-t540 * t500 + t536 * t580) * t549 + t712 * t425) - g(2) * (-pkin(5) * t674 - (t536 * t499 + t540 * t578) * t546 + (-t536 * t500 - t540 * t580) * t549 + t712 * t431) - g(3) * (pkin(5) * t668 - (-pkin(9) * t648 - t501 * t538) * t546 + (pkin(3) * t648 + t502 * t538) * t549 + t712 * t461) + (-t441 * t545 + t470) * t400 + t635 * t374 + t637 * t368) * MDP(23); MDP(8) * t609 + MDP(9) * t608 + qJDD(4) * MDP(10) + (qJD(4) * t418 - t439 * t626 - t581 - t704) * MDP(11) + (-t664 + (-qJD(3) * t439 - t404) * t548 + t565) * MDP(12) + (-t520 * t655 - t670) * MDP(13) + ((-t443 + t658) * t547 + (-t444 + t657) * t544) * MDP(14) + ((-t497 * t545 + t520 * t639) * qJD(3) + t573) * MDP(15) + ((t495 * t545 - t520 * t641) * qJD(3) + t572) * MDP(16) + t520 * MDP(17) * t626 + (-t375 * t626 - pkin(4) * t444 - t418 * t495 + t486 * t520 + (-t520 * t701 + t559) * t544 + (pkin(10) * qJD(5) * t520 - t367 - t704) * t547) * MDP(18) + (t376 * t626 + pkin(4) * t443 + t678 - t418 * t497 + (-pkin(10) * t617 + t636) * t520 + t559 * t547 + t602) * MDP(19) + (-t368 * t626 - t406 * t495 - t444 * t525 - t492 * t511 - t633 * t520 + (-t400 * t624 + (t400 + t694) * qJD(5)) * t544 + t714 * t547) * MDP(20) + (t679 - t406 * t497 + t443 * t525 - t492 * t512 + t634 * t520 + (pkin(5) * t656 + t400 * t547) * qJD(5) + (t374 * t545 - t400 * t639) * qJD(3) + t602) * MDP(21) + (-t443 * t511 - t444 * t512 - t633 * t497 - t634 * t495 + (t368 * t520 + t359) * t547 + (-t358 + t677) * t544 - t565) * MDP(22) + (t359 * t512 - t358 * t511 - t362 * t525 - g(1) * (-t392 * t685 + t395 * t525) - g(2) * (t393 * t525 - t397 * t685) - g(3) * (t433 * t685 - t525 * t703) + (pkin(5) * t617 - t406) * t400 + t634 * t374 + t633 * t368) * MDP(23) + (-MDP(6) * t545 * t548 + MDP(7) * t629) * t551; t497 * t495 * MDP(13) + (-t491 + t696) * MDP(14) + (-t443 - t658) * MDP(15) + (-t444 - t657) * MDP(16) + t492 * MDP(17) + (-t376 * t520 - t412 * t497 + t552) * MDP(18) + (-t375 * t520 + t412 * t495 + t555) * MDP(19) + (0.2e1 * t695 + t683 - t677 + (-t400 + t589) * t497 + t552) * MDP(20) + (-pkin(5) * t696 + t682 - t373 * t520 + (qJD(6) + t400) * t495 + t555) * MDP(21) + (pkin(5) * t443 - t495 * t638) * MDP(22) + (t638 * t374 + (-t400 * t497 + t358 + t697) * pkin(5)) * MDP(23); (t577 - t657) * MDP(20) + (t585 + t658) * MDP(21) + (-t491 - t696) * MDP(22) + (t368 * t497 + t374 * t495 - t714) * MDP(23) + (MDP(20) * t706 - MDP(21) * t592) * t544;];
tau = t1;
