% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:01:31
% EndTime: 2019-03-09 06:01:41
% DurationCPUTime: 6.40s
% Computational Cost: add. (5360->502), mult. (11256->647), div. (0->0), fcn. (7611->14), ass. (0->229)
t552 = sin(pkin(10));
t530 = pkin(1) * t552 + pkin(7);
t509 = t530 * qJD(1);
t556 = sin(qJ(3));
t559 = cos(qJ(3));
t458 = qJD(2) * t559 - t556 * t509;
t598 = pkin(3) * t556 - pkin(8) * t559;
t495 = t598 * qJD(1);
t558 = cos(qJ(4));
t475 = t558 * t495;
t555 = sin(qJ(4));
t661 = t558 * t559;
t590 = pkin(4) * t556 - pkin(9) * t661;
t685 = pkin(8) + pkin(9);
t620 = qJD(4) * t685;
t710 = qJD(1) * t590 - t458 * t555 + t558 * t620 + t475;
t641 = qJD(1) * t559;
t619 = t555 * t641;
t649 = t558 * t458 + t555 * t495;
t709 = pkin(9) * t619 - t555 * t620 - t649;
t633 = qJD(4) * t556;
t611 = qJD(1) * t633;
t630 = qJD(1) * qJD(3);
t612 = t559 * t630;
t627 = qJDD(1) * t556;
t702 = qJD(3) * qJD(4) + t612 + t627;
t426 = (qJDD(3) - t611) * t555 + t702 * t558;
t636 = qJD(3) * t558;
t642 = qJD(1) * t556;
t486 = -t555 * t642 + t636;
t638 = qJD(3) * t555;
t487 = t558 * t642 + t638;
t554 = sin(qJ(5));
t600 = t702 * t555 + t558 * t611;
t578 = t558 * qJDD(3) - t600;
t684 = cos(qJ(5));
t613 = t684 * qJD(5);
t631 = qJD(5) * t554;
t376 = -t684 * t426 - t486 * t613 + t487 * t631 - t554 * t578;
t585 = t554 * t486 + t487 * t684;
t377 = qJD(5) * t585 + t554 * t426 - t684 * t578;
t432 = -t684 * t486 + t487 * t554;
t429 = t432 ^ 2;
t542 = t559 * qJDD(1);
t482 = t556 * t630 + qJDD(4) - t542;
t477 = qJDD(5) + t482;
t528 = -qJD(4) + t641;
t519 = -qJD(5) + t528;
t686 = t585 ^ 2;
t708 = t477 * MDP(23) + (-t519 * t585 - t377) * MDP(22) + t432 * MDP(19) * t585 + (-t432 * t519 - t376) * MDP(21) + (-t429 + t686) * MDP(20);
t707 = t432 * qJ(6);
t665 = t554 * t555;
t584 = t558 * t684 - t665;
t691 = qJD(4) + qJD(5);
t692 = t684 * qJD(4) + t613;
t651 = -t558 * t692 + t584 * t641 + t665 * t691;
t489 = t554 * t558 + t555 * t684;
t438 = t691 * t489;
t650 = -t489 * t641 + t438;
t504 = t530 * qJDD(1);
t706 = qJD(2) * qJD(3) + t504;
t635 = qJD(3) * t559;
t621 = t509 * t635 + t706 * t556;
t675 = qJDD(3) * pkin(3);
t589 = t621 - t675;
t626 = t559 * qJDD(2);
t414 = t589 - t626;
t548 = qJ(1) + pkin(10);
t539 = sin(t548);
t540 = cos(t548);
t597 = g(1) * t540 + g(2) * t539;
t586 = t597 * t556;
t573 = g(3) * t559 - t586;
t705 = -qJD(4) * pkin(8) * t528 + t414 + t573;
t617 = t555 * t635;
t632 = qJD(4) * t558;
t704 = t556 * t632 + t617;
t703 = qJD(3) * t458;
t448 = -qJD(3) * pkin(3) - t458;
t427 = -pkin(4) * t486 + t448;
t551 = qJ(4) + qJ(5);
t544 = sin(t551);
t545 = cos(t551);
t666 = t545 * t559;
t444 = -t539 * t666 + t540 * t544;
t446 = t539 * t544 + t540 * t666;
t459 = t556 * qJD(2) + t559 * t509;
t449 = qJD(3) * pkin(8) + t459;
t553 = cos(pkin(10));
t531 = -pkin(1) * t553 - pkin(2);
t478 = -pkin(3) * t559 - pkin(8) * t556 + t531;
t450 = t478 * qJD(1);
t406 = t449 * t558 + t450 * t555;
t413 = qJDD(3) * pkin(8) + qJDD(2) * t556 + t504 * t559 + t703;
t498 = t598 * qJD(3);
t428 = qJD(1) * t498 + qJDD(1) * t478;
t423 = t558 * t428;
t363 = pkin(4) * t482 - pkin(9) * t426 - qJD(4) * t406 - t413 * t555 + t423;
t622 = t558 * t413 + t555 * t428 + t450 * t632;
t634 = qJD(4) * t555;
t581 = -t449 * t634 + t622;
t368 = pkin(9) * t578 + t581;
t405 = -t449 * t555 + t558 * t450;
t398 = -pkin(9) * t487 + t405;
t391 = -pkin(4) * t528 + t398;
t399 = pkin(9) * t486 + t406;
t601 = -t554 * t363 - t684 * t368 - t391 * t613 + t399 * t631;
t667 = t545 * t556;
t701 = g(1) * t446 - g(2) * t444 + g(3) * t667 + t427 * t432 + t601;
t699 = qJ(6) * t585;
t392 = pkin(5) * t432 + qJD(6) + t427;
t698 = t392 * t585;
t660 = qJDD(2) - g(3);
t697 = t559 * t660;
t463 = t558 * t478;
t662 = t556 * t558;
t669 = t530 * t555;
t415 = -pkin(9) * t662 + t463 + (-pkin(4) - t669) * t559;
t493 = t530 * t661;
t645 = t555 * t478 + t493;
t664 = t555 * t556;
t425 = -pkin(9) * t664 + t645;
t652 = t554 * t415 + t684 * t425;
t695 = -t459 + (-t619 + t634) * pkin(4);
t512 = t685 * t555;
t513 = t685 * t558;
t646 = -t554 * t512 + t684 * t513;
t694 = qJD(5) * t646 + t709 * t554 + t710 * t684;
t693 = -t512 * t613 - t513 * t631 - t710 * t554 + t709 * t684;
t668 = t544 * t559;
t443 = t539 * t668 + t540 * t545;
t445 = t539 * t545 - t540 * t668;
t679 = g(3) * t556;
t690 = -g(1) * t445 + g(2) * t443 + t544 * t679;
t397 = t684 * t399;
t372 = t554 * t391 + t397;
t569 = -qJD(5) * t372 + t684 * t363 - t554 * t368;
t689 = -t427 * t585 + t569 + t690;
t616 = t555 * t633;
t576 = t558 * t635 - t616;
t688 = t482 * t662 - t528 * t576;
t546 = t558 * pkin(4);
t677 = pkin(3) + t546;
t499 = pkin(4) * t555 + pkin(5) * t544;
t676 = pkin(7) + t499;
t674 = t426 * t555;
t673 = t486 * t528;
t672 = t487 * t528;
t671 = t499 * t559;
t670 = t528 * t558;
t395 = t554 * t399;
t663 = t555 * t559;
t371 = t684 * t391 - t395;
t366 = t371 - t699;
t364 = -pkin(5) * t519 + t366;
t659 = -t366 + t364;
t658 = -t650 * qJ(6) + qJD(6) * t584 + t693;
t657 = -pkin(5) * t642 + t651 * qJ(6) - t489 * qJD(6) - t694;
t599 = t684 * t635;
t400 = t438 * t556 + t554 * t617 - t558 * t599;
t461 = t584 * t556;
t656 = -t461 * t377 + t400 * t432;
t401 = t555 * t599 - t554 * t616 - t631 * t664 + (t554 * t635 + t556 * t692) * t558;
t460 = t489 * t556;
t655 = t401 * t519 - t460 * t477;
t654 = t684 * t398 - t395;
t648 = t478 * t632 + t555 * t498;
t637 = qJD(3) * t556;
t647 = t558 * t498 + t637 * t669;
t466 = pkin(4) * t664 + t556 * t530;
t500 = pkin(5) * t545 + t546;
t549 = t556 ^ 2;
t644 = -t559 ^ 2 + t549;
t510 = qJD(1) * t531;
t639 = qJD(3) * t486;
t439 = t704 * pkin(4) + t530 * t635;
t618 = t528 * t638;
t610 = t376 * t559 + t585 * t637;
t609 = -t398 * t554 - t397;
t607 = t684 * t415 - t425 * t554;
t605 = -t426 * t559 + t487 * t637;
t604 = t528 * t530 + t449;
t603 = -t684 * t512 - t513 * t554;
t602 = -qJD(4) * t450 - t413;
t596 = g(1) * t539 - g(2) * t540;
t557 = sin(qJ(1));
t560 = cos(qJ(1));
t595 = g(1) * t557 - g(2) * t560;
t594 = -t376 * t460 + t401 * t585;
t593 = -t400 * t519 - t461 * t477;
t494 = pkin(3) + t500;
t547 = -qJ(6) - t685;
t592 = t494 * t559 - t547 * t556;
t588 = t595 * pkin(1);
t587 = pkin(2) + t592;
t583 = t377 * t559 - t432 * t637;
t582 = -t482 * t555 + t528 * t632;
t385 = t590 * qJD(3) + (-t493 + (pkin(9) * t556 - t478) * t555) * qJD(4) + t647;
t388 = (-t556 * t636 - t634 * t559) * t530 - t704 * pkin(9) + t648;
t580 = t554 * t385 + t684 * t388 + t415 * t613 - t425 * t631;
t577 = -qJD(1) * t510 + t597;
t575 = -pkin(8) * t482 - t448 * t528;
t574 = 0.2e1 * qJD(3) * t510 - qJDD(3) * t530;
t561 = qJD(3) ^ 2;
t571 = -0.2e1 * qJDD(1) * t531 - t530 * t561 + t596;
t568 = -qJD(5) * t652 + t684 * t385 - t554 * t388;
t566 = -pkin(4) * t578 + t589;
t564 = t377 * pkin(5) + qJDD(6) + t566;
t537 = pkin(4) * t684 + pkin(5);
t503 = qJDD(3) * t559 - t556 * t561;
t502 = qJDD(3) * t556 + t559 * t561;
t456 = t539 * t555 + t540 * t661;
t455 = t539 * t558 - t540 * t663;
t454 = -t539 * t661 + t540 * t555;
t453 = t539 * t663 + t540 * t558;
t418 = qJ(6) * t584 + t646;
t417 = -qJ(6) * t489 + t603;
t387 = t566 - t626;
t379 = -qJ(6) * t460 + t652;
t378 = -pkin(5) * t559 - qJ(6) * t461 + t607;
t370 = t654 - t699;
t369 = t609 + t707;
t367 = t372 - t707;
t360 = t564 - t626;
t359 = -qJ(6) * t401 - qJD(6) * t460 + t580;
t358 = pkin(5) * t637 + t400 * qJ(6) - t461 * qJD(6) + t568;
t357 = -qJ(6) * t377 - qJD(6) * t432 - t601;
t356 = t477 * pkin(5) + t376 * qJ(6) - qJD(6) * t585 + t569;
t1 = [(qJDD(1) * t549 + 0.2e1 * t556 * t612) * MDP(5) + (-t593 + t610) * MDP(21) + ((t552 ^ 2 + t553 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t588) * MDP(4) + (t556 * t574 + t559 * t571) * MDP(10) + (-t556 * t571 + t559 * t574) * MDP(11) + t595 * MDP(2) + (t605 + t688) * MDP(14) + (t648 * t528 - t645 * t482 - g(1) * t453 - g(2) * t455 + (-t604 * t634 + (t448 * t558 + t487 * t530) * qJD(3) + t622) * t559 + (-t448 * t634 + t414 * t558 + t530 * t426 + (-t530 * t670 - t406) * qJD(3)) * t556) * MDP(18) + ((-t578 + t618) * t559 + (t582 + t639) * t556) * MDP(15) + (t426 * t662 + t487 * t576) * MDP(12) + (g(1) * t560 + g(2) * t557) * MDP(3) + (-t356 * t461 - t357 * t460 - t358 * t585 - t359 * t432 + t364 * t400 - t367 * t401 + t376 * t378 - t377 * t379 + t556 * t596) * MDP(26) + (-g(1) * t443 - g(2) * t445 - t372 * t637 - t466 * t376 + t387 * t461 - t427 * t400 + t439 * t585 - t477 * t652 + t519 * t580 - t559 * t601) * MDP(25) + (-t376 * t461 - t400 * t585) * MDP(19) + (t583 + t655) * MDP(22) + (-t594 + t656) * MDP(20) + ((t486 * t558 - t487 * t555) * t635 + (t558 * t578 - t674 + (-t486 * t555 - t487 * t558) * qJD(4)) * t556) * MDP(13) + (-t482 * t559 - t528 * t637) * MDP(16) + (-t477 * t559 - t519 * t637) * MDP(23) + (-g(1) * t444 - g(2) * t446 + t371 * t637 + t466 * t377 + t387 * t460 + t427 * t401 + t439 * t432 + t477 * t607 - t519 * t568 - t559 * t569) * MDP(24) + t502 * MDP(7) + t503 * MDP(8) + (t357 * t379 + t367 * t359 + t356 * t378 + t364 * t358 + t360 * (pkin(5) * t460 + t466) + t392 * (pkin(5) * t401 + t439) + t588 + (-g(1) * t676 - g(2) * t587) * t540 + (g(1) * t587 - g(2) * t676) * t539) * MDP(27) + qJDD(1) * MDP(1) + 0.2e1 * (t542 * t556 - t630 * t644) * MDP(6) + (-(-t478 * t634 + t647) * t528 + t463 * t482 - g(1) * t454 - g(2) * t456 + (-t530 * t639 - t423 + t604 * t632 + (qJD(3) * t448 - t482 * t530 - t602) * t555) * t559 + (t405 * qJD(3) + t414 * t555 + t448 * t632 - t530 * t578) * t556) * MDP(17); t660 * MDP(4) + t503 * MDP(10) - t502 * MDP(11) + ((t578 + t618) * t559 + (t582 - t639) * t556) * MDP(17) + (t605 - t688) * MDP(18) + (-t583 + t655) * MDP(24) + (t593 + t610) * MDP(25) + (t594 + t656) * MDP(26) + (-t356 * t460 + t357 * t461 - t360 * t559 - t364 * t401 - t367 * t400 + t392 * t637 - g(3)) * MDP(27); MDP(7) * t627 + MDP(8) * t542 + qJDD(3) * MDP(9) + (qJD(3) * t459 + t556 * t577 - t621 + t697) * MDP(10) + (t703 + (qJD(3) * t509 - t660) * t556 + (t577 - t706) * t559) * MDP(11) + (-t487 * t670 + t674) * MDP(12) + ((t426 - t673) * t558 + (t578 + t672) * t555) * MDP(13) + ((-t487 * t556 + t528 * t661) * qJD(1) - t582) * MDP(14) + (t528 * t634 + t482 * t558 + (-t486 * t556 - t528 * t663) * qJD(1)) * MDP(15) + (-pkin(3) * t600 + t475 * t528 + t459 * t486 + (-t458 * t528 + t575) * t555 + (t675 - t705) * t558) * MDP(17) + (-pkin(3) * t426 - t459 * t487 - t649 * t528 + t705 * t555 + t575 * t558) * MDP(18) + (-t376 * t489 - t585 * t651) * MDP(19) + (-t376 * t584 - t377 * t489 + t432 * t651 - t585 * t650) * MDP(20) + (t477 * t489 + t519 * t651) * MDP(21) + (t477 * t584 + t519 * t650) * MDP(22) + (-g(3) * t666 - t377 * t677 - t387 * t584 + t427 * t650 + t432 * t695 + t477 * t603 + t519 * t694 + t597 * t667) * MDP(24) + (t376 * t677 + t387 * t489 - t651 * t427 - t646 * t477 + t519 * t693 + t573 * t544 + t585 * t695) * MDP(25) + (-t356 * t489 + t357 * t584 + t364 * t651 - t367 * t650 + t376 * t417 - t377 * t418 - t432 * t658 - t559 * t597 - t585 * t657 - t679) * MDP(26) + (t357 * t418 + t356 * t417 + t360 * (-pkin(5) * t584 - t677) - g(3) * t592 + (pkin(5) * t650 + t695) * t392 + t658 * t367 + t657 * t364 + t597 * (t494 * t556 + t547 * t559)) * MDP(27) + (t528 * MDP(16) - t405 * MDP(17) + t406 * MDP(18) - MDP(21) * t585 + t432 * MDP(22) + t519 * MDP(23) - t371 * MDP(24) + t372 * MDP(25)) * t642 + (-MDP(5) * t556 * t559 + MDP(6) * t644) * qJD(1) ^ 2; -t487 * t486 * MDP(12) + (-t486 ^ 2 + t487 ^ 2) * MDP(13) + (t426 + t673) * MDP(14) + (t578 - t672) * MDP(15) + t482 * MDP(16) + (-t449 * t632 - g(1) * t455 + g(2) * t453 - t406 * t528 - t448 * t487 + t423 + (t602 + t679) * t555) * MDP(17) + (g(1) * t456 - g(2) * t454 + g(3) * t662 - t405 * t528 - t448 * t486 - t581) * MDP(18) + (t609 * t519 + (-t487 * t432 + t477 * t684 + t519 * t631) * pkin(4) + t689) * MDP(24) + (-t654 * t519 + (-t554 * t477 - t487 * t585 + t519 * t613) * pkin(4) + t701) * MDP(25) + (-t364 * t432 + t367 * t585 + t369 * t585 + t370 * t432 + t537 * t376 + (-t377 * t554 + (-t432 * t684 + t554 * t585) * qJD(5)) * pkin(4)) * MDP(26) + (t356 * t537 - t367 * t370 - t364 * t369 - pkin(5) * t698 - g(1) * (t500 * t539 - t540 * t671) - g(2) * (-t500 * t540 - t539 * t671) + t499 * t679 + (t357 * t554 - t392 * t487 + (-t364 * t554 + t367 * t684) * qJD(5)) * pkin(4)) * MDP(27) + t708; (-t372 * t519 + t689) * MDP(24) + (-t371 * t519 + t701) * MDP(25) + (pkin(5) * t376 - t432 * t659) * MDP(26) + (t659 * t367 + (t356 + t690 - t698) * pkin(5)) * MDP(27) + t708; (-t429 - t686) * MDP(26) + (t364 * t585 + t367 * t432 + t564 - t586 - t697) * MDP(27);];
tau  = t1;
