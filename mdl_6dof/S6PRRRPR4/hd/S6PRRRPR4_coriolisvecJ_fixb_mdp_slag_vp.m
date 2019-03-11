% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:20:41
% EndTime: 2019-03-08 23:20:54
% DurationCPUTime: 7.67s
% Computational Cost: add. (5016->476), mult. (12707->673), div. (0->0), fcn. (9658->12), ass. (0->224)
t557 = cos(qJ(6));
t554 = sin(qJ(4));
t558 = cos(qJ(4));
t614 = t558 * qJD(3);
t555 = sin(qJ(3));
t624 = qJD(2) * t555;
t518 = t554 * t624 - t614;
t622 = qJD(3) * t554;
t520 = t558 * t624 + t622;
t549 = sin(pkin(12));
t551 = cos(pkin(12));
t575 = -t518 * t551 - t520 * t549;
t647 = t557 * t575;
t455 = t518 * t549 - t520 * t551;
t553 = sin(qJ(6));
t663 = t455 * t553;
t398 = t647 + t663;
t559 = cos(qJ(3));
t623 = qJD(2) * t559;
t539 = -qJD(4) + t623;
t532 = -qJD(6) + t539;
t666 = t398 * t532;
t613 = qJD(2) * qJD(3);
t594 = t559 * t613;
t618 = qJD(4) * t555;
t597 = t554 * t618;
t612 = qJD(3) * qJD(4);
t477 = -qJD(2) * t597 + (t594 + t612) * t558;
t556 = sin(qJ(2));
t550 = sin(pkin(6));
t627 = qJD(1) * t550;
t606 = t556 * t627;
t524 = qJD(2) * pkin(8) + t606;
t552 = cos(pkin(6));
t653 = t552 * t555;
t537 = qJD(1) * t653;
t484 = t524 * t559 + t537;
t474 = qJD(3) * pkin(9) + t484;
t527 = -pkin(3) * t559 - pkin(9) * t555 - pkin(2);
t560 = cos(qJ(2));
t605 = t560 * t627;
t485 = qJD(2) * t527 - t605;
t660 = t485 * t554;
t428 = t474 * t558 + t660;
t625 = qJD(2) * t550;
t602 = t560 * t625;
t621 = qJD(3) * t555;
t626 = qJD(1) * t559;
t447 = -t524 * t621 + (qJD(3) * t552 + t602) * t626;
t584 = pkin(3) * t555 - pkin(9) * t559;
t522 = t584 * qJD(3);
t482 = (t522 + t606) * qJD(2);
t587 = t447 * t554 - t482 * t558;
t563 = -qJD(4) * t428 - t587;
t595 = t555 * t613;
t365 = pkin(4) * t595 - qJ(5) * t477 - qJD(5) * t520 + t563;
t617 = qJD(4) * t558;
t596 = t555 * t617;
t620 = qJD(3) * t559;
t600 = t554 * t620;
t678 = t596 + t600;
t478 = qJD(2) * t678 + t554 * t612;
t619 = qJD(4) * t554;
t566 = t447 * t558 - t474 * t619 + t482 * t554 + t485 * t617;
t367 = -qJ(5) * t478 - qJD(5) * t518 + t566;
t353 = t365 * t549 + t367 * t551;
t422 = -t477 * t549 - t478 * t551;
t349 = pkin(10) * t422 + t353;
t483 = -t524 * t555 + t552 * t626;
t473 = -qJD(3) * pkin(3) - t483;
t442 = pkin(4) * t518 + qJD(5) + t473;
t391 = -pkin(5) * t575 + t442;
t352 = t365 * t551 - t367 * t549;
t423 = t477 * t551 - t478 * t549;
t348 = pkin(5) * t595 - pkin(10) * t423 + t352;
t426 = -t474 * t554 + t485 * t558;
t407 = -qJ(5) * t520 + t426;
t389 = -pkin(4) * t539 + t407;
t408 = -qJ(5) * t518 + t428;
t654 = t551 * t408;
t369 = t389 * t549 + t654;
t676 = pkin(10) * t575;
t360 = t369 + t676;
t615 = qJD(6) * t553;
t589 = t553 * t348 - t360 * t615;
t691 = -t557 * t349 - t391 * t398 - t589;
t677 = -t455 * t557 + t553 * t575;
t690 = MDP(25) * t595 + (-t398 ^ 2 + t677 ^ 2) * MDP(22) - t398 * MDP(21) * t677;
t645 = t559 * t560;
t689 = -(t554 * t556 + t558 * t645) * t627 + t554 * t522 + t527 * t617;
t670 = pkin(8) * t554;
t688 = t558 * t522 + t621 * t670 - (-t554 * t645 + t556 * t558) * t627;
t667 = t677 * t532;
t646 = t558 * t559;
t541 = pkin(8) * t646;
t572 = pkin(4) * t555 - qJ(5) * t646;
t616 = qJD(5) * t558;
t686 = -t555 * t616 + t572 * qJD(3) + (-t541 + (qJ(5) * t555 - t527) * t554) * qJD(4) + t688;
t650 = t555 * t558;
t685 = (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t650 + (-qJD(5) * t555 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t559) * t554 + t689;
t521 = t584 * qJD(2);
t586 = -t483 * t554 + t521 * t558;
t669 = -qJ(5) - pkin(9);
t591 = qJD(4) * t669;
t684 = -qJD(2) * t572 - qJD(5) * t554 + t558 * t591 - t586;
t601 = t554 * t623;
t635 = t483 * t558 + t521 * t554;
t683 = -qJ(5) * t601 - t554 * t591 - t616 + t635;
t590 = t348 * t557 - t553 * t349;
t682 = -t391 * t677 + t590;
t681 = pkin(10) * t455;
t512 = t549 * t558 + t551 * t554;
t567 = t512 * t559;
t680 = qJD(2) * t567 - qJD(4) * t512;
t574 = t549 * t554 - t551 * t558;
t679 = t539 * t574;
t675 = MDP(5) * t555;
t547 = t555 ^ 2;
t674 = MDP(6) * (-t559 ^ 2 + t547);
t641 = -t549 * t685 + t551 * t686;
t640 = t549 * t686 + t551 * t685;
t637 = t549 * t683 + t551 * t684;
t636 = t549 * t684 - t551 * t683;
t672 = -t484 + (-t601 + t619) * pkin(4);
t588 = -t422 * t557 + t423 * t553;
t359 = qJD(6) * t677 + t588;
t671 = pkin(4) * t549;
t668 = qJD(2) * pkin(2);
t585 = t555 * t602;
t448 = qJD(1) * t585 + qJD(3) * t537 + t524 * t620;
t665 = t448 * t554;
t664 = t448 * t558;
t662 = t473 * t554;
t661 = t477 * t554;
t659 = t518 * t539;
t658 = t520 * t539;
t657 = t539 * t558;
t402 = t549 * t408;
t656 = t550 * t556;
t655 = t550 * t560;
t652 = t554 * t555;
t651 = t554 * t559;
t561 = qJD(3) ^ 2;
t649 = t555 * t561;
t368 = t389 * t551 - t402;
t355 = -pkin(5) * t539 + t368 + t681;
t648 = t557 * t355;
t644 = t559 * t561;
t599 = t559 * t614;
t439 = t512 * t618 + t549 * t600 - t551 * t599;
t643 = -pkin(5) * t621 - pkin(10) * t439 - t641;
t438 = -qJD(3) * t567 + t574 * t618;
t642 = pkin(10) * t438 + t640;
t576 = -t512 * t553 - t557 * t574;
t639 = qJD(6) * t576 + t553 * t680 + t557 * t679;
t454 = t512 * t557 - t553 * t574;
t638 = qJD(6) * t454 + t553 * t679 - t557 * t680;
t375 = t407 * t551 - t402;
t514 = t558 * t527;
t459 = -qJ(5) * t650 + t514 + (-pkin(4) - t670) * t559;
t631 = t527 * t554 + t541;
t466 = -qJ(5) * t652 + t631;
t410 = t459 * t549 + t466 * t551;
t634 = -pkin(5) * t680 + t672;
t528 = t669 * t554;
t529 = t669 * t558;
t471 = t528 * t549 - t529 * t551;
t630 = pkin(4) * t652 + pkin(8) * t555;
t610 = qJD(6) * t647 + t422 * t553 + t423 * t557;
t608 = pkin(4) * t678 + pkin(8) * t620;
t607 = -pkin(4) * t558 - pkin(3);
t603 = t556 * t625;
t598 = t539 * t619;
t592 = MDP(16) * t621;
t374 = -t407 * t549 - t654;
t409 = t459 * t551 - t466 * t549;
t470 = t528 * t551 + t529 * t549;
t525 = -t605 - t668;
t583 = -t525 - t605;
t444 = -pkin(10) * t574 + t471;
t582 = pkin(5) * t624 + pkin(10) * t679 + qJD(6) * t444 - t637;
t443 = -pkin(10) * t512 + t470;
t581 = pkin(10) * t680 + qJD(6) * t443 + t636;
t351 = t553 * t355 + t557 * t360;
t495 = t574 * t555;
t383 = -pkin(5) * t559 + pkin(10) * t495 + t409;
t494 = t512 * t555;
t384 = -pkin(10) * t494 + t410;
t580 = t383 * t553 + t384 * t557;
t503 = t559 * t656 + t653;
t464 = -t503 * t554 - t558 * t655;
t569 = -t503 * t558 + t554 * t655;
t412 = t464 * t551 + t549 * t569;
t413 = t464 * t549 - t551 * t569;
t579 = t412 * t557 - t413 * t553;
t578 = t412 * t553 + t413 * t557;
t577 = -t494 * t557 + t495 * t553;
t437 = -t494 * t553 - t495 * t557;
t573 = qJD(2) * t547 - t539 * t559;
t414 = pkin(4) * t478 + t448;
t542 = pkin(4) * t551 + pkin(5);
t571 = t542 * t553 + t557 * t671;
t570 = t542 * t557 - t553 * t671;
t502 = -t552 * t559 + t555 * t656;
t358 = t455 * t615 + t610;
t564 = qJD(3) * (-t583 - t668);
t562 = qJD(2) ^ 2;
t489 = pkin(5) * t574 + t607;
t463 = qJD(3) * t503 + t585;
t462 = -qJD(3) * t502 + t559 * t602;
t460 = pkin(5) * t494 + t630;
t432 = pkin(4) * t520 - pkin(5) * t455;
t411 = -pkin(5) * t438 + t608;
t401 = qJD(4) * t464 + t462 * t558 + t554 * t603;
t400 = qJD(4) * t569 - t462 * t554 + t558 * t603;
t380 = -pkin(5) * t422 + t414;
t379 = qJD(6) * t437 - t438 * t557 - t439 * t553;
t378 = qJD(6) * t577 + t438 * t553 - t439 * t557;
t373 = t400 * t549 + t401 * t551;
t371 = t400 * t551 - t401 * t549;
t362 = t375 + t681;
t361 = t374 - t676;
t350 = -t360 * t553 + t648;
t1 = [(-t400 * t539 + t463 * t518 + t478 * t502) * MDP(17) + (t401 * t539 + t463 * t520 + t477 * t502) * MDP(18) + (t371 * t455 + t373 * t575 - t412 * t423 + t413 * t422) * MDP(19) + (t352 * t412 + t353 * t413 + t368 * t371 + t369 * t373 + t414 * t502 + t442 * t463) * MDP(20) + (-(-qJD(6) * t578 + t371 * t557 - t373 * t553) * t532 - t463 * t398 + t502 * t359) * MDP(26) + ((qJD(6) * t579 + t371 * t553 + t373 * t557) * t532 + t463 * t677 + t502 * t358) * MDP(27) + (-t463 * MDP(10) - t462 * MDP(11) + (MDP(17) * t464 + MDP(18) * t569 + MDP(26) * t579 - MDP(27) * t578) * t624) * qJD(3) + ((-MDP(10) * t555 - MDP(11) * t559) * t560 * t613 + (-t560 * MDP(4) + (-MDP(10) * t559 + MDP(11) * t555 - MDP(3)) * t556) * t562) * t550; 0.2e1 * t594 * t675 - 0.2e1 * t613 * t674 + MDP(7) * t644 - MDP(8) * t649 + (-pkin(8) * t644 + t555 * t564) * MDP(10) + (pkin(8) * t649 + t559 * t564) * MDP(11) + (t477 * t650 + (-t597 + t599) * t520) * MDP(12) + ((-t518 * t558 - t520 * t554) * t620 + (-t661 - t478 * t558 + (t518 * t554 - t520 * t558) * qJD(4)) * t555) * MDP(13) + (t539 * t597 - t477 * t559 + (t520 * t555 + t558 * t573) * qJD(3)) * MDP(14) + (t539 * t596 + t478 * t559 + (-t518 * t555 - t554 * t573) * qJD(3)) * MDP(15) + (-t539 - t623) * t592 + ((t527 * t619 - t688) * t539 + ((pkin(8) * t518 + t662) * qJD(3) + (t660 + (pkin(8) * t539 + t474) * t558) * qJD(4) + t587) * t559 + (-t518 * t605 + t473 * t617 + pkin(8) * t478 + t665 + ((-pkin(8) * t651 + t514) * qJD(2) + t426) * qJD(3)) * t555) * MDP(17) + (t689 * t539 + (t473 * t614 + (qJD(3) * t520 - t598) * pkin(8) + t566) * t559 + (-t520 * t605 - t473 * t619 + pkin(8) * t477 + t664 + (-pkin(8) * t657 - qJD(2) * t631 - t428) * qJD(3)) * t555) * MDP(18) + (t352 * t495 - t353 * t494 + t368 * t439 + t369 * t438 - t409 * t423 + t410 * t422 + t455 * t641 + t575 * t640) * MDP(19) + (t353 * t410 + t352 * t409 + t414 * t630 + (-t555 * t605 + t608) * t442 + t640 * t369 + t641 * t368) * MDP(20) + (t358 * t437 + t378 * t677) * MDP(21) + (t358 * t577 - t359 * t437 + t378 * t398 - t379 * t677) * MDP(22) + (-t358 * t559 - t378 * t532 + (qJD(2) * t437 + t677) * t621) * MDP(23) + (t359 * t559 + t379 * t532 + (qJD(2) * t577 + t398) * t621) * MDP(24) + (-t532 - t623) * MDP(25) * t621 + (-t590 * t559 - t411 * t398 + t460 * t359 - t380 * t577 + t391 * t379 + (t553 * t642 + t557 * t643) * t532 + (t351 * t559 + t532 * t580) * qJD(6) + (t398 * t605 + ((t383 * t557 - t384 * t553) * qJD(2) + t350) * qJD(3)) * t555) * MDP(26) + (((qJD(6) * t355 + t349) * t557 + t589) * t559 + t411 * t677 + t460 * t358 + t380 * t437 + t391 * t378 + ((qJD(6) * t383 + t642) * t557 + (-qJD(6) * t384 - t643) * t553) * t532 + (-t677 * t605 + (-qJD(2) * t580 - t351) * qJD(3)) * t555) * MDP(27); (qJD(3) * t484 - t448) * MDP(10) + t583 * t623 * MDP(11) + (-t520 * t657 + t661) * MDP(12) + ((t477 + t659) * t558 + (-t478 + t658) * t554) * MDP(13) + (-t539 * t617 + (t539 * t646 + (-t520 + t622) * t555) * qJD(2)) * MDP(14) + (t598 + (-t539 * t651 + (t518 + t614) * t555) * qJD(2)) * MDP(15) + (-pkin(3) * t478 - t664 + t586 * t539 - t484 * t518 + (pkin(9) * t657 + t662) * qJD(4) + (-t426 * t555 + (-pkin(9) * t621 - t473 * t559) * t554) * qJD(2)) * MDP(17) + (-pkin(3) * t477 + t665 - t635 * t539 - t484 * t520 + (-pkin(9) * t539 * t554 + t473 * t558) * qJD(4) + (-t473 * t646 + (-pkin(9) * t614 + t428) * t555) * qJD(2)) * MDP(18) + (-t352 * t512 - t353 * t574 - t368 * t679 + t369 * t680 + t422 * t471 - t423 * t470 + t455 * t637 + t636 * t575) * MDP(19) + (t352 * t470 + t353 * t471 + t637 * t368 + t636 * t369 + t414 * t607 + t442 * t672) * MDP(20) + (t358 * t454 + t639 * t677) * MDP(21) + (t358 * t576 - t359 * t454 + t398 * t639 - t638 * t677) * MDP(22) + (t489 * t359 - t380 * t576 + t638 * t391 - t398 * t634) * MDP(26) + (t489 * t358 + t380 * t454 + t639 * t391 + t634 * t677) * MDP(27) + (-t639 * MDP(23) + t638 * MDP(24) + (t553 * t581 + t557 * t582) * MDP(26) + (-t553 * t582 + t557 * t581) * MDP(27)) * t532 + (-t525 * MDP(10) + t539 * MDP(16) + (qJD(3) * t454 - t677) * MDP(23) + (qJD(3) * t576 - t398) * MDP(24) + t532 * MDP(25) + ((t443 * t557 - t444 * t553) * qJD(3) - t350) * MDP(26) + (-(t443 * t553 + t444 * t557) * qJD(3) + t351) * MDP(27)) * t624 + (-t559 * t675 + t674) * t562; t520 * t518 * MDP(12) + (-t518 ^ 2 + t520 ^ 2) * MDP(13) + (t477 - t659) * MDP(14) + (-t478 - t658) * MDP(15) + qJD(2) * t592 + (-t428 * t539 - t473 * t520 + t563) * MDP(17) + (-t426 * t539 + t473 * t518 - t566) * MDP(18) + ((t422 * t549 - t423 * t551) * pkin(4) + (t368 - t375) * t575 + (-t369 - t374) * t455) * MDP(19) + (-t368 * t374 - t369 * t375 + (t352 * t551 + t353 * t549 - t442 * t520) * pkin(4)) * MDP(20) + (t358 + t666) * MDP(23) + (-t359 - t667) * MDP(24) + (t570 * t595 + (t361 * t557 - t362 * t553) * t532 + t432 * t398 + (t532 * t571 - t351) * qJD(6) + t682) * MDP(26) + (-t571 * t595 - (t361 * t553 + t362 * t557) * t532 - t432 * t677 + (t532 * t570 - t648) * qJD(6) + t691) * MDP(27) + t690; (-t455 ^ 2 - t575 ^ 2) * MDP(19) + (-t368 * t455 - t369 * t575 + t414) * MDP(20) + (t359 - t667) * MDP(26) + (t358 - t666) * MDP(27); (t610 + t666) * MDP(23) + (-t588 - t667) * MDP(24) + (-t351 * t532 + t682) * MDP(26) + (-t350 * t532 + t691) * MDP(27) + (MDP(23) * t663 - MDP(24) * t677 - MDP(26) * t351 - MDP(27) * t648) * qJD(6) + t690;];
tauc  = t1;
