% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRP5_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:16:36
% EndTime: 2019-03-09 03:16:48
% DurationCPUTime: 10.54s
% Computational Cost: add. (8404->548), mult. (19821->671), div. (0->0), fcn. (15443->14), ass. (0->223)
t552 = sin(pkin(10));
t554 = cos(pkin(10));
t553 = sin(pkin(9));
t555 = cos(pkin(9));
t559 = sin(qJ(3));
t684 = cos(qJ(3));
t518 = t684 * t553 + t559 * t555;
t690 = t518 * qJD(1);
t481 = qJD(3) * t552 + t554 * t690;
t558 = sin(qJ(5));
t683 = cos(qJ(5));
t494 = t552 * t690;
t639 = qJD(3) * t554;
t697 = -t494 + t639;
t585 = t683 * t697;
t431 = -t481 * t558 + t585;
t719 = t431 ^ 2;
t625 = t684 * t555;
t590 = -t559 * t553 + t625;
t710 = t590 * qJD(1);
t498 = qJD(5) - t710;
t718 = t431 * t498;
t717 = MDP(24) + MDP(26);
t716 = -MDP(25) + MDP(28);
t600 = t558 * t697;
t711 = t481 * t683 + t600;
t685 = t711 ^ 2;
t715 = t498 * t711;
t633 = qJD(1) * qJD(2);
t677 = pkin(7) + qJ(2);
t686 = qJDD(1) * t677 + t633;
t492 = t686 * t553;
t493 = t686 * t555;
t525 = t677 * t553;
t519 = qJD(1) * t525;
t527 = t677 * t555;
t520 = qJD(1) * t527;
t622 = qJD(3) * t684;
t638 = qJD(3) * t559;
t584 = -t492 * t684 - t559 * t493 + t519 * t638 - t520 * t622;
t673 = qJDD(3) * pkin(3);
t418 = qJDD(4) - t584 - t673;
t551 = pkin(9) + qJ(3);
t546 = cos(t551);
t537 = g(3) * t546;
t544 = sin(t551);
t560 = sin(qJ(1));
t561 = cos(qJ(1));
t609 = g(1) * t561 + g(2) * t560;
t576 = t609 * t544 - t537;
t572 = -t418 + t576;
t509 = t518 * qJD(3);
t632 = qJDD(1) * t553;
t603 = -qJDD(1) * t625 + t559 * t632;
t468 = qJD(1) * t509 + t603;
t465 = qJDD(5) + t468;
t517 = t552 * t683 + t558 * t554;
t586 = -t558 * t552 + t683 * t554;
t621 = qJD(5) * t683;
t636 = qJD(5) * t558;
t693 = -t552 * t636 + t554 * t621;
t704 = t586 * t710 - t693;
t605 = t517 * t465 - t498 * t704;
t676 = pkin(8) + qJ(4);
t524 = t676 * t552;
t526 = t676 * t554;
t478 = -t558 * t524 + t526 * t683;
t550 = pkin(10) + qJ(5);
t543 = sin(t550);
t712 = t478 * t465 + t543 * t576;
t507 = t517 * qJD(5);
t641 = -t517 * t710 + t507;
t606 = t586 * t465 - t498 * t641;
t709 = t481 * t710;
t464 = pkin(3) * t690 - qJ(4) * t710;
t591 = t684 * t519 + t559 * t520;
t420 = t552 * t464 - t554 * t591;
t705 = -qJD(4) * t554 + t420;
t674 = qJDD(1) * pkin(1);
t542 = qJDD(2) - t674;
t695 = g(1) * t560 - g(2) * t561;
t599 = -t542 + t695;
t703 = qJD(3) * t690;
t656 = t544 * t561;
t657 = t544 * t560;
t702 = g(1) * t656 + g(2) * t657 - t537;
t541 = t555 * pkin(2) + pkin(1);
t467 = -pkin(3) * t590 - qJ(4) * t518 - t541;
t479 = -t559 * t525 + t527 * t684;
t422 = t554 * t467 - t479 * t552;
t658 = t518 * t554;
t406 = -pkin(4) * t590 - pkin(8) * t658 + t422;
t423 = t552 * t467 + t554 * t479;
t659 = t518 * t552;
t414 = -pkin(8) * t659 + t423;
t700 = t558 * t406 + t683 * t414;
t419 = t554 * t464 + t552 * t591;
t662 = t710 * t554;
t399 = pkin(4) * t690 - pkin(8) * t662 + t419;
t663 = t710 * t552;
t409 = -pkin(8) * t663 + t420;
t587 = -t524 * t683 - t558 * t526;
t699 = -qJD(4) * t586 - qJD(5) * t587 + t558 * t399 + t683 * t409;
t698 = -qJD(4) * t517 - qJD(5) * t478 - t399 * t683 + t558 * t409;
t696 = 0.2e1 * t710;
t692 = -t684 * t525 - t559 * t527;
t508 = t590 * qJD(3);
t578 = t518 * qJDD(1);
t568 = qJD(1) * t508 + t578;
t631 = qJDD(1) * t555;
t407 = -pkin(2) * t631 + t468 * pkin(3) - qJ(4) * t568 - qJD(4) * t690 + t542;
t592 = -t559 * t492 + t493 * t684;
t415 = qJDD(3) * qJ(4) + (qJD(4) - t591) * qJD(3) + t592;
t379 = t554 * t407 - t552 * t415;
t380 = t552 * t407 + t554 * t415;
t691 = -t379 * t552 + t380 * t554;
t689 = qJ(2) * qJDD(1);
t567 = t552 * qJDD(3) + t554 * t568;
t460 = t552 * t568;
t616 = qJDD(3) * t554 - t460;
t383 = -qJD(5) * t585 + t481 * t636 - t558 * t616 - t683 * t567;
t688 = -t383 * t586 - t641 * t711;
t678 = g(3) * t544;
t575 = -t609 * t546 - t678;
t436 = pkin(3) * t509 - qJ(4) * t508 - qJD(4) * t518;
t448 = t590 * qJD(2) + qJD(3) * t692;
t404 = t554 * t436 - t448 * t552;
t660 = t508 * t554;
t387 = pkin(4) * t509 - pkin(8) * t660 + t404;
t405 = t552 * t436 + t554 * t448;
t661 = t508 * t552;
t394 = -pkin(8) * t661 + t405;
t687 = -qJD(5) * t700 + t387 * t683 - t558 * t394;
t499 = t710 ^ 2;
t682 = pkin(5) * t465;
t675 = qJ(6) * t465;
t523 = -qJD(1) * t541 + qJD(2);
t445 = -pkin(3) * t710 - qJ(4) * t690 + t523;
t472 = -t559 * t519 + t684 * t520;
t463 = qJD(3) * qJ(4) + t472;
t412 = t554 * t445 - t463 * t552;
t390 = -pkin(4) * t710 - pkin(8) * t481 + t412;
t413 = t552 * t445 + t554 * t463;
t398 = pkin(8) * t697 + t413;
t370 = t558 * t390 + t398 * t683;
t672 = t370 * t498;
t668 = t431 * t690;
t667 = t711 * t431;
t666 = t711 * t690;
t665 = t468 * t552;
t655 = t546 * t560;
t654 = t546 * t561;
t653 = t554 * t468;
t545 = cos(t550);
t649 = t560 * t545;
t648 = t561 * t543;
t434 = pkin(4) * t663 + t472;
t645 = pkin(5) * t641 + qJ(6) * t704 - qJD(6) * t517 - t434;
t644 = -qJ(6) * t690 - t699;
t643 = pkin(5) * t690 - t698;
t640 = t553 ^ 2 + t555 ^ 2;
t635 = t472 * qJD(3);
t369 = t390 * t683 - t558 * t398;
t634 = qJD(6) - t369;
t540 = pkin(4) * t554 + pkin(3);
t620 = pkin(4) * t552 + t677;
t618 = t640 * qJD(1) ^ 2;
t367 = t468 * pkin(4) - pkin(8) * t567 + t379;
t373 = pkin(8) * t616 + t380;
t615 = t683 * t367 - t558 * t373 - t390 * t636 - t398 * t621;
t449 = qJD(2) * t518 - t525 * t638 + t527 * t622;
t614 = 0.2e1 * t640;
t612 = -g(1) * t657 + g(2) * t656;
t485 = t543 * t655 + t545 * t561;
t487 = t546 * t648 - t649;
t611 = -g(1) * t485 + g(2) * t487;
t486 = t546 * t649 - t648;
t488 = t543 * t560 + t545 * t654;
t610 = g(1) * t486 - g(2) * t488;
t566 = t558 * t567 - t683 * t616;
t384 = qJD(5) * t711 + t566;
t607 = -t517 * t384 - t431 * t704;
t604 = pkin(3) * t546 + qJ(4) * t544;
t424 = pkin(4) * t661 + t449;
t602 = -t412 * t552 + t413 * t554;
t601 = t540 * t546 + t544 * t676;
t450 = pkin(4) * t659 - t692;
t598 = pkin(5) * t545 + qJ(6) * t543 + t540;
t595 = t695 * t546;
t593 = t406 * t683 - t558 * t414;
t583 = t558 * t367 + t683 * t373 + t390 * t621 - t398 * t636;
t582 = t558 * t387 + t683 * t394 + t406 * t621 - t414 * t636;
t581 = t587 * t465 + t545 * t702;
t461 = -qJD(3) * pkin(3) + qJD(4) + t591;
t573 = g(1) * t487 + g(2) * t485 + t543 * t678 + t615;
t571 = t614 * t633 - t609;
t425 = -pkin(4) * t697 + t461;
t381 = -pkin(5) * t431 - qJ(6) * t711 + t425;
t570 = t381 * t711 + qJDD(6) - t573;
t569 = -g(1) * t488 - g(2) * t486 - t545 * t678 + t583;
t391 = -pkin(4) * t616 + t418;
t361 = t384 * pkin(5) + t383 * qJ(6) - qJD(6) * t711 + t391;
t565 = -qJD(5) * t600 - t481 * t621 - t566;
t528 = t561 * t541;
t521 = -qJDD(1) * t541 + qJDD(2);
t466 = -pkin(5) * t586 - qJ(6) * t517 - t540;
t457 = t586 * t518;
t456 = t517 * t518;
t417 = t508 * t517 + t518 * t693;
t416 = t507 * t518 - t508 * t586;
t396 = pkin(5) * t711 - qJ(6) * t431;
t393 = t456 * pkin(5) - t457 * qJ(6) + t450;
t378 = pkin(5) * t590 - t593;
t377 = -qJ(6) * t590 + t700;
t376 = -t383 - t718;
t368 = t498 * qJ(6) + t370;
t366 = -t498 * pkin(5) + t634;
t363 = pkin(5) * t417 + qJ(6) * t416 - qJD(6) * t457 + t424;
t362 = -t509 * pkin(5) - t687;
t360 = qJ(6) * t509 - qJD(6) * t590 + t582;
t359 = qJDD(6) - t615 - t682;
t358 = qJD(6) * t498 + t583 + t675;
t1 = [(-t358 * t590 + t360 * t498 - t361 * t457 - t363 * t711 + t368 * t509 + t377 * t465 + t381 * t416 + t383 * t393 - t611) * MDP(28) + (t383 * t590 - t416 * t498 + t457 * t465 + t509 * t711) * MDP(21) + (-t370 * t509 - t450 * t383 + t391 * t457 - t425 * t416 + t424 * t711 - t465 * t700 - t498 * t582 + t583 * t590 + t611) * MDP(25) + (-t383 * t457 - t416 * t711) * MDP(19) + (t555 * MDP(4) - t553 * MDP(5)) * (t599 + t674) + (-t518 * t468 + t508 * t710 - t509 * t690 + t568 * t590) * MDP(9) + (t405 * t710 - t423 * t468 + t380 * t590 - t413 * t509 + t449 * t481 - t692 * t567 + t418 * t658 + t461 * t660 - g(1) * (t552 * t655 + t554 * t561) - g(2) * (-t552 * t654 + t554 * t560)) * MDP(16) + (-t404 * t710 + t422 * t468 - t379 * t590 + t412 * t509 - t449 * t697 + t692 * t616 + t554 * t595 + (t418 * t518 + t461 * t508 - t609) * t552) * MDP(15) + (t508 * t690 + t518 * t568) * MDP(8) + (-t465 * t590 + t498 * t509) * MDP(23) + (-qJD(3) * t509 + qJDD(3) * t590) * MDP(11) + (-t379 * t658 - t380 * t659 - t404 * t481 + t405 * t697 - t412 * t660 - t413 * t661 - t422 * t567 + t423 * t616 - t612) * MDP(17) + t695 * MDP(2) + (t369 * t509 + t450 * t384 + t391 * t456 + t425 * t417 - t424 * t431 + t593 * t465 + t498 * t687 - t590 * t615 + t610) * MDP(24) + (-t358 * t456 + t359 * t457 + t360 * t431 + t362 * t711 - t366 * t416 - t368 * t417 - t377 * t384 - t378 * t383 - t612) * MDP(27) + (t383 * t456 - t384 * t457 - t416 * t431 - t417 * t711) * MDP(20) + (t359 * t590 + t361 * t456 - t362 * t498 - t363 * t431 - t366 * t509 - t378 * t465 + t381 * t417 + t384 * t393 + t610) * MDP(26) + (t384 * t590 - t417 * t498 + t431 * t509 - t456 * t465) * MDP(22) + (t358 * t377 + t368 * t360 + t361 * t393 + t381 * t363 + t359 * t378 + t366 * t362 - g(1) * (-pkin(5) * t486 - qJ(6) * t485) - g(2) * (pkin(5) * t488 + qJ(6) * t487 + t528) + (-g(1) * t620 - g(2) * t601) * t561 + (-g(1) * (-t541 - t601) - g(2) * t620) * t560) * MDP(29) + t609 * MDP(3) + (-t448 * qJD(3) - t479 * qJDD(3) + t523 * t508 + t521 * t518 - t541 * t568 + t612) * MDP(14) + (pkin(1) * t599 + (t640 * t689 + t571) * qJ(2)) * MDP(7) + (t614 * t689 + t571) * MDP(6) + qJDD(1) * MDP(1) + (-qJD(3) * t449 + qJDD(3) * t692 - t468 * t541 + t509 * t523 - t521 * t590 + t595) * MDP(13) + (-g(2) * t528 + t379 * t422 + t380 * t423 + t412 * t404 + t413 * t405 - t418 * t692 + t461 * t449 + (-g(1) * t677 - g(2) * t604) * t561 + (-g(1) * (-t541 - t604) - g(2) * t677) * t560) * MDP(18) + (qJD(3) * t508 + qJDD(3) * t518) * MDP(10); -MDP(4) * t631 + MDP(5) * t632 - MDP(6) * t618 + (-qJ(2) * t618 - t599) * MDP(7) + (t603 + 0.2e1 * t703) * MDP(13) + (qJD(3) * t696 + t578) * MDP(14) + (-t499 * t552 + t690 * t697 + t653) * MDP(15) + (-t481 * t690 - t499 * t554 - t665) * MDP(16) + ((-qJDD(1) * t658 + t494 * t710 - t639 * t696) * t554 + (-t460 - t709) * t552) * MDP(17) + (t379 * t554 + t380 * t552 - t461 * t690 - t602 * t710 - t695) * MDP(18) + (t607 - t688) * MDP(27) + (t358 * t517 - t359 * t586 + t366 * t641 - t368 * t704 - t381 * t690 - t695) * MDP(29) + t716 * (t605 + t666) + t717 * (t606 + t668); -t499 * MDP(9) + t578 * MDP(10) - t603 * MDP(11) + qJDD(3) * MDP(12) + (t576 + t584 + t635) * MDP(13) + (-t523 * t710 - t575 - t592) * MDP(14) + (-qJ(4) * t665 - pkin(3) * t460 - t472 * t494 + t419 * t710 + (t572 + t635 + t673) * t554) * MDP(15) + (-pkin(3) * t567 - qJ(4) * t653 - t461 * t662 - t472 * t481 - t705 * t710) * MDP(16) + (qJ(4) * t554 * t616 - g(1) * t654 - g(2) * t655 + t412 * t662 + t413 * t663 + t419 * t481 - t705 * t697 - t678 + t691) * MDP(17) + (-t412 * t419 - t413 * t420 - t461 * t472 + t602 * qJD(4) + t572 * pkin(3) + (t575 + t691) * qJ(4)) * MDP(18) + (-t383 * t517 - t704 * t711) * MDP(19) + (t607 + t688) * MDP(20) + (t605 - t666) * MDP(21) + (t606 - t668) * MDP(22) + (-t540 * t384 - t391 * t586 + t641 * t425 + t431 * t434 + t498 * t698 + t581) * MDP(24) + (t540 * t383 + t391 * t517 - t425 * t704 - t434 * t711 + t498 * t699 - t712) * MDP(25) + (-t361 * t586 + t381 * t641 + t384 * t466 - t431 * t645 - t498 * t643 + t581) * MDP(26) + (t358 * t586 + t359 * t517 - t366 * t704 - t368 * t641 + t383 * t587 - t384 * t478 + t431 * t644 + t643 * t711 + t575) * MDP(27) + (-t361 * t517 + t381 * t704 + t383 * t466 + t498 * t644 - t645 * t711 + t712) * MDP(28) + (t358 * t478 - t359 * t587 + t361 * t466 + t645 * t381 + t644 * t368 + t643 * t366 + (-g(3) * t598 - t609 * t676) * t546 + (-g(3) * t676 + t598 * t609) * t544) * MDP(29) + (-(-qJD(4) + t461) * t710 * MDP(15) + (t418 - t702) * MDP(16) + (qJ(4) * t567 + qJD(4) * t481) * MDP(17)) * t552 + (-MDP(13) * t523 - MDP(15) * t412 + MDP(16) * t413 - MDP(23) * t498 - MDP(24) * t369 + MDP(25) * t370 + MDP(26) * t366 - MDP(28) * t368 - MDP(8) * t710 + MDP(9) * t690) * t690; (-t616 - t709) * MDP(15) + (-t697 * t710 + t567) * MDP(16) + (-t481 ^ 2 - t697 ^ 2) * MDP(17) + (t412 * t481 - t413 * t697 - t572) * MDP(18) + (-t685 - t719) * MDP(27) + (-t366 * t711 - t368 * t431 + t361 - t576) * MDP(29) + t716 * (t383 - t718) + t717 * (-t565 + t715); -MDP(19) * t667 + (t685 - t719) * MDP(20) + t376 * MDP(21) + (t565 + t715) * MDP(22) + t465 * MDP(23) + (-t425 * t711 + t573 + t672) * MDP(24) + (t369 * t498 - t425 * t431 - t569) * MDP(25) + (t396 * t431 - t570 + t672 + 0.2e1 * t682) * MDP(26) + (pkin(5) * t383 - qJ(6) * t384 + (t368 - t370) * t711 - (t366 - t634) * t431) * MDP(27) + (0.2e1 * t675 + t381 * t431 + t396 * t711 + (0.2e1 * qJD(6) - t369) * t498 + t569) * MDP(28) + (t358 * qJ(6) - t359 * pkin(5) - t381 * t396 - t366 * t370 - g(1) * (-pkin(5) * t487 + qJ(6) * t488) - g(2) * (-pkin(5) * t485 + qJ(6) * t486) - (-pkin(5) * t543 + qJ(6) * t545) * t678 + t634 * t368) * MDP(29); (-qJDD(5) - t603 - t667 - t703) * MDP(26) + t376 * MDP(27) + (-t498 ^ 2 - t685) * MDP(28) + (-t368 * t498 + t570 - t682) * MDP(29);];
tau  = t1;
