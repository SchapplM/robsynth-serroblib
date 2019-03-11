% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRPRR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:49:48
% EndTime: 2019-03-09 03:49:58
% DurationCPUTime: 7.47s
% Computational Cost: add. (4894->466), mult. (11454->573), div. (0->0), fcn. (9058->12), ass. (0->207)
t511 = qJDD(3) - qJDD(5);
t515 = qJD(3) - qJD(5);
t524 = sin(qJ(6));
t528 = cos(qJ(6));
t520 = sin(pkin(10));
t526 = sin(qJ(3));
t610 = qJD(3) * t526;
t591 = t520 * t610;
t521 = cos(pkin(10));
t645 = cos(qJ(3));
t593 = t645 * t521;
t572 = qJD(1) * t593;
t588 = qJDD(1) * t645;
t597 = qJDD(1) * t526;
t594 = qJD(3) * t572 + t520 * t588 + t521 * t597;
t416 = qJD(1) * t591 - t594;
t474 = t520 * t645 + t526 * t521;
t466 = t474 * qJD(3);
t563 = t520 * t597 - t521 * t588;
t417 = qJD(1) * t466 + t563;
t620 = t520 * t526;
t592 = qJD(1) * t620;
t461 = -t572 + t592;
t463 = t474 * qJD(1);
t525 = sin(qJ(5));
t529 = cos(qJ(5));
t608 = qJD(5) * t529;
t609 = qJD(5) * t525;
t574 = t529 * t416 - t525 * t417 - t461 * t608 + t463 * t609;
t604 = qJD(6) * t528;
t596 = -t524 * t511 - t515 * t604 - t528 * t574;
t605 = qJD(6) * t524;
t654 = t461 * t525 + t529 * t463;
t346 = -t605 * t654 + t596;
t394 = -t515 * t524 + t528 * t654;
t499 = t528 * t511;
t347 = t394 * qJD(6) - t524 * t574 + t499;
t364 = qJD(5) * t654 - t416 * t525 - t529 * t417;
t360 = qJDD(6) + t364;
t618 = t528 * t360;
t656 = -t529 * t461 + t463 * t525;
t661 = qJD(6) + t656;
t585 = t605 * t661 - t618;
t567 = t524 * t656 * t661 + t585;
t628 = t656 * t515;
t629 = t654 * t515;
t631 = t394 * t654;
t632 = t394 * t661;
t392 = t528 * t515 + t524 * t654;
t633 = t392 * t654;
t634 = t392 * t661;
t636 = t346 * t524;
t619 = t524 * t360;
t671 = t661 * t528;
t667 = t661 * t671 + t619;
t679 = -((t347 + t632) * t524 - (t346 - t634) * t528) * MDP(27) + (t394 * t671 + t636) * MDP(26) + (-t631 + t667) * MDP(28) - (t574 + t628) * MDP(21) - (t364 + t629) * MDP(22) - (t567 - t633) * MDP(29) + (t654 ^ 2 - t656 ^ 2) * MDP(20) - t511 * MDP(23) + (MDP(19) * t656 - MDP(30) * t661) * t654;
t639 = pkin(7) + qJ(2);
t484 = t639 * t521;
t476 = qJD(1) * t484;
t458 = t526 * t476;
t483 = t639 * t520;
t475 = qJD(1) * t483;
t420 = -t645 * t475 - t458;
t601 = qJD(4) - t420;
t665 = pkin(5) * t654;
t502 = pkin(2) * t521 + pkin(1);
t480 = -qJD(1) * t502 + qJD(2);
t395 = t461 * pkin(3) - t463 * qJ(4) + t480;
t377 = -pkin(4) * t461 - t395;
t348 = pkin(5) * t656 - pkin(9) * t654 + t377;
t531 = -pkin(3) - pkin(4);
t662 = -pkin(8) * t463 + t601;
t383 = qJD(3) * t531 + t662;
t421 = -t526 * t475 + t645 * t476;
t391 = pkin(8) * t461 + t421;
t518 = qJD(3) * qJ(4);
t384 = t391 + t518;
t362 = t383 * t525 + t384 * t529;
t357 = -pkin(9) * t515 + t362;
t341 = t348 * t528 - t357 * t524;
t664 = t341 * t654;
t342 = t348 * t524 + t357 * t528;
t663 = t342 * t654;
t519 = qJDD(1) * pkin(1);
t506 = qJDD(2) - t519;
t530 = cos(qJ(1));
t510 = g(2) * t530;
t527 = sin(qJ(1));
t655 = g(1) * t527 - t510;
t553 = -t506 + t655;
t514 = pkin(10) + qJ(3);
t507 = sin(t514);
t508 = cos(t514);
t557 = t507 * t525 + t508 * t529;
t433 = t557 * t530;
t599 = qJD(1) * qJD(2);
t649 = qJDD(1) * t639 + t599;
t436 = t649 * t520;
t437 = t649 * t521;
t590 = qJD(3) * t645;
t575 = t645 * t436 + t526 * t437 - t475 * t610 + t476 * t590;
t556 = qJDD(4) + t575;
t354 = pkin(8) * t416 + qJDD(3) * t531 + t556;
t516 = qJDD(3) * qJ(4);
t517 = qJD(3) * qJD(4);
t595 = -t526 * t436 + t645 * t437 - t475 * t590;
t370 = -t476 * t610 + t516 + t517 + t595;
t355 = pkin(8) * t417 + t370;
t549 = t525 * t354 + t529 * t355 + t383 * t608 - t384 * t609;
t431 = t557 * t527;
t624 = t507 * t529;
t568 = g(2) * t431 + g(3) * (-t508 * t525 + t624);
t659 = g(1) * t433 + t377 * t656 - t549 + t568;
t622 = t508 * t527;
t430 = t525 * t622 - t527 * t624;
t621 = t508 * t530;
t623 = t507 * t530;
t432 = t525 * t621 - t529 * t623;
t548 = g(1) * t432 + g(2) * t430 + g(3) * t557;
t552 = -t529 * t354 + t525 * t355 + t383 * t609 + t384 * t608;
t658 = -t377 * t654 + t548 - t552;
t425 = -t526 * t483 + t645 * t484;
t614 = t529 * qJ(4) + t525 * t531;
t640 = g(2) * t527;
t641 = g(1) * t530;
t570 = t640 + t641;
t653 = qJ(2) * qJDD(1);
t413 = t463 * pkin(3) + t461 * qJ(4);
t385 = -pkin(4) * t463 - t413;
t478 = -pkin(9) + t614;
t338 = pkin(5) * t511 + t552;
t547 = -t338 + t548;
t652 = t661 * (-pkin(9) * t656 + qJD(6) * t478 + t385 - t665) + t547;
t651 = t661 * (pkin(9) * t661 + t665) - t547;
t396 = -t483 * t590 + qJD(2) * t593 + (-qJD(2) * t520 - qJD(3) * t484) * t526;
t650 = -qJD(3) * t396 - qJDD(3) * t425 - t507 * t655;
t378 = pkin(8) * t466 + t396;
t397 = t474 * qJD(2) + qJD(3) * t425;
t465 = -t521 * t590 + t591;
t379 = t465 * pkin(8) + t397;
t424 = t483 * t645 + t526 * t484;
t399 = -t474 * pkin(8) + t424;
t473 = -t593 + t620;
t400 = pkin(8) * t473 + t425;
t560 = t399 * t529 - t400 * t525;
t343 = qJD(5) * t560 + t378 * t529 + t379 * t525;
t361 = t383 * t529 - t384 * t525;
t356 = pkin(5) * t515 - t361;
t414 = t473 * pkin(3) - t474 * qJ(4) - t502;
t389 = -pkin(4) * t473 - t414;
t419 = t473 * t525 + t474 * t529;
t558 = t529 * t473 - t474 * t525;
t358 = -pkin(5) * t558 - pkin(9) * t419 + t389;
t369 = t399 * t525 + t400 * t529;
t374 = qJD(5) * t558 - t465 * t529 + t466 * t525;
t580 = -pkin(9) * t511 + qJD(6) * t348 + t549;
t648 = t338 * t419 + t356 * t374 - t369 * t360 - (qJD(6) * t358 + t343) * t661 + t558 * t580 + t641;
t647 = g(3) * t507 + qJD(3) * (t420 + t458) + t508 * t570 - t595;
t646 = t463 ^ 2;
t643 = g(1) * t431;
t638 = MDP(5) * t520;
t637 = qJDD(3) * pkin(3);
t635 = t356 * t419;
t627 = t461 * t463;
t561 = -qJ(4) * t525 + t529 * t531;
t616 = qJD(5) * t561 - t391 * t525 + t529 * t662;
t615 = t614 * qJD(5) + t391 * t529 + t525 * t662;
t613 = t520 ^ 2 + t521 ^ 2;
t611 = qJD(3) * t421;
t607 = qJD(6) * t357;
t606 = qJD(6) * t654;
t598 = qJDD(1) * t521;
t388 = t466 * pkin(3) + t465 * qJ(4) - t474 * qJD(4);
t479 = -pkin(2) * t598 + t506;
t586 = t613 * qJD(1) ^ 2;
t577 = t515 ^ 2;
t576 = t529 * t515;
t573 = 0.2e1 * t613;
t566 = pkin(3) * t508 + qJ(4) * t507;
t564 = t358 * t360 + t643;
t562 = -t607 - t510;
t376 = -pkin(4) * t466 - t388;
t551 = t502 + t566;
t367 = t417 * pkin(3) + t416 * qJ(4) - t463 * qJD(4) + t479;
t550 = t374 * t528 - t419 * t605;
t545 = t568 - t580;
t544 = g(1) * t623 - g(3) * t508 + t507 * t640 - t575;
t349 = -pkin(4) * t417 - t367;
t543 = -pkin(9) * t360 + (t356 + t361) * t661;
t542 = g(1) * t622 - g(2) * t621 - qJD(3) * t397 - qJDD(3) * t424;
t539 = -t478 * t360 + (-t356 - t616) * t661;
t538 = t573 * t599 - t570;
t537 = t395 * t463 + qJDD(4) - t544;
t477 = pkin(5) - t561;
t455 = t461 ^ 2;
t423 = t433 * t528 - t524 * t527;
t422 = -t433 * t524 - t527 * t528;
t411 = t518 + t421;
t410 = -qJD(3) * pkin(3) + t601;
t387 = (t461 - t592) * qJD(3) + t594;
t375 = qJD(5) * t419 - t465 * t525 - t529 * t466;
t371 = t556 - t637;
t345 = pkin(5) * t375 - pkin(9) * t374 + t376;
t344 = qJD(5) * t369 + t378 * t525 - t379 * t529;
t340 = pkin(5) * t364 + pkin(9) * t574 + t349;
t339 = t528 * t340;
t1 = [(-g(2) * t422 - t342 * t375 + t344 * t394 - t560 * t346 + (-(-qJD(6) * t369 + t345) * t661 + (t340 - t607) * t558 - qJD(6) * t635 - t564) * t524 + t648 * t528) * MDP(32) + (-g(2) * t423 - t339 * t558 + t341 * t375 + t344 * t392 - t560 * t347 + (t345 * t661 + (t357 * t558 - t369 * t661 + t635) * qJD(6) + t564) * t528 + t648 * t524) * MDP(31) + (-t360 * t558 + t375 * t661) * MDP(30) + (-t346 * t558 + t375 * t394 + t419 * t618 + t550 * t661) * MDP(28) + (-t419 * t619 + t347 * t558 - t375 * t392 + (-t374 * t524 - t419 * t604) * t661) * MDP(29) + (t367 * t473 + t388 * t461 + t395 * t466 + t414 * t417 + t542) * MDP(15) + (-t417 * t502 + t466 * t480 + t473 * t479 + t542) * MDP(13) + (-g(1) * t430 + g(2) * t432 + t343 * t515 + t349 * t419 + t369 * t511 + t374 * t377 + t376 * t654 - t389 * t574) * MDP(25) + (t374 * t654 - t419 * t574) * MDP(19) + (-t364 * t419 - t374 * t656 - t375 * t654 - t558 * t574) * MDP(20) + (-g(2) * t433 + t344 * t515 - t349 * t558 + t364 * t389 + t375 * t377 + t376 * t656 - t511 * t560 + t643) * MDP(24) + (t375 * t515 - t511 * t558) * MDP(22) + (t367 * t414 + t370 * t425 + t371 * t424 + t395 * t388 + t411 * t396 + t410 * t397 + (-g(1) * t639 - g(2) * t551) * t530 + (g(1) * t551 - g(2) * t639) * t527) * MDP(18) + t655 * MDP(2) + (-qJD(3) * t466 - qJDD(3) * t473) * MDP(11) + (t416 * t473 - t417 * t474 + t461 * t465 - t463 * t466) * MDP(9) + (-t416 * t474 - t463 * t465) * MDP(8) + (-qJD(3) * t465 + qJDD(3) * t474) * MDP(10) + (-t367 * t474 - t388 * t463 + t395 * t465 + t414 * t416 - t650) * MDP(17) + (t416 * t502 - t465 * t480 + t474 * t479 + t650) * MDP(14) + qJDD(1) * MDP(1) + (pkin(1) * t553 + (t613 * t653 + t538) * qJ(2)) * MDP(7) + (t573 * t653 + t538) * MDP(6) + (-t374 * t515 - t419 * t511) * MDP(21) + (-t370 * t473 + t371 * t474 - t396 * t461 + t397 * t463 - t410 * t465 - t411 * t466 - t416 * t424 - t417 * t425 - t570) * MDP(16) + (MDP(4) * t521 - t638) * (t553 + t519) + (t346 * t419 * t528 + t550 * t394) * MDP(26) + t570 * MDP(3) + ((-t392 * t528 - t394 * t524) * t374 + (-t636 - t347 * t528 + (t392 * t524 - t394 * t528) * qJD(6)) * t419) * MDP(27); -MDP(4) * t598 + qJDD(1) * t638 - MDP(6) * t586 + (-qJ(2) * t586 - t553) * MDP(7) + (-t455 - t646) * MDP(16) + (-t410 * t463 + t411 * t461 + t367 - t655) * MDP(18) + (-t364 + t629) * MDP(24) + (t574 - t628) * MDP(25) + (t567 + t633) * MDP(31) + (t631 + t667) * MDP(32) + (MDP(13) + MDP(15)) * (0.2e1 * qJD(3) * t463 + t563) + (-MDP(14) + MDP(17)) * ((t461 + t592) * qJD(3) - t594); (-t385 * t656 - t511 * t561 + t515 * t615 - t658) * MDP(24) + (-t385 * t654 + t511 * t614 + t515 * t616 - t659) * MDP(25) + (-t371 * pkin(3) - g(3) * t566 + t370 * qJ(4) - t395 * t413 - t410 * t421 + t411 * t601 + t570 * (pkin(3) * t507 - qJ(4) * t508)) * MDP(18) + (t477 * t346 + t615 * t394 + t524 * t652 + t539 * t528 - t663) * MDP(32) + (t477 * t347 + t615 * t392 + t539 * t524 - t528 * t652 + t664) * MDP(31) + (-t395 * t461 + t413 * t463 + 0.2e1 * t516 + 0.2e1 * t517 - t647) * MDP(17) + (t461 * t480 + t647) * MDP(14) + qJDD(3) * MDP(12) + t387 * MDP(10) + (pkin(3) * t416 - qJ(4) * t417 + (t411 - t421) * t463 + (t410 - t601) * t461) * MDP(16) + (-t463 * t480 + t544 + t611) * MDP(13) - t563 * MDP(11) + MDP(8) * t627 + (-t413 * t461 - t537 + t611 + 0.2e1 * t637) * MDP(15) + (-t455 + t646) * MDP(9) - t679; (-qJDD(3) + t627) * MDP(15) + t387 * MDP(16) + (-qJD(3) ^ 2 - t646) * MDP(17) + (-qJD(3) * t411 + t537 - t637) * MDP(18) + (-t463 * t656 - t511 * t529 - t525 * t577) * MDP(24) + (-t463 * t654 + t511 * t525 - t529 * t577) * MDP(25) + (-t529 * t347 + (-t463 * t528 + t524 * t576) * t661 + (-t392 * t515 - t604 * t661 - t619) * t525) * MDP(31) + (-t529 * t346 + (t463 * t524 + t528 * t576) * t661 + (-t394 * t515 + t585) * t525) * MDP(32); (-t362 * t515 + t658) * MDP(24) + (-t361 * t515 + t659) * MDP(25) + (-pkin(5) * t347 - t362 * t392 + t543 * t524 - t528 * t651 - t664) * MDP(31) + (-pkin(5) * t346 - t362 * t394 + t524 * t651 + t543 * t528 + t663) * MDP(32) + t679; t394 * t392 * MDP(26) + (-t392 ^ 2 + t394 ^ 2) * MDP(27) + (t596 + t634) * MDP(28) + (-t499 + t632) * MDP(29) + t360 * MDP(30) + (-g(1) * t422 + t342 * t661 - t356 * t394 + t339) * MDP(31) + (g(1) * t423 + t341 * t661 + t356 * t392) * MDP(32) + (-MDP(29) * t606 + MDP(31) * t562 + MDP(32) * t545) * t528 + (-MDP(28) * t606 + (qJD(6) * t515 + t574) * MDP(29) + t545 * MDP(31) + (-t340 - t562) * MDP(32)) * t524;];
tau  = t1;
