% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRR10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR10_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR10_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:07:43
% EndTime: 2019-03-09 04:08:04
% DurationCPUTime: 11.63s
% Computational Cost: add. (25327->664), mult. (52037->894), div. (0->0), fcn. (56270->8), ass. (0->342)
t436 = sin(pkin(10));
t437 = cos(pkin(10));
t611 = sin(qJ(5));
t612 = cos(qJ(5));
t400 = -t611 * t436 + t612 * t437;
t438 = sin(qJ(6));
t440 = cos(qJ(6));
t463 = t436 * t612 + t437 * t611;
t332 = t400 * t438 + t440 * t463;
t487 = t440 * t400 - t438 * t463;
t512 = Ifges(7,5) * t487 - Ifges(7,6) * t332;
t600 = pkin(8) + qJ(4);
t411 = t600 * t436;
t413 = t600 * t437;
t343 = -t411 * t611 + t413 * t612;
t288 = t400 * pkin(9) + t343;
t342 = -t612 * t411 - t413 * t611;
t464 = -pkin(9) * t463 + t342;
t169 = t288 * t440 + t438 * t464;
t668 = -t288 * t438 + t440 * t464;
t727 = -t169 * mrSges(7,1) - t668 * mrSges(7,2);
t31 = t512 + t727;
t730 = t31 * qJD(6);
t441 = cos(qJ(3));
t376 = t463 * t441;
t378 = t400 * t441;
t284 = -t376 * t438 + t378 * t440;
t439 = sin(qJ(3));
t227 = mrSges(7,1) * t439 - t284 * mrSges(7,3);
t656 = -t227 / 0.2e1;
t610 = pkin(3) * t439;
t410 = -qJ(4) * t441 + qJ(2) + t610;
t396 = t437 * t410;
t442 = -pkin(1) - pkin(7);
t494 = -t436 * t442 + pkin(4);
t523 = t437 * t441;
t323 = -pkin(8) * t523 + t439 * t494 + t396;
t428 = t439 * t442;
t349 = t436 * t410 + t437 * t428;
t526 = t436 * t441;
t336 = -pkin(8) * t526 + t349;
t198 = t612 * t323 - t336 * t611;
t155 = -t378 * pkin(9) + t198;
t131 = t439 * pkin(5) + t155;
t199 = t323 * t611 + t336 * t612;
t156 = -t376 * pkin(9) + t199;
t549 = t156 * t440;
t66 = t131 * t438 + t549;
t72 = -t155 * t438 - t549;
t722 = t66 + t72;
t613 = t441 / 0.2e1;
t517 = t441 * t442;
t398 = pkin(4) * t526 - t517;
t334 = t376 * pkin(5) + t398;
t489 = -t440 * t376 - t378 * t438;
t702 = t284 * mrSges(7,1);
t708 = t489 * mrSges(7,2) + t702;
t729 = t334 * t708;
t375 = t463 * t439;
t377 = t400 * t439;
t278 = t440 * t375 + t377 * t438;
t282 = -t375 * t438 + t440 * t377;
t46 = -t282 * mrSges(7,1) + t278 * mrSges(7,2);
t728 = t46 * qJD(6);
t595 = Ifges(7,4) * t332;
t203 = Ifges(7,2) * t487 + t595;
t634 = t332 / 0.2e1;
t320 = Ifges(7,4) * t487;
t204 = Ifges(7,1) * t332 + t320;
t674 = -Ifges(7,2) * t332 + t204 + t320;
t690 = -t487 / 0.2e1;
t705 = -t332 / 0.2e1;
t711 = Ifges(7,1) * t487 - t595;
t726 = t203 * t634 + t674 * t690 + t711 * t705;
t550 = t156 * t438;
t65 = t131 * t440 - t550;
t660 = m(7) / 0.2e1;
t73 = t155 * t440 - t550;
t725 = (-t65 + t73) * t660 + t656;
t596 = Ifges(7,4) * t284;
t145 = Ifges(7,2) * t489 + t439 * Ifges(7,6) + t596;
t425 = -pkin(4) * t437 - pkin(3);
t362 = -pkin(5) * t400 + t425;
t225 = -mrSges(7,2) * t439 + mrSges(7,3) * t489;
t657 = t225 / 0.2e1;
t701 = t332 * mrSges(7,1);
t709 = t487 * mrSges(7,2) + t701;
t710 = Ifges(7,1) * t489 - t596;
t261 = Ifges(7,4) * t489;
t147 = Ifges(7,1) * t284 + t439 * Ifges(7,5) + t261;
t697 = -Ifges(7,2) * t284 + t261;
t715 = t147 + t697;
t724 = t708 * t362 / 0.2e1 + t709 * t334 / 0.2e1 + t715 * t487 / 0.4e1 + t668 * t657 + t674 * t489 / 0.4e1 + (t710 / 0.4e1 - t145 / 0.4e1) * t332 + (t711 / 0.4e1 - t203 / 0.4e1) * t284;
t652 = t278 / 0.2e1;
t723 = -t282 / 0.2e1;
t645 = t282 / 0.2e1;
t641 = -t701 / 0.2e1;
t503 = t702 / 0.2e1;
t721 = mrSges(7,3) * t278;
t718 = t169 * t332;
t658 = m(7) * pkin(5);
t506 = -t658 / 0.2e1;
t570 = t332 * mrSges(7,2);
t575 = t487 * mrSges(7,1);
t513 = t575 / 0.2e1 - t570 / 0.2e1;
t560 = t463 * mrSges(6,2);
t562 = t400 * mrSges(6,1);
t694 = t332 * t438 + t440 * t487;
t714 = -t562 / 0.2e1 + t560 / 0.2e1 + t694 * t506 - t513;
t577 = t284 * mrSges(7,2);
t584 = t489 * mrSges(7,1);
t516 = t584 / 0.2e1 - t577 / 0.2e1;
t565 = t378 * mrSges(6,2);
t567 = t376 * mrSges(6,1);
t518 = t440 * t489;
t521 = t438 * t284;
t712 = t518 + t521;
t713 = t567 / 0.2e1 + t565 / 0.2e1 + t712 * t506 - t516;
t683 = Ifges(7,5) * t489;
t704 = Ifges(7,6) * t284;
t515 = t683 - t704;
t469 = Ifges(7,5) * t723 + Ifges(7,6) * t652;
t414 = pkin(3) * t441 + qJ(4) * t439;
t404 = t437 * t414;
t524 = t437 * t439;
t328 = pkin(8) * t524 + t441 * t494 + t404;
t351 = t436 * t414 + t437 * t517;
t527 = t436 * t439;
t341 = pkin(8) * t527 + t351;
t205 = t612 * t328 - t341 * t611;
t607 = t441 * pkin(5);
t152 = t377 * pkin(9) + t205 + t607;
t206 = t611 * t328 + t612 * t341;
t174 = pkin(9) * t375 + t206;
t74 = t152 * t440 - t174 * t438;
t75 = t152 * t438 + t174 * t440;
t457 = -t469 + t75 * mrSges(7,2) / 0.2e1 - t74 * mrSges(7,1) / 0.2e1;
t486 = Ifges(7,3) * t613 - t457;
t497 = t704 / 0.2e1 - t683 / 0.2e1;
t706 = 0.2e1 * mrSges(7,2);
t688 = -t278 / 0.2e1;
t699 = t225 * t688;
t350 = -t436 * t517 + t404;
t472 = -t350 * t436 + t351 * t437;
t696 = -t284 * t66 - t489 * t65;
t695 = t668 * t487 + t718;
t689 = t489 / 0.2e1;
t434 = t436 ^ 2;
t435 = t437 ^ 2;
t508 = t434 + t435;
t687 = mrSges(5,3) * t508;
t289 = t378 * mrSges(6,1) - t376 * mrSges(6,2);
t676 = t708 + t289;
t337 = mrSges(6,1) * t463 + t400 * mrSges(6,2);
t675 = t709 + t337;
t394 = Ifges(6,4) * t400;
t340 = Ifges(6,1) * t463 + t394;
t673 = -Ifges(6,2) * t463 + t340 + t394;
t672 = Ifges(6,5) * t400 - Ifges(6,6) * t463 + t512;
t671 = -Ifges(6,5) * t376 - Ifges(6,6) * t378 + t515;
t348 = -t428 * t436 + t396;
t406 = -t439 * mrSges(5,2) - mrSges(5,3) * t526;
t408 = t439 * mrSges(5,1) - mrSges(5,3) * t523;
t669 = -m(5) * (t348 * t437 + t349 * t436) - t436 * t406 - t437 * t408;
t563 = t378 * Ifges(6,4);
t272 = -t376 * Ifges(6,2) + t439 * Ifges(6,6) + t563;
t366 = Ifges(6,4) * t376;
t274 = t378 * Ifges(6,1) + t439 * Ifges(6,5) - t366;
t292 = -Ifges(6,2) * t378 - t366;
t293 = -Ifges(6,1) * t376 - t563;
t597 = Ifges(6,4) * t463;
t339 = Ifges(6,2) * t400 + t597;
t345 = -mrSges(6,2) * t439 - t376 * mrSges(6,3);
t564 = t378 * mrSges(6,3);
t347 = mrSges(6,1) * t439 - t564;
t484 = Ifges(6,1) * t400 - t597;
t535 = t343 * t378;
t536 = t342 * t376;
t615 = t439 / 0.4e1;
t638 = t487 / 0.2e1;
t648 = -t489 / 0.2e1;
t654 = -t284 / 0.2e1;
t667 = (t536 / 0.2e1 - t535 / 0.2e1) * mrSges(6,3) + (t292 + t274) * t400 / 0.4e1 + (-t339 / 0.4e1 + t484 / 0.4e1) * t378 + (t638 * t73 + t648 * t668 + t65 * t690 + t705 * t722) * mrSges(7,3) + t722 * t668 * t660 - (t272 / 0.4e1 - t293 / 0.4e1) * t463 + t425 * t289 / 0.2e1 + t398 * t337 / 0.2e1 + t672 * t615 + t342 * t345 / 0.2e1 - t343 * t347 / 0.2e1 - t673 * t376 / 0.4e1 + (t654 * mrSges(7,3) + t725) * t169 + t724;
t161 = t577 - t584;
t202 = t570 - t575;
t621 = t463 / 0.2e1;
t625 = t378 / 0.2e1;
t666 = (t334 * t463 + t362 * t378) * t660 + t161 * t621 + t202 * t625;
t665 = 2 * qJD(3);
t664 = m(5) / 0.2e1;
t663 = -m(6) / 0.2e1;
t662 = m(6) / 0.2e1;
t661 = -m(7) / 0.2e1;
t643 = t284 / 0.2e1;
t630 = -t375 / 0.2e1;
t629 = t375 / 0.2e1;
t628 = -t376 / 0.2e1;
t626 = -t377 / 0.2e1;
t623 = t400 / 0.2e1;
t620 = -t463 / 0.2e1;
t619 = -t436 / 0.2e1;
t618 = t436 / 0.2e1;
t617 = t437 / 0.2e1;
t616 = t439 / 0.2e1;
t614 = -t441 / 0.2e1;
t609 = pkin(5) * t378;
t608 = pkin(5) * t463;
t606 = t65 * mrSges(7,2);
t605 = t66 * mrSges(7,1);
t604 = t72 * mrSges(7,1);
t603 = t73 * mrSges(7,2);
t599 = Ifges(5,4) * t436;
t598 = Ifges(5,4) * t437;
t594 = Ifges(5,2) * t436;
t592 = pkin(5) * qJD(5);
t585 = t278 * mrSges(7,1);
t579 = t282 * mrSges(7,2);
t144 = -Ifges(7,4) * t282 + Ifges(7,2) * t278 + Ifges(7,6) * t441;
t146 = -Ifges(7,1) * t282 + Ifges(7,4) * t278 + Ifges(7,5) * t441;
t160 = -t579 - t585;
t224 = -mrSges(7,2) * t441 + t721;
t226 = mrSges(7,1) * t441 + mrSges(7,3) * t282;
t271 = -Ifges(6,4) * t377 + Ifges(6,2) * t375 + Ifges(6,6) * t441;
t273 = -Ifges(6,1) * t377 + Ifges(6,4) * t375 + Ifges(6,5) * t441;
t566 = t377 * mrSges(6,2);
t568 = t375 * mrSges(6,1);
t290 = -t566 - t568;
t291 = t565 + t567;
t397 = -pkin(4) * t527 + t428;
t333 = -pkin(5) * t375 + t397;
t344 = -mrSges(6,2) * t441 + mrSges(6,3) * t375;
t346 = mrSges(6,1) * t441 + t377 * mrSges(6,3);
t373 = Ifges(5,6) * t441 + (t594 - t598) * t439;
t374 = Ifges(5,5) * t441 + (-Ifges(5,1) * t437 + t599) * t439;
t556 = t437 * mrSges(5,2);
t558 = t436 * mrSges(5,1);
t485 = t556 + t558;
t387 = t485 * t439;
t388 = t485 * t441;
t405 = -mrSges(5,2) * t441 + mrSges(5,3) * t527;
t407 = mrSges(5,1) * t441 + mrSges(5,3) * t524;
t470 = Ifges(6,5) * t626 + Ifges(6,6) * t629;
t555 = t437 * Ifges(5,5);
t557 = t436 * Ifges(5,6);
t3 = m(7) * (t333 * t334 + t65 * t74 + t66 * t75) + m(6) * (t198 * t205 + t199 * t206 + t397 * t398) + t145 * t652 + m(5) * (t348 * t350 + t349 * t351) + t146 * t643 + t273 * t625 + t274 * t626 + t271 * t628 + t272 * t629 + (Ifges(7,5) * t643 + Ifges(7,6) * t689 + Ifges(6,5) * t625 + Ifges(6,6) * t628 + qJ(2) * mrSges(4,1) + t373 * t619 + t374 * t617 + t442 * t387 + (-Ifges(4,4) + t555 / 0.2e1 - t557 / 0.2e1) * t441) * t441 + t144 * t689 + t349 * t405 + t351 * t406 + t348 * t407 + t350 * t408 + t397 * t291 + t398 * t290 + t199 * t344 + t206 * t345 + t198 * t346 + t205 * t347 + t333 * t161 + t334 * t160 + t75 * t225 + t65 * t226 + t74 * t227 + t66 * t224 + (-qJ(2) * mrSges(4,2) + t442 * t388 + (Ifges(4,4) - t555 + t557) * t439 + (-Ifges(4,1) + Ifges(4,2) + Ifges(5,3) + Ifges(6,3) + Ifges(7,3) - Ifges(5,1) * t435 / 0.2e1 - m(5) * t442 ^ 2 + (t598 - t594 / 0.2e1) * t436) * t441 + t469 + t470) * t439 + t147 * t723;
t576 = t3 * qJD(1);
t573 = t487 * mrSges(7,3);
t569 = t332 * mrSges(7,3);
t561 = t400 * mrSges(6,3);
t559 = t463 * mrSges(6,3);
t6 = m(7) * (t334 * t609 + t65 * t72 + t66 * t73) + t161 * t609 + t729 + t710 * t643 + t145 * t654 + t73 * t225 + t72 * t227 + t293 * t625 - t378 * t272 / 0.2e1 + t398 * t289 + t198 * t345 + (-t347 - t564) * t199 + t696 * mrSges(7,3) - (t274 / 0.2e1 + t292 / 0.2e1 - t198 * mrSges(6,3)) * t376 + t671 * t616 + t715 * t689;
t554 = t6 * qJD(1);
t9 = t729 + t515 * t616 + t65 * t225 - t66 * t227 + (t710 / 0.2e1 - t145 / 0.2e1 - t66 * mrSges(7,3)) * t284 + (t147 / 0.2e1 + t697 / 0.2e1 - t65 * mrSges(7,3)) * t489;
t553 = t9 * qJD(1);
t552 = -mrSges(5,1) * t437 + mrSges(5,2) * t436 - mrSges(4,1);
t52 = (t278 * t284 + t282 * t489) * t660 + (t375 * t378 - t376 * t377) * t662;
t551 = qJD(1) * t52;
t548 = t668 * t489;
t547 = t169 * t284;
t541 = t282 * t284;
t543 = t278 * t489;
t451 = (-t541 / 0.2e1 + t543 / 0.2e1) * mrSges(7,3) + t699 + t227 * t723 + t708 * t614;
t17 = t451 - t513;
t546 = t17 * qJD(1);
t26 = t439 * mrSges(4,1) + t441 * mrSges(4,2) + t332 * t225 + t487 * t227 + t463 * t345 + t400 * t347 + mrSges(3,3) + (m(4) + m(3)) * qJ(2) + m(7) * (t332 * t66 + t487 * t65) + m(6) * (t198 * t400 + t199 * t463) - t669;
t545 = t26 * qJD(1);
t544 = t284 * t440;
t542 = t489 * t438;
t538 = t332 * t440;
t537 = t487 * t438;
t532 = t375 * t376;
t530 = t377 * t378;
t528 = t436 * t407;
t525 = t437 * t405;
t522 = t438 * t224;
t520 = t439 * t441;
t519 = t440 * t226;
t507 = qJD(3) * t439;
t505 = t658 / 0.2e1;
t504 = m(5) * t616;
t498 = t709 * t614;
t491 = t508 * qJ(4);
t44 = m(7) * (-t520 + t541 - t543) + m(6) * (-t520 + t530 + t532) + m(5) * (-0.1e1 + t508) * t520;
t465 = t406 * t617 + t408 * t619;
t473 = -t348 * t436 + t349 * t437;
t444 = (-t160 / 0.2e1 - t290 / 0.2e1 + t387 / 0.2e1 + t465) * t441 + (t161 / 0.2e1 + t291 / 0.2e1 + t525 / 0.2e1 - t528 / 0.2e1 + t388 / 0.2e1) * t439 + (t473 * t441 + (t472 - 0.2e1 * t517) * t439) * t664 + (-t198 * t376 + t199 * t378 - t205 * t375 + t206 * t377 - t397 * t441 + t398 * t439) * t662 + (-t278 * t74 + t282 * t75 - t333 * t441 + t334 * t439 - t696) * t660 + t226 * t688 + t227 * t689 + t224 * t645 + t225 * t643 + t346 * t630 + t347 * t628 + t377 * t344 / 0.2e1 + t345 * t625;
t453 = (t342 * t400 + t343 * t463) * t662 + t695 * t660;
t8 = t444 - t453;
t480 = t8 * qJD(1) + t44 * qJD(2);
t447 = (-t278 * t648 - t284 * t645) * mrSges(7,3) + (-t530 / 0.2e1 - t532 / 0.2e1) * mrSges(6,3) + (-t722 * t278 - t378 * t607) * t660 + t699 + t345 * t630 + t347 * t626 + t725 * t282;
t11 = (-t708 / 0.2e1 - t289 / 0.2e1) * t441 + t447 + t714;
t479 = t11 * qJD(1);
t24 = t489 * t225 - t284 * t227 - t376 * t345 - t378 * t347 + m(7) * (-t284 * t65 + t489 * t66) + m(6) * (-t198 * t378 - t199 * t376) + t669 * t441;
t478 = qJD(1) * t24 + qJD(2) * t52;
t55 = (t625 + t544 / 0.2e1 - t542 / 0.2e1) * t658 + t676;
t64 = (t621 + t538 / 0.2e1 - t537 / 0.2e1) * t658 + t675;
t477 = qJD(1) * t55 + qJD(3) * t64;
t77 = t689 * t706 + 0.2e1 * t503;
t82 = t690 * t706 + 0.2e1 * t641;
t476 = qJD(1) * t77 - qJD(3) * t82;
t466 = m(7) * (t438 * t75 + t440 * t74);
t456 = -t205 * mrSges(6,1) / 0.2e1 + t206 * mrSges(6,2) / 0.2e1 - t470;
t1 = (-t522 / 0.2e1 - t519 / 0.2e1 - t466 / 0.2e1 + t666) * pkin(5) + (-Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1) * t441 + t456 + t457 + t667;
t12 = -t202 * t608 - t425 * t337 + t484 * t620 + t339 * t621 - t673 * t400 / 0.2e1 + (t695 - t718) * mrSges(7,3) + (-m(7) * t608 - t709) * t362 - t668 * t573 + t726;
t448 = -t463 * t607 * t660 - (t690 + t638) * t721;
t20 = (-t709 / 0.2e1 - t337 / 0.2e1) * t441 + t448 + t713;
t462 = t1 * qJD(1) + t20 * qJD(2) - t12 * qJD(3);
t23 = -t362 * t709 + t726;
t29 = t498 - t516;
t446 = (-t547 / 0.2e1 - t548 / 0.2e1) * mrSges(7,3) + t169 * t656 + t512 * t615 + t724;
t5 = t446 - t486;
t461 = t5 * qJD(1) + t29 * qJD(2) - t23 * qJD(3);
t445 = (-t284 * t705 + t489 * t638) * mrSges(7,3) + (-t376 * t623 + t378 * t621) * mrSges(6,3) + t473 * t664 + (-t198 * t463 + t199 * t400 - t342 * t378 - t343 * t376) * t662 + (t169 * t489 - t284 * t668 - t332 * t65 + t487 * t66) * t660 + t227 * t705 + t487 * t657 + t345 * t623 + t347 * t620 + t465;
t450 = t397 * t663 + t333 * t661 + t585 / 0.2e1 + t579 / 0.2e1 + t568 / 0.2e1 + t566 / 0.2e1;
t14 = t445 + (-m(5) * t442 / 0.2e1 + t556 / 0.2e1 + t558 / 0.2e1) * t439 + t450;
t36 = -t332 * t569 - t487 * t573 - m(7) * (t169 * t487 - t332 * t668) - m(6) * (-t342 * t463 + t343 * t400) - m(5) * t491 + (-t400 ^ 2 - t463 ^ 2) * mrSges(6,3) - t687;
t452 = (t278 * t332 + t282 * t487) * t660 + (t375 * t463 + t377 * t400) * t662;
t50 = (t661 + t663 + (t434 / 0.2e1 + t435 / 0.2e1 - 0.1e1 / 0.2e1) * m(5)) * t439 + t452;
t460 = qJD(1) * t14 + qJD(2) * t50 - qJD(3) * t36;
t449 = (t440 * t657 + t438 * t656 + (-t521 / 0.2e1 - t518 / 0.2e1) * mrSges(7,3)) * pkin(5) - t497;
t16 = (-t65 / 0.2e1 + t73 / 0.2e1) * mrSges(7,2) + (-t66 / 0.2e1 - t72 / 0.2e1) * mrSges(7,1) + t449 + t497;
t409 = (mrSges(7,1) * t438 + mrSges(7,2) * t440) * pkin(5);
t45 = (t652 + t688) * mrSges(7,2) + (t723 + t645) * mrSges(7,1);
t458 = -qJD(1) * t16 - qJD(2) * t45 + qJD(5) * t409;
t402 = t409 * qJD(6);
t338 = t560 - t562;
t130 = -t463 * t506 + (t537 - t538) * t505;
t104 = (t378 + t542 - t544) * t505;
t83 = t701 / 0.2e1 + t641;
t78 = t503 - t702 / 0.2e1;
t49 = t508 * t504 + t452 + t504 + (m(6) + m(7)) * t616;
t48 = t52 * qJD(4);
t30 = t498 + t516;
t19 = t614 * t675 + t448 - t713;
t18 = t451 + t513;
t15 = -t606 / 0.2e1 - t605 / 0.2e1 - t603 / 0.2e1 + t604 / 0.2e1 + t449 - t497;
t13 = t445 - mrSges(5,2) * t524 / 0.2e1 - mrSges(5,1) * t527 / 0.2e1 + t442 * t504 - t450;
t10 = t614 * t676 + t447 - t714;
t7 = t444 + t453;
t4 = t446 + t486;
t2 = Ifges(6,3) * t613 - t456 + t486 + t666 * pkin(5) + (t466 + t522 + t519) * pkin(5) / 0.2e1 + t667;
t21 = [qJD(2) * t26 + qJD(3) * t3 + qJD(4) * t24 + qJD(5) * t6 + qJD(6) * t9, t545 + 0.2e1 * ((-t278 * t487 + t282 * t332) * t660 + (-t375 * t400 + t377 * t463) * t662) * qJD(2) + t7 * qJD(3) + t48 + t10 * qJD(5) + t18 * qJD(6), t576 + t7 * qJD(2) + t13 * qJD(4) + t2 * qJD(5) + t4 * qJD(6) + (-Ifges(4,5) - t437 * (Ifges(5,1) * t436 + t598) / 0.2e1 + (Ifges(5,2) * t437 + t599) * t618 + t552 * t442) * t507 + ((t169 * t75 + t333 * t362 + t668 * t74) * t660 + (t205 * t342 + t206 * t343 + t397 * t425) * t662 + (-pkin(3) * t428 + qJ(4) * t472) * t664) * t665 + (t472 * mrSges(5,3) + (t525 - t528) * qJ(4) + t203 * t652 + t146 * t634 + t144 * t638 + t273 * t621 + t271 * t623 + t340 * t626 + t339 * t629 + t373 * t617 + t374 * t618 + (Ifges(5,5) * t436 + Ifges(6,5) * t463 + Ifges(7,5) * t332 + Ifges(5,6) * t437 + Ifges(6,6) * t400 + Ifges(7,6) * t487) * t613 + t668 * t226 - Ifges(4,6) * t441 + t425 * t290 + t397 * t338 + pkin(3) * t387 + t362 * t160 + t343 * t344 + t342 * t346 + t333 * t202 + t169 * t224 + t75 * t573 - t205 * t559 + t206 * t561 - t74 * t569 - mrSges(4,2) * t517 + t204 * t723) * qJD(3), qJD(3) * t13 + qJD(5) * t104 + qJD(6) * t78 + t478, t554 + t10 * qJD(2) + t2 * qJD(3) + t104 * qJD(4) + (-t199 * mrSges(6,1) - t198 * mrSges(6,2) - t603 + t604 + t671) * qJD(5) + t15 * qJD(6) + (m(7) * (t438 * t73 + t440 * t72) - t712 * mrSges(7,3)) * t592, t553 + t18 * qJD(2) + t4 * qJD(3) + t78 * qJD(4) + t15 * qJD(5) + (t515 - t605 - t606) * qJD(6); qJD(3) * t8 + qJD(5) * t11 + qJD(6) * t17 + t48 - t545, t44 * qJD(3), t49 * qJD(4) + t19 * qJD(5) + t30 * qJD(6) + (t202 + t338 + t552) * t507 + ((t362 * t439 + t547 + t548) * t660 + (t425 * t439 + t535 - t536) * t662 + (t441 * t491 - t610) * t664) * t665 + t480 + (-t489 * t569 + t284 * t573 + t376 * t559 + t378 * t561 + (-mrSges(4,2) + t687) * t441) * qJD(3), qJD(3) * t49 + t551, t19 * qJD(3) + (-t377 * mrSges(6,1) + t375 * mrSges(6,2) + (-t278 * t438 - t282 * t440) * t658 + t46) * qJD(5) + t728 + t479, t30 * qJD(3) + t46 * qJD(5) + t546 + t728; -qJD(2) * t8 + qJD(4) * t14 + qJD(5) * t1 + qJD(6) * t5 - t576, qJD(4) * t50 + qJD(5) * t20 + qJD(6) * t29 - t480, -qJD(4) * t36 - qJD(5) * t12 - qJD(6) * t23, qJD(5) * t130 + qJD(6) * t83 + t460, t130 * qJD(4) + (-t343 * mrSges(6,1) - t342 * mrSges(6,2) + t672 + t727) * qJD(5) + t730 + (m(7) * (-t169 * t440 + t438 * t668) - t694 * mrSges(7,3)) * t592 + t462, t83 * qJD(4) + t31 * qJD(5) + t461 + t730; -qJD(3) * t14 + qJD(5) * t55 + qJD(6) * t77 - t478, -qJD(3) * t50 - t551, qJD(5) * t64 - qJD(6) * t82 - t460, 0, t477, t476; -qJD(2) * t11 - qJD(3) * t1 - qJD(4) * t55 + qJD(6) * t16 - t554, -qJD(3) * t20 + qJD(6) * t45 - t479, -qJD(4) * t64 - t462, -t477, -t402, -t402 - t458; -qJD(2) * t17 - qJD(3) * t5 - qJD(4) * t77 - qJD(5) * t16 - t553, -t29 * qJD(3) - t45 * qJD(5) - t546, qJD(4) * t82 - t461, -t476, t458, 0;];
Cq  = t21;
