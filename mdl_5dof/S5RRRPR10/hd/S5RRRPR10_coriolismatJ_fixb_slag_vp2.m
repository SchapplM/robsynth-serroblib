% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR10_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR10_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:26:47
% EndTime: 2019-12-31 21:27:14
% DurationCPUTime: 13.06s
% Computational Cost: add. (24547->701), mult. (59252->1027), div. (0->0), fcn. (65284->10), ass. (0->367)
t480 = sin(pkin(5));
t485 = cos(qJ(2));
t612 = t480 * t485;
t481 = cos(pkin(5));
t483 = sin(qJ(2));
t613 = t480 * t483;
t674 = sin(qJ(3));
t675 = cos(qJ(3));
t429 = t481 * t675 - t613 * t674;
t430 = t481 * t674 + t613 * t675;
t629 = sin(pkin(10));
t630 = cos(pkin(10));
t339 = -t630 * t429 + t430 * t629;
t648 = t339 * mrSges(5,3);
t302 = mrSges(5,2) * t612 - t648;
t742 = -t302 / 0.2e1;
t741 = -t339 / 0.2e1;
t739 = t339 / 0.2e1;
t574 = t629 * pkin(3);
t469 = t574 + pkin(9);
t738 = m(6) * t469 + mrSges(6,3);
t462 = pkin(7) * t613;
t673 = pkin(1) * t485;
t431 = t481 * t673 - t462;
t432 = (pkin(2) * t483 - pkin(8) * t485) * t480;
t350 = -t431 * t674 + t675 * t432;
t352 = t675 * t431 + t674 * t432;
t737 = -t674 * t350 + t675 * t352;
t669 = Ifges(4,4) * t430;
t736 = -Ifges(4,6) * t612 + t669;
t734 = Ifges(5,3) + Ifges(4,3);
t512 = t429 * t629 + t430 * t630;
t649 = t512 * mrSges(5,3);
t482 = sin(qJ(5));
t478 = t482 ^ 2;
t484 = cos(qJ(5));
t479 = t484 ^ 2;
t601 = t478 + t479;
t733 = t469 * t601;
t286 = -t482 * t512 - t484 * t612;
t287 = -t482 * t612 + t484 * t512;
t167 = -mrSges(6,1) * t286 + mrSges(6,2) * t287;
t303 = -mrSges(5,1) * t612 - t649;
t732 = t167 - t303;
t731 = Ifges(6,5) * t484 - Ifges(6,6) * t482;
t474 = Ifges(4,5) * t675;
t547 = -Ifges(4,6) * t674 + t474;
t475 = Ifges(6,4) * t484;
t455 = t482 * Ifges(6,1) + t475;
t730 = -Ifges(6,2) * t482 + t475;
t476 = Ifges(4,4) * t675;
t454 = -Ifges(4,2) * t674 + t476;
t643 = t430 * mrSges(4,3);
t531 = -mrSges(4,1) * t612 - t643;
t729 = t531 + t643;
t589 = Ifges(4,4) * t674;
t453 = Ifges(4,2) * t675 + t589;
t457 = Ifges(4,1) * t674 + t476;
t598 = t675 / 0.2e1;
t728 = -t674 * t453 / 0.2e1 + t457 * t598;
t727 = -t675 * mrSges(4,1) + t674 * mrSges(4,2);
t438 = t629 * t674 - t630 * t675;
t440 = -t629 * t675 - t630 * t674;
t726 = -Ifges(5,5) * t438 + Ifges(5,6) * t440 + t547;
t723 = -t167 / 0.2e1 + t303 / 0.2e1;
t722 = mrSges(6,3) * t601;
t714 = m(5) * pkin(3);
t721 = t629 * t714 - mrSges(5,2);
t595 = t674 * pkin(8);
t519 = -qJ(4) * t674 - t595;
t472 = t675 * qJ(4);
t477 = t675 * pkin(8);
t602 = t477 + t472;
t395 = -t630 * t519 + t602 * t629;
t720 = t629 * t519 + t602 * t630;
t433 = t481 * t483 * pkin(1) + pkin(7) * t612;
t719 = t433 * mrSges(3,1) + t431 * mrSges(3,2);
t446 = -mrSges(6,1) * t484 + mrSges(6,2) * t482;
t576 = t630 * pkin(3);
t470 = -t576 - pkin(4);
t718 = m(6) * t470 - t630 * t714 - mrSges(5,1) + t446;
t717 = t480 ^ 2;
t716 = m(5) / 0.2e1;
t715 = m(6) / 0.2e1;
t713 = -mrSges(6,1) / 0.2e1;
t651 = t287 * Ifges(6,4);
t119 = t286 * Ifges(6,2) + t339 * Ifges(6,6) + t651;
t712 = -t119 / 0.2e1;
t412 = pkin(8) * t481 + t433;
t413 = (-pkin(2) * t485 - pkin(8) * t483 - pkin(1)) * t480;
t310 = -t412 * t674 + t675 * t413;
t256 = -t430 * qJ(4) + t310;
t232 = -pkin(3) * t612 + t256;
t311 = t412 * t675 + t413 * t674;
t257 = t429 * qJ(4) + t311;
t558 = t629 * t257;
t129 = t232 * t630 - t558;
t121 = pkin(4) * t612 - t129;
t711 = t121 / 0.2e1;
t653 = t286 * mrSges(6,3);
t191 = -mrSges(6,2) * t339 + t653;
t710 = -t191 / 0.2e1;
t652 = t287 * mrSges(6,3);
t192 = mrSges(6,1) * t339 - t652;
t709 = -t192 / 0.2e1;
t471 = -pkin(3) * t675 - pkin(2);
t361 = t438 * pkin(4) + t440 * pkin(9) + t471;
t240 = t361 * t484 - t482 * t720;
t708 = t240 / 0.2e1;
t596 = t674 * pkin(3);
t362 = -t440 * pkin(4) + t438 * pkin(9) + t596;
t245 = t362 * t482 - t395 * t484;
t707 = t245 / 0.2e1;
t706 = -t286 / 0.2e1;
t705 = -t287 / 0.2e1;
t704 = t512 / 0.2e1;
t701 = t339 / 0.4e1;
t636 = t484 * mrSges(6,2);
t640 = t482 * mrSges(6,1);
t447 = t636 + t640;
t356 = t447 * t440;
t697 = -t356 / 0.2e1;
t617 = t440 * t482;
t592 = mrSges(6,3) * t617;
t364 = -mrSges(6,2) * t438 + t592;
t696 = t364 / 0.2e1;
t695 = t395 / 0.2e1;
t403 = t440 * t612;
t694 = t403 / 0.2e1;
t693 = -t403 / 0.2e1;
t404 = t438 * t612;
t692 = -t404 / 0.2e1;
t690 = -t438 / 0.2e1;
t689 = -t438 / 0.4e1;
t688 = t438 / 0.4e1;
t687 = -t440 / 0.2e1;
t686 = -t440 / 0.4e1;
t685 = -t446 / 0.2e1;
t684 = t447 / 0.2e1;
t683 = t731 / 0.4e1;
t682 = -t469 / 0.2e1;
t681 = t470 / 0.2e1;
t680 = t482 / 0.2e1;
t678 = -t484 / 0.2e1;
t677 = -t484 / 0.4e1;
t676 = t484 / 0.2e1;
t672 = pkin(3) * t430;
t671 = Ifges(3,4) * t483;
t670 = Ifges(3,4) * t485;
t668 = Ifges(5,4) * t440;
t667 = Ifges(5,5) * t404;
t666 = Ifges(6,5) * t287;
t351 = -t404 * t484 + t482 * t613;
t665 = Ifges(6,5) * t351;
t664 = Ifges(6,5) * t438;
t662 = Ifges(5,6) * t403;
t661 = Ifges(6,6) * t286;
t349 = t404 * t482 + t484 * t613;
t660 = Ifges(6,6) * t349;
t659 = Ifges(6,6) * t438;
t657 = Ifges(6,3) * t512;
t656 = Ifges(6,3) * t403;
t655 = Ifges(6,3) * t440;
t241 = t361 * t482 + t484 * t720;
t654 = t241 * mrSges(6,3);
t247 = t630 * t257;
t130 = t629 * t232 + t247;
t133 = t256 * t629 + t247;
t134 = t256 * t630 - t558;
t137 = Ifges(6,6) * t512 - t339 * t730;
t638 = t482 * Ifges(6,4);
t456 = Ifges(6,1) * t484 - t638;
t138 = Ifges(6,5) * t512 - t339 * t456;
t206 = t447 * t339;
t624 = t339 * t482;
t209 = -mrSges(6,2) * t512 + mrSges(6,3) * t624;
t623 = t339 * t484;
t210 = mrSges(6,1) * t512 + mrSges(6,3) * t623;
t211 = mrSges(5,1) * t339 + mrSges(5,2) * t512;
t524 = t462 + (-pkin(2) - t673) * t481;
t353 = -t429 * pkin(3) + t524;
t511 = -Ifges(5,1) * t512 + Ifges(5,4) * t339 + Ifges(5,5) * t612;
t527 = Ifges(5,4) * t512 - Ifges(5,6) * t612;
t644 = t429 * mrSges(4,3);
t529 = mrSges(4,2) * t612 + t644;
t545 = mrSges(4,1) * t430 + mrSges(4,2) * t429;
t426 = Ifges(4,4) * t429;
t549 = Ifges(4,5) * t612 - t426;
t557 = mrSges(5,1) * t512 - mrSges(5,2) * t339;
t122 = -pkin(9) * t612 + t130;
t149 = t339 * pkin(4) - pkin(9) * t512 + t353;
t57 = -t122 * t482 + t149 * t484;
t58 = t122 * t484 + t149 * t482;
t279 = Ifges(6,4) * t286;
t120 = Ifges(6,1) * t287 + Ifges(6,5) * t339 + t279;
t607 = t484 * t120;
t611 = t482 * t119;
t187 = pkin(4) * t512 + pkin(9) * t339 + t672;
t64 = -t134 * t482 + t187 * t484;
t65 = t134 * t484 + t187 * t482;
t3 = t130 * t649 - t310 * t529 - t524 * t545 - t58 * t209 - t57 * t210 + t121 * t206 - t64 * t192 - t353 * t557 - t134 * t302 + t137 * t706 + t138 * t705 - t65 * t191 - m(5) * (t130 * t134 + t353 * t672) - m(6) * (t57 * t64 + t58 * t65) + (-t666 / 0.2e1 - t661 / 0.2e1 + (-Ifges(5,2) - Ifges(6,3)) * t339 + t527) * t512 - (t731 * t741 - t607 / 0.2e1 + t611 / 0.2e1 + t129 * mrSges(5,3) + t511) * t339 + (t310 * mrSges(4,3) + t549) * t429 + (-pkin(3) * t211 + (-Ifges(4,1) + Ifges(4,2)) * t429 + t736) * t430 + t729 * t311 + (m(5) * t129 - m(6) * t121 - t732) * t133;
t650 = t3 * qJD(1);
t276 = (pkin(3) * t483 - t472 * t485) * t480 + t350;
t578 = t485 * t674;
t552 = t480 * t578;
t309 = -qJ(4) * t552 + t352;
t162 = t276 * t630 - t309 * t629;
t154 = -pkin(4) * t613 - t162;
t163 = t629 * t276 + t630 * t309;
t400 = pkin(3) * t552 + t433;
t498 = t480 * (Ifges(4,6) * t483 + t454 * t485);
t458 = Ifges(4,1) * t675 - t589;
t499 = t480 * (Ifges(4,5) * t483 + t458 * t485);
t500 = t674 * (Ifges(4,2) * t429 + t736);
t448 = mrSges(4,1) * t674 + mrSges(4,2) * t675;
t502 = t448 * t612;
t503 = t675 * (Ifges(4,1) * t430 - t549);
t506 = t480 * (-mrSges(4,2) * t483 - mrSges(4,3) * t578);
t507 = t480 * (-mrSges(4,3) * t485 * t675 + mrSges(4,1) * t483);
t508 = -Ifges(5,4) * t404 + Ifges(5,2) * t403 + Ifges(5,6) * t613;
t509 = -Ifges(5,2) * t339 + t527;
t510 = -Ifges(5,1) * t404 + Ifges(5,4) * t403 + Ifges(5,5) * t613;
t515 = -t656 + t660 + t665;
t516 = Ifges(6,3) * t339 + t661 + t666;
t517 = Ifges(6,4) * t351 + Ifges(6,2) * t349 - Ifges(6,6) * t403;
t518 = Ifges(6,1) * t351 + Ifges(6,4) * t349 - Ifges(6,5) * t403;
t528 = -mrSges(5,2) * t613 + mrSges(5,3) * t403;
t530 = mrSges(5,1) * t613 + mrSges(5,3) * t404;
t541 = mrSges(6,2) * t403 + mrSges(6,3) * t349;
t542 = -mrSges(6,1) * t403 - mrSges(6,3) * t351;
t543 = -t349 * mrSges(6,1) + t351 * mrSges(6,2);
t645 = t404 * mrSges(5,2);
t646 = t403 * mrSges(5,1);
t544 = -t645 - t646;
t567 = -t612 / 0.2e1;
t568 = t613 / 0.2e1;
t569 = -t613 / 0.2e1;
t155 = pkin(9) * t613 + t163;
t218 = -pkin(4) * t403 + pkin(9) * t404 + t400;
t75 = -t155 * t482 + t218 * t484;
t76 = t155 * t484 + t218 * t482;
t4 = ((Ifges(3,2) * t485 + t671) * t568 + (Ifges(3,1) * t483 + t670) * t567) * t480 - (t485 * (-Ifges(3,2) * t483 + t670) + t483 * (Ifges(3,1) * t485 - t671)) * t717 / 0.2e1 - t512 * t510 / 0.2e1 - t429 * t498 / 0.2e1 - t430 * t499 / 0.2e1 + t508 * t739 + t515 * t741 - t433 * (-mrSges(4,1) * t429 + mrSges(4,2) * t430) - t400 * t211 - t351 * t120 / 0.2e1 - t163 * t302 - t162 * t303 - t76 * t191 - t75 * t192 - t154 * t167 + t518 * t705 + t517 * t706 + t349 * t712 + t511 * t692 + t509 * t693 + t516 * t694 + (Ifges(5,3) * t613 + t500 + t662 - t667) * t612 / 0.2e1 - m(5) * (t129 * t162 + t130 * t163 + t353 * t400) - m(6) * (t121 * t154 + t57 * t75 + t58 * t76) - t353 * t544 - t58 * t541 - t57 * t542 - t121 * t543 + (Ifges(4,5) * t430 + Ifges(5,5) * t512 + Ifges(4,6) * t429 - Ifges(5,6) * t339 - t734 * t612) * t569 + (pkin(1) * (mrSges(3,1) * t483 + mrSges(3,2) * t485) + t485 * (Ifges(4,3) * t483 + t547 * t485) / 0.2e1) * t717 + (Ifges(3,5) * t567 + Ifges(3,6) * t568 - (Ifges(3,5) * t485 - Ifges(3,6) * t483) * t480 / 0.2e1 + t719) * t481 - t130 * t528 - t352 * t529 - t129 * t530 - t350 * t531 - m(4) * (t310 * t350 + t311 * t352 + t433 * t524) - t524 * t502 + t503 * t567 - t311 * t506 - t310 * t507;
t647 = t4 * qJD(1);
t642 = t438 * mrSges(5,3);
t641 = t440 * mrSges(5,3);
t634 = t57 * t482;
t633 = t58 * t484;
t166 = mrSges(6,1) * t287 + mrSges(6,2) * t286;
t168 = Ifges(6,5) * t286 - Ifges(6,6) * t287;
t169 = -Ifges(6,2) * t287 + t279;
t170 = Ifges(6,1) * t286 - t651;
t7 = t121 * t166 + t168 * t739 + t57 * t191 - t58 * t192 + (-t58 * mrSges(6,3) + t712 + t170 / 0.2e1) * t287 + (-t57 * mrSges(6,3) + t169 / 0.2e1 + t120 / 0.2e1) * t286;
t632 = t7 * qJD(1);
t631 = t75 * t482;
t538 = -t633 + t634;
t606 = t484 * t191;
t610 = t482 * t192;
t13 = t732 * t512 - (t302 + t606 - t610) * t339 + m(6) * (t121 * t512 + t339 * t538) + m(5) * (-t129 * t512 - t130 * t339);
t628 = qJD(1) * t13;
t625 = t512 * t395;
t622 = t395 * t133;
t620 = t395 * t440;
t619 = t438 * t482;
t618 = t438 * t484;
t616 = t440 * t484;
t615 = t469 * t482;
t614 = t469 * t484;
t298 = -t440 * t730 + t659;
t609 = t482 * t298;
t366 = mrSges(6,1) * t438 + mrSges(6,3) * t616;
t608 = t482 * t366;
t300 = -t440 * t456 + t664;
t605 = t484 * t300;
t604 = t484 * t364;
t599 = t714 / 0.2e1;
t597 = t674 / 0.2e1;
t590 = t672 / 0.2e1;
t582 = -t659 / 0.2e1;
t581 = -t652 / 0.2e1;
t580 = -t649 / 0.2e1;
t579 = mrSges(6,3) * t678;
t573 = t624 / 0.2e1;
t572 = -t623 / 0.2e1;
t571 = t619 / 0.2e1;
t570 = -t618 / 0.2e1;
t564 = t120 / 0.4e1 + t169 / 0.4e1;
t563 = t170 / 0.4e1 - t119 / 0.4e1;
t358 = t455 * t440;
t562 = -t298 / 0.4e1 + t358 / 0.4e1;
t451 = t484 * Ifges(6,2) + t638;
t357 = t451 * t440;
t561 = t300 / 0.4e1 + t357 / 0.4e1;
t560 = t455 / 0.4e1 + t730 / 0.4e1;
t559 = -t456 / 0.4e1 + t451 / 0.4e1;
t370 = -t440 * mrSges(5,1) - t438 * mrSges(5,2);
t556 = mrSges(5,3) * t576;
t555 = mrSges(5,3) * t574;
t553 = t596 / 0.2e1;
t548 = t710 + t653 / 0.2e1;
t546 = t581 + t709;
t449 = t482 * Ifges(6,5) + t484 * Ifges(6,6);
t244 = t362 * t484 + t395 * t482;
t295 = -t438 * t731 - t655;
t296 = Ifges(6,3) * t438 - t440 * t731;
t297 = -t440 * Ifges(6,6) - t438 * t730;
t299 = -t440 * Ifges(6,5) - t438 * t456;
t355 = t447 * t438;
t363 = mrSges(6,2) * t440 + mrSges(6,3) * t619;
t365 = -mrSges(6,1) * t440 + mrSges(6,3) * t618;
t371 = mrSges(5,1) * t438 - mrSges(5,2) * t440;
t437 = Ifges(5,4) * t438;
t372 = t440 * Ifges(5,2) - t437;
t373 = -Ifges(5,2) * t438 - t668;
t374 = -t438 * Ifges(5,1) + t668;
t375 = -Ifges(5,1) * t440 - t437;
t12 = -pkin(2) * t448 - t720 * t356 - t395 * t355 + t241 * t363 + t245 * t364 + t240 * t365 + t244 * t366 + t454 * t598 + t458 * t597 + t371 * t596 + m(6) * (t240 * t244 + t241 * t245 + t395 * t720) + (t373 / 0.2e1 - t374 / 0.2e1 - t296 / 0.2e1 + t299 * t678 + t297 * t680) * t440 + (-t375 / 0.2e1 - t372 / 0.2e1 + t295 / 0.2e1 - t605 / 0.2e1 + t609 / 0.2e1) * t438 + t728 + (m(5) * t596 + t370) * t471;
t486 = (-t339 * t731 + t511 + t611 + t657) * t688 + (-t648 - t206) * t695 - (-t609 / 0.4e1 + Ifges(5,4) * t689 + Ifges(5,1) * t686) * t339 + (t457 + t454) * t429 / 0.4e1 + t503 / 0.4e1 + (-t641 / 0.2e1 + t697) * t133 + (t644 / 0.2e1 - t529 / 0.2e1) * t595 + (t240 * t64 + t241 * t65 + t244 * t57 + t245 * t58 + t622) * t715 - t500 / 0.4e1 - (t605 + t375 + t372) * t339 / 0.4e1 + (-t395 * t130 + t622 + (t353 * t674 + t430 * t471) * pkin(3)) * t716 + t395 * t742 + (t121 * t715 + t580 - t723 + (-t129 + t134) * t716) * t720 + t675 * (-Ifges(4,2) * t430 + t426) / 0.4e1 + t674 * (Ifges(4,1) * t429 - t669) / 0.4e1 + t130 * t641 / 0.2e1 - t138 * t616 / 0.4e1 + t137 * t617 / 0.4e1 + t371 * t590 - t430 * t453 / 0.4e1 + t430 * t458 / 0.4e1 + t471 * t557 / 0.2e1 + t58 * t363 / 0.2e1 + t57 * t365 / 0.2e1 + t64 * t366 / 0.2e1 + t353 * t370 / 0.2e1 + t286 * t297 / 0.4e1 + t287 * t299 / 0.4e1 + t244 * t192 / 0.2e1 + t241 * t209 / 0.2e1 + t191 * t707 + t210 * t708 - t355 * t711 + t607 * t689 + t65 * t696 + t295 * t701 + (-t134 / 0.2e1 + t129 / 0.2e1) * t642 - pkin(2) * t545 / 0.2e1 - t729 * t477 / 0.2e1 - t726 * t612 / 0.4e1 + t524 * t448 / 0.2e1 + t211 * t553 + t516 * t686 + t440 * t509 / 0.4e1 + (-t373 / 0.4e1 + t296 / 0.4e1 - Ifges(5,2) * t689 - Ifges(5,4) * t686 + t374 / 0.4e1) * t512;
t487 = t667 / 0.2e1 - t662 / 0.2e1 - (t162 * t630 + t163 * t629) * t714 / 0.2e1 + mrSges(6,3) * t631 / 0.2e1 - m(6) * (t154 * t470 + (t76 * t484 - t631) * t469) / 0.2e1 - t541 * t614 / 0.2e1 + t542 * t615 / 0.2e1 - t530 * t576 / 0.2e1 - t528 * t574 / 0.2e1 + t76 * t579 - t349 * t451 / 0.4e1 - t351 * t455 / 0.4e1 + t403 * t449 / 0.4e1 + t352 * mrSges(4,2) / 0.2e1 - t350 * mrSges(4,1) / 0.2e1 - t162 * mrSges(5,1) / 0.2e1 + t163 * mrSges(5,2) / 0.2e1 + Ifges(4,6) * t552 / 0.2e1 - t470 * t543 / 0.2e1 + t567 * t474 + t517 * t677 + t154 * t685 - t482 * t518 / 0.4e1 + t734 * t569;
t2 = t487 + t486;
t537 = t2 * qJD(1) + t12 * qJD(2);
t354 = t440 * t446;
t23 = -t395 * t354 + t241 * t366 + ((-t300 / 0.2e1 - t357 / 0.2e1 - t664 / 0.2e1) * t482 + (t358 / 0.2e1 - t298 / 0.2e1 - t654 + t582) * t484) * t440 + (-t364 + t592) * t240;
t488 = (-t240 * mrSges(6,3) / 0.2e1 + t561) * t286 + (-t654 / 0.2e1 + t562) * t287 + (t449 * t701 + t170 * t677 + t484 * t119 / 0.4e1 + (-t634 / 0.2e1 + t633 / 0.2e1) * mrSges(6,3) + (t169 + t120) * t482 / 0.4e1) * t440 + t354 * t711 + t191 * t708 + t241 * t709 + t166 * t695 + t168 * t688 + t57 * t696 - t58 * t366 / 0.2e1;
t497 = t665 / 0.2e1 + t660 / 0.2e1 - t656 / 0.2e1 + t75 * mrSges(6,1) / 0.2e1 - t76 * mrSges(6,2) / 0.2e1;
t5 = t488 - t497;
t536 = t5 * qJD(1) - t23 * qJD(2);
t522 = t604 / 0.2e1 - t608 / 0.2e1;
t523 = t610 / 0.2e1 - t606 / 0.2e1;
t532 = t240 * t482 - t241 * t484;
t489 = -t522 * t339 + (t580 + t723) * t440 + (t648 / 0.2e1 + t742 + t523) * t438 + (t129 * t440 - t130 * t438 - t339 * t720 + t625) * t716 + (-t121 * t440 + t339 * t532 + t438 * t538 + t625) * t715 + t512 * t697;
t492 = t400 * t716 + (t482 * t76 + t484 * t75) * t715 - t646 / 0.2e1 - t645 / 0.2e1 + t541 * t680 + t542 * t676;
t11 = t489 - t492;
t33 = (t356 + t641) * t440 + (-t604 + t608 + t642) * t438 + m(6) * (t438 * t532 - t620) + m(5) * (-t438 * t720 - t620);
t535 = -qJD(1) * t11 - qJD(2) * t33;
t491 = (-t339 * t733 + t470 * t512) * t715 + t446 * t704 + (-t339 * t629 - t512 * t630) * t599 + t741 * t722;
t495 = (t482 * t65 + t484 * t64) * t715 + t209 * t680 + t210 * t676 + m(5) * t590;
t18 = t491 - t495 - t557;
t490 = (-t438 * t733 - t440 * t470) * t715 + t440 * t685 + (-t438 * t629 + t440 * t630) * t599 + t690 * t722;
t494 = (t244 * t484 + t245 * t482) * t715 + t363 * t680 + t365 * t676 + m(5) * t553;
t35 = -t370 + t490 - t494;
t534 = qJD(1) * t18 + qJD(2) * t35;
t24 = (mrSges(6,2) * t739 + t548) * t484 + (mrSges(6,1) * t739 - t546) * t482;
t525 = t636 / 0.2e1 + t640 / 0.2e1;
t514 = t525 * t438;
t59 = t514 - t522;
t533 = qJD(1) * t24 + qJD(2) * t59;
t526 = mrSges(6,2) * t707 + t244 * t713;
t521 = t364 * t682 + t562;
t520 = t366 * t682 + t561;
t496 = t559 * t484 + t560 * t482 + (t479 / 0.2e1 + t478 / 0.2e1) * t469 * mrSges(6,3);
t504 = t354 * t681 + t395 * t684 + t438 * t683;
t22 = (t664 / 0.2e1 + t520) * t484 + (t582 + t521) * t482 + (Ifges(6,3) / 0.2e1 + t496) * t440 + t504 + t526;
t234 = t470 * t447 + (t455 / 0.2e1 + t730 / 0.2e1) * t484 + (t456 / 0.2e1 - t451 / 0.2e1) * t482;
t493 = t121 * t684 + t166 * t681 + t286 * t560 - t287 * t559 + t339 * t683;
t505 = -t657 / 0.2e1 + t64 * t713 + t65 * mrSges(6,2) / 0.2e1;
t9 = (Ifges(6,5) * t739 + t469 * t546 + t564) * t484 + (Ifges(6,6) * t741 + t469 * t548 + t563) * t482 + t493 + t505;
t513 = t9 * qJD(1) + t22 * qJD(2) + t234 * qJD(3);
t60 = t514 + t522;
t37 = t490 + t494;
t25 = t286 * t579 + t339 * t525 + t482 * t581 - t523;
t21 = Ifges(6,5) * t570 + Ifges(6,6) * t571 - t655 / 0.2e1 + t520 * t484 + t521 * t482 + t496 * t440 + t504 - t526;
t19 = t491 + t495;
t10 = t489 + t492;
t8 = Ifges(6,5) * t572 + Ifges(6,6) * t573 + t564 * t484 + t563 * t482 + (t482 * t710 + t192 * t678 + (t286 * t680 + t287 * t678) * mrSges(6,3)) * t469 + t493 - t505;
t6 = t488 + t497;
t1 = -t487 + t486;
t14 = [-qJD(2) * t4 - qJD(3) * t3 + qJD(4) * t13 + qJD(5) * t7, t1 * qJD(3) + t10 * qJD(4) + t6 * qJD(5) - t647 + (m(4) * (-pkin(2) * t433 + t737 * pkin(8)) + t737 * mrSges(4,3) + m(5) * (t163 * t720 + t400 * t471) + t720 * t528 - t719 + m(6) * (t240 * t75 + t241 * t76) - t163 * t642 - t518 * t616 / 0.2e1 + t517 * t617 / 0.2e1 - Ifges(3,6) * t613 - t507 * t595 + t400 * t371 + t76 * t364 + t75 * t366 + t351 * t300 / 0.2e1 - t154 * t356 + t349 * t298 / 0.2e1 + t508 * t690 + t375 * t692 + t296 * t693 + t373 * t694 + (Ifges(4,5) * t674 - Ifges(5,5) * t440 + Ifges(4,6) * t675 - Ifges(5,6) * t438) * t568 + (-m(5) * t162 + m(6) * t154 - t530 + t543) * t395 + (Ifges(3,5) + t728) * t612 + t471 * t544 + t241 * t541 + t240 * t542 + t433 * t727 + t506 * t477 + t499 * t597 + t498 * t598 + t162 * t641 + t510 * t687 - pkin(2) * t502 + t438 * t515 / 0.2e1) * qJD(2), -t650 + t1 * qJD(2) + (Ifges(4,5) * t429 - Ifges(4,6) * t430 + t209 * t614 - t210 * t615 - t470 * t206 + t451 * t573 + t455 * t572 + t449 * t704 + t138 * t680 + t137 * t676 - t310 * mrSges(4,2) - t311 * mrSges(4,1) - (Ifges(5,5) - t556) * t339 + (-Ifges(5,6) - t555) * t512 + t721 * t134 + t718 * t133 + t738 * (-t482 * t64 + t484 * t65)) * qJD(3) + t19 * qJD(4) + t8 * qJD(5), qJD(2) * t10 + qJD(3) * t19 + qJD(5) * t25 + t628, t632 + t6 * qJD(2) + t8 * qJD(3) + t25 * qJD(4) + (-mrSges(6,1) * t58 - mrSges(6,2) * t57 + t168) * qJD(5); qJD(3) * t2 + qJD(4) * t11 + qJD(5) * t5 + t647, qJD(3) * t12 + qJD(4) * t33 - qJD(5) * t23, (t727 * pkin(8) + t297 * t676 + t299 * t680 - t470 * t355 + t363 * t614 - t365 * t615 + t718 * t720 - t721 * t395 + t438 * t556 + t440 * t555 + t449 * t687 + t451 * t571 + t455 * t570 + t726 + t738 * (-t244 * t482 + t245 * t484)) * qJD(3) + t37 * qJD(4) + t21 * qJD(5) + t537, qJD(3) * t37 + qJD(5) * t60 - t535, t21 * qJD(3) + t60 * qJD(4) + (-mrSges(6,1) * t241 - mrSges(6,2) * t240 + t440 * t449) * qJD(5) + t536; -qJD(2) * t2 + qJD(4) * t18 + qJD(5) * t9 + t650, qJD(4) * t35 + qJD(5) * t22 - t537, t234 * qJD(5), t534, (t446 * t469 + t731) * qJD(5) + t513; -qJD(2) * t11 - qJD(3) * t18 - qJD(5) * t24 - t628, -qJD(3) * t35 - qJD(5) * t59 + t535, -t534, 0, -qJD(5) * t447 - t533; -qJD(2) * t5 - qJD(3) * t9 + qJD(4) * t24 - t632, -qJD(3) * t22 + qJD(4) * t59 - t536, -t513, t533, 0;];
Cq = t14;
