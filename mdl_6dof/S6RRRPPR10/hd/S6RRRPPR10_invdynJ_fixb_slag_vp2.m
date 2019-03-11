% Calculate vector of inverse dynamics joint torques for
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPPR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR10_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR10_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR10_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:21:30
% EndTime: 2019-03-09 16:22:53
% DurationCPUTime: 54.66s
% Computational Cost: add. (17211->1044), mult. (41118->1384), div. (0->0), fcn. (32367->14), ass. (0->462)
t383 = cos(qJ(2));
t373 = sin(pkin(6));
t505 = qJD(1) * t373;
t482 = t383 * t505;
t333 = -qJD(3) + t482;
t379 = sin(qJ(2));
t529 = cos(pkin(6));
t461 = t529 * qJD(1);
t447 = pkin(1) * t461;
t287 = pkin(8) * t482 + t379 * t447;
t378 = sin(qJ(3));
t450 = t378 * t482;
t500 = qJD(3) * t378;
t705 = -qJD(4) * t378 - t287 + (-t450 + t500) * pkin(3);
t504 = qJD(1) * t379;
t483 = t373 * t504;
t284 = -pkin(8) * t483 + t383 * t447;
t411 = t373 * (pkin(2) * t379 - pkin(9) * t383);
t285 = qJD(1) * t411;
t382 = cos(qJ(3));
t196 = -t378 * t284 + t382 * t285;
t557 = pkin(3) + qJ(5);
t475 = t557 * t379;
t499 = qJD(3) * t382;
t508 = t382 * t383;
t594 = pkin(4) + pkin(9);
t704 = -(pkin(4) * t508 - t475) * t505 + t196 + t594 * t499;
t703 = qJD(5) * t382 - t705 + t333 * (-qJ(4) * t382 + qJ(5) * t378);
t702 = Ifges(4,4) + Ifges(5,6);
t372 = sin(pkin(11));
t374 = cos(pkin(11));
t660 = t703 * t372 + t374 * t704;
t659 = t372 * t704 - t703 * t374;
t377 = sin(qJ(6));
t381 = cos(qJ(6));
t640 = -t372 * t377 + t374 * t381;
t298 = t640 * qJD(6);
t416 = t461 + qJD(2);
t481 = t382 * t504;
t256 = t373 * t481 + t378 * t416;
t656 = t640 * t256;
t647 = t298 + t656;
t398 = qJD(3) * t416;
t498 = qJD(1) * qJD(2);
t403 = qJDD(1) * t379 + t383 * t498;
t458 = t529 * qJDD(1);
t413 = t458 + qJDD(2);
t167 = t378 * t398 - t382 * t413 + (qJD(3) * t481 + t378 * t403) * t373;
t290 = (-qJDD(1) * t383 + t379 * t498) * t373;
t278 = qJDD(3) + t290;
t117 = t167 * t374 - t278 * t372;
t589 = t117 / 0.2e1;
t118 = t167 * t372 + t278 * t374;
t588 = t118 / 0.2e1;
t392 = t403 * t373;
t516 = t373 * t379;
t492 = t378 * t516;
t448 = qJD(3) * t492;
t166 = qJD(1) * t448 - t378 * t413 + (-t392 - t398) * t382;
t584 = -t166 / 0.2e1;
t582 = -t167 / 0.2e1;
t571 = t278 / 0.2e1;
t672 = Ifges(4,5) - Ifges(5,4);
t671 = Ifges(4,6) - Ifges(5,5);
t670 = Ifges(4,3) + Ifges(5,1);
t669 = Ifges(6,3) + Ifges(4,1);
t509 = t378 * t383;
t409 = t372 * t509 + t374 * t379;
t244 = t409 * t505;
t449 = t382 * t482;
t518 = t372 * t378;
t701 = -pkin(5) * t449 + t244 * pkin(10) + (pkin(5) * t382 - pkin(10) * t518) * qJD(3) + t660;
t410 = -t372 * t379 + t374 * t509;
t243 = t410 * t505;
t700 = t659 + (t374 * t500 - t243) * pkin(10);
t419 = t381 * t372 + t377 * t374;
t299 = t419 * qJD(6);
t404 = t419 * t256;
t646 = -t299 - t404;
t692 = -pkin(8) * t373 * t498 + pkin(1) * t458;
t497 = qJDD(1) * t373;
t694 = pkin(8) * t497 + qJD(2) * t447;
t206 = t379 * t692 + t383 * t694;
t185 = pkin(9) * t413 + t206;
t192 = t290 * pkin(2) + (-qJDD(1) * pkin(1) - pkin(9) * t403) * t373;
t241 = pkin(9) * t416 + t287;
t248 = (-pkin(2) * t383 - pkin(9) * t379 - pkin(1)) * t505;
t65 = -t378 * t185 + t192 * t382 - t241 * t499 - t248 * t500;
t399 = qJDD(4) - t65;
t31 = -pkin(4) * t166 + qJD(5) * t333 - t278 * t557 + t399;
t255 = t378 * t483 - t382 * t416;
t207 = -t379 * t694 + t383 * t692;
t186 = -pkin(2) * t413 - t207;
t384 = t166 * qJ(4) - t256 * qJD(4) + t186;
t33 = t255 * qJD(5) + t167 * t557 + t384;
t10 = t374 * t31 - t33 * t372;
t11 = t372 * t31 + t374 * t33;
t583 = t166 / 0.2e1;
t157 = qJDD(6) - t166;
t585 = t157 / 0.2e1;
t204 = t255 * t372 - t333 * t374;
t459 = t374 * t255 + t333 * t372;
t113 = t204 * t381 + t377 * t459;
t38 = -qJD(6) * t113 + t117 * t381 - t118 * t377;
t604 = t38 / 0.2e1;
t686 = -t204 * t377 + t381 * t459;
t37 = qJD(6) * t686 + t117 * t377 + t118 * t381;
t605 = t37 / 0.2e1;
t240 = -pkin(2) * t416 - t284;
t385 = -t256 * qJ(4) + t240;
t102 = t255 * t557 + t385;
t158 = t241 * t378 - t382 * t248;
t412 = pkin(4) * t256 + t158;
t97 = t333 * t557 + qJD(4) + t412;
t46 = -t102 * t372 + t374 * t97;
t30 = pkin(5) * t256 - pkin(10) * t204 + t46;
t47 = t374 * t102 + t372 * t97;
t40 = pkin(10) * t459 + t47;
t12 = t30 * t381 - t377 * t40;
t5 = -pkin(5) * t166 - pkin(10) * t118 + t10;
t6 = pkin(10) * t117 + t11;
t1 = qJD(6) * t12 + t377 * t5 + t381 * t6;
t13 = t30 * t377 + t381 * t40;
t2 = -qJD(6) * t13 - t377 * t6 + t381 * t5;
t628 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t696 = -t278 / 0.2e1;
t7 = Ifges(7,5) * t37 + Ifges(7,6) * t38 + Ifges(7,3) * t157;
t699 = t628 + t10 * mrSges(6,1) - t11 * mrSges(6,2) + Ifges(5,4) * t696 + 0.2e1 * Ifges(6,5) * t588 + Ifges(7,5) * t605 + 0.2e1 * Ifges(6,6) * t589 + Ifges(7,6) * t604 + Ifges(7,3) * t585 + 0.2e1 * t584 * t669 + (-t583 + t584) * Ifges(5,2) + t7 / 0.2e1 + t702 * t582 + (t672 + Ifges(4,5)) * t571;
t133 = t255 * pkin(3) + t385;
t698 = t240 * mrSges(4,2) - t133 * mrSges(5,3);
t138 = pkin(3) * t333 + qJD(4) + t158;
t697 = t138 * mrSges(5,1) + t12 * mrSges(7,1) - t13 * mrSges(7,2) + t158 * mrSges(4,3);
t581 = t167 / 0.2e1;
t431 = mrSges(5,2) * t382 - mrSges(5,3) * t378;
t371 = pkin(11) + qJ(6);
t366 = sin(t371);
t367 = cos(t371);
t433 = t366 * mrSges(7,1) + t367 * mrSges(7,2);
t437 = mrSges(4,1) * t382 - mrSges(4,2) * t378;
t564 = pkin(5) * t372;
t375 = -pkin(10) - qJ(5);
t678 = m(7) * t375;
t695 = -t378 * (m(7) * t564 + t433) + t382 * (-mrSges(7,3) + t678) + t431 - t437;
t360 = qJ(4) + t564;
t435 = mrSges(6,1) * t372 + mrSges(6,2) * t374;
t693 = -m(7) * t360 - t435;
t159 = t382 * t241 + t378 * t248;
t322 = t333 * qJ(4);
t139 = t322 - t159;
t690 = t139 * mrSges(5,1) - t159 * mrSges(4,3);
t436 = t374 * mrSges(6,1) - t372 * mrSges(6,2);
t364 = pkin(5) * t374 + pkin(4);
t556 = pkin(9) + t364;
t674 = -mrSges(5,1) - mrSges(4,3);
t613 = m(6) * t594 + m(7) * t556 - mrSges(3,2) + t436 - t674;
t630 = -m(6) * qJ(5) - mrSges(6,3);
t621 = mrSges(4,1) - mrSges(5,2) - t630 - t678;
t636 = -m(7) - m(5) - m(6);
t689 = -pkin(3) * t636 + t621;
t685 = Ifges(4,6) * t696 + 0.2e1 * Ifges(5,3) * t581 + t702 * t583 + (t581 - t582) * Ifges(4,2) + (-t671 + Ifges(5,5)) * t571;
t380 = sin(qJ(1));
t566 = cos(qJ(1));
t439 = t529 * t566;
t303 = t379 * t439 + t380 * t383;
t484 = t373 * t566;
t227 = t303 * t378 + t382 * t484;
t463 = t379 * t529;
t305 = -t380 * t463 + t383 * t566;
t514 = t373 * t382;
t231 = t305 * t378 - t380 * t514;
t300 = -t382 * t529 + t492;
t401 = -g(1) * t231 - g(2) * t227 - g(3) * t300;
t684 = -t1 * t419 - t12 * t646 - t13 * t647 - t2 * t640 - t401;
t249 = Ifges(5,6) * t255;
t149 = -t333 * Ifges(5,4) - t256 * Ifges(5,2) + t249;
t251 = qJD(6) + t256;
t576 = t251 / 0.2e1;
t590 = t113 / 0.2e1;
t592 = t686 / 0.2e1;
t250 = Ifges(4,4) * t255;
t663 = t459 * Ifges(6,6);
t665 = t204 * Ifges(6,5);
t624 = t663 + t665;
t633 = -t333 * Ifges(4,5) + t113 * Ifges(7,5) + t686 * Ifges(7,6) + t251 * Ifges(7,3) + t256 * t669 - t250 + t624;
t683 = Ifges(7,5) * t590 + Ifges(7,6) * t592 + Ifges(7,3) * t576 - t149 / 0.2e1 + t633 / 0.2e1 + t697;
t682 = -t290 / 0.2e1;
t681 = t392 / 0.2e1;
t680 = t413 / 0.2e1;
t673 = mrSges(4,2) - mrSges(5,3);
t527 = qJ(4) * t378;
t466 = -pkin(2) - t527;
t311 = -t382 * t557 + t466;
t336 = t594 * t378;
t315 = t374 * t336;
t195 = pkin(5) * t378 + t315 + (pkin(10) * t382 - t311) * t372;
t213 = t374 * t311 + t372 * t336;
t510 = t374 * t382;
t205 = -pkin(10) * t510 + t213;
t106 = t195 * t381 - t205 * t377;
t668 = qJD(6) * t106 + t377 * t701 + t381 * t700;
t107 = t195 * t377 + t205 * t381;
t667 = -qJD(6) * t107 - t377 * t700 + t381 * t701;
t555 = -pkin(10) - t557;
t320 = t555 * t372;
t321 = t555 * t374;
t217 = t320 * t381 + t321 * t377;
t123 = -pkin(4) * t255 + t159;
t116 = t374 * t123;
t528 = qJ(4) * t255;
t137 = t256 * t557 + t528;
t48 = -pkin(5) * t255 + t116 + (-pkin(10) * t256 - t137) * t372;
t525 = t256 * t374;
t69 = t372 * t123 + t374 * t137;
t56 = pkin(10) * t525 + t69;
t662 = -qJD(5) * t640 - qJD(6) * t217 + t377 * t56 - t381 * t48;
t216 = -t320 * t377 + t321 * t381;
t661 = -qJD(5) * t419 + qJD(6) * t216 - t377 * t48 - t381 * t56;
t197 = t382 * t284 + t378 * t285;
t176 = -qJ(4) * t483 - t197;
t144 = -pkin(4) * t450 - t176;
t655 = t243 * pkin(5) - t556 * t500 - t144;
t654 = -t594 * t500 - t144;
t653 = -t255 * t671 + t256 * t672 - t333 * t670;
t160 = t243 * t381 - t244 * t377;
t199 = t299 * t382 + t500 * t640;
t652 = t160 - t199;
t161 = t243 * t377 + t244 * t381;
t279 = t640 * t382;
t198 = -qJD(6) * t279 + t419 * t500;
t651 = t161 - t198;
t302 = t379 * t380 - t383 * t439;
t521 = t302 * t382;
t650 = -pkin(3) * t521 - t302 * t527;
t462 = t383 * t529;
t304 = t379 * t566 + t380 * t462;
t520 = t304 * t382;
t649 = -pkin(3) * t520 - t304 * t527;
t648 = (t449 - t499) * qJ(4) + t705;
t309 = pkin(1) * t462 - pkin(8) * t516;
t645 = -m(6) * qJ(4) - t433 + t693;
t644 = -t133 * (-mrSges(5,2) * t378 - mrSges(5,3) * t382) - t240 * (mrSges(4,1) * t378 + mrSges(4,2) * t382);
t643 = -t378 * t671 + t382 * t672;
t641 = -t166 * t672 - t167 * t671 + t278 * t670;
t64 = t382 * t185 + t378 * t192 - t241 * t500 + t248 * t499;
t639 = -t378 * t65 + t382 * t64;
t53 = -t278 * qJ(4) + t333 * qJD(4) - t64;
t58 = -pkin(3) * t278 + t399;
t638 = t378 * t58 - t382 * t53;
t76 = mrSges(6,2) * t166 + mrSges(6,3) * t117;
t77 = -mrSges(6,1) * t166 - mrSges(6,3) * t118;
t637 = t372 * t76 + t374 * t77;
t421 = t10 * t374 + t11 * t372;
t554 = Ifges(3,4) * t379;
t609 = t373 ^ 2;
t634 = (t379 * (t383 * Ifges(3,1) - t554) / 0.2e1 - pkin(1) * (mrSges(3,1) * t379 + mrSges(3,2) * t383)) * t609;
t119 = -mrSges(6,1) * t459 + mrSges(6,2) * t204;
t210 = mrSges(5,1) * t255 + mrSges(5,3) * t333;
t59 = -mrSges(7,1) * t686 + mrSges(7,2) * t113;
t632 = t119 + t59 - t210;
t103 = qJD(5) + t123 - t322;
t631 = -t103 * t436 + t690;
t455 = mrSges(3,3) * t483;
t627 = -m(4) * t240 + mrSges(3,1) * t416 - mrSges(4,1) * t255 - mrSges(4,2) * t256 - t455;
t626 = t46 * mrSges(6,1) - t47 * mrSges(6,2);
t625 = t240 * mrSges(4,1) - t133 * mrSges(5,2);
t623 = -m(5) * qJ(4) + t645 + t673;
t619 = -mrSges(7,3) - t621;
t512 = t374 * t378;
t618 = mrSges(6,1) * t518 + mrSges(6,2) * t512 + mrSges(3,1) - t695;
t617 = t673 + t693;
t616 = -t65 * mrSges(4,1) + t64 * mrSges(4,2) - t58 * mrSges(5,2) + t53 * mrSges(5,3);
t434 = t367 * mrSges(7,1) - t366 * mrSges(7,2);
t614 = -t434 - t613;
t577 = -t251 / 0.2e1;
t579 = -t204 / 0.2e1;
t580 = -t459 / 0.2e1;
t591 = -t113 / 0.2e1;
t593 = -t686 / 0.2e1;
t611 = -Ifges(6,5) * t579 - Ifges(7,5) * t591 - Ifges(6,6) * t580 - Ifges(7,6) * t593 - Ifges(7,3) * t577 + t626;
t607 = Ifges(7,4) * t605 + Ifges(7,2) * t604 + Ifges(7,6) * t585;
t606 = Ifges(7,1) * t605 + Ifges(7,4) * t604 + Ifges(7,5) * t585;
t42 = Ifges(6,4) * t118 + Ifges(6,2) * t117 - Ifges(6,6) * t166;
t603 = -t42 / 0.2e1;
t43 = Ifges(6,1) * t118 + Ifges(6,4) * t117 - Ifges(6,5) * t166;
t602 = t43 / 0.2e1;
t548 = Ifges(7,4) * t113;
t51 = Ifges(7,2) * t686 + t251 * Ifges(7,6) + t548;
t601 = -t51 / 0.2e1;
t600 = t51 / 0.2e1;
t109 = Ifges(7,4) * t686;
t52 = t113 * Ifges(7,1) + t251 * Ifges(7,5) + t109;
t599 = -t52 / 0.2e1;
t598 = t52 / 0.2e1;
t95 = t204 * Ifges(6,4) + Ifges(6,2) * t459 + Ifges(6,6) * t256;
t596 = -t95 / 0.2e1;
t96 = Ifges(6,1) * t204 + Ifges(6,4) * t459 + Ifges(6,5) * t256;
t595 = -t96 / 0.2e1;
t537 = t256 * Ifges(4,4);
t150 = -t255 * Ifges(4,2) - t333 * Ifges(4,6) + t537;
t586 = -t150 / 0.2e1;
t575 = -t255 / 0.2e1;
t574 = t255 / 0.2e1;
t573 = -t256 / 0.2e1;
t572 = t256 / 0.2e1;
t569 = -t333 / 0.2e1;
t568 = t333 / 0.2e1;
t565 = pkin(1) * t373;
t563 = pkin(9) * t302;
t562 = pkin(9) * t304;
t368 = t382 * pkin(9);
t558 = qJD(3) / 0.2e1;
t513 = t373 * t383;
t310 = pkin(1) * t463 + pkin(8) * t513;
t274 = pkin(9) * t529 + t310;
t507 = pkin(2) * t513 + pkin(9) * t516;
t275 = -t507 - t565;
t286 = qJD(2) * t411;
t288 = t309 * qJD(2);
t105 = -t274 * t499 - t275 * t500 + t382 * t286 - t378 * t288;
t503 = qJD(2) * t383;
t479 = t373 * t503;
t226 = -t448 + (qJD(3) * t529 + t479) * t382;
t74 = t226 * pkin(4) + (-qJD(2) * t475 + qJD(5) * t383) * t373 - t105;
t301 = t378 * t529 + t379 * t514;
t225 = qJD(3) * t301 + t378 * t479;
t289 = t310 * qJD(2);
t386 = -t226 * qJ(4) - t301 * qJD(4) + t289;
t75 = t300 * qJD(5) + t225 * t557 + t386;
t27 = t372 * t74 + t374 * t75;
t553 = Ifges(3,4) * t383;
t552 = Ifges(4,4) * t378;
t551 = Ifges(4,4) * t382;
t550 = Ifges(6,4) * t372;
t549 = Ifges(6,4) * t374;
t547 = Ifges(5,6) * t378;
t546 = Ifges(5,6) * t382;
t536 = t256 * Ifges(5,6);
t526 = t256 * t372;
t524 = t300 * qJ(5);
t523 = t302 * t366;
t522 = t302 * t367;
t517 = t372 * t382;
t515 = t373 * t380;
t190 = -t378 * t274 + t382 * t275;
t355 = pkin(3) * t513;
t174 = -t190 + t355;
t128 = t301 * pkin(4) + qJ(5) * t513 + t174;
t273 = -pkin(2) * t529 - t309;
t291 = t300 * pkin(3);
t460 = t301 * qJ(4) - t291;
t172 = t273 - t460;
t134 = t172 + t524;
t68 = t372 * t128 + t374 * t134;
t191 = t382 * t274 + t378 * t275;
t337 = t382 * pkin(4) + t368;
t506 = t566 * pkin(1) + pkin(8) * t515;
t501 = qJD(3) * t256;
t496 = pkin(9) * t500;
t495 = pkin(9) * t499;
t493 = qJ(4) * t513;
t292 = t302 * pkin(2);
t489 = -t292 + t650;
t294 = t304 * pkin(2);
t488 = -t294 + t649;
t487 = Ifges(3,5) * t392 - Ifges(3,6) * t290 + Ifges(3,3) * t413;
t486 = t305 * pkin(2) + t506;
t480 = qJD(2) * t516;
t478 = t516 / 0.2e1;
t14 = -t38 * mrSges(7,1) + t37 * mrSges(7,2);
t474 = -t505 / 0.2e1;
t473 = t505 / 0.2e1;
t471 = t501 / 0.2e1;
t467 = -pkin(1) * t380 + pkin(8) * t484;
t465 = pkin(9) * t303 - t292;
t464 = pkin(9) * t305 - t294;
t26 = -t372 * t75 + t374 * t74;
t127 = -t166 * mrSges(5,1) + t278 * mrSges(5,2);
t60 = -t117 * mrSges(6,1) + t118 * mrSges(6,2);
t66 = t374 * t128 - t134 * t372;
t228 = t303 * t382 - t378 * t484;
t454 = mrSges(3,3) * t482;
t232 = t305 * t382 + t378 * t515;
t453 = t232 * pkin(3) + t486;
t444 = t383 * t474;
t443 = t383 * t473;
t441 = -t303 * pkin(2) + t467;
t438 = mrSges(4,1) * t300 + mrSges(4,2) * t301;
t432 = -t300 * mrSges(5,2) - t301 * mrSges(5,3);
t430 = Ifges(4,1) * t382 - t552;
t429 = Ifges(6,1) * t372 + t549;
t428 = -Ifges(4,2) * t378 + t551;
t426 = Ifges(6,2) * t374 + t550;
t424 = Ifges(6,5) * t372 + Ifges(6,6) * t374;
t423 = -Ifges(5,2) * t382 + t547;
t422 = Ifges(5,3) * t378 - t546;
t224 = t300 * t372 - t374 * t513;
t49 = pkin(5) * t301 - pkin(10) * t224 + t66;
t223 = t300 * t374 + t372 * t513;
t54 = pkin(10) * t223 + t68;
t15 = -t377 * t54 + t381 * t49;
t16 = t377 * t49 + t381 * t54;
t140 = t223 * t381 - t224 * t377;
t141 = t223 * t377 + t224 * t381;
t415 = -pkin(3) * t228 + t441;
t173 = t493 - t191;
t405 = qJ(4) * t231 + t453;
t104 = -t274 * t500 + t275 * t499 + t378 * t286 + t382 * t288;
t400 = -g(1) * t232 - g(2) * t228 - g(3) * t301;
t136 = -t300 * pkin(4) - t173;
t393 = -qJ(4) * t227 + t415;
t389 = Ifges(3,6) * t529 + (t383 * Ifges(3,2) + t554) * t373;
t92 = -qJ(4) * t480 + qJD(4) * t513 - t104;
t388 = t373 * t416 * (Ifges(3,5) * t383 - Ifges(3,6) * t379);
t39 = -pkin(4) * t167 + qJDD(5) - t53;
t78 = -t225 * pkin(4) - t92;
t347 = Ifges(3,4) * t482;
t327 = -pkin(3) * t382 + t466;
t308 = (-mrSges(3,1) * t383 + mrSges(3,2) * t379) * t373;
t297 = pkin(5) * t510 + t337;
t283 = -mrSges(3,2) * t416 + t454;
t280 = t419 * t382;
t237 = Ifges(3,1) * t483 + Ifges(3,5) * t416 + t347;
t236 = Ifges(3,6) * qJD(2) + qJD(1) * t389;
t212 = -t311 * t372 + t315;
t211 = mrSges(5,1) * t256 - mrSges(5,2) * t333;
t209 = -mrSges(4,1) * t333 - mrSges(4,3) * t256;
t208 = mrSges(4,2) * t333 - mrSges(4,3) * t255;
t194 = t225 * t372 + t374 * t480;
t193 = t225 * t374 - t372 * t480;
t189 = -mrSges(5,2) * t255 - mrSges(5,3) * t256;
t187 = pkin(3) * t256 + t528;
t177 = -pkin(3) * t483 - t196;
t170 = t231 * t366 + t304 * t367;
t169 = t231 * t367 - t304 * t366;
t147 = -t333 * Ifges(5,5) + t255 * Ifges(5,3) - t536;
t143 = mrSges(6,1) * t256 - mrSges(6,3) * t204;
t142 = -mrSges(6,2) * t256 + mrSges(6,3) * t459;
t126 = mrSges(5,1) * t167 - mrSges(5,3) * t278;
t125 = -mrSges(4,2) * t278 - mrSges(4,3) * t167;
t124 = mrSges(4,1) * t278 + mrSges(4,3) * t166;
t101 = t225 * pkin(3) + t386;
t99 = -pkin(3) * t480 - t105;
t98 = -t223 * pkin(5) + t136;
t93 = -t256 * t364 - t158;
t89 = mrSges(7,1) * t251 - mrSges(7,3) * t113;
t88 = -mrSges(7,2) * t251 + mrSges(7,3) * t686;
t87 = mrSges(4,1) * t167 - mrSges(4,2) * t166;
t86 = -mrSges(5,2) * t167 + mrSges(5,3) * t166;
t83 = -pkin(5) * t459 + t103;
t67 = -t137 * t372 + t116;
t63 = -qJD(6) * t141 + t193 * t381 - t194 * t377;
t62 = qJD(6) * t140 + t193 * t377 + t194 * t381;
t57 = -t193 * pkin(5) + t78;
t55 = t167 * pkin(3) + t384;
t25 = -mrSges(7,2) * t157 + mrSges(7,3) * t38;
t24 = mrSges(7,1) * t157 - mrSges(7,3) * t37;
t21 = -pkin(5) * t117 + t39;
t20 = pkin(10) * t193 + t27;
t19 = pkin(5) * t226 - pkin(10) * t194 + t26;
t4 = -qJD(6) * t16 + t19 * t381 - t20 * t377;
t3 = qJD(6) * t15 + t19 * t377 + t20 * t381;
t8 = [(Ifges(7,1) * t62 + Ifges(7,4) * t63) * t590 + (Ifges(7,1) * t141 + Ifges(7,4) * t140) * t605 + (Ifges(6,1) * t224 + Ifges(6,4) * t223) * t588 + t204 * (Ifges(6,1) * t194 + Ifges(6,4) * t193) / 0.2e1 + (t206 * mrSges(3,3) + Ifges(3,4) * t681 - Ifges(5,4) * t583 - Ifges(4,5) * t584 - Ifges(5,5) * t581 + Ifges(3,2) * t682 + Ifges(3,6) * t680 - Ifges(4,6) * t582 - t571 * t670 + t616 - t641 / 0.2e1) * t513 + (t1 * t140 - t12 * t62 + t13 * t63 - t141 * t2) * mrSges(7,3) + (-t10 * t224 + t11 * t223 + t193 * t47 - t194 * t46) * mrSges(6,3) + m(3) * (pkin(1) ^ 2 * qJDD(1) * t609 + t206 * t310 + t207 * t309 + t287 * t288) + (t207 * t529 - t290 * t565 + t309 * t413) * mrSges(3,1) + (Ifges(3,1) * t392 - Ifges(3,4) * t290 + Ifges(3,5) * t413) * t478 + t459 * (Ifges(6,4) * t194 + Ifges(6,2) * t193) / 0.2e1 + (t609 * qJD(1) * (-Ifges(3,2) * t379 + t553) + t373 * t237) * t503 / 0.2e1 + (t663 / 0.2e1 + t665 / 0.2e1 + t669 * t572 + t672 * t569 - Ifges(5,2) * t573 - Ifges(5,6) * t574 + Ifges(4,4) * t575 + t626 + t683 + t698) * t226 + (t58 * mrSges(5,1) - t65 * mrSges(4,3) + Ifges(4,4) * t582 - Ifges(5,6) * t581 + t699) * t301 + (-m(6) * t405 - mrSges(2,1) * t566 + t380 * mrSges(2,2) - m(5) * (t405 + t562) - m(4) * (t486 + t562) - m(7) * t453 - t170 * mrSges(7,1) - t169 * mrSges(7,2) - m(3) * t506 - t305 * mrSges(3,1) - mrSges(3,3) * t515 + t617 * t231 - t613 * t304 + t619 * t232) * g(2) + t78 * t119 + (Ifges(6,4) * t224 + Ifges(6,2) * t223) * t589 + (-t207 * t516 - t284 * t479 - t287 * t480 - t290 * t310 - t309 * t392) * mrSges(3,3) + (-Ifges(4,4) * t572 + Ifges(5,6) * t573 + Ifges(5,3) * t574 - Ifges(4,2) * t575 + t586 + t147 / 0.2e1 - t671 * t569 + t625 + t690) * t225 + (t653 * t478 + t388 / 0.2e1) * qJD(2) + (Ifges(3,3) * t529 + (Ifges(3,5) * t379 + Ifges(3,6) * t383) * t373) * t680 + (Ifges(3,5) * t529 + (t379 * Ifges(3,1) + t553) * t373) * t681 + t389 * t682 + m(4) * (t104 * t159 - t105 * t158 + t186 * t273 + t190 * t65 + t191 * t64) + m(5) * (t101 * t133 + t138 * t99 + t139 * t92 + t172 * t55 + t173 * t53 + t174 * t58) + m(7) * (t1 * t16 + t12 * t4 + t13 * t3 + t15 * t2 + t21 * t98 + t57 * t83) + m(6) * (t10 * t66 + t103 * t78 + t11 * t68 + t136 * t39 + t26 * t46 + t27 * t47) + (Ifges(7,4) * t62 + Ifges(7,2) * t63) * t592 + (Ifges(7,4) * t141 + Ifges(7,2) * t140) * t604 + t186 * t438 + t55 * t432 + (t53 * mrSges(5,1) - t64 * mrSges(4,3) - Ifges(4,4) * t584 + Ifges(5,6) * t583 + t685) * t300 + (Ifges(7,5) * t62 + Ifges(7,6) * t63) * t576 + (Ifges(7,5) * t141 + Ifges(7,6) * t140) * t585 + t21 * (-mrSges(7,1) * t140 + mrSges(7,2) * t141) + t27 * t142 + t26 * t143 + t136 * t60 - pkin(1) * t308 * t497 + (-t206 * t529 - t310 * t413 - t392 * t565) * mrSges(3,2) + (-m(3) * t467 + t303 * mrSges(3,1) - mrSges(3,3) * t484 - m(6) * t393 + t380 * mrSges(2,1) + mrSges(2,2) * t566 - m(5) * (t393 - t563) - m(4) * (t441 - t563) - m(7) * t415 + t522 * mrSges(7,1) - t523 * mrSges(7,2) - (-t433 + t617) * t227 + t613 * t302 - t619 * t228) * g(1) + t98 * t14 + t3 * t88 + t4 * t89 + t83 * (-mrSges(7,1) * t63 + mrSges(7,2) * t62) + t68 * t76 + t66 * t77 + t57 * t59 + Ifges(2,3) * qJDD(1) + t15 * t24 + t16 * t25 + (Ifges(6,5) * t194 + Ifges(6,6) * t193) * t572 + (Ifges(6,5) * t224 + Ifges(6,6) * t223) * t584 + (t138 * mrSges(5,2) - t158 * mrSges(4,1) + Ifges(4,5) * t572 + Ifges(5,4) * t573 + Ifges(5,5) * t574 + Ifges(4,6) * t575 - t236 / 0.2e1 - t139 * mrSges(5,3) - t159 * mrSges(4,2) + t670 * t569) * t480 + t529 * t487 / 0.2e1 + t288 * t283 + t172 * t86 + t173 * t126 + t174 * t127 + t101 * t189 + t190 * t124 + t191 * t125 + t193 * t95 / 0.2e1 + t194 * t96 / 0.2e1 + t103 * (-mrSges(6,1) * t193 + mrSges(6,2) * t194) + t104 * t208 + t105 * t209 + t92 * t210 + t99 * t211 + (-m(3) * t284 - t627) * t289 + t223 * t42 / 0.2e1 + t39 * (-mrSges(6,1) * t223 + mrSges(6,2) * t224) + t62 * t598 + t63 * t600 + t224 * t602 + t141 * t606 + t140 * t607 + t634 * t498 + t273 * t87; (-Ifges(3,2) * t483 + t382 * t633 + t237 + t347) * t444 + (Ifges(6,4) * t244 + Ifges(6,2) * t243) * t580 + (t256 * (Ifges(5,4) * t379 + t383 * t423) + t255 * (Ifges(4,6) * t379 + t383 * t428) + t379 * t236 + (t379 * t670 + t383 * t643) * t333) * t473 + (-t12 * t449 + t21 * t279 + t652 * t83) * mrSges(7,1) + t551 * t584 - t546 * t583 + (t256 * (Ifges(4,5) * t379 + t383 * t430) + t255 * (Ifges(5,5) * t379 + t383 * t422) + t653 * t379) * t474 - t547 * t581 + t552 * t582 + (t125 - t126) * t368 + (Ifges(6,3) * t573 - t611) * t449 + (Ifges(6,5) * t244 + Ifges(6,6) * t243) * t573 + (-Ifges(7,4) * t280 - Ifges(7,2) * t279) * t604 + (-Ifges(7,1) * t280 - Ifges(7,4) * t279) * t605 + (-Ifges(7,5) * t280 - Ifges(7,6) * t279) * t585 + (-t1 * t279 + t12 * t651 - t13 * t652 + t2 * t280) * mrSges(7,3) + (Ifges(7,4) * t198 + Ifges(7,2) * t199) * t592 + (Ifges(7,4) * t161 + Ifges(7,2) * t160) * t593 + (Ifges(6,1) * t244 + Ifges(6,4) * t243) * t579 + (t372 * t96 + t374 * t95 + t147) * t500 / 0.2e1 + (-t283 + t454) * t284 + (t10 * t517 - t11 * t510 - t243 * t47 + t244 * t46) * mrSges(6,3) + (t643 * t569 + (-t428 / 0.2e1 + t422 / 0.2e1) * t255 + t47 * (-mrSges(6,2) * t382 + mrSges(6,3) * t512) + t46 * (mrSges(6,1) * t382 - mrSges(6,3) * t518) - t644 + ((t158 * t382 - t159 * t378) * m(4) + m(5) * (t138 * t382 + t139 * t378)) * pkin(9)) * qJD(3) + ((t127 - t124) * pkin(9) + t147 * t444 + t150 * t443 + t424 * t471 + (t204 * t429 + t426 * t459) * t558 + t636 * t493 * g(3) + t699) * t378 + t430 * t471 + (Ifges(7,1) * t198 + Ifges(7,4) * t199) * t590 + (Ifges(7,1) * t161 + Ifges(7,4) * t160) * t591 + (t495 - t177) * t211 + t683 * t499 + (t158 * (mrSges(4,1) * t379 - mrSges(4,3) * t508) - t138 * (mrSges(5,1) * t508 + mrSges(5,2) * t379) - t159 * (-mrSges(4,2) * t379 - mrSges(4,3) * t509) - t139 * (mrSges(5,1) * t509 - mrSges(5,3) * t379)) * t505 + t487 - t186 * t437 + t55 * t431 + (Ifges(6,3) * t471 + t149 * t443 + t39 * t436 - t424 * t584 - t426 * t589 - t429 * t588 + t558 * t624 - t685) * t382 + (t13 * t449 - t21 * t280 - t651 * t83) * mrSges(7,2) + (-m(4) * t507 + t308 + t636 * (t382 * t355 + t507) + (-t409 * mrSges(6,1) - t410 * mrSges(6,2) + t630 * t508 + t695 * t383 + (-m(6) * pkin(4) - m(7) * t364 - t434 + t674) * t379) * t373) * g(3) + (Ifges(7,5) * t198 + Ifges(7,6) * t199) * t576 + (Ifges(7,5) * t161 + Ifges(7,6) * t160) * t577 + (-t388 / 0.2e1 - t634 * qJD(1)) * qJD(1) + (-pkin(2) * t186 + pkin(9) * t639 + t158 * t196 - t159 * t197) * m(4) + t106 * t24 + t107 * t25 - pkin(2) * t87 + (-t495 - t196) * t209 + (-t496 - t197) * t208 + (t496 - t176) * t210 + t337 * t60 + t327 * t86 + t667 * t89 + t668 * t88 + (t1 * t107 + t106 * t2 + t12 * t667 + t13 * t668 + t21 * t297 + t655 * t83) * m(7) + t297 * t14 + t660 * t143 + (t10 * t212 + t103 * t654 + t11 * t213 + t337 * t39 + t46 * t660 + t47 * t659) * m(6) + t659 * t142 - t206 * mrSges(3,2) + t207 * mrSges(3,1) + (t455 + t627) * t287 + t212 * t77 + t213 * t76 - t423 * t501 / 0.2e1 + t244 * t595 + t243 * t596 + t198 * t598 + t161 * t599 + t199 * t600 + t160 * t601 + t510 * t603 - t280 * t606 - t279 * t607 + (t586 + t631) * t500 - t103 * (-mrSges(6,1) * t243 + mrSges(6,2) * t244) + t638 * mrSges(5,1) + t639 * mrSges(4,3) + t644 * t482 + (pkin(9) * t638 + t133 * t648 - t138 * t177 - t139 * t176 + t327 * t55) * m(5) + t648 * t189 + (-m(5) * (t464 + t649) - m(4) * t464 - m(7) * t488 - m(6) * (-qJ(5) * t520 + t488) + mrSges(6,3) * t520 + t614 * t305 + t618 * t304) * g(1) + (-m(5) * (t465 + t650) - m(4) * t465 - m(7) * t489 - m(6) * (-qJ(5) * t521 + t489) + mrSges(6,3) * t521 + t614 * t303 + t618 * t302) * g(2) + t654 * t119 + t655 * t59 - t43 * t517 / 0.2e1; (Ifges(7,5) * t404 + Ifges(7,6) * t656) * t577 + (Ifges(7,1) * t404 + Ifges(7,4) * t656) * t591 + (Ifges(7,4) * t404 + Ifges(7,2) * t656) * t593 + (Ifges(7,1) * t640 - Ifges(7,4) * t419) * t605 + t640 * t606 + (Ifges(7,5) * t640 - Ifges(7,6) * t419) * t585 + t21 * (mrSges(7,1) * t419 + mrSges(7,2) * t640) + (Ifges(7,4) * t640 - Ifges(7,2) * t419) * t604 + (-Ifges(7,5) * t299 - Ifges(7,6) * t298) * t576 + (-Ifges(7,1) * t299 - Ifges(7,4) * t298) * t590 + (-Ifges(7,4) * t299 - Ifges(7,2) * t298) * t592 + t404 * t599 + (Ifges(5,2) * t572 - t568 * t672 - t573 * t669 + t611 + t697 + t698) * t255 - t616 + (t208 - t210) * t158 + (t227 * t689 + t228 * t623) * g(2) + (t231 * t689 + t232 * t623) * g(1) + (t149 + t249) * t575 + (Ifges(6,5) * t374 - Ifges(6,6) * t372) * t584 + (Ifges(6,1) * t374 - t550) * t588 + t641 + (t147 - t537) * t573 + (t150 + t536) * t572 + t39 * t435 + t684 * mrSges(7,3) + (t209 - t211) * t159 - t69 * t142 - t67 * t143 - pkin(3) * t127 + t412 * t119 - t419 * t607 + (t103 * t412 - t46 * t67 - t47 * t69 + qJ(4) * t39 - t421 * t557 + (-t372 * t47 - t374 * t46) * qJD(5)) * m(6) - t637 * t557 - t93 * t59 + (-t142 * t372 - t143 * t374) * qJD(5) + (-t126 + t60) * qJ(4) + t360 * t14 + (-Ifges(4,2) * t574 + Ifges(5,3) * t575 + t424 * t573 + t426 * t580 + t429 * t579 - t568 * t671 - t625 - t631) * t256 + t661 * t88 + t662 * t89 + (t1 * t217 + t12 * t662 + t13 * t661 + t2 * t216 + t21 * t360 - t83 * t93) * m(7) + (-pkin(3) * t58 - qJ(4) * t53 - t133 * t187 - t138 * t159 - t139 * t158) * m(5) + (-Ifges(6,2) * t372 + t549) * t589 - t187 * t189 + t216 * t24 + t217 * t25 + t526 * t595 + t525 * t596 - t299 * t598 + t374 * t602 + t372 * t603 + (-m(5) * t139 + m(6) * t103 + m(7) * t83 + t632) * qJD(4) + (-t250 + t633) * t574 + (t46 * t526 - t47 * t525 - t421) * mrSges(6,3) + (-m(5) * t460 + t432 + t438 - m(7) * (t300 * t375 - t291) - m(6) * (-t291 - t524) + t300 * mrSges(6,3) + t645 * t301) * g(3) + t647 * t601 + (mrSges(7,1) * t647 + mrSges(7,2) * t646) * t83; t640 * t24 + t419 * t25 + t646 * t89 + t647 * t88 + t632 * t333 + (t142 * t374 - t143 * t372 + t189) * t256 + t127 + (t333 * t83 - t684) * m(7) + (t401 + t103 * t333 - (t372 * t46 - t374 * t47) * t256 + t421) * m(6) + (t133 * t256 - t139 * t333 + t401 + t58) * m(5) + t637; t113 * t89 - t686 * t88 - t459 * t142 + t204 * t143 + t14 + t60 + (t113 * t12 - t13 * t686 + t21 + t400) * m(7) + (t204 * t46 - t459 * t47 + t39 + t400) * m(6); -t83 * (mrSges(7,1) * t113 + mrSges(7,2) * t686) + (Ifges(7,1) * t686 - t548) * t591 + t51 * t590 + (Ifges(7,5) * t686 - Ifges(7,6) * t113) * t577 - t12 * t88 + t13 * t89 - g(1) * (mrSges(7,1) * t169 - mrSges(7,2) * t170) - g(2) * ((t227 * t367 - t523) * mrSges(7,1) + (-t227 * t366 - t522) * mrSges(7,2)) - g(3) * ((t300 * t367 + t366 * t513) * mrSges(7,1) + (-t300 * t366 + t367 * t513) * mrSges(7,2)) + (t113 * t13 + t12 * t686) * mrSges(7,3) + t7 + (-Ifges(7,2) * t113 + t109 + t52) * t593 + t628;];
tau  = t8;
