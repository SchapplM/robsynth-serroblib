% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:53:55
% EndTime: 2019-03-09 16:55:27
% DurationCPUTime: 58.58s
% Computational Cost: add. (22765->984), mult. (55333->1289), div. (0->0), fcn. (44515->14), ass. (0->424)
t390 = cos(qJ(2));
t386 = sin(qJ(2));
t379 = sin(pkin(6));
t499 = qJD(1) * t379;
t469 = t386 * t499;
t381 = cos(pkin(6));
t498 = qJD(1) * t381;
t487 = pkin(1) * t498;
t308 = -pkin(8) * t469 + t390 * t487;
t409 = (pkin(2) * t386 - pkin(9) * t390) * t379;
t309 = qJD(1) * t409;
t385 = sin(qJ(3));
t389 = cos(qJ(3));
t223 = -t385 * t308 + t389 * t309;
t383 = -qJ(4) - pkin(9);
t450 = qJD(3) * t383;
t506 = t389 * t390;
t720 = -(pkin(3) * t386 - qJ(4) * t506) * t499 - t223 - qJD(4) * t385 + t389 * t450;
t224 = t389 * t308 + t385 * t309;
t468 = t390 * t499;
t440 = t385 * t468;
t719 = -qJ(4) * t440 - qJD(4) * t389 - t385 * t450 + t224;
t682 = Ifges(6,4) + Ifges(7,4);
t378 = sin(pkin(11));
t380 = cos(pkin(11));
t332 = t378 * t389 + t380 * t385;
t257 = t332 * t468;
t318 = t332 * qJD(3);
t718 = -t257 + t318;
t342 = qJD(3) - t468;
t581 = -t342 / 0.2e1;
t361 = qJD(2) + t498;
t282 = t361 * t389 - t385 * t469;
t283 = t361 * t385 + t389 * t469;
t413 = t282 * t378 + t380 * t283;
t592 = -t413 / 0.2e1;
t448 = t380 * t282 - t283 * t378;
t594 = -t448 / 0.2e1;
t717 = -Ifges(5,4) * t592 - Ifges(5,2) * t594 - Ifges(5,6) * t581;
t683 = Ifges(6,1) + Ifges(7,1);
t681 = Ifges(6,5) + Ifges(7,5);
t680 = Ifges(6,2) + Ifges(7,2);
t679 = Ifges(6,6) + Ifges(7,6);
t656 = t378 * t720 - t719 * t380;
t311 = pkin(8) * t468 + t386 * t487;
t495 = qJD(3) * t385;
t653 = -t311 + (-t440 + t495) * pkin(3);
t264 = -t361 * pkin(2) - t308;
t200 = -t282 * pkin(3) + qJD(4) + t264;
t199 = qJD(5) - t448;
t384 = sin(qJ(5));
t388 = cos(qJ(5));
t169 = t342 * t384 + t388 * t413;
t265 = pkin(9) * t361 + t311;
t297 = (-pkin(2) * t390 - pkin(9) * t386 - pkin(1)) * t379;
t269 = qJD(1) * t297;
t186 = -t265 * t385 + t389 * t269;
t158 = -qJ(4) * t283 + t186;
t147 = pkin(3) * t342 + t158;
t187 = t265 * t389 + t269 * t385;
t159 = qJ(4) * t282 + t187;
t512 = t380 * t159;
t85 = t378 * t147 + t512;
t80 = pkin(10) * t342 + t85;
t96 = -pkin(4) * t448 - pkin(10) * t413 + t200;
t30 = -t384 * t80 + t388 * t96;
t25 = -qJ(6) * t169 + t30;
t21 = pkin(5) * t199 + t25;
t168 = t342 * t388 - t384 * t413;
t31 = t384 * t96 + t388 * t80;
t26 = qJ(6) * t168 + t31;
t716 = t200 * mrSges(5,1) + t30 * mrSges(6,1) + t21 * mrSges(7,1) - t31 * mrSges(6,2) - t26 * mrSges(7,2) - t85 * mrSges(5,3) - t717;
t715 = pkin(10) * t469 - t656;
t331 = t378 * t385 - t380 * t389;
t258 = t331 * t468;
t319 = t331 * qJD(3);
t714 = t653 + (-t258 + t319) * pkin(10) + t718 * pkin(4);
t713 = t682 * t168;
t712 = t682 * t169;
t580 = t342 / 0.2e1;
t591 = t413 / 0.2e1;
t593 = t448 / 0.2e1;
t596 = t199 / 0.2e1;
t605 = t169 / 0.2e1;
t607 = t168 / 0.2e1;
t678 = Ifges(6,3) + Ifges(7,3);
t670 = t168 * t679 + t169 * t681 + t199 * t678;
t711 = -Ifges(5,4) * t591 - Ifges(5,2) * t593 - Ifges(5,6) * t580 + t596 * t678 + t605 * t681 + t607 * t679 + t670 / 0.2e1 + t716;
t372 = pkin(5) * t388 + pkin(4);
t709 = -m(6) * pkin(4) - m(7) * t372 - mrSges(5,1);
t382 = -qJ(6) - pkin(10);
t398 = -m(6) * pkin(10) + m(7) * t382 + mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t708 = Ifges(4,3) + Ifges(5,3);
t669 = t168 * t680 + t199 * t679 + t712;
t668 = t169 * t683 + t199 * t681 + t713;
t658 = t719 * t378 + t380 * t720;
t228 = t258 * t384 + t388 * t469;
t492 = qJD(5) * t388;
t465 = t332 * t492;
t528 = t319 * t384;
t405 = t465 - t528;
t655 = t228 + t405;
t493 = qJD(5) * t384;
t530 = t448 * t384;
t707 = t493 - t530;
t141 = Ifges(5,1) * t413 + Ifges(5,4) * t448 + t342 * Ifges(5,5);
t154 = t378 * t159;
t84 = t147 * t380 - t154;
t636 = -t200 * mrSges(5,2) + t84 * mrSges(5,3);
t706 = -t636 + Ifges(5,1) * t591 + Ifges(5,4) * t593 + Ifges(5,5) * t580 + t141 / 0.2e1;
t491 = qJD(1) * qJD(2);
t315 = (qJDD(1) * t386 + t390 * t491) * t379;
t489 = qJDD(1) * t381;
t360 = qJDD(2) + t489;
t188 = qJD(3) * t282 + t315 * t389 + t360 * t385;
t189 = -qJD(3) * t283 - t315 * t385 + t360 * t389;
t136 = -t188 * t378 + t189 * t380;
t135 = qJDD(5) - t136;
t314 = (-qJDD(1) * t390 + t386 * t491) * t379;
t301 = qJDD(3) + t314;
t490 = qJDD(1) * t379;
t696 = pkin(8) * t490 + qJD(2) * t487;
t697 = -pkin(8) * t379 * t491 + pkin(1) * t489;
t235 = t386 * t697 + t390 * t696;
t209 = pkin(9) * t360 + t235;
t214 = -pkin(1) * t490 + pkin(2) * t314 - pkin(9) * t315;
t94 = -qJD(3) * t187 - t209 * t385 + t389 * t214;
t55 = pkin(3) * t301 - qJ(4) * t188 - qJD(4) * t283 + t94;
t494 = qJD(3) * t389;
t93 = t389 * t209 + t385 * t214 - t265 * t495 + t269 * t494;
t59 = qJ(4) * t189 + qJD(4) * t282 + t93;
t20 = t378 * t55 + t380 * t59;
t18 = pkin(10) * t301 + t20;
t137 = t188 * t380 + t189 * t378;
t236 = -t386 * t696 + t390 * t697;
t210 = -t360 * pkin(2) - t236;
t143 = -t189 * pkin(3) + qJDD(4) + t210;
t34 = -t136 * pkin(4) - t137 * pkin(10) + t143;
t4 = -qJD(5) * t31 - t18 * t384 + t388 * t34;
t70 = qJD(5) * t168 + t137 * t388 + t301 * t384;
t1 = pkin(5) * t135 - qJ(6) * t70 - qJD(6) * t169 + t4;
t584 = t301 / 0.2e1;
t613 = t137 / 0.2e1;
t614 = t136 / 0.2e1;
t615 = t135 / 0.2e1;
t71 = -qJD(5) * t169 - t137 * t384 + t301 * t388;
t622 = t71 / 0.2e1;
t623 = t70 / 0.2e1;
t3 = t388 * t18 + t384 * t34 + t96 * t492 - t493 * t80;
t2 = qJ(6) * t71 + qJD(6) * t168 + t3;
t632 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t677 = t135 * t678 + t679 * t71 + t681 * t70;
t704 = t677 / 0.2e1 + t143 * mrSges(5,1) + t1 * mrSges(7,1) - t20 * mrSges(5,3) + t615 * t678 + t622 * t679 + t623 * t681 + (-t584 - t301 / 0.2e1) * Ifges(5,6) + (-t614 - t136 / 0.2e1) * Ifges(5,2) + (-t613 - t137 / 0.2e1) * Ifges(5,4) + t632;
t626 = m(7) * pkin(5);
t675 = t681 * t135 + t682 * t71 + t683 * t70;
t702 = t675 / 0.2e1;
t676 = t135 * t679 + t680 * t71 + t682 * t70;
t701 = t676 / 0.2e1;
t516 = t379 * t386;
t320 = t381 * t389 - t385 * t516;
t410 = t320 * pkin(3);
t699 = t384 * t715 + t388 * t714;
t373 = pkin(3) * t389 + pkin(2);
t238 = pkin(4) * t331 - pkin(10) * t332 - t373;
t698 = t238 * t492 + t384 * t714 - t388 * t715;
t657 = pkin(4) * t469 - t658;
t576 = cos(qJ(1));
t470 = t576 * t390;
t387 = sin(qJ(1));
t509 = t386 * t387;
t325 = -t381 * t509 + t470;
t514 = t379 * t389;
t254 = -t325 * t385 + t387 * t514;
t695 = t682 * t388;
t694 = t682 * t384;
t597 = -t199 / 0.2e1;
t606 = -t169 / 0.2e1;
t608 = -t168 / 0.2e1;
t693 = t597 * t678 + t606 * t681 + t608 * t679 - t716 + t717;
t377 = qJ(3) + pkin(11);
t374 = sin(t377);
t375 = cos(t377);
t432 = -mrSges(4,1) * t389 + mrSges(4,2) * t385;
t692 = -m(4) * pkin(2) + t374 * t398 + t375 * t709 + t432;
t471 = t576 * t386;
t508 = t387 * t390;
t323 = t381 * t471 + t508;
t472 = t379 * t576;
t244 = t323 * t375 - t374 * t472;
t322 = -t381 * t470 + t509;
t691 = t244 * t384 - t322 * t388;
t690 = -t244 * t388 - t322 * t384;
t19 = -t378 * t59 + t380 * t55;
t687 = t143 * mrSges(5,2) - t19 * mrSges(5,3) + 0.2e1 * Ifges(5,1) * t613 + 0.2e1 * Ifges(5,4) * t614 + 0.2e1 * Ifges(5,5) * t584;
t601 = t188 / 0.2e1;
t600 = t189 / 0.2e1;
t559 = -mrSges(7,1) - mrSges(6,1);
t684 = mrSges(6,2) + mrSges(7,2);
t229 = -t258 * t388 + t384 * t469;
t346 = t383 * t385;
t347 = t383 * t389;
t263 = t346 * t378 - t347 * t380;
t250 = t388 * t263;
t411 = qJ(6) * t319 - qJD(6) * t332;
t674 = qJ(6) * t229 + t411 * t388 + (-t250 + (qJ(6) * t332 - t238) * t384) * qJD(5) + t699 + t718 * pkin(5);
t673 = (-qJD(5) * t263 + t411) * t384 + t698 + (-t228 - t465) * qJ(6);
t167 = t384 * t238 + t250;
t672 = -qJD(5) * t167 + t699;
t671 = -t263 * t493 + t698;
t667 = t235 * mrSges(3,2);
t666 = mrSges(7,1) + t626;
t665 = -m(4) * pkin(9) - mrSges(4,3) - mrSges(5,3);
t664 = pkin(5) * t655 + t657;
t574 = pkin(3) * t283;
t120 = pkin(4) * t413 - pkin(10) * t448 + t574;
t89 = t158 * t380 - t154;
t42 = t384 * t120 + t388 * t89;
t573 = pkin(3) * t378;
t370 = pkin(10) + t573;
t505 = qJ(6) + t370;
t446 = qJD(5) * t505;
t663 = qJ(6) * t530 + qJD(6) * t388 - t384 * t446 - t42;
t41 = t388 * t120 - t384 * t89;
t529 = t448 * t388;
t662 = -pkin(5) * t413 + qJ(6) * t529 - qJD(6) * t384 - t388 * t446 - t41;
t88 = t158 * t378 + t512;
t661 = pkin(5) * t707 - t88;
t110 = -mrSges(6,1) * t168 + mrSges(6,2) * t169;
t176 = mrSges(5,1) * t342 - mrSges(5,3) * t413;
t660 = t176 - t110;
t445 = mrSges(3,3) * t469;
t659 = -mrSges(3,1) * t361 - mrSges(4,1) * t282 + mrSges(4,2) * t283 + t445;
t527 = t319 * t388;
t404 = t332 * t493 + t527;
t654 = t229 + t404;
t513 = t379 * t390;
t328 = t381 * t386 * pkin(1) + pkin(8) * t513;
t296 = pkin(9) * t381 + t328;
t213 = t389 * t296 + t385 * t297;
t362 = pkin(8) * t516;
t575 = pkin(1) * t390;
t327 = t381 * t575 - t362;
t427 = mrSges(7,1) * t384 + mrSges(7,2) * t388;
t429 = mrSges(6,1) * t384 + mrSges(6,2) * t388;
t79 = -pkin(4) * t342 - t84;
t54 = -pkin(5) * t168 + qJD(6) + t79;
t651 = t54 * t427 + t79 * t429;
t650 = -t384 * t679 + t388 * t681;
t649 = -t384 * t680 + t695;
t648 = t388 * t683 - t694;
t396 = t323 * t385 + t389 * t472;
t647 = Ifges(4,5) * t188 + Ifges(5,5) * t137 + Ifges(4,6) * t189 + Ifges(5,6) * t136 + t301 * t708;
t646 = -t492 + t529;
t644 = -t385 * t94 + t389 * t93;
t643 = t3 * t388 - t384 * t4;
t642 = m(6) + m(5) + m(7);
t641 = -mrSges(6,1) - t666;
t640 = mrSges(3,2) + t665;
t428 = -mrSges(7,1) * t388 + mrSges(7,2) * t384;
t430 = -mrSges(6,1) * t388 + mrSges(6,2) * t384;
t639 = -t428 - t430 - t709;
t635 = -t94 * mrSges(4,1) - t19 * mrSges(5,1) + t93 * mrSges(4,2) + t20 * mrSges(5,2);
t634 = -t384 * t626 + t665;
t633 = mrSges(3,1) - t692;
t515 = t379 * t387;
t248 = t325 * t375 + t374 * t515;
t294 = t374 * t381 + t375 * t516;
t631 = -g(1) * t248 - g(2) * t244 - g(3) * t294;
t628 = t379 ^ 2;
t627 = m(5) * pkin(3);
t621 = pkin(1) * mrSges(3,1);
t620 = pkin(1) * mrSges(3,2);
t619 = Ifges(4,4) * t601 + Ifges(4,2) * t600 + Ifges(4,6) * t584;
t618 = Ifges(4,1) * t601 + Ifges(4,4) * t600 + Ifges(4,5) * t584;
t610 = -t141 / 0.2e1;
t537 = t283 * Ifges(4,4);
t180 = t282 * Ifges(4,2) + t342 * Ifges(4,6) + t537;
t603 = t180 / 0.2e1;
t272 = Ifges(4,4) * t282;
t181 = t283 * Ifges(4,1) + t342 * Ifges(4,5) + t272;
t602 = t181 / 0.2e1;
t585 = t283 / 0.2e1;
t579 = t381 / 0.2e1;
t572 = pkin(3) * t380;
t571 = pkin(5) * t169;
t557 = mrSges(6,3) * t168;
t556 = mrSges(6,3) * t169;
t555 = mrSges(7,3) * t168;
t554 = mrSges(7,3) * t169;
t553 = Ifges(3,4) * t386;
t552 = Ifges(3,4) * t390;
t551 = Ifges(4,4) * t385;
t550 = Ifges(4,4) * t389;
t545 = Ifges(3,6) * t361;
t542 = t186 * mrSges(4,3);
t541 = t187 * mrSges(4,3);
t540 = t448 * Ifges(5,6);
t539 = t413 * Ifges(5,5);
t538 = t282 * Ifges(4,6);
t536 = t283 * Ifges(4,5);
t535 = t361 * Ifges(3,5);
t310 = qJD(2) * t409;
t312 = t327 * qJD(2);
t149 = -qJD(3) * t213 + t389 * t310 - t312 * t385;
t496 = qJD(2) * t379;
t466 = t390 * t496;
t252 = qJD(3) * t320 + t389 * t466;
t321 = t381 * t385 + t386 * t514;
t467 = t386 * t496;
t101 = pkin(3) * t467 - qJ(4) * t252 - qJD(4) * t321 + t149;
t148 = -t296 * t495 + t297 * t494 + t385 * t310 + t389 * t312;
t251 = -qJD(3) * t321 - t385 * t466;
t108 = qJ(4) * t251 + qJD(4) * t320 + t148;
t47 = t378 * t101 + t380 * t108;
t212 = -t385 * t296 + t389 * t297;
t163 = -pkin(3) * t513 - t321 * qJ(4) + t212;
t177 = qJ(4) * t320 + t213;
t105 = t378 * t163 + t380 * t177;
t100 = -pkin(10) * t513 + t105;
t231 = -t380 * t320 + t321 * t378;
t232 = t320 * t378 + t321 * t380;
t295 = t362 + (-pkin(2) - t575) * t381;
t237 = t295 - t410;
t138 = t231 * pkin(4) - t232 * pkin(10) + t237;
t49 = t388 * t100 + t384 * t138;
t524 = t323 * t384;
t522 = t325 * t384;
t520 = t332 * t384;
t519 = t332 * t388;
t518 = t375 * t384;
t517 = t375 * t388;
t510 = t384 * t390;
t507 = t388 * t390;
t500 = t576 * pkin(1) + pkin(8) * t515;
t481 = t379 * t510;
t479 = t385 * t515;
t477 = t379 * t507;
t473 = Ifges(3,5) * t315 - Ifges(3,6) * t314 + Ifges(3,3) * t360;
t23 = -t71 * mrSges(7,1) + t70 * mrSges(7,2);
t453 = -t493 / 0.2e1;
t451 = -pkin(1) * t387 + pkin(8) * t472;
t58 = -t136 * mrSges(5,1) + t137 * mrSges(5,2);
t48 = -t100 * t384 + t388 * t138;
t46 = t101 * t380 - t378 * t108;
t104 = t163 * t380 - t378 * t177;
t166 = t388 * t238 - t263 * t384;
t350 = t385 * t472;
t447 = -t323 * t389 + t350;
t262 = -t380 * t346 - t347 * t378;
t444 = mrSges(3,3) * t468;
t435 = t254 * pkin(3);
t99 = pkin(4) * t513 - t104;
t433 = mrSges(4,1) * t320 - mrSges(4,2) * t321;
t426 = Ifges(4,1) * t389 - t551;
t423 = Ifges(3,2) * t390 + t553;
t422 = -Ifges(4,2) * t385 + t550;
t419 = Ifges(4,5) * t389 - Ifges(4,6) * t385;
t324 = t381 * t508 + t471;
t414 = pkin(3) * t479 - t324 * t383 + t325 * t373 + t500;
t194 = -t248 * t384 + t324 * t388;
t17 = -pkin(4) * t301 - t19;
t197 = -t384 * t232 - t477;
t408 = -t388 * t232 + t481;
t403 = pkin(3) * t350 + t322 * t383 - t323 * t373 + t451;
t44 = pkin(10) * t467 + t47;
t173 = -t380 * t251 + t252 * t378;
t174 = t251 * t378 + t252 * t380;
t313 = t328 * qJD(2);
t208 = -t251 * pkin(3) + t313;
t83 = t173 * pkin(4) - t174 * pkin(10) + t208;
t7 = -t100 * t493 + t138 * t492 + t384 * t83 + t388 * t44;
t243 = t323 * t374 + t375 * t472;
t395 = -mrSges(3,2) - t634;
t394 = -g(1) * t324 - g(2) * t322 + g(3) * t513;
t43 = -pkin(4) * t467 - t46;
t393 = t396 * pkin(3);
t8 = -qJD(5) * t49 - t384 * t44 + t388 * t83;
t371 = -pkin(4) - t572;
t355 = Ifges(3,4) * t468;
t338 = -t372 - t572;
t330 = t505 * t388;
t329 = t505 * t384;
t326 = (-mrSges(3,1) * t390 + mrSges(3,2) * t386) * t379;
t307 = -t361 * mrSges(3,2) + t444;
t293 = t374 * t516 - t381 * t375;
t260 = Ifges(3,1) * t469 + t355 + t535;
t259 = t423 * t499 + t545;
t255 = t325 * t389 + t479;
t247 = t325 * t374 - t375 * t515;
t240 = mrSges(4,1) * t342 - mrSges(4,3) * t283;
t239 = -mrSges(4,2) * t342 + mrSges(4,3) * t282;
t218 = pkin(5) * t520 + t262;
t195 = t248 * t388 + t324 * t384;
t179 = t342 * Ifges(4,3) + t536 + t538;
t175 = -mrSges(5,2) * t342 + mrSges(5,3) * t448;
t161 = -mrSges(4,2) * t301 + mrSges(4,3) * t189;
t160 = mrSges(4,1) * t301 - mrSges(4,3) * t188;
t152 = -qJ(6) * t520 + t167;
t145 = pkin(5) * t331 - qJ(6) * t519 + t166;
t144 = -mrSges(5,1) * t448 + mrSges(5,2) * t413;
t142 = -mrSges(4,1) * t189 + mrSges(4,2) * t188;
t139 = t342 * Ifges(5,3) + t539 + t540;
t126 = mrSges(6,1) * t199 - t556;
t125 = mrSges(7,1) * t199 - t554;
t124 = -mrSges(6,2) * t199 + t557;
t123 = -mrSges(7,2) * t199 + t555;
t115 = qJD(5) * t408 - t384 * t174 + t388 * t467;
t114 = qJD(5) * t197 + t388 * t174 + t384 * t467;
t109 = -mrSges(7,1) * t168 + mrSges(7,2) * t169;
t107 = mrSges(5,1) * t301 - mrSges(5,3) * t137;
t106 = -mrSges(5,2) * t301 + mrSges(5,3) * t136;
t72 = -pkin(5) * t197 + t99;
t39 = qJ(6) * t197 + t49;
t38 = -mrSges(6,2) * t135 + mrSges(6,3) * t71;
t37 = -mrSges(7,2) * t135 + mrSges(7,3) * t71;
t36 = mrSges(6,1) * t135 - mrSges(6,3) * t70;
t35 = mrSges(7,1) * t135 - mrSges(7,3) * t70;
t28 = pkin(5) * t231 + qJ(6) * t408 + t48;
t24 = -mrSges(6,1) * t71 + mrSges(6,2) * t70;
t22 = -pkin(5) * t115 + t43;
t9 = -pkin(5) * t71 + qJDD(6) + t17;
t6 = qJ(6) * t115 + qJD(6) * t197 + t7;
t5 = pkin(5) * t173 - qJ(6) * t114 + qJD(6) * t408 + t8;
t10 = [-t381 * t667 + t659 * t313 + t252 * t602 + t251 * t603 + (-m(3) * t451 + t323 * mrSges(3,1) - mrSges(3,3) * t472 + t387 * mrSges(2,1) + t576 * mrSges(2,2) - m(4) * (-pkin(2) * t323 + t451) - t447 * mrSges(4,1) - t396 * mrSges(4,2) - m(7) * (-t244 * t372 + t403) - m(5) * t403 + t244 * mrSges(5,1) - m(6) * (-pkin(4) * t244 + t403) + t559 * t690 - t684 * t691 + t395 * t322 - t398 * t243) * g(1) + (-t576 * mrSges(2,1) - m(3) * t500 - t325 * mrSges(3,1) - m(7) * (t248 * t372 + t414) - m(5) * t414 - t248 * mrSges(5,1) - m(6) * (pkin(4) * t248 + t414) - m(4) * (pkin(2) * t325 + t500) - t255 * mrSges(4,1) - t254 * mrSges(4,2) + (-mrSges(3,3) * t379 + mrSges(2,2)) * t387 + t559 * t195 - t684 * t194 - t395 * t324 + t398 * t247) * g(2) + (Ifges(4,1) * t321 + Ifges(4,4) * t320) * t601 - t210 * t433 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t628 + t235 * t328 + t236 * t327 - t308 * t313 + t311 * t312) + (-t327 * mrSges(3,3) + Ifges(3,5) * t579 + (t386 * Ifges(3,1) + t552 - t620) * t379) * t315 + ((-t308 * mrSges(3,3) + t260 / 0.2e1 + t535 / 0.2e1 + (-t620 + t552 / 0.2e1) * t499) * t390 + (-t311 * mrSges(3,3) + t539 / 0.2e1 + t540 / 0.2e1 + t84 * mrSges(5,1) - t85 * mrSges(5,2) + t536 / 0.2e1 + t538 / 0.2e1 - t187 * mrSges(4,2) + t186 * mrSges(4,1) - t259 / 0.2e1 + t179 / 0.2e1 + t139 / 0.2e1 - t545 / 0.2e1 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t342 + (-t621 - t553 / 0.2e1 + (Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1) * t390) * t499) * t386) * t496 - (t328 * mrSges(3,3) + Ifges(3,6) * t579 + (t423 + t621) * t379) * t314 + (t327 * mrSges(3,1) - t328 * mrSges(3,2) + Ifges(3,3) * t579 + (Ifges(3,5) * t386 + Ifges(3,6) * t390) * t379) * t360 + t704 * t231 + (Ifges(4,4) * t321 + Ifges(4,2) * t320) * t600 + m(4) * (t148 * t187 + t149 * t186 + t210 * t295 + t212 * t94 + t213 * t93 + t264 * t313) + m(5) * (t104 * t19 + t105 * t20 + t143 * t237 + t200 * t208 + t46 * t84 + t47 * t85) + m(7) * (t1 * t28 + t2 * t39 + t21 * t5 + t22 * t54 + t26 * t6 + t72 * t9) + m(6) * (t17 * t99 + t3 * t49 + t30 * t8 + t31 * t7 + t4 * t48 + t43 * t79) + (-t647 / 0.2e1 + mrSges(3,3) * t235 - Ifges(4,5) * t601 - Ifges(5,5) * t613 - Ifges(4,6) * t600 - Ifges(5,6) * t614 - t708 * t584 + t635) * t513 + t706 * t174 + t282 * (Ifges(4,4) * t252 + Ifges(4,2) * t251) / 0.2e1 + (-t21 * mrSges(7,3) - t30 * mrSges(6,3) + t668 / 0.2e1 + t681 * t596 + t682 * t607 + t683 * t605 + t54 * mrSges(7,2) + t79 * mrSges(6,2)) * t114 + t236 * (mrSges(3,1) * t381 - mrSges(3,3) * t516) + t473 * t579 + (Ifges(4,5) * t252 + Ifges(4,6) * t251) * t580 + (Ifges(4,1) * t252 + Ifges(4,4) * t251) * t585 + (t26 * mrSges(7,3) + t31 * mrSges(6,3) + t669 / 0.2e1 + t679 * t596 + t680 * t607 + t682 * t605 - t54 * mrSges(7,1) - t79 * mrSges(6,1)) * t115 + (-t186 * t252 + t187 * t251 + t320 * t93 - t321 * t94) * mrSges(4,3) + (-pkin(1) * t326 * t379 + Ifges(2,3)) * qJDD(1) + t711 * t173 + (-t17 * mrSges(6,1) - t9 * mrSges(7,1) + t3 * mrSges(6,3) + t2 * mrSges(7,3) + t679 * t615 + t680 * t622 + t682 * t623 + t701) * t197 + t687 * t232 + (t4 * mrSges(6,3) + t1 * mrSges(7,3) - t17 * mrSges(6,2) - t9 * mrSges(7,2) - t675 / 0.2e1 - t681 * t615 - t682 * t622 - t683 * t623) * t408 + (Ifges(4,5) * t321 + Ifges(4,6) * t320) * t584 + t312 * t307 + t295 * t142 + t264 * (-mrSges(4,1) * t251 + mrSges(4,2) * t252) + t237 * t58 + t148 * t239 + t149 * t240 + t208 * t144 + t212 * t160 + t213 * t161 + t47 * t175 + t46 * t176 + t6 * t123 + t7 * t124 + t5 * t125 + t8 * t126 + t321 * t618 + t320 * t619 + t28 * t35 + t39 * t37 + t48 * t36 + t49 * t38 + t72 * t23 + t99 * t24 + t105 * t106 + t104 * t107 + t22 * t109 + t43 * t110; (t693 - t670 / 0.2e1) * t257 + (-t541 - t180 / 0.2e1) * t495 + (t602 - t542) * t494 + t664 * t109 + (m(4) * ((-t186 * t389 - t187 * t385) * qJD(3) + t644) - t240 * t494 - t239 * t495 - t385 * t160 + t389 * t161) * pkin(9) + t644 * mrSges(4,3) + (t200 * t258 + t469 * t85) * mrSges(5,2) - t84 * (mrSges(5,1) * t469 + mrSges(5,3) * t258) + (-Ifges(5,1) * t258 + Ifges(5,5) * t469) * t592 + (-Ifges(5,4) * t258 + Ifges(5,6) * t469) * t594 + (-Ifges(5,5) * t258 + Ifges(5,3) * t469) * t581 + (t282 * t422 + t283 * t426 + t342 * t419) * qJD(3) / 0.2e1 + (t24 - t107) * t262 + t653 * t144 + (mrSges(7,1) * t655 - mrSges(7,2) * t654) * t54 + (mrSges(6,1) * t655 - mrSges(6,2) * t654) * t79 + (-t3 * t520 + t30 * t654 - t31 * t655 - t4 * t519) * mrSges(6,3) + (-t1 * t519 - t2 * t520 + t21 * t654 - t26 * t655) * mrSges(7,3) + t656 * t175 + t657 * t110 + (-t143 * t373 - t19 * t262 + t20 * t263 + t200 * t653 + t656 * t85 + t658 * t84) * m(5) + t658 * t176 + (-m(4) * t264 + t445 - t659) * t311 + t519 * t702 + (Ifges(4,2) * t389 + t551) * t600 + (Ifges(4,1) * t385 + t550) * t601 + t440 * t603 - ((-Ifges(3,2) * t469 + t389 * t181 + t260 + t355) * t390 + (t179 + t139) * t386 + t342 * (Ifges(4,3) * t386 + t390 * t419) + t283 * (Ifges(4,5) * t386 + t390 * t426) + t282 * (Ifges(4,6) * t386 + t390 * t422) + t361 * (Ifges(3,5) * t390 - Ifges(3,6) * t386)) * t499 / 0.2e1 + t259 * t469 / 0.2e1 + t210 * t432 + (-t386 * (Ifges(3,1) * t390 - t553) / 0.2e1 + pkin(1) * (mrSges(3,1) * t386 + mrSges(3,2) * t390)) * qJD(1) ^ 2 * t628 + (-t186 * (mrSges(4,1) * t386 - mrSges(4,3) * t506) - t187 * (-mrSges(4,3) * t385 * t390 - mrSges(4,2) * t386)) * t499 + t704 * t331 - t258 * t610 + t668 * (t332 * t453 - t527 / 0.2e1 - t229 / 0.2e1) + t669 * (t528 / 0.2e1 - t228 / 0.2e1 - t465 / 0.2e1) - t706 * t319 + t473 + (-pkin(2) * t210 - t186 * t223 - t187 * t224) * m(4) + (Ifges(4,5) * t385 + Ifges(4,6) * t389) * t584 - t667 + t711 * t318 + (t444 - t307) * t308 + (-t642 * t373 * t513 + t326 + (t559 * (t375 * t507 + t384 * t386) - t684 * (-t375 * t510 + t386 * t388) + t692 * t390 + (t383 * t642 + t634) * t386) * t379) * g(3) + (t17 * t429 + t427 * t9 + t615 * t650 + t622 * t649 + t623 * t648 + t687) * t332 - t373 * t58 + t263 * t106 - t224 * t239 - t223 * t240 + t236 * mrSges(3,1) + t218 * t23 + t167 * t38 + t166 * t36 + t152 * t37 + t145 * t35 - pkin(2) * t142 + t385 * t618 + t389 * t619 + t671 * t124 + t672 * t126 + (t166 * t4 + t167 * t3 + t17 * t262 + t30 * t672 + t31 * t671 + t657 * t79) * m(6) + t673 * t123 + t674 * t125 + (t1 * t145 + t152 * t2 + t21 * t674 + t218 * t9 + t26 * t673 + t54 * t664) * m(7) - t676 * t520 / 0.2e1 + (t228 * t679 + t229 * t681) * t597 + (-t404 * t681 - t405 * t679) * t596 + (t228 * t680 + t229 * t682) * t608 + (-t404 * t682 - t405 * t680) * t607 + (t228 * t682 + t229 * t683) * t606 + (-t404 * t683 - t405 * t682) * t605 + (-t524 * t626 - t642 * (-t322 * t373 - t323 * t383) + t559 * (-t322 * t517 + t524) - t684 * (t322 * t518 + t323 * t388) + t640 * t323 + t633 * t322) * g(2) + (-t522 * t626 - t642 * (-t324 * t373 - t325 * t383) + t559 * (-t324 * t517 + t522) - t684 * (t324 * t518 + t325 * t388) + t640 * t325 + t633 * t324) * g(1) + t342 * t264 * (mrSges(4,1) * t385 + mrSges(4,2) * t389); t647 + (-t1 * t329 + t2 * t330 + t21 * t662 + t26 * t663 + t338 * t9 + t54 * t661) * m(7) + t663 * t123 + t661 * t109 + t662 * t125 + (m(6) * ((-t30 * t388 - t31 * t384) * qJD(5) + t643) - t126 * t492 - t124 * t493 - t384 * t36 + t388 * t38) * t370 + (-m(6) * (pkin(10) * t294 + t410) - t433 - m(7) * (-t294 * t382 + t410) - m(5) * t410 + t294 * mrSges(5,2) + t639 * t293) * g(3) + (-m(6) * (pkin(10) * t248 + t435) - m(7) * (-t248 * t382 + t435) - m(5) * t435 + t248 * mrSges(5,2) - mrSges(4,1) * t254 + mrSges(4,2) * t255 + t639 * t247) * g(1) - (-Ifges(4,2) * t283 + t181 + t272) * t282 / 0.2e1 + (Ifges(5,1) * t592 + Ifges(5,4) * t594 + Ifges(5,5) * t581 + t597 * t650 + t606 * t648 + t608 * t649 + t610 + t636 - t651) * t448 + t651 * qJD(5) + (m(5) * t84 + t660) * t88 + (t30 * t646 - t31 * t707 + t631 + t643) * mrSges(6,3) + (-t1 * t384 + t2 * t388 + t21 * t646 - t26 * t707 + t631) * mrSges(7,3) + t388 * t701 + t384 * t702 + t670 * t592 - t283 * (Ifges(4,1) * t282 - t537) / 0.2e1 + t107 * t572 + t106 * t573 + t17 * t430 - t635 + (t530 / 0.2e1 + t453) * t669 - t144 * t574 + t9 * t428 + t283 * t541 + t282 * t542 + (t17 * t371 - t30 * t41 - t31 * t42 - t79 * t88) * m(6) + t693 * t413 + (t388 * t680 + t694) * t622 + (t384 * t683 + t695) * t623 + (Ifges(4,5) * t282 - Ifges(4,6) * t283) * t581 + t180 * t585 + (t492 / 0.2e1 - t529 / 0.2e1) * t668 + (-m(6) * (t244 * pkin(10) - t393) - m(7) * (-t244 * t382 - t393) + t244 * mrSges(5,2) - mrSges(4,2) * t447 + (t627 + mrSges(4,1)) * t396 + t639 * t243) * g(2) + (t168 * t649 + t169 * t648 + t199 * t650) * qJD(5) / 0.2e1 + t371 * t24 + t338 * t23 - t329 * t35 + t330 * t37 - t264 * (mrSges(4,1) * t283 + mrSges(4,2) * t282) - t186 * t239 + t187 * t240 - t89 * t175 - t42 * t124 - t41 * t126 + (t384 * t681 + t388 * t679) * t615 + (t19 * t380 + t20 * t378) * t627 - m(5) * (t200 * t574 + t85 * t89); -t448 * t175 + (-t109 + t660) * t413 + (t35 + t36 + t199 * (t123 + t124)) * t388 + (t37 + t38 - t199 * (t125 + t126)) * t384 + t58 + (t1 * t388 + t2 * t384 - t413 * t54 + t394 + t199 * (-t21 * t384 + t26 * t388)) * m(7) + (-t413 * t79 + t3 * t384 + t4 * t388 + t394 + t199 * (-t30 * t384 + t31 * t388)) * m(6) + (t413 * t84 - t448 * t85 + t143 + t394) * m(5); (-t641 * t691 - t684 * t690) * g(2) + (-t684 * (-t294 * t388 + t481) + t641 * (-t294 * t384 - t477)) * g(3) + t632 + t21 * t555 - t109 * t571 + (-t169 * t680 + t668 + t713) * t608 + (t194 * t641 + t195 * t684) * g(1) + (t168 * t683 - t712) * t606 + (t168 * t681 - t169 * t679) * t597 + t666 * t1 + t669 * t605 + (-m(7) * t571 - mrSges(7,1) * t169 - mrSges(7,2) * t168) * t54 + (t556 + t126) * t31 + (t557 - t124) * t30 + (t554 - m(7) * (-t21 + t25) + t125) * t26 - t79 * (mrSges(6,1) * t169 + mrSges(6,2) * t168) + pkin(5) * t35 - t25 * t123 + t677; -t168 * t123 + t169 * t125 + (-g(1) * t247 - g(2) * t243 - g(3) * t293 - t168 * t26 + t21 * t169 + t9) * m(7) + t23;];
tau  = t10;
