% Calculate vector of inverse dynamics joint torques for
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:38:22
% EndTime: 2019-03-09 15:40:01
% DurationCPUTime: 64.32s
% Computational Cost: add. (28367->1037), mult. (69026->1418), div. (0->0), fcn. (56364->18), ass. (0->446)
t422 = cos(qJ(2));
t418 = sin(qJ(2));
t410 = sin(pkin(6));
t505 = qJD(1) * t410;
t479 = t418 * t505;
t413 = cos(pkin(6));
t504 = qJD(1) * t413;
t494 = pkin(1) * t504;
t330 = -pkin(8) * t479 + t422 * t494;
t440 = t410 * (pkin(2) * t418 - pkin(9) * t422);
t331 = qJD(1) * t440;
t417 = sin(qJ(3));
t421 = cos(qJ(3));
t249 = -t417 * t330 + t421 * t331;
t414 = -qJ(4) - pkin(9);
t470 = qJD(3) * t414;
t509 = t421 * t422;
t709 = -(pkin(3) * t418 - qJ(4) * t509) * t505 - t249 - qJD(4) * t417 + t421 * t470;
t250 = t421 * t330 + t417 * t331;
t478 = t422 * t505;
t462 = t417 * t478;
t708 = -qJ(4) * t462 - qJD(4) * t421 - t417 * t470 + t250;
t406 = pkin(12) + qJ(6);
t401 = sin(t406);
t403 = cos(t406);
t411 = cos(pkin(12));
t398 = pkin(5) * t411 + pkin(4);
t408 = sin(pkin(12));
t455 = -mrSges(6,1) * t411 + mrSges(6,2) * t408;
t425 = -m(6) * pkin(4) - m(7) * t398 - mrSges(5,1) + t455;
t707 = -mrSges(7,1) * t403 + mrSges(7,2) * t401 + t425;
t409 = sin(pkin(11));
t412 = cos(pkin(11));
t358 = t409 * t421 + t412 * t417;
t283 = t358 * t478;
t340 = t358 * qJD(3);
t706 = -t283 + t340;
t369 = qJD(3) - t478;
t568 = -t369 / 0.2e1;
t388 = qJD(2) + t504;
t304 = t388 * t421 - t417 * t479;
t305 = t388 * t417 + t421 * t479;
t443 = t304 * t409 + t412 * t305;
t577 = -t443 / 0.2e1;
t234 = -t412 * t304 + t305 * t409;
t578 = t234 / 0.2e1;
t705 = Ifges(5,4) * t577 + Ifges(5,2) * t578 + Ifges(5,6) * t568;
t630 = mrSges(5,2) - m(6) * qJ(5) - mrSges(6,3) + m(7) * (-pkin(10) - qJ(5));
t652 = t409 * t709 - t708 * t412;
t333 = pkin(8) * t478 + t418 * t494;
t500 = qJD(3) * t417;
t651 = -t333 + (-t462 + t500) * pkin(3);
t430 = -mrSges(7,3) + t630;
t704 = qJ(5) * t479 - t652;
t356 = t409 * t417 - t412 * t421;
t284 = t356 * t478;
t341 = t356 * qJD(3);
t703 = -qJD(5) * t358 + t651 + (-t284 + t341) * qJ(5) + t706 * pkin(4);
t291 = pkin(9) * t388 + t333;
t319 = (-pkin(2) * t422 - pkin(9) * t418 - pkin(1)) * t410;
t295 = qJD(1) * t319;
t221 = -t291 * t417 + t421 * t295;
t185 = -qJ(4) * t305 + t221;
t171 = pkin(3) * t369 + t185;
t222 = t291 * t421 + t295 * t417;
t186 = qJ(4) * t304 + t222;
t513 = t412 * t186;
t102 = t409 * t171 + t513;
t198 = t369 * t408 + t411 * t443;
t100 = qJ(5) * t369 + t102;
t290 = -t388 * pkin(2) - t330;
t232 = -t304 * pkin(3) + qJD(4) + t290;
t111 = t234 * pkin(4) - qJ(5) * t443 + t232;
t52 = -t100 * t408 + t411 * t111;
t35 = pkin(5) * t234 - pkin(10) * t198 + t52;
t469 = t411 * t369 - t408 * t443;
t53 = t411 * t100 + t408 * t111;
t41 = pkin(10) * t469 + t53;
t416 = sin(qJ(6));
t420 = cos(qJ(6));
t12 = t35 * t420 - t41 * t416;
t13 = t35 * t416 + t41 * t420;
t702 = mrSges(5,1) * t232 + mrSges(6,1) * t52 + mrSges(7,1) * t12 - mrSges(6,2) * t53 - mrSges(7,2) * t13 - mrSges(5,3) * t102 + t705;
t252 = -t284 * t411 + t408 * t479;
t523 = t341 * t411;
t701 = t252 + t523;
t251 = t284 * t408 + t411 * t479;
t644 = -t408 * t341 + t251;
t567 = t369 / 0.2e1;
t576 = t443 / 0.2e1;
t579 = -t234 / 0.2e1;
t231 = qJD(6) + t234;
t581 = t231 / 0.2e1;
t588 = t198 / 0.2e1;
t590 = t469 / 0.2e1;
t125 = t198 * t420 + t416 * t469;
t602 = t125 / 0.2e1;
t671 = -t198 * t416 + t420 * t469;
t604 = t671 / 0.2e1;
t700 = -Ifges(5,4) * t576 + Ifges(6,5) * t588 + Ifges(7,5) * t602 - Ifges(5,2) * t579 - Ifges(5,6) * t567 + Ifges(6,6) * t590 + Ifges(7,6) * t604 + Ifges(6,3) * t578 + Ifges(7,3) * t581 + t702;
t407 = qJ(3) + pkin(11);
t402 = sin(t407);
t404 = cos(t407);
t457 = -mrSges(4,1) * t421 + mrSges(4,2) * t417;
t699 = -m(4) * pkin(2) + t430 * t402 + t404 * t707 + t457;
t498 = qJD(1) * qJD(2);
t337 = (qJDD(1) * t418 + t422 * t498) * t410;
t496 = qJDD(1) * t413;
t387 = qJDD(2) + t496;
t226 = qJD(3) * t304 + t337 * t421 + t387 * t417;
t227 = -qJD(3) * t305 - t337 * t417 + t387 * t421;
t156 = t226 * t412 + t227 * t409;
t336 = (-qJDD(1) * t422 + t418 * t498) * t410;
t323 = qJDD(3) + t336;
t118 = -t156 * t408 + t323 * t411;
t607 = t118 / 0.2e1;
t119 = t156 * t411 + t323 * t408;
t606 = t119 / 0.2e1;
t155 = t226 * t409 - t412 * t227;
t597 = t155 / 0.2e1;
t698 = Ifges(4,3) + Ifges(5,3);
t664 = t704 * t408 + t703 * t411;
t663 = t703 * t408 - t704 * t411;
t497 = qJDD(1) * t410;
t675 = pkin(8) * t497 + qJD(2) * t494;
t676 = -pkin(8) * t410 * t498 + pkin(1) * t496;
t261 = t418 * t676 + t422 * t675;
t241 = pkin(9) * t387 + t261;
t247 = -pkin(1) * t497 + pkin(2) * t336 - pkin(9) * t337;
t110 = -qJD(3) * t222 - t241 * t417 + t421 * t247;
t79 = pkin(3) * t323 - qJ(4) * t226 - qJD(4) * t305 + t110;
t499 = qJD(3) * t421;
t109 = t421 * t241 + t417 * t247 - t291 * t500 + t295 * t499;
t83 = qJ(4) * t227 + qJD(4) * t304 + t109;
t32 = t409 * t79 + t412 * t83;
t29 = qJ(5) * t323 + qJD(5) * t369 + t32;
t262 = -t418 * t675 + t422 * t676;
t242 = -t387 * pkin(2) - t262;
t167 = -t227 * pkin(3) + qJDD(4) + t242;
t48 = t155 * pkin(4) - t156 * qJ(5) - qJD(5) * t443 + t167;
t10 = -t29 * t408 + t411 * t48;
t11 = t411 * t29 + t408 * t48;
t571 = t323 / 0.2e1;
t596 = t156 / 0.2e1;
t598 = -t155 / 0.2e1;
t153 = qJDD(6) + t155;
t599 = t153 / 0.2e1;
t40 = -qJD(6) * t125 + t118 * t420 - t119 * t416;
t622 = t40 / 0.2e1;
t39 = qJD(6) * t671 + t118 * t416 + t119 * t420;
t623 = t39 / 0.2e1;
t5 = pkin(5) * t155 - pkin(10) * t119 + t10;
t6 = pkin(10) * t118 + t11;
t1 = qJD(6) * t12 + t416 * t5 + t420 * t6;
t2 = -qJD(6) * t13 - t416 * t6 + t420 * t5;
t634 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t7 = Ifges(7,5) * t39 + Ifges(7,6) * t40 + Ifges(7,3) * t153;
t694 = mrSges(5,1) * t167 + mrSges(6,1) * t10 + 0.2e1 * Ifges(6,5) * t606 + Ifges(7,5) * t623 + 0.2e1 * Ifges(6,6) * t607 + Ifges(7,6) * t622 + 0.2e1 * Ifges(6,3) * t597 + Ifges(7,3) * t599 + t7 / 0.2e1 - mrSges(6,2) * t11 - mrSges(5,3) * t32 + (-t323 / 0.2e1 - t571) * Ifges(5,6) + (t597 - t598) * Ifges(5,2) + (-t156 / 0.2e1 - t596) * Ifges(5,4) + t634;
t31 = -t409 * t83 + t412 * t79;
t693 = t167 * mrSges(5,2) - mrSges(5,3) * t31 + 0.2e1 * Ifges(5,1) * t596 + 0.2e1 * Ifges(5,4) * t598 + 0.2e1 * Ifges(5,5) * t571;
t691 = pkin(5) * t706 + t701 * pkin(10) + t664;
t690 = -pkin(10) * t644 + t663;
t654 = t708 * t409 + t412 * t709;
t182 = t409 * t186;
t101 = t171 * t412 - t182;
t685 = t232 * mrSges(5,2) - t101 * mrSges(5,3);
t159 = Ifges(5,1) * t443 - t234 * Ifges(5,4) + t369 * Ifges(5,5);
t684 = Ifges(5,1) * t576 + Ifges(5,4) * t579 + Ifges(5,5) * t567 + t159 / 0.2e1;
t681 = m(7) * pkin(5);
t517 = t410 * t418;
t344 = t413 * t421 - t417 * t517;
t679 = t344 * pkin(3);
t678 = t198 * Ifges(6,5) + t125 * Ifges(7,5) + Ifges(6,6) * t469 + Ifges(7,6) * t671 + t234 * Ifges(6,3) + t231 * Ifges(7,3);
t677 = t261 * mrSges(3,2);
t653 = pkin(4) * t479 - t654;
t674 = Ifges(4,5) * t226 + Ifges(5,5) * t156 + Ifges(4,6) * t227 - Ifges(5,6) * t155 + t323 * t698;
t565 = cos(qJ(1));
t480 = t565 * t422;
t419 = sin(qJ(1));
t511 = t418 * t419;
t349 = -t413 * t511 + t480;
t515 = t410 * t421;
t280 = -t349 * t417 + t419 * t515;
t582 = -t231 / 0.2e1;
t589 = -t198 / 0.2e1;
t591 = -t469 / 0.2e1;
t603 = -t125 / 0.2e1;
t605 = -t671 / 0.2e1;
t673 = Ifges(6,5) * t589 + Ifges(7,5) * t603 + Ifges(6,6) * t591 + Ifges(7,6) * t605 + Ifges(6,3) * t579 + Ifges(7,3) * t582 - t702 - t705;
t453 = t401 * mrSges(7,1) + t403 * mrSges(7,2);
t530 = t411 * mrSges(6,2);
t454 = t408 * mrSges(6,1) + t530;
t665 = -m(4) * pkin(9) - mrSges(4,3) - mrSges(5,3);
t672 = -t408 * t681 - t453 - t454 + t665;
t670 = -t110 * mrSges(4,1) - t31 * mrSges(5,1) + t109 * mrSges(4,2) + t32 * mrSges(5,2);
t669 = -m(7) - m(6);
t584 = t226 / 0.2e1;
t583 = t227 / 0.2e1;
t400 = pkin(3) * t421 + pkin(2);
t266 = pkin(4) * t356 - qJ(5) * t358 - t400;
t373 = t414 * t417;
t374 = t414 * t421;
t289 = t373 * t409 - t374 * t412;
t193 = t411 * t266 - t289 * t408;
t519 = t358 * t411;
t169 = pkin(5) * t356 - pkin(10) * t519 + t193;
t194 = t408 * t266 + t411 * t289;
t520 = t358 * t408;
t178 = -pkin(10) * t520 + t194;
t94 = t169 * t416 + t178 * t420;
t667 = -qJD(6) * t94 - t416 * t690 + t420 * t691;
t93 = t169 * t420 - t178 * t416;
t666 = qJD(6) * t93 + t416 * t691 + t420 * t690;
t561 = pkin(3) * t409;
t396 = qJ(5) + t561;
t555 = pkin(10) + t396;
t351 = t555 * t408;
t352 = t555 * t411;
t263 = -t351 * t420 - t352 * t416;
t442 = t408 * t416 - t411 * t420;
t524 = t234 * t411;
t106 = t185 * t412 - t182;
t562 = pkin(3) * t305;
t141 = pkin(4) * t443 + qJ(5) * t234 + t562;
t60 = -t106 * t408 + t411 * t141;
t45 = pkin(5) * t443 + pkin(10) * t524 + t60;
t525 = t234 * t408;
t61 = t411 * t106 + t408 * t141;
t51 = pkin(10) * t525 + t61;
t662 = -qJD(5) * t442 + qJD(6) * t263 - t416 * t45 - t420 * t51;
t264 = -t351 * t416 + t352 * t420;
t359 = t408 * t420 + t411 * t416;
t661 = -qJD(5) * t359 - qJD(6) * t264 + t416 * t51 - t420 * t45;
t172 = t251 * t420 - t252 * t416;
t342 = t442 * qJD(6);
t180 = t341 * t359 + t342 * t358;
t658 = t172 - t180;
t173 = t251 * t416 + t252 * t420;
t343 = t359 * qJD(6);
t179 = t341 * t442 - t343 * t358;
t657 = t173 - t179;
t132 = -mrSges(6,1) * t469 + mrSges(6,2) * t198;
t205 = mrSges(5,1) * t369 - mrSges(5,3) * t443;
t656 = t205 - t132;
t655 = pkin(5) * t644 + t653;
t514 = t410 * t422;
t354 = t413 * t418 * pkin(1) + pkin(8) * t514;
t318 = pkin(9) * t413 + t354;
t245 = t421 * t318 + t417 * t319;
t467 = mrSges(3,3) * t479;
t650 = -mrSges(3,1) * t388 - mrSges(4,1) * t304 + mrSges(4,2) * t305 + t467;
t162 = t442 * t234;
t649 = t342 + t162;
t161 = t359 * t234;
t648 = t343 + t161;
t389 = pkin(8) * t517;
t563 = pkin(1) * t422;
t353 = t413 * t563 - t389;
t481 = t565 * t418;
t510 = t419 * t422;
t347 = t413 * t481 + t510;
t482 = t410 * t565;
t429 = t347 * t417 + t421 * t482;
t145 = -mrSges(6,2) * t234 + mrSges(6,3) * t469;
t146 = mrSges(6,1) * t234 - mrSges(6,3) * t198;
t643 = t145 * t411 - t146 * t408;
t642 = t109 * t421 - t110 * t417;
t641 = -t10 * t408 + t11 * t411;
t638 = m(5) - t669;
t635 = -t408 * (mrSges(6,1) + t681) + mrSges(3,2) - t530 + t665;
t99 = -pkin(4) * t369 + qJD(5) - t101;
t633 = t99 * t454 + t685;
t632 = mrSges(3,2) + t672;
t631 = mrSges(3,1) - t699;
t627 = t410 ^ 2;
t626 = Ifges(7,4) * t623 + Ifges(7,2) * t622 + Ifges(7,6) * t599;
t625 = Ifges(7,1) * t623 + Ifges(7,4) * t622 + Ifges(7,5) * t599;
t624 = m(5) * pkin(3);
t43 = t119 * Ifges(6,4) + t118 * Ifges(6,2) + t155 * Ifges(6,6);
t621 = t43 / 0.2e1;
t620 = Ifges(6,1) * t606 + Ifges(6,4) * t607 + Ifges(6,5) * t597;
t545 = Ifges(7,4) * t125;
t55 = Ifges(7,2) * t671 + Ifges(7,6) * t231 + t545;
t619 = -t55 / 0.2e1;
t618 = t55 / 0.2e1;
t120 = Ifges(7,4) * t671;
t56 = Ifges(7,1) * t125 + Ifges(7,5) * t231 + t120;
t617 = -t56 / 0.2e1;
t616 = t56 / 0.2e1;
t97 = t198 * Ifges(6,4) + Ifges(6,2) * t469 + t234 * Ifges(6,6);
t612 = t97 / 0.2e1;
t98 = t198 * Ifges(6,1) + Ifges(6,4) * t469 + t234 * Ifges(6,5);
t611 = -t98 / 0.2e1;
t610 = t98 / 0.2e1;
t609 = pkin(1) * mrSges(3,1);
t608 = pkin(1) * mrSges(3,2);
t601 = Ifges(4,4) * t584 + Ifges(4,2) * t583 + Ifges(4,6) * t571;
t600 = Ifges(4,1) * t584 + Ifges(4,4) * t583 + Ifges(4,5) * t571;
t593 = -t159 / 0.2e1;
t533 = t305 * Ifges(4,4);
t215 = t304 * Ifges(4,2) + t369 * Ifges(4,6) + t533;
t586 = t215 / 0.2e1;
t296 = Ifges(4,4) * t304;
t216 = t305 * Ifges(4,1) + t369 * Ifges(4,5) + t296;
t585 = t216 / 0.2e1;
t572 = t305 / 0.2e1;
t566 = t413 / 0.2e1;
t560 = pkin(3) * t412;
t501 = qJD(2) * t418;
t332 = qJD(2) * t440;
t334 = t353 * qJD(2);
t175 = -qJD(3) * t245 + t421 * t332 - t334 * t417;
t502 = qJD(2) * t410;
t476 = t422 * t502;
t278 = qJD(3) * t344 + t421 * t476;
t345 = t413 * t417 + t418 * t515;
t477 = t410 * t501;
t117 = pkin(3) * t477 - qJ(4) * t278 - qJD(4) * t345 + t175;
t174 = -t318 * t500 + t319 * t499 + t417 * t332 + t421 * t334;
t277 = -qJD(3) * t345 - t417 * t476;
t133 = qJ(4) * t277 + qJD(4) * t344 + t174;
t66 = t409 * t117 + t412 * t133;
t59 = (qJ(5) * t501 - qJD(5) * t422) * t410 + t66;
t202 = -t412 * t277 + t278 * t409;
t203 = t277 * t409 + t278 * t412;
t335 = t354 * qJD(2);
t240 = -t277 * pkin(3) + t335;
t258 = t344 * t409 + t345 * t412;
t88 = t202 * pkin(4) - t203 * qJ(5) - t258 * qJD(5) + t240;
t27 = t408 * t88 + t411 * t59;
t551 = Ifges(3,4) * t418;
t550 = Ifges(3,4) * t422;
t549 = Ifges(4,4) * t417;
t548 = Ifges(4,4) * t421;
t547 = Ifges(6,4) * t408;
t546 = Ifges(6,4) * t411;
t544 = Ifges(3,6) * t388;
t538 = t221 * mrSges(4,3);
t537 = t222 * mrSges(4,3);
t536 = t234 * Ifges(5,6);
t535 = t443 * Ifges(5,5);
t534 = t304 * Ifges(4,6);
t532 = t305 * Ifges(4,5);
t531 = t388 * Ifges(3,5);
t516 = t410 * t419;
t244 = -t417 * t318 + t421 * t319;
t192 = -pkin(3) * t514 - t345 * qJ(4) + t244;
t206 = qJ(4) * t344 + t245;
t129 = t409 * t192 + t412 * t206;
t114 = -qJ(5) * t514 + t129;
t257 = -t412 * t344 + t345 * t409;
t317 = t389 + (-pkin(2) - t563) * t413;
t265 = t317 - t679;
t154 = t257 * pkin(4) - t258 * qJ(5) + t265;
t71 = t411 * t114 + t408 * t154;
t506 = t565 * pkin(1) + pkin(8) * t516;
t487 = t417 * t516;
t483 = Ifges(3,5) * t337 - Ifges(3,6) * t336 + Ifges(3,3) * t387;
t14 = -t40 * mrSges(7,1) + t39 * mrSges(7,2);
t471 = -pkin(1) * t419 + pkin(8) * t482;
t26 = -t408 * t59 + t411 * t88;
t82 = t155 * mrSges(5,1) + t156 * mrSges(5,2);
t64 = -t118 * mrSges(6,1) + t119 * mrSges(6,2);
t70 = -t114 * t408 + t411 * t154;
t65 = t117 * t412 - t409 * t133;
t105 = t185 * t409 + t513;
t128 = t192 * t412 - t409 * t206;
t270 = t347 * t404 - t402 * t482;
t377 = t417 * t482;
t468 = -t347 * t421 + t377;
t288 = -t412 * t373 - t374 * t409;
t466 = mrSges(3,3) * t478;
t115 = pkin(4) * t514 - t128;
t458 = mrSges(4,1) * t344 - mrSges(4,2) * t345;
t452 = Ifges(4,1) * t421 - t549;
t451 = Ifges(6,1) * t411 - t547;
t450 = t422 * Ifges(3,2) + t551;
t449 = -Ifges(4,2) * t417 + t548;
t448 = -Ifges(6,2) * t408 + t546;
t447 = Ifges(4,5) * t421 - Ifges(4,6) * t417;
t446 = Ifges(6,5) * t411 - Ifges(6,6) * t408;
t445 = -t408 * t52 + t411 * t53;
t230 = t411 * t258 - t408 * t514;
t50 = pkin(5) * t257 - pkin(10) * t230 + t70;
t229 = -t408 * t258 - t411 * t514;
t57 = pkin(10) * t229 + t71;
t17 = -t416 * t57 + t420 * t50;
t18 = t416 * t50 + t420 * t57;
t163 = t229 * t420 - t230 * t416;
t164 = t229 * t416 + t230 * t420;
t269 = t347 * t402 + t404 * t482;
t273 = t349 * t402 - t404 * t516;
t315 = t402 * t517 - t413 * t404;
t433 = -g(1) * t273 - g(2) * t269 - g(3) * t315;
t30 = -pkin(4) * t323 + qJDD(5) - t31;
t346 = -t413 * t480 + t511;
t348 = t413 * t510 + t481;
t427 = -g(1) * t348 - g(2) * t346 + g(3) * t514;
t62 = -pkin(4) * t477 - t65;
t399 = -pkin(4) - t560;
t382 = Ifges(3,4) * t478;
t365 = -t398 - t560;
t350 = (-mrSges(3,1) * t422 + mrSges(3,2) * t418) * t410;
t329 = -t388 * mrSges(3,2) + t466;
t316 = t402 * t413 + t404 * t517;
t286 = Ifges(3,1) * t479 + t382 + t531;
t285 = t450 * t505 + t544;
t281 = t349 * t421 + t487;
t274 = t349 * t404 + t402 * t516;
t268 = mrSges(4,1) * t369 - mrSges(4,3) * t305;
t267 = -mrSges(4,2) * t369 + mrSges(4,3) * t304;
t256 = t442 * t358;
t255 = t359 * t358;
t248 = pkin(5) * t520 + t288;
t214 = t369 * Ifges(4,3) + t532 + t534;
t212 = t274 * t403 + t348 * t401;
t211 = -t274 * t401 + t348 * t403;
t204 = -mrSges(5,2) * t369 - mrSges(5,3) * t234;
t190 = -mrSges(4,2) * t323 + mrSges(4,3) * t227;
t189 = mrSges(4,1) * t323 - mrSges(4,3) * t226;
t188 = t203 * t411 + t408 * t477;
t187 = -t203 * t408 + t411 * t477;
t168 = mrSges(5,1) * t234 + mrSges(5,2) * t443;
t160 = -mrSges(4,1) * t227 + mrSges(4,2) * t226;
t157 = t369 * Ifges(5,3) + t535 - t536;
t131 = mrSges(5,1) * t323 - mrSges(5,3) * t156;
t130 = -mrSges(5,2) * t323 - mrSges(5,3) * t155;
t92 = mrSges(7,1) * t231 - mrSges(7,3) * t125;
t91 = -mrSges(7,2) * t231 + mrSges(7,3) * t671;
t90 = -pkin(5) * t229 + t115;
t89 = -pkin(5) * t525 + t105;
t78 = -pkin(5) * t469 + t99;
t74 = mrSges(6,1) * t155 - mrSges(6,3) * t119;
t73 = -mrSges(6,2) * t155 + mrSges(6,3) * t118;
t69 = -qJD(6) * t164 + t187 * t420 - t188 * t416;
t68 = qJD(6) * t163 + t187 * t416 + t188 * t420;
t67 = -mrSges(7,1) * t671 + mrSges(7,2) * t125;
t49 = -pkin(5) * t187 + t62;
t25 = -mrSges(7,2) * t153 + mrSges(7,3) * t40;
t24 = mrSges(7,1) * t153 - mrSges(7,3) * t39;
t21 = pkin(10) * t187 + t27;
t20 = -pkin(5) * t118 + t30;
t19 = pkin(5) * t202 - pkin(10) * t188 + t26;
t4 = -qJD(6) * t18 + t19 * t420 - t21 * t416;
t3 = qJD(6) * t17 + t19 * t416 + t21 * t420;
t8 = [(t678 / 0.2e1 + t700) * t202 + (-t674 / 0.2e1 + mrSges(3,3) * t261 - Ifges(4,5) * t584 - Ifges(5,5) * t596 - Ifges(4,6) * t583 - Ifges(5,6) * t598 - t698 * t571 + t670) * t514 + (t684 + t685) * t203 + (Ifges(7,4) * t164 + Ifges(7,2) * t163) * t622 + (Ifges(7,4) * t68 + Ifges(7,2) * t69) * t604 + (t353 * mrSges(3,1) - t354 * mrSges(3,2) + Ifges(3,3) * t566 + (Ifges(3,5) * t418 + Ifges(3,6) * t422) * t410) * t387 + (-t353 * mrSges(3,3) + Ifges(3,5) * t566 + (t418 * Ifges(3,1) + t550 - t608) * t410) * t337 + ((-t330 * mrSges(3,3) + t286 / 0.2e1 + t531 / 0.2e1 + (-t608 + t550 / 0.2e1) * t505) * t422 + (-t333 * mrSges(3,3) + t214 / 0.2e1 + t157 / 0.2e1 - t285 / 0.2e1 + t534 / 0.2e1 + t532 / 0.2e1 + t221 * mrSges(4,1) - t222 * mrSges(4,2) + t535 / 0.2e1 - t536 / 0.2e1 - t102 * mrSges(5,2) + t101 * mrSges(5,1) - t544 / 0.2e1 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t369 + (-t609 - t551 / 0.2e1 + (Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1) * t422) * t505) * t418) * t502 - (t354 * mrSges(3,3) + Ifges(3,6) * t566 + (t450 + t609) * t410) * t336 + (t1 * t163 - t12 * t68 + t13 * t69 - t164 * t2) * mrSges(7,3) + (Ifges(6,5) * t188 + Ifges(6,6) * t187) * t578 + (Ifges(6,5) * t230 + Ifges(6,6) * t229) * t597 + (Ifges(4,1) * t345 + Ifges(4,4) * t344) * t584 + (Ifges(7,5) * t68 + Ifges(7,6) * t69) * t581 + (Ifges(7,5) * t164 + Ifges(7,6) * t163) * t599 + (-t10 * t230 + t11 * t229 + t187 * t53 - t188 * t52) * mrSges(6,3) + (Ifges(4,4) * t345 + Ifges(4,2) * t344) * t583 + (-pkin(1) * t350 * t410 + Ifges(2,3)) * qJDD(1) + t693 * t258 + t694 * t257 + t650 * t335 + m(6) * (t10 * t70 + t11 * t71 + t115 * t30 + t26 * t52 + t27 * t53 + t62 * t99) + m(7) * (t1 * t18 + t12 * t4 + t13 * t3 + t17 * t2 + t20 * t90 + t49 * t78) + m(5) * (t101 * t65 + t102 * t66 + t128 * t31 + t129 * t32 + t167 * t265 + t232 * t240) + m(4) * (t109 * t245 + t110 * t244 + t174 * t222 + t175 * t221 + t242 * t317 + t290 * t335) - t242 * t458 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t627 + t261 * t354 + t262 * t353 - t330 * t335 + t333 * t334) + (Ifges(7,1) * t164 + Ifges(7,4) * t163) * t623 + (Ifges(7,1) * t68 + Ifges(7,4) * t69) * t602 + (-t565 * mrSges(2,1) - m(4) * (pkin(2) * t349 + t506) - t281 * mrSges(4,1) - t280 * mrSges(4,2) - t212 * mrSges(7,1) - t211 * mrSges(7,2) - m(3) * t506 - t349 * mrSges(3,1) + (-mrSges(3,3) * t410 + mrSges(2,2)) * t419 + t425 * t274 + t635 * t348 + t430 * t273 - t638 * (pkin(3) * t487 - t348 * t414 + t349 * t400 + t506)) * g(2) + (Ifges(4,5) * t345 + Ifges(4,6) * t344) * t571 + (Ifges(6,4) * t188 + Ifges(6,2) * t187) * t590 + (Ifges(6,1) * t230 + Ifges(6,4) * t229) * t606 + (Ifges(6,1) * t188 + Ifges(6,4) * t187) * t588 + t262 * (mrSges(3,1) * t413 - mrSges(3,3) * t517) + (t109 * t344 - t110 * t345 - t221 * t278 + t222 * t277) * mrSges(4,3) + t334 * t329 + t317 * t160 + t290 * (-mrSges(4,1) * t277 + mrSges(4,2) * t278) + t174 * t267 + t175 * t268 + t265 * t82 + t244 * t189 + t245 * t190 + t240 * t168 + t30 * (-mrSges(6,1) * t229 + mrSges(6,2) * t230) + t66 * t204 + t65 * t205 + t99 * (-mrSges(6,1) * t187 + mrSges(6,2) * t188) + t20 * (-mrSges(7,1) * t163 + mrSges(7,2) * t164) + t27 * t145 + t26 * t146 + t62 * t132 + t129 * t130 + t128 * t131 + t115 * t64 + t90 * t14 + t3 * t91 + t4 * t92 + t78 * (-mrSges(7,1) * t69 + mrSges(7,2) * t68) + t70 * t74 + t71 * t73 + t49 * t67 - t413 * t677 + t163 * t626 + t164 * t625 + t68 * t616 + t69 * t618 + t230 * t620 + t229 * t621 + t187 * t612 + t188 * t610 + (-m(3) * t471 + t347 * mrSges(3,1) - mrSges(3,3) * t482 + t419 * mrSges(2,1) + t565 * mrSges(2,2) - m(4) * (-pkin(2) * t347 + t471) - t468 * mrSges(4,1) - t429 * mrSges(4,2) - t707 * t270 + (t453 - t635) * t346 - t430 * t269 + t638 * (-pkin(3) * t377 - t346 * t414 + t347 * t400 - t471)) * g(1) + t17 * t24 + t18 * t25 + t483 * t566 + (Ifges(4,5) * t278 + Ifges(4,6) * t277) * t567 + (Ifges(4,1) * t278 + Ifges(4,4) * t277) * t572 + t278 * t585 + t277 * t586 + t345 * t600 + t344 * t601 + t304 * (Ifges(4,4) * t278 + Ifges(4,2) * t277) / 0.2e1 + (Ifges(6,4) * t230 + Ifges(6,2) * t229) * t607; (-t10 * t519 - t11 * t520 + t52 * t701 - t644 * t53) * mrSges(6,3) + t700 * t340 + (t350 - t638 * t400 * t514 + (t699 * t422 + (t414 * t638 + t672) * t418) * t410) * g(3) - (t446 * t578 + t448 * t590 + t451 * t588 + t633 + t684) * t341 + (t466 - t329) * t330 + (-t538 + t585) * t499 + (Ifges(7,4) * t179 + Ifges(7,2) * t180) * t604 + (Ifges(7,4) * t173 + Ifges(7,2) * t172) * t605 + (Ifges(7,5) * t179 + Ifges(7,6) * t180) * t581 + (Ifges(7,5) * t173 + Ifges(7,6) * t172) * t582 + (Ifges(6,1) * t252 + Ifges(6,4) * t251) * t589 + (Ifges(7,1) * t179 + Ifges(7,4) * t180) * t602 + (Ifges(7,1) * t173 + Ifges(7,4) * t172) * t603 - (t305 * (Ifges(4,5) * t418 + t422 * t452) + t304 * (Ifges(4,6) * t418 + t422 * t449) + t388 * (Ifges(3,5) * t422 - Ifges(3,6) * t418) + t369 * (Ifges(4,3) * t418 + t422 * t447) + (-Ifges(3,2) * t479 + t421 * t216 + t286 + t382) * t422 + (t214 + t157) * t418) * t505 / 0.2e1 + (-t537 - t215 / 0.2e1) * t500 + (Ifges(6,5) * t252 + Ifges(6,6) * t251) * t579 + (-Ifges(5,4) * t284 + Ifges(5,6) * t479) * t578 + (t102 * t479 + t232 * t284) * mrSges(5,2) + (Ifges(6,4) * t252 + Ifges(6,2) * t251) * t591 - t43 * t520 / 0.2e1 + (m(4) * ((-t221 * t421 - t222 * t417) * qJD(3) + t642) - t268 * t499 - t267 * t500 + t421 * t190 - t417 * t189) * pkin(9) + t642 * mrSges(4,3) - t644 * t97 / 0.2e1 + (-t222 * (-mrSges(4,3) * t417 * t422 - mrSges(4,2) * t418) - t221 * (mrSges(4,1) * t418 - mrSges(4,3) * t509)) * t505 + (pkin(1) * (mrSges(3,1) * t418 + mrSges(3,2) * t422) - t418 * (Ifges(3,1) * t422 - t551) / 0.2e1) * qJD(1) ^ 2 * t627 + t666 * t91 + (t1 * t94 + t12 * t667 + t13 * t666 + t2 * t93 + t20 * t248 + t655 * t78) * m(7) + t667 * t92 + (t30 * t454 + t446 * t597 + t448 * t607 + t451 * t606 + t693) * t358 + t694 * t356 + (-Ifges(7,4) * t256 - Ifges(7,2) * t255) * t622 + (-Ifges(7,5) * t256 - Ifges(7,6) * t255) * t599 + (-Ifges(7,1) * t256 - Ifges(7,4) * t255) * t623 + t20 * (mrSges(7,1) * t255 - mrSges(7,2) * t256) + (-Ifges(5,1) * t284 + Ifges(5,5) * t479) * t577 + (-Ifges(5,5) * t284 + Ifges(5,3) * t479) * t568 - t101 * (mrSges(5,1) * t479 + mrSges(5,3) * t284) + (-m(4) * t290 + t467 - t650) * t333 + t651 * t168 + t652 * t204 + t242 * t457 + (t64 - t131) * t288 + (t304 * t449 + t305 * t452 + t369 * t447) * qJD(3) / 0.2e1 - t400 * t82 + t483 + t653 * t132 + (t101 * t654 + t102 * t652 - t167 * t400 + t232 * t651 - t288 * t31 + t289 * t32) * m(5) + t654 * t205 + t655 * t67 + (mrSges(7,1) * t658 - mrSges(7,2) * t657) * t78 + (-t1 * t255 + t12 * t657 - t13 * t658 + t2 * t256) * mrSges(7,3) - t677 + t285 * t479 / 0.2e1 + (-t638 * (-t348 * t400 - t349 * t414) + t632 * t349 + t631 * t348) * g(1) + (-t638 * (-t346 * t400 - t347 * t414) + t632 * t347 + t631 * t346) * g(2) + t369 * t290 * (mrSges(4,1) * t417 + mrSges(4,2) * t421) + t289 * t130 - t250 * t267 - t249 * t268 + t262 * mrSges(3,1) - t99 * (-mrSges(6,1) * t251 + mrSges(6,2) * t252) + t248 * t14 + t193 * t74 + t194 * t73 + (-pkin(2) * t242 - t221 * t249 - t222 * t250) * m(4) - pkin(2) * t160 + t663 * t145 + t664 * t146 + (t10 * t193 + t11 * t194 + t288 * t30 + t52 * t664 + t53 * t663 + t653 * t99) * m(6) + t94 * t25 + t93 * t24 + t678 * (t340 / 0.2e1 - t283 / 0.2e1) + t673 * t283 - t256 * t625 - t255 * t626 + t173 * t617 + t180 * t618 + t172 * t619 + t519 * t620 + t252 * t611 + t179 * t616 - t523 * t610 + (Ifges(4,5) * t417 + Ifges(4,6) * t421) * t571 + (Ifges(4,2) * t421 + t549) * t583 + (Ifges(4,1) * t417 + t548) * t584 + t462 * t586 - t284 * t593 + t417 * t600 + t421 * t601; (mrSges(4,2) * t281 - t707 * t273 + t630 * t274 + (-pkin(3) * t638 - mrSges(4,1)) * t280) * g(1) + (-mrSges(4,2) * t468 - t707 * t269 + t630 * t270 + (-pkin(3) * t669 + mrSges(4,1) + t624) * t429) * g(2) + (-t315 * t707 + t630 * t316 - t638 * t679 - t458) * g(3) + (Ifges(7,5) * t162 + Ifges(7,6) * t161) * t582 + (-Ifges(7,5) * t342 - Ifges(7,6) * t343) * t581 - (Ifges(5,1) * t577 + Ifges(5,4) * t578 + Ifges(5,5) * t568 + t446 * t579 + t448 * t591 + t451 * t589 + t593 - t633) * t234 + t20 * (mrSges(7,1) * t442 + mrSges(7,2) * t359) + (Ifges(7,4) * t359 - Ifges(7,2) * t442) * t622 + (Ifges(7,1) * t359 - Ifges(7,4) * t442) * t623 + (Ifges(7,5) * t359 - Ifges(7,6) * t442) * t599 + (-g(1) * t274 - g(2) * t270 - g(3) * t316 - t1 * t442 + t12 * t649 - t13 * t648 - t2 * t359) * mrSges(7,3) - t442 * t626 - t168 * t562 + (-t408 * t74 + t411 * t73) * t396 + (-Ifges(7,1) * t342 - Ifges(7,4) * t343) * t602 + (-Ifges(7,4) * t342 - Ifges(7,2) * t343) * t604 + (-t52 * t524 - t525 * t53 + t641) * mrSges(6,3) + (t445 * qJD(5) - t105 * t99 + t30 * t399 + t396 * t641 - t52 * t60 - t53 * t61) * m(6) + t643 * qJD(5) + (Ifges(7,4) * t162 + Ifges(7,2) * t161) * t605 + (Ifges(7,1) * t162 + Ifges(7,4) * t161) * t603 - m(5) * (t102 * t106 + t232 * t562) + (mrSges(7,1) * t648 - mrSges(7,2) * t649) * t78 + t30 * t455 - t670 - (-Ifges(4,2) * t305 + t216 + t296) * t304 / 0.2e1 + t399 * t64 + (m(5) * t101 + t656) * t105 + t365 * t14 - t305 * (Ifges(4,1) * t304 - t533) / 0.2e1 - t290 * (mrSges(4,1) * t305 + mrSges(4,2) * t304) - t221 * t267 + t222 * t268 + t263 * t24 + t264 * t25 - t106 * t204 - t61 * t145 - t60 * t146 + t661 * t92 + t662 * t91 + (t1 * t264 + t12 * t661 + t13 * t662 + t2 * t263 + t20 * t365 - t78 * t89) * m(7) - t89 * t67 + t678 * t577 + t673 * t443 + (t31 * t412 + t32 * t409) * t624 + t359 * t625 + t162 * t617 - t343 * t618 + t161 * t619 + t408 * t620 + t411 * t621 - t524 * t611 - t525 * t612 - t342 * t616 + (Ifges(6,1) * t408 + t546) * t606 + (Ifges(6,2) * t411 + t547) * t607 + t305 * t537 + t304 * t538 + t674 + t131 * t560 + t130 * t561 + (Ifges(4,5) * t304 - Ifges(4,6) * t305) * t568 + t215 * t572 + (Ifges(6,5) * t408 + Ifges(6,6) * t411) * t597; -t442 * t24 + t359 * t25 + t408 * t73 + t411 * t74 - t648 * t92 - t649 * t91 - (-t204 - t643) * t234 + (-t67 + t656) * t443 + t82 + (t1 * t359 - t12 * t648 - t13 * t649 - t2 * t442 - t443 * t78 + t427) * m(7) + (t10 * t411 + t11 * t408 + t234 * t445 - t443 * t99 + t427) * m(6) + (t101 * t443 + t102 * t234 + t167 + t427) * m(5); t125 * t92 - t671 * t91 - t469 * t145 + t198 * t146 + t14 + t64 + (t12 * t125 - t13 * t671 + t20 + t433) * m(7) + (t198 * t52 - t469 * t53 + t30 + t433) * m(6); -t78 * (mrSges(7,1) * t125 + mrSges(7,2) * t671) + (Ifges(7,1) * t671 - t545) * t603 + t55 * t602 + (Ifges(7,5) * t671 - Ifges(7,6) * t125) * t582 - t12 * t91 + t13 * t92 - g(1) * (mrSges(7,1) * t211 - mrSges(7,2) * t212) - g(2) * ((-t270 * t401 + t346 * t403) * mrSges(7,1) + (-t270 * t403 - t346 * t401) * mrSges(7,2)) - g(3) * ((-t316 * t401 - t403 * t514) * mrSges(7,1) + (-t316 * t403 + t401 * t514) * mrSges(7,2)) + (t12 * t671 + t125 * t13) * mrSges(7,3) + t7 + (-Ifges(7,2) * t125 + t120 + t56) * t605 + t634;];
tau  = t8;
