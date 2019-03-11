% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR12_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR12_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR12_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR12_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR12_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:38:13
% EndTime: 2019-03-09 19:40:06
% DurationCPUTime: 70.88s
% Computational Cost: add. (37654->1200), mult. (89842->1650), div. (0->0), fcn. (73827->18), ass. (0->476)
t455 = sin(qJ(2));
t577 = cos(pkin(6));
t534 = pkin(1) * t577;
t435 = t455 * t534;
t454 = sin(qJ(3));
t458 = cos(qJ(3));
t488 = pkin(3) * t454 - qJ(4) * t458;
t449 = sin(pkin(6));
t459 = cos(qJ(2));
t563 = t449 * t459;
t748 = (t435 + (pkin(8) + t488) * t563) * qJD(1) - qJD(3) * t488 + qJD(4) * t454;
t515 = t577 * qJD(1);
t506 = pkin(1) * t515;
t553 = qJD(1) * t449;
t530 = t455 * t553;
t366 = -pkin(8) * t530 + t459 * t506;
t479 = (pkin(2) * t455 - pkin(9) * t459) * t449;
t367 = qJD(1) * t479;
t278 = t458 * t366 + t454 * t367;
t252 = qJ(4) * t530 + t278;
t448 = sin(pkin(12));
t450 = cos(pkin(12));
t549 = qJD(3) * t454;
t541 = pkin(9) * t549;
t688 = -t748 * t450 + (t252 + t541) * t448;
t747 = t450 * t252 + t448 * t748;
t560 = t458 * t459;
t331 = (t448 * t455 + t450 * t560) * t553;
t529 = t459 * t553;
t507 = t454 * t529;
t561 = t450 * t458;
t746 = -pkin(4) * t507 + t331 * pkin(10) + (pkin(4) * t454 - pkin(10) * t561) * qJD(3) + t688;
t330 = (-t448 * t560 + t450 * t455) * t553;
t562 = t450 * t454;
t565 = t448 * t458;
t745 = pkin(10) * t330 - (-pkin(9) * t562 - pkin(10) * t565) * qJD(3) + t747;
t453 = sin(qJ(5));
t457 = cos(qJ(5));
t240 = t330 * t457 - t331 * t453;
t486 = t448 * t453 - t450 * t457;
t379 = t486 * qJD(5);
t399 = t448 * t457 + t450 * t453;
t548 = qJD(3) * t458;
t282 = t379 * t454 - t399 * t548;
t685 = t240 - t282;
t241 = t330 * t453 + t331 * t457;
t380 = t399 * qJD(5);
t281 = -t380 * t454 - t486 * t548;
t684 = t241 - t281;
t484 = t515 + qJD(2);
t344 = t454 * t530 - t458 * t484;
t242 = t399 * t344;
t683 = t242 + t380;
t243 = t486 * t344;
t682 = t243 + t379;
t681 = t507 - t549;
t405 = -pkin(3) * t458 - qJ(4) * t454 - pkin(2);
t393 = t450 * t405;
t301 = -pkin(10) * t562 + t393 + (-pkin(9) * t448 - pkin(4)) * t458;
t343 = pkin(9) * t561 + t448 * t405;
t566 = t448 * t454;
t313 = -pkin(10) * t566 + t343;
t544 = qJD(5) * t457;
t545 = qJD(5) * t453;
t697 = t301 * t544 - t313 * t545 + t453 * t746 - t745 * t457;
t210 = t453 * t301 + t457 * t313;
t696 = -qJD(5) * t210 + t745 * t453 + t457 * t746;
t555 = pkin(8) * t563 + t435;
t357 = t577 * pkin(9) + t555;
t328 = qJD(2) * pkin(9) + qJD(1) * t357;
t336 = (-pkin(2) * t459 - pkin(9) * t455 - pkin(1)) * t553;
t237 = -t454 * t328 + t336 * t458;
t345 = t454 * t484 + t458 * t530;
t263 = pkin(3) * t345 + qJ(4) * t344;
t167 = -t237 * t448 + t450 * t263;
t571 = t344 * t450;
t130 = pkin(4) * t345 + pkin(10) * t571 + t167;
t168 = t450 * t237 + t448 * t263;
t572 = t344 * t448;
t149 = pkin(10) * t572 + t168;
t596 = pkin(10) + qJ(4);
t406 = t596 * t448;
t407 = t596 * t450;
t323 = -t453 * t406 + t457 * t407;
t695 = -t399 * qJD(4) - qJD(5) * t323 - t457 * t130 + t149 * t453;
t546 = qJD(4) * t450;
t547 = qJD(4) * t448;
t694 = -t406 * t544 + (-t149 + t546) * t457 + (-qJD(5) * t407 - t130 - t547) * t453;
t447 = pkin(12) + qJ(5);
t443 = qJ(6) + t447;
t437 = sin(t443);
t438 = cos(t443);
t441 = sin(t447);
t442 = cos(t447);
t744 = -mrSges(6,1) * t442 - mrSges(7,1) * t438 + mrSges(6,2) * t441 + mrSges(7,2) * t437;
t439 = pkin(4) * t450 + pkin(3);
t498 = -mrSges(5,1) * t450 + mrSges(5,2) * t448;
t743 = -m(5) * pkin(3) + t498 - m(7) * (pkin(5) * t442 + t439) - m(6) * t439;
t499 = t458 * mrSges(4,1) - t454 * mrSges(4,2);
t670 = m(7) * (-pkin(11) - t596) - mrSges(7,3) - m(6) * t596 - mrSges(6,3) - m(5) * qJ(4) - mrSges(5,3);
t671 = -t743 - t744;
t742 = t454 * t670 - t458 * t671 - t499;
t411 = qJD(3) - t529;
t287 = t345 * t450 + t411 * t448;
t513 = -t345 * t448 + t450 * t411;
t191 = t287 * t457 + t453 * t513;
t452 = sin(qJ(6));
t456 = cos(qJ(6));
t720 = -t287 * t453 + t457 * t513;
t741 = -t191 * t452 + t456 * t720;
t118 = t191 * t456 + t452 * t720;
t715 = -m(5) - m(4);
t740 = -pkin(5) * t681 + pkin(11) * t684 + t696;
t739 = pkin(11) * t685 - t697;
t497 = t448 * mrSges(5,1) + t450 * mrSges(5,2);
t738 = mrSges(4,3) + t497;
t737 = -pkin(11) * t683 + t694;
t736 = -pkin(5) * t345 + pkin(11) * t682 + t695;
t277 = -t454 * t366 + t367 * t458;
t253 = -pkin(3) * t530 - t277;
t440 = pkin(9) * t548;
t731 = t440 - t253;
t730 = -t441 * mrSges(6,1) - t437 * mrSges(7,1) - t442 * mrSges(6,2) - t438 * mrSges(7,2);
t327 = -pkin(2) * t484 - t366;
t206 = t344 * pkin(3) - t345 * qJ(4) + t327;
t238 = t458 * t328 + t454 * t336;
t211 = qJ(4) * t411 + t238;
t137 = t450 * t206 - t211 * t448;
t138 = t448 * t206 + t450 * t211;
t338 = qJD(5) + t344;
t106 = pkin(10) * t513 + t138;
t93 = pkin(4) * t344 - pkin(10) * t287 + t137;
t52 = -t106 * t453 + t457 * t93;
t46 = -pkin(11) * t191 + t52;
t45 = pkin(5) * t338 + t46;
t53 = t106 * t457 + t453 * t93;
t47 = pkin(11) * t720 + t53;
t579 = t452 * t47;
t16 = t45 * t456 - t579;
t578 = t456 * t47;
t17 = t45 * t452 + t578;
t583 = t238 * mrSges(4,3);
t699 = t411 * Ifges(4,6);
t729 = -t583 - t699 / 0.2e1 + t327 * mrSges(4,1) + t137 * mrSges(5,1) + t52 * mrSges(6,1) + t16 * mrSges(7,1) - t138 * mrSges(5,2) - t53 * mrSges(6,2) - t17 * mrSges(7,2);
t728 = Ifges(5,6) * t513;
t727 = t287 * Ifges(5,5);
t726 = mrSges(4,2) + t670;
t601 = pkin(4) * t448;
t686 = pkin(4) * t330 + t548 * t601 + t731;
t542 = qJDD(1) * t449;
t725 = pkin(8) * t542 + qJD(2) * t506;
t512 = t577 * qJDD(1);
t543 = qJD(1) * qJD(2);
t724 = -pkin(8) * t449 * t543 + pkin(1) * t512;
t714 = -m(6) - m(7);
t672 = -t714 - t715;
t723 = pkin(2) * t672 + mrSges(3,1) - t742;
t402 = pkin(5) * t441 + t601;
t662 = -mrSges(3,2) + m(7) * (pkin(9) + t402) + m(6) * (pkin(9) + t601) + t738;
t373 = (qJDD(1) * t455 + t459 * t543) * t449;
t425 = t512 + qJDD(2);
t550 = qJD(3) * t344;
t244 = t458 * t373 + t454 * t425 - t550;
t372 = (-qJDD(1) * t459 + t455 * t543) * t449;
t361 = qJDD(3) + t372;
t197 = t244 * t450 + t361 * t448;
t245 = qJD(3) * t345 + t454 * t373 - t458 * t425;
t290 = t455 * t724 + t459 * t725;
t261 = pkin(9) * t425 + t290;
t540 = pkin(1) * t542;
t269 = pkin(2) * t372 - pkin(9) * t373 - t540;
t133 = t458 * t261 + t454 * t269 - t328 * t549 + t336 * t548;
t105 = qJ(4) * t361 + qJD(4) * t411 + t133;
t291 = -t455 * t725 + t459 * t724;
t262 = -t425 * pkin(2) - t291;
t111 = t245 * pkin(3) - t244 * qJ(4) - t345 * qJD(4) + t262;
t60 = -t105 * t448 + t450 * t111;
t43 = pkin(4) * t245 - pkin(10) * t197 + t60;
t196 = -t244 * t448 + t361 * t450;
t61 = t450 * t105 + t448 * t111;
t51 = pkin(10) * t196 + t61;
t13 = -qJD(5) * t53 + t457 * t43 - t453 * t51;
t236 = qJDD(5) + t245;
t83 = qJD(5) * t720 + t196 * t453 + t197 * t457;
t6 = pkin(5) * t236 - pkin(11) * t83 + t13;
t12 = -t106 * t545 + t453 * t43 + t457 * t51 + t93 * t544;
t84 = -qJD(5) * t191 + t196 * t457 - t197 * t453;
t7 = pkin(11) * t84 + t12;
t2 = qJD(6) * t16 + t452 * t6 + t456 * t7;
t3 = -qJD(6) * t17 - t452 * t7 + t456 * t6;
t719 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t718 = t13 * mrSges(6,1) - t12 * mrSges(6,2);
t582 = t345 * Ifges(4,4);
t227 = -t344 * Ifges(4,2) + t582 + t699;
t614 = t338 / 0.2e1;
t332 = qJD(6) + t338;
t616 = t332 / 0.2e1;
t629 = t191 / 0.2e1;
t631 = t720 / 0.2e1;
t636 = t118 / 0.2e1;
t638 = t741 / 0.2e1;
t673 = t191 * Ifges(6,5) + t118 * Ifges(7,5) + Ifges(6,6) * t720 + Ifges(7,6) * t741 + t344 * Ifges(5,3) + t338 * Ifges(6,3) + t332 * Ifges(7,3) + t727 + t728;
t716 = -t227 / 0.2e1 + Ifges(6,5) * t629 + Ifges(7,5) * t636 + Ifges(6,6) * t631 + Ifges(7,6) * t638 + Ifges(6,3) * t614 + Ifges(7,3) * t616 + t673 / 0.2e1;
t31 = qJD(6) * t741 + t452 * t84 + t456 * t83;
t655 = t31 / 0.2e1;
t32 = -qJD(6) * t118 - t452 * t83 + t456 * t84;
t654 = t32 / 0.2e1;
t647 = t83 / 0.2e1;
t646 = t84 / 0.2e1;
t628 = t196 / 0.2e1;
t627 = t197 / 0.2e1;
t225 = qJDD(6) + t236;
t626 = t225 / 0.2e1;
t624 = t236 / 0.2e1;
t623 = t244 / 0.2e1;
t622 = -t245 / 0.2e1;
t621 = t245 / 0.2e1;
t609 = t361 / 0.2e1;
t209 = t457 * t301 - t313 * t453;
t363 = t486 * t454;
t181 = -pkin(5) * t458 + pkin(11) * t363 + t209;
t362 = t399 * t454;
t186 = -pkin(11) * t362 + t210;
t107 = t181 * t456 - t186 * t452;
t708 = qJD(6) * t107 + t452 * t740 - t456 * t739;
t108 = t181 * t452 + t186 * t456;
t707 = -qJD(6) * t108 + t452 * t739 + t456 * t740;
t322 = -t457 * t406 - t407 * t453;
t279 = -pkin(11) * t399 + t322;
t280 = -pkin(11) * t486 + t323;
t182 = t279 * t456 - t280 * t452;
t706 = qJD(6) * t182 + t452 * t736 + t456 * t737;
t183 = t279 * t452 + t280 * t456;
t705 = -qJD(6) * t183 - t452 * t737 + t456 * t736;
t700 = t411 * Ifges(4,5);
t657 = m(7) * pkin(5);
t698 = -mrSges(6,1) - t657;
t564 = t449 * t455;
t388 = -pkin(8) * t564 + t459 * t534;
t356 = -pkin(2) * t577 - t388;
t381 = t454 * t564 - t458 * t577;
t382 = t454 * t577 + t458 * t564;
t249 = t381 * pkin(3) - t382 * qJ(4) + t356;
t556 = pkin(2) * t563 + pkin(9) * t564;
t602 = pkin(1) * t449;
t358 = -t556 - t602;
t266 = t458 * t357 + t454 * t358;
t250 = -qJ(4) * t563 + t266;
t162 = t450 * t249 - t250 * t448;
t310 = t382 * t450 - t448 * t563;
t125 = pkin(4) * t381 - pkin(10) * t310 + t162;
t163 = t448 * t249 + t450 * t250;
t309 = -t382 * t448 - t450 * t563;
t139 = pkin(10) * t309 + t163;
t69 = t453 * t125 + t457 * t139;
t270 = -t362 * t456 + t363 * t452;
t151 = qJD(6) * t270 + t281 * t456 + t282 * t452;
t159 = t240 * t452 + t241 * t456;
t691 = t151 - t159;
t271 = -t362 * t452 - t363 * t456;
t152 = -qJD(6) * t271 - t281 * t452 + t282 * t456;
t158 = t240 * t456 - t241 * t452;
t690 = t152 - t158;
t689 = pkin(5) * t685 + t686;
t687 = -t450 * t541 - t747;
t173 = t287 * Ifges(5,4) + Ifges(5,2) * t513 + Ifges(5,6) * t344;
t634 = -t173 / 0.2e1;
t680 = -t237 * mrSges(4,3) + t448 * t634;
t134 = -t454 * t261 + t269 * t458 - t328 * t548 - t336 * t549;
t112 = -pkin(3) * t361 + qJDD(4) - t134;
t208 = -pkin(3) * t411 + qJD(4) - t237;
t679 = t112 * t454 + t208 * t548;
t678 = t133 * t458 - t134 * t454;
t677 = -t448 * t60 + t450 * t61;
t195 = -pkin(4) * t572 + t238;
t676 = pkin(5) * t683 - t195;
t37 = Ifges(6,5) * t83 + Ifges(6,6) * t84 + Ifges(6,3) * t236;
t8 = Ifges(7,5) * t31 + Ifges(7,6) * t32 + Ifges(7,3) * t225;
t675 = t197 * Ifges(5,5) + t196 * Ifges(5,6) + t245 * Ifges(5,3) + t37 + t8;
t593 = Ifges(3,4) * t455;
t660 = t449 ^ 2;
t674 = (pkin(1) * (mrSges(3,1) * t455 + mrSges(3,2) * t459) - t455 * (Ifges(3,1) * t459 - t593) / 0.2e1) * t660;
t669 = mrSges(4,1) + t671;
t510 = mrSges(3,3) * t530;
t666 = -m(4) * t327 + mrSges(3,1) * t484 - mrSges(4,1) * t344 - mrSges(4,2) * t345 - t510;
t664 = pkin(9) * t715 - t662 + t730;
t661 = -mrSges(4,1) + t743;
t658 = Ifges(7,4) * t655 + Ifges(7,2) * t654 + Ifges(7,6) * t626;
t656 = Ifges(7,1) * t655 + Ifges(7,4) * t654 + Ifges(7,5) * t626;
t653 = Ifges(6,4) * t647 + Ifges(6,2) * t646 + Ifges(6,6) * t624;
t652 = Ifges(6,1) * t647 + Ifges(6,4) * t646 + Ifges(6,5) * t624;
t586 = Ifges(7,4) * t118;
t58 = Ifges(7,2) * t741 + t332 * Ifges(7,6) + t586;
t651 = -t58 / 0.2e1;
t650 = t58 / 0.2e1;
t113 = Ifges(7,4) * t741;
t59 = t118 * Ifges(7,1) + t332 * Ifges(7,5) + t113;
t649 = -t59 / 0.2e1;
t648 = t59 / 0.2e1;
t87 = t197 * Ifges(5,4) + t196 * Ifges(5,2) + t245 * Ifges(5,6);
t645 = t87 / 0.2e1;
t644 = Ifges(5,1) * t627 + Ifges(5,4) * t628 + Ifges(5,5) * t621;
t587 = Ifges(6,4) * t191;
t101 = Ifges(6,2) * t720 + t338 * Ifges(6,6) + t587;
t643 = -t101 / 0.2e1;
t642 = t101 / 0.2e1;
t188 = Ifges(6,4) * t720;
t102 = t191 * Ifges(6,1) + t338 * Ifges(6,5) + t188;
t641 = -t102 / 0.2e1;
t640 = t102 / 0.2e1;
t639 = -t741 / 0.2e1;
t637 = -t118 / 0.2e1;
t635 = Ifges(4,1) * t623 + Ifges(4,4) * t622 + Ifges(4,5) * t609;
t174 = t287 * Ifges(5,1) + Ifges(5,4) * t513 + Ifges(5,5) * t344;
t633 = t174 / 0.2e1;
t632 = -t720 / 0.2e1;
t630 = -t191 / 0.2e1;
t620 = -t513 / 0.2e1;
t619 = -t287 / 0.2e1;
t617 = -t332 / 0.2e1;
t615 = -t338 / 0.2e1;
t613 = -t344 / 0.2e1;
t612 = t344 / 0.2e1;
t610 = t345 / 0.2e1;
t606 = cos(qJ(1));
t605 = sin(qJ(1));
t604 = mrSges(7,3) * t16;
t603 = mrSges(7,3) * t17;
t600 = pkin(5) * t191;
t444 = t454 * pkin(9);
t595 = mrSges(6,3) * t720;
t594 = mrSges(6,3) * t191;
t592 = Ifges(3,4) * t459;
t591 = Ifges(4,4) * t454;
t590 = Ifges(4,4) * t458;
t589 = Ifges(5,4) * t448;
t588 = Ifges(5,4) * t450;
t502 = t577 * t606;
t383 = t455 * t605 - t459 * t502;
t570 = t383 * t437;
t569 = t383 * t438;
t568 = t383 * t441;
t567 = t383 * t442;
t368 = qJD(2) * t479;
t370 = t388 * qJD(2);
t184 = -t357 * t549 + t358 * t548 + t454 * t368 + t458 * t370;
t552 = qJD(2) * t455;
t171 = (qJ(4) * t552 - qJD(4) * t459) * t449 + t184;
t551 = qJD(2) * t459;
t527 = t449 * t551;
t311 = qJD(3) * t382 + t454 * t527;
t312 = -qJD(3) * t381 + t458 * t527;
t371 = t555 * qJD(2);
t180 = t311 * pkin(3) - t312 * qJ(4) - t382 * qJD(4) + t371;
t95 = t450 * t171 + t448 * t180;
t384 = t455 * t502 + t459 * t605;
t533 = t449 * t606;
t315 = t384 * t458 - t454 * t533;
t559 = (-t315 * t437 + t569) * mrSges(7,1) + (-t315 * t438 - t570) * mrSges(7,2);
t501 = t577 * t605;
t386 = -t455 * t501 + t459 * t606;
t532 = t449 * t605;
t319 = t386 * t458 + t454 * t532;
t385 = t455 * t606 + t459 * t501;
t232 = -t319 * t437 + t385 * t438;
t233 = t319 * t438 + t385 * t437;
t558 = t232 * mrSges(7,1) - t233 * mrSges(7,2);
t557 = (-t382 * t437 - t438 * t563) * mrSges(7,1) + (-t382 * t438 + t437 * t563) * mrSges(7,2);
t401 = pkin(4) * t566 + t444;
t554 = t606 * pkin(1) + pkin(8) * t532;
t538 = Ifges(4,5) * t244 - Ifges(4,6) * t245 + Ifges(4,3) * t361;
t537 = Ifges(3,5) * t373 - Ifges(3,6) * t372 + Ifges(3,3) * t425;
t536 = t386 * pkin(2) + t554;
t528 = t449 * t552;
t524 = t564 / 0.2e1;
t44 = -t84 * mrSges(6,1) + t83 * mrSges(6,2);
t11 = -t32 * mrSges(7,1) + t31 * mrSges(7,2);
t516 = t548 / 0.2e1;
t121 = -t196 * mrSges(5,1) + t197 * mrSges(5,2);
t68 = t457 * t125 - t139 * t453;
t94 = -t171 * t448 + t450 * t180;
t265 = -t454 * t357 + t358 * t458;
t509 = mrSges(3,3) * t529;
t503 = -pkin(1) * t605 + pkin(8) * t533;
t251 = pkin(3) * t563 - t265;
t500 = mrSges(4,1) * t381 + mrSges(4,2) * t382;
t494 = Ifges(4,1) * t458 - t591;
t493 = Ifges(5,1) * t450 - t589;
t492 = -Ifges(4,2) * t454 + t590;
t491 = -Ifges(5,2) * t448 + t588;
t490 = Ifges(4,5) * t458 - Ifges(4,6) * t454;
t489 = Ifges(5,5) * t450 - Ifges(5,6) * t448;
t213 = t309 * t453 + t310 * t457;
t54 = pkin(5) * t381 - pkin(11) * t213 + t68;
t212 = t309 * t457 - t310 * t453;
t56 = pkin(11) * t212 + t69;
t22 = -t452 * t56 + t456 * t54;
t23 = t452 * t54 + t456 * t56;
t147 = t212 * t456 - t213 * t452;
t148 = t212 * t452 + t213 * t456;
t246 = -t319 * t441 + t385 * t442;
t299 = -t399 * t452 - t456 * t486;
t300 = t399 * t456 - t452 * t486;
t481 = -t384 * pkin(2) + t503;
t185 = -t357 * t548 - t358 * t549 + t368 * t458 - t454 * t370;
t480 = t719 + t8;
t477 = t327 * (mrSges(4,1) * t454 + mrSges(4,2) * t458);
t276 = t312 * t450 + t448 * t528;
t75 = pkin(4) * t311 - pkin(10) * t276 + t94;
t275 = -t312 * t448 + t450 * t528;
t85 = pkin(10) * t275 + t95;
t20 = t125 * t544 - t139 * t545 + t453 * t75 + t457 * t85;
t314 = t384 * t454 + t458 * t533;
t318 = t386 * t454 - t458 * t532;
t473 = -g(1) * t318 - g(2) * t314 - g(3) * t381;
t194 = -pkin(4) * t309 + t251;
t170 = -pkin(4) * t513 + t208;
t179 = -pkin(3) * t528 - t185;
t21 = -qJD(5) * t69 - t453 * t85 + t457 * t75;
t464 = Ifges(3,6) * t577 + (t459 * Ifges(3,2) + t593) * t449;
t78 = -pkin(4) * t196 + t112;
t463 = t449 * t484 * (Ifges(3,5) * t459 - Ifges(3,6) * t455);
t140 = -pkin(4) * t275 + t179;
t420 = Ifges(3,4) * t529;
t387 = (-mrSges(3,1) * t459 + mrSges(3,2) * t455) * t449;
t369 = t555 * qJD(1);
t365 = -mrSges(3,2) * t484 + t509;
t352 = pkin(5) * t486 - t439;
t342 = -pkin(9) * t565 + t393;
t337 = Ifges(4,4) * t344;
t325 = Ifges(3,1) * t530 + Ifges(3,5) * t484 + t420;
t324 = Ifges(3,6) * qJD(2) + qJD(1) * t464;
t303 = pkin(5) * t362 + t401;
t298 = mrSges(4,1) * t411 - mrSges(4,3) * t345;
t297 = -mrSges(4,2) * t411 - mrSges(4,3) * t344;
t247 = t319 * t442 + t385 * t441;
t228 = t345 * Ifges(4,1) - t337 + t700;
t226 = t345 * Ifges(4,5) - t344 * Ifges(4,6) + t411 * Ifges(4,3);
t215 = mrSges(5,1) * t344 - mrSges(5,3) * t287;
t214 = -mrSges(5,2) * t344 + mrSges(5,3) * t513;
t205 = -qJD(6) * t300 + t379 * t452 - t380 * t456;
t204 = qJD(6) * t299 - t379 * t456 - t380 * t452;
t200 = -mrSges(4,2) * t361 - mrSges(4,3) * t245;
t199 = mrSges(4,1) * t361 - mrSges(4,3) * t244;
t198 = -mrSges(5,1) * t513 + mrSges(5,2) * t287;
t166 = mrSges(6,1) * t338 - t594;
t165 = -mrSges(6,2) * t338 + t595;
t164 = mrSges(4,1) * t245 + mrSges(4,2) * t244;
t161 = t242 * t452 + t243 * t456;
t160 = t242 * t456 - t243 * t452;
t153 = t244 * Ifges(4,4) - t245 * Ifges(4,2) + t361 * Ifges(4,6);
t142 = mrSges(5,1) * t245 - mrSges(5,3) * t197;
t141 = -mrSges(5,2) * t245 + mrSges(5,3) * t196;
t135 = -pkin(5) * t212 + t194;
t129 = -qJD(5) * t213 + t275 * t457 - t276 * t453;
t128 = qJD(5) * t212 + t275 * t453 + t276 * t457;
t120 = -mrSges(6,1) * t720 + mrSges(6,2) * t191;
t103 = -pkin(5) * t720 + t170;
t98 = mrSges(7,1) * t332 - mrSges(7,3) * t118;
t97 = -mrSges(7,2) * t332 + mrSges(7,3) * t741;
t70 = -pkin(5) * t129 + t140;
t67 = -mrSges(6,2) * t236 + mrSges(6,3) * t84;
t66 = mrSges(6,1) * t236 - mrSges(6,3) * t83;
t63 = -mrSges(7,1) * t741 + mrSges(7,2) * t118;
t49 = -qJD(6) * t148 - t128 * t452 + t129 * t456;
t48 = qJD(6) * t147 + t128 * t456 + t129 * t452;
t42 = -pkin(5) * t84 + t78;
t25 = -mrSges(7,2) * t225 + mrSges(7,3) * t32;
t24 = mrSges(7,1) * t225 - mrSges(7,3) * t31;
t19 = t456 * t46 - t579;
t18 = -t452 * t46 - t578;
t15 = pkin(11) * t129 + t20;
t14 = pkin(5) * t311 - pkin(11) * t128 + t21;
t5 = -qJD(6) * t23 + t14 * t456 - t15 * t452;
t4 = qJD(6) * t22 + t14 * t452 + t15 * t456;
t1 = [(Ifges(7,4) * t48 + Ifges(7,2) * t49) * t638 + (Ifges(7,4) * t148 + Ifges(7,2) * t147) * t654 + t276 * t633 + t382 * t635 - t324 * t528 / 0.2e1 + t237 * (mrSges(4,1) * t528 - mrSges(4,3) * t312) - t387 * t540 + (-m(3) * t554 - mrSges(2,1) * t606 - t386 * mrSges(3,1) - t247 * mrSges(6,1) - t233 * mrSges(7,1) + mrSges(2,2) * t605 - t246 * mrSges(6,2) - t232 * mrSges(7,2) - mrSges(3,3) * t532 + t714 * t536 + t715 * (pkin(9) * t385 + t536) + t661 * t319 - t662 * t385 + t726 * t318) * g(2) + (Ifges(5,4) * t310 + Ifges(5,2) * t309) * t628 + t513 * (Ifges(5,4) * t276 + Ifges(5,2) * t275) / 0.2e1 + t128 * t640 + t129 * t642 + t310 * t644 + t309 * t645 + t48 * t648 + t49 * t650 + (Ifges(5,1) * t310 + Ifges(5,4) * t309) * t627 + t287 * (Ifges(5,1) * t276 + Ifges(5,4) * t275) / 0.2e1 + (t660 * qJD(1) * (-Ifges(3,2) * t455 + t592) + t449 * t325) * t551 / 0.2e1 + (-t290 * t577 - t373 * t602 - t425 * t555) * mrSges(3,2) + m(3) * (pkin(1) ^ 2 * qJDD(1) * t660 + t290 * t555 + t291 * t388 + t369 * t370) + t411 * (Ifges(4,5) * t312 + Ifges(4,3) * t528) / 0.2e1 + (Ifges(4,5) * t382 - Ifges(4,3) * t563) * t609 - t372 * t464 / 0.2e1 - t674 * t543 + (Ifges(6,4) * t128 + Ifges(6,2) * t129) * t631 + (Ifges(6,4) * t213 + Ifges(6,2) * t212) * t646 + (t727 / 0.2e1 + t728 / 0.2e1 + t716 - Ifges(4,2) * t613 - Ifges(4,4) * t610 + Ifges(5,3) * t612 + t729) * t311 + (Ifges(7,1) * t48 + Ifges(7,4) * t49) * t636 + (Ifges(7,1) * t148 + Ifges(7,4) * t147) * t655 + (Ifges(5,5) * t627 + Ifges(5,6) * t628 + Ifges(6,6) * t646 + Ifges(6,5) * t647 + Ifges(5,3) * t621 - Ifges(4,2) * t622 - Ifges(4,4) * t623 + Ifges(6,3) * t624 + Ifges(7,3) * t626 - Ifges(4,6) * t609 + t675 / 0.2e1 - t133 * mrSges(4,3) + Ifges(7,6) * t654 + Ifges(7,5) * t655 - t61 * mrSges(5,2) + t60 * mrSges(5,1) - t153 / 0.2e1 + t718 + t719) * t381 + (Ifges(7,5) * t48 + Ifges(7,6) * t49) * t616 + (Ifges(7,5) * t148 + Ifges(7,6) * t147) * t626 + (t463 / 0.2e1 + t226 * t524) * qJD(2) + (t133 * t563 - t238 * t528 + t312 * t327) * mrSges(4,2) - t538 * t563 / 0.2e1 + t373 * (Ifges(3,5) * t577 + (t455 * Ifges(3,1) + t592) * t449) / 0.2e1 + (Ifges(5,5) * t310 + Ifges(5,6) * t309) * t621 + (Ifges(5,5) * t276 + Ifges(5,6) * t275) * t612 + (t147 * t2 - t148 * t3 - t16 * t48 + t17 * t49) * mrSges(7,3) + (t12 * t212 - t128 * t52 + t129 * t53 - t13 * t213) * mrSges(6,3) + (-m(3) * t366 - t666) * t371 + (-t137 * t276 + t138 * t275 + t309 * t61 - t310 * t60) * mrSges(5,3) + m(7) * (t103 * t70 + t135 * t42 + t16 * t5 + t17 * t4 + t2 * t23 + t22 * t3) + m(6) * (t12 * t69 + t13 * t68 + t140 * t170 + t194 * t78 + t20 * t53 + t21 * t52) + m(5) * (t112 * t251 + t137 * t94 + t138 * t95 + t162 * t60 + t163 * t61 + t179 * t208) + (Ifges(4,4) * t312 + Ifges(4,6) * t528) * t613 + (Ifges(4,4) * t382 - Ifges(4,6) * t563) * t622 + t134 * (-mrSges(4,1) * t563 - t382 * mrSges(4,3)) + t194 * t44 + t179 * t198 + (Ifges(6,1) * t128 + Ifges(6,4) * t129) * t629 + (Ifges(6,1) * t213 + Ifges(6,4) * t212) * t647 + t262 * t500 + t577 * t537 / 0.2e1 + t425 * (Ifges(3,3) * t577 + (Ifges(3,5) * t455 + Ifges(3,6) * t459) * t449) / 0.2e1 + t213 * t652 + t212 * t653 + t148 * t656 + t147 * t658 + t20 * t165 + t21 * t166 + t170 * (-mrSges(6,1) * t129 + mrSges(6,2) * t128) + t162 * t142 + t163 * t141 + t42 * (-mrSges(7,1) * t147 + mrSges(7,2) * t148) + t140 * t120 + t135 * t11 + Ifges(2,3) * qJDD(1) + (Ifges(4,1) * t382 - Ifges(4,5) * t563) * t623 + (Ifges(4,1) * t312 + Ifges(4,5) * t528) * t610 + t103 * (-mrSges(7,1) * t49 + mrSges(7,2) * t48) + t4 * t97 + t5 * t98 + t70 * t63 + t68 * t66 + t69 * t67 + t22 * t24 + t23 * t25 + m(4) * (t133 * t266 + t134 * t265 + t184 * t238 + t185 * t237 + t262 * t356) + (t290 * t563 - t291 * t564 - t366 * t527 - t369 * t528 - t372 * t555 - t373 * t388) * mrSges(3,3) + (t291 * t577 - t372 * t602 + t388 * t425) * mrSges(3,1) + (Ifges(3,4) * t373 - Ifges(3,2) * t372 + Ifges(3,6) * t425) * t563 / 0.2e1 + (Ifges(3,1) * t373 - Ifges(3,4) * t372 + Ifges(3,5) * t425) * t524 + (Ifges(6,5) * t128 + Ifges(6,6) * t129) * t614 + (Ifges(6,5) * t213 + Ifges(6,6) * t212) * t624 + t78 * (-mrSges(6,1) * t212 + mrSges(6,2) * t213) + t95 * t214 + t94 * t215 + (-m(3) * t503 + mrSges(2,1) * t605 + t384 * mrSges(3,1) + t568 * mrSges(6,1) + t570 * mrSges(7,1) + mrSges(2,2) * t606 + t567 * mrSges(6,2) + t569 * mrSges(7,2) - mrSges(3,3) * t533 + t714 * t481 + t715 * (-t383 * pkin(9) + t481) - (t661 + t744) * t315 + t662 * t383 - t726 * t314) * g(1) + t251 * t121 + t265 * t199 + t266 * t200 + t275 * t173 / 0.2e1 + t208 * (-mrSges(5,1) * t275 + mrSges(5,2) * t276) + t184 * t297 + t185 * t298 + t112 * (-mrSges(5,1) * t309 + mrSges(5,2) * t310) + t312 * t228 / 0.2e1 + t356 * t164 + t370 * t365; (-Ifges(5,5) * t458 + t454 * t493) * t627 + (-Ifges(5,6) * t458 + t454 * t491) * t628 + (Ifges(6,1) * t241 + Ifges(6,4) * t240 + Ifges(6,5) * t507) * t630 + (Ifges(6,4) * t241 + Ifges(6,2) * t240 + Ifges(6,6) * t507) * t632 + t330 * t634 + t454 * t635 + (Ifges(7,1) * t159 + Ifges(7,4) * t158 + Ifges(7,5) * t507) * t637 + (-t237 * t277 - t238 * t278 - pkin(2) * t262 + ((-t237 * t458 - t238 * t454) * qJD(3) + t678) * pkin(9)) * m(4) + t678 * mrSges(4,3) - t477 * t529 + (-t237 * (mrSges(4,1) * t455 - mrSges(4,3) * t560) - t238 * (-mrSges(4,3) * t454 * t459 - mrSges(4,2) * t455)) * t553 + (t509 - t365) * t366 + t696 * t166 + (t12 * t210 + t13 * t209 + t170 * t686 + t401 * t78 + t52 * t696 + t53 * t697) * m(6) + t697 * t165 + t228 * t516 - ((-Ifges(3,2) * t530 + t458 * t228 + t454 * t673 + t325 + t420) * t459 + t411 * (Ifges(4,3) * t455 + t459 * t490) + t345 * (Ifges(4,5) * t455 + t459 * t494) + t455 * t226) * t553 / 0.2e1 + (Ifges(7,4) * t159 + Ifges(7,2) * t158 + Ifges(7,6) * t507) * t639 + t281 * t640 + t241 * t641 + t282 * t642 + t240 * t643 + t562 * t644 + t151 * t648 + t159 * t649 + t152 * t650 + t158 * t651 + t227 * t507 / 0.2e1 - t138 * (-mrSges(5,2) * t507 + t330 * mrSges(5,3)) - t137 * (mrSges(5,1) * t507 - t331 * mrSges(5,3)) + (t344 * (Ifges(4,6) * t455 + t459 * t492) + t455 * t324) * t553 / 0.2e1 + (t287 * (Ifges(5,5) * t454 + t458 * t493) + t513 * (Ifges(5,6) * t454 + t458 * t491) + t411 * t490 + t345 * t494) * qJD(3) / 0.2e1 + (Ifges(5,5) * t331 + Ifges(5,6) * t330 + Ifges(5,3) * t507) * t613 + (Ifges(6,5) * t241 + Ifges(6,6) * t240 + Ifges(6,3) * t507) * t615 + (Ifges(7,5) * t159 + Ifges(7,6) * t158 + Ifges(7,3) * t507) * t617 + (Ifges(5,1) * t331 + Ifges(5,4) * t330 + Ifges(5,5) * t507) * t619 + (Ifges(5,4) * t331 + Ifges(5,2) * t330 + Ifges(5,6) * t507) * t620 + (-Ifges(5,3) * t458 + t454 * t489) * t621 + (Ifges(4,2) * t458 + t591) * t622 + (Ifges(4,1) * t454 + t590) * t623 + (Ifges(7,5) * t271 + Ifges(7,6) * t270 - Ifges(7,3) * t458) * t626 + (Ifges(4,5) * t454 + Ifges(4,6) * t458) * t609 - t675 * t458 / 0.2e1 + t707 * t98 + (t103 * t689 + t107 * t3 + t108 * t2 + t16 * t707 + t17 * t708 + t303 * t42) * m(7) + t708 * t97 + (Ifges(7,1) * t151 + Ifges(7,4) * t152) * t636 + (Ifges(6,4) * t281 + Ifges(6,2) * t282) * t631 + (t450 * t516 - t331 / 0.2e1) * t174 + (t121 - t199) * t444 + (Ifges(6,5) * t281 + Ifges(6,6) * t282) * t614 + (Ifges(6,1) * t281 + Ifges(6,4) * t282) * t629 + t60 * (-mrSges(5,1) * t458 - mrSges(5,3) * t562) + (-t463 / 0.2e1 + t674 * qJD(1)) * qJD(1) + (Ifges(5,3) * t454 + t458 * t489) * t550 / 0.2e1 - t492 * t550 / 0.2e1 + (t510 + t666) * t369 + (Ifges(7,4) * t151 + Ifges(7,2) * t152) * t638 - t87 * t566 / 0.2e1 + t61 * (mrSges(5,2) * t458 - mrSges(5,3) * t566) + (-t440 - t277) * t298 + t689 * t63 + (mrSges(7,2) * t681 + mrSges(7,3) * t690) * t17 + (-mrSges(7,1) * t690 + mrSges(7,2) * t691) * t103 + (-mrSges(7,1) * t681 - mrSges(7,3) * t691) * t16 + t679 * t497 + t680 * t548 + (t137 * (mrSges(5,1) * t454 - mrSges(5,3) * t561) + t477 + t138 * (-mrSges(5,2) * t454 - mrSges(5,3) * t565)) * qJD(3) + (-t541 - t278) * t297 + t209 * t66 + t210 * t67 - t262 * t499 - t363 * t652 - t362 * t653 + (Ifges(7,4) * t271 + Ifges(7,2) * t270 - Ifges(7,6) * t458) * t654 + (Ifges(7,1) * t271 + Ifges(7,4) * t270 - Ifges(7,5) * t458) * t655 + t271 * t656 + t270 * t658 - pkin(2) * t164 + (Ifges(7,5) * t151 + Ifges(7,6) * t152) * t616 + t108 * t25 + t107 * t24 + t2 * (mrSges(7,2) * t458 + mrSges(7,3) * t270) + t3 * (-mrSges(7,1) * t458 - mrSges(7,3) * t271) + t458 * t153 / 0.2e1 + t12 * (mrSges(6,2) * t458 - mrSges(6,3) * t362) + (-Ifges(6,4) * t363 - Ifges(6,2) * t362 - Ifges(6,6) * t458) * t646 + t13 * (-mrSges(6,1) * t458 + mrSges(6,3) * t363) + (-Ifges(6,1) * t363 - Ifges(6,4) * t362 - Ifges(6,5) * t458) * t647 + (-Ifges(6,5) * t363 - Ifges(6,6) * t362 - Ifges(6,3) * t458) * t624 + t78 * (mrSges(6,1) * t362 - mrSges(6,2) * t363) + (-mrSges(6,1) * t681 + mrSges(6,3) * t684) * t52 + (mrSges(6,1) * t685 - mrSges(6,2) * t684) * t170 + (mrSges(6,2) * t681 - mrSges(6,3) * t685) * t53 + t686 * t120 + t687 * t214 + (pkin(9) * t679 + t137 * t688 + t138 * t687 - t208 * t253 + t342 * t60 + t343 * t61) * m(5) + t688 * t215 + t458 * pkin(9) * t200 + t731 * t198 + (t387 - t672 * t556 + (t742 * t459 + (-m(6) * t601 - m(7) * t402 + t730 - t738) * t455) * t449) * g(3) + t537 + (t383 * t723 + t384 * t664) * g(2) + (t385 * t723 + t386 * t664) * g(1) + t42 * (-mrSges(7,1) * t270 + mrSges(7,2) * t271) - t290 * mrSges(3,2) + t291 * mrSges(3,1) + t303 * t11 - t208 * (-mrSges(5,1) * t330 + mrSges(5,2) * t331) + t342 * t142 + t343 * t141 + (-t583 + t716) * t549 + t401 * t44; (Ifges(7,5) * t300 + Ifges(7,6) * t299) * t626 + (Ifges(5,1) * t448 + t588) * t627 + (Ifges(5,2) * t450 + t589) * t628 + t571 * t633 + t676 * t63 + (-t137 * t571 - t138 * t572 + t677) * mrSges(5,3) + (-t137 * t167 - t138 * t168 - t208 * t238 - pkin(3) * t112 + (-t137 * t448 + t138 * t450) * qJD(4) + t677 * qJ(4)) * m(5) + (-t489 * t613 - t493 * t619 - t491 * t620 + t208 * t497 + t327 * mrSges(4,2) + t700 / 0.2e1 + t680) * t344 + t694 * t165 + (t12 * t323 + t13 * t322 - t170 * t195 - t439 * t78 + t52 * t695 + t53 * t694) * m(6) + t695 * t166 + (Ifges(7,4) * t636 + Ifges(7,2) * t638 + Ifges(7,6) * t616 + t603 + t650) * t205 + (t381 * t671 + t382 * t670 + t500) * g(3) + (Ifges(7,1) * t636 + Ifges(7,4) * t638 + Ifges(7,5) * t616 - t604 + t648) * t204 + (t318 * t669 + t319 * t726) * g(1) + (t314 * t669 + t315 * t726) * g(2) + (Ifges(6,4) * t399 - Ifges(6,2) * t486) * t646 + (Ifges(6,1) * t399 - Ifges(6,4) * t486) * t647 + (Ifges(6,5) * t399 - Ifges(6,6) * t486) * t624 + t78 * (mrSges(6,1) * t486 + mrSges(6,2) * t399) - t486 * t653 + (-t12 * t486 - t13 * t399 + t52 * t682 - t53 * t683) * mrSges(6,3) + (t16 * t161 - t160 * t17 + t2 * t299 - t3 * t300) * mrSges(7,3) + (Ifges(6,1) * t243 + Ifges(6,4) * t242) * t630 - t379 * t640 + t243 * t641 - t380 * t642 + t242 * t643 + t448 * t644 + t450 * t645 + t161 * t649 + t160 * t651 + t399 * t652 + (Ifges(6,4) * t243 + Ifges(6,2) * t242) * t632 + (Ifges(5,5) * t448 + Ifges(5,6) * t450) * t621 + t227 * t610 - (-Ifges(4,1) * t344 - t582 + t673) * t345 / 0.2e1 + (-Ifges(6,4) * t379 - Ifges(6,2) * t380) * t631 + (-Ifges(6,5) * t379 - Ifges(6,6) * t380) * t614 + (-Ifges(6,1) * t379 - Ifges(6,4) * t380) * t629 + t705 * t98 + (t103 * t676 + t16 * t705 + t17 * t706 + t182 * t3 + t183 * t2 + t352 * t42) * m(7) + t706 * t97 + (t228 - t337) * t612 + (t141 * t450 - t142 * t448) * qJ(4) + (Ifges(5,5) * t619 + Ifges(6,5) * t630 + Ifges(7,5) * t637 - Ifges(4,2) * t612 + Ifges(5,6) * t620 + Ifges(6,6) * t632 + Ifges(7,6) * t639 + Ifges(5,3) * t613 + Ifges(6,3) * t615 + Ifges(7,3) * t617 - t729) * t345 + (-t198 + t298) * t238 + (Ifges(6,5) * t243 + Ifges(6,6) * t242) * t615 + (Ifges(7,5) * t161 + Ifges(7,6) * t160) * t617 + (-t547 - t167) * t215 + (Ifges(7,4) * t161 + Ifges(7,2) * t160) * t639 + (t546 - t168) * t214 + (Ifges(7,1) * t161 + Ifges(7,4) * t160) * t637 - t195 * t120 + t112 * t498 + t538 + (Ifges(7,4) * t300 + Ifges(7,2) * t299) * t654 + (Ifges(7,1) * t300 + Ifges(7,4) * t299) * t655 + t300 * t656 + t299 * t658 + t182 * t24 + t183 * t25 - t133 * mrSges(4,2) + t134 * mrSges(4,1) - pkin(3) * t121 + ((-t161 + t204) * mrSges(7,2) + (t160 - t205) * mrSges(7,1)) * t103 - t439 * t44 + (mrSges(6,1) * t683 - mrSges(6,2) * t682) * t170 - t237 * t297 + t42 * (-mrSges(7,1) * t299 + mrSges(7,2) * t300) + t322 * t66 + t323 * t67 + t352 * t11; t118 * t98 - t741 * t97 - t720 * t165 + t191 * t166 - t513 * t214 + t287 * t215 + t11 + t121 + t44 + (t118 * t16 - t17 * t741 + t42 + t473) * m(7) + (t191 * t52 - t53 * t720 + t473 + t78) * m(6) + (t137 * t287 - t138 * t513 + t112 + t473) * m(5); t101 * t629 + (Ifges(6,1) * t720 - t587) * t630 + (-(-t382 * t442 + t441 * t563) * mrSges(6,2) - t557 + t698 * (-t382 * t441 - t442 * t563)) * g(3) + (-(-t315 * t442 - t568) * mrSges(6,2) - t559 + t698 * (-t315 * t441 + t567)) * g(2) + (mrSges(6,2) * t247 + t246 * t698 - t558) * g(1) - (mrSges(7,1) * t103 + Ifges(7,4) * t637 + Ifges(7,2) * t639 + Ifges(7,6) * t617 - t603 + t651) * t118 + (-mrSges(7,2) * t103 + Ifges(7,1) * t637 + Ifges(7,4) * t639 + Ifges(7,5) * t617 + t604 + t649) * t741 + (Ifges(6,5) * t720 - Ifges(6,6) * t191) * t615 + (-Ifges(6,2) * t191 + t102 + t188) * t632 + t480 - t63 * t600 - m(7) * (t103 * t600 + t16 * t18 + t17 * t19) + (t594 + t166) * t53 - t170 * (mrSges(6,1) * t191 + mrSges(6,2) * t720) + t718 + t37 + (t2 * t452 + t3 * t456 + (-t16 * t452 + t17 * t456) * qJD(6)) * t657 - t19 * t97 - t18 * t98 + ((-t452 * t98 + t456 * t97) * qJD(6) + t24 * t456 + t25 * t452) * pkin(5) + (t595 - t165) * t52; -t103 * (mrSges(7,1) * t118 + mrSges(7,2) * t741) + (Ifges(7,1) * t741 - t586) * t637 + t58 * t636 + (Ifges(7,5) * t741 - Ifges(7,6) * t118) * t617 - t16 * t97 + t17 * t98 - g(1) * t558 - g(2) * t559 - g(3) * t557 + (t118 * t17 + t16 * t741) * mrSges(7,3) + t480 + (-Ifges(7,2) * t118 + t113 + t59) * t639;];
tau  = t1;
