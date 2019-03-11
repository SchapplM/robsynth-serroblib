% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:27:51
% EndTime: 2019-03-09 13:29:07
% DurationCPUTime: 43.73s
% Computational Cost: add. (31907->1013), mult. (86883->1389), div. (0->0), fcn. (71937->16), ass. (0->467)
t394 = cos(pkin(6));
t398 = sin(qJ(2));
t546 = t394 * t398;
t379 = pkin(1) * t546;
t393 = sin(pkin(6));
t403 = cos(qJ(2));
t548 = t393 * t403;
t597 = pkin(8) + qJ(3);
t301 = (t548 * t597 + t379) * qJD(1);
t563 = sin(pkin(12));
t289 = t563 * t301;
t606 = pkin(1) * t394;
t381 = t403 * t606;
t373 = qJD(1) * t381;
t487 = t597 * t398;
t467 = t393 * t487;
t300 = -qJD(1) * t467 + t373;
t564 = cos(pkin(12));
t199 = t300 * t564 - t289;
t529 = qJD(1) * t393;
t460 = t563 * t529;
t461 = t564 * t529;
t314 = -t398 * t460 + t403 * t461;
t315 = -t398 * t461 - t403 * t460;
t493 = t398 * t529;
t472 = pkin(2) * t493;
t232 = -pkin(3) * t315 - pkin(9) * t314 + t472;
t397 = sin(qJ(4));
t402 = cos(qJ(4));
t134 = -t199 * t397 + t232 * t402;
t494 = t563 * pkin(2);
t383 = t494 + pkin(9);
t596 = pkin(10) + t383;
t482 = qJD(4) * t596;
t555 = t314 * t402;
t714 = pkin(4) * t315 + pkin(10) * t555 - t402 * t482 - t134;
t135 = t199 * t402 + t232 * t397;
t556 = t314 * t397;
t713 = -pkin(10) * t556 + t397 * t482 + t135;
t308 = qJD(4) - t314;
t306 = qJD(5) + t308;
t615 = t306 / 0.2e1;
t376 = qJD(1) * t394 + qJD(2);
t272 = t315 * t397 + t376 * t402;
t396 = sin(qJ(5));
t401 = cos(qJ(5));
t437 = t315 * t402 - t376 * t397;
t439 = t272 * t396 - t401 * t437;
t622 = t439 / 0.2e1;
t712 = Ifges(6,4) * t622 + Ifges(6,6) * t615;
t476 = t272 * t401 + t396 * t437;
t166 = qJD(6) - t476;
t627 = -t166 / 0.2e1;
t395 = sin(qJ(6));
t400 = cos(qJ(6));
t144 = t306 * t395 + t400 * t439;
t633 = -t144 / 0.2e1;
t143 = t306 * t400 - t395 * t439;
t635 = -t143 / 0.2e1;
t711 = Ifges(7,5) * t633 + Ifges(7,6) * t635 + Ifges(7,3) * t627;
t706 = mrSges(6,2) - mrSges(7,3);
t624 = t476 / 0.2e1;
t710 = Ifges(6,2) * t624;
t455 = -mrSges(7,1) * t400 + mrSges(7,2) * t395;
t705 = mrSges(6,1) - t455;
t347 = t396 * t397 - t401 * t402;
t225 = t347 * t314;
t669 = qJD(4) + qJD(5);
t292 = t669 * t347;
t709 = t225 - t292;
t348 = t396 * t402 + t397 * t401;
t224 = t348 * t314;
t293 = t669 * t348;
t678 = t293 - t224;
t656 = m(7) * pkin(5);
t708 = t656 + t705;
t655 = m(7) * pkin(11);
t473 = -t655 + t706;
t343 = t596 * t397;
t344 = t596 * t402;
t436 = -t343 * t401 - t344 * t396;
t687 = qJD(5) * t436 + t396 * t714 - t401 * t713;
t477 = t564 * t301;
t198 = t300 * t563 + t477;
t525 = qJD(4) * t397;
t672 = -t198 + (t525 - t556) * pkin(4);
t276 = pkin(2) * t376 + t300;
t183 = t276 * t564 - t289;
t175 = -pkin(3) * t376 - t183;
t141 = -pkin(4) * t272 + t175;
t184 = t276 * t563 + t477;
t176 = pkin(9) * t376 + t184;
t389 = pkin(2) * t403 + pkin(1);
t335 = -t389 * t529 + qJD(3);
t195 = -pkin(3) * t314 + pkin(9) * t315 + t335;
t128 = t176 * t402 + t195 * t397;
t110 = pkin(10) * t272 + t128;
t539 = t401 * t110;
t127 = -t176 * t397 + t195 * t402;
t109 = pkin(10) * t437 + t127;
t97 = pkin(4) * t308 + t109;
t55 = t396 * t97 + t539;
t53 = pkin(11) * t306 + t55;
t86 = -pkin(5) * t476 - pkin(11) * t439 + t141;
t25 = -t395 * t53 + t400 * t86;
t26 = t395 * t86 + t400 * t53;
t707 = t711 + t712;
t702 = -mrSges(6,1) * t141 - mrSges(7,1) * t25 + mrSges(7,2) * t26 + mrSges(6,3) * t55 + t707 + t710;
t704 = pkin(11) * t315 + t687;
t703 = pkin(5) * t678 - pkin(11) * t709 + t672;
t519 = qJD(1) * qJD(2);
t330 = (qJDD(1) * t403 - t398 * t519) * t393;
t331 = (qJDD(1) * t398 + t403 * t519) * t393;
t252 = t330 * t563 + t331 * t564;
t517 = qJDD(1) * t394;
t375 = qJDD(2) + t517;
t154 = qJD(4) * t272 + t252 * t402 + t375 * t397;
t251 = t330 * t564 - t331 * t563;
t250 = qJDD(4) - t251;
t530 = pkin(8) * t548 + t379;
t329 = t530 * qJD(2);
t508 = pkin(1) * t517;
t371 = t403 * t508;
t550 = t393 * t398;
t490 = qJD(3) * t550;
t518 = qJDD(1) * t393;
t507 = pkin(8) * t518;
t180 = -t398 * t507 + pkin(2) * t375 - qJ(3) * t331 + t371 + (-t329 - t490) * qJD(1);
t515 = qJD(2) * t606;
t468 = qJD(1) * t515;
t496 = t398 * t508 + (t468 + t507) * t403;
t526 = qJD(3) * t403;
t527 = qJD(2) * t398;
t188 = qJ(3) * t330 + (-pkin(8) * t527 + t526) * t529 + t496;
t126 = t180 * t563 + t188 * t564;
t120 = pkin(9) * t375 + t126;
t296 = -pkin(1) * t518 - pkin(2) * t330 + qJDD(3);
t140 = -pkin(3) * t251 - pkin(9) * t252 + t296;
t51 = -qJD(4) * t128 - t120 * t397 + t140 * t402;
t34 = pkin(4) * t250 - pkin(10) * t154 + t51;
t155 = qJD(4) * t437 - t252 * t397 + t375 * t402;
t524 = qJD(4) * t402;
t50 = t120 * t402 + t140 * t397 - t176 * t525 + t195 * t524;
t41 = pkin(10) * t155 + t50;
t522 = qJD(5) * t401;
t523 = qJD(5) * t396;
t10 = -t110 * t523 + t34 * t396 + t401 * t41 + t522 * t97;
t242 = qJDD(5) + t250;
t84 = qJD(5) * t476 + t154 * t401 + t155 * t396;
t47 = qJD(6) * t143 + t242 * t395 + t400 * t84;
t48 = -qJD(6) * t144 + t242 * t400 - t395 * t84;
t85 = -qJD(5) * t439 - t154 * t396 + t155 * t401;
t83 = qJDD(6) - t85;
t14 = Ifges(7,5) * t47 + Ifges(7,6) * t48 + Ifges(7,3) * t83;
t621 = t242 / 0.2e1;
t642 = t85 / 0.2e1;
t643 = t84 / 0.2e1;
t644 = t83 / 0.2e1;
t648 = t48 / 0.2e1;
t649 = t47 / 0.2e1;
t125 = t180 * t564 - t188 * t563;
t119 = -pkin(3) * t375 - t125;
t89 = -pkin(4) * t155 + t119;
t22 = -pkin(5) * t85 - pkin(11) * t84 + t89;
t7 = pkin(11) * t242 + t10;
t2 = qJD(6) * t25 + t22 * t395 + t400 * t7;
t3 = -qJD(6) * t26 + t22 * t400 - t395 * t7;
t662 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t701 = t662 + mrSges(6,1) * t89 - mrSges(6,3) * t10 + Ifges(7,5) * t649 + Ifges(7,6) * t648 + Ifges(7,3) * t644 + t14 / 0.2e1 + (-t621 - t242 / 0.2e1) * Ifges(6,6) + (-t642 - t85 / 0.2e1) * Ifges(6,2) + (-t643 - t84 / 0.2e1) * Ifges(6,4);
t691 = m(7) + m(6);
t692 = -m(5) - m(4);
t511 = t691 - t692;
t626 = t166 / 0.2e1;
t632 = t144 / 0.2e1;
t634 = t143 / 0.2e1;
t699 = Ifges(7,5) * t632 + Ifges(7,6) * t634 + Ifges(7,3) * t626 - t702 - t710 - t712;
t414 = t398 * t564 + t403 * t563;
t323 = t414 * t393;
t294 = -t323 * t397 + t394 * t402;
t463 = t294 * pkin(4);
t275 = -t343 * t396 + t344 * t401;
t686 = -qJD(5) * t275 + t396 * t713 + t401 * t714;
t146 = mrSges(6,1) * t306 - mrSges(6,3) * t439;
t95 = -mrSges(7,1) * t143 + mrSges(7,2) * t144;
t683 = t95 - t146;
t345 = t398 * t563 - t403 * t564;
t322 = t345 * t393;
t698 = -Ifges(5,5) * t437 + Ifges(6,5) * t439 + Ifges(5,6) * t272 + Ifges(6,6) * t476 + Ifges(5,3) * t308 + Ifges(6,3) * t306;
t392 = qJ(4) + qJ(5);
t390 = sin(t392);
t391 = cos(t392);
t457 = -mrSges(5,1) * t402 + mrSges(5,2) * t397;
t418 = m(5) * pkin(3) - t457;
t697 = -t390 * t473 + t391 * t708 + t418;
t404 = cos(qJ(1));
t538 = t403 * t404;
t399 = sin(qJ(1));
t543 = t398 * t399;
t696 = t394 * t538 - t543;
t531 = t414 * t394;
t262 = -t345 * t404 - t399 * t531;
t549 = t393 * t399;
t237 = -t262 * t397 + t402 * t549;
t616 = -t306 / 0.2e1;
t623 = -t439 / 0.2e1;
t625 = -t476 / 0.2e1;
t695 = -Ifges(6,4) * t623 - Ifges(6,2) * t625 - Ifges(6,6) * t616 + t702 + t711;
t650 = Ifges(6,1) * t643 + Ifges(6,4) * t642 + Ifges(6,5) * t621;
t629 = t154 / 0.2e1;
t628 = t155 / 0.2e1;
t620 = t250 / 0.2e1;
t21 = -mrSges(7,1) * t48 + mrSges(7,2) * t47;
t66 = mrSges(6,1) * t242 - mrSges(6,3) * t84;
t690 = -t66 + t21;
t495 = t564 * pkin(2);
t384 = -t495 - pkin(3);
t602 = t402 * pkin(4);
t356 = t384 - t602;
t266 = pkin(5) * t347 - pkin(11) * t348 + t356;
t164 = t266 * t400 - t275 * t395;
t689 = qJD(6) * t164 + t395 * t703 + t400 * t704;
t165 = t266 * t395 + t275 * t400;
t688 = -qJD(6) * t165 - t395 * t704 + t400 * t703;
t685 = -pkin(5) * t315 - t686;
t297 = pkin(2) * t394 + t381 - t467;
t309 = qJ(3) * t548 + t530;
t219 = t297 * t563 + t309 * t564;
t197 = pkin(9) * t394 + t219;
t377 = pkin(2) * t548;
t484 = pkin(9) * t323 + t377;
t607 = pkin(1) * t393;
t236 = pkin(3) * t322 - t484 - t607;
t136 = -t197 * t397 + t236 * t402;
t295 = t323 * t402 + t394 * t397;
t114 = pkin(4) * t322 - pkin(10) * t295 + t136;
t137 = t197 * t402 + t236 * t397;
t123 = pkin(10) * t294 + t137;
t682 = t114 * t396 + t123 * t401;
t161 = t225 * t395 - t315 * t400;
t520 = qJD(6) * t400;
t429 = -t292 * t395 + t348 * t520;
t681 = t161 + t429;
t162 = -t225 * t400 - t315 * t395;
t521 = qJD(6) * t395;
t428 = t292 * t400 + t348 * t521;
t680 = t162 + t428;
t590 = mrSges(4,3) * t315;
t679 = mrSges(4,1) * t376 + mrSges(5,1) * t272 + mrSges(5,2) * t437 + t590;
t190 = t294 * t396 + t295 * t401;
t312 = t322 * t400;
t158 = -t190 * t395 + t312;
t285 = -t323 * t390 + t391 * t394;
t286 = t323 * t391 + t390 * t394;
t676 = -t285 * t705 + t286 * t706;
t230 = t262 * t390 - t391 * t549;
t231 = t262 * t391 + t390 * t549;
t675 = t230 * t705 + t231 * t706;
t257 = t345 * t399 - t404 * t531;
t547 = t393 * t404;
t226 = t257 * t390 - t391 * t547;
t227 = -t257 * t391 - t390 * t547;
t674 = -t226 * t705 + t227 * t706;
t454 = mrSges(7,1) * t395 + mrSges(7,2) * t400;
t667 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t459 = m(5) * pkin(9) - t667;
t673 = -t454 - t459;
t540 = t399 * t403;
t542 = t398 * t404;
t338 = -t394 * t540 - t542;
t671 = -t397 * t51 + t402 * t50;
t23 = mrSges(7,1) * t83 - mrSges(7,3) * t47;
t24 = -mrSges(7,2) * t83 + mrSges(7,3) * t48;
t670 = -t23 * t395 + t24 * t400;
t587 = Ifges(3,4) * t398;
t668 = pkin(1) * (mrSges(3,1) * t398 + mrSges(3,2) * t403) - t398 * (Ifges(3,1) * t403 - t587) / 0.2e1;
t100 = -mrSges(7,2) * t166 + mrSges(7,3) * t143;
t101 = mrSges(7,1) * t166 - mrSges(7,3) * t144;
t145 = -mrSges(6,2) * t306 + mrSges(6,3) * t476;
t666 = t100 * t400 - t101 * t395 + t145;
t560 = t110 * t396;
t54 = t401 * t97 - t560;
t665 = mrSges(6,2) * t141 - mrSges(6,3) * t54;
t11 = -qJD(5) * t55 + t34 * t401 - t396 * t41;
t317 = qJD(2) * t322;
t217 = qJD(4) * t294 - t317 * t402;
t316 = qJD(2) * t323;
t374 = t403 * t515;
t281 = t374 + (-qJD(2) * t487 + t526) * t393;
t488 = t597 * t393;
t282 = -t490 + (-t403 * t488 - t379) * qJD(2);
t179 = t281 * t564 + t282 * t563;
t491 = t393 * t527;
t233 = pkin(2) * t491 + pkin(3) * t316 + pkin(9) * t317;
t88 = -qJD(4) * t137 - t179 * t397 + t233 * t402;
t63 = pkin(4) * t316 - pkin(10) * t217 + t88;
t216 = -qJD(4) * t295 + t317 * t397;
t87 = t179 * t402 - t197 * t525 + t233 * t397 + t236 * t524;
t71 = pkin(10) * t216 + t87;
t20 = -qJD(5) * t682 - t396 * t71 + t401 * t63;
t664 = -mrSges(4,1) - t697;
t471 = pkin(8) * t491;
t268 = -qJD(1) * t471 + t496;
t269 = -pkin(8) * t331 - t398 * t468 + t371;
t663 = t269 * mrSges(3,1) + t125 * mrSges(4,1) - t268 * mrSges(3,2) - t126 * mrSges(4,2) + Ifges(3,5) * t331 + Ifges(4,5) * t252 + Ifges(3,6) * t330;
t117 = pkin(5) * t439 - pkin(11) * t476;
t52 = -pkin(5) * t306 - t54;
t430 = t52 * t454;
t142 = Ifges(7,4) * t143;
t78 = Ifges(7,1) * t144 + Ifges(7,5) * t166 + t142;
t566 = t400 * t78;
t608 = t395 / 0.2e1;
t163 = Ifges(6,4) * t476;
t573 = t306 * Ifges(6,5);
t105 = Ifges(6,1) * t439 + t163 + t573;
t637 = -t105 / 0.2e1;
t582 = Ifges(7,4) * t144;
t77 = Ifges(7,2) * t143 + Ifges(7,6) * t166 + t582;
t659 = -t430 - t566 / 0.2e1 + t77 * t608 + t637 - t665;
t15 = Ifges(7,4) * t47 + Ifges(7,2) * t48 + Ifges(7,6) * t83;
t653 = t15 / 0.2e1;
t16 = Ifges(7,1) * t47 + Ifges(7,4) * t48 + Ifges(7,5) * t83;
t652 = t16 / 0.2e1;
t641 = Ifges(5,4) * t629 + Ifges(5,2) * t628 + Ifges(5,6) * t620;
t640 = Ifges(5,1) * t629 + Ifges(5,4) * t628 + Ifges(5,5) * t620;
t636 = t105 / 0.2e1;
t574 = t437 * Ifges(5,4);
t150 = t272 * Ifges(5,2) + t308 * Ifges(5,6) - t574;
t631 = t150 / 0.2e1;
t267 = Ifges(5,4) * t272;
t151 = -Ifges(5,1) * t437 + t308 * Ifges(5,5) + t267;
t630 = t151 / 0.2e1;
t619 = -t272 / 0.2e1;
t618 = t437 / 0.2e1;
t617 = -t437 / 0.2e1;
t614 = -t308 / 0.2e1;
t612 = -t315 / 0.2e1;
t605 = pkin(4) * t437;
t604 = pkin(4) * t396;
t603 = pkin(4) * t401;
t595 = mrSges(4,1) * t322;
t593 = mrSges(5,2) * t402;
t591 = mrSges(4,3) * t314;
t589 = mrSges(7,3) * t395;
t588 = mrSges(7,3) * t400;
t586 = Ifges(4,4) * t315;
t585 = Ifges(5,4) * t397;
t584 = Ifges(5,4) * t402;
t581 = Ifges(7,4) * t395;
t580 = Ifges(7,4) * t400;
t579 = t127 * mrSges(5,3);
t578 = t128 * mrSges(5,3);
t571 = t376 * Ifges(3,5);
t570 = t376 * Ifges(3,6);
t554 = t322 * t395;
t552 = t348 * t395;
t551 = t348 * t400;
t516 = Ifges(6,5) * t84 + Ifges(6,6) * t85 + Ifges(6,3) * t242;
t510 = -m(3) * pkin(1) - mrSges(2,1);
t506 = t397 * t549;
t368 = t397 * t547;
t500 = t566 / 0.2e1;
t499 = Ifges(5,5) * t154 + Ifges(5,6) * t155 + Ifges(5,3) * t250;
t492 = t403 * t529;
t485 = -t521 / 0.2e1;
t483 = -pkin(5) * t230 + pkin(11) * t231;
t333 = pkin(2) * t546 - t488;
t475 = -t333 * t399 + t389 * t404;
t474 = -m(3) * pkin(8) - mrSges(3,3) - mrSges(4,3);
t470 = mrSges(3,3) * t493;
t469 = mrSges(3,3) * t492;
t465 = t696 * pkin(2);
t464 = t237 * pkin(4);
t458 = mrSges(5,1) * t294 - mrSges(5,2) * t295;
t453 = Ifges(5,1) * t402 - t585;
t452 = Ifges(7,1) * t400 - t581;
t451 = -Ifges(5,2) * t397 + t584;
t450 = -Ifges(7,2) * t395 + t580;
t449 = Ifges(5,5) * t402 - Ifges(5,6) * t397;
t448 = Ifges(7,5) * t400 - Ifges(7,6) * t395;
t447 = t25 * t400 + t26 * t395;
t446 = -t25 * t395 + t26 * t400;
t59 = pkin(11) * t322 + t682;
t218 = t297 * t564 - t309 * t563;
t196 = -pkin(3) * t394 - t218;
t157 = t196 - t463;
t438 = t294 * t401 - t295 * t396;
t94 = -pkin(5) * t438 - pkin(11) * t190 + t157;
t30 = t395 * t94 + t400 * t59;
t29 = -t395 * t59 + t400 * t94;
t178 = t281 * t563 - t282 * t564;
t64 = t114 * t401 - t123 * t396;
t159 = t190 * t400 + t554;
t193 = -mrSges(5,2) * t308 + mrSges(5,3) * t272;
t194 = mrSges(5,1) * t308 + mrSges(5,3) * t437;
t440 = t193 * t402 - t194 * t397;
t431 = t257 * t397 - t402 * t547;
t426 = t143 * t450;
t425 = t144 * t452;
t424 = t166 * t448;
t19 = t114 * t522 - t123 * t523 + t396 * t63 + t401 * t71;
t412 = t394 * t345;
t261 = t399 * t412 - t404 * t414;
t388 = pkin(3) + t602;
t405 = -pkin(10) - pkin(9);
t420 = pkin(4) * t506 + t261 * t405 + t262 * t388 + t475;
t131 = -pkin(4) * t216 + t178;
t416 = t431 * pkin(4);
t413 = -qJD(6) * t447 - t3 * t395;
t409 = t2 * t400 + t413;
t408 = (-t100 * t395 - t101 * t400) * qJD(6) + t670;
t8 = -pkin(5) * t242 - t11;
t406 = t11 * mrSges(6,1) - t10 * mrSges(6,2) + t16 * t608 + t2 * t588 + t400 * t653 + t8 * t455 + (Ifges(7,1) * t395 + t580) * t649 + (Ifges(7,2) * t400 + t581) * t648 + t516 + t77 * t485 + (Ifges(7,5) * t395 + Ifges(7,6) * t400) * t644 + (t430 + t500) * qJD(6) + (t426 + t425 + t424) * qJD(6) / 0.2e1;
t387 = -pkin(5) - t603;
t372 = Ifges(3,4) * t492;
t366 = Ifges(3,3) * t375;
t365 = Ifges(4,3) * t375;
t349 = -t377 - t607;
t341 = -pkin(8) * t550 + t381;
t340 = (-mrSges(3,1) * t403 + mrSges(3,2) * t398) * t393;
t339 = -t394 * t543 + t538;
t337 = -t394 * t542 - t540;
t328 = t374 - t471;
t327 = t530 * qJD(1);
t326 = -pkin(8) * t493 + t373;
t325 = -mrSges(3,2) * t376 + t469;
t324 = mrSges(3,1) * t376 - t470;
t307 = Ifges(4,4) * t314;
t299 = Ifges(3,1) * t493 + t372 + t571;
t298 = t570 + (Ifges(3,2) * t403 + t587) * t529;
t283 = -mrSges(4,2) * t376 + t591;
t280 = t285 * pkin(5);
t258 = -t399 * t414 - t404 * t412;
t248 = Ifges(4,6) * t251;
t247 = t252 * mrSges(4,2);
t239 = -mrSges(4,1) * t314 - mrSges(4,2) * t315;
t238 = t262 * t402 + t506;
t235 = mrSges(4,1) * t375 - mrSges(4,3) * t252;
t234 = -mrSges(4,2) * t375 + mrSges(4,3) * t251;
t223 = -t315 * Ifges(4,1) + t376 * Ifges(4,5) + t307;
t222 = Ifges(4,2) * t314 + Ifges(4,6) * t376 - t586;
t214 = t226 * pkin(5);
t148 = t231 * t400 - t261 * t395;
t147 = -t231 * t395 - t261 * t400;
t130 = -mrSges(5,2) * t250 + mrSges(5,3) * t155;
t129 = mrSges(5,1) * t250 - mrSges(5,3) * t154;
t116 = -mrSges(6,1) * t476 + mrSges(6,2) * t439;
t108 = qJD(5) * t190 - t216 * t401 + t217 * t396;
t107 = qJD(5) * t438 + t216 * t396 + t217 * t401;
t102 = -mrSges(5,1) * t155 + mrSges(5,2) * t154;
t99 = t117 - t605;
t75 = -qJD(6) * t159 - t107 * t395 + t316 * t400;
t74 = qJD(6) * t158 + t107 * t400 + t316 * t395;
t67 = -mrSges(6,2) * t242 + mrSges(6,3) * t85;
t58 = -pkin(5) * t322 - t64;
t57 = t109 * t401 - t560;
t56 = t109 * t396 + t539;
t42 = pkin(5) * t108 - pkin(11) * t107 + t131;
t37 = -mrSges(6,1) * t85 + mrSges(6,2) * t84;
t36 = t117 * t395 + t400 * t54;
t35 = t117 * t400 - t395 * t54;
t32 = t395 * t99 + t400 * t57;
t31 = -t395 * t57 + t400 * t99;
t18 = -pkin(5) * t316 - t20;
t17 = pkin(11) * t316 + t19;
t5 = -qJD(6) * t30 - t17 * t395 + t400 * t42;
t4 = qJD(6) * t29 + t17 * t400 + t395 * t42;
t1 = [t190 * t650 + t159 * t652 + t158 * t653 + t216 * t631 + t107 * t636 + t295 * t640 + t294 * t641 + (Ifges(5,1) * t217 + Ifges(5,4) * t216 + Ifges(5,5) * t316) * t617 + (Ifges(5,5) * t295 + Ifges(5,6) * t294 + Ifges(5,3) * t322) * t620 + (Ifges(5,4) * t295 + Ifges(5,2) * t294 + Ifges(5,6) * t322) * t628 + (Ifges(5,1) * t295 + Ifges(5,4) * t294 + Ifges(5,5) * t322) * t629 + t217 * t630 + t296 * t595 + (t516 + t499) * t322 / 0.2e1 + m(6) * (t10 * t682 + t11 * t64 + t131 * t141 + t157 * t89 + t19 * t55 + t20 * t54) + t682 * t67 + (-t10 * t322 + t107 * t141 + t190 * t89 - t316 * t55) * mrSges(6,2) + (t474 * g(2) * t399 + (mrSges(3,1) * t330 - mrSges(3,2) * t331 + (m(3) * t607 - t340) * qJDD(1)) * pkin(1) + (t474 - t593) * g(1) * t404 + (mrSges(3,3) * t268 + Ifges(3,4) * t331 + Ifges(3,2) * t330 + Ifges(3,6) * t375) * t403 + (-mrSges(3,3) * t269 + Ifges(3,1) * t331 + Ifges(3,4) * t330 + Ifges(3,5) * t375) * t398 + ((t571 / 0.2e1 - t326 * mrSges(3,3) + t299 / 0.2e1) * t403 + (-t570 / 0.2e1 - t327 * mrSges(3,3) - t298 / 0.2e1 + (m(4) * t335 + t239) * pkin(2)) * t398 + (t403 * (Ifges(3,4) * t403 - Ifges(3,2) * t398) / 0.2e1 - t668) * t529) * qJD(2)) * t393 - t679 * t178 - t701 * t438 + (-Ifges(4,1) * t317 - Ifges(4,4) * t316) * t612 + (-t126 * t322 + t183 * t317 - t184 * t316) * mrSges(4,3) + t376 * (-Ifges(4,5) * t317 - Ifges(4,6) * t316) / 0.2e1 + t335 * (mrSges(4,1) * t316 - mrSges(4,2) * t317) + t314 * (-Ifges(4,4) * t317 - Ifges(4,2) * t316) / 0.2e1 + m(3) * (t268 * t530 + t269 * t341 - t326 * t329 + t327 * t328) + (Ifges(6,4) * t190 + Ifges(6,6) * t322) * t642 + (Ifges(6,4) * t107 + Ifges(6,6) * t316) * t624 + (t341 * mrSges(3,1) - t530 * mrSges(3,2) - Ifges(4,6) * t322 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t394) * t375 + (t330 * t530 - t331 * t341) * mrSges(3,3) + (Ifges(7,5) * t159 + Ifges(7,6) * t158) * t644 + (Ifges(7,5) * t74 + Ifges(7,6) * t75) * t626 + (mrSges(4,2) * t296 - mrSges(4,3) * t125 + Ifges(4,1) * t252 + Ifges(4,4) * t251 + Ifges(4,5) * t375) * t323 + (Ifges(7,4) * t159 + Ifges(7,2) * t158) * t648 + (Ifges(7,4) * t74 + Ifges(7,2) * t75) * t634 + t699 * t108 + m(5) * (t119 * t196 + t127 * t88 + t128 * t87 + t136 * t51 + t137 * t50 + t175 * t178) + m(7) * (t18 * t52 + t2 * t30 + t25 * t5 + t26 * t4 + t29 * t3 + t58 * t8) + m(4) * (t125 * t218 + t126 * t219 - t178 * t183 + t179 * t184 + t296 * t349) + (-Ifges(4,2) * t322 + Ifges(4,6) * t394 / 0.2e1 - t349 * mrSges(4,1)) * t251 + (t158 * t2 - t159 * t3 - t25 * t74 + t26 * t75) * mrSges(7,3) - t252 * Ifges(4,4) * t322 + (Ifges(6,5) * t107 + Ifges(6,3) * t316) * t615 + (Ifges(6,5) * t190 + Ifges(6,3) * t322) * t621 + (t248 / 0.2e1 + t365 / 0.2e1 + t366 / 0.2e1 + t663) * t394 + (Ifges(7,1) * t159 + Ifges(7,4) * t158) * t649 + (Ifges(7,1) * t74 + Ifges(7,4) * t75) * t632 + (-m(4) * t475 - t262 * mrSges(4,1) + mrSges(2,2) * t399 - m(7) * (pkin(5) * t231 + t420) - t148 * mrSges(7,1) - t147 * mrSges(7,2) - m(6) * t420 - t231 * mrSges(6,1) - m(5) * (pkin(3) * t262 + t475) - t238 * mrSges(5,1) - t237 * mrSges(5,2) - t339 * mrSges(3,1) - t338 * mrSges(3,2) + t510 * t404 + t473 * t230 + t459 * t261) * g(2) + (Ifges(6,1) * t190 + Ifges(6,5) * t322) * t643 + (Ifges(6,1) * t107 + Ifges(6,5) * t316) * t622 + Ifges(2,3) * qJDD(1) + t349 * t247 + t328 * t325 - t329 * t324 + t50 * (-mrSges(5,2) * t322 + mrSges(5,3) * t294) + t51 * (mrSges(5,1) * t322 - mrSges(5,3) * t295) + t11 * (mrSges(6,1) * t322 - mrSges(6,3) * t190) - t316 * t222 / 0.2e1 - t317 * t223 / 0.2e1 + t308 * (Ifges(5,5) * t217 + Ifges(5,6) * t216 + Ifges(5,3) * t316) / 0.2e1 + t272 * (Ifges(5,4) * t217 + Ifges(5,2) * t216 + Ifges(5,6) * t316) / 0.2e1 + t127 * (mrSges(5,1) * t316 - mrSges(5,3) * t217) + t128 * (-mrSges(5,2) * t316 + mrSges(5,3) * t216) + t54 * (mrSges(6,1) * t316 - mrSges(6,3) * t107) + t179 * t283 + t219 * t234 + t218 * t235 + t175 * (-mrSges(5,1) * t216 + mrSges(5,2) * t217) + t87 * t193 + t88 * t194 + t196 * t102 + t8 * (-mrSges(7,1) * t158 + mrSges(7,2) * t159) + t157 * t37 + t19 * t145 + t20 * t146 + t136 * t129 + t137 * t130 + t131 * t116 + t4 * t100 - t119 * t458 + t5 * t101 + t18 * t95 + (-t368 * mrSges(5,1) - t337 * mrSges(3,1) + t696 * mrSges(3,2) + (t333 * t511 + mrSges(2,2)) * t404 + (t389 * t511 - t510) * t399 - (mrSges(4,1) + t418) * t257 + t473 * t226 + t708 * t227 + t673 * t258 + t691 * (-pkin(4) * t368 - t257 * t388 + t258 * t405)) * g(1) + t29 * t23 + t30 * t24 + t698 * t316 / 0.2e1 + t58 * t21 + t64 * t66 + t52 * (-mrSges(7,1) * t75 + mrSges(7,2) * t74) + t75 * t77 / 0.2e1 + t74 * t78 / 0.2e1; (t630 - t579) * t524 + (-t150 / 0.2e1 - t578) * t525 + t551 * t652 + t556 * t631 - t225 * t637 + t397 * t640 + t402 * t641 + t222 * t612 + (-Ifges(5,3) * t315 + t314 * t449) * t614 + (-Ifges(5,5) * t315 + t314 * t453) * t618 + (-Ifges(5,6) * t315 + t314 * t451) * t619 + (Ifges(5,5) * t397 + Ifges(5,6) * t402) * t620 + (Ifges(5,2) * t402 + t585) * t628 + (Ifges(5,1) * t397 + t584) * t629 + t183 * t591 + (-Ifges(7,5) * t428 - Ifges(7,6) * t429) * t626 + (Ifges(7,5) * t162 + Ifges(7,6) * t161) * t627 + t234 * t494 + t235 * t495 + (t272 * t451 + t308 * t449 - t437 * t453) * qJD(4) / 0.2e1 + (t10 * t275 + t11 * t436 + t141 * t672 + t356 * t89 + t54 * t686 + t55 * t687) * m(6) + (t164 * t3 + t165 * t2 + t25 * t688 + t26 * t689 - t436 * t8 + t52 * t685) * m(7) - t690 * t436 + (t89 * mrSges(6,2) - t11 * mrSges(6,3) + t448 * t644 + t450 * t648 + t452 * t649 + t454 * t8 + t485 * t78 + 0.2e1 * t650) * t348 + (-Ifges(7,1) * t428 - Ifges(7,4) * t429) * t632 + (Ifges(7,1) * t162 + Ifges(7,4) * t161) * t633 + (mrSges(3,2) * t339 - t691 * (t261 * t388 - t262 * t405) + t664 * t261 + t673 * t262 + (-pkin(2) * t511 - mrSges(3,1)) * t338) * g(1) + t701 * t347 + (t141 * t225 - t315 * t55) * mrSges(6,2) + (-Ifges(6,1) * t225 - Ifges(6,5) * t315) * t623 + (-Ifges(6,5) * t225 - Ifges(6,3) * t315) * t616 - t54 * (-mrSges(6,1) * t315 + mrSges(6,3) * t225) + (-Ifges(6,4) * t225 - Ifges(6,6) * t315) * t625 - (Ifges(4,2) * t315 + t223 + t307) * t314 / 0.2e1 + (t470 + t324) * t327 + (-Ifges(7,4) * t428 - Ifges(7,2) * t429) * t634 + (Ifges(7,4) * t162 + Ifges(7,2) * t161) * t635 + (t469 - t325) * t326 - ((-Ifges(3,2) * t493 + t299 + t372) * t403 + t376 * (Ifges(3,5) * t403 - Ifges(3,6) * t398)) * t529 / 0.2e1 + (-mrSges(3,1) * t696 - mrSges(3,2) * t337 + t692 * t465 - t691 * (t257 * t405 + t258 * t388 + t465) + t664 * t258 - t673 * t257) * g(2) + t699 * t293 + t127 * mrSges(5,1) * t315 - t128 * mrSges(5,2) * t315 + ((t125 * t564 + t126 * t563) * pkin(2) + t183 * t198 - t184 * t199 - t335 * t472) * m(4) + t365 + t366 - t184 * t590 + (t119 * t384 - t127 * t134 - t128 * t135 - t175 * t198) * m(5) + t298 * t493 / 0.2e1 - (Ifges(6,1) * t622 + Ifges(6,4) * t624 + Ifges(6,5) * t615 + t500 + t636 + t665) * t292 + t248 - t151 * t555 / 0.2e1 - t15 * t552 / 0.2e1 - t239 * t472 - t376 * (Ifges(4,5) * t314 + Ifges(4,6) * t315) / 0.2e1 + t384 * t102 + t356 * t37 - t335 * (-mrSges(4,1) * t315 + mrSges(4,2) * t314) + t275 * t67 - t199 * t283 - t135 * t193 - t134 * t194 - t162 * t78 / 0.2e1 + t164 * t23 + t165 * t24 + t119 * t457 + t668 * qJD(1) ^ 2 * t393 ^ 2 + t308 * t175 * (mrSges(5,1) * t397 + t593) + (-t194 * t524 - t193 * t525 - t397 * t129 + m(5) * ((-t127 * t402 - t128 * t397) * qJD(4) + t671) + t402 * t130) * t383 + (t127 * t555 + t128 * t556 + t671) * mrSges(5,3) + t672 * t116 + t679 * t198 + t663 - t681 * t77 / 0.2e1 + (mrSges(7,1) * t681 - mrSges(7,2) * t680) * t52 + (-t2 * t552 + t25 * t680 - t26 * t681 - t3 * t551) * mrSges(7,3) + t695 * t224 + (-m(4) * t377 - m(5) * t484 + t340 + t595 - t691 * (-t322 * t388 - t323 * t405 + t377) + t697 * t322 + (-t454 + t667) * t323) * g(3) + (Ifges(4,1) * t314 + t586 + t698) * t315 / 0.2e1 + t685 * t95 + t686 * t146 + t687 * t145 + t688 * t101 + t689 * t100; (t67 + t408) * t348 + t690 * t347 - t666 * t292 + (t116 - t679) * t315 + t247 + t402 * t129 + t397 * t130 - t251 * mrSges(4,1) + t225 * t145 - t161 * t101 - t162 * t100 + (-t283 - t440) * t314 + t440 * qJD(4) + t683 * t678 + (-t161 * t25 - t162 * t26 - t292 * t446 + t347 * t8 + t348 * t409 + t52 * t678) * m(7) + (t10 * t348 - t11 * t347 + t141 * t315 - t678 * t54 + t55 * t709) * m(6) + (t175 * t315 + t397 * t50 + t402 * t51 + t308 * (-t127 * t397 + t128 * t402)) * m(5) + (-t183 * t315 - t184 * t314 + t296) * m(4) + (-t394 * g(3) + (-g(1) * t399 + g(2) * t404) * t393) * t511; t150 * t617 + t66 * t603 + t67 * t604 + t272 * t579 - t437 * t578 + (Ifges(5,5) * t272 + Ifges(5,6) * t437) * t614 + (Ifges(5,2) * t437 + t151 + t267) * t619 - t175 * (-mrSges(5,1) * t437 + mrSges(5,2) * t272) + (Ifges(6,1) * t623 + Ifges(6,4) * t625 + Ifges(6,5) * t616 + t25 * t588 + t26 * t589 + t448 * t627 + t450 * t635 + t452 * t633 + t659) * t476 + (-t25 * t31 - t26 * t32 - t52 * t56 + t387 * t8 + (t396 * t52 + t401 * t446) * qJD(5) * pkin(4)) * m(7) + (-m(6) * t416 - t431 * mrSges(5,1) - (t257 * t402 + t368) * mrSges(5,2) - m(7) * (pkin(11) * t227 + t214 + t416) + t674) * g(2) + (Ifges(5,1) * t272 + t574) * t618 + ((t10 * t396 + t11 * t401 + (-t396 * t54 + t401 * t55) * qJD(5)) * pkin(4) + t141 * t605 + t54 * t56 - t55 * t57) * m(6) + (-t25 * t520 - t26 * t521) * mrSges(7,3) - t3 * t589 + t406 + t116 * t605 + t499 + t387 * t21 - t127 * t193 + t128 * t194 - t57 * t145 - t32 * t100 - t31 * t101 + t666 * pkin(4) * t522 + (m(7) * t409 - t100 * t521 - t101 * t520 + t670) * (pkin(11) + t604) + (-m(7) * (t464 + t483) - m(6) * t464 - mrSges(5,1) * t237 + mrSges(5,2) * t238 + t675) * g(1) + (-m(6) * t463 - m(7) * (pkin(11) * t286 + t280 + t463) - t458 + t676) * g(3) + t695 * t439 - t50 * mrSges(5,2) + t51 * mrSges(5,1) + t683 * (pkin(4) * t523 - t56); (-m(7) * t214 + t674) * g(2) + (-m(7) * t483 + t675) * g(1) + (-m(7) * t280 + t676) * g(3) + (-t424 / 0.2e1 - t425 / 0.2e1 - t426 / 0.2e1 - t573 / 0.2e1 - t163 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t439 + t447 * mrSges(7,3) + t659) * t476 + t413 * mrSges(7,3) + t409 * t655 + ((-g(2) * t227 - g(3) * t286) * m(7) + t408) * pkin(11) - t683 * t55 - t8 * t656 + (t702 + t707) * t439 + t406 - pkin(5) * t21 - m(7) * (t25 * t35 + t26 * t36 + t52 * t55) - t54 * t145 - t36 * t100 - t35 * t101; -t52 * (mrSges(7,1) * t144 + mrSges(7,2) * t143) + (Ifges(7,1) * t143 - t582) * t633 + t77 * t632 + (Ifges(7,5) * t143 - Ifges(7,6) * t144) * t627 - t25 * t100 + t26 * t101 - g(1) * (mrSges(7,1) * t147 - mrSges(7,2) * t148) - g(2) * ((-t227 * t395 - t258 * t400) * mrSges(7,1) + (-t227 * t400 + t258 * t395) * mrSges(7,2)) - g(3) * ((-t286 * t395 + t312) * mrSges(7,1) + (-t286 * t400 - t554) * mrSges(7,2)) + (t143 * t25 + t144 * t26) * mrSges(7,3) + t14 + (-Ifges(7,2) * t144 + t142 + t78) * t635 + t662;];
tau  = t1;
