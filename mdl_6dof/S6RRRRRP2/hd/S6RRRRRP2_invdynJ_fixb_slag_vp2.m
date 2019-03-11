% Calculate vector of inverse dynamics joint torques for
% S6RRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:00:42
% EndTime: 2019-03-10 01:01:26
% DurationCPUTime: 27.34s
% Computational Cost: add. (22911->808), mult. (52502->1007), div. (0->0), fcn. (38837->14), ass. (0->392)
t362 = sin(qJ(5));
t367 = cos(qJ(5));
t364 = sin(qJ(3));
t365 = sin(qJ(2));
t369 = cos(qJ(3));
t370 = cos(qJ(2));
t293 = -t364 * t365 + t369 * t370;
t275 = t293 * qJD(1);
t294 = t364 * t370 + t365 * t369;
t276 = t294 * qJD(1);
t363 = sin(qJ(4));
t368 = cos(qJ(4));
t401 = t275 * t363 + t368 * t276;
t359 = qJD(2) + qJD(3);
t437 = qJD(4) + t359;
t204 = t362 * t401 - t367 * t437;
t434 = t368 * t275 - t276 * t363;
t223 = qJD(5) - t434;
t205 = t362 * t437 + t367 * t401;
t529 = Ifges(6,4) * t205;
t112 = -Ifges(6,2) * t204 + Ifges(6,6) * t223 + t529;
t676 = -t112 / 0.2e1;
t657 = Ifges(6,1) + Ifges(7,1);
t629 = Ifges(6,5) + Ifges(7,4);
t357 = t370 * pkin(2);
t346 = t357 + pkin(1);
t323 = t346 * qJD(1);
t248 = -pkin(3) * t275 - t323;
t624 = t437 * Ifges(5,5);
t675 = t248 * mrSges(5,2) + t624 / 0.2e1 + t362 * t676;
t471 = qJD(1) * qJD(2);
t303 = qJDD(1) * t370 - t365 * t471;
t304 = qJDD(1) * t365 + t370 * t471;
t387 = t293 * qJD(3);
t208 = qJD(1) * t387 + t303 * t364 + t304 * t369;
t388 = t294 * qJD(3);
t209 = -qJD(1) * t388 + t303 * t369 - t304 * t364;
t118 = qJD(4) * t434 + t208 * t368 + t209 * t363;
t358 = qJDD(2) + qJDD(3);
t353 = qJDD(4) + t358;
t474 = qJD(5) * t204;
t71 = t367 * t118 + t362 * t353 - t474;
t585 = t71 / 0.2e1;
t72 = qJD(5) * t205 + t362 * t118 - t367 * t353;
t583 = t72 / 0.2e1;
t119 = -qJD(4) * t401 - t208 * t363 + t209 * t368;
t117 = qJDD(5) - t119;
t582 = t117 / 0.2e1;
t674 = -mrSges(6,3) - mrSges(7,2);
t673 = -Ifges(6,4) + Ifges(7,5);
t628 = Ifges(7,2) + Ifges(6,3);
t656 = Ifges(6,6) - Ifges(7,6);
t199 = Ifges(6,4) * t204;
t526 = Ifges(7,5) * t204;
t651 = t205 * t657 + t223 * t629 - t199 + t526;
t472 = qJD(5) * t367;
t618 = t434 * t367;
t672 = t472 - t618;
t473 = qJD(5) * t362;
t619 = t434 * t362;
t671 = -t473 + t619;
t573 = -t434 / 0.2e1;
t576 = -t223 / 0.2e1;
t578 = -t205 / 0.2e1;
t579 = t204 / 0.2e1;
t580 = -t204 / 0.2e1;
t372 = -pkin(8) - pkin(7);
t325 = t372 * t370;
t300 = qJD(1) * t325;
t277 = t364 * t300;
t324 = t372 * t365;
t299 = qJD(1) * t324;
t284 = qJD(2) * pkin(2) + t299;
t233 = t369 * t284 + t277;
t269 = t276 * pkin(9);
t202 = t233 - t269;
t189 = pkin(3) * t359 + t202;
t280 = t369 * t300;
t234 = t284 * t364 - t280;
t552 = pkin(9) * t275;
t203 = t234 + t552;
t192 = t368 * t203;
t129 = t363 * t189 + t192;
t124 = pkin(10) * t437 + t129;
t136 = -pkin(4) * t434 - pkin(10) * t401 + t248;
t53 = -t124 * t362 + t136 * t367;
t42 = -pkin(5) * t223 + qJD(6) - t53;
t54 = t124 * t367 + t136 * t362;
t43 = qJ(6) * t223 + t54;
t623 = t437 * Ifges(5,6);
t587 = -t248 * mrSges(5,1) - t53 * mrSges(6,1) + t42 * mrSges(7,1) + t54 * mrSges(6,2) - t43 * mrSges(7,3) + t623 / 0.2e1;
t670 = -Ifges(5,2) * t573 + Ifges(6,6) * t579 + Ifges(7,6) * t580 + t576 * t628 + t578 * t629 + t587;
t524 = Ifges(7,5) * t367;
t411 = Ifges(7,3) * t362 + t524;
t527 = Ifges(6,4) * t367;
t415 = -Ifges(6,2) * t362 + t527;
t198 = Ifges(7,5) * t205;
t109 = Ifges(7,6) * t223 + Ifges(7,3) * t204 + t198;
t493 = t362 * t109;
t566 = -t367 / 0.2e1;
t525 = Ifges(7,5) * t362;
t528 = Ifges(6,4) * t362;
t606 = t367 * t657 + t525 - t528;
t607 = -t362 * t656 + t367 * t629;
t633 = -t401 / 0.2e1;
t191 = t363 * t203;
t128 = t368 * t189 - t191;
t123 = -pkin(4) * t437 - t128;
t419 = t362 * mrSges(7,1) - t367 * mrSges(7,3);
t421 = mrSges(6,1) * t362 + mrSges(6,2) * t367;
t55 = t204 * pkin(5) - t205 * qJ(6) + t123;
t644 = t123 * t421 + t55 * t419;
t669 = Ifges(5,1) * t633 + t411 * t580 + t415 * t579 - t493 / 0.2e1 - t644 + t651 * t566 + t578 * t606 + t576 * t607 - t675;
t668 = t629 * t582 + t673 * t583 + t585 * t657;
t632 = mrSges(6,1) + mrSges(7,1);
t631 = mrSges(6,2) - mrSges(7,3);
t539 = mrSges(7,2) * t204;
t150 = mrSges(7,3) * t223 - t539;
t535 = mrSges(6,3) * t204;
t151 = -mrSges(6,2) * t223 - t535;
t667 = -t151 - t150;
t534 = mrSges(6,3) * t205;
t152 = mrSges(6,1) * t223 - t534;
t538 = mrSges(7,2) * t205;
t153 = -mrSges(7,1) * t223 + t538;
t666 = t153 - t152;
t410 = pkin(5) * t367 + qJ(6) * t362;
t306 = -pkin(4) - t410;
t361 = qJ(2) + qJ(3);
t356 = qJ(4) + t361;
t341 = sin(t356);
t420 = -t367 * mrSges(7,1) - t362 * mrSges(7,3);
t541 = mrSges(6,1) * t367;
t665 = (-m(7) * t306 - t420 + t541) * t341;
t520 = t401 * mrSges(5,3);
t615 = mrSges(5,1) * t437 - mrSges(6,1) * t204 - mrSges(6,2) * t205 - t520;
t591 = -m(6) * t123 + t615;
t664 = -m(5) * t128 - t591;
t530 = Ifges(5,4) * t401;
t626 = Ifges(5,2) * t434;
t162 = t530 + t623 + t626;
t220 = Ifges(5,4) * t434;
t163 = Ifges(5,1) * t401 + t220 + t624;
t19 = Ifges(7,5) * t71 + Ifges(7,6) * t117 + Ifges(7,3) * t72;
t20 = Ifges(6,4) * t71 - Ifges(6,2) * t72 + Ifges(6,6) * t117;
t290 = t304 * pkin(7);
t245 = qJDD(2) * pkin(2) - pkin(8) * t304 - t290;
t289 = t303 * pkin(7);
t246 = pkin(8) * t303 + t289;
t477 = qJD(3) * t369;
t478 = qJD(3) * t364;
t147 = t364 * t245 + t369 * t246 + t284 * t477 + t300 * t478;
t106 = pkin(9) * t209 + t147;
t475 = qJD(4) * t368;
t476 = qJD(4) * t363;
t148 = -qJD(3) * t234 + t369 * t245 - t246 * t364;
t96 = pkin(3) * t358 - pkin(9) * t208 + t148;
t29 = -t363 * t106 - t189 * t476 - t203 * t475 + t368 * t96;
t26 = -pkin(4) * t353 - t29;
t28 = t368 * t106 + t189 * t475 - t203 * t476 + t363 * t96;
t439 = t472 / 0.2e1;
t440 = -t473 / 0.2e1;
t449 = t493 / 0.2e1;
t571 = t401 / 0.2e1;
t584 = -t72 / 0.2e1;
t617 = -t204 * t656 + t205 * t629 + t223 * t628;
t25 = pkin(10) * t353 + t28;
t518 = qJDD(1) * pkin(1);
t266 = -pkin(2) * t303 - t518;
t176 = -pkin(3) * t209 + t266;
t37 = -pkin(4) * t119 - pkin(10) * t118 + t176;
t6 = -t124 * t473 + t136 * t472 + t367 * t25 + t362 * t37;
t2 = qJ(6) * t117 + qJD(6) * t223 + t6;
t7 = -qJD(5) * t54 - t25 * t362 + t367 * t37;
t4 = -pkin(5) * t117 + qJDD(6) - t7;
t641 = t2 * t367 + t362 * t4;
t642 = -t362 * t7 + t367 * t6;
t9 = pkin(5) * t72 - qJ(6) * t71 - qJD(6) * t205 + t26;
t662 = (t205 * t606 + t223 * t607) * qJD(5) / 0.2e1 + t362 * t668 + t26 * (mrSges(6,2) * t362 - t541) + (t220 + t163) * t573 + (-t415 / 0.2e1 + t411 / 0.2e1) * t474 + (-Ifges(7,3) * t367 + t525) * t583 + (Ifges(6,2) * t367 + t528) * t584 + t367 * t20 / 0.2e1 + Ifges(5,3) * t353 + t9 * t420 - t28 * mrSges(5,2) + t29 * mrSges(5,1) + (t449 + t644) * qJD(5) + t651 * t439 + t112 * t440 + t19 * t566 + t162 * t571 + (t362 * t629 + t367 * t656) * t582 + (t362 * t657 - t524 + t527) * t585 + Ifges(5,5) * t118 + Ifges(5,6) * t119 + (-t53 * t672 + t54 * t671 + t642) * mrSges(6,3) + (t42 * t672 + t43 * t671 + t641) * mrSges(7,2) + (-t530 + t617) * t633;
t132 = t202 * t363 + t192;
t650 = -pkin(3) * t476 + t132;
t240 = -t299 * t364 + t280;
t210 = t240 - t552;
t241 = t369 * t299 + t277;
t211 = -t269 + t241;
t561 = pkin(2) * t369;
t345 = pkin(3) + t561;
t489 = t364 * t368;
t649 = t368 * t210 - t211 * t363 + t345 * t476 + (t364 * t475 + (t363 * t369 + t489) * qJD(3)) * pkin(2);
t268 = pkin(5) * t473 - qJ(6) * t472 - qJD(6) * t362;
t648 = -pkin(5) * t619 + qJ(6) * t618 + t268;
t366 = sin(qJ(1));
t502 = t341 * t362;
t465 = mrSges(6,2) * t502;
t342 = cos(t356);
t499 = t342 * t366;
t647 = -t366 * t465 + t499 * t674;
t371 = cos(qJ(1));
t497 = t342 * t371;
t646 = -t371 * t465 + t497 * t674;
t645 = t674 * t341;
t354 = sin(t361);
t355 = cos(t361);
t422 = mrSges(5,1) * t341 + mrSges(5,2) * t342;
t643 = mrSges(4,1) * t354 + mrSges(4,2) * t355 + t422;
t601 = g(1) * t371 + g(2) * t366;
t33 = -mrSges(7,2) * t72 + mrSges(7,3) * t117;
t34 = mrSges(6,1) * t117 - mrSges(6,3) * t71;
t35 = -t117 * mrSges(7,1) + t71 * mrSges(7,2);
t36 = -mrSges(6,2) * t117 - mrSges(6,3) * t72;
t640 = t367 * (t36 + t33) + t362 * (-t34 + t35) + m(6) * ((-t362 * t54 - t367 * t53) * qJD(5) + t642) + m(7) * ((-t362 * t43 + t367 * t42) * qJD(5) + t641) + t667 * t473 + t666 * t472;
t635 = m(6) + m(7);
t634 = t303 / 0.2e1;
t565 = t370 / 0.2e1;
t555 = pkin(5) * t401;
t625 = t370 * Ifges(3,2);
t622 = t648 + t649;
t621 = t648 - t650;
t409 = pkin(5) * t362 - qJ(6) * t367;
t620 = -t409 * t434 - t129 + t268;
t219 = t401 * qJ(6);
t239 = t293 * t363 + t294 * t368;
t255 = -pkin(3) * t293 - t346;
t400 = t368 * t293 - t294 * t363;
t160 = -pkin(4) * t400 - pkin(10) * t239 + t255;
t249 = t369 * t324 + t325 * t364;
t221 = -pkin(9) * t294 + t249;
t250 = t364 * t324 - t369 * t325;
t222 = pkin(9) * t293 + t250;
t165 = t221 * t363 + t222 * t368;
t616 = t362 * t160 + t367 * t165;
t614 = t368 * t221 - t222 * t363;
t141 = t210 * t363 + t211 * t368;
t490 = t363 * t364;
t229 = t345 * t475 + (-t364 * t476 + (t368 * t369 - t490) * qJD(3)) * pkin(2);
t613 = t229 - t141;
t611 = t342 * mrSges(5,1) - t341 * mrSges(5,2);
t483 = t342 * pkin(4) + t341 * pkin(10);
t436 = t355 * mrSges(4,1) - mrSges(4,2) * t354;
t556 = pkin(4) * t341;
t559 = pkin(3) * t354;
t610 = m(7) * t559 - m(6) * (-t556 - t559) + t665;
t609 = t371 * t665 + t646;
t608 = t366 * t665 + t647;
t604 = t117 * t628 + t629 * t71 - t656 * t72;
t480 = qJD(1) * t370;
t481 = qJD(1) * t365;
t553 = pkin(7) * t370;
t554 = pkin(7) * t365;
t603 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t481) * t553 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t480) * t554;
t602 = t289 * t370 + t290 * t365;
t498 = t342 * t367;
t500 = t342 * t362;
t599 = -t498 * t632 + t500 * t631 - t611 + t645;
t596 = -m(3) * pkin(7) + m(4) * t372 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t458 = qJD(2) * t372;
t301 = t365 * t458;
t302 = t370 * t458;
t183 = t369 * t301 + t364 * t302 + t324 * t477 + t325 * t478;
t243 = -qJD(2) * t294 - t388;
t157 = pkin(9) * t243 + t183;
t184 = -qJD(3) * t250 - t301 * t364 + t369 * t302;
t242 = qJD(2) * t293 + t387;
t158 = -pkin(9) * t242 + t184;
t40 = qJD(4) * t614 + t157 * t368 + t158 * t363;
t143 = qJD(4) * t400 + t242 * t368 + t243 * t363;
t144 = qJD(4) * t239 + t242 * t363 - t368 * t243;
t479 = qJD(2) * t365;
t351 = pkin(2) * t479;
t228 = -pkin(3) * t243 + t351;
t49 = pkin(4) * t144 - pkin(10) * t143 + t228;
t13 = -qJD(5) * t616 - t362 * t40 + t367 * t49;
t595 = m(7) * pkin(5) + t632;
t322 = -mrSges(3,1) * t370 + mrSges(3,2) * t365;
t594 = -m(3) * pkin(1) - m(4) * t346 - mrSges(2,1) + t322 - t436 - t611;
t593 = -t436 + t599;
t592 = m(7) * qJ(6) - t631;
t175 = pkin(4) * t401 - pkin(10) * t434;
t589 = t7 * mrSges(6,1) - t4 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3);
t577 = t205 / 0.2e1;
t575 = t223 / 0.2e1;
t568 = t276 / 0.2e1;
t315 = pkin(10) * t499;
t564 = m(7) * t315;
t318 = pkin(10) * t497;
t563 = m(7) * t318;
t562 = pkin(2) * t365;
t560 = pkin(3) * t276;
t338 = pkin(3) * t355;
t558 = pkin(3) * t363;
t557 = pkin(3) * t368;
t547 = g(3) * t341;
t537 = mrSges(4,3) * t275;
t536 = mrSges(4,3) * t276;
t533 = Ifges(3,4) * t365;
t532 = Ifges(3,4) * t370;
t531 = Ifges(4,4) * t276;
t523 = t128 * mrSges(5,3);
t522 = t129 * mrSges(5,3);
t521 = t434 * mrSges(5,3);
t517 = t143 * t367;
t511 = t229 * t362;
t510 = t229 * t367;
t508 = t239 * t367;
t501 = t341 * t371;
t491 = t362 * t366;
t488 = t366 * t367;
t487 = t367 * t371;
t486 = t371 * t362;
t74 = t367 * t128 + t362 * t175;
t133 = t202 * t368 - t191;
t149 = t175 + t560;
t63 = t367 * t133 + t362 * t149;
t349 = pkin(2) * t481;
t142 = t149 + t349;
t65 = t367 * t141 + t362 * t142;
t271 = pkin(2) * t489 + t363 * t345;
t482 = t338 + t357;
t463 = pkin(3) * t475;
t459 = t338 + t483;
t457 = t239 * t472;
t438 = t471 / 0.2e1;
t298 = pkin(1) + t482;
t360 = -pkin(9) + t372;
t433 = t371 * t298 - t360 * t366;
t270 = -pkin(2) * t490 + t345 * t368;
t432 = t362 * t463;
t431 = t367 * t463;
t430 = pkin(5) * t498 + qJ(6) * t500 + t483;
t428 = -t366 * t556 + t315;
t427 = -pkin(4) * t501 + t318;
t425 = mrSges(3,1) * t365 + mrSges(3,2) * t370;
t416 = t533 + t625;
t413 = Ifges(3,5) * t370 - Ifges(3,6) * t365;
t406 = t362 * t42 + t367 * t43;
t405 = -t362 * t53 + t367 * t54;
t404 = t338 + t430;
t73 = -t128 * t362 + t175 * t367;
t62 = -t133 * t362 + t149 * t367;
t64 = -t141 * t362 + t142 * t367;
t84 = t160 * t367 - t165 * t362;
t396 = pkin(1) * t425;
t394 = t143 * t362 + t457;
t393 = t239 * t473 - t517;
t391 = t365 * (Ifges(3,1) * t370 - t533);
t12 = t160 * t472 - t165 * t473 + t362 * t49 + t367 * t40;
t41 = qJD(4) * t165 + t157 * t363 - t368 * t158;
t217 = Ifges(4,2) * t275 + Ifges(4,6) * t359 + t531;
t267 = Ifges(4,4) * t275;
t218 = t276 * Ifges(4,1) + t359 * Ifges(4,5) + t267;
t373 = (t522 + t670) * t401 + (t523 + t669) * t434 - (-Ifges(4,2) * t276 + t218 + t267) * t275 / 0.2e1 - t276 * (Ifges(4,1) * t275 - t531) / 0.2e1 + t323 * (mrSges(4,1) * t276 + mrSges(4,2) * t275) + Ifges(4,3) * t358 - t359 * (Ifges(4,5) * t275 - Ifges(4,6) * t276) / 0.2e1 + Ifges(4,5) * t208 + Ifges(4,6) * t209 + t148 * mrSges(4,1) - t147 * mrSges(4,2) + t234 * t536 + t233 * t537 + t217 * t568 + t662;
t348 = Ifges(3,4) * t480;
t344 = -pkin(4) - t557;
t305 = -t559 - t562;
t287 = t371 * t305;
t286 = t366 * t305;
t285 = t306 - t557;
t274 = Ifges(3,1) * t481 + Ifges(3,5) * qJD(2) + t348;
t273 = Ifges(3,6) * qJD(2) + qJD(1) * t416;
t264 = -pkin(4) - t270;
t262 = t342 * t487 + t491;
t261 = t342 * t486 - t488;
t260 = t342 * t488 - t486;
t259 = t342 * t491 + t487;
t253 = mrSges(4,1) * t359 - t536;
t252 = -mrSges(4,2) * t359 + t537;
t251 = t349 + t560;
t247 = -t270 + t306;
t232 = -mrSges(4,1) * t275 + mrSges(4,2) * t276;
t212 = -mrSges(5,2) * t437 + t521;
t188 = -mrSges(4,2) * t358 + mrSges(4,3) * t209;
t187 = mrSges(4,1) * t358 - mrSges(4,3) * t208;
t174 = -mrSges(5,1) * t434 + mrSges(5,2) * t401;
t138 = mrSges(7,1) * t204 - mrSges(7,3) * t205;
t137 = pkin(5) * t205 + qJ(6) * t204;
t99 = -mrSges(5,2) * t353 + mrSges(5,3) * t119;
t98 = mrSges(5,1) * t353 - mrSges(5,3) * t118;
t92 = t239 * t409 - t614;
t61 = pkin(5) * t400 - t84;
t60 = -qJ(6) * t400 + t616;
t51 = -t73 - t555;
t50 = t74 + t219;
t47 = -t64 - t555;
t46 = t219 + t65;
t45 = -t62 - t555;
t44 = t219 + t63;
t31 = mrSges(6,1) * t72 + mrSges(6,2) * t71;
t30 = mrSges(7,1) * t72 - mrSges(7,3) * t71;
t14 = t409 * t143 + (qJD(5) * t410 - qJD(6) * t367) * t239 + t41;
t11 = -pkin(5) * t144 - t13;
t10 = qJ(6) * t144 - qJD(6) * t400 + t12;
t1 = [(-mrSges(3,1) * t554 - mrSges(3,2) * t553 + 0.2e1 * Ifges(3,6) * t565) * qJDD(2) + (Ifges(3,1) * t304 + Ifges(3,4) * t634 + Ifges(3,5) * qJDD(2) - t438 * t625) * t365 + m(4) * (t147 * t250 + t148 * t249 + t183 * t234 + t184 * t233 - t266 * t346 - t323 * t351) + t457 * t676 + (t303 * t553 + t304 * t554 + t602) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t602) + (t274 * t565 + t413 * qJD(2) / 0.2e1 - t603) * qJD(2) + m(6) * (t12 * t54 + t13 * t53 + t6 * t616 + t7 * t84) + t616 * t36 + (t176 * mrSges(5,2) - t29 * mrSges(5,3) + Ifges(5,1) * t118 + Ifges(5,4) * t119 + Ifges(5,5) * t353 + t109 * t439 + t26 * t421 + t411 * t583 + t415 * t584 + t419 * t9 + t582 * t607 + t585 * t606 + (-t20 / 0.2e1 + t19 / 0.2e1 - t6 * mrSges(6,3) - t2 * mrSges(7,2)) * t362) * t239 + (-t657 * t393 + t394 * t673) * t577 + (-mrSges(4,1) * t266 + mrSges(4,3) * t147 + Ifges(4,4) * t208 + Ifges(4,2) * t209 + Ifges(4,6) * t358) * t293 + (t595 * t260 + t592 * t259 + (m(5) * t298 - t635 * (-t298 - t483) - t594 - t645) * t366 + (t596 + (m(5) + t635) * t360) * t371) * g(1) + (mrSges(4,2) * t266 - mrSges(4,3) * t148 + Ifges(4,1) * t208 + Ifges(4,4) * t209 + Ifges(4,5) * t358) * t294 + (-t393 * t629 - t394 * t656) * t575 + (t393 * t53 - t394 * t54 - t508 * t7) * mrSges(6,3) + t304 * t532 / 0.2e1 + (-Ifges(6,4) * t393 - Ifges(6,2) * t394) * t580 + (-m(5) * t433 + t674 * t501 - t635 * (pkin(4) * t497 + pkin(10) * t501 + t433) - t595 * t262 - t592 * t261 + t594 * t371 + t596 * t366) * g(2) - (-m(5) * t29 + m(6) * t26 + t31 - t98) * t614 + (-Ifges(7,5) * t393 + Ifges(7,3) * t394) * t579 + t232 * t351 + (-t393 * t42 - t394 * t43 + t4 * t508) * mrSges(7,2) + (Ifges(3,4) * t304 + Ifges(3,2) * t303) * t565 - t396 * t471 + (t370 * t532 + t391) * t438 + m(7) * (t10 * t43 + t11 * t42 + t14 * t55 + t2 * t60 + t4 * t61 + t9 * t92) - t273 * t479 / 0.2e1 + t508 * t668 + (-t604 / 0.2e1 - t176 * mrSges(5,1) + t28 * mrSges(5,3) + Ifges(5,4) * t118 + Ifges(5,2) * t119 + Ifges(5,6) * t353 - Ifges(6,6) * t584 - Ifges(7,6) * t583 - t582 * t628 - t585 * t629 - t589) * t400 + t416 * t634 - t322 * t518 + (-t523 + t449 + t220 / 0.2e1 + t163 / 0.2e1 + Ifges(5,1) * t571 + t675) * t143 + (t617 / 0.2e1 - t522 + Ifges(6,6) * t580 + Ifges(7,6) * t579 - t626 / 0.2e1 - t162 / 0.2e1 - Ifges(5,4) * t571 + t629 * t577 + t628 * t575 - t587) * t144 + (-t233 * t242 + t234 * t243) * mrSges(4,3) + m(5) * (t129 * t40 + t165 * t28 + t176 * t255 + t228 * t248) + Ifges(2,3) * qJDD(1) + t651 * (t239 * t440 + t517 / 0.2e1) + t664 * t41 + t359 * (Ifges(4,5) * t242 + Ifges(4,6) * t243) / 0.2e1 - t346 * (-mrSges(4,1) * t209 + mrSges(4,2) * t208) - t323 * (-mrSges(4,1) * t243 + mrSges(4,2) * t242) - pkin(1) * (-mrSges(3,1) * t303 + mrSges(3,2) * t304) + t275 * (Ifges(4,4) * t242 + Ifges(4,2) * t243) / 0.2e1 + t255 * (-mrSges(5,1) * t119 + mrSges(5,2) * t118) + t243 * t217 / 0.2e1 + t249 * t187 + t250 * t188 + t183 * t252 + t184 * t253 + t242 * t218 / 0.2e1 + t55 * (mrSges(7,1) * t394 + mrSges(7,3) * t393) + t123 * (mrSges(6,1) * t394 - mrSges(6,2) * t393) + t228 * t174 + t40 * t212 + t165 * t99 + t10 * t150 + t12 * t151 + t13 * t152 + t11 * t153 + t14 * t138 + t60 * t33 + t61 * t35 + (Ifges(4,1) * t242 + Ifges(4,4) * t243) * t568 + t84 * t34 + t92 * t30; (-t511 - t64) * t152 + (-t65 + t510) * t151 + (-t46 + t510) * t150 + (-t47 + t511) * t153 - m(4) * (t233 * t240 + t234 * t241 - t323 * t349) + (t603 + (t396 - t391 / 0.2e1) * qJD(1)) * qJD(1) + (-m(5) * t482 - m(6) * (t357 + t459) - m(7) * (t357 + t404) + t322 - m(4) * t357 + t593) * g(3) + (t229 * t405 + t26 * t264 - t53 * t64 - t54 * t65 - g(2) * (t286 + t428) - g(1) * (t287 + t427)) * m(6) + (t129 * t613 - t248 * t251 + t270 * t29 + t271 * t28) * m(5) + t608 * g(2) + t609 * g(1) - t232 * t349 - t413 * t471 / 0.2e1 + t613 * t212 + t622 * t138 + (t229 * t406 + t247 * t9 - t42 * t47 - t43 * t46 - g(1) * (t287 + t318) - g(2) * (t286 + t315) + t622 * t55) * m(7) - (-Ifges(3,2) * t481 + t274 + t348) * t480 / 0.2e1 + t273 * t481 / 0.2e1 + t664 * t649 + Ifges(3,3) * qJDD(2) + Ifges(3,6) * t303 + Ifges(3,5) * t304 - t289 * mrSges(3,2) - t290 * mrSges(3,1) + t270 * t98 + t271 * t99 + t264 * t31 + t247 * t30 - t251 * t174 - t241 * t252 - t240 * t253 + t640 * (pkin(10) + t271) + t601 * (m(4) * t562 - m(5) * t305 + t425 + t643) + t187 * t561 + (t252 * t477 - t253 * t478 + m(4) * (t147 * t364 + t148 * t369 + (-t233 * t364 + t234 * t369) * qJD(3)) + t364 * t188) * pkin(2) + t373; (t463 - t133) * t212 + (t432 - t45) * t153 + (-t62 - t432) * t152 + (t431 - t63) * t151 + (t431 - t44) * t150 + (-m(5) * t338 - m(6) * t459 - m(7) * t404 + t593) * g(3) + (t366 * t610 - t564 + t647) * g(2) + (t371 * t610 - t563 + t646) * g(1) - t174 * t560 + t621 * t138 + (t285 * t9 + t406 * t463 - t42 * t45 - t43 * t44 + t55 * t621) * m(7) + (t26 * t344 + (t123 * t363 + t368 * t405) * qJD(4) * pkin(3) - g(1) * t318 - g(2) * t315 - t123 * t132 - t53 * t62 - t54 * t63) * m(6) + t344 * t31 + t285 * t30 - t233 * t252 + t234 * t253 + t98 * t557 + t99 * t558 + (t128 * t132 - t129 * t133 - t248 * t560 + (t28 * t363 + t29 * t368 + (-t128 * t363 + t129 * t368) * qJD(4)) * pkin(3)) * m(5) + t373 + t650 * t615 + (m(5) * t559 + t643) * t601 + t640 * (pkin(10) + t558); t669 * t434 + t601 * t422 + (t520 + t591) * t129 + (-t564 + t608) * g(2) + (-t563 + t609) * g(1) + (t306 * t9 - t42 * t51 - t43 * t50 + t55 * t620) * m(7) + t620 * t138 + (-t212 + t521) * t128 + (-m(6) * t483 - m(7) * t430 + t599) * g(3) + t670 * t401 + (-pkin(4) * t26 - g(1) * t427 - g(2) * t428 - t53 * t73 - t54 * t74) * m(6) + t306 * t30 + t640 * pkin(10) - t50 * t150 - t74 * t151 - t73 * t152 - t51 * t153 - pkin(4) * t31 + t662; t589 + (t419 + t421) * t547 + (-t204 * t629 - t205 * t656) * t576 + (t259 * t632 + t260 * t631) * g(2) + (t261 * t632 + t262 * t631) * g(1) + (-t204 * t657 + t109 + t198 - t529) * t578 + (-Ifges(6,2) * t205 - t199 + t651) * t579 + (-m(7) * t42 + t534 - t666) * t54 + t604 + (Ifges(7,3) * t205 - t526) * t580 + t112 * t577 + (t409 * t547 - g(2) * (-pkin(5) * t259 + qJ(6) * t260) - g(1) * (-pkin(5) * t261 + qJ(6) * t262) - pkin(5) * t4 + qJ(6) * t2 + qJD(6) * t43 - t137 * t55) * m(7) - t55 * (mrSges(7,1) * t205 + mrSges(7,3) * t204) - t123 * (mrSges(6,1) * t205 - mrSges(6,2) * t204) + qJD(6) * t150 - t137 * t138 + qJ(6) * t33 - pkin(5) * t35 + (-m(7) * t43 - t535 + t667) * t53 + t43 * t538 + t42 * t539; t205 * t138 - t223 * t150 + (-g(1) * t261 - g(2) * t259 - g(3) * t502 + t55 * t205 - t43 * t223 + t4) * m(7) + t35;];
tau  = t1;
