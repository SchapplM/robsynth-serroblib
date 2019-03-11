% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:14:17
% EndTime: 2019-03-09 17:15:21
% DurationCPUTime: 38.59s
% Computational Cost: add. (9866->867), mult. (21166->1075), div. (0->0), fcn. (14080->8), ass. (0->386)
t657 = Ifges(4,1) + Ifges(5,1);
t656 = Ifges(5,4) + Ifges(4,5);
t654 = Ifges(4,6) - Ifges(5,6);
t326 = sin(qJ(2));
t313 = t326 * pkin(8);
t330 = cos(qJ(2));
t318 = t330 * pkin(2);
t424 = -pkin(1) - t318;
t357 = t424 - t313;
t232 = t357 * qJD(1);
t445 = qJD(1) * t330;
t307 = pkin(7) * t445;
t274 = qJD(2) * pkin(8) + t307;
t325 = sin(qJ(3));
t329 = cos(qJ(3));
t159 = t329 * t232 - t325 * t274;
t573 = qJD(4) - t159;
t665 = -mrSges(6,3) - mrSges(7,3);
t621 = Ifges(6,4) + Ifges(7,4);
t547 = pkin(8) - pkin(9);
t276 = t547 * t329;
t423 = -pkin(7) * t325 - pkin(3);
t458 = t329 * t330;
t344 = -pkin(9) * t458 + (-pkin(4) + t423) * t326;
t387 = pkin(2) * t326 - pkin(8) * t330;
t260 = t387 * qJD(1);
t468 = t260 * t329;
t664 = -qJD(1) * t344 + qJD(3) * t276 + t468;
t224 = t325 * t260;
t446 = qJD(1) * t326;
t296 = qJ(4) * t446;
t441 = qJD(3) * t325;
t461 = t326 * t329;
t462 = t325 * t330;
t663 = t224 + t296 + (-pkin(7) * t461 + pkin(9) * t462) * qJD(1) + t547 * t441;
t323 = -qJ(6) - pkin(9);
t561 = -m(7) * t323 - t665;
t658 = mrSges(5,2) + mrSges(4,3);
t662 = m(6) * pkin(9) + t561 - t658;
t420 = t329 * t446;
t251 = qJD(2) * t325 + t420;
t661 = -pkin(9) * t251 + t573;
t288 = qJD(3) - t445;
t282 = qJD(5) - t288;
t522 = t282 / 0.2e1;
t422 = t325 * t446;
t436 = t329 * qJD(2);
t250 = t422 - t436;
t324 = sin(qJ(5));
t328 = cos(qJ(5));
t360 = t250 * t324 + t251 * t328;
t538 = t360 / 0.2e1;
t155 = t250 * t328 - t251 * t324;
t541 = t155 / 0.2e1;
t617 = Ifges(7,3) + Ifges(6,3);
t618 = Ifges(6,6) + Ifges(7,6);
t620 = Ifges(6,5) + Ifges(7,5);
t660 = t617 * t522 + t620 * t538 + t618 * t541;
t435 = qJD(1) * qJD(2);
t265 = qJDD(1) * t326 + t330 * t435;
t442 = qJD(3) * t250;
t148 = qJDD(2) * t325 + t265 * t329 - t442;
t546 = t148 / 0.2e1;
t149 = qJD(3) * t251 - t329 * qJDD(2) + t265 * t325;
t544 = t149 / 0.2e1;
t264 = t330 * qJDD(1) - t326 * t435;
t248 = qJDD(3) - t264;
t530 = t248 / 0.2e1;
t528 = t250 / 0.2e1;
t659 = -t251 / 0.2e1;
t521 = -t288 / 0.2e1;
t623 = Ifges(6,1) + Ifges(7,1);
t655 = Ifges(5,2) + Ifges(4,3);
t619 = Ifges(6,2) + Ifges(7,2);
t463 = t325 * t328;
t466 = t324 * t329;
t359 = -t463 + t466;
t571 = qJD(3) - qJD(5);
t164 = t571 * t359;
t206 = t359 * t330;
t190 = qJD(1) * t206;
t653 = t190 - t164;
t590 = -t325 * t654 + t329 * t656;
t493 = Ifges(5,5) * t325;
t495 = Ifges(4,4) * t325;
t588 = t329 * t657 + t493 - t495;
t380 = t329 * mrSges(5,1) + t325 * mrSges(5,3);
t382 = mrSges(4,1) * t329 - mrSges(4,2) * t325;
t652 = -t382 - t380;
t439 = qJD(3) * t329;
t651 = -t325 * qJD(4) - t307 + (t329 * t445 - t439) * qJ(4);
t327 = sin(qJ(1));
t331 = cos(qJ(1));
t575 = g(1) * t331 + g(2) * t327;
t545 = -t149 / 0.2e1;
t629 = t264 / 0.2e1;
t649 = t265 / 0.2e1;
t639 = t621 * t155;
t608 = t282 * t620 + t360 * t623 + t639;
t648 = t608 / 0.2e1;
t637 = t621 * t360;
t609 = t155 * t619 + t282 * t618 + t637;
t647 = t609 / 0.2e1;
t646 = t656 * t530 + (-Ifges(4,4) + Ifges(5,5)) * t544 + t657 * t546;
t231 = qJDD(5) - t248;
t40 = qJD(5) * t155 + t148 * t328 + t149 * t324;
t41 = -qJD(5) * t360 - t148 * t324 + t149 * t328;
t645 = -t619 * t41 / 0.2e1 - t621 * t40 / 0.2e1 - t618 * t231 / 0.2e1;
t275 = t547 * t325;
t173 = t324 * t275 + t328 * t276;
t611 = -qJD(5) * t173 + t324 * t663 + t328 * t664;
t437 = qJD(5) * t328;
t438 = qJD(5) * t324;
t610 = t275 * t437 - t276 * t438 + t324 * t664 - t328 * t663;
t160 = t325 * t232 + t329 * t274;
t106 = pkin(9) * t250 + t160;
t332 = -pkin(3) - pkin(4);
t267 = t328 * qJ(4) + t324 * t332;
t602 = -qJD(5) * t267 - t328 * t106 - t324 * t661;
t266 = -qJ(4) * t324 + t328 * t332;
t601 = qJD(5) * t266 - t324 * t106 + t328 * t661;
t642 = qJ(6) * t155;
t456 = t331 * t325;
t221 = -t327 * t329 + t330 * t456;
t457 = t330 * t331;
t222 = t325 * t327 + t329 * t457;
t361 = t221 * t324 + t222 * t328;
t625 = mrSges(6,2) + mrSges(7,2);
t641 = t625 * t361;
t241 = Ifges(4,4) * t250;
t483 = t250 * Ifges(5,5);
t640 = t251 * t657 + t288 * t656 - t241 + t483;
t421 = t325 * t445;
t427 = t332 * t325;
t597 = qJD(3) * t427 - t332 * t421 - t651;
t204 = t324 * t461 - t326 * t463;
t467 = t324 * t325;
t358 = t328 * t329 + t467;
t599 = t358 * t326;
t398 = t204 * mrSges(7,1) + mrSges(7,2) * t599;
t399 = t204 * mrSges(6,1) + mrSges(6,2) * t599;
t638 = -t398 - t399;
t81 = t288 * t332 + t661;
t278 = t288 * qJ(4);
t88 = t106 + t278;
t31 = -t324 * t88 + t328 * t81;
t600 = qJ(6) * t360;
t17 = t31 - t600;
t16 = pkin(5) * t282 + t17;
t306 = pkin(7) * t446;
t273 = -qJD(2) * pkin(2) + t306;
t123 = t250 * pkin(3) - t251 * qJ(4) + t273;
t96 = -pkin(4) * t250 - t123;
t58 = -pkin(5) * t155 + qJD(6) + t96;
t636 = t96 * mrSges(6,2) + t58 * mrSges(7,2) - mrSges(6,3) * t31 - mrSges(7,3) * t16 + t648;
t32 = t324 * t81 + t328 * t88;
t18 = t32 + t642;
t497 = Ifges(3,4) * t326;
t606 = t330 * Ifges(3,2);
t374 = t497 + t606;
t634 = t31 * mrSges(6,1) + t16 * mrSges(7,1) - t32 * mrSges(6,2) - t18 * mrSges(7,2) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t374 / 0.2e1 + t655 * t521 + t656 * t659 + t654 * t528 + t660;
t552 = m(7) * pkin(5);
t632 = -m(6) - m(5);
t631 = -m(7) - m(6);
t539 = -t360 / 0.2e1;
t628 = pkin(5) * t360;
t443 = qJD(2) * t330;
t430 = pkin(7) * t443;
t627 = mrSges(7,1) + mrSges(6,1);
t626 = -mrSges(4,2) + mrSges(5,3);
t624 = -mrSges(3,3) + mrSges(2,2);
t616 = t231 * t620 + t40 * t623 + t41 * t621;
t29 = -mrSges(7,2) * t231 + mrSges(7,3) * t41;
t30 = -mrSges(6,2) * t231 + mrSges(6,3) * t41;
t615 = t29 + t30;
t614 = qJ(6) * t653 - qJD(6) * t358 + t610;
t163 = t571 * t358;
t207 = t358 * t330;
t191 = qJD(1) * t207;
t613 = pkin(5) * t446 + qJD(6) * t359 + t611 + (-t163 + t191) * qJ(6);
t607 = -pkin(5) * t653 + t597;
t605 = -t642 + t602;
t604 = -t600 + t601;
t384 = mrSges(3,1) * t330 - mrSges(3,2) * t326;
t603 = -t384 - mrSges(2,1);
t351 = t326 * t359;
t449 = t318 + t313;
t269 = -pkin(1) - t449;
t289 = pkin(7) * t462;
t317 = t330 * pkin(3);
t507 = pkin(9) * t326;
t147 = pkin(4) * t330 + t289 + t317 + (-t269 - t507) * t329;
t290 = pkin(7) * t458;
t189 = t325 * t269 + t290;
t174 = -qJ(4) * t330 + t189;
t464 = t325 * t326;
t158 = pkin(9) * t464 + t174;
t72 = t324 * t147 + t328 * t158;
t596 = (-t421 + t441) * pkin(3) + t651;
t595 = -qJD(2) * mrSges(3,1) + mrSges(4,1) * t250 + mrSges(4,2) * t251 + mrSges(3,3) * t446;
t594 = t326 * t655 + t330 * t590;
t593 = t326 * t656 + t330 * t588;
t477 = t329 * mrSges(5,3);
t379 = t325 * mrSges(5,1) - t477;
t381 = mrSges(4,1) * t325 + mrSges(4,2) * t329;
t592 = -t123 * t379 - t273 * t381;
t591 = t325 * t656 + t329 * t654;
t492 = Ifges(5,5) * t329;
t494 = Ifges(4,4) * t329;
t589 = t325 * t657 - t492 + t494;
t585 = t148 * t656 - t149 * t654 + t248 * t655;
t584 = t231 * t617 + t40 * t620 + t41 * t618;
t240 = Ifges(5,5) * t251;
t124 = t288 * Ifges(5,6) + t250 * Ifges(5,3) + t240;
t304 = Ifges(3,4) * t445;
t583 = Ifges(3,1) * t446 + Ifges(3,5) * qJD(2) + t325 * t124 + t304;
t246 = t264 * pkin(7);
t247 = t265 * pkin(7);
t582 = t246 * t330 + t247 * t326;
t459 = t327 * t330;
t219 = t325 * t459 + t329 * t331;
t220 = t327 * t458 - t456;
t581 = t219 * t328 - t220 * t324;
t580 = t658 * t326;
t473 = qJDD(1) * pkin(1);
t162 = -pkin(2) * t264 - pkin(8) * t265 - t473;
t213 = qJDD(2) * pkin(8) + t246;
t59 = t325 * t162 + t329 * t213 + t232 * t439 - t274 * t441;
t60 = t162 * t329 - t325 * t213 - t232 * t441 - t274 * t439;
t579 = -t325 * t60 + t329 * t59;
t42 = t248 * qJ(4) + t288 * qJD(4) + t59;
t349 = qJDD(4) - t60;
t47 = -pkin(3) * t248 + t349;
t578 = t325 * t47 + t329 * t42;
t117 = t278 + t160;
t577 = -t117 * mrSges(5,2) - t160 * mrSges(4,3);
t115 = -pkin(3) * t288 + t573;
t576 = t115 * mrSges(5,2) - t159 * mrSges(4,3);
t572 = -m(5) + t631;
t570 = -t552 - t627;
t312 = t325 * qJ(4);
t268 = -t329 * pkin(3) - pkin(2) - t312;
t568 = -Ifges(5,5) * t148 / 0.2e1 - Ifges(5,6) * t248 / 0.2e1 + Ifges(4,4) * t546 + Ifges(4,6) * t530 + (Ifges(5,3) + Ifges(4,2)) * t545;
t472 = t219 * t324;
t130 = t220 * t328 + t472;
t564 = t625 * t130;
t401 = pkin(5) * t324 + qJ(4);
t562 = -m(7) * t401 - t626;
t301 = pkin(5) * t328 + pkin(4);
t559 = m(6) * pkin(4) + m(7) * t301 + mrSges(4,1) + mrSges(5,1);
t558 = -t96 * mrSges(6,1) - t58 * mrSges(7,1) + mrSges(6,3) * t32 + mrSges(7,3) * t18;
t21 = -pkin(9) * t148 + t248 * t332 + t349;
t23 = pkin(9) * t149 + t42;
t4 = -qJD(5) * t32 + t328 * t21 - t23 * t324;
t1 = pkin(5) * t231 - qJ(6) * t40 - qJD(6) * t360 + t4;
t3 = t324 * t21 + t328 * t23 + t81 * t437 - t438 * t88;
t2 = qJ(6) * t41 + qJD(6) * t155 + t3;
t556 = t4 * mrSges(6,1) + t1 * mrSges(7,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t555 = -t60 * mrSges(4,1) + t47 * mrSges(5,1) + t59 * mrSges(4,2) - t42 * mrSges(5,3);
t551 = t40 / 0.2e1;
t550 = t41 / 0.2e1;
t542 = -t155 / 0.2e1;
t531 = t231 / 0.2e1;
t529 = -t250 / 0.2e1;
t526 = t251 / 0.2e1;
t523 = -t282 / 0.2e1;
t516 = m(7) * t326;
t510 = pkin(7) * t326;
t504 = -qJD(1) / 0.2e1;
t503 = qJD(3) / 0.2e1;
t502 = -pkin(3) - t301;
t501 = mrSges(6,3) * t155;
t500 = mrSges(6,3) * t360;
t499 = mrSges(7,3) * t155;
t498 = mrSges(7,3) * t360;
t496 = Ifges(3,4) * t330;
t482 = t251 * Ifges(4,4);
t460 = t326 * t331;
t263 = t387 * qJD(2);
t453 = t325 * t263 + t269 * t439;
t167 = t251 * pkin(3) + t250 * qJ(4);
t419 = t330 * t436;
t451 = qJ(4) * t419 + qJD(4) * t461;
t448 = t331 * pkin(1) + t327 * pkin(7);
t444 = qJD(2) * t326;
t440 = qJD(3) * t326;
t431 = pkin(7) * t444;
t429 = pkin(8) * t441;
t428 = pkin(8) * t439;
t214 = -qJDD(2) * pkin(2) + t247;
t418 = t325 * t440;
t127 = -t250 * Ifges(4,2) + t288 * Ifges(4,6) + t482;
t417 = -t325 * t127 / 0.2e1;
t14 = -t41 * mrSges(7,1) + t40 * mrSges(7,2);
t409 = -t445 / 0.2e1;
t406 = t443 / 0.2e1;
t403 = -t440 / 0.2e1;
t402 = t439 / 0.2e1;
t400 = t435 / 0.2e1;
t92 = -t248 * mrSges(5,1) + t148 * mrSges(5,2);
t71 = t328 * t147 - t158 * t324;
t172 = t328 * t275 - t276 * t324;
t188 = t269 * t329 - t289;
t239 = t329 * pkin(4) - t268;
t393 = pkin(2) * t457 + pkin(8) * t460 + t448;
t389 = -pkin(7) + t427;
t112 = -pkin(4) * t251 - t167;
t388 = t423 * t326;
t385 = qJD(3) * t290 - t263 * t329 + t269 * t441;
t383 = mrSges(3,1) * t326 + mrSges(3,2) * t330;
t373 = -Ifges(4,2) * t325 + t494;
t372 = Ifges(4,2) * t329 + t495;
t369 = Ifges(3,5) * t330 - Ifges(3,6) * t326;
t366 = Ifges(5,3) * t325 + t492;
t365 = -Ifges(5,3) * t329 + t493;
t362 = t222 * pkin(3) + t393;
t138 = t221 * t328 - t222 * t324;
t177 = -pkin(7) * t420 + t224;
t44 = t149 * pkin(3) - t148 * qJ(4) - t251 * qJD(4) + t214;
t356 = pkin(1) * t383;
t355 = t329 * t332 - t312;
t352 = t326 * (Ifges(3,1) * t330 - t497);
t75 = pkin(9) * t418 + qJD(2) * t344 + t385;
t298 = qJ(4) * t444;
t76 = t298 + (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t461 + (-qJD(4) + (-pkin(7) * qJD(3) + pkin(9) * qJD(2)) * t325) * t330 + t453;
t12 = t147 * t437 - t158 * t438 + t324 * t75 + t328 * t76;
t285 = qJ(4) * t461;
t171 = t326 * t389 + t285;
t348 = -t418 + t419;
t347 = t325 * t443 + t326 * t439;
t345 = -g(1) * t221 - g(2) * t219 - g(3) * t464;
t24 = -pkin(4) * t149 - t44;
t340 = Ifges(4,6) * t326 + t330 * t373;
t339 = Ifges(5,6) * t326 + t330 * t366;
t13 = -qJD(5) * t72 - t324 * t76 + t328 * t75;
t101 = (-t326 * t436 - t330 * t441) * pkin(7) + t453;
t336 = t556 + t584;
t82 = t355 * t440 + t389 * t443 + t451;
t319 = t331 * pkin(7);
t271 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t445;
t259 = -pkin(5) + t266;
t234 = t381 * t326;
t211 = t221 * pkin(3);
t209 = t219 * pkin(3);
t198 = -t285 + (pkin(3) * t325 + pkin(7)) * t326;
t181 = -mrSges(5,2) * t250 + mrSges(5,3) * t288;
t180 = -mrSges(5,1) * t288 + mrSges(5,2) * t251;
t179 = mrSges(4,1) * t288 - mrSges(4,3) * t251;
t178 = -mrSges(4,2) * t288 - mrSges(4,3) * t250;
t176 = pkin(7) * t422 + t468;
t175 = -t188 + t317;
t170 = pkin(5) * t358 + t239;
t168 = mrSges(5,1) * t250 - mrSges(5,3) * t251;
t166 = qJD(1) * t388 - t468;
t165 = t177 + t296;
t119 = -qJ(6) * t358 + t173;
t118 = qJ(6) * t359 + t172;
t111 = pkin(5) * t204 + t171;
t110 = mrSges(6,1) * t282 - t500;
t109 = mrSges(7,1) * t282 - t498;
t108 = -mrSges(6,2) * t282 + t501;
t107 = -mrSges(7,2) * t282 + t499;
t102 = t325 * t431 - t385;
t100 = pkin(3) * t347 + qJ(4) * t418 + t430 - t451;
t95 = qJD(2) * t388 + t385;
t94 = -mrSges(5,2) * t149 + mrSges(5,3) * t248;
t93 = -mrSges(4,2) * t248 - mrSges(4,3) * t149;
t91 = mrSges(4,1) * t248 - mrSges(4,3) * t148;
t87 = -qJD(4) * t330 + t101 + t298;
t84 = qJD(2) * t207 + t351 * t571;
t83 = -qJD(2) * t206 + t163 * t326;
t78 = -mrSges(6,1) * t155 + mrSges(6,2) * t360;
t77 = -mrSges(7,1) * t155 + mrSges(7,2) * t360;
t74 = mrSges(4,1) * t149 + mrSges(4,2) * t148;
t73 = mrSges(5,1) * t149 - mrSges(5,3) * t148;
t67 = t112 - t628;
t51 = -qJ(6) * t204 + t72;
t48 = pkin(5) * t330 - qJ(6) * t599 + t71;
t33 = -pkin(5) * t83 + t82;
t28 = mrSges(6,1) * t231 - mrSges(6,3) * t40;
t27 = mrSges(7,1) * t231 - mrSges(7,3) * t40;
t15 = -mrSges(6,1) * t41 + mrSges(6,2) * t40;
t7 = -pkin(5) * t41 + qJDD(6) + t24;
t6 = qJ(6) * t83 - qJD(6) * t204 + t12;
t5 = -pkin(5) * t444 - qJ(6) * t84 - qJD(6) * t599 + t13;
t8 = [(qJD(2) * t339 - t365 * t440) * t528 + (qJD(2) * t340 - t372 * t440) * t529 + t417 * t443 + (t472 * t552 + t572 * (-t220 * pkin(3) - qJ(4) * t219 + t319) + t624 * t331 + (-m(4) - m(3)) * t319 + t559 * t220 + t626 * t219 + t627 * t130 + t625 * t581 + (m(3) * pkin(1) + t631 * t424 + (-m(4) - m(5)) * t357 + (-m(7) * (-pkin(8) - t323) + m(6) * t547 + t665) * t326 + t580 - t603) * t327) * g(1) + (t496 * t400 + Ifges(3,6) * qJDD(2) - Ifges(5,6) * t544 - Ifges(4,6) * t545 + pkin(7) * (-qJDD(2) * mrSges(3,2) + mrSges(3,3) * t264) + t620 * t551 + t618 * t550 - t656 * t546 + t617 * t531 - t655 * t530 + t555 + t556 + Ifges(3,4) * t649 + Ifges(3,2) * t629 + t584 / 0.2e1 - t585 / 0.2e1) * t330 + (t159 * mrSges(4,1) - t115 * mrSges(5,1) - t160 * mrSges(4,2) + t117 * mrSges(5,3) - t634 - t660) * t444 + (-qJDD(2) * mrSges(3,1) + t74) * t510 + t640 * (t325 * t403 + t329 * t406) + t204 * t645 + t461 * t646 + t83 * t647 + t84 * t648 + t496 * t649 + (m(4) * t214 * pkin(7) + Ifges(3,1) * t265 + Ifges(3,4) * t629 + Ifges(3,5) * qJDD(2) + t124 * t402 + t366 * t544 + t373 * t545 + t44 * t379 - t400 * t606 + t530 * t590 + t546 * t588) * t326 + (t621 * t83 + t623 * t84) * t538 + t48 * t27 + t51 * t29 + (-t159 * t348 - t160 * t347 - t461 * t60) * mrSges(4,3) + (t115 * t348 - t117 * t347 + t461 * t47) * mrSges(5,2) - t271 * t431 + t352 * t400 + (-t42 * mrSges(5,2) - t59 * mrSges(4,3) - t568) * t464 + (t618 * t83 + t620 * t84) * t522 + (t619 * t83 + t621 * t84) * t541 + t273 * (mrSges(4,1) * t347 + mrSges(4,2) * t348) + t123 * (mrSges(5,1) * t347 - mrSges(5,3) * t348) + (-t204 * t619 + t599 * t621) * t550 + (-t204 * t621 + t599 * t623) * t551 + (-t1 * t599 - t16 * t84 + t18 * t83 - t2 * t204) * mrSges(7,3) + t616 * t599 / 0.2e1 + (-t204 * t618 + t599 * t620) * t531 + (-t204 * t3 - t31 * t84 + t32 * t83 - t4 * t599) * mrSges(6,3) + t384 * t473 - t356 * t435 + t374 * t629 + (-m(3) * t448 - m(4) * t393 - m(7) * t362 + t632 * (qJ(4) * t221 + t362) + t603 * t331 + t624 * t327 - t559 * t222 + t562 * t221 - t627 * t361 - t625 * t138 + t662 * t460) * g(2) + (t265 * t510 + t582) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t582) + t583 * t406 + t7 * t398 + t24 * t399 + Ifges(2,3) * qJDD(1) + m(4) * (t101 * t160 + t102 * t159 + t188 * t60 + t189 * t59 + t273 * t430) + qJD(2) ^ 2 * t369 / 0.2e1 + (qJD(2) * t593 - t440 * t589) * t526 + (qJD(2) * t594 - t440 * t591) * t288 / 0.2e1 + t595 * t430 - pkin(1) * (-mrSges(3,1) * t264 + mrSges(3,2) * t265) + t214 * t234 + m(7) * (t1 * t48 + t111 * t7 + t16 * t5 + t18 * t6 + t2 * t51 + t33 * t58) + m(6) * (t12 * t32 + t13 * t31 + t171 * t24 + t3 * t72 + t4 * t71 + t82 * t96) + m(5) * (t100 * t123 + t115 * t95 + t117 * t87 + t174 * t42 + t175 * t47 + t198 * t44) + t71 * t28 + t72 * t30 + t33 * t77 + t82 * t78 + t58 * (-mrSges(7,1) * t83 + mrSges(7,2) * t84) + t329 * t127 * t403 + t96 * (-mrSges(6,1) * t83 + mrSges(6,2) * t84) + t6 * t107 + t12 * t108 + t5 * t109 + t13 * t110 + t111 * t14 + t100 * t168 + t171 * t15 + t174 * t94 + t175 * t92 + t101 * t178 + t102 * t179 + t95 * t180 + t87 * t181 + t188 * t91 + t189 * t93 + t198 * t73; (-pkin(2) * t214 - t159 * t176 - t160 * t177 - t273 * t307) * m(4) + (-t115 * t166 - t117 * t165 + t596 * t123 + t268 * t44) * m(5) + (-Ifges(3,2) * t409 - t523 * t617 - t539 * t620 - t542 * t618 + t634) * t446 + (t304 + t583) * t409 + (-t373 / 0.2e1 + t366 / 0.2e1) * t442 + (t428 - t166) * t180 + (-t429 - t177) * t178 + t640 * (t329 * t409 + t402) + t358 * t645 + t325 * t646 + t190 * t647 + ((t340 / 0.2e1 - t339 / 0.2e1) * t250 - t160 * (-mrSges(4,2) * t326 - mrSges(4,3) * t462) - t117 * (-mrSges(5,2) * t462 + mrSges(5,3) * t326) - t159 * (mrSges(4,1) * t326 - mrSges(4,3) * t458) - t115 * (-mrSges(5,1) * t326 + mrSges(5,2) * t458) + (t356 - t352 / 0.2e1) * qJD(1)) * qJD(1) + (-t428 - t176) * t179 + t613 * t109 + (t1 * t118 + t119 * t2 + t16 * t613 + t170 * t7 + t18 * t614 + t58 * t607) * m(7) + t614 * t107 + (t503 * t588 + t504 * t593) * t251 - t58 * (mrSges(7,1) * t190 + mrSges(7,2) * t191) - t96 * (mrSges(6,1) * t190 + mrSges(6,2) * t191) + (-t429 - t165) * t181 + t568 * t329 + t127 * t421 / 0.2e1 - t616 * t359 / 0.2e1 + t7 * (mrSges(7,1) * t358 - mrSges(7,2) * t359) + t24 * (mrSges(6,1) * t358 - mrSges(6,2) * t359) + (-t358 * t621 - t359 * t623) * t551 + (-t358 * t618 - t359 * t620) * t531 + (-t358 * t619 - t359 * t621) * t550 + (-t190 * t621 + t191 * t623) * t539 + (-t190 * t618 + t191 * t620) * t523 + (-t190 * t619 + t191 * t621) * t542 + (t1 * t359 + t16 * t191 + t18 * t190 - t2 * t358) * mrSges(7,3) + (t190 * t32 + t191 * t31 - t3 * t358 + t359 * t4) * mrSges(6,3) - t369 * t435 / 0.2e1 + (t383 - (-t325 * t401 + t329 * t502 - pkin(2)) * t516 + (-m(5) * t268 - m(6) * (-pkin(2) + t355) + m(4) * pkin(2) - t652) * t326 + t662 * t330 - t351 * t625 + t599 * t627) * t575 + t576 * t439 + (t124 / 0.2e1 + t577) * t441 + t578 * mrSges(5,2) + t579 * mrSges(4,3) + (t503 * t590 + t504 * t594) * t288 - (t618 * t522 + t621 * t538 + t619 * t541 + t558 + t647) * t164 + (t620 * t522 + t623 * t538 + t621 * t541 + t636) * t163 + ((t92 - t91) * t325 + (t94 + t93) * t329 + ((-t159 * t329 - t160 * t325) * qJD(3) + t579) * m(4) + ((t115 * t329 - t117 * t325) * qJD(3) + t578) * m(5) + (g(1) * t457 + g(2) * t459) * (-m(4) + t572)) * pkin(8) + (-m(4) * t449 - m(6) * (pkin(4) * t458 - t507) - t384 + t572 * (pkin(3) * t458 + qJ(4) * t462 + t449) + (-m(7) * (pkin(5) * t467 + t301 * t329) + t652) * t330 + t561 * t326 - t627 * t207 + t625 * t206 - t580) * g(3) + Ifges(3,3) * qJDD(2) - t214 * t382 - t44 * t380 + t589 * t546 + t591 * t530 + (t417 - t592) * qJD(3) + t592 * t445 - t595 * t307 + t596 * t168 + t597 * t78 + Ifges(3,5) * t265 + t268 * t73 + Ifges(3,6) * t264 - t246 * mrSges(3,2) - t247 * mrSges(3,1) + t239 * t15 + t607 * t77 - t608 * t191 / 0.2e1 + t610 * t108 + t611 * t110 + (t172 * t4 + t173 * t3 + t239 * t24 + t31 * t611 + t32 * t610 + t597 * t96) * m(6) + t365 * t544 + t372 * t545 - pkin(2) * t74 + t271 * t306 + t118 * t27 + t119 * t29 + t170 * t14 + t172 * t28 + t173 * t30; t127 * t526 + (t273 * mrSges(4,2) - t123 * mrSges(5,3) - t521 * t656 + t576) * t250 + (t179 - t180) * t160 + (-t178 - t181) * t159 + (-(pkin(5) * t466 + t325 * t502) * t516 + t234 + (-m(6) * t427 - t477 - (-m(5) * pkin(3) - mrSges(5,1)) * t325) * t326 + t572 * t285 + t638) * g(3) + (m(7) * t209 + t632 * (qJ(4) * t220 - t209) + t562 * t220 + t559 * t219 + t627 * t581 - t564) * g(2) - t555 + (-Ifges(4,2) * t251 - t241 + t640) * t528 + (m(7) * t211 + t632 * (qJ(4) * t222 - t211) + t562 * t222 + t559 * t221 + t627 * t138 - t641) * g(1) - t483 * t529 + t601 * t108 + (-t112 * t96 + t266 * t4 + t267 * t3 + t31 * t602 + t32 * t601) * m(6) + t602 * t110 + t604 * t107 + t605 * t109 + (t1 * t259 + t16 * t605 + t18 * t604 + t2 * t267 - t58 * t67) * m(7) + (-pkin(3) * t47 + qJ(4) * t42 - t115 * t160 + t117 * t573 - t123 * t167) * m(5) - t336 + t615 * t267 + (t618 * t523 + t621 * t539 + t619 * t542 - t558) * t360 + t585 + t609 * t539 + (-t620 * t523 - t623 * t539 - t621 * t542 + t636) * t155 + (-t273 * mrSges(4,1) - t123 * mrSges(5,1) + Ifges(5,3) * t529 - t521 * t654 - t577) * t251 + t266 * t28 + t259 * t27 - t67 * t77 - pkin(3) * t92 + qJ(4) * t94 - t112 * t78 + (-t250 * t657 + t124 + t240 - t482) * t659 - t167 * t168 + qJD(4) * t181; -t288 * t181 + (t168 - t77 - t78) * t251 + (t27 + t28 + t282 * (t107 + t108)) * t328 + (-t282 * (t109 + t110) + t615) * t324 + t92 + (t1 * t328 + t2 * t324 - t251 * t58 + t345 + t282 * (-t16 * t324 + t18 * t328)) * m(7) + (-t251 * t96 + t3 * t324 + t328 * t4 + t345 + t282 * (-t31 * t324 + t32 * t328)) * m(6) + (-t117 * t288 + t123 * t251 + t345 + t47) * m(5); t16 * t499 + pkin(5) * t27 + t336 - t17 * t107 + t1 * t552 - t58 * (mrSges(7,1) * t360 + mrSges(7,2) * t155) - t96 * (mrSges(6,1) * t360 + mrSges(6,2) * t155) + (t155 * t623 - t637) * t539 + t609 * t538 + (t155 * t620 - t360 * t618) * t523 + (-m(7) * t58 - t77) * t628 + (t500 + t110) * t32 + (t501 - t108) * t31 + (t498 - m(7) * (-t16 + t17) + t109) * t18 + (pkin(5) * t359 * t516 - t638) * g(3) + (t570 * t581 + t564) * g(2) + (t138 * t570 + t641) * g(1) + (-t360 * t619 + t608 + t639) * t542; -t155 * t107 + t360 * t109 + (-g(3) * t330 - t18 * t155 + t16 * t360 + t326 * t575 + t7) * m(7) + t14;];
tau  = t8;
