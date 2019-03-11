% Calculate vector of inverse dynamics joint torques for
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPPR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR6_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR6_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR6_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR6_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:47:57
% EndTime: 2019-03-09 15:49:11
% DurationCPUTime: 49.74s
% Computational Cost: add. (18425->966), mult. (44771->1264), div. (0->0), fcn. (35724->14), ass. (0->436)
t363 = cos(qJ(2));
t359 = sin(qJ(2));
t353 = sin(pkin(6));
t471 = qJD(1) * t353;
t439 = t359 * t471;
t355 = cos(pkin(6));
t470 = qJD(1) * t355;
t458 = pkin(1) * t470;
t277 = -pkin(8) * t439 + t363 * t458;
t386 = (pkin(2) * t359 - pkin(9) * t363) * t353;
t278 = qJD(1) * t386;
t358 = sin(qJ(3));
t362 = cos(qJ(3));
t190 = -t358 * t277 + t362 * t278;
t356 = -qJ(4) - pkin(9);
t429 = qJD(3) * t356;
t475 = t362 * t363;
t674 = -(pkin(3) * t359 - qJ(4) * t475) * t471 - t190 - qJD(4) * t358 + t362 * t429;
t191 = t362 * t277 + t358 * t278;
t438 = t363 * t471;
t419 = t358 * t438;
t673 = -qJ(4) * t419 - qJD(4) * t362 - t358 * t429 + t191;
t352 = sin(pkin(11));
t354 = cos(pkin(11));
t482 = t354 * t362;
t224 = -t352 * t419 + t438 * t482;
t465 = qJD(3) * t362;
t466 = qJD(3) * t358;
t288 = -t352 * t466 + t354 * t465;
t672 = t224 - t288;
t334 = qJD(2) + t470;
t249 = t334 * t362 - t358 * t439;
t250 = t334 * t358 + t362 * t439;
t393 = t249 * t352 + t354 * t250;
t551 = t393 / 0.2e1;
t176 = -t354 * t249 + t250 * t352;
t555 = -t176 / 0.2e1;
t671 = -Ifges(6,2) * t551 - Ifges(6,6) * t555;
t311 = -qJD(3) + t438;
t542 = t311 / 0.2e1;
t650 = Ifges(6,1) + Ifges(5,3);
t670 = Ifges(4,3) + t650;
t299 = t352 * t362 + t354 * t358;
t223 = t299 * t438;
t287 = t299 * qJD(3);
t657 = t223 - t287;
t280 = pkin(8) * t438 + t359 * t458;
t613 = -t280 + (-t419 + t466) * pkin(3);
t357 = sin(qJ(6));
t361 = cos(qJ(6));
t231 = pkin(9) * t334 + t280;
t266 = (-pkin(2) * t363 - pkin(9) * t359 - pkin(1)) * t353;
t235 = qJD(1) * t266;
t163 = -t231 * t358 + t362 * t235;
t132 = -qJ(4) * t250 + t163;
t125 = -pkin(3) * t311 + t132;
t164 = t231 * t362 + t235 * t358;
t133 = qJ(4) * t249 + t164;
t129 = t352 * t133;
t65 = t125 * t354 - t129;
t413 = qJD(5) - t65;
t580 = pkin(4) + pkin(10);
t629 = pkin(5) * t393;
t40 = t311 * t580 + t413 + t629;
t230 = -t334 * pkin(2) - t277;
t174 = -t249 * pkin(3) + qJD(4) + t230;
t365 = -qJ(5) * t393 + t174;
t56 = t176 * t580 + t365;
t15 = -t357 * t56 + t361 * t40;
t16 = t357 * t40 + t361 * t56;
t61 = pkin(4) * t311 + t413;
t74 = t176 * pkin(4) + t365;
t669 = -mrSges(6,1) * t61 - mrSges(7,1) * t15 - mrSges(5,2) * t174 + mrSges(7,2) * t16 + mrSges(5,3) * t65 + mrSges(6,3) * t74 - Ifges(6,4) * t542 + t671;
t543 = -t311 / 0.2e1;
t552 = -t393 / 0.2e1;
t554 = t176 / 0.2e1;
t483 = t354 * t133;
t66 = t352 * t125 + t483;
t62 = qJ(5) * t311 - t66;
t668 = mrSges(6,1) * t62 - mrSges(5,3) * t66 + mrSges(5,1) * t174 - mrSges(6,2) * t74 + Ifges(6,5) * t543 + Ifges(5,6) * t542 + (Ifges(5,2) + Ifges(6,3)) * t554 + (Ifges(5,4) + Ifges(6,6)) * t552;
t631 = m(7) + m(6);
t462 = qJD(1) * qJD(2);
t284 = (qJDD(1) * t359 + t363 * t462) * t353;
t460 = qJDD(1) * t355;
t333 = qJDD(2) + t460;
t166 = qJD(3) * t249 + t284 * t362 + t333 * t358;
t167 = -qJD(3) * t250 - t284 * t358 + t333 * t362;
t109 = t166 * t352 - t354 * t167;
t576 = -t109 / 0.2e1;
t110 = t166 * t354 + t167 * t352;
t573 = t110 / 0.2e1;
t283 = (-qJDD(1) * t363 + t359 * t462) * t353;
t270 = qJDD(3) + t283;
t546 = t270 / 0.2e1;
t627 = -Ifges(6,4) + Ifges(5,5);
t626 = Ifges(6,5) - Ifges(5,6);
t623 = t673 * t352 + t354 * t674;
t667 = qJ(5) * t672 - qJD(5) * t299 + t613;
t628 = mrSges(5,2) - mrSges(6,3);
t392 = -qJ(5) * t631 + t628;
t408 = mrSges(7,1) * t357 + mrSges(7,2) * t361;
t602 = t392 - t408;
t574 = -t110 / 0.2e1;
t666 = t573 - t574;
t575 = t109 / 0.2e1;
t665 = t575 - t576;
t173 = qJD(6) + t393;
t556 = t173 / 0.2e1;
t141 = t176 * t357 - t311 * t361;
t563 = t141 / 0.2e1;
t140 = t176 * t361 + t311 * t357;
t565 = t140 / 0.2e1;
t664 = Ifges(5,1) * t551 + Ifges(5,4) * t555 + Ifges(7,5) * t563 - Ifges(6,2) * t552 - Ifges(6,6) * t554 + Ifges(7,6) * t565 + Ifges(7,3) * t556 + t543 * t627 - t669;
t662 = -Ifges(5,4) * t551 - Ifges(5,2) * t555 + Ifges(6,6) * t552 + Ifges(6,3) * t554 + t668;
t521 = -mrSges(6,2) + mrSges(5,1);
t646 = Ifges(5,1) * t393 - t176 * Ifges(5,4) - t311 * Ifges(5,5) + t141 * Ifges(7,5) + t140 * Ifges(7,6) + t173 * Ifges(7,3);
t487 = t353 * t359;
t427 = t580 * t487;
t660 = -pkin(5) * t672 + qJD(1) * t427 - t623;
t620 = t352 * t674 - t673 * t354;
t659 = t580 * t657 - t667;
t535 = pkin(1) * t355;
t457 = qJD(2) * t535;
t461 = qJDD(1) * t353;
t644 = pkin(8) * t461 + qJD(1) * t457;
t645 = -pkin(8) * t353 * t462 + pkin(1) * t460;
t200 = t359 * t645 + t363 * t644;
t184 = pkin(9) * t333 + t200;
t189 = -pkin(1) * t461 + pkin(2) * t283 - pkin(9) * t284;
t73 = -qJD(3) * t164 - t184 * t358 + t362 * t189;
t39 = pkin(3) * t270 - qJ(4) * t166 - qJD(4) * t250 + t73;
t72 = t362 * t184 + t358 * t189 - t231 * t466 + t235 * t465;
t48 = qJ(4) * t167 + qJD(4) * t249 + t72;
t13 = -t352 * t48 + t354 * t39;
t415 = qJDD(5) - t13;
t11 = -pkin(4) * t270 + t415;
t201 = -t359 * t644 + t363 * t645;
t185 = -t333 * pkin(2) - t201;
t119 = -t167 * pkin(3) + qJDD(4) + t185;
t364 = -t110 * qJ(5) - qJD(5) * t393 + t119;
t20 = t109 * pkin(4) + t364;
t107 = qJDD(6) + t110;
t577 = t107 / 0.2e1;
t55 = -qJD(6) * t141 + t109 * t361 - t270 * t357;
t587 = t55 / 0.2e1;
t54 = qJD(6) * t140 + t109 * t357 + t270 * t361;
t588 = t54 / 0.2e1;
t12 = t109 * t580 + t364;
t5 = pkin(5) * t110 - t270 * t580 + t415;
t1 = qJD(6) * t15 + t12 * t361 + t357 * t5;
t2 = -qJD(6) * t16 - t12 * t357 + t361 * t5;
t605 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t630 = -t270 / 0.2e1;
t7 = Ifges(7,5) * t54 + Ifges(7,6) * t55 + Ifges(7,3) * t107;
t653 = t605 + mrSges(6,1) * t11 + mrSges(5,2) * t119 - mrSges(5,3) * t13 - mrSges(6,3) * t20 + 0.2e1 * Ifges(5,1) * t573 + 0.2e1 * Ifges(5,4) * t576 + Ifges(6,4) * t630 + Ifges(7,5) * t588 + Ifges(7,6) * t587 + Ifges(7,3) * t577 - t665 * Ifges(6,6) + t666 * Ifges(6,2) + t7 / 0.2e1 + (t627 + Ifges(5,5)) * t546;
t651 = pkin(5) * t176;
t289 = t355 * t362 - t358 * t487;
t391 = t289 * pkin(3);
t540 = cos(qJ(1));
t441 = t540 * t359;
t360 = sin(qJ(1));
t477 = t360 * t363;
t292 = t355 * t441 + t477;
t442 = t353 * t540;
t373 = t292 * t358 + t362 * t442;
t369 = t373 * pkin(3);
t647 = t200 * mrSges(3,2);
t621 = qJ(5) * t439 - t620;
t440 = t540 * t363;
t478 = t359 * t360;
t294 = -t355 * t478 + t440;
t485 = t353 * t362;
t220 = -t294 * t358 + t360 * t485;
t641 = Ifges(4,5) * t166 + Ifges(4,6) * t167 + t626 * t109 + t627 * t110 + t270 * t670;
t640 = -Ifges(5,4) * t552 - Ifges(5,2) * t554 + Ifges(6,6) * t551 + Ifges(6,3) * t555 + t542 * t626 - t668;
t557 = -t173 / 0.2e1;
t564 = -t141 / 0.2e1;
t566 = -t140 / 0.2e1;
t639 = Ifges(5,1) * t552 + Ifges(5,4) * t554 + Ifges(7,5) * t564 + Ifges(7,6) * t566 + Ifges(7,3) * t557 + t542 * t627 + t669 + t671;
t14 = t352 * t39 + t354 * t48;
t10 = -qJ(5) * t270 + qJD(5) * t311 - t14;
t637 = -t73 * mrSges(4,1) - t13 * mrSges(5,1) + t72 * mrSges(4,2) + t14 * mrSges(5,2) - t11 * mrSges(6,2) + t10 * mrSges(6,3);
t634 = mrSges(5,1) * t119 + mrSges(6,1) * t10 - mrSges(6,2) * t20 - mrSges(5,3) * t14 - Ifges(5,4) * t666 + Ifges(6,5) * t546 + Ifges(5,2) * t665 + Ifges(5,6) * t630 + 0.2e1 * Ifges(6,6) * t574 + 0.2e1 * Ifges(6,3) * t575;
t632 = m(7) * pkin(10);
t559 = t166 / 0.2e1;
t558 = t167 / 0.2e1;
t298 = t352 * t358 - t482;
t346 = pkin(3) * t362 + pkin(2);
t389 = -qJ(5) * t299 - t346;
t179 = t298 * t580 + t389;
t315 = t356 * t358;
t316 = t356 * t362;
t228 = -t354 * t315 - t316 * t352;
t194 = pkin(5) * t299 + t228;
t122 = -t179 * t357 + t194 * t361;
t625 = qJD(6) * t122 + t357 * t660 - t361 * t659;
t123 = t179 * t361 + t194 * t357;
t624 = -qJD(6) * t123 + t357 * t659 + t361 * t660;
t147 = mrSges(6,1) * t176 + mrSges(6,3) * t311;
t88 = -mrSges(7,1) * t140 + mrSges(7,2) * t141;
t501 = -t147 + t88;
t622 = pkin(4) * t439 - t623;
t619 = pkin(5) * t657 - t621;
t618 = -pkin(4) * t657 + t667;
t148 = mrSges(6,1) * t393 - mrSges(6,2) * t311;
t150 = -mrSges(5,1) * t311 - mrSges(5,3) * t393;
t617 = -t148 + t150;
t426 = mrSges(3,3) * t439;
t616 = -mrSges(3,1) * t334 - mrSges(4,1) * t249 + mrSges(4,2) * t250 + t426;
t196 = t223 * t361 - t357 * t439;
t464 = qJD(6) * t357;
t497 = t287 * t361;
t383 = t298 * t464 - t497;
t615 = t196 + t383;
t197 = t223 * t357 + t361 * t439;
t463 = qJD(6) * t361;
t480 = t357 * t287;
t384 = t298 * t463 + t480;
t614 = t197 - t384;
t484 = t353 * t363;
t297 = pkin(8) * t484 + t359 * t535;
t265 = pkin(9) * t355 + t297;
t188 = t362 * t265 + t358 * t266;
t290 = t355 * t358 + t359 * t485;
t198 = -t354 * t289 + t290 * t352;
t479 = t357 * t363;
t319 = t353 * t479;
t171 = t198 * t361 + t319;
t335 = pkin(8) * t487;
t534 = pkin(1) * t363;
t296 = t355 * t534 - t335;
t25 = mrSges(7,1) * t107 - mrSges(7,3) * t54;
t26 = -mrSges(7,2) * t107 + mrSges(7,3) * t55;
t612 = t361 * t25 + t357 * t26;
t611 = -t358 * t73 + t362 * t72;
t609 = -mrSges(7,3) - t632;
t351 = qJ(3) + pkin(11);
t348 = sin(t351);
t349 = cos(t351);
t411 = -mrSges(4,1) * t362 + mrSges(4,2) * t358;
t608 = -m(4) * pkin(2) + t348 * t628 - t349 * t521 + t411;
t279 = qJD(2) * t386;
t281 = t296 * qJD(2);
t127 = -qJD(3) * t188 + t362 * t279 - t281 * t358;
t468 = qJD(2) * t353;
t436 = t363 * t468;
t218 = qJD(3) * t289 + t362 * t436;
t467 = qJD(2) * t359;
t437 = t353 * t467;
t78 = pkin(3) * t437 - qJ(4) * t218 - qJD(4) * t290 + t127;
t126 = -t265 * t466 + t266 * t465 + t358 * t279 + t362 * t281;
t217 = -qJD(3) * t290 - t358 * t436;
t87 = qJ(4) * t217 + qJD(4) * t289 + t126;
t30 = t352 * t78 + t354 * t87;
t27 = -t353 * (qJ(5) * t467 - qJD(5) * t363) - t30;
t409 = mrSges(7,1) * t361 - mrSges(7,2) * t357;
t41 = -t62 - t651;
t512 = t141 * Ifges(7,4);
t59 = t140 * Ifges(7,2) + t173 * Ifges(7,6) + t512;
t586 = -t59 / 0.2e1;
t607 = t361 * t586 + t41 * t409;
t606 = t521 + t632;
t603 = -m(4) * pkin(9) - m(7) * pkin(5) - mrSges(6,1) - mrSges(4,3) - mrSges(5,3);
t210 = t292 * t348 + t349 * t442;
t486 = t353 * t360;
t214 = t294 * t348 - t349 * t486;
t262 = t348 * t487 - t355 * t349;
t601 = g(1) * t214 + g(2) * t210 + g(3) * t262;
t600 = t348 * t408 + mrSges(3,1) - t608;
t599 = -t65 * mrSges(5,1) + t66 * mrSges(5,2) - t61 * mrSges(6,2) + t62 * mrSges(6,3);
t597 = -t409 + t603;
t596 = mrSges(3,2) + t597;
t593 = t353 ^ 2;
t592 = Ifges(7,1) * t588 + Ifges(7,4) * t587 + Ifges(7,5) * t577;
t585 = t59 / 0.2e1;
t139 = Ifges(7,4) * t140;
t60 = t141 * Ifges(7,1) + t173 * Ifges(7,5) + t139;
t584 = -t60 / 0.2e1;
t583 = t60 / 0.2e1;
t582 = Ifges(4,4) * t559 + Ifges(4,2) * t558 + Ifges(4,6) * t546;
t581 = Ifges(4,1) * t559 + Ifges(4,4) * t558 + Ifges(4,5) * t546;
t579 = pkin(1) * mrSges(3,1);
t578 = pkin(1) * mrSges(3,2);
t508 = t250 * Ifges(4,4);
t157 = t249 * Ifges(4,2) - t311 * Ifges(4,6) + t508;
t561 = t157 / 0.2e1;
t236 = Ifges(4,4) * t249;
t158 = t250 * Ifges(4,1) - t311 * Ifges(4,5) + t236;
t560 = t158 / 0.2e1;
t547 = t250 / 0.2e1;
t541 = t355 / 0.2e1;
t533 = pkin(3) * t250;
t532 = pkin(3) * t352;
t531 = pkin(3) * t354;
t530 = t1 * t357;
t520 = mrSges(7,3) * t361;
t519 = Ifges(3,4) * t359;
t518 = Ifges(3,4) * t363;
t517 = Ifges(4,4) * t358;
t516 = Ifges(4,4) * t362;
t515 = Ifges(7,4) * t357;
t514 = Ifges(7,4) * t361;
t513 = Ifges(3,6) * t334;
t511 = t163 * mrSges(4,3);
t510 = t164 * mrSges(4,3);
t509 = t249 * Ifges(4,6);
t507 = t250 * Ifges(4,5);
t506 = t334 * Ifges(3,5);
t500 = qJ(5) * t348;
t499 = t393 * t357;
t291 = -t355 * t440 + t478;
t496 = t291 * t349;
t293 = t355 * t477 + t441;
t494 = t293 * t349;
t492 = t298 * t357;
t491 = t298 * t361;
t476 = t361 * t363;
t187 = -t358 * t265 + t362 * t266;
t138 = -pkin(3) * t484 - t290 * qJ(4) + t187;
t151 = qJ(4) * t289 + t188;
t82 = t352 * t138 + t354 * t151;
t474 = -t291 * t346 - t292 * t356;
t473 = -t293 * t346 - t294 * t356;
t282 = pkin(8) * t436 + t359 * t457;
t472 = t540 * pkin(1) + pkin(8) * t486;
t337 = pkin(4) * t484;
t453 = qJ(5) * t484;
t451 = t358 * t486;
t449 = t353 * t476;
t444 = Ifges(3,5) * t284 - Ifges(3,6) * t283 + Ifges(3,3) * t333;
t345 = -pkin(4) - t531;
t433 = -t471 / 0.2e1;
t431 = -t464 / 0.2e1;
t430 = -pkin(1) * t360 + pkin(8) * t442;
t29 = -t352 * t87 + t354 * t78;
t84 = t110 * mrSges(6,1) + t270 * mrSges(6,2);
t68 = t132 * t352 + t483;
t81 = t354 * t138 - t352 * t151;
t211 = t292 * t349 - t348 * t442;
t320 = t358 * t442;
t428 = -t292 * t362 + t320;
t425 = mrSges(3,3) * t438;
t422 = -pkin(4) * t496 - t291 * t500 + t474;
t421 = -pkin(4) * t494 - t293 * t500 + t473;
t417 = t363 * t433;
t416 = t220 * pkin(3);
t183 = -pkin(3) * t217 + t282;
t77 = t337 - t81;
t414 = -t521 + t609;
t412 = mrSges(4,1) * t289 - mrSges(4,2) * t290;
t406 = Ifges(4,1) * t362 - t517;
t405 = Ifges(7,1) * t357 + t514;
t404 = Ifges(3,2) * t363 + t519;
t403 = -Ifges(4,2) * t358 + t516;
t402 = Ifges(7,2) * t361 + t515;
t401 = Ifges(4,5) * t362 - Ifges(4,6) * t358;
t400 = Ifges(7,5) * t357 + Ifges(7,6) * t361;
t399 = qJ(5) * t176 + t533;
t397 = t15 * t357 - t16 * t361;
t199 = t289 * t352 + t290 * t354;
t57 = t199 * pkin(5) + pkin(10) * t484 + t77;
t264 = t335 + (-pkin(2) - t534) * t355;
t202 = t264 - t391;
t367 = -t199 * qJ(5) + t202;
t71 = t198 * t580 + t367;
t23 = -t357 * t71 + t361 * t57;
t24 = t357 * t57 + t361 * t71;
t98 = -mrSges(7,2) * t173 + mrSges(7,3) * t140;
t99 = mrSges(7,1) * t173 - mrSges(7,3) * t141;
t396 = -t357 * t99 + t361 * t98;
t395 = -t357 * t98 - t361 * t99;
t394 = pkin(3) * t451 - t293 * t356 + t294 * t346 + t472;
t69 = t132 * t354 - t129;
t229 = t315 * t352 - t316 * t354;
t76 = t453 - t82;
t385 = -t357 * t198 + t449;
t380 = pkin(3) * t320 + t291 * t356 - t292 * t346 + t430;
t379 = t230 * (mrSges(4,1) * t358 + mrSges(4,2) * t362);
t371 = mrSges(3,2) + t603;
t370 = -g(1) * t293 - g(2) * t291 + g(3) * t484;
t146 = t217 * t352 + t218 * t354;
t368 = -qJ(5) * t146 - qJD(5) * t199 + t183;
t366 = -qJD(6) * t397 + t2 * t361 + t530;
t343 = qJ(5) + t532;
t325 = Ifges(3,4) * t438;
t304 = t346 * t484;
t295 = (-mrSges(3,1) * t363 + mrSges(3,2) * t359) * t353;
t276 = -t334 * mrSges(3,2) + t425;
t226 = Ifges(3,1) * t439 + t325 + t506;
t225 = t404 * t471 + t513;
t221 = t294 * t362 + t451;
t215 = t294 * t349 + t348 * t486;
t205 = -mrSges(4,1) * t311 - mrSges(4,3) * t250;
t204 = mrSges(4,2) * t311 + mrSges(4,3) * t249;
t203 = pkin(4) * t298 + t389;
t195 = -pkin(5) * t298 + t229;
t169 = t214 * t357 + t293 * t361;
t168 = t214 * t361 - t293 * t357;
t156 = -t311 * Ifges(4,3) + t507 + t509;
t149 = mrSges(5,2) * t311 - mrSges(5,3) * t176;
t145 = -t354 * t217 + t218 * t352;
t135 = -mrSges(4,2) * t270 + mrSges(4,3) * t167;
t134 = mrSges(4,1) * t270 - mrSges(4,3) * t166;
t121 = -mrSges(6,2) * t176 - mrSges(6,3) * t393;
t120 = mrSges(5,1) * t176 + mrSges(5,2) * t393;
t117 = -mrSges(4,1) * t167 + mrSges(4,2) * t166;
t115 = -t311 * Ifges(6,1) - Ifges(6,4) * t393 + t176 * Ifges(6,5);
t112 = Ifges(5,5) * t393 - t176 * Ifges(5,6) - t311 * Ifges(5,3);
t108 = t198 * pkin(4) + t367;
t102 = t110 * mrSges(6,3);
t101 = t110 * mrSges(5,2);
t95 = pkin(4) * t393 + t399;
t92 = qJD(6) * t171 + t357 * t145 + t361 * t437;
t91 = qJD(6) * t385 + t361 * t145 - t357 * t437;
t86 = mrSges(5,1) * t270 - mrSges(5,3) * t110;
t85 = -mrSges(5,2) * t270 - mrSges(5,3) * t109;
t83 = mrSges(6,1) * t109 - mrSges(6,3) * t270;
t64 = t393 * t580 + t399;
t63 = -t198 * pkin(5) - t76;
t51 = t69 - t629;
t50 = t68 - t651;
t49 = pkin(4) * t145 + t368;
t45 = -t109 * mrSges(6,2) - t102;
t44 = t109 * mrSges(5,1) + t101;
t33 = t145 * t580 + t368;
t28 = -pkin(4) * t437 - t29;
t22 = -t145 * pkin(5) - t27;
t21 = pkin(5) * t146 - qJD(2) * t427 - t29;
t19 = t357 * t50 + t361 * t64;
t18 = -t357 * t64 + t361 * t50;
t17 = -mrSges(7,1) * t55 + mrSges(7,2) * t54;
t8 = Ifges(7,4) * t54 + Ifges(7,2) * t55 + Ifges(7,6) * t107;
t6 = -pkin(5) * t109 - t10;
t4 = -qJD(6) * t24 + t21 * t361 - t33 * t357;
t3 = qJD(6) * t23 + t21 * t357 + t33 * t361;
t9 = [(Ifges(7,4) * t92 + Ifges(7,2) * t91) * t565 + (Ifges(7,1) * t92 + Ifges(7,4) * t91) * t563 + t634 * t198 + (Ifges(4,4) * t290 + Ifges(4,2) * t289) * t558 + (Ifges(4,1) * t290 + Ifges(4,4) * t289) * t559 + t662 * t145 + (t646 / 0.2e1 + t664) * t146 + (-m(5) * t380 - m(4) * (-pkin(2) * t292 + t430) - t428 * mrSges(4,1) - t373 * mrSges(4,2) - m(3) * t430 + t292 * mrSges(3,1) - mrSges(3,3) * t442 + t360 * mrSges(2,1) + t540 * mrSges(2,2) - t414 * t211 - t602 * t210 + (-t371 + t409) * t291 + t631 * (pkin(4) * t211 - t380)) * g(1) + (-Ifges(7,4) * t385 + Ifges(7,2) * t171) * t587 + (-Ifges(7,1) * t385 + Ifges(7,4) * t171) * t588 - t385 * t592 + (-Ifges(7,5) * t385 + Ifges(7,6) * t171) * t577 + (t1 * t171 - t15 * t92 + t16 * t91 + t2 * t385) * mrSges(7,3) + t6 * (-mrSges(7,1) * t171 - mrSges(7,2) * t385) + m(7) * (t1 * t24 + t15 * t4 + t16 * t3 + t2 * t23 + t22 * t41 + t6 * t63) + m(5) * (t119 * t202 + t13 * t81 + t14 * t82 + t174 * t183 + t29 * t65 + t30 * t66) + m(4) * (t126 * t164 + t127 * t163 + t185 * t264 + t187 * t73 + t188 * t72 + t230 * t282) + t290 * t581 + t289 * t582 + t92 * t583 + t91 * t585 + (Ifges(4,1) * t218 + Ifges(4,4) * t217) * t547 + t218 * t560 + t217 * t561 + t653 * t199 + (Ifges(3,5) * t541 - t296 * mrSges(3,3) + (t359 * Ifges(3,1) + t518 - t578) * t353) * t284 - (Ifges(3,6) * t541 + t297 * mrSges(3,3) + (t404 + t579) * t353) * t283 + (Ifges(4,5) * t290 + Ifges(4,6) * t289 + t198 * t626) * t546 + (Ifges(3,3) * t541 + t296 * mrSges(3,1) - t297 * mrSges(3,2) + (Ifges(3,5) * t359 + Ifges(3,6) * t363) * t353) * t333 + t249 * (Ifges(4,4) * t218 + Ifges(4,2) * t217) / 0.2e1 + (Ifges(7,5) * t92 + Ifges(7,6) * t91) * t556 + t444 * t541 + (-t163 * t218 + t164 * t217 + t289 * t72 - t290 * t73) * mrSges(4,3) + m(6) * (t10 * t76 + t108 * t20 + t11 * t77 + t27 * t62 + t28 * t61 + t49 * t74) + t23 * t25 + t24 * t26 + (Ifges(4,5) * t218 + Ifges(4,6) * t217 + t145 * t626) * t543 - t185 * t412 + t281 * t276 + (-m(3) * t472 - t294 * mrSges(3,1) - m(4) * (pkin(2) * t294 + t472) - t221 * mrSges(4,1) - t220 * mrSges(4,2) - t169 * mrSges(7,1) - t168 * mrSges(7,2) - m(5) * t394 - t540 * mrSges(2,1) + (-mrSges(3,3) * t353 + mrSges(2,2)) * t360 + t392 * t214 + t414 * t215 + t371 * t293 - t631 * (t215 * pkin(4) + t394)) * g(2) + t264 * t117 + t230 * (-mrSges(4,1) * t217 + mrSges(4,2) * t218) + t201 * (mrSges(3,1) * t355 - mrSges(3,3) * t487) + (-pkin(1) * t295 * t353 + Ifges(2,3)) * qJDD(1) + t127 * t205 + t202 * t44 + t126 * t204 + t183 * t120 + t187 * t134 + t188 * t135 + t63 * t17 + ((t226 / 0.2e1 - t277 * mrSges(3,3) + t506 / 0.2e1 + (-t578 + t518 / 0.2e1) * t471) * t363 + (-t225 / 0.2e1 + t156 / 0.2e1 + t112 / 0.2e1 + t115 / 0.2e1 + t507 / 0.2e1 - t164 * mrSges(4,2) + t163 * mrSges(4,1) + t509 / 0.2e1 - t280 * mrSges(3,3) - t513 / 0.2e1 + (Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1) * t393 + (Ifges(6,5) / 0.2e1 - Ifges(5,6) / 0.2e1) * t176 + (-Ifges(4,3) / 0.2e1 - Ifges(5,3) / 0.2e1 - Ifges(6,1) / 0.2e1) * t311 + (-t579 - t519 / 0.2e1 + (Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1) * t363) * t471 - t599) * t359) * t468 + t76 * t83 + t77 * t84 + t82 * t85 + t81 * t86 - t355 * t647 + t22 * t88 + t41 * (-mrSges(7,1) * t91 + mrSges(7,2) * t92) + t3 * t98 + t4 * t99 + t108 * t45 + (-t641 / 0.2e1 + mrSges(3,3) * t200 - Ifges(6,4) * t574 - Ifges(4,5) * t559 - Ifges(5,5) * t573 - Ifges(6,5) * t575 - Ifges(4,6) * t558 - Ifges(5,6) * t576 - t670 * t546 + t637) * t484 + t49 * t121 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t593 + t200 * t297 + t201 * t296 - t277 * t282 + t280 * t281) + t27 * t147 + t28 * t148 + t30 * t149 + t29 * t150 + t171 * t8 / 0.2e1 + t616 * t282; (Ifges(7,1) * t384 - Ifges(7,4) * t383) * t563 + (Ifges(7,1) * t197 + Ifges(7,4) * t196) * t564 + (t560 - t511) * t465 + (-t164 * (-mrSges(4,3) * t358 * t363 - mrSges(4,2) * t359) - t163 * (mrSges(4,1) * t359 - mrSges(4,3) * t475)) * t471 + (t400 * t577 + t402 * t587 + t405 * t588 - t409 * t6 + t431 * t59 + t546 * t626 + t634) * t298 + (-pkin(2) * t185 - t163 * t190 - t164 * t191) * m(4) + t646 * (-t224 / 0.2e1 + t288 / 0.2e1) + (t362 * t158 + t226 + t325) * t417 + (t84 - t86) * t228 + (t425 - t276) * t277 + (t543 * t626 + t662) * t287 + t664 * t288 + (-t359 * (Ifges(3,1) * t363 - t519) / 0.2e1 + pkin(1) * (mrSges(3,1) * t359 + mrSges(3,2) * t363)) * qJD(1) ^ 2 * t593 + (t401 * t543 + t379) * qJD(3) + (Ifges(7,5) * t384 - Ifges(7,6) * t383) * t556 + (Ifges(7,5) * t197 + Ifges(7,6) * t196) * t557 + t358 * t581 + t362 * t582 + t480 * t583 + t197 * t584 + t497 * t585 + t196 * t586 + t492 * t592 + (t250 * (Ifges(4,5) * t359 + t363 * t406) + t249 * (Ifges(4,6) * t359 + t363 * t403) + t334 * (Ifges(3,5) * t363 - Ifges(3,6) * t359) + (t112 + t156 + t115) * t359) * t433 + (Ifges(4,2) * t362 + t517) * t558 + (Ifges(4,1) * t358 + t516) * t559 + t653 * t299 + t639 * t224 + t640 * t223 - t647 + (t85 - t83) * t229 + (Ifges(4,5) * t358 + Ifges(4,6) * t362) * t546 + t444 + (Ifges(7,4) * t384 - Ifges(7,2) * t383) * t565 + (Ifges(7,4) * t197 + Ifges(7,2) * t196) * t566 - t379 * t438 + (-m(5) * t304 + t295 - t631 * (t349 * t337 + t348 * t453 + t304) + ((-mrSges(7,1) * t479 - mrSges(7,2) * t476) * t348 + ((m(5) + t631) * t356 + t597) * t359 + (t349 * t609 + t608) * t363) * t353) * g(3) + t185 * t411 - t346 * t44 + t624 * t99 + t625 * t98 + (t1 * t123 + t122 * t2 + t15 * t624 + t16 * t625 + t195 * t6 + t41 * t619) * m(7) + (-t157 / 0.2e1 - t510) * t466 + t201 * mrSges(3,1) + t203 * t45 - t191 * t204 - t190 * t205 + t619 * t88 + t620 * t149 + t621 * t147 + t622 * t148 + (-t10 * t229 + t11 * t228 + t20 * t203 + t61 * t622 + t618 * t74 + t62 * t621) * m(6) + t623 * t150 + (-t119 * t346 - t13 * t228 + t14 * t229 + t174 * t613 + t620 * t66 + t623 * t65) * m(5) + t195 * t17 + t419 * t561 + (-m(7) * (-pkin(10) * t496 + t422) + mrSges(7,3) * t496 - m(6) * t422 - m(5) * t474 + t596 * t292 + t600 * t291) * g(2) + (-m(7) * (-pkin(10) * t494 + t421) + mrSges(7,3) * t494 - m(6) * t421 - m(5) * t473 + t596 * t294 + t600 * t293) * g(1) + (t311 * (Ifges(4,3) * t359 + t363 * t401) + t359 * t225) * t471 / 0.2e1 + (qJD(6) * t60 + t8) * t491 / 0.2e1 + (t249 * t403 + t250 * t406) * qJD(3) / 0.2e1 - pkin(2) * t117 + (Ifges(6,4) * t551 + Ifges(5,5) * t552 + Ifges(6,5) * t555 - Ifges(3,2) * t417 + Ifges(5,6) * t554 + t542 * t650 + t599) * t439 + t122 * t25 + t123 * t26 + (-t205 * t465 - t204 * t466 + m(4) * ((-t163 * t362 - t164 * t358) * qJD(3) + t611) - t358 * t134 + t362 * t135) * pkin(9) + t611 * mrSges(4,3) + t613 * t120 + (mrSges(7,1) * t615 - mrSges(7,2) * t614) * t41 + (t1 * t491 + t15 * t614 - t16 * t615 - t2 * t492) * mrSges(7,3) + (-m(4) * t230 + t426 - t616) * t280 + t618 * t121; -t637 + t641 + (-t174 * t533 + t65 * t68 - t66 * t69 + (t13 * t354 + t14 * t352) * pkin(3)) * m(5) + t499 * t584 + (-Ifges(7,2) * t357 + t514) * t587 + (Ifges(7,1) * t361 - t515) * t588 + t361 * t592 + (-t15 * t18 - t16 * t19 + t343 * t6 + (-t51 + qJD(5)) * t41) * m(7) + t157 * t547 - t120 * t533 - t639 * t176 + (-t16 * t520 + t400 * t557 + t402 * t566 + t405 * t564 + t607 + t640) * t393 + t250 * t510 + t249 * t511 - t2 * t520 + (-t10 * t343 + t11 * t345 - t61 * t68 - t74 * t95 + (-qJD(5) + t69) * t62) * m(6) + t60 * t431 + (t147 - t149) * t69 + t86 * t531 + t85 * t532 + (Ifges(4,5) * t249 - Ifges(4,6) * t250) * t542 - t250 * (Ifges(4,1) * t249 - t508) / 0.2e1 + t6 * t408 - t357 * t8 / 0.2e1 + t345 * t84 + (-m(5) * t391 - t412 - t631 * (-t262 * pkin(4) + t391) + t606 * t262 + t602 * (t348 * t355 + t349 * t487)) * g(3) + (-m(5) * t416 - mrSges(4,1) * t220 + mrSges(4,2) * t221 - t631 * (-t214 * pkin(4) + t416) + t606 * t214 + t602 * t215) * g(1) - t230 * (mrSges(4,1) * t250 + mrSges(4,2) * t249) - t163 * t204 + t164 * t205 + t501 * qJD(5) + (Ifges(7,5) * t361 - Ifges(7,6) * t357) * t577 + (-t16 * t463 - t530 + (t464 + t499) * t15 + t601) * mrSges(7,3) - (t140 * t402 + t141 * t405 + t173 * t400) * qJD(6) / 0.2e1 - (-Ifges(4,2) * t250 + t158 + t236) * t249 / 0.2e1 - t51 * t88 + (t17 - t83) * t343 - t19 * t98 - t18 * t99 - t95 * t121 + t646 * t554 + t607 * qJD(6) + (m(7) * t366 + t463 * t98 - t464 * t99 + t612) * (-pkin(10) + t345) + t617 * t68 + (m(5) * t369 + mrSges(4,1) * t373 - mrSges(4,2) * t428 + t606 * t210 + t602 * t211 + t631 * (t210 * pkin(4) + t369)) * g(2); -t357 * t25 + t361 * t26 + t101 - t102 + t521 * t109 + t395 * qJD(6) - (-t149 - t501) * t176 + (t395 + t617) * t393 + (t1 * t361 + t176 * t41 - t2 * t357 + t370 - t173 * (t15 * t361 + t16 * t357)) * m(7) + (-t176 * t62 - t393 * t61 + t20 + t370) * m(6) + (t176 * t66 + t393 * t65 + t119 + t370) * m(5); t501 * t311 + t396 * qJD(6) + (t121 + t396) * t393 + t84 + (t311 * t41 - t393 * t397 + t366 - t601) * m(7) + (-t311 * t62 + t393 * t74 + t11 - t601) * m(6) + t612; -t41 * (mrSges(7,1) * t141 + mrSges(7,2) * t140) + (Ifges(7,1) * t140 - t512) * t564 + t59 * t563 + (Ifges(7,5) * t140 - Ifges(7,6) * t141) * t557 - t15 * t98 + t16 * t99 - g(1) * (mrSges(7,1) * t168 - mrSges(7,2) * t169) - g(2) * ((t210 * t361 - t291 * t357) * mrSges(7,1) + (-t210 * t357 - t291 * t361) * mrSges(7,2)) - g(3) * ((t262 * t361 + t319) * mrSges(7,1) + (-t262 * t357 + t449) * mrSges(7,2)) + (t140 * t15 + t141 * t16) * mrSges(7,3) + t7 + (-Ifges(7,2) * t141 + t139 + t60) * t566 + t605;];
tau  = t9;
