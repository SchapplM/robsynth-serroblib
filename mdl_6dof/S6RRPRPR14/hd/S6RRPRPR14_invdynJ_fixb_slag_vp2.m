% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR14_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR14_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR14_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR14_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:32:46
% EndTime: 2019-03-09 11:33:34
% DurationCPUTime: 34.46s
% Computational Cost: add. (10472->949), mult. (24777->1173), div. (0->0), fcn. (18280->10), ass. (0->435)
t289 = sin(pkin(6));
t293 = sin(qJ(2));
t455 = t289 * t293;
t400 = qJD(2) * t455;
t297 = cos(qJ(2));
t453 = t289 * t297;
t216 = qJD(1) * t400 - qJDD(1) * t453;
t290 = cos(pkin(6));
t421 = qJDD(1) * t290;
t270 = qJDD(2) + t421;
t292 = sin(qJ(4));
t296 = cos(qJ(4));
t436 = qJD(1) * t290;
t272 = qJD(2) + t436;
t437 = qJD(1) * t289;
t401 = t297 * t437;
t185 = t272 * t292 + t296 * t401;
t431 = qJD(4) * t185;
t104 = -t292 * t216 - t296 * t270 + t431;
t516 = -t104 / 0.2e1;
t596 = Ifges(5,1) * t516;
t240 = t292 * t401;
t428 = qJD(4) * t296;
t105 = -qJD(4) * t240 - t296 * t216 + t270 * t292 + t272 * t428;
t514 = -t105 / 0.2e1;
t595 = Ifges(5,4) * t514;
t419 = pkin(1) * t436;
t267 = t297 * t419;
t594 = qJD(3) - t267;
t433 = qJD(2) * t297;
t217 = (qJD(1) * t433 + qJDD(1) * t293) * t289;
t201 = qJDD(4) + t217;
t501 = t201 / 0.2e1;
t429 = qJD(4) * t292;
t593 = pkin(4) * t428 + qJ(5) * t429 + t594;
t102 = qJDD(6) - t104;
t402 = t293 * t437;
t247 = qJD(4) + t402;
t291 = sin(qJ(6));
t295 = cos(qJ(6));
t118 = t185 * t295 - t247 * t291;
t41 = qJD(6) * t118 + t105 * t291 + t201 * t295;
t119 = t185 * t291 + t247 * t295;
t42 = -qJD(6) * t119 + t105 * t295 - t201 * t291;
t7 = Ifges(7,5) * t41 + Ifges(7,6) * t42 + Ifges(7,3) * t102;
t592 = t596 + t595 + Ifges(5,5) * t501 + t7 / 0.2e1;
t591 = -mrSges(3,1) + mrSges(4,2);
t519 = pkin(3) + pkin(8);
t372 = -qJ(5) * t292 - t519;
t518 = pkin(4) + pkin(10);
t590 = (pkin(10) * qJD(4) - qJD(5)) * t296 - (-t296 * t518 + t372) * t402 + t593;
t520 = pkin(2) + pkin(9);
t486 = pkin(5) + t520;
t383 = qJD(4) * t486;
t450 = t292 * t293;
t262 = pkin(2) * t402;
t338 = pkin(9) * t293 - qJ(3) * t297;
t170 = t338 * t437 + t262;
t265 = pkin(8) * t401;
t213 = t293 * t419 + t265;
t172 = pkin(3) * t401 + t213;
t88 = -t292 * t170 + t172 * t296;
t589 = -t292 * t383 - (pkin(5) * t450 - t297 * t518) * t437 + t88;
t387 = -qJ(3) * t293 - pkin(1);
t136 = (-t297 * t520 + t387) * t437;
t308 = -t272 * t520 + t402 * t519 + t594;
t58 = t136 * t292 - t296 * t308;
t53 = -pkin(4) * t247 + qJD(5) + t58;
t186 = t272 * t296 - t240;
t328 = pkin(5) * t186 + t58;
t577 = qJD(5) + t328;
t38 = -t247 * t518 + t577;
t254 = t272 * qJ(3);
t127 = t254 + t172;
t324 = -qJ(5) * t186 + t127;
t52 = t185 * t518 + t324;
t16 = -t291 * t52 + t295 * t38;
t17 = t291 * t38 + t295 * t52;
t539 = t16 * mrSges(7,1) - t17 * mrSges(7,2);
t182 = Ifges(5,4) * t185;
t184 = qJD(6) + t186;
t558 = t186 * Ifges(5,1) + t247 * Ifges(5,5) + t119 * Ifges(7,5) + t118 * Ifges(7,6) + t184 * Ifges(7,3) - t182;
t181 = Ifges(6,6) * t185;
t94 = t247 * Ifges(6,4) - t186 * Ifges(6,2) + t181;
t588 = -t539 - t558 / 0.2e1 - t53 * mrSges(6,1) - t58 * mrSges(5,3) + t94 / 0.2e1;
t513 = t105 / 0.2e1;
t515 = t104 / 0.2e1;
t562 = Ifges(6,5) - Ifges(5,6);
t568 = -t201 / 0.2e1;
t582 = Ifges(6,6) * t515;
t587 = Ifges(5,4) * t515 + Ifges(5,6) * t568 + 0.2e1 * Ifges(6,3) * t513 + t582 + (t513 - t514) * Ifges(5,2) + (t562 + Ifges(6,5)) * t501;
t21 = mrSges(7,1) * t102 - mrSges(7,3) * t41;
t22 = -mrSges(7,2) * t102 + mrSges(7,3) * t42;
t85 = -mrSges(7,2) * t184 + mrSges(7,3) * t118;
t86 = mrSges(7,1) * t184 - mrSges(7,3) * t119;
t550 = -t291 * t86 + t295 * t85;
t586 = -qJD(6) * t550 - t295 * t21 - t291 * t22;
t517 = t102 / 0.2e1;
t525 = t42 / 0.2e1;
t526 = t41 / 0.2e1;
t495 = pkin(1) * t290;
t418 = qJD(2) * t495;
t375 = qJD(1) * t418;
t408 = pkin(1) * t421;
t120 = -pkin(8) * t216 + t293 * t408 + t297 * t375;
t103 = -t270 * qJ(3) - t272 * qJD(3) - t120;
t74 = -pkin(3) * t216 - t103;
t301 = qJ(5) * t104 - qJD(5) * t186 + t74;
t11 = t105 * t518 + t301;
t432 = qJD(3) * t293;
t309 = -qJ(3) * t217 + (-pkin(1) * qJDD(1) - qJD(1) * t432) * t289;
t580 = qJD(4) * t308 + t216 * t520 + t309;
t273 = pkin(8) * t455;
t121 = -qJD(2) * t265 - qJDD(1) * t273 - t293 * t375 + t297 * t408;
t313 = qJDD(3) - t121;
t73 = pkin(3) * t217 - t270 * t520 + t313;
t15 = -t136 * t428 - t292 * t580 + t296 * t73;
t323 = qJDD(5) - t15;
t5 = -pkin(5) * t104 - t201 * t518 + t323;
t1 = qJD(6) * t16 + t11 * t295 + t291 * t5;
t2 = -qJD(6) * t17 - t11 * t291 + t295 * t5;
t541 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t565 = Ifges(6,4) - Ifges(5,5);
t585 = -t501 * t565 + t541 + t596 + Ifges(6,4) * t568 + Ifges(7,5) * t526 + Ifges(6,6) * t514 + Ifges(7,6) * t525 + Ifges(7,3) * t517 + (-t515 + t516) * Ifges(6,2);
t337 = t16 * t291 - t17 * t295;
t491 = t1 * t291;
t307 = -t337 * qJD(6) + t2 * t295 + t491;
t490 = t185 * pkin(5);
t59 = t296 * t136 + t292 * t308;
t55 = -t247 * qJ(5) - t59;
t43 = -t55 - t490;
t579 = t247 * t43 - t307;
t355 = mrSges(7,1) * t291 + mrSges(7,2) * t295;
t318 = m(7) * qJ(5) + t355;
t358 = mrSges(5,1) * t292 + mrSges(5,2) * t296;
t578 = t318 * t296 - t358;
t566 = mrSges(6,2) - mrSges(5,1);
t569 = m(7) + m(6);
t570 = m(7) * pkin(10);
t576 = pkin(4) * t569 - t566 + t570;
t60 = pkin(4) * t185 + t324;
t575 = t127 * mrSges(5,1) - t60 * mrSges(6,2);
t574 = -t127 * mrSges(5,2) + t60 * mrSges(6,3);
t567 = -mrSges(6,1) - mrSges(5,3);
t488 = mrSges(4,3) - mrSges(3,2);
t564 = Ifges(3,5) - Ifges(4,4);
t563 = Ifges(4,5) - Ifges(3,6);
t287 = t292 * pkin(4);
t382 = -qJ(5) * t296 + qJ(3);
t237 = t287 + t382;
t492 = pkin(10) * t292;
t232 = t237 + t492;
t239 = t486 * t296;
t134 = -t232 * t291 + t239 * t295;
t561 = qJD(6) * t134 + t291 * t589 + t295 * t590;
t135 = t232 * t295 + t239 * t291;
t560 = -qJD(6) * t135 - t291 * t590 + t295 * t589;
t68 = -mrSges(5,2) * t201 - mrSges(5,3) * t105;
t69 = mrSges(6,1) * t105 - mrSges(6,3) * t201;
t559 = -t69 + t68;
t485 = mrSges(6,1) * t185;
t128 = -mrSges(6,3) * t247 + t485;
t61 = -mrSges(7,1) * t118 + mrSges(7,2) * t119;
t458 = t128 - t61;
t446 = t293 * t296;
t89 = t296 * t170 + t292 * t172;
t556 = -t296 * t383 - (pkin(5) * t446 + qJ(5) * t297) * t437 - t89;
t484 = mrSges(6,1) * t186;
t129 = mrSges(6,2) * t247 + t484;
t482 = mrSges(5,3) * t186;
t131 = mrSges(5,1) * t247 - t482;
t441 = t129 - t131;
t179 = (-t291 * t297 - t295 * t446) * t437;
t425 = qJD(6) * t292;
t555 = t291 * t425 - t295 * t428 + t179;
t180 = (-t291 * t446 + t295 * t297) * t437;
t424 = qJD(6) * t295;
t554 = -t291 * t428 - t292 * t424 + t180;
t553 = -qJD(5) * t296 - (-pkin(4) * t296 + t372) * t402 + t593;
t115 = mrSges(5,1) * t185 + mrSges(5,2) * t186;
t378 = mrSges(4,1) * t401;
t208 = -mrSges(4,3) * t272 - t378;
t552 = -t208 + t115;
t377 = mrSges(3,3) * t402;
t379 = mrSges(4,1) * t402;
t551 = t272 * t591 + t377 + t379;
t14 = -t136 * t429 + t292 * t73 + t296 * t580;
t548 = t14 * t292 + t15 * t296;
t10 = -qJ(5) * t201 - qJD(5) * t247 - t14;
t12 = -pkin(4) * t201 + t323;
t547 = -t10 * t292 - t12 * t296;
t546 = -m(6) * qJ(5) - mrSges(6,3);
t545 = -t567 - t591;
t356 = t295 * mrSges(7,1) - t291 * mrSges(7,2);
t477 = Ifges(7,4) * t119;
t50 = Ifges(7,2) * t118 + Ifges(7,6) * t184 + t477;
t460 = t295 * t50;
t544 = t43 * t356 - t460 / 0.2e1;
t422 = m(5) + t569;
t410 = m(4) + t422;
t543 = qJ(3) * t410 + t488;
t264 = pkin(2) * t400;
t133 = t264 + (qJD(2) * t338 - t432) * t289;
t494 = pkin(1) * t297;
t403 = -pkin(2) - t494;
t139 = pkin(3) * t455 + t273 + (-pkin(9) + t403) * t290;
t439 = pkin(2) * t453 + qJ(3) * t455;
t496 = pkin(1) * t289;
t193 = -t439 - t496;
t276 = pkin(9) * t453;
t164 = t193 - t276;
t282 = t293 * t495;
t173 = (t453 * t519 + t282) * qJD(2);
t322 = -t296 * t133 - t139 * t428 + t164 * t429 - t292 * t173;
t29 = -t289 * (qJ(5) * t433 + qJD(5) * t293) + t322;
t540 = -m(7) * (-pkin(5) - pkin(9)) + t356 + t545;
t70 = -t104 * mrSges(6,1) + t201 * mrSges(6,2);
t538 = t586 - t70;
t537 = t127 * (mrSges(5,1) * t296 - mrSges(5,2) * t292) - t60 * (mrSges(6,2) * t296 - mrSges(6,3) * t292);
t536 = -mrSges(5,2) + t318 - t546;
t535 = pkin(9) * t422 + t545;
t294 = sin(qJ(1));
t444 = t294 * t297;
t298 = cos(qJ(1));
t445 = t293 * t298;
t227 = t290 * t444 + t445;
t454 = t289 * t294;
t156 = -t227 * t296 + t292 * t454;
t442 = t297 * t298;
t447 = t293 * t294;
t225 = -t290 * t442 + t447;
t452 = t289 * t298;
t158 = t225 * t296 + t292 * t452;
t223 = t290 * t292 + t296 * t453;
t534 = g(1) * t156 - g(2) * t158 + g(3) * t223;
t462 = t292 * mrSges(6,2);
t351 = -t296 * mrSges(6,3) - t462;
t533 = -m(6) * t382 - m(7) * t492 - t292 * mrSges(7,3) - t351 - t488 + t578 + (-m(4) - m(7) - m(5)) * qJ(3);
t507 = -t184 / 0.2e1;
t510 = -t119 / 0.2e1;
t512 = -t118 / 0.2e1;
t532 = -Ifges(7,5) * t510 - Ifges(7,6) * t512 - Ifges(7,3) * t507 + t539;
t531 = Ifges(7,1) * t526 + Ifges(7,4) * t525 + Ifges(7,5) * t517;
t530 = m(7) * pkin(5);
t117 = Ifges(7,4) * t118;
t51 = Ifges(7,1) * t119 + Ifges(7,5) * t184 + t117;
t524 = -t51 / 0.2e1;
t466 = t186 * Ifges(5,4);
t95 = -t185 * Ifges(5,2) + t247 * Ifges(5,6) + t466;
t522 = -t95 / 0.2e1;
t511 = t118 / 0.2e1;
t509 = t119 / 0.2e1;
t506 = t184 / 0.2e1;
t505 = -t185 / 0.2e1;
t504 = t185 / 0.2e1;
t503 = -t186 / 0.2e1;
t502 = t186 / 0.2e1;
t499 = -t247 / 0.2e1;
t498 = t247 / 0.2e1;
t493 = pkin(10) * t223;
t487 = Ifges(3,4) + Ifges(4,6);
t483 = mrSges(5,3) * t185;
t481 = mrSges(7,3) * t295;
t480 = Ifges(3,4) * t293;
t479 = Ifges(5,4) * t292;
t478 = Ifges(5,4) * t296;
t476 = Ifges(7,4) * t291;
t475 = Ifges(7,4) * t295;
t474 = Ifges(4,6) * t293;
t473 = Ifges(4,6) * t297;
t472 = Ifges(6,6) * t292;
t471 = Ifges(6,6) * t296;
t465 = t186 * Ifges(6,6);
t457 = qJ(5) * t185;
t456 = t186 * t291;
t451 = t291 * t292;
t449 = t292 * t295;
t84 = t292 * t139 + t296 * t164;
t231 = pkin(8) * t453 + t282;
t438 = t298 * pkin(1) + pkin(8) * t454;
t434 = qJD(1) ^ 2 * t289 ^ 2;
t427 = qJD(4) * t520;
t426 = qJD(6) * t291;
t420 = pkin(4) * t455;
t417 = Ifges(6,1) / 0.2e1 + Ifges(5,3) / 0.2e1;
t416 = Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1;
t415 = Ifges(3,5) / 0.2e1 - Ifges(4,4) / 0.2e1;
t414 = Ifges(6,5) / 0.2e1 - Ifges(5,6) / 0.2e1;
t413 = -Ifges(3,6) / 0.2e1 + Ifges(4,5) / 0.2e1;
t412 = Ifges(4,6) / 0.2e1 + Ifges(3,4) / 0.2e1;
t409 = mrSges(7,3) + t570;
t407 = t292 * t453;
t228 = -t290 * t447 + t442;
t405 = t228 * pkin(2) + t438;
t404 = t276 + t439;
t192 = -t290 * qJ(3) - t231;
t399 = t289 * t433;
t398 = t292 * t427;
t397 = t296 * t427;
t392 = t437 / 0.2e1;
t388 = -t426 / 0.2e1;
t386 = -pkin(1) * t294 + pkin(8) * t452;
t219 = t225 * pkin(2);
t385 = -pkin(9) * t225 - t219;
t221 = t227 * pkin(2);
t384 = -pkin(9) * t227 - t221;
t151 = t217 * mrSges(4,1) + t270 * mrSges(4,2);
t218 = t223 * pkin(4);
t224 = t290 * t296 - t407;
t381 = qJ(5) * t224 - t218;
t83 = t139 * t296 - t292 * t164;
t380 = t519 * t455;
t376 = mrSges(3,3) * t401;
t374 = pkin(3) * t454 + t405;
t163 = pkin(3) * t453 - t192;
t371 = t292 * t402;
t366 = t293 * t392;
t226 = t290 * t445 + t444;
t364 = -t226 * pkin(2) + t386;
t361 = -t409 + t566;
t360 = mrSges(5,1) * t223 + mrSges(5,2) * t224;
t154 = t223 * t295 - t291 * t455;
t155 = t223 * t291 + t295 * t455;
t357 = mrSges(7,1) * t154 - mrSges(7,2) * t155;
t354 = mrSges(4,2) * t297 - mrSges(4,3) * t293;
t353 = -t223 * mrSges(6,2) - t224 * mrSges(6,3);
t350 = Ifges(5,1) * t292 + t478;
t349 = Ifges(7,1) * t295 - t476;
t348 = Ifges(7,1) * t291 + t475;
t347 = Ifges(5,2) * t296 + t479;
t346 = Ifges(6,4) * t292 + Ifges(6,5) * t296;
t345 = -Ifges(7,2) * t291 + t475;
t344 = Ifges(7,2) * t295 + t476;
t343 = Ifges(5,5) * t292 + Ifges(5,6) * t296;
t342 = Ifges(7,5) * t295 - Ifges(7,6) * t291;
t341 = Ifges(7,5) * t291 + Ifges(7,6) * t295;
t340 = Ifges(6,2) * t292 + t471;
t339 = Ifges(6,3) * t296 + t472;
t54 = pkin(5) * t224 - t455 * t518 - t83;
t87 = t163 - t381;
t64 = t87 + t493;
t23 = -t291 * t64 + t295 * t54;
t24 = t291 * t54 + t295 * t64;
t334 = t292 * t53 - t296 * t55;
t333 = t292 * t58 + t296 * t59;
t330 = -qJ(5) * t569 + mrSges(5,2) - mrSges(6,3);
t329 = pkin(3) * t452 + t364;
t212 = pkin(8) * t402 - t267;
t268 = t297 * t418;
t214 = -pkin(8) * t400 + t268;
t78 = -qJ(5) * t455 - t84;
t37 = -t292 * t133 - t139 * t429 - t164 * t428 + t173 * t296;
t327 = -t225 * t292 + t296 * t452;
t321 = t356 + t530;
t285 = t290 * qJD(3);
t138 = -qJD(2) * t380 + t268 + t285;
t312 = t15 * mrSges(5,1) - t14 * mrSges(5,2) + t12 * mrSges(6,2) - t10 * mrSges(6,3);
t112 = -pkin(2) * t270 + t313;
t311 = t121 * mrSges(3,1) - t120 * mrSges(3,2) + t112 * mrSges(4,2) - t103 * mrSges(4,3);
t306 = qJD(4) * t334 + t547;
t305 = qJD(4) * t333 + t548;
t152 = qJD(4) * t223 - t292 * t400;
t304 = qJ(5) * t152 - qJD(5) * t224 + t138;
t261 = Ifges(3,4) * t401;
t253 = Ifges(4,1) * t270;
t252 = Ifges(3,3) * t270;
t238 = t486 * t292;
t230 = t290 * t494 - t273;
t229 = (-mrSges(3,1) * t297 + mrSges(3,2) * t293) * t289;
t215 = t231 * qJD(2);
t211 = -qJ(3) * t401 + t262;
t210 = t354 * t437;
t207 = -mrSges(3,2) * t272 + t376;
t200 = Ifges(4,4) * t217;
t199 = Ifges(3,5) * t217;
t198 = Ifges(4,5) * t216;
t197 = Ifges(3,6) * t216;
t196 = t290 * t403 + t273;
t195 = t228 * t287;
t194 = t226 * t287;
t191 = Ifges(6,1) * t201;
t190 = Ifges(5,3) * t201;
t183 = -t214 - t285;
t177 = t272 * t295 - t291 * t371;
t176 = -t272 * t291 - t295 * t371;
t175 = (-pkin(2) * t297 + t387) * t437;
t174 = t264 + (-qJ(3) * t433 - t432) * t289;
t171 = -qJD(1) * t380 + t267;
t169 = -t254 - t213;
t168 = t272 * Ifges(4,4) + (-Ifges(4,2) * t293 - t473) * t437;
t167 = t272 * Ifges(4,5) + (-t297 * Ifges(4,3) - t474) * t437;
t166 = Ifges(3,1) * t402 + t272 * Ifges(3,5) + t261;
t165 = t272 * Ifges(3,6) + (t297 * Ifges(3,2) + t480) * t437;
t162 = -pkin(2) * t272 + qJD(3) + t212;
t157 = t227 * t292 + t296 * t454;
t153 = -qJD(4) * t407 + t290 * t428 - t296 * t400;
t150 = mrSges(4,1) * t216 - mrSges(4,3) * t270;
t130 = -mrSges(5,2) * t247 - t483;
t116 = -mrSges(6,2) * t185 - mrSges(6,3) * t186;
t114 = pkin(4) * t186 + t457;
t108 = t156 * t291 + t228 * t295;
t107 = t156 * t295 - t228 * t291;
t106 = pkin(2) * t216 + t309;
t101 = Ifges(6,4) * t104;
t100 = Ifges(5,5) * t104;
t99 = Ifges(6,5) * t105;
t98 = Ifges(5,6) * t105;
t96 = t247 * Ifges(6,1) - t186 * Ifges(6,4) + t185 * Ifges(6,5);
t93 = t186 * Ifges(5,5) - t185 * Ifges(5,6) + t247 * Ifges(5,3);
t92 = t247 * Ifges(6,5) + t185 * Ifges(6,3) - t465;
t82 = -pkin(4) * t401 - t88;
t81 = -qJ(5) * t401 - t89;
t80 = t186 * t518 + t457;
t79 = -t83 - t420;
t76 = qJD(6) * t154 + t153 * t291 + t295 * t399;
t75 = -qJD(6) * t155 + t153 * t295 - t291 * t399;
t67 = mrSges(5,1) * t201 + mrSges(5,3) * t104;
t56 = -pkin(5) * t223 - t78;
t48 = t59 - t490;
t46 = pkin(4) * t153 + t304;
t45 = mrSges(5,1) * t105 - mrSges(5,2) * t104;
t44 = -mrSges(6,2) * t105 + mrSges(6,3) * t104;
t31 = t153 * t518 + t304;
t30 = -pkin(4) * t399 - t37;
t26 = -pkin(5) * t153 - t29;
t25 = -pkin(5) * t152 - t399 * t518 - t37;
t20 = t291 * t48 + t295 * t80;
t19 = -t291 * t80 + t295 * t48;
t18 = pkin(4) * t105 + t301;
t13 = -mrSges(7,1) * t42 + mrSges(7,2) * t41;
t8 = t41 * Ifges(7,4) + t42 * Ifges(7,2) + t102 * Ifges(7,6);
t6 = -pkin(5) * t105 - t10;
t4 = -qJD(6) * t24 + t25 * t295 - t291 * t31;
t3 = qJD(6) * t23 + t25 * t291 + t295 * t31;
t9 = [(Ifges(7,1) * t76 + Ifges(7,4) * t75) * t509 + (Ifges(7,1) * t155 + Ifges(7,4) * t154) * t526 + (-Ifges(5,1) * t502 - Ifges(5,4) * t505 - Ifges(7,5) * t509 + Ifges(6,2) * t503 + Ifges(6,6) * t504 - Ifges(7,6) * t511 - Ifges(7,3) * t506 + t498 * t565 + t574 + t588) * t152 + (Ifges(7,5) * t76 + Ifges(7,6) * t75) * t506 + (Ifges(7,5) * t155 + Ifges(7,6) * t154) * t517 + (mrSges(6,1) * t10 - mrSges(5,3) * t14 - Ifges(5,4) * t516 + t582 + t587) * t223 + m(5) * (t127 * t138 + t14 * t84 + t15 * t83 + t163 * t74 - t322 * t59 - t37 * t58) - t322 * t130 + (-Ifges(5,2) * t505 + t522 - Ifges(5,4) * t502 + Ifges(6,6) * t503 + Ifges(6,3) * t504 - t59 * mrSges(5,3) + t92 / 0.2e1 + t55 * mrSges(6,1) + t562 * t498 + t575) * t153 + t18 * t353 - t6 * t357 + t74 * t360 + (t199 / 0.2e1 - t197 / 0.2e1 + t252 / 0.2e1 + t253 / 0.2e1 - t200 / 0.2e1 + t198 / 0.2e1 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t270 + t415 * t217 + t413 * t216 + t311) * t290 + t155 * t531 + (t1 * t154 - t155 * t2 - t16 * t76 + t17 * t75) * mrSges(7,3) + m(3) * (t120 * t231 + t121 * t230 + t212 * t215 + t213 * t214) + Ifges(2,3) * qJDD(1) + t230 * (mrSges(3,1) * t270 - mrSges(3,3) * t217) + t231 * (-mrSges(3,2) * t270 - mrSges(3,3) * t216) + t193 * (-mrSges(4,2) * t216 - mrSges(4,3) * t217) + t214 * t207 + t183 * t208 + t174 * t210 + t192 * t150 + t196 * t151 + t163 * t45 + t154 * t8 / 0.2e1 + t138 * t115 + t29 * t128 + t30 * t129 + t37 * t131 + t46 * t116 + t83 * t67 + t84 * t68 + t3 * t85 + t4 * t86 + t87 * t44 + t43 * (-mrSges(7,1) * t75 + mrSges(7,2) * t76) + t76 * t51 / 0.2e1 + t78 * t69 + t79 * t70 + t75 * t50 / 0.2e1 + t26 * t61 + t56 * t13 + t23 * t21 + t24 * t22 + (t12 * mrSges(6,1) - t15 * mrSges(5,3) - Ifges(6,6) * t513 + t585 + t592 + t595) * t224 + t551 * t215 + (Ifges(7,4) * t76 + Ifges(7,2) * t75) * t511 + (Ifges(7,4) * t155 + Ifges(7,2) * t154) * t525 + ((-mrSges(3,1) * t216 - mrSges(3,2) * t217 + (m(3) * t496 - t229) * qJDD(1)) * pkin(1) + (-t103 * mrSges(4,1) + t106 * mrSges(4,2) + t120 * mrSges(3,3) - t563 * t270 + t487 * t217 + (-Ifges(3,2) - Ifges(4,3)) * t216) * t297 + (t190 / 0.2e1 + t191 / 0.2e1 + t101 / 0.2e1 + t99 / 0.2e1 - t100 / 0.2e1 - t98 / 0.2e1 - t121 * mrSges(3,3) + t112 * mrSges(4,1) - t106 * mrSges(4,3) + t564 * t270 + (Ifges(3,1) + Ifges(4,2)) * t217 - t487 * t216 + t417 * t201 + t414 * t105 + t416 * t104 + t312) * t293 + ((-t175 * mrSges(4,2) + t169 * mrSges(4,1) - t213 * mrSges(3,3) - t165 / 0.2e1 + t167 / 0.2e1 + t413 * t272 + (-pkin(1) * mrSges(3,1) - t293 * t412) * t437) * t293 + (-t175 * mrSges(4,3) + t212 * mrSges(3,3) + t162 * mrSges(4,1) - t58 * mrSges(5,1) - t55 * mrSges(6,3) - t59 * mrSges(5,2) + t53 * mrSges(6,2) + t166 / 0.2e1 - t168 / 0.2e1 + t93 / 0.2e1 + t96 / 0.2e1 + t415 * t272 + t417 * t247 - t416 * t186 + t414 * t185 + ((-pkin(1) * mrSges(3,2) + t297 * t412) * t289 + (-Ifges(4,3) / 0.2e1 + Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1 + Ifges(4,2) / 0.2e1) * t455) * qJD(1)) * t297) * qJD(2) + (mrSges(3,3) + mrSges(4,1)) * (-g(1) * t298 - g(2) * t294)) * t289 + (-m(4) * t364 - m(3) * t386 + mrSges(2,1) * t294 + mrSges(2,2) * t298 - m(5) * t329 + t543 * t225 + t361 * t327 + (t330 - t355) * t158 + (t321 + t535) * t226 + t569 * (-pkin(4) * t327 - t329)) * g(1) + (-m(4) * t405 - m(3) * t438 - mrSges(2,1) * t298 + mrSges(2,2) * t294 - m(5) * t374 - t108 * mrSges(7,1) - t107 * mrSges(7,2) - t543 * t227 + t330 * t156 + t361 * t157 + (-t530 - t535) * t228 - t569 * (t157 * pkin(4) + t374)) * g(2) + m(7) * (t1 * t24 + t16 * t4 + t17 * t3 + t2 * t23 + t26 * t43 + t56 * t6) + m(6) * (t10 * t78 + t12 * t79 + t18 * t87 + t29 * t55 + t30 * t53 + t46 * t60) + m(4) * (t103 * t192 + t106 * t193 + t112 * t196 + t162 * t215 + t169 * t183 + t174 * t175); (t74 * qJ(3) + t58 * t88 - t59 * t89 + (-t171 + qJD(3)) * t127) * m(5) + (t18 * t237 - t53 * t82 - t55 * t81 + t553 * t60) * m(6) + (-m(4) * t169 + t207 - t208 - t376) * t212 - t532 * t371 + (-t293 * (Ifges(3,1) * t297 - t480) / 0.2e1 + pkin(1) * (mrSges(3,1) * t293 + mrSges(3,2) * t297)) * t434 + (t366 * t92 + t585) * t296 + (-t55 * (-mrSges(6,1) * t446 - mrSges(6,3) * t297) - t59 * (-mrSges(5,2) * t297 + mrSges(5,3) * t446) + t58 * (mrSges(5,1) * t297 - mrSges(5,3) * t450) - t53 * (mrSges(6,1) * t450 + mrSges(6,2) * t297) - t175 * (-mrSges(4,2) * t293 - mrSges(4,3) * t297)) * t437 + (-m(6) * (t195 + t384) - m(5) * t384 - m(7) * (t195 - t221) + m(4) * t221 + t533 * t228 + t540 * t227) * g(1) + (-m(6) * (t194 + t385) - m(5) * t385 - m(7) * (t194 - t219) + m(4) * t219 + t533 * t226 + t540 * t225) * g(2) + t537 * t402 + (Ifges(7,1) * t180 + Ifges(7,4) * t179) * t510 + (t428 * t55 - t547) * mrSges(6,1) + (-t428 * t59 - t548) * mrSges(5,3) + (-pkin(2) * t112 - qJ(3) * t103 - qJD(3) * t169 - t175 * t211) * m(4) + t588 * t429 - t471 * t513 + t478 * t514 - t162 * t378 + t472 * t515 - t479 * t516 + t252 + t253 + (t341 * t517 + t344 * t525 + t348 * t526 - t356 * t6 + t366 * t94 + t388 * t50 + t587) * t292 + (-(t67 - t70) * t296 - t559 * t292 - m(5) * t305 - m(6) * t306) * t520 + (-t397 - t89) * t130 + (t397 - t81) * t128 + (-t398 - t82) * t129 + (t398 - t88) * t131 + (-m(4) * t439 - m(5) * t404 + t229 - t569 * (t292 * t420 + t404) + (t354 + (-t321 + t567) * t297 + (-t292 * t409 - t296 * t546 + t462 + t578) * t293) * t289) * g(3) + t198 + t199 - t200 - t197 + t18 * t351 + t74 * t358 + (Ifges(7,4) * t180 + Ifges(7,2) * t179) * t512 + t428 * t522 + t165 * t366 + (t186 * (Ifges(6,4) * t297 - t293 * t340) + t185 * (Ifges(5,6) * t297 + t293 * t347) + t297 * t168) * t392 + t180 * t524 + t451 * t531 - t169 * t379 + ((-Ifges(7,3) * t292 + t296 * t341) * t506 + (-Ifges(7,5) * t292 + t296 * t348) * t509 + (-Ifges(7,6) * t292 + t296 * t344) * t511 + t346 * t498 + t343 * t499 + t537 + (-t350 / 0.2e1 - t340 / 0.2e1) * t186) * qJD(4) + t296 * t592 + (-t150 + t45) * qJ(3) + (t342 * t506 + t345 * t511 + t349 * t509) * t425 + (t291 * t51 + t460 + t92) * t428 / 0.2e1 + (t297 * (Ifges(4,3) * t293 - t473) + t293 * (-Ifges(4,2) * t297 + t474)) * t434 / 0.2e1 + (qJD(6) * t51 + t8) * t449 / 0.2e1 - ((t292 * t558 + t296 * t95 + t167) * t293 + t186 * (Ifges(5,5) * t297 + t293 * t350) + t185 * (Ifges(6,5) * t297 - t293 * t339) + (-Ifges(3,2) * t402 + t166 + t261 + t93 + t96) * t297 + ((Ifges(6,1) + Ifges(5,3)) * t297 + (t343 - t346) * t293) * t247 + (t293 * t563 + t297 * t564) * t272) * t437 / 0.2e1 + (t347 / 0.2e1 + t339 / 0.2e1) * t431 + t237 * t44 - t238 * t13 - t211 * t210 - t179 * t50 / 0.2e1 - t171 * t115 - pkin(2) * t151 + t134 * t21 + t135 * t22 + (-m(4) * t162 + t377 - t551) * t213 + t552 * qJD(3) + t553 * t116 + (mrSges(7,1) * t555 - mrSges(7,2) * t554) * t43 + (t1 * t449 + t16 * t554 - t17 * t555 - t2 * t451) * mrSges(7,3) + t556 * t61 + t560 * t86 + t561 * t85 + (t1 * t135 + t134 * t2 + t560 * t16 + t561 * t17 - t238 * t6 + t556 * t43) * m(7) + (Ifges(7,5) * t180 + Ifges(7,6) * t179) * t507 + t311; -t176 * t86 - t177 * t85 + (-t116 - t552) * t272 + t210 * t402 + (t13 + t441 * t402 + (t291 * t85 + t295 * t86 + t441) * qJD(4) + t559) * t292 + (t67 + t247 * (t130 - t458) + t538) * t296 + t151 + (-g(1) * t227 - g(2) * t225 + g(3) * t453) * t410 + (-t16 * t176 - t17 * t177 + (t6 + (t16 * t295 + t17 * t291) * qJD(4)) * t292 + t579 * t296) * m(7) + (-t272 * t60 + t334 * t402 + t306) * m(6) + (-t127 * t272 + t333 * t402 + t305) * m(5) + (t169 * t272 + t175 * t402 + t112) * m(4); (t92 - t466) * t503 + (t94 + t181) * t505 + (-t17 * t424 - t491 + (t426 + t456) * t16 + t534) * mrSges(7,3) + t544 * qJD(6) + (qJ(5) * t6 - t16 * t19 - t17 * t20 + t43 * t577) * m(7) + t312 + (-m(7) * t307 + t586) * t518 + (t95 + t465) * t502 + (t360 - m(6) * t381 + t353 - m(7) * (-t218 - t493) - t318 * t224) * g(3) + t328 * t61 + (-Ifges(5,1) * t503 + Ifges(6,2) * t502 + t499 * t565 + t532 - t574) * t185 + (Ifges(6,3) * t505 - t17 * t481 + t341 * t507 + t344 * t512 + t348 * t510 + t499 * t562 + t544 - t575) * t186 - t55 * t484 + t53 * t485 + t190 + t191 + t6 * t355 + t51 * t388 + t342 * t517 + t456 * t524 + t345 * t525 + t349 * t526 + t295 * t531 - t2 * t481 + (-t69 + t13) * qJ(5) + (-m(6) * t55 - t128 + t130 + t483) * t58 + t101 + t99 - t100 - t98 - (t118 * t344 + t119 * t348 + t184 * t341) * qJD(6) / 0.2e1 + (-pkin(4) * t12 - qJ(5) * t10 - qJD(5) * t55 - t114 * t60) * m(6) + (t156 * t576 - t157 * t536) * g(1) + (-t158 * t576 + t327 * t536) * g(2) - t291 * t8 / 0.2e1 - t114 * t116 - t20 * t85 - t19 * t86 - pkin(4) * t70 + (-m(6) * t53 - t441 + t482) * t59 - t458 * qJD(5) + (-Ifges(5,2) * t186 - t182 + t558) * t504; t458 * t247 + (t116 + t550) * t186 + (-t186 * t337 - t534 - t579) * m(7) + (t186 * t60 + t247 * t55 + t12 - t534) * m(6) - t538; -t43 * (mrSges(7,1) * t119 + mrSges(7,2) * t118) + (Ifges(7,1) * t118 - t477) * t510 + t50 * t509 + (Ifges(7,5) * t118 - Ifges(7,6) * t119) * t507 - t16 * t85 + t17 * t86 - g(1) * (mrSges(7,1) * t107 - mrSges(7,2) * t108) - g(2) * ((-t158 * t295 - t226 * t291) * mrSges(7,1) + (t158 * t291 - t226 * t295) * mrSges(7,2)) - g(3) * t357 + (t118 * t16 + t119 * t17) * mrSges(7,3) + t7 + (-Ifges(7,2) * t119 + t117 + t51) * t512 + t541;];
tau  = t9;
