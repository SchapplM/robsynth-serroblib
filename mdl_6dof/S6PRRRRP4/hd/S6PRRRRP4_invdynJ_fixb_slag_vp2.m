% Calculate vector of inverse dynamics joint torques for
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:13:10
% EndTime: 2019-03-09 00:14:01
% DurationCPUTime: 30.84s
% Computational Cost: add. (10277->759), mult. (23134->1003), div. (0->0), fcn. (17264->14), ass. (0->366)
t307 = qJ(4) + qJ(5);
t304 = sin(t307);
t305 = cos(t307);
t598 = mrSges(6,2) - mrSges(7,3);
t528 = -m(7) * qJ(6) + t598;
t599 = mrSges(6,1) + mrSges(7,1);
t605 = m(7) * pkin(5) + t599;
t608 = -t304 * t528 + t305 * t605;
t567 = -Ifges(6,4) + Ifges(7,5);
t607 = t567 + Ifges(7,5);
t313 = sin(qJ(3));
t317 = cos(qJ(3));
t410 = qJD(2) * qJD(3);
t261 = qJDD(2) * t313 + t317 * t410;
t312 = sin(qJ(4));
t316 = cos(qJ(4));
t417 = qJD(3) * t316;
t421 = qJD(2) * t313;
t339 = t312 * t421 - t417;
t331 = t339 * qJD(4);
t142 = qJDD(3) * t312 + t261 * t316 - t331;
t252 = qJD(3) * t312 + t316 * t421;
t143 = -qJD(4) * t252 + qJDD(3) * t316 - t261 * t312;
t311 = sin(qJ(5));
t315 = cos(qJ(5));
t154 = t311 * t252 + t315 * t339;
t46 = -qJD(5) * t154 + t315 * t142 + t311 * t143;
t513 = t46 / 0.2e1;
t324 = t315 * t252 - t311 * t339;
t47 = qJD(5) * t324 + t311 * t142 - t315 * t143;
t511 = t47 / 0.2e1;
t260 = qJDD(2) * t317 - t313 * t410;
t249 = qJDD(4) - t260;
t243 = qJDD(5) + t249;
t494 = t243 / 0.2e1;
t314 = sin(qJ(2));
t309 = sin(pkin(6));
t424 = qJD(1) * t309;
t394 = t314 * t424;
t263 = qJD(2) * pkin(8) + t394;
t310 = cos(pkin(6));
t423 = qJD(1) * t310;
t196 = -t313 * t263 + t317 * t423;
t364 = pkin(3) * t313 - pkin(9) * t317;
t257 = t364 * qJD(2);
t130 = -t196 * t312 + t316 * t257;
t435 = t316 * t317;
t344 = pkin(4) * t313 - pkin(10) * t435;
t319 = -pkin(10) - pkin(9);
t395 = qJD(4) * t319;
t602 = -qJD(2) * t344 + t316 * t395 - t130;
t131 = t316 * t196 + t312 * t257;
t419 = qJD(2) * t317;
t390 = t312 * t419;
t601 = -pkin(10) * t390 - t312 * t395 + t131;
t181 = -qJD(3) * pkin(3) - t196;
t132 = pkin(4) * t339 + t181;
t293 = qJD(4) - t419;
t282 = qJD(5) + t293;
t197 = t317 * t263 + t313 * t423;
t182 = qJD(3) * pkin(9) + t197;
t267 = -pkin(3) * t317 - pkin(9) * t313 - pkin(2);
t318 = cos(qJ(2));
t393 = t318 * t424;
t199 = qJD(2) * t267 - t393;
t109 = t316 * t182 + t312 * t199;
t72 = -pkin(10) * t339 + t109;
t461 = t311 * t72;
t108 = -t182 * t312 + t316 * t199;
t71 = -pkin(10) * t252 + t108;
t61 = pkin(4) * t293 + t71;
t29 = t315 * t61 - t461;
t551 = qJD(6) - t29;
t26 = -pkin(5) * t282 + t551;
t48 = t154 * pkin(5) - qJ(6) * t324 + t132;
t146 = Ifges(6,4) * t154;
t464 = Ifges(7,5) * t154;
t566 = Ifges(7,4) + Ifges(6,5);
t568 = Ifges(6,1) + Ifges(7,1);
t560 = t282 * t566 + t324 * t568 - t146 + t464;
t600 = mrSges(6,2) * t132 - mrSges(7,3) * t48 + t26 * mrSges(7,2) - t29 * mrSges(6,3) + t560 / 0.2e1;
t565 = -Ifges(6,6) + Ifges(7,6);
t564 = Ifges(6,3) + Ifges(7,2);
t489 = t282 / 0.2e1;
t498 = t324 / 0.2e1;
t501 = t154 / 0.2e1;
t502 = -t154 / 0.2e1;
t597 = Ifges(6,4) * t502 + Ifges(7,5) * t501 + t489 * t566 + t498 * t568 + t600;
t459 = t315 * t72;
t30 = t311 * t61 + t459;
t28 = qJ(6) * t282 + t30;
t145 = Ifges(7,5) * t324;
t55 = t282 * Ifges(7,6) + t154 * Ifges(7,3) + t145;
t465 = Ifges(6,4) * t324;
t58 = -t154 * Ifges(6,2) + t282 * Ifges(6,6) + t465;
t596 = -t28 * mrSges(7,2) - t30 * mrSges(6,3) + t132 * mrSges(6,1) + t48 * mrSges(7,1) + t55 / 0.2e1 - t58 / 0.2e1;
t422 = qJD(2) * t309;
t384 = qJD(1) * t422;
t278 = t314 * t384;
t409 = qJDD(1) * t309;
t220 = t318 * t409 - t278;
t205 = -qJDD(2) * pkin(2) - t220;
t119 = -pkin(3) * t260 - pkin(9) * t261 + t205;
t408 = qJDD(1) * t310;
t418 = qJD(3) * t313;
t279 = t318 * t384;
t221 = t314 * t409 + t279;
t590 = qJDD(2) * pkin(8) + qJD(3) * t423 + t221;
t101 = -t263 * t418 + t313 * t408 + t317 * t590;
t84 = qJDD(3) * pkin(9) + t101;
t25 = -qJD(4) * t109 + t316 * t119 - t312 * t84;
t15 = pkin(4) * t249 - pkin(10) * t142 + t25;
t413 = qJD(4) * t316;
t415 = qJD(4) * t312;
t24 = t312 * t119 - t182 * t415 + t199 * t413 + t316 * t84;
t21 = pkin(10) * t143 + t24;
t6 = -qJD(5) * t30 + t15 * t315 - t21 * t311;
t3 = -pkin(5) * t243 + qJDD(6) - t6;
t416 = qJD(3) * t317;
t102 = -t263 * t416 - t313 * t590 + t317 * t408;
t85 = -qJDD(3) * pkin(3) - t102;
t49 = -pkin(4) * t143 + t85;
t512 = -t47 / 0.2e1;
t7 = pkin(5) * t47 - qJ(6) * t46 - qJD(6) * t324 + t49;
t595 = mrSges(6,2) * t49 + mrSges(7,2) * t3 - mrSges(6,3) * t6 - mrSges(7,3) * t7 + Ifges(6,4) * t512 + 0.2e1 * t494 * t566 + t511 * t607 + 0.2e1 * t513 * t568;
t575 = -m(7) - m(6);
t593 = -m(4) + t575;
t74 = mrSges(7,1) * t154 - mrSges(7,3) * t324;
t75 = mrSges(6,1) * t154 + mrSges(6,2) * t324;
t592 = t74 + t75;
t276 = t319 * t312;
t277 = t319 * t316;
t345 = t315 * t276 + t277 * t311;
t557 = qJD(5) * t345 + t311 * t602 - t315 * t601;
t185 = t276 * t311 - t277 * t315;
t556 = -qJD(5) * t185 + t311 * t601 + t315 * t602;
t32 = t315 * t71 - t461;
t411 = qJD(5) * t315;
t591 = pkin(4) * t411 - t32;
t471 = mrSges(6,3) * t324;
t128 = mrSges(6,1) * t282 - t471;
t129 = -mrSges(7,1) * t282 + mrSges(7,2) * t324;
t550 = t128 - t129;
t549 = -t197 + (-t390 + t415) * pkin(4);
t362 = -mrSges(5,1) * t316 + mrSges(5,2) * t312;
t337 = m(5) * pkin(3) - t362;
t589 = mrSges(4,1) + t337 + t608;
t361 = t312 * mrSges(5,1) + t316 * mrSges(5,2);
t588 = -t361 + mrSges(3,2) - mrSges(4,3);
t387 = t312 * t416;
t333 = t313 * t413 + t387;
t458 = cos(pkin(11));
t373 = t458 * t318;
t308 = sin(pkin(11));
t446 = t308 * t314;
t229 = -t310 * t446 + t373;
t171 = t308 * t309 * t313 + t229 * t317;
t374 = t458 * t314;
t445 = t308 * t318;
t228 = t310 * t445 + t374;
t587 = -t171 * t312 + t228 * t316;
t227 = t310 * t374 + t445;
t375 = t309 * t458;
t169 = t227 * t317 - t313 * t375;
t226 = -t310 * t373 + t446;
t586 = -t169 * t312 + t226 * t316;
t583 = -Ifges(6,2) * t502 + Ifges(7,3) * t501 + t489 * t565 + t498 * t567 + t596;
t490 = -t282 / 0.2e1;
t499 = -t324 / 0.2e1;
t582 = -Ifges(6,2) * t501 + Ifges(7,3) * t502 + t490 * t565 + t499 * t567 - t596;
t73 = pkin(5) * t324 + qJ(6) * t154;
t579 = -Ifges(6,4) * t501 - Ifges(7,5) * t502 - t490 * t566 - t499 * t568 + t600;
t578 = pkin(4) * t312 * t575 - m(5) * pkin(8) - t304 * t605 - t305 * t528 + t588;
t301 = pkin(4) * t316 + pkin(3);
t363 = mrSges(4,1) * t317 - mrSges(4,2) * t313;
t401 = m(5) * pkin(9) + mrSges(5,3);
t533 = t313 * t401 + t317 * t337 + mrSges(3,1) + t363;
t569 = mrSges(6,3) + mrSges(7,2);
t577 = t533 + (t575 * t319 + t569) * t313 + (-t575 * t301 + t608) * t317;
t576 = -m(4) - m(5);
t504 = t142 / 0.2e1;
t503 = t143 / 0.2e1;
t493 = t249 / 0.2e1;
t574 = t260 / 0.2e1;
t573 = t261 / 0.2e1;
t572 = t339 / 0.2e1;
t33 = mrSges(6,1) * t243 - mrSges(6,3) * t46;
t34 = -t243 * mrSges(7,1) + t46 * mrSges(7,2);
t562 = t34 - t33;
t35 = -mrSges(6,2) * t243 - mrSges(6,3) * t47;
t36 = -mrSges(7,2) * t47 + mrSges(7,3) * t243;
t561 = t35 + t36;
t253 = t311 * t312 - t315 * t316;
t537 = qJD(4) + qJD(5);
t165 = t537 * t253;
t254 = t311 * t316 + t312 * t315;
t166 = t537 * t254;
t209 = t254 * t419;
t338 = t253 * t317;
t210 = qJD(2) * t338;
t559 = -qJD(6) * t254 + t549 + (t165 - t210) * qJ(6) + (t166 - t209) * pkin(5);
t558 = -qJ(6) * t421 + t557;
t555 = pkin(5) * t421 - t556;
t64 = -mrSges(5,1) * t143 + mrSges(5,2) * t142;
t554 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t261 + t64;
t553 = qJD(6) + t591;
t216 = t254 * t313;
t126 = -mrSges(7,2) * t154 + mrSges(7,3) * t282;
t472 = mrSges(6,3) * t154;
t127 = -mrSges(6,2) * t282 - t472;
t433 = t126 + t127;
t251 = t316 * t267;
t438 = t313 * t316;
t162 = -pkin(10) * t438 + t251 + (-pkin(8) * t312 - pkin(4)) * t317;
t295 = pkin(8) * t435;
t208 = t312 * t267 + t295;
t441 = t312 * t313;
t176 = -pkin(10) * t441 + t208;
t548 = t311 * t162 + t315 * t176;
t434 = t317 * t318;
t192 = (-t312 * t434 + t314 * t316) * t424;
t259 = t364 * qJD(3);
t402 = pkin(8) * t418;
t426 = t316 * t259 + t312 * t402;
t547 = -qJD(4) * t208 - t192 + t426;
t120 = t312 * t259 + t267 * t413 + (-t313 * t417 - t317 * t415) * pkin(8);
t440 = t312 * t314;
t193 = (t316 * t434 + t440) * t424;
t546 = -t193 + t120;
t248 = Ifges(5,4) * t339;
t136 = Ifges(5,1) * t252 + Ifges(5,5) * t293 - t248;
t302 = Ifges(4,4) * t419;
t545 = Ifges(4,1) * t421 + Ifges(4,5) * qJD(3) + t316 * t136 + t302;
t400 = mrSges(4,3) * t421;
t544 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t339 + t252 * mrSges(5,2) + t400;
t414 = qJD(4) * t313;
t334 = -t312 * t414 + t316 * t416;
t543 = t243 * t564 + t46 * t566 + t47 * t565;
t443 = t309 * t317;
t234 = t310 * t313 + t314 * t443;
t442 = t309 * t318;
t157 = t234 * t304 + t305 * t442;
t397 = t304 * t442;
t158 = t234 * t305 - t397;
t542 = t157 * t599 + t598 * t158;
t541 = t101 * t317 - t102 * t313;
t540 = t24 * t316 - t25 * t312;
t105 = t171 * t304 - t228 * t305;
t106 = t171 * t305 + t228 * t304;
t539 = t105 * t599 + t598 * t106;
t103 = t169 * t304 - t226 * t305;
t104 = t169 * t305 + t226 * t304;
t538 = t103 * t599 + t598 * t104;
t536 = Ifges(5,5) * t252 - Ifges(5,6) * t339 + Ifges(5,3) * t293 + t154 * t565 + t282 * t564 + t324 * t566;
t531 = -t401 + mrSges(4,2) - t569;
t76 = t344 * qJD(3) + (-t295 + (pkin(10) * t313 - t267) * t312) * qJD(4) + t426;
t83 = -pkin(10) * t333 + t120;
t23 = -qJD(5) * t548 - t311 * t83 + t315 * t76;
t527 = -m(5) * t181 - t544;
t526 = -t25 * mrSges(5,1) + t24 * mrSges(5,2);
t412 = qJD(5) * t311;
t5 = t311 * t15 + t315 * t21 + t61 * t411 - t412 * t72;
t2 = qJ(6) * t243 + qJD(6) * t282 + t5;
t523 = -t6 * mrSges(6,1) + t3 * mrSges(7,1) + t5 * mrSges(6,2) - t2 * mrSges(7,3);
t470 = Ifges(4,4) * t313;
t356 = t317 * Ifges(4,2) + t470;
t522 = t26 * mrSges(7,1) + t30 * mrSges(6,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(2) * t356 / 0.2e1 - t28 * mrSges(7,3) - t29 * mrSges(6,1);
t516 = mrSges(6,1) * t49 + mrSges(7,1) * t7 - mrSges(7,2) * t2 - mrSges(6,3) * t5 + 0.2e1 * Ifges(7,3) * t511 - t46 * Ifges(6,4) / 0.2e1 - t243 * Ifges(6,6) / 0.2e1 + t607 * t513 + (t565 + Ifges(7,6)) * t494 + (-t512 + t511) * Ifges(6,2);
t320 = qJD(2) ^ 2;
t510 = Ifges(5,1) * t504 + Ifges(5,4) * t503 + Ifges(5,5) * t493;
t492 = t252 / 0.2e1;
t487 = pkin(4) * t252;
t486 = pkin(4) * t311;
t484 = pkin(4) * t315;
t306 = t313 * pkin(8);
t479 = -qJD(2) / 0.2e1;
t478 = qJD(4) / 0.2e1;
t473 = mrSges(5,3) * t252;
t469 = Ifges(4,4) * t317;
t468 = Ifges(5,4) * t252;
t467 = Ifges(5,4) * t312;
t466 = Ifges(5,4) * t316;
t460 = t313 * t85;
t444 = t309 * t314;
t439 = t312 * t317;
t425 = pkin(2) * t442 + pkin(8) * t444;
t262 = pkin(4) * t441 + t306;
t420 = qJD(2) * t314;
t303 = pkin(8) * t416;
t399 = mrSges(4,3) * t419;
t398 = t313 * t442;
t396 = Ifges(5,5) * t142 + Ifges(5,6) * t143 + Ifges(5,3) * t249;
t198 = pkin(4) * t333 + t303;
t392 = t309 * t420;
t391 = t318 * t422;
t135 = -Ifges(5,2) * t339 + Ifges(5,6) * t293 + t468;
t385 = -t312 * t135 / 0.2e1;
t372 = t410 / 0.2e1;
t371 = -t103 * pkin(5) + qJ(6) * t104;
t370 = -t105 * pkin(5) + qJ(6) * t106;
t369 = -t157 * pkin(5) + qJ(6) * t158;
t366 = t586 * pkin(4);
t365 = t587 * pkin(4);
t358 = Ifges(5,1) * t316 - t467;
t357 = Ifges(5,1) * t312 + t466;
t355 = -Ifges(5,2) * t312 + t466;
t354 = Ifges(5,2) * t316 + t467;
t353 = Ifges(4,5) * t317 - Ifges(4,6) * t313;
t352 = Ifges(5,5) * t316 - Ifges(5,6) * t312;
t351 = Ifges(5,5) * t312 + Ifges(5,6) * t316;
t77 = t162 * t315 - t176 * t311;
t174 = -t234 * t312 - t316 * t442;
t343 = -t234 * t316 + t312 * t442;
t346 = t315 * t174 + t311 * t343;
t81 = t174 * t311 - t315 * t343;
t233 = -t310 * t317 + t313 * t444;
t342 = t181 * t361;
t264 = -qJD(2) * pkin(2) - t393;
t341 = t264 * (mrSges(4,1) * t313 + mrSges(4,2) * t317);
t340 = t313 * (Ifges(4,1) * t317 - t470);
t22 = t162 * t411 - t176 * t412 + t311 * t76 + t315 * t83;
t335 = t174 * pkin(4);
t332 = t339 * mrSges(5,3);
t327 = Ifges(5,5) * t313 + t317 * t358;
t326 = Ifges(5,6) * t313 + t317 * t355;
t325 = Ifges(5,3) * t313 + t317 * t352;
t323 = -t523 + t543;
t300 = -pkin(5) - t484;
t297 = qJ(6) + t486;
t272 = -qJD(3) * mrSges(4,2) + t399;
t256 = t363 * qJD(2);
t222 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t260;
t219 = t228 * pkin(2);
t218 = t226 * pkin(2);
t217 = t253 * t313;
t207 = -pkin(8) * t439 + t251;
t195 = mrSges(5,1) * t293 - t473;
t194 = -t293 * mrSges(5,2) - t332;
t177 = -mrSges(4,1) * t260 + mrSges(4,2) * t261;
t173 = -qJD(3) * t233 + t317 * t391;
t172 = qJD(3) * t234 + t313 * t391;
t170 = -t229 * t313 + t308 * t443;
t168 = -t227 * t313 - t317 * t375;
t147 = pkin(5) * t253 - qJ(6) * t254 - t301;
t115 = pkin(5) * t216 + qJ(6) * t217 + t262;
t114 = -mrSges(5,2) * t249 + mrSges(5,3) * t143;
t113 = mrSges(5,1) * t249 - mrSges(5,3) * t142;
t111 = t192 * t311 + t193 * t315;
t110 = -t315 * t192 + t193 * t311;
t87 = -t412 * t441 + (t438 * t537 + t387) * t315 + t334 * t311;
t86 = -qJD(3) * t338 - t216 * t537;
t70 = pkin(5) * t317 - t77;
t69 = -qJ(6) * t317 + t548;
t66 = qJD(4) * t174 + t173 * t316 + t312 * t392;
t65 = qJD(4) * t343 - t173 * t312 + t316 * t392;
t53 = t487 + t73;
t50 = Ifges(5,4) * t142 + Ifges(5,2) * t143 + Ifges(5,6) * t249;
t31 = t311 * t71 + t459;
t27 = pkin(5) * t87 - qJ(6) * t86 + qJD(6) * t217 + t198;
t19 = -pkin(5) * t418 - t23;
t18 = qJ(6) * t418 - qJD(6) * t317 + t22;
t17 = qJD(5) * t81 + t311 * t66 - t315 * t65;
t16 = qJD(5) * t346 + t311 * t65 + t315 * t66;
t13 = mrSges(6,1) * t47 + mrSges(6,2) * t46;
t12 = mrSges(7,1) * t47 - mrSges(7,3) * t46;
t1 = [m(2) * qJDD(1) + t174 * t113 - t343 * t114 + t173 * t272 + t66 * t194 + t65 * t195 + t234 * t222 + t561 * t81 - t562 * t346 - t550 * t17 + t433 * t16 + (t12 + t13 + t554) * t233 + (t544 + t592) * t172 + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t320 - t177) * t318 + (-mrSges(3,1) * t320 - mrSges(3,2) * qJDD(2) - qJD(2) * t256) * t314) * t309 + (-m(2) - m(3) + t575 + t576) * g(3) + m(4) * (t101 * t234 - t102 * t233 - t172 * t196 + t173 * t197 + (-t205 * t318 + t264 * t420) * t309) + m(3) * (qJDD(1) * t310 ^ 2 + (t220 * t318 + t221 * t314) * t309) + m(7) * (t16 * t28 + t17 * t26 + t172 * t48 + t2 * t81 + t233 * t7 - t3 * t346) + m(6) * (t132 * t172 + t16 * t30 - t17 * t29 + t233 * t49 + t346 * t6 + t5 * t81) + m(5) * (t108 * t65 + t109 * t66 + t172 * t181 + t174 * t25 + t233 * t85 - t24 * t343); (t108 * mrSges(5,1) - t109 * mrSges(5,2) + Ifges(6,6) * t502 + Ifges(7,6) * t501 + t489 * t564 + t498 * t566 - t522 - t197 * mrSges(4,3) + t536 / 0.2e1) * t418 + (t279 - t221) * mrSges(3,2) + (-t108 * t334 - t109 * t333 - t24 * t441 - t25 * t438) * mrSges(5,3) + qJD(3) ^ 2 * t353 / 0.2e1 + (-m(6) * t132 - m(7) * t48 + t527 - t592) * t313 * t393 - t339 * (qJD(3) * t326 - t354 * t414) / 0.2e1 + t293 * (qJD(3) * t325 - t351 * t414) / 0.2e1 - t272 * t402 - t205 * t363 + t256 * t394 + (t115 * t7 + t2 * t69 + t27 * t48 + t3 * t70 + (-t111 + t18) * t28 + (-t110 + t19) * t26) * m(7) - (t135 * t316 + t136 * t312) * t414 / 0.2e1 + t181 * (mrSges(5,1) * t333 + mrSges(5,2) * t334) + (-t196 * t416 + t541) * mrSges(4,3) + t340 * t372 + (t278 + t220) * mrSges(3,1) + t438 * t510 + (t576 * t425 - t569 * t398 + t575 * (-t319 * t398 + t425) + t528 * (-t305 * t444 + t317 * t397) + (t575 * (pkin(4) * t440 + t301 * t434) - t605 * (t304 * t314 + t305 * t434) - t533 * t318 + t588 * t314) * t309) * g(3) + (qJD(3) * t327 - t357 * t414) * t492 + Ifges(3,3) * qJDD(2) + qJDD(3) * (Ifges(4,5) * t313 + Ifges(4,6) * t317) + (-Ifges(5,6) * t503 - Ifges(5,5) * t504 - Ifges(7,6) * t511 + (-Ifges(4,2) * t313 + t469) * t372 - t272 * t393 + Ifges(4,4) * t573 + Ifges(4,2) * t574 - Ifges(5,3) * t493 + pkin(8) * t222 - Ifges(6,6) * t512 - t566 * t513 - t564 * t494 + t523 + t526) * t317 + (Ifges(4,1) * t261 + Ifges(4,4) * t574 + t352 * t493 + t355 * t503 + t358 * t504) * t313 + t554 * t306 + t262 * t13 - t433 * t111 + t544 * t303 + t545 * t416 / 0.2e1 + t546 * t194 + t547 * t195 + (t207 * t25 + t208 * t24 + (t181 * t416 + t460) * pkin(8) + t546 * t109 + t547 * t108) * m(5) + t550 * t110 + (-pkin(2) * t205 + ((-t196 * t317 - t197 * t313) * qJD(3) + t541) * pkin(8) - (t264 * t314 + (-t196 * t313 + t197 * t317) * t318) * t424) * m(4) - (t396 + t543) * t317 / 0.2e1 + t207 * t113 + t208 * t114 + t198 * t75 - pkin(2) * t177 + t18 * t126 + t22 * t127 + t23 * t128 + t516 * t216 + t19 * t129 + t115 * t12 + t77 * t33 + t69 * t36 + t70 * t34 + t27 * t74 - t50 * t441 / 0.2e1 - t595 * t217 + t597 * t86 + t583 * t87 + t469 * t573 + t356 * t574 + qJD(3) * t341 + (m(5) * t219 + t593 * (pkin(8) * t229 - t219) + t578 * t229 + t577 * t228) * g(1) + (m(5) * t218 + t593 * (pkin(8) * t227 - t218) + t578 * t227 + t577 * t226) * g(2) + (t132 * t198 + t262 * t49 + t5 * t548 + t6 * t77 + (-t111 + t22) * t30 + (t110 + t23) * t29) * m(6) + t548 * t35 + t385 * t416 + t361 * t460; (t325 * t479 + t352 * t478) * t293 + (t327 * t479 + t358 * t478) * t252 + (-t272 + t399) * t196 - t355 * t331 / 0.2e1 + t136 * t413 / 0.2e1 - t353 * t410 / 0.2e1 - t320 * t340 / 0.2e1 + t85 * t362 + (t575 * (t170 * t301 - t171 * t319) + t531 * t171 - t589 * t170) * g(1) + (t575 * (t168 * t301 - t169 * t319) + t531 * t169 - t589 * t168) * g(2) + t135 * t390 / 0.2e1 + t354 * t503 + t357 * t504 + t312 * t510 + (-pkin(3) * t85 - t108 * t130 - t109 * t131) * m(5) + (t385 + t342) * qJD(4) + t351 * t493 + Ifges(4,3) * qJDD(3) + (-t108 * (mrSges(5,1) * t313 - mrSges(5,3) * t435) - t109 * (-mrSges(5,2) * t313 - mrSges(5,3) * t439) + t326 * t572 - t341) * qJD(2) + t316 * t50 / 0.2e1 + (Ifges(6,6) * t501 + Ifges(7,6) * t502 + t490 * t564 + t499 * t566 + t522) * t421 + t555 * t129 + t556 * t128 + t557 * t127 + t558 * t126 + t559 * t74 + t561 * t185 - t301 * t13 + Ifges(4,6) * t260 + Ifges(4,5) * t261 - (-Ifges(4,2) * t421 + t302 + t545) * t419 / 0.2e1 + t549 * t75 + (m(5) * ((-t108 * t316 - t109 * t312) * qJD(4) + t540) - t195 * t413 - t194 * t415 + t316 * t114 - t312 * t113) * pkin(9) + (-t108 * t413 - t109 * t415 + t540) * mrSges(5,3) - t536 * t421 / 0.2e1 - t131 * t194 - t130 * t195 + (t400 + t527) * t197 + t147 * t12 + t516 * t253 - t101 * mrSges(4,2) + t102 * mrSges(4,1) - pkin(3) * t64 + t595 * t254 + t579 * t210 - t597 * t165 + t582 * t209 + t583 * t166 + (t575 * (-t233 * t301 - t234 * t319) + t531 * t234 + t589 * t233) * g(3) + (t132 * t549 + t185 * t5 + t29 * t556 + t30 * t557 - t301 * t49 + t345 * t6) * m(6) + (t147 * t7 + t185 * t2 + t26 * t555 + t28 * t558 - t3 * t345 + t48 * t559) * m(7) - t562 * t345 - t342 * t419; (t2 * t297 + t28 * t553 + t3 * t300 - t48 * t53) * m(7) - t75 * t487 - t526 + (-Ifges(5,2) * t252 + t136 - t248) * t572 + (-t132 * t487 + t29 * t31 - t30 * t32 + (t311 * t5 + t315 * t6 + (-t29 * t311 + t30 * t315) * qJD(5)) * pkin(4)) * m(6) + t323 + (t195 + t473) * t109 + t396 - t293 * (-Ifges(5,5) * t339 - Ifges(5,6) * t252) / 0.2e1 - t181 * (t252 * mrSges(5,1) - mrSges(5,2) * t339) + (-t194 - t332) * t108 - t252 * (-Ifges(5,1) * t339 - t468) / 0.2e1 + t33 * t484 + t35 * t486 + t135 * t492 + t553 * t126 + t300 * t34 + t297 * t36 - t53 * t74 + t579 * t154 + t582 * t324 + (-t586 * mrSges(5,1) - (-t169 * t316 - t226 * t312) * mrSges(5,2) - m(7) * (t366 + t371) - m(6) * t366 + t538) * g(2) + (-t587 * mrSges(5,1) - (-t171 * t316 - t228 * t312) * mrSges(5,2) - m(7) * (t365 + t370) - m(6) * t365 + t539) * g(1) + (-m(7) * t26 + t550) * (-pkin(4) * t412 + t31) + t591 * t127 + (-m(7) * (t335 + t369) - m(6) * t335 - mrSges(5,1) * t174 - mrSges(5,2) * t343 + t542) * g(3); t542 * g(3) + t323 + t539 * g(1) + t538 * g(2) + t58 * t498 + (Ifges(7,3) * t324 - t464) * t502 + (t550 + t471) * t30 + (-t433 - t472) * t29 + (t154 * t26 + t28 * t324) * mrSges(7,2) - t48 * (mrSges(7,1) * t324 + mrSges(7,3) * t154) - t132 * (mrSges(6,1) * t324 - mrSges(6,2) * t154) + qJD(6) * t126 - t73 * t74 - pkin(5) * t34 + qJ(6) * t36 + (-t154 * t566 + t324 * t565) * t490 + (-pkin(5) * t3 - t370 * g(1) - t371 * g(2) - t369 * g(3) + qJ(6) * t2 - t26 * t30 + t28 * t551 - t48 * t73) * m(7) + (-Ifges(6,2) * t324 - t146 + t560) * t501 + (-t154 * t568 + t145 - t465 + t55) * t499; -t282 * t126 + t324 * t74 + (-g(1) * t105 - g(2) * t103 - g(3) * t157 - t28 * t282 + t324 * t48 + t3) * m(7) + t34;];
tau  = t1;
