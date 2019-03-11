% Calculate vector of inverse dynamics joint torques for
% S6RRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:52:57
% EndTime: 2019-03-09 21:53:33
% DurationCPUTime: 21.22s
% Computational Cost: add. (30148->799), mult. (72306->1055), div. (0->0), fcn. (55077->18), ass. (0->381)
t365 = qJ(2) + qJ(3);
t360 = qJ(4) + t365;
t344 = pkin(11) + t360;
t334 = sin(t344);
t335 = cos(t344);
t368 = sin(qJ(6));
t506 = mrSges(7,2) * t368;
t600 = t334 * t506 + t335 * (m(7) * pkin(10) + mrSges(7,3));
t363 = qJD(2) + qJD(3);
t356 = qJD(4) + t363;
t376 = cos(qJ(2));
t378 = -pkin(8) - pkin(7);
t327 = t378 * t376;
t310 = qJD(1) * t327;
t370 = sin(qJ(3));
t290 = t370 * t310;
t371 = sin(qJ(2));
t326 = t378 * t371;
t309 = qJD(1) * t326;
t296 = qJD(2) * pkin(2) + t309;
t375 = cos(qJ(3));
t242 = t375 * t296 + t290;
t304 = t370 * t376 + t371 * t375;
t289 = t304 * qJD(1);
t283 = t289 * pkin(9);
t203 = t242 - t283;
t192 = pkin(3) * t363 + t203;
t293 = t375 * t310;
t243 = t296 * t370 - t293;
t303 = -t370 * t371 + t375 * t376;
t288 = t303 * qJD(1);
t513 = pkin(9) * t288;
t204 = t243 + t513;
t369 = sin(qJ(4));
t193 = t369 * t204;
t374 = cos(qJ(4));
t137 = t374 * t192 - t193;
t235 = t288 * t369 + t289 * t374;
t227 = qJ(5) * t235;
t115 = t137 - t227;
t109 = pkin(4) * t356 + t115;
t367 = cos(pkin(11));
t195 = t374 * t204;
t138 = t192 * t369 + t195;
t425 = t374 * t288 - t289 * t369;
t482 = qJ(5) * t425;
t116 = t138 + t482;
t366 = sin(pkin(11));
t480 = t116 * t366;
t58 = t109 * t367 - t480;
t56 = -pkin(5) * t356 - t58;
t373 = cos(qJ(6));
t589 = t367 * t235 + t366 * t425;
t153 = t356 * t373 - t368 * t589;
t154 = t356 * t368 + t373 * t589;
t572 = mrSges(6,1) * t356 + mrSges(7,1) * t153 - mrSges(7,2) * t154 - mrSges(6,3) * t589;
t587 = -m(6) * t58 + m(7) * t56 - t572;
t590 = -t235 * t366 + t367 * t425;
t599 = pkin(5) * t589 - pkin(10) * t590;
t454 = qJD(1) * qJD(2);
t313 = qJDD(1) * t376 - t371 * t454;
t314 = qJDD(1) * t371 + t376 * t454;
t395 = t303 * qJD(3);
t207 = qJD(1) * t395 + t313 * t370 + t314 * t375;
t396 = t304 * qJD(3);
t208 = -qJD(1) * t396 + t313 * t375 - t314 * t370;
t127 = qJD(4) * t425 + t207 * t374 + t208 * t369;
t362 = qJDD(2) + qJDD(3);
t355 = qJDD(4) + t362;
t302 = t314 * pkin(7);
t253 = qJDD(2) * pkin(2) - pkin(8) * t314 - t302;
t301 = t313 * pkin(7);
t254 = pkin(8) * t313 + t301;
t150 = -qJD(3) * t243 + t375 * t253 - t254 * t370;
t114 = pkin(3) * t362 - pkin(9) * t207 + t150;
t459 = qJD(3) * t375;
t460 = qJD(3) * t370;
t149 = t370 * t253 + t375 * t254 + t296 * t459 + t310 * t460;
t122 = pkin(9) * t208 + t149;
t39 = -qJD(4) * t138 + t374 * t114 - t122 * t369;
t23 = pkin(4) * t355 - qJ(5) * t127 - qJD(5) * t235 + t39;
t128 = -qJD(4) * t235 - t207 * t369 + t208 * t374;
t457 = qJD(4) * t374;
t458 = qJD(4) * t369;
t38 = t369 * t114 + t374 * t122 + t192 * t457 - t204 * t458;
t25 = qJ(5) * t128 + qJD(5) * t425 + t38;
t10 = t23 * t367 - t25 * t366;
t101 = Ifges(6,4) * t589 + Ifges(6,2) * t590 + Ifges(6,6) * t356;
t102 = Ifges(6,1) * t589 + Ifges(6,4) * t590 + t356 * Ifges(6,5);
t11 = t366 * t23 + t367 * t25;
t73 = t127 * t367 + t128 * t366;
t47 = qJD(6) * t153 + t355 * t368 + t373 * t73;
t48 = -qJD(6) * t154 + t355 * t373 - t368 * t73;
t72 = -t127 * t366 + t128 * t367;
t71 = qJDD(6) - t72;
t14 = Ifges(7,4) * t47 + Ifges(7,2) * t48 + Ifges(7,6) * t71;
t15 = Ifges(7,1) * t47 + Ifges(7,4) * t48 + Ifges(7,5) * t71;
t167 = qJD(6) - t590;
t361 = t376 * pkin(2);
t349 = t361 + pkin(1);
t325 = t349 * qJD(1);
t255 = -pkin(3) * t288 - t325;
t180 = -pkin(4) * t425 + qJD(5) + t255;
t481 = qJDD(1) * pkin(1);
t281 = -pkin(2) * t313 - t481;
t181 = -pkin(3) * t208 + t281;
t92 = -pkin(4) * t128 + qJDD(5) + t181;
t19 = -pkin(5) * t72 - pkin(10) * t73 + t92;
t112 = t367 * t116;
t59 = t366 * t109 + t112;
t57 = pkin(10) * t356 + t59;
t89 = -pkin(5) * t590 - pkin(10) * t589 + t180;
t26 = -t368 * t57 + t373 * t89;
t8 = pkin(10) * t355 + t11;
t2 = qJD(6) * t26 + t19 * t368 + t373 * t8;
t27 = t368 * t89 + t373 * t57;
t3 = -qJD(6) * t27 + t19 * t373 - t368 * t8;
t414 = mrSges(7,1) * t368 + mrSges(7,2) * t373;
t398 = t56 * t414;
t409 = Ifges(7,5) * t373 - Ifges(7,6) * t368;
t496 = Ifges(7,4) * t373;
t411 = -Ifges(7,2) * t368 + t496;
t497 = Ifges(7,4) * t368;
t413 = Ifges(7,1) * t373 - t497;
t456 = qJD(6) * t368;
t432 = -t456 / 0.2e1;
t152 = Ifges(7,4) * t153;
t80 = Ifges(7,1) * t154 + Ifges(7,5) * t167 + t152;
t484 = t373 * t80;
t442 = t484 / 0.2e1;
t455 = qJD(6) * t373;
t501 = mrSges(7,3) * t373;
t502 = mrSges(7,3) * t368;
t510 = t58 * mrSges(6,3);
t526 = mrSges(6,3) * t59;
t528 = t368 / 0.2e1;
t530 = -t356 / 0.2e1;
t537 = -t589 / 0.2e1;
t538 = -t590 / 0.2e1;
t539 = -t167 / 0.2e1;
t541 = -t154 / 0.2e1;
t542 = -t153 / 0.2e1;
t547 = t71 / 0.2e1;
t548 = t48 / 0.2e1;
t549 = t47 / 0.2e1;
t579 = t27 * mrSges(7,2);
t580 = t26 * mrSges(7,1);
t508 = mrSges(7,1) * t373;
t595 = t506 - t508;
t7 = -pkin(5) * t355 - t10;
t78 = Ifges(7,5) * t154 + Ifges(7,6) * t153 + Ifges(7,3) * t167;
t493 = t154 * Ifges(7,4);
t79 = t153 * Ifges(7,2) + t167 * Ifges(7,6) + t493;
t598 = (Ifges(7,5) * t368 + Ifges(7,6) * t373) * t547 + (Ifges(7,2) * t373 + t497) * t548 + t15 * t528 + t2 * t501 + (-t26 * t455 - t27 * t456) * mrSges(7,3) + t7 * t595 - t3 * t502 - t11 * mrSges(6,2) + t10 * mrSges(6,1) + t373 * t14 / 0.2e1 + Ifges(5,5) * t127 + Ifges(5,6) * t128 + (Ifges(7,1) * t368 + t496) * t549 - t38 * mrSges(5,2) + t39 * mrSges(5,1) + t79 * t432 + Ifges(6,6) * t72 + Ifges(6,5) * t73 + (Ifges(5,3) + Ifges(6,3)) * t355 + (t398 + t442) * qJD(6) + (t153 * t411 + t154 * t413 + t167 * t409) * qJD(6) / 0.2e1 + (-t180 * mrSges(6,2) + Ifges(6,1) * t537 + Ifges(6,4) * t538 + Ifges(6,5) * t530 + t26 * t501 + t27 * t502 + t409 * t539 + t411 * t542 + t413 * t541 - t398 - t484 / 0.2e1 + t79 * t528 + t510 - t102 / 0.2e1) * t590 + (t579 - t580 - t180 * mrSges(6,1) - Ifges(6,4) * t537 + Ifges(7,5) * t541 - Ifges(6,2) * t538 - Ifges(6,6) * t530 + Ifges(7,6) * t542 + Ifges(7,3) * t539 + t526 + t101 / 0.2e1 - t78 / 0.2e1) * t589;
t520 = pkin(4) * t235;
t228 = Ifges(5,4) * t425;
t164 = Ifges(5,1) * t235 + Ifges(5,5) * t356 + t228;
t597 = t164 + t228;
t345 = sin(t360);
t346 = cos(t360);
t596 = mrSges(5,1) * t345 + mrSges(6,1) * t334 + mrSges(5,2) * t346 + mrSges(6,2) * t335;
t372 = sin(qJ(1));
t377 = cos(qJ(1));
t594 = g(1) * t377 + g(2) * t372;
t593 = -t346 * mrSges(5,1) - t335 * mrSges(6,1) + t345 * mrSges(5,2) + (mrSges(6,2) - mrSges(7,3)) * t334;
t357 = sin(t365);
t358 = cos(t365);
t592 = mrSges(4,1) * t357 + mrSges(4,2) * t358 + t596;
t583 = -m(7) - m(6);
t498 = Ifges(5,4) * t235;
t163 = Ifges(5,2) * t425 + t356 * Ifges(5,6) + t498;
t582 = t163 / 0.2e1;
t581 = t313 / 0.2e1;
t531 = t355 / 0.2e1;
t527 = t376 / 0.2e1;
t536 = -t425 / 0.2e1;
t578 = mrSges(5,2) * t425;
t577 = Ifges(5,1) * t425;
t247 = t303 * t369 + t304 * t374;
t576 = Ifges(5,5) * t247;
t575 = Ifges(5,5) * t425;
t246 = t303 * t374 - t304 * t369;
t574 = Ifges(5,6) * t246;
t573 = t376 * Ifges(3,2);
t158 = -mrSges(6,2) * t356 + mrSges(6,3) * t590;
t98 = -mrSges(7,2) * t167 + mrSges(7,3) * t153;
t99 = mrSges(7,1) * t167 - mrSges(7,3) * t154;
t407 = -t368 * t99 + t373 * t98;
t567 = -t158 - t407;
t256 = t375 * t326 + t327 * t370;
t229 = -pkin(9) * t304 + t256;
t257 = t370 * t326 - t375 * t327;
t230 = pkin(9) * t303 + t257;
t166 = t369 * t229 + t374 * t230;
t248 = -t309 * t370 + t293;
t209 = t248 - t513;
t249 = t375 * t309 + t290;
t210 = -t283 + t249;
t145 = t369 * t209 + t374 * t210;
t524 = pkin(2) * t375;
t348 = pkin(3) + t524;
t468 = t369 * t370;
t238 = t348 * t457 + (-t370 * t458 + (t374 * t375 - t468) * qJD(3)) * pkin(2);
t566 = t238 - t145;
t144 = t374 * t209 - t210 * t369;
t467 = t370 * t374;
t239 = -t348 * t458 + (-t370 * t457 + (-t369 * t375 - t467) * qJD(3)) * pkin(2);
t565 = t239 - t144;
t421 = t335 * pkin(5) + t334 * pkin(10);
t430 = t358 * mrSges(4,1) - mrSges(4,2) * t357;
t179 = t246 * t366 + t247 * t367;
t250 = qJD(2) * t303 + t395;
t251 = -qJD(2) * t304 - t396;
t146 = qJD(4) * t246 + t250 * t374 + t251 * t369;
t147 = -qJD(4) * t247 - t250 * t369 + t251 * t374;
t94 = t146 * t367 + t147 * t366;
t403 = t179 * t455 + t368 * t94;
t462 = qJD(1) * t376;
t463 = qJD(1) * t371;
t514 = pkin(7) * t376;
t515 = pkin(7) * t371;
t563 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t463) * t514 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t462) * t515;
t562 = t301 * t376 + t302 * t371;
t560 = 0.2e1 * t531;
t559 = t335 * t595 + t593;
t557 = -t430 + t559;
t555 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t364 = -pkin(9) + t378;
t554 = -m(3) * pkin(7) + m(4) * t378 + m(5) * t364 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3);
t323 = -mrSges(3,1) * t376 + mrSges(3,2) * t371;
t343 = pkin(3) * t358;
t464 = t343 + t361;
t553 = mrSges(2,1) + m(5) * (pkin(1) + t464) + m(4) * t349 + t430 + m(3) * pkin(1) - t323 - t593;
t20 = mrSges(7,1) * t71 - mrSges(7,3) * t47;
t21 = -mrSges(7,2) * t71 + mrSges(7,3) * t48;
t552 = m(7) * (t2 * t373 - t3 * t368 + (-t26 * t373 - t27 * t368) * qJD(6)) - t99 * t455 - t98 * t456 + t373 * t21 - t368 * t20;
t408 = -t26 * t368 + t27 * t373;
t551 = m(6) * t59 + m(7) * t408 - t567;
t540 = t154 / 0.2e1;
t535 = -t235 / 0.2e1;
t534 = t235 / 0.2e1;
t532 = t289 / 0.2e1;
t525 = pkin(2) * t371;
t523 = pkin(3) * t289;
t522 = pkin(3) * t357;
t521 = pkin(3) * t374;
t519 = pkin(4) * t345;
t336 = pkin(4) * t346;
t518 = pkin(4) * t366;
t517 = pkin(4) * t367;
t516 = pkin(5) * t334;
t505 = mrSges(4,3) * t288;
t504 = mrSges(5,3) * t137;
t503 = mrSges(5,3) * t425;
t500 = Ifges(3,4) * t371;
t499 = Ifges(3,4) * t376;
t495 = pkin(3) * qJD(4);
t494 = t138 * mrSges(5,3);
t492 = t289 * mrSges(4,3);
t491 = t289 * Ifges(4,4);
t479 = t179 * t368;
t478 = t179 * t373;
t472 = t366 * t369;
t471 = t367 * t369;
t470 = t368 * t372;
t469 = t368 * t377;
t466 = t372 * t373;
t465 = t373 * t377;
t142 = t374 * t203 - t193;
t284 = -pkin(2) * t468 + t374 * t348;
t278 = pkin(4) + t284;
t285 = pkin(2) * t467 + t348 * t369;
t220 = t366 * t278 + t367 * t285;
t347 = pkin(4) + t521;
t280 = pkin(3) * t471 + t366 * t347;
t461 = qJD(2) * t371;
t452 = Ifges(7,5) * t47 + Ifges(7,6) * t48 + Ifges(7,3) * t71;
t353 = pkin(2) * t461;
t448 = t334 * t508;
t352 = pkin(2) * t463;
t441 = t336 + t421;
t440 = t336 + t464;
t439 = qJD(2) * t378;
t435 = -t72 * mrSges(6,1) + t73 * mrSges(6,2);
t236 = -pkin(3) * t251 + t353;
t431 = t454 / 0.2e1;
t141 = -t203 * t369 - t195;
t165 = t374 * t229 - t230 * t369;
t424 = t343 + t441;
t264 = -pkin(3) * t303 - t349;
t423 = t600 * t372;
t422 = t600 * t377;
t189 = t523 + t520;
t299 = -t519 - t522;
t420 = -g(1) * t372 + g(2) * t377;
t419 = mrSges(3,1) * t371 + mrSges(3,2) * t376;
t412 = t500 + t573;
t410 = Ifges(3,5) * t376 - Ifges(3,6) * t371;
t134 = qJ(5) * t246 + t166;
t399 = -qJ(5) * t247 + t165;
t88 = t367 * t134 + t366 * t399;
t178 = -t367 * t246 + t247 * t366;
t196 = -pkin(4) * t246 + t264;
t96 = pkin(5) * t178 - pkin(10) * t179 + t196;
t36 = t368 * t96 + t373 * t88;
t35 = -t368 * t88 + t373 * t96;
t219 = t278 * t367 - t285 * t366;
t130 = -pkin(4) * t147 + t236;
t279 = -pkin(3) * t472 + t347 * t367;
t404 = pkin(1) * t419;
t402 = t179 * t456 - t373 * t94;
t401 = t141 - t482;
t400 = t144 - t482;
t397 = t371 * (Ifges(3,1) * t376 - t500);
t311 = t371 * t439;
t312 = t376 * t439;
t185 = t375 * t311 + t370 * t312 + t326 * t459 + t327 * t460;
t160 = pkin(9) * t251 + t185;
t186 = -qJD(3) * t257 - t311 * t370 + t375 * t312;
t161 = -pkin(9) * t250 + t186;
t81 = t374 * t160 + t369 * t161 + t229 * t457 - t230 * t458;
t272 = t299 - t525;
t385 = m(7) * (t272 - t516) - t448;
t384 = m(7) * (t299 - t516) - t448;
t91 = t189 + t599;
t383 = m(7) * (-t516 - t519) - t448;
t82 = -qJD(4) * t166 - t160 * t369 + t374 * t161;
t381 = -qJ(5) * t146 - qJD(5) * t247 + t82;
t225 = t288 * Ifges(4,2) + t363 * Ifges(4,6) + t491;
t282 = Ifges(4,4) * t288;
t226 = t289 * Ifges(4,1) + t363 * Ifges(4,5) + t282;
t379 = t598 + t577 * t535 + t597 * t536 + t225 * t532 + t242 * t505 + t243 * t492 + t325 * (mrSges(4,1) * t289 + mrSges(4,2) * t288) + t425 * t504 + t575 * t530 - t255 * t578 - t289 * (Ifges(4,1) * t288 - t491) / 0.2e1 + Ifges(4,3) * t362 - t363 * (Ifges(4,5) * t288 - Ifges(4,6) * t289) / 0.2e1 + Ifges(4,6) * t208 + Ifges(4,5) * t207 - t149 * mrSges(4,2) + t150 * mrSges(4,1) - (-Ifges(4,2) * t289 + t226 + t282) * t288 / 0.2e1 + (-t255 * mrSges(5,1) - Ifges(5,4) * t535 - Ifges(5,2) * t536 - Ifges(5,6) * t530 + t494 + t582) * t235;
t359 = -qJ(5) + t364;
t351 = Ifges(3,4) * t462;
t342 = -pkin(5) - t517;
t287 = Ifges(3,1) * t463 + Ifges(3,5) * qJD(2) + t351;
t286 = Ifges(3,6) * qJD(2) + qJD(1) * t412;
t273 = -pkin(5) - t279;
t270 = pkin(1) + t440;
t268 = t335 * t465 + t470;
t267 = -t335 * t469 + t466;
t266 = -t335 * t466 + t469;
t265 = t335 * t470 + t465;
t260 = mrSges(4,1) * t363 - t492;
t259 = -mrSges(4,2) * t363 + t505;
t258 = t352 + t523;
t241 = -mrSges(4,1) * t288 + mrSges(4,2) * t289;
t214 = -pkin(5) - t219;
t212 = mrSges(5,1) * t356 - mrSges(5,3) * t235;
t211 = -mrSges(5,2) * t356 + t503;
t191 = -mrSges(4,2) * t362 + mrSges(4,3) * t208;
t190 = mrSges(4,1) * t362 - mrSges(4,3) * t207;
t182 = t189 + t352;
t175 = -mrSges(5,1) * t425 + mrSges(5,2) * t235;
t129 = -t227 + t145;
t123 = -t227 + t142;
t119 = -mrSges(5,2) * t355 + mrSges(5,3) * t128;
t118 = mrSges(5,1) * t355 - mrSges(5,3) * t127;
t107 = -mrSges(6,1) * t590 + mrSges(6,2) * t589;
t95 = t520 + t599;
t93 = t146 * t366 - t367 * t147;
t90 = t352 + t91;
t75 = t367 * t129 + t366 * t400;
t66 = mrSges(6,1) * t355 - mrSges(6,3) * t73;
t65 = -mrSges(6,2) * t355 + mrSges(6,3) * t72;
t64 = t367 * t123 + t366 * t401;
t61 = t115 * t367 - t480;
t60 = t115 * t366 + t112;
t43 = qJ(5) * t147 + qJD(5) * t246 + t81;
t34 = t368 * t95 + t373 * t61;
t33 = -t368 * t61 + t373 * t95;
t32 = t368 * t90 + t373 * t75;
t31 = -t368 * t75 + t373 * t90;
t30 = t368 * t91 + t373 * t64;
t29 = -t368 * t64 + t373 * t91;
t28 = pkin(5) * t93 - pkin(10) * t94 + t130;
t18 = -mrSges(7,1) * t48 + mrSges(7,2) * t47;
t17 = t366 * t381 + t367 * t43;
t5 = -qJD(6) * t36 - t17 * t368 + t28 * t373;
t4 = qJD(6) * t35 + t17 * t373 + t28 * t368;
t1 = [t587 * (t366 * t43 - t367 * t381) + (t376 * t499 + t397) * t431 + t589 * (Ifges(6,1) * t94 - Ifges(6,4) * t93) / 0.2e1 + t590 * (Ifges(6,4) * t94 - Ifges(6,2) * t93) / 0.2e1 + (t313 * t514 + t314 * t515 + t562) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t562) + (t287 * t527 + t410 * qJD(2) / 0.2e1 - t563) * qJD(2) - t403 * t79 / 0.2e1 + t56 * (mrSges(7,1) * t403 - mrSges(7,2) * t402) + t153 * (-Ifges(7,4) * t402 - Ifges(7,2) * t403 + Ifges(7,6) * t93) / 0.2e1 + t167 * (-Ifges(7,5) * t402 - Ifges(7,6) * t403 + Ifges(7,3) * t93) / 0.2e1 + (Ifges(5,1) * t146 + Ifges(5,4) * t147) * t534 + (-Ifges(7,1) * t402 - Ifges(7,4) * t403 + Ifges(7,5) * t93) * t540 + (Ifges(3,4) * t314 + Ifges(3,2) * t313) * t527 + (-mrSges(4,1) * t281 + mrSges(4,3) * t149 + Ifges(4,4) * t207 + Ifges(4,2) * t208 + Ifges(4,6) * t362) * t303 + (Ifges(4,1) * t250 + Ifges(4,4) * t251) * t532 + t147 * t494 + m(6) * (t11 * t88 + t130 * t180 + t17 * t59 + t196 * t92) + m(7) * (t2 * t36 + t26 * t5 + t27 * t4 + t3 * t35) + (mrSges(4,2) * t281 - mrSges(4,3) * t150 + Ifges(4,1) * t207 + Ifges(4,4) * t208 + Ifges(4,5) * t362) * t304 + (-t242 * t250 + t243 * t251) * mrSges(4,3) + t314 * t499 / 0.2e1 + (t92 * mrSges(6,2) - t10 * mrSges(6,3) + Ifges(6,1) * t73 + Ifges(6,4) * t72 + Ifges(6,5) * t560 + t409 * t547 + t411 * t548 + t413 * t549 + t7 * t414 + t80 * t432) * t179 + (Ifges(7,3) * t547 + Ifges(7,6) * t548 + Ifges(7,5) * t549 + t452 / 0.2e1 - Ifges(6,4) * t73 - Ifges(6,2) * t72 - t11 * mrSges(6,3) + t92 * mrSges(6,1) - t560 * Ifges(6,6) + t555) * t178 + (-m(6) * t10 + m(7) * t7 + t18 - t66) * (t134 * t366 - t367 * t399) + (-mrSges(5,1) * t264 + Ifges(5,4) * t247 + Ifges(5,2) * t246) * t128 + (t246 * t38 - t247 * t39) * mrSges(5,3) + m(5) * (t137 * t82 + t138 * t81 + t165 * t39 + t166 * t38 + t181 * t264 + t236 * t255) + (Ifges(5,5) * t146 + Ifges(6,5) * t94 + Ifges(5,6) * t147 - Ifges(6,6) * t93) * t356 / 0.2e1 + t425 * (Ifges(5,4) * t146 + Ifges(5,2) * t147) / 0.2e1 - t94 * t510 - t146 * t504 + Ifges(2,3) * qJDD(1) - t323 * t481 + t15 * t478 / 0.2e1 - t14 * t479 / 0.2e1 - t286 * t461 / 0.2e1 - t404 * t454 + t363 * (Ifges(4,5) * t250 + Ifges(4,6) * t251) / 0.2e1 - t349 * (-mrSges(4,1) * t208 + mrSges(4,2) * t207) - t325 * (-mrSges(4,1) * t251 + mrSges(4,2) * t250) - pkin(1) * (-mrSges(3,1) * t313 + mrSges(3,2) * t314) + t196 * t435 + t288 * (Ifges(4,4) * t250 + Ifges(4,2) * t251) / 0.2e1 + t251 * t225 / 0.2e1 + t255 * (-mrSges(5,1) * t147 + mrSges(5,2) * t146) + t256 * t190 + t257 * t191 + t185 * t259 + t186 * t260 + t181 * (-mrSges(5,1) * t246 + mrSges(5,2) * t247) + t250 * t226 / 0.2e1 + t236 * t175 + t81 * t211 + t82 * t212 + t180 * (mrSges(6,1) * t93 + mrSges(6,2) * t94) + t146 * t164 / 0.2e1 + t165 * t118 + t166 * t119 + t17 * t158 + t130 * t107 - t93 * t101 / 0.2e1 + t94 * t102 / 0.2e1 + t4 * t98 + t5 * t99 + t93 * t78 / 0.2e1 + m(4) * (t149 * t257 + t150 * t256 + t185 * t243 + t186 * t242 - t281 * t349 - t325 * t353) + t241 * t353 + t35 * t20 + t36 * t21 + (t574 / 0.2e1 + t576 / 0.2e1) * t355 + (t574 + t576) * t531 - t93 * t526 - t93 * t579 + t94 * t442 + t93 * t580 + t412 * t581 + t147 * t582 + (-mrSges(3,1) * t515 - mrSges(3,2) * t514 + 0.2e1 * Ifges(3,6) * t527) * qJDD(2) + (Ifges(3,1) * t314 + Ifges(3,4) * t581 + Ifges(3,5) * qJDD(2) - t431 * t573) * t371 + (-t268 * mrSges(7,1) - t267 * mrSges(7,2) + t583 * (t377 * t270 - t359 * t372) + t554 * t372 + (-m(7) * t421 - t553) * t377) * g(2) + (-t266 * mrSges(7,1) - t265 * mrSges(7,2) + (-t359 * t583 + t554) * t377 + (-m(7) * (-t270 - t421) + m(6) * t270 + t553) * t372) * g(1) + t88 * t65 + (mrSges(5,2) * t264 + Ifges(5,1) * t247 + Ifges(5,4) * t246) * t127 + (-t2 * t479 + t26 * t402 - t27 * t403 - t3 * t478) * mrSges(7,3); (t563 + (t404 - t397 / 0.2e1) * qJD(1)) * qJD(1) + t587 * ((-t239 + t400) * t367 + (-t129 + t238) * t366) + t190 * t524 + (t214 * t7 - t26 * t31 - t27 * t32) * m(7) + (t10 * t219 + t11 * t220 - t180 * t182 - t59 * t75) * m(6) + t551 * (t238 * t367 + t239 * t366) + t552 * (pkin(10) + t220) - (-Ifges(3,2) * t463 + t287 + t351) * t462 / 0.2e1 + t594 * (m(4) * t525 - m(6) * t272 - m(5) * (-t522 - t525) + t419 + t592) + (-m(4) * t361 - m(5) * t464 - m(6) * t440 + t323 - m(7) * (t361 + t424) + t557) * g(3) + (m(4) * (t149 * t370 + t150 * t375 + (-t242 * t370 + t243 * t375) * qJD(3)) - t260 * t460 + t370 * t191 + t259 * t459) * pkin(2) + t379 + t286 * t463 / 0.2e1 - t410 * t454 / 0.2e1 - t241 * t352 + Ifges(3,6) * t313 + Ifges(3,5) * t314 - t301 * mrSges(3,2) - t302 * mrSges(3,1) + Ifges(3,3) * qJDD(2) + t284 * t118 + t285 * t119 - t258 * t175 - t249 * t259 - t248 * t260 + t214 * t18 + t219 * t66 + t220 * t65 - t182 * t107 - t75 * t158 - t32 * t98 - t31 * t99 - m(4) * (t242 * t248 + t243 * t249 - t325 * t352) + t565 * t212 + t566 * t211 + (t137 * t565 + t138 * t566 - t255 * t258 + t284 * t39 + t285 * t38) * m(5) - g(1) * (t377 * t385 + t422) - g(2) * (t372 * t385 + t423); (-t26 * t29 - t27 * t30 + t273 * t7) * m(7) + (t10 * t279 + t11 * t280 - t180 * t189 - t59 * t64) * m(6) + t587 * (-t123 * t366 + t367 * t401 + (t366 * t374 + t471) * t495) + t551 * (t367 * t374 - t472) * t495 + t118 * t521 + t552 * (pkin(10) + t280) - t175 * t523 - m(5) * (t137 * t141 + t138 * t142 + t255 * t523) + t594 * (m(5) * t522 - m(6) * t299 + t592) + (-m(5) * t343 - m(6) * (t336 + t343) - m(7) * t424 + t557) * g(3) + (m(5) * (t369 * t38 + t374 * t39 + (-t137 * t369 + t138 * t374) * qJD(4)) - t212 * t458 + t369 * t119 + t211 * t457) * pkin(3) + t379 + t279 * t66 + t280 * t65 + t273 * t18 - t242 * t259 + t243 * t260 - t142 * t211 - t141 * t212 - t189 * t107 - t64 * t158 - t30 * t98 - t29 * t99 - g(1) * (t377 * t384 + t422) - g(2) * (t372 * t384 + t423); t598 + ((t10 * t367 + t11 * t366) * pkin(4) - t180 * t520 + t58 * t60 - t59 * t61) * m(6) + t163 * t534 + (m(6) * t519 + t596) * t594 + (-Ifges(5,2) * t235 + t597) * t536 + t66 * t517 + t65 * t518 + t235 * t494 + t552 * (pkin(10) + t518) + (-m(6) * t336 - m(7) * t441 + t559) * g(3) + (-t26 * t33 - t27 * t34 + t342 * t7 - t56 * t60) * m(7) + (t503 - t211) * t137 - t107 * t520 + t342 * t18 + t138 * t212 - t61 * t158 - t34 * t98 - t33 * t99 + t572 * t60 + (-Ifges(5,6) * t235 + t575) * t530 - g(1) * (t377 * t383 + t422) - g(2) * (t372 * t383 + t423) + (-t498 + t577) * t535 - t255 * (mrSges(5,1) * t235 + t578); t373 * t20 + t368 * t21 + t572 * t589 + t407 * qJD(6) + t567 * t590 + t435 + (t167 * t408 + t2 * t368 + t3 * t373 - t56 * t589 + t420) * m(7) + (t58 * t589 - t59 * t590 + t420 + t92) * m(6); -t56 * (mrSges(7,1) * t154 + mrSges(7,2) * t153) + (Ifges(7,1) * t153 - t493) * t541 + t79 * t540 + (Ifges(7,5) * t153 - Ifges(7,6) * t154) * t539 - t26 * t98 + t27 * t99 - g(1) * (mrSges(7,1) * t267 - mrSges(7,2) * t268) - g(2) * (-mrSges(7,1) * t265 + mrSges(7,2) * t266) + g(3) * t414 * t334 + (t153 * t26 + t154 * t27) * mrSges(7,3) + t452 + (-Ifges(7,2) * t154 + t152 + t80) * t542 + t555;];
tau  = t1;
