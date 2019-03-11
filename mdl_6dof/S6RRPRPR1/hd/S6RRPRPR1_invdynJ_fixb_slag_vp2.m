% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:07:55
% EndTime: 2019-03-09 10:08:35
% DurationCPUTime: 25.19s
% Computational Cost: add. (20424->802), mult. (47717->1051), div. (0->0), fcn. (37132->18), ass. (0->377)
t346 = sin(pkin(11));
t348 = cos(pkin(11));
t351 = sin(qJ(6));
t355 = cos(qJ(6));
t382 = t346 * t351 - t348 * t355;
t263 = t382 * qJD(6);
t347 = sin(pkin(10));
t353 = sin(qJ(2));
t356 = cos(qJ(2));
t457 = cos(pkin(10));
t281 = -t347 * t353 + t356 * t457;
t367 = qJD(1) * t281;
t489 = cos(qJ(4));
t246 = t489 * t367;
t407 = t457 * t353;
t430 = qJD(1) * t356;
t260 = -qJD(1) * t407 - t347 * t430;
t352 = sin(qJ(4));
t205 = t260 * t352 + t246;
t570 = t382 * t205;
t584 = t263 - t570;
t284 = t346 * t355 + t348 * t351;
t264 = t284 * qJD(6);
t571 = t284 * t205;
t583 = t264 - t571;
t343 = pkin(11) + qJ(6);
t334 = sin(t343);
t336 = cos(t343);
t392 = -mrSges(6,1) * t348 + mrSges(6,2) * t346;
t582 = -mrSges(7,1) * t336 + mrSges(7,2) * t334 + t392;
t198 = qJD(6) - t205;
t344 = qJD(2) + qJD(4);
t364 = t352 * t367;
t361 = -t260 * t489 + t364;
t184 = t344 * t346 + t348 * t361;
t402 = t348 * t344 - t346 * t361;
t138 = t184 * t355 + t351 * t402;
t466 = Ifges(7,4) * t138;
t563 = -t184 * t351 + t355 * t402;
t63 = Ifges(7,2) * t563 + Ifges(7,6) * t198 + t466;
t514 = t63 / 0.2e1;
t131 = Ifges(7,4) * t563;
t64 = Ifges(7,1) * t138 + Ifges(7,5) * t198 + t131;
t512 = t64 / 0.2e1;
t573 = t205 * t346;
t580 = pkin(5) * t573;
t579 = pkin(9) * t573;
t578 = -mrSges(6,3) - mrSges(7,3);
t345 = qJ(2) + pkin(10);
t338 = qJ(4) + t345;
t324 = sin(t338);
t325 = cos(t338);
t478 = t348 * pkin(5);
t326 = pkin(4) + t478;
t350 = -pkin(9) - qJ(5);
t577 = m(7) * t325 * t350 + (m(7) * t326 - t582) * t324;
t572 = t205 * t348;
t576 = pkin(5) * t361 - pkin(9) * t572;
t425 = qJD(1) * qJD(2);
t411 = t353 * t425;
t424 = qJDD(1) * t356;
t293 = -t411 + t424;
t294 = qJDD(1) * t353 + t356 * t425;
t221 = t293 * t457 - t347 * t294;
t222 = t347 * t293 + t294 * t457;
t426 = qJD(4) * t352;
t128 = qJD(4) * t246 + t352 * t221 + t222 * t489 + t260 * t426;
t412 = qJD(4) * t489;
t129 = qJD(4) * t364 - t489 * t221 + t352 * t222 - t260 * t412;
t349 = -qJ(3) - pkin(7);
t311 = t349 * t356;
t292 = qJD(1) * t311;
t265 = t347 * t292;
t309 = t349 * t353;
t291 = qJD(1) * t309;
t274 = qJD(2) * pkin(2) + t291;
t212 = t457 * t274 + t265;
t483 = pkin(8) * t260;
t178 = qJD(2) * pkin(3) + t212 + t483;
t408 = t457 * t292;
t213 = t347 * t274 - t408;
t366 = pkin(8) * t367;
t180 = t366 + t213;
t117 = t352 * t178 + t180 * t489;
t109 = t344 * qJ(5) + t117;
t223 = -qJD(1) * pkin(1) - pkin(2) * t430 - pkin(3) * t367 + qJD(3);
t132 = -pkin(4) * t205 - qJ(5) * t361 + t223;
t70 = -t109 * t346 + t348 * t132;
t51 = -pkin(5) * t205 - pkin(9) * t184 + t70;
t71 = t348 * t109 + t346 * t132;
t57 = pkin(9) * t402 + t71;
t17 = -t351 * t57 + t355 * t51;
t18 = t351 * t51 + t355 * t57;
t341 = qJDD(2) + qJDD(4);
t106 = -t128 * t346 + t341 * t348;
t278 = t294 * pkin(7);
t428 = qJD(3) * t353;
t210 = qJDD(2) * pkin(2) - qJ(3) * t294 - qJD(1) * t428 - t278;
t329 = pkin(7) * t424;
t429 = qJD(2) * t353;
t418 = pkin(7) * t429;
t427 = qJD(3) * t356;
t217 = qJ(3) * t293 + t329 + (-t418 + t427) * qJD(1);
t165 = t457 * t210 - t217 * t347;
t123 = qJDD(2) * pkin(3) - pkin(8) * t222 + t165;
t166 = t347 * t210 + t457 * t217;
t133 = pkin(8) * t221 + t166;
t49 = t352 * t123 + t489 * t133 + t178 * t412 - t180 * t426;
t46 = qJ(5) * t341 + qJD(5) * t344 + t49;
t455 = qJDD(1) * pkin(1);
t251 = -pkin(2) * t293 + qJDD(3) - t455;
t181 = -pkin(3) * t221 + t251;
t54 = pkin(4) * t129 - qJ(5) * t128 - qJD(5) * t361 + t181;
t15 = t346 * t54 + t348 * t46;
t11 = pkin(9) * t106 + t15;
t107 = t128 * t348 + t341 * t346;
t14 = -t346 * t46 + t348 * t54;
t6 = pkin(5) * t129 - pkin(9) * t107 + t14;
t2 = qJD(6) * t17 + t11 * t355 + t351 * t6;
t3 = -qJD(6) * t18 - t11 * t351 + t355 * t6;
t50 = t123 * t489 - t352 * t133 - t178 * t426 - t180 * t412;
t47 = -t341 * pkin(4) + qJDD(5) - t50;
t30 = -t106 * pkin(5) + t47;
t354 = sin(qJ(1));
t357 = cos(qJ(1));
t385 = -t14 * t346 + t15 * t348;
t44 = Ifges(6,4) * t107 + Ifges(6,2) * t106 + Ifges(6,6) * t129;
t441 = t325 * t357;
t442 = t325 * t354;
t467 = Ifges(6,4) * t348;
t468 = Ifges(6,4) * t346;
t490 = t348 / 0.2e1;
t499 = t198 / 0.2e1;
t500 = -t198 / 0.2e1;
t504 = t138 / 0.2e1;
t505 = -t138 / 0.2e1;
t506 = t563 / 0.2e1;
t507 = -t563 / 0.2e1;
t508 = t129 / 0.2e1;
t127 = qJDD(6) + t129;
t509 = t127 / 0.2e1;
t510 = t107 / 0.2e1;
t511 = t106 / 0.2e1;
t516 = Ifges(6,1) * t510 + Ifges(6,4) * t511 + Ifges(6,5) * t508;
t42 = -qJD(6) * t138 + t106 * t355 - t107 * t351;
t517 = t42 / 0.2e1;
t41 = qJD(6) * t563 + t106 * t351 + t107 * t355;
t518 = t41 / 0.2e1;
t519 = Ifges(7,1) * t518 + Ifges(7,4) * t517 + Ifges(7,5) * t509;
t520 = Ifges(7,4) * t518 + Ifges(7,2) * t517 + Ifges(7,6) * t509;
t116 = t178 * t489 - t352 * t180;
t108 = -t344 * pkin(4) + qJD(5) - t116;
t89 = -pkin(5) * t402 + t108;
t575 = (Ifges(6,5) * t346 + Ifges(6,6) * t348) * t508 + (Ifges(6,1) * t346 + t467) * t510 + (Ifges(6,2) * t348 + t468) * t511 + t346 * t516 + t284 * t519 + t44 * t490 + (-Ifges(7,4) * t263 - Ifges(7,2) * t264) * t506 + (-Ifges(7,5) * t263 - Ifges(7,6) * t264) * t499 + (-Ifges(7,1) * t263 - Ifges(7,4) * t264) * t504 + Ifges(5,5) * t128 - Ifges(5,6) * t129 + (t572 * t70 + t573 * t71 + t385) * mrSges(6,3) - t382 * t520 + (Ifges(7,5) * t284 - Ifges(7,6) * t382) * t509 + (Ifges(7,4) * t284 - Ifges(7,2) * t382) * t517 + (Ifges(7,1) * t284 - Ifges(7,4) * t382) * t518 + t30 * (mrSges(7,1) * t382 + mrSges(7,2) * t284) + t47 * t392 + (t354 * t577 + t442 * t578) * g(2) + (t357 * t577 + t441 * t578) * g(1) + Ifges(5,3) * t341 - t583 * t514 + (t17 * t584 - t583 * t18 - t2 * t382 - t284 * t3) * mrSges(7,3) - t584 * t512 + (t583 * mrSges(7,1) - mrSges(7,2) * t584) * t89 + (-Ifges(7,4) * t570 - Ifges(7,2) * t571) * t507 + (-Ifges(7,1) * t570 - Ifges(7,4) * t571) * t505 + (-Ifges(7,5) * t570 - Ifges(7,6) * t571) * t500 - t49 * mrSges(5,2) + t50 * mrSges(5,1);
t538 = -t324 * t350 + t325 * t326;
t574 = m(7) * t538;
t310 = -t356 * mrSges(3,1) + t353 * mrSges(3,2);
t335 = sin(t345);
t337 = cos(t345);
t569 = -t337 * mrSges(4,1) + t335 * mrSges(4,2) + t310;
t568 = -t325 * mrSges(5,1) + (mrSges(5,2) - mrSges(7,3)) * t324;
t533 = g(1) * t357 + g(2) * t354;
t565 = t223 * mrSges(5,2) - t116 * mrSges(5,3);
t391 = mrSges(6,1) * t346 + mrSges(6,2) * t348;
t377 = t108 * t391;
t97 = t184 * Ifges(6,4) + Ifges(6,2) * t402 - Ifges(6,6) * t205;
t98 = t184 * Ifges(6,1) + Ifges(6,4) * t402 - Ifges(6,5) * t205;
t564 = t98 * t490 + t377 - t346 * t97 / 0.2e1;
t160 = pkin(4) * t361 - qJ(5) * t205;
t561 = t293 / 0.2e1;
t560 = t294 / 0.2e1;
t495 = -t361 / 0.2e1;
t477 = qJD(2) / 0.2e1;
t340 = t356 * pkin(2);
t328 = t340 + pkin(1);
t550 = t402 * Ifges(6,6);
t554 = t184 * Ifges(6,5);
t557 = t138 * Ifges(7,5) + Ifges(7,6) * t563 - Ifges(6,3) * t205 + t198 * Ifges(7,3) + t550 + t554;
t283 = t347 * t356 + t407;
t556 = Ifges(4,5) * t283;
t555 = Ifges(4,6) * t281;
t552 = t344 * Ifges(5,5);
t551 = t344 * Ifges(5,6);
t458 = qJDD(2) / 0.2e1;
t414 = t457 * pkin(2);
t327 = t414 + pkin(3);
t488 = pkin(2) * t347;
t256 = t352 * t327 + t489 * t488;
t252 = qJ(5) + t256;
t226 = (-pkin(9) - t252) * t346;
t339 = t348 * pkin(9);
t446 = t252 * t348;
t227 = t339 + t446;
t170 = t226 * t355 - t227 * t351;
t320 = t352 * t488;
t244 = -qJD(4) * t320 + t327 * t412;
t239 = qJD(5) + t244;
t219 = -t347 * t291 + t408;
t186 = t219 - t366;
t220 = t457 * t291 + t265;
t187 = t220 + t483;
t141 = t352 * t186 + t187 * t489;
t431 = qJD(1) * t353;
t331 = pkin(2) * t431;
t230 = -pkin(3) * t260 + t331;
t142 = t160 + t230;
t77 = -t141 * t346 + t348 * t142;
t55 = t576 + t77;
t78 = t348 * t141 + t346 * t142;
t61 = t78 - t579;
t549 = qJD(6) * t170 - t239 * t382 - t351 * t55 - t355 * t61;
t171 = t226 * t351 + t227 * t355;
t548 = -qJD(6) * t171 - t239 * t284 + t351 * t61 - t355 * t55;
t305 = t350 * t346;
t456 = qJ(5) * t348;
t306 = t339 + t456;
t224 = t305 * t355 - t306 * t351;
t85 = -t116 * t346 + t348 * t160;
t58 = t576 + t85;
t86 = t348 * t116 + t346 * t160;
t66 = t86 - t579;
t547 = -qJD(5) * t382 + qJD(6) * t224 - t351 * t58 - t355 * t66;
t225 = t305 * t351 + t306 * t355;
t546 = -qJD(5) * t284 - qJD(6) * t225 + t351 * t66 - t355 * t58;
t543 = mrSges(5,1) * t344 + mrSges(6,1) * t402 - mrSges(6,2) * t184 - mrSges(5,3) * t361;
t228 = t457 * t309 + t311 * t347;
t195 = -pkin(8) * t283 + t228;
t229 = t347 * t309 - t457 * t311;
t196 = pkin(8) * t281 + t229;
t542 = t489 * t195 - t352 * t196;
t541 = t244 - t141;
t146 = mrSges(6,2) * t205 + mrSges(6,3) * t402;
t147 = -mrSges(6,1) * t205 - mrSges(6,3) * t184;
t536 = t348 * t146 - t346 * t147;
t277 = -pkin(7) * t411 + t329;
t535 = t277 * t356 + t278 * t353;
t534 = -g(1) * t354 + g(2) * t357;
t315 = t324 * mrSges(6,3);
t532 = t325 * t582 - t315 + t568;
t67 = -t106 * mrSges(6,1) + t107 * mrSges(6,2);
t530 = m(6) * t47 + t67;
t529 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t342 = -pkin(8) + t349;
t528 = -m(3) * pkin(7) + m(4) * t349 + m(6) * t342 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - t391;
t527 = -m(5) * t116 + m(6) * t108 - t543;
t526 = -m(4) * t328 - (m(6) * pkin(4) - t392) * t325 - mrSges(2,1) - m(3) * pkin(1) + t568 + t569;
t524 = t223 * mrSges(5,1) + t70 * mrSges(6,1) + t17 * mrSges(7,1) - t71 * mrSges(6,2) - t18 * mrSges(7,2) - t117 * mrSges(5,3);
t386 = Ifges(6,5) * t348 - Ifges(6,6) * t346;
t388 = -Ifges(6,2) * t346 + t467;
t390 = Ifges(6,1) * t348 - t468;
t491 = -t344 / 0.2e1;
t497 = t205 / 0.2e1;
t501 = -t184 / 0.2e1;
t502 = -t402 / 0.2e1;
t523 = Ifges(5,1) * t495 + Ifges(5,5) * t491 + t386 * t497 + t388 * t502 + t390 * t501 - t565;
t496 = -t205 / 0.2e1;
t522 = Ifges(6,5) * t501 + Ifges(7,5) * t505 - Ifges(5,2) * t496 - Ifges(5,6) * t491 + Ifges(6,6) * t502 + Ifges(7,6) * t507 + Ifges(6,3) * t497 + Ifges(7,3) * t500 - t524;
t494 = t361 / 0.2e1;
t492 = -t260 / 0.2e1;
t487 = pkin(2) * t353;
t486 = pkin(4) * t324;
t485 = pkin(5) * t346;
t484 = pkin(7) * t356;
t409 = qJD(2) * t349;
t257 = t353 * t409 + t427;
t258 = t356 * t409 - t428;
t193 = -t257 * t347 + t457 * t258;
t262 = t281 * qJD(2);
t174 = -pkin(8) * t262 + t193;
t194 = t457 * t257 + t347 * t258;
t261 = t283 * qJD(2);
t175 = -pkin(8) * t261 + t194;
t79 = qJD(4) * t542 + t352 * t174 + t489 * t175;
t378 = t281 * t489 - t352 * t283;
t167 = qJD(4) * t378 - t352 * t261 + t262 * t489;
t216 = t352 * t281 + t283 * t489;
t168 = qJD(4) * t216 + t261 * t489 + t352 * t262;
t332 = pkin(2) * t429;
t231 = pkin(3) * t261 + t332;
t83 = pkin(4) * t168 - qJ(5) * t167 - qJD(5) * t216 + t231;
t32 = t346 * t83 + t348 * t79;
t472 = Ifges(3,4) * t353;
t471 = Ifges(3,4) * t356;
t470 = Ifges(4,4) * t260;
t469 = Ifges(5,4) * t361;
t461 = t213 * mrSges(4,3);
t73 = mrSges(6,1) * t129 - mrSges(6,3) * t107;
t459 = t346 * t73;
t454 = t167 * t346;
t453 = t167 * t348;
t448 = t216 * t346;
t447 = t216 * t348;
t313 = t324 * qJ(5);
t440 = t334 * t357;
t439 = t336 * t357;
t435 = t354 * t334;
t434 = t354 * t336;
t240 = -pkin(3) * t281 - t328;
t153 = -pkin(4) * t378 - qJ(5) * t216 + t240;
t155 = t352 * t195 + t196 * t489;
t88 = t346 * t153 + t348 * t155;
t433 = t325 * pkin(4) + t313;
t432 = pkin(3) * t337 + t340;
t423 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t431) * t484;
t422 = Ifges(7,5) * t41 + Ifges(7,6) * t42 + Ifges(7,3) * t127;
t12 = -t42 * mrSges(7,1) + t41 * mrSges(7,2);
t31 = -t346 * t79 + t348 * t83;
t405 = -t221 * mrSges(4,1) + t222 * mrSges(4,2);
t404 = t129 * mrSges(5,1) + t128 * mrSges(5,2);
t87 = t348 * t153 - t155 * t346;
t140 = -t489 * t186 + t187 * t352;
t399 = qJ(5) * t442 - t354 * t486;
t398 = qJ(5) * t441 - t357 * t486;
t255 = t327 * t489 - t320;
t396 = mrSges(3,1) * t353 + mrSges(3,2) * t356;
t393 = mrSges(5,1) * t324 + mrSges(5,2) * t325;
t389 = t356 * Ifges(3,2) + t472;
t387 = Ifges(3,5) * t356 - Ifges(3,6) * t353;
t384 = -t346 * t70 + t348 * t71;
t65 = -pkin(5) * t378 - pkin(9) * t447 + t87;
t74 = -pkin(9) * t448 + t88;
t26 = -t351 * t74 + t355 * t65;
t27 = t351 * t65 + t355 * t74;
t253 = -pkin(4) - t255;
t379 = pkin(1) * t396;
t376 = t353 * (Ifges(3,1) * t356 - t472);
t365 = mrSges(4,3) * t367;
t80 = qJD(4) * t155 - t489 * t174 + t352 * t175;
t330 = Ifges(3,4) * t430;
t308 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t430;
t298 = -qJD(1) * t328 + qJD(3);
t295 = -pkin(3) * t335 - t487;
t290 = pkin(1) + t432;
t276 = t357 * t295;
t275 = t354 * t295;
t271 = Ifges(3,1) * t431 + Ifges(3,5) * qJD(2) + t330;
t270 = Ifges(3,6) * qJD(2) + qJD(1) * t389;
t269 = t357 * t290;
t254 = Ifges(4,4) * t367;
t238 = qJD(2) * mrSges(4,1) + mrSges(4,3) * t260;
t237 = -qJD(2) * mrSges(4,2) + t365;
t236 = t253 - t478;
t235 = t325 * t439 + t435;
t234 = -t325 * t440 + t434;
t233 = -t325 * t434 + t440;
t232 = t325 * t435 + t439;
t211 = -mrSges(4,1) * t367 - t260 * mrSges(4,2);
t208 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t222;
t207 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t221;
t202 = -t260 * Ifges(4,1) + Ifges(4,5) * qJD(2) + t254;
t201 = Ifges(4,2) * t367 + Ifges(4,6) * qJD(2) - t470;
t197 = Ifges(5,4) * t205;
t188 = -mrSges(5,2) * t344 + mrSges(5,3) * t205;
t163 = t382 * t216;
t162 = t284 * t216;
t161 = -mrSges(5,1) * t205 + mrSges(5,2) * t361;
t157 = Ifges(5,1) * t361 + t197 + t552;
t156 = Ifges(5,2) * t205 + t469 + t551;
t114 = -mrSges(5,2) * t341 - mrSges(5,3) * t129;
t113 = mrSges(5,1) * t341 - mrSges(5,3) * t128;
t102 = pkin(5) * t448 - t542;
t93 = t140 + t580;
t92 = mrSges(7,1) * t198 - mrSges(7,3) * t138;
t91 = -mrSges(7,2) * t198 + mrSges(7,3) * t563;
t90 = t117 + t580;
t76 = -mrSges(7,1) * t563 + mrSges(7,2) * t138;
t72 = -mrSges(6,2) * t129 + mrSges(6,3) * t106;
t69 = -t167 * t284 + t216 * t263;
t68 = -t167 * t382 - t216 * t264;
t56 = pkin(5) * t454 + t80;
t29 = -mrSges(7,2) * t127 + mrSges(7,3) * t42;
t28 = mrSges(7,1) * t127 - mrSges(7,3) * t41;
t24 = -pkin(9) * t454 + t32;
t19 = pkin(5) * t168 - pkin(9) * t453 + t31;
t5 = -qJD(6) * t27 + t19 * t355 - t24 * t351;
t4 = qJD(6) * t26 + t19 * t351 + t24 * t355;
t1 = [(-t14 * t447 - t15 * t448 - t453 * t70 - t454 * t71) * mrSges(6,3) + (t181 * mrSges(5,2) - t50 * mrSges(5,3) + Ifges(5,1) * t128 - Ifges(5,4) * t129 + Ifges(5,5) * t341 + t386 * t508 + t388 * t511 + t390 * t510 + t391 * t47) * t216 + (Ifges(4,4) * t262 - Ifges(4,2) * t261) * t367 / 0.2e1 + (Ifges(4,1) * t262 - Ifges(4,4) * t261) * t492 + (Ifges(4,5) * t262 - Ifges(4,6) * t261) * t477 + t298 * (mrSges(4,1) * t261 + mrSges(4,2) * t262) + (-Ifges(7,5) * t163 - Ifges(7,6) * t162) * t509 + (-Ifges(7,4) * t163 - Ifges(7,2) * t162) * t517 + (-Ifges(7,1) * t163 - Ifges(7,4) * t162) * t518 + t30 * (mrSges(7,1) * t162 - mrSges(7,2) * t163) + (-mrSges(3,2) * t484 + t555 / 0.2e1 + t556 / 0.2e1 + Ifges(3,6) * t356 / 0.2e1) * qJDD(2) + (t555 + t556) * t458 + (-t162 * t2 + t163 * t3 - t17 * t68 + t18 * t69) * mrSges(7,3) + (Ifges(3,4) * t560 + Ifges(3,2) * t561 + Ifges(3,6) * t458 + t271 * t477) * t356 + (-pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t294) + Ifges(3,1) * t294 + Ifges(3,4) * t561 + 0.2e1 * Ifges(3,5) * t458) * t353 + (Ifges(7,5) * t68 + Ifges(7,6) * t69) * t499 + t262 * t202 / 0.2e1 + t68 * t512 + t69 * t514 + t447 * t516 - t163 * t519 - t162 * t520 + (-t165 * t283 + t166 * t281 - t212 * t262) * mrSges(4,3) - t261 * t201 / 0.2e1 + (t387 * t477 - t423) * qJD(2) + (Ifges(7,4) * t68 + Ifges(7,2) * t69) * t506 + (Ifges(5,1) * t494 + t386 * t496 + Ifges(5,4) * t497 + t157 / 0.2e1 + t402 * t388 / 0.2e1 + t184 * t390 / 0.2e1 + t552 / 0.2e1 + t564 + t565) * t167 + t211 * t332 + (Ifges(7,1) * t68 + Ifges(7,4) * t69) * t504 + (t293 * t484 + t535) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t535) + (t356 * (-Ifges(3,2) * t353 + t471) + t376) * t425 / 0.2e1 + t527 * t80 + t471 * t560 + t389 * t561 - t261 * t461 + Ifges(2,3) * qJDD(1) + t231 * t161 - t310 * t455 + t194 * t237 + t193 * t238 + t228 * t208 + t229 * t207 - t44 * t448 / 0.2e1 + t79 * t188 - t270 * t429 / 0.2e1 - t379 * t425 + t155 * t114 + m(4) * (t165 * t228 + t166 * t229 + t193 * t212 + t194 * t213 - t251 * t328 + t298 * t332) - t308 * t418 + t32 * t146 + t31 * t147 + (t557 / 0.2e1 + Ifges(7,6) * t506 - Ifges(5,4) * t494 + Ifges(6,3) * t496 - Ifges(5,2) * t497 + Ifges(7,3) * t499 + Ifges(7,5) * t504 - t156 / 0.2e1 + t550 / 0.2e1 + t554 / 0.2e1 - t551 / 0.2e1 + t524) * t168 + t102 * t12 + t4 * t91 + t5 * t92 + t87 * t73 + t88 * t72 + t89 * (-mrSges(7,1) * t69 + mrSges(7,2) * t68) + t56 * t76 + m(6) * (t14 * t87 + t15 * t88 + t31 * t70 + t32 * t71) + m(5) * (t117 * t79 + t155 * t49 + t181 * t240 + t223 * t231) + t240 * t404 - t328 * t405 + m(7) * (t102 * t30 + t17 * t5 + t18 * t4 + t2 * t27 + t26 * t3 + t56 * t89) - (t181 * mrSges(5,1) + t14 * mrSges(6,1) - t15 * mrSges(6,2) - t49 * mrSges(5,3) - Ifges(5,4) * t128 + Ifges(6,5) * t510 + Ifges(7,5) * t518 + Ifges(5,2) * t129 - Ifges(5,6) * t341 + Ifges(6,6) * t511 + Ifges(7,6) * t517 + Ifges(6,3) * t508 + Ifges(7,3) * t509 + t529) * t378 - (Ifges(6,5) * t107 + Ifges(6,6) * t106 + Ifges(6,3) * t129 + t422) * t378 / 0.2e1 + (-t233 * mrSges(7,1) - t232 * mrSges(7,2) + (-m(7) * (-t342 + t485) + m(5) * t342 + t528) * t357 + (-m(7) * (-t290 - t538) - m(6) * (-t290 - t313) + t315 + m(5) * t290 - t526) * t354) * g(1) + t251 * (-mrSges(4,1) * t281 + mrSges(4,2) * t283) + t221 * (Ifges(4,4) * t283 + Ifges(4,2) * t281) + t222 * (Ifges(4,1) * t283 + Ifges(4,4) * t281) - pkin(1) * (-mrSges(3,1) * t293 + mrSges(3,2) * t294) + (-m(6) * t269 - t235 * mrSges(7,1) - t234 * mrSges(7,2) + (-m(7) - m(5)) * (-t354 * t342 + t269) + (-m(7) * t485 + t528) * t354 + (-(m(6) * qJ(5) + mrSges(6,3)) * t324 - t574 + t526) * t357) * g(2) + t26 * t28 + t27 * t29 - (-m(5) * t50 - t113 + t530) * t542; (t270 / 0.2e1 + pkin(7) * t308) * t431 - t298 * (-t260 * mrSges(4,1) + mrSges(4,2) * t367) + (-t212 * t219 - t213 * t220 - t298 * t331 + (t165 * t457 + t166 * t347) * pkin(2)) * m(4) + (-m(4) * t340 - m(6) * (t432 + t433) - m(7) * (t538 + t432) - m(5) * t432 + t532 + t569) * g(3) + (-t108 * t140 - t70 * t77 - t71 * t78 + t239 * t384 + t252 * t385 + t253 * t47 - g(1) * (t276 + t398) - g(2) * (t275 + t399)) * m(6) + t548 * t92 + (-g(1) * t276 - g(2) * t275 + t548 * t17 + t170 * t3 + t171 * t2 + t549 * t18 + t236 * t30 - t89 * t93) * m(7) + t549 * t91 + t212 * t365 + (-Ifges(5,4) * t495 + t156 / 0.2e1 + t522) * t361 + (Ifges(5,4) * t496 - t377 - t157 / 0.2e1 + t523) * t205 + t207 * t488 + t201 * t492 + t255 * t113 + t253 * t67 + t533 * (m(4) * t487 - m(5) * t295 + mrSges(4,1) * t335 + mrSges(4,2) * t337 + t393 + t396) + t575 + t543 * t140 + t536 * t239 + t541 * t188 - (Ifges(4,2) * t260 + t202 + t254) * t367 / 0.2e1 - (-Ifges(3,2) * t431 + t271 + t330) * t430 / 0.2e1 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + t557 * t495 + t260 * (Ifges(4,1) * t367 + t470) / 0.2e1 - t260 * t461 - t252 * t459 + t236 * t12 - t220 * t237 - t219 * t238 - t230 * t161 + Ifges(4,6) * t221 + Ifges(4,5) * t222 - t387 * t425 / 0.2e1 + t170 * t28 + t171 * t29 + t165 * mrSges(4,1) - t166 * mrSges(4,2) - t78 * t146 - t77 * t147 - t211 * t331 - t93 * t76 - t277 * mrSges(3,2) - t278 * mrSges(3,1) + Ifges(3,6) * t293 + Ifges(3,5) * t294 - qJD(2) * (Ifges(4,5) * t367 + Ifges(4,6) * t260) / 0.2e1 + t208 * t414 + t72 * t446 + (t116 * t140 + t117 * t541 - t223 * t230 + t255 * t50) * m(5) + (t114 + t49 * m(5) + (m(7) * t89 + t527 + t76) * qJD(4)) * t256 - t98 * t572 / 0.2e1 + t97 * t573 / 0.2e1 + (t423 + (-t376 / 0.2e1 + t379) * qJD(1)) * qJD(1); -t237 * t367 - t260 * t238 - t382 * t28 + t284 * t29 + t346 * t72 + t348 * t73 + t404 + t405 - t583 * t92 - t584 * t91 + (-t188 - t536) * t205 + (-t76 + t543) * t361 + (-t17 * t583 - t18 * t584 + t2 * t284 - t3 * t382 - t361 * t89 + t534) * m(7) + (-t108 * t361 + t14 * t348 + t15 * t346 - t205 * t384 + t534) * m(6) + (t116 * t361 - t117 * t205 + t181 + t534) * m(5) + (-t212 * t260 - t213 * t367 + t251 + t534) * m(4); t533 * t393 - (-t523 + t564) * t205 + t546 * t92 + (t17 * t546 + t18 * t547 + t2 * t225 + t224 * t3 - t30 * t326 - t89 * t90) * m(7) + t547 * t91 + t522 * t361 + t156 * t494 + (-t469 + t557) * t495 + t575 + t543 * t117 + (-pkin(4) * t47 - g(1) * t398 - g(2) * t399 + qJ(5) * t385 + qJD(5) * t384 - t108 * t117 - t70 * t85 - t71 * t86) * m(6) + t536 * qJD(5) - qJ(5) * t459 + t224 * t28 + t225 * t29 - t116 * t188 - t86 * t146 - t85 * t147 - t90 * t76 - pkin(4) * t67 + (t157 + t197) * t496 + (-m(6) * t433 + t532 - t574) * g(3) - t326 * t12 + t72 * t456; t138 * t92 - t563 * t91 - t402 * t146 + t184 * t147 - m(6) * (-t184 * t70 + t402 * t71) + t12 + t530 + (g(3) * t325 - t324 * t533) * (m(6) + m(7)) + (t138 * t17 - t18 * t563 + t30) * m(7); -t89 * (mrSges(7,1) * t138 + mrSges(7,2) * t563) + (Ifges(7,1) * t563 - t466) * t505 + t63 * t504 + (Ifges(7,5) * t563 - Ifges(7,6) * t138) * t500 - t17 * t91 + t18 * t92 - g(1) * (mrSges(7,1) * t234 - mrSges(7,2) * t235) - g(2) * (-mrSges(7,1) * t232 + mrSges(7,2) * t233) - g(3) * (-mrSges(7,1) * t334 - mrSges(7,2) * t336) * t324 + (t138 * t18 + t17 * t563) * mrSges(7,3) + t422 + (-Ifges(7,2) * t138 + t131 + t64) * t507 + t529;];
tau  = t1;
