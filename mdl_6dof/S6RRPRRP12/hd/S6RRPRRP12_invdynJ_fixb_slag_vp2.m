% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP12_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP12_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP12_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP12_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP12_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP12_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:50:25
% EndTime: 2019-03-09 12:51:19
% DurationCPUTime: 34.32s
% Computational Cost: add. (10085->806), mult. (20392->1021), div. (0->0), fcn. (12726->10), ass. (0->383)
t597 = mrSges(6,2) - mrSges(7,3);
t587 = m(7) * qJ(6) - t597;
t589 = m(7) * pkin(5) + mrSges(7,1);
t531 = -mrSges(6,1) - t589;
t289 = sin(qJ(5));
t290 = sin(qJ(4));
t293 = cos(qJ(4));
t406 = qJD(5) * t289;
t409 = qJD(4) * t293;
t488 = cos(qJ(5));
t377 = t488 * qJD(5);
t545 = t488 * qJD(4) + t377;
t126 = -t289 * t409 - t290 * t545 - t293 * t406;
t291 = sin(qJ(2));
t324 = -t289 * t293 - t290 * t488;
t314 = t291 * t324;
t150 = qJD(1) * t314;
t595 = t126 + t150;
t410 = qJD(4) * t290;
t127 = -t289 * t410 - t290 * t406 + t293 * t545;
t273 = t291 * qJD(1);
t381 = t488 * t293;
t358 = t291 * t381;
t438 = t289 * t290;
t149 = qJD(1) * t358 - t273 * t438;
t594 = t149 + t127;
t294 = cos(qJ(2));
t414 = qJD(1) * t294;
t268 = pkin(7) * t414;
t207 = pkin(3) * t414 + t268;
t287 = qJD(2) * qJ(3);
t178 = t287 + t207;
t405 = t290 * qJD(2);
t319 = t293 * t414 + t405;
t129 = pkin(4) * t319 + t178;
t412 = qJD(2) * t293;
t318 = t290 * t414 - t412;
t121 = -t289 * t318 + t488 * t319;
t301 = -t289 * t319 - t318 * t488;
t47 = t121 * pkin(5) - qJ(6) * t301 + t129;
t118 = Ifges(6,4) * t121;
t253 = t273 + qJD(4);
t239 = qJD(5) + t253;
t461 = Ifges(7,5) * t121;
t566 = Ifges(7,4) + Ifges(6,5);
t569 = Ifges(6,1) + Ifges(7,1);
t558 = t239 * t566 + t301 * t569 - t118 + t461;
t599 = mrSges(6,2) * t129 - mrSges(7,3) * t47 + t558 / 0.2e1;
t572 = -m(6) - m(7);
t598 = mrSges(6,1) + mrSges(7,1);
t596 = mrSges(7,2) + mrSges(6,3);
t568 = -Ifges(4,4) + Ifges(3,5);
t567 = Ifges(6,4) - Ifges(7,5);
t565 = Ifges(4,5) - Ifges(3,6);
t564 = Ifges(7,2) + Ifges(6,3);
t563 = Ifges(6,6) - Ifges(7,6);
t348 = t294 * mrSges(4,2) - t291 * mrSges(4,3);
t353 = mrSges(3,1) * t294 - mrSges(3,2) * t291;
t593 = t348 - t353;
t490 = t239 / 0.2e1;
t499 = t301 / 0.2e1;
t502 = t121 / 0.2e1;
t503 = -t121 / 0.2e1;
t592 = Ifges(6,4) * t503 + Ifges(7,5) * t502 + t566 * t490 + t569 * t499 + t599;
t276 = t291 * qJ(3);
t368 = -pkin(1) - t276;
t507 = pkin(2) + pkin(8);
t140 = (-t294 * t507 + t368) * qJD(1);
t266 = pkin(7) * t273;
t206 = -pkin(3) * t273 - t266;
t536 = qJD(3) - t206;
t155 = -qJD(2) * t507 + t536;
t91 = t293 * t140 + t290 * t155;
t79 = -pkin(9) * t319 + t91;
t452 = t289 * t79;
t90 = -t140 * t290 + t293 * t155;
t78 = pkin(9) * t318 + t90;
t70 = pkin(4) * t253 + t78;
t27 = t488 * t70 - t452;
t23 = -t239 * pkin(5) + qJD(6) - t27;
t591 = mrSges(7,2) * t23 - mrSges(6,3) * t27;
t590 = qJD(3) + t266;
t288 = qJ(4) + qJ(5);
t274 = sin(t288);
t275 = cos(t288);
t350 = t290 * mrSges(5,1) + t293 * mrSges(5,2);
t588 = t531 * t274 + t275 * t587 - t350;
t404 = qJD(1) * qJD(2);
t210 = -t294 * qJDD(1) + t291 * t404;
t312 = t319 * qJD(4);
t113 = qJDD(2) * t293 + t210 * t290 - t312;
t114 = qJD(4) * t318 - qJDD(2) * t290 + t210 * t293;
t40 = -qJD(5) * t121 + t113 * t488 + t289 * t114;
t516 = t40 / 0.2e1;
t41 = qJD(5) * t301 + t289 * t113 - t114 * t488;
t514 = t41 / 0.2e1;
t211 = qJDD(1) * t291 + t294 * t404;
t198 = qJDD(4) + t211;
t189 = qJDD(5) + t198;
t495 = t189 / 0.2e1;
t391 = t488 * t79;
t33 = t289 * t78 + t391;
t584 = pkin(4) * t406 - t33;
t485 = pkin(4) * t293;
t262 = pkin(3) + t485;
t548 = pkin(4) * t409 + t262 * t273 + t590;
t200 = -t381 + t438;
t28 = t289 * t70 + t391;
t196 = t211 * pkin(7);
t369 = qJDD(3) + t196;
t116 = pkin(3) * t211 - qJDD(2) * t507 + t369;
t444 = qJDD(1) * pkin(1);
t310 = -qJ(3) * t211 - qJD(3) * t273 - t444;
t88 = t210 * t507 + t310;
t26 = -qJD(4) * t91 + t293 * t116 - t290 * t88;
t19 = pkin(4) * t198 - pkin(9) * t113 + t26;
t25 = t290 * t116 - t140 * t410 + t155 * t409 + t293 * t88;
t21 = pkin(9) * t114 + t25;
t5 = t289 * t19 + t488 * t21 + t70 * t377 - t406 * t79;
t6 = -qJD(5) * t28 + t19 * t488 - t289 * t21;
t583 = t200 * t6 - t27 * t595 - t594 * t28 + t324 * t5;
t2 = qJ(6) * t189 + qJD(6) * t239 + t5;
t24 = t239 * qJ(6) + t28;
t3 = -t189 * pkin(5) + qJDD(6) - t6;
t582 = t2 * t324 - t200 * t3 + t23 * t595 - t594 * t24;
t581 = -t294 * t596 + t593;
t491 = -t239 / 0.2e1;
t500 = -t301 / 0.2e1;
t580 = Ifges(6,4) * t502 + Ifges(7,5) * t503 + t491 * t566 + t500 * t569 - t599;
t71 = pkin(5) * t301 + qJ(6) * t121;
t313 = t319 * mrSges(5,3);
t135 = -t253 * mrSges(5,2) - t313;
t470 = mrSges(5,3) * t318;
t136 = mrSges(5,1) * t253 + t470;
t334 = t293 * t135 - t290 * t136;
t85 = mrSges(5,1) * t198 - mrSges(5,3) * t113;
t86 = -mrSges(5,2) * t198 + mrSges(5,3) * t114;
t579 = t334 * qJD(4) + t290 * t86 + t293 * t85;
t117 = Ifges(7,5) * t301;
t59 = Ifges(7,6) * t239 + Ifges(7,3) * t121 + t117;
t462 = Ifges(6,4) * t301;
t62 = -Ifges(6,2) * t121 + Ifges(6,6) * t239 + t462;
t521 = -t62 / 0.2e1 + t59 / 0.2e1 + mrSges(6,1) * t129 + mrSges(7,1) * t47;
t515 = -t41 / 0.2e1;
t195 = t210 * pkin(7);
t153 = -qJDD(2) * qJ(3) - qJD(2) * qJD(3) + t195;
t119 = -pkin(3) * t210 - t153;
t69 = -pkin(4) * t114 + t119;
t7 = pkin(5) * t41 - qJ(6) * t40 - qJD(6) * t301 + t69;
t578 = 0.2e1 * Ifges(7,3) * t514 + t69 * mrSges(6,1) + t7 * mrSges(7,1) - t40 * Ifges(6,4) / 0.2e1 - t189 * Ifges(6,6) / 0.2e1 + (t514 - t515) * Ifges(6,2) + (-t567 + Ifges(7,5)) * t516 + (-t563 + Ifges(7,6)) * t495;
t577 = Ifges(6,2) * t503 - Ifges(7,3) * t502 + t490 * t563 + t499 * t567;
t576 = -Ifges(6,2) * t502 + Ifges(7,3) * t503 - t491 * t563 - t500 * t567;
t459 = Ifges(4,6) * t294;
t338 = -t291 * Ifges(4,2) - t459;
t575 = t23 * mrSges(7,1) + t28 * mrSges(6,2) + Ifges(4,4) * qJD(2) / 0.2e1 + qJD(1) * t338 / 0.2e1 - t24 * mrSges(7,3) - t27 * mrSges(6,1);
t573 = -m(4) - m(5);
t571 = t319 / 0.2e1;
t562 = t189 * t566 + t40 * t569 - t41 * t567;
t29 = mrSges(6,1) * t189 - mrSges(6,3) * t40;
t30 = -t189 * mrSges(7,1) + t40 * mrSges(7,2);
t561 = t30 - t29;
t31 = -mrSges(6,2) * t189 - mrSges(6,3) * t41;
t32 = -mrSges(7,2) * t41 + mrSges(7,3) * t189;
t560 = t31 + t32;
t559 = pkin(5) * t594 - qJ(6) * t595 + qJD(6) * t200 + t548;
t472 = mrSges(7,2) * t121;
t94 = mrSges(7,3) * t239 - t472;
t469 = mrSges(6,3) * t121;
t95 = -mrSges(6,2) * t239 - t469;
t474 = t94 + t95;
t468 = mrSges(6,3) * t301;
t96 = mrSges(6,1) * t239 - t468;
t471 = mrSges(7,2) * t301;
t97 = -mrSges(7,1) * t239 + t471;
t473 = t96 - t97;
t281 = t294 * pkin(2);
t416 = t281 + t276;
t355 = pkin(8) * t294 + t416;
t191 = -pkin(1) - t355;
t506 = pkin(3) + pkin(7);
t230 = t506 * t291;
t202 = t290 * t230;
t125 = t293 * t191 + t202;
t426 = t293 * t294;
t107 = -pkin(9) * t426 + t125;
t203 = t293 * t230;
t367 = pkin(9) * t294 - t191;
t99 = pkin(4) * t291 + t290 * t367 + t203;
t553 = t488 * t107 + t289 * t99;
t163 = t324 * t294;
t128 = mrSges(5,1) * t319 - mrSges(5,2) * t318;
t389 = mrSges(4,1) * t414;
t222 = -qJD(2) * mrSges(4,3) - t389;
t552 = -t222 + t128;
t551 = qJD(2) * mrSges(3,2) - mrSges(3,3) * t414 + t222;
t390 = mrSges(4,1) * t273;
t550 = mrSges(3,3) * t273 + t390 + (-mrSges(3,1) + mrSges(4,2)) * qJD(2);
t549 = (t274 * t587 + t275 * t589) * t294;
t460 = Ifges(4,6) * t291;
t547 = t291 * (-Ifges(4,2) * t294 + t460) + t294 * (Ifges(4,3) * t291 - t459);
t546 = t291 * t565 + t294 * t568;
t544 = t189 * t564 + t40 * t566 - t41 * t563;
t295 = cos(qJ(1));
t292 = sin(qJ(1));
t431 = t291 * t292;
t158 = t274 * t295 + t275 * t431;
t159 = -t274 * t431 + t275 * t295;
t543 = -t158 * t598 - t159 * t597;
t430 = t291 * t295;
t156 = t274 * t292 - t275 * t430;
t157 = t274 * t430 + t275 * t292;
t542 = t156 * t598 + t157 * t597;
t541 = -t195 * t294 + t196 * t291;
t160 = -qJDD(2) * pkin(2) + t369;
t540 = -t153 * t294 + t160 * t291;
t538 = t25 * t290 + t26 * t293;
t537 = g(1) * t295 + g(2) * t292;
t535 = qJD(4) + qJD(5);
t465 = Ifges(5,4) * t318;
t103 = -Ifges(5,2) * t319 + Ifges(5,6) * t253 - t465;
t192 = Ifges(5,4) * t319;
t104 = -Ifges(5,1) * t318 + t253 * Ifges(5,5) - t192;
t337 = -t294 * Ifges(4,3) - t460;
t533 = Ifges(4,5) * qJD(2) + qJD(1) * t337 + t293 * t103 + t290 * t104;
t265 = Ifges(3,4) * t414;
t532 = Ifges(3,1) * t273 + Ifges(3,5) * qJD(2) - Ifges(5,5) * t318 - Ifges(5,6) * t319 + Ifges(5,3) * t253 - t121 * t563 + t239 * t564 + t301 * t566 + t265;
t435 = t290 * t291;
t332 = pkin(4) * t294 - pkin(9) * t435;
t445 = qJ(3) * t294;
t336 = pkin(8) * t291 - t445;
t413 = qJD(2) * t291;
t361 = pkin(2) * t413 - qJD(3) * t291;
t137 = qJD(2) * t336 + t361;
t411 = qJD(2) * t294;
t209 = t506 * t411;
t362 = -t137 * t290 + t293 * t209;
t48 = t332 * qJD(2) + (t293 * t367 - t202) * qJD(4) + t362;
t408 = qJD(4) * t294;
t379 = t290 * t408;
t316 = t291 * t412 + t379;
t66 = t293 * t137 - t191 * t410 + t290 * t209 + t230 * t409;
t51 = pkin(9) * t316 + t66;
t11 = -qJD(5) * t553 - t289 * t51 + t48 * t488;
t357 = -m(5) * t507 - mrSges(5,3);
t448 = t294 * mrSges(4,3);
t530 = -t448 + t588 * t294 + (-mrSges(4,2) - t357 + (m(4) - t572) * pkin(2) + t596) * t291;
t526 = -m(5) * pkin(3) - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t520 = mrSges(7,2) * t24 + mrSges(6,3) * t28 - t521;
t513 = -t113 * Ifges(5,4) / 0.2e1 - t114 * Ifges(5,2) / 0.2e1 - t198 * Ifges(5,6) / 0.2e1;
t505 = t113 / 0.2e1;
t504 = t114 / 0.2e1;
t494 = t198 / 0.2e1;
t493 = -t318 / 0.2e1;
t487 = pkin(4) * t318;
t486 = pkin(4) * t289;
t484 = pkin(7) * t291;
t279 = t294 * pkin(7);
t477 = -qJD(1) / 0.2e1;
t476 = -qJD(4) / 0.2e1;
t475 = pkin(9) + t507;
t267 = pkin(2) * t273;
t161 = qJD(1) * t336 + t267;
t110 = -t161 * t290 + t293 * t207;
t87 = qJD(1) * t332 + t110;
t111 = t293 * t161 + t290 * t207;
t93 = pkin(9) * t273 * t293 + t111;
t45 = t289 * t87 + t488 * t93;
t467 = Ifges(3,4) * t291;
t466 = Ifges(3,4) * t294;
t464 = Ifges(5,4) * t290;
t463 = Ifges(5,4) * t293;
t439 = t275 * t294;
t434 = t290 * t292;
t433 = t290 * t294;
t432 = t290 * t295;
t429 = t292 * t293;
t425 = t293 * t295;
t424 = t294 * t295;
t296 = -pkin(9) - pkin(8);
t423 = t294 * t296;
t422 = t295 * t296;
t259 = t290 * pkin(4) + qJ(3);
t256 = pkin(4) * t432;
t387 = t291 * t429;
t417 = pkin(4) * t387 + t256;
t231 = t294 * pkin(3) + t279;
t415 = t295 * pkin(1) + t292 * pkin(7);
t399 = pkin(4) * t434;
t254 = pkin(4) * t435;
t396 = t488 * pkin(4);
t392 = -t572 - t573;
t388 = t290 * t430;
t386 = t291 * t425;
t385 = Ifges(5,5) * t113 + Ifges(5,6) * t114 + Ifges(5,3) * t198;
t255 = pkin(4) * t426;
t173 = t255 + t231;
t380 = t291 * t405;
t214 = t475 * t293;
t371 = -t409 / 0.2e1;
t366 = -t404 / 0.2e1;
t166 = t211 * mrSges(4,1) + qJDD(2) * mrSges(4,2);
t364 = -t156 * pkin(5) + qJ(6) * t157;
t363 = t158 * pkin(5) - qJ(6) * t159;
t360 = pkin(4) * t377;
t359 = pkin(2) * t424 + qJ(3) * t430 + t415;
t356 = pkin(4) * t386 - t399;
t354 = t475 * t410;
t352 = mrSges(3,1) * t291 + mrSges(3,2) * t294;
t351 = mrSges(5,1) * t293 - mrSges(5,2) * t290;
t347 = Ifges(5,1) * t293 - t464;
t346 = Ifges(5,1) * t290 + t463;
t345 = t294 * Ifges(3,2) + t467;
t343 = -Ifges(5,2) * t290 + t463;
t342 = Ifges(5,2) * t293 + t464;
t340 = Ifges(5,5) * t293 - Ifges(5,6) * t290;
t339 = Ifges(5,5) * t290 + Ifges(5,6) * t293;
t335 = t90 * t290 - t91 * t293;
t212 = -qJD(2) * pkin(2) + t590;
t220 = -t268 - t287;
t333 = t212 * t294 + t220 * t291;
t331 = t368 - t281;
t44 = -t289 * t93 + t488 * t87;
t327 = pkin(1) * t352;
t55 = -t289 * t107 + t488 * t99;
t213 = t475 * t290;
t131 = -t213 * t488 - t289 * t214;
t325 = t289 * t213 - t214 * t488;
t10 = -t107 * t406 + t289 * t48 + t99 * t377 + t488 * t51;
t179 = t331 * qJD(1);
t323 = t179 * (-mrSges(4,2) * t291 - t448);
t322 = t291 * (Ifges(3,1) * t294 - t467);
t315 = -t293 * t408 + t380;
t309 = Ifges(5,5) * t294 + t291 * t346;
t308 = Ifges(5,6) * t294 + t291 * t342;
t307 = Ifges(5,3) * t294 + t291 * t339;
t305 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3) + t544;
t132 = -pkin(4) * t379 + (-pkin(7) - t262) * t413;
t302 = -qJD(4) * t335 + t538;
t282 = t295 * pkin(7);
t261 = -t396 - pkin(5);
t258 = qJ(6) + t486;
t250 = qJ(3) * t424;
t249 = t292 * t445;
t240 = t360 + qJD(6);
t216 = -pkin(1) - t416;
t208 = t506 * t413;
t205 = -qJ(3) * t414 + t267;
t204 = t348 * qJD(1);
t197 = qJD(4) * t214;
t190 = t351 * t294;
t183 = -t290 * t431 + t425;
t182 = t387 + t432;
t181 = t388 + t429;
t180 = t386 - t434;
t174 = Ifges(3,6) * qJD(2) + qJD(1) * t345;
t167 = -qJ(3) * t411 + t361;
t165 = mrSges(4,1) * t210 - qJDD(2) * mrSges(4,3);
t162 = t289 * t433 - t294 * t381;
t124 = -t191 * t290 + t203;
t115 = -pkin(5) * t324 + qJ(6) * t200 + t259;
t112 = pkin(2) * t210 + t310;
t82 = -pkin(5) * t162 - qJ(6) * t163 + t173;
t81 = qJD(2) * t358 - t163 * t535 - t289 * t380;
t80 = t200 * t294 * t535 - qJD(2) * t314;
t75 = qJD(5) * t131 - t289 * t197 - t354 * t488;
t74 = qJD(5) * t325 - t197 * t488 + t289 * t354;
t73 = mrSges(6,1) * t121 + mrSges(6,2) * t301;
t72 = mrSges(7,1) * t121 - mrSges(7,3) * t301;
t67 = -qJD(4) * t125 + t362;
t65 = -mrSges(5,1) * t114 + mrSges(5,2) * t113;
t58 = -t487 + t71;
t54 = t113 * Ifges(5,1) + t114 * Ifges(5,4) + t198 * Ifges(5,5);
t52 = -t291 * pkin(5) - t55;
t50 = qJ(6) * t291 + t553;
t43 = -pkin(5) * t414 - t44;
t42 = qJ(6) * t414 + t45;
t34 = t488 * t78 - t452;
t22 = -pkin(5) * t81 - qJ(6) * t80 - qJD(6) * t163 + t132;
t17 = mrSges(6,1) * t41 + mrSges(6,2) * t40;
t16 = mrSges(7,1) * t41 - mrSges(7,3) * t40;
t9 = -pkin(5) * t411 - t11;
t8 = qJ(6) * t411 + qJD(6) * t291 + t10;
t1 = [(-m(3) * t415 - t181 * mrSges(5,1) - t180 * mrSges(5,2) + t573 * t359 + t572 * (pkin(4) * t388 + t292 * t262 - t294 * t422 + t359) + t531 * t157 - t587 * t156 + (-m(5) * pkin(8) - mrSges(5,3) - t596) * t424 + (-mrSges(2,1) + t593) * t295 + t526 * t292) * g(2) + (-t183 * mrSges(5,1) + t182 * mrSges(5,2) + t572 * (t295 * t262 + t292 * t423 + t282) + t531 * t159 - t587 * t158 + (-m(3) + t573) * t282 + t526 * t295 + (m(3) * pkin(1) - m(4) * t331 - m(5) * t368 - t294 * t357 + mrSges(2,1) + t572 * (t331 - t254) - t581) * t292) * g(1) + m(4) * (t112 * t216 + t167 * t179 + (qJD(2) * t333 + t540) * pkin(7)) + m(5) * (t119 * t231 + t124 * t26 + t125 * t25 - t178 * t208 + t66 * t91 + t67 * t90) + (t163 * t69 - t291 * t5) * mrSges(6,2) + (-t210 * t279 + t211 * t484 + t541) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t541) + t546 * qJD(2) ^ 2 / 0.2e1 + t547 * t366 + (Ifges(7,5) * t163 + Ifges(7,6) * t291) * t514 + t103 * t379 / 0.2e1 + (t212 * mrSges(4,1) + Ifges(7,6) * t502 + t550 * pkin(7) + Ifges(6,6) * t503 + t532 / 0.2e1 - t91 * mrSges(5,2) + t90 * mrSges(5,1) + t564 * t490 + t566 * t499 - t575) * t411 + m(6) * (t10 * t28 + t11 * t27 + t129 * t132 + t173 * t69 + t5 * t553 + t55 * t6) + t553 * t31 + t353 * t444 + (-t163 * t7 + t2 * t291) * mrSges(7,3) + (t220 * mrSges(4,1) + t551 * pkin(7) + t533 / 0.2e1 - t174 / 0.2e1) * t413 + (-qJDD(2) * mrSges(3,2) - t165) * t279 + (t591 + t592) * t80 - t210 * t345 / 0.2e1 + t112 * t348 + t210 * t337 / 0.2e1 - t211 * t338 / 0.2e1 + m(7) * (t2 * t50 + t22 * t47 + t23 * t9 + t24 * t8 + t3 * t52 + t7 * t82) + t540 * mrSges(4,1) + qJD(2) * t323 + (-Ifges(3,4) * t210 + Ifges(3,5) * qJDD(2) + t385 + t544) * t291 / 0.2e1 + t211 * t466 / 0.2e1 + t178 * (-mrSges(5,1) * t316 + mrSges(5,2) * t315) + (t294 * (-Ifges(3,2) * t291 + t466) + t322) * t404 / 0.2e1 + (Ifges(6,4) * t163 + Ifges(6,6) * t291) * t515 + t294 * t104 * t371 + (t520 + t577) * t81 + (t2 * mrSges(7,2) + t5 * mrSges(6,3) - t578) * t162 + t211 * t291 * Ifges(3,1) + t91 * mrSges(5,3) * t316 + Ifges(2,3) * qJDD(1) - t90 * mrSges(5,3) * t315 - t294 * (Ifges(4,5) * qJDD(2) - Ifges(4,6) * t211 + Ifges(4,3) * t210) / 0.2e1 + t294 * (Ifges(3,4) * t211 - Ifges(3,2) * t210 + Ifges(3,6) * qJDD(2)) / 0.2e1 + (-qJDD(2) * mrSges(3,1) + t166) * t484 + t3 * (-t291 * mrSges(7,1) + mrSges(7,2) * t163) - t291 * (Ifges(4,4) * qJDD(2) - Ifges(4,2) * t211 + Ifges(4,6) * t210) / 0.2e1 + t6 * (t291 * mrSges(6,1) - mrSges(6,3) * t163) + (Ifges(5,6) * t291 - t294 * t342) * t504 + (Ifges(5,5) * t291 - t294 * t346) * t505 + t426 * t513 + (qJD(2) * t309 - t347 * t408) * t493 + (Ifges(5,3) * t291 - t294 * t339) * t494 + t231 * t65 + t216 * (-mrSges(4,2) * t210 - mrSges(4,3) * t211) - t208 * t128 - pkin(1) * (mrSges(3,1) * t210 + mrSges(3,2) * t211) + t167 * t204 + t119 * t190 - t54 * t433 / 0.2e1 + t26 * (mrSges(5,1) * t291 + mrSges(5,3) * t433) + t173 * t17 + t25 * (-mrSges(5,2) * t291 - mrSges(5,3) * t426) + t66 * t135 + t67 * t136 + t132 * t73 + t124 * t85 + t125 * t86 - t319 * (qJD(2) * t308 - t343 * t408) / 0.2e1 + t253 * (qJD(2) * t307 - t340 * t408) / 0.2e1 + t562 * t163 / 0.2e1 + (t163 * t566 + t291 * t564) * t495 + (t291 * t568 - t294 * t565) * qJDD(2) / 0.2e1 + (t163 * t569 + t291 * t566) * t516 + t50 * t32 + t52 * t30 + t55 * t29 + t22 * t72 + t82 * t16 - t327 * t404 + t8 * t94 + t10 * t95 + t11 * t96 + t9 * t97; (t119 * qJ(3) - t110 * t90 - t111 * t91 + t178 * t536) * m(5) + t537 * t352 + (-t409 * t91 + t410 * t90 - t538) * mrSges(5,3) - t550 * t268 - t551 * t266 + t552 * qJD(3) + t546 * t366 + t548 * t73 + (-t562 / 0.2e1 - mrSges(6,2) * t69 + mrSges(7,3) * t7 - Ifges(6,4) * t515 - Ifges(7,5) * t514 - t566 * t495 - t569 * t516) * t200 + t582 * mrSges(7,2) + t583 * mrSges(6,3) + (Ifges(6,6) * t502 + Ifges(7,6) * t503 + t564 * t491 + t566 * t500 + t575) * t414 + (t325 * t6 + t131 * t5 + t259 * t69 + (-t45 + t74) * t28 + (-t44 - t75) * t27 + t548 * t129) * m(6) + (t115 * t7 - t325 * t3 + t131 * t2 + t559 * t47 + (-t42 + t74) * t24 + (-t43 + t75) * t23) * m(7) - t561 * t325 - (t309 * t477 + t346 * t476) * t318 + (-m(5) * t302 - t579) * t507 + (t178 * t351 + t307 * t477 + t339 * t476) * t253 + (-pkin(2) * t160 - qJ(3) * t153 - qJD(3) * t220 - t179 * t205) * m(4) + t592 * t126 + (-m(4) * t416 - m(5) * t355 - t294 * mrSges(5,3) + t572 * (t254 + t416 - t423) + t588 * t291 + t581) * g(3) + t119 * t350 + t342 * t312 / 0.2e1 + (Ifges(4,1) + Ifges(3,3)) * qJDD(2) + (t521 - t576) * t149 - t533 * t273 / 0.2e1 - t580 * t150 + (-m(4) * t333 * pkin(7) - t90 * (mrSges(5,1) * t294 - mrSges(5,3) * t435) - t323 - t91 * (mrSges(5,3) * t291 * t293 - mrSges(5,2) * t294) + t308 * t571 + (t547 / 0.2e1 + t327 - t322 / 0.2e1) * qJD(1)) * qJD(1) - (-Ifges(3,2) * t273 + t265 + t532) * t414 / 0.2e1 + (-t165 + t65) * qJ(3) + (t521 - t577) * t127 - t578 * t324 + t293 * t54 / 0.2e1 + t343 * t504 + t347 * t505 + t290 * t513 + t340 * t494 + t259 * t17 - t206 * t128 - t205 * t204 + t195 * mrSges(3,2) - t196 * mrSges(3,1) - pkin(2) * t166 + t160 * mrSges(4,2) - t153 * mrSges(4,3) - t111 * t135 - t110 * t136 - t104 * t410 / 0.2e1 + t115 * t16 - t473 * t75 + t474 * t74 + t559 * t72 + t560 * t131 + t565 * t210 + t568 * t211 - t212 * t389 - t220 * t390 + (t572 * (t294 * t256 + t291 * t422 + t250) + t573 * t250 + t530 * t295) * g(1) + (t572 * (t294 * t399 + t296 * t431 + t249) + t573 * t249 + t530 * t292) * g(2) - t42 * t94 - t45 * t95 - t44 * t96 - t43 * t97 + t174 * t273 / 0.2e1 + t103 * t371; ((t204 + t334) * qJD(1) - t537 * t392) * t291 - t560 * t324 + (-t72 - t73 - t552) * qJD(2) + t166 + t561 * t200 + t392 * t294 * g(3) + t594 * t474 + t595 * t473 + (-qJD(2) * t47 - t582) * m(7) + (-qJD(2) * t129 - t583) * m(6) + (-qJD(2) * t178 - t273 * t335 + t302) * m(5) + (qJD(2) * t220 + t179 * t273 + t160) * m(4) + t579; (-t313 - t135) * t90 + (-t470 + t136) * t91 + ((t488 * t6 + t289 * t5 + (-t27 * t289 + t28 * t488) * qJD(5)) * pkin(4) + t129 * t487 + t27 * t33 - t28 * t34) * m(6) + (m(7) * t255 - (-m(6) * t485 - mrSges(6,1) * t275) * t294 + t190 + t549) * g(3) + (-mrSges(5,1) * t180 + mrSges(5,2) * t181 - m(7) * (t356 + t364) - m(6) * t356 + t542) * g(1) + (-mrSges(5,1) * t182 - mrSges(5,2) * t183 - m(7) * (t363 + t417) - m(6) * t417 + t543) * g(2) + (t2 * t258 + t261 * t3 - t47 * t58 + (t240 - t34) * t24 + t584 * t23) * m(7) - t584 * t473 + t318 * (-Ifges(5,1) * t319 + t465) / 0.2e1 - t253 * (-Ifges(5,5) * t319 + Ifges(5,6) * t318) / 0.2e1 - t178 * (-mrSges(5,1) * t318 - mrSges(5,2) * t319) + (Ifges(5,2) * t318 + t104 - t192) * t571 + (-t580 + t591) * t121 + t95 * t360 + (t520 + t576) * t301 + t305 + t73 * t487 + t103 * t493 + t31 * t486 + t258 * t32 + t261 * t30 + t240 * t94 - t474 * t34 + t385 - t25 * mrSges(5,2) + t26 * mrSges(5,1) + t29 * t396 - t58 * t72; t305 + (Ifges(7,3) * t301 - t461) * t503 + t62 * t499 + t24 * t471 + t23 * t472 - t129 * (mrSges(6,1) * t301 - mrSges(6,2) * t121) - t47 * (mrSges(7,1) * t301 + mrSges(7,3) * t121) - pkin(5) * t30 + qJ(6) * t32 - t71 * t72 + qJD(6) * t94 + (-t121 * t566 - t301 * t563) * t491 + (mrSges(6,1) * t439 + t549) * g(3) + t543 * g(2) + t542 * g(1) + (-pkin(5) * t3 - g(1) * t364 - g(2) * t363 + qJ(6) * t2 + qJD(6) * t24 - t47 * t71) * m(7) + (-Ifges(6,2) * t301 - t118 + t558) * t502 + (-t121 * t569 + t117 - t462 + t59) * t500 + (-m(7) * t23 + t468 + t473) * t28 + (-m(7) * t24 - t469 - t474) * t27; t301 * t72 - t239 * t94 + (-g(1) * t156 + g(2) * t158 - g(3) * t439 - t24 * t239 + t301 * t47 + t3) * m(7) + t30;];
tau  = t1;
