% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PPRRRP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_coriolismatJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:52:35
% EndTime: 2019-03-08 18:52:45
% DurationCPUTime: 6.44s
% Computational Cost: add. (10939->539), mult. (31034->780), div. (0->0), fcn. (34329->12), ass. (0->294)
t347 = cos(qJ(5));
t570 = Ifges(6,6) + Ifges(7,6);
t582 = t570 * t347;
t346 = sin(qJ(3));
t478 = sin(pkin(6));
t422 = sin(pkin(12)) * t478;
t505 = cos(qJ(3));
t343 = sin(pkin(7));
t405 = cos(pkin(12)) * t478;
t479 = cos(pkin(7));
t480 = cos(pkin(6));
t557 = t480 * t343 + t479 * t405;
t173 = t346 * t422 - t557 * t505;
t174 = t557 * t346 + t505 * t422;
t344 = sin(qJ(5));
t348 = cos(qJ(4));
t465 = t344 * t348;
t102 = t173 * t465 + t174 * t347;
t462 = t347 * t348;
t103 = -t173 * t462 + t174 * t344;
t497 = -qJ(6) - pkin(10);
t300 = t497 * t344;
t304 = t497 * t347;
t327 = -pkin(5) * t347 - pkin(4);
t345 = sin(qJ(4));
t475 = t173 * t345;
t476 = t103 * t347;
t477 = t102 * t344;
t572 = mrSges(7,3) + mrSges(6,3);
t581 = t572 * (t477 / 0.2e1 - t476 / 0.2e1) - m(6) * (pkin(4) * t475 + (t476 - t477) * pkin(10)) / 0.2e1 - m(7) * (t102 * t300 - t103 * t304 - t327 * t475) / 0.2e1;
t298 = -pkin(4) * t348 - t345 * pkin(10) - pkin(3);
t281 = t347 * t298;
t464 = t345 * t347;
t406 = -qJ(6) * t464 + t281;
t172 = (-pkin(9) * t344 - pkin(5)) * t348 + t406;
t216 = pkin(9) * t462 + t344 * t298;
t466 = t344 * t345;
t180 = -qJ(6) * t466 + t216;
t455 = pkin(9) * t465;
t215 = t281 - t455;
t331 = t344 * mrSges(7,1);
t332 = t347 * mrSges(7,2);
t558 = t332 + t331;
t263 = t558 * t348;
t305 = t344 * mrSges(6,1) + t347 * mrSges(6,2);
t264 = t305 * t348;
t501 = pkin(5) * t344;
t447 = pkin(9) + t501;
t294 = t447 * t348;
t537 = m(7) / 0.2e1;
t539 = m(6) / 0.2e1;
t580 = (pkin(9) * t348 + t215 * t344 - t216 * t347) * t539 + (t172 * t344 - t180 * t347 + t294) * t537 + t264 / 0.2e1 + t263 / 0.2e1;
t366 = -t343 * t405 + t479 * t480;
t127 = t174 * t345 - t348 * t366;
t128 = t174 * t348 + t345 * t366;
t314 = t345 * pkin(4) - pkin(10) * t348;
t220 = pkin(9) * t466 + t347 * t314;
t175 = t345 * pkin(5) - qJ(6) * t462 + t220;
t221 = -pkin(9) * t464 + t344 * t314;
t183 = -qJ(6) * t465 + t221;
t293 = t447 * t345;
t286 = mrSges(6,2) * t348 - mrSges(6,3) * t466;
t522 = -t286 / 0.2e1;
t482 = t348 * mrSges(7,2);
t285 = -mrSges(7,3) * t466 + t482;
t524 = -t285 / 0.2e1;
t430 = t524 + t522;
t290 = -mrSges(6,1) * t348 - mrSges(6,3) * t464;
t515 = t290 / 0.2e1;
t484 = t348 * mrSges(7,1);
t289 = -mrSges(7,3) * t464 - t484;
t517 = t289 / 0.2e1;
t352 = t430 * t347 + (t515 + t517) * t344 + t580;
t292 = t345 * mrSges(6,1) - mrSges(6,3) * t462;
t291 = t345 * mrSges(7,1) - mrSges(7,3) * t462;
t514 = t291 / 0.2e1;
t428 = t514 + t292 / 0.2e1;
t287 = -t345 * mrSges(7,2) - mrSges(7,3) * t465;
t288 = -t345 * mrSges(6,2) - mrSges(6,3) * t465;
t429 = t287 / 0.2e1 + t288 / 0.2e1;
t262 = t305 * t345;
t261 = t558 * t345;
t528 = t261 / 0.2e1;
t431 = t528 + t262 / 0.2e1;
t499 = pkin(9) * t345;
t83 = -t128 * t344 + t173 * t347;
t84 = t128 * t347 + t173 * t344;
t579 = t431 * t128 + t429 * t84 + t428 * t83 + (t128 * t499 + t220 * t83 + t221 * t84) * t539 + (t128 * t293 + t175 * t83 + t183 * t84) * t537 + t352 * t127;
t493 = Ifges(7,4) * t344;
t401 = Ifges(7,1) * t347 - t493;
t232 = -t348 * Ifges(7,5) + t345 * t401;
t494 = Ifges(6,4) * t344;
t402 = Ifges(6,1) * t347 - t494;
t234 = -t348 * Ifges(6,5) + t345 * t402;
t449 = Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1;
t408 = t449 * t348;
t578 = -t408 + t232 / 0.2e1 + t234 / 0.2e1;
t468 = t343 * t346;
t255 = t345 * t479 + t348 * t468;
t446 = t343 * t505;
t177 = -t344 * t255 - t347 * t446;
t178 = t255 * t347 - t344 * t446;
t483 = t348 * mrSges(5,2);
t306 = t345 * mrSges(5,1) + t483;
t390 = -t306 * t446 / 0.2e1;
t576 = (t220 * t177 + t221 * t178 + t255 * t499) * t539 + (t175 * t177 + t178 * t183 + t255 * t293) * t537 + t390 + t255 * t431 + t178 * t429 + t177 * t428;
t457 = m(6) / 0.4e1 + m(7) / 0.4e1;
t575 = 0.4e1 * t457;
t574 = mrSges(7,3) / 0.2e1;
t573 = t558 / 0.2e1;
t504 = m(7) * t293;
t571 = Ifges(6,5) + Ifges(7,5);
t569 = Ifges(6,3) + Ifges(7,3);
t535 = m(7) * pkin(5);
t568 = -mrSges(7,1) - t535;
t487 = t345 * mrSges(5,2);
t567 = -t348 * mrSges(5,1) - mrSges(4,1) + t487;
t562 = -t261 - t262;
t561 = t285 + t286;
t560 = t289 + t290;
t301 = -mrSges(7,1) * t347 + mrSges(7,2) * t344;
t496 = mrSges(6,1) * t347;
t404 = -mrSges(6,2) * t344 + t496;
t559 = t301 - t404;
t552 = -t220 * t344 + t221 * t347;
t335 = Ifges(7,4) * t347;
t551 = -t344 * Ifges(7,1) - t335;
t336 = Ifges(6,4) * t347;
t550 = -t344 * Ifges(6,1) - t336;
t549 = t285 / 0.2e1 + t286 / 0.2e1;
t547 = t332 / 0.2e1 + t331 / 0.2e1;
t452 = mrSges(7,1) / 0.2e1 + mrSges(6,1) / 0.2e1;
t458 = t535 / 0.2e1;
t388 = t458 + t452;
t533 = mrSges(7,2) / 0.2e1;
t534 = mrSges(6,2) / 0.2e1;
t451 = t533 + t534;
t409 = t451 * t347;
t510 = t305 / 0.2e1;
t546 = t409 + t510 + t573 + (t388 + t458) * t344;
t260 = t404 * t345;
t456 = pkin(5) * t464;
t259 = mrSges(7,1) * t464 - mrSges(7,2) * t466;
t531 = t259 / 0.2e1;
t545 = t456 * t537 + t260 / 0.2e1 + t531;
t179 = t406 - t455;
t461 = -t172 + t179;
t516 = -t290 / 0.2e1;
t518 = -t289 / 0.2e1;
t544 = t461 * t537 + t516 + t518;
t542 = 0.2e1 * m(7);
t541 = 2 * qJD(4);
t540 = m(5) / 0.2e1;
t532 = -t175 / 0.2e1;
t530 = -t260 / 0.2e1;
t512 = t300 / 0.2e1;
t511 = t301 / 0.2e1;
t509 = t306 / 0.2e1;
t508 = t344 / 0.2e1;
t507 = t345 / 0.2e1;
t506 = t347 / 0.2e1;
t503 = m(7) * t327;
t502 = m(7) * t345;
t341 = t348 ^ 2;
t500 = pkin(9) * t341;
t495 = mrSges(7,3) * t347;
t490 = t344 * t83;
t489 = t344 * t84;
t486 = t347 * t83;
t485 = t347 * t84;
t481 = -t404 - mrSges(5,1);
t474 = t177 * t347;
t473 = t178 * t344;
t444 = t348 * t505;
t213 = (-t344 * t444 + t346 * t347) * t343;
t472 = t213 * t344;
t214 = (t344 * t346 + t347 * t444) * t343;
t471 = t214 * t347;
t467 = t344 * t177;
t463 = t347 * t178;
t338 = t344 ^ 2;
t340 = t347 ^ 2;
t459 = t340 + t338;
t454 = t294 * t537;
t453 = t502 / 0.2e1;
t450 = t574 + mrSges(6,3) / 0.2e1;
t448 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t445 = t345 * t505;
t441 = -t475 / 0.2e1;
t436 = t466 / 0.2e1;
t434 = -t464 / 0.2e1;
t307 = t347 * Ifges(7,2) + t493;
t308 = t347 * Ifges(6,2) + t494;
t427 = -t307 / 0.2e1 - t308 / 0.2e1;
t426 = -t551 / 0.2e1 - t550 / 0.2e1;
t425 = m(7) * t461;
t424 = pkin(10) * t459;
t423 = qJD(5) * (-mrSges(6,2) - mrSges(7,2));
t420 = m(7) * t456;
t339 = t345 ^ 2;
t416 = t339 * t446;
t415 = t343 * t445;
t254 = t345 * t468 - t348 * t479;
t414 = t254 * t445;
t410 = t446 / 0.2e1;
t407 = qJD(4) * (t301 + t481);
t400 = -Ifges(6,2) * t344 + t336;
t399 = -Ifges(7,2) * t344 + t335;
t398 = t485 - t490;
t12 = (t102 * t83 + t103 * t84 - t127 * t475) * t575 + m(5) * (-t127 * t345 - t128 * t348 + t174) * t173;
t8 = ((-t254 * t345 - t255 * t348) * t173 + (t127 * t445 + t128 * t444 + t173 * t346 - t174 * t505) * t343) * t540 + (t537 + t539) * (t177 * t102 + t178 * t103 + t213 * t83 + t214 * t84 + (t127 * t446 - t173 * t254) * t345);
t397 = t12 * qJD(1) + t8 * qJD(2);
t386 = t128 - t398;
t15 = t386 * t127 * t575;
t393 = t463 - t467;
t383 = t255 - t393;
t9 = 0.2e1 * (t127 * t383 + t254 * t386) * t457;
t396 = t15 * qJD(1) + t9 * qJD(2);
t38 = t383 * t254 * t575;
t395 = t9 * qJD(1) + t38 * qJD(2);
t39 = m(5) * (t255 * t444 - t346 * t446 + t414) * t343 + (t177 * t213 + t178 * t214 + t343 * t414) * t575;
t394 = t8 * qJD(1) + t39 * qJD(2);
t392 = t300 * t344 + t304 * t347;
t391 = t345 * t410;
t389 = qJD(5) * (-mrSges(6,1) + t568);
t205 = -t259 - t420;
t269 = -m(7) * t501 - t558;
t387 = qJD(3) * t205 + qJD(4) * t269;
t385 = t425 / 0.2e1 + t518;
t349 = (-pkin(4) * t415 + (t471 - t472) * pkin(10)) * t539 + (t300 * t213 - t304 * t214 + t327 * t415) * t537 + t390 + t559 * t391 + t572 * (-t472 / 0.2e1 + t471 / 0.2e1);
t11 = -t349 + (t560 * t508 - t561 * t347 / 0.2e1 + t580) * t254 + t576;
t229 = Ifges(7,6) * t345 + t348 * t399;
t231 = Ifges(6,6) * t345 + t348 * t400;
t233 = Ifges(7,5) * t345 + t348 * t401;
t235 = Ifges(6,5) * t345 + t348 * t402;
t228 = -t348 * Ifges(7,6) + t345 * t399;
t230 = -t348 * Ifges(6,6) + t345 * t400;
t373 = -t228 / 0.2e1 - t230 / 0.2e1 + t448 * t348;
t17 = -pkin(3) * t306 + t172 * t291 + t175 * t289 + t180 * t287 + t183 * t285 + t215 * t292 + t216 * t288 + t220 * t290 + t221 * t286 + t294 * t261 + t293 * t263 + m(7) * (t172 * t175 + t180 * t183 + t293 * t294) + m(6) * (t215 * t220 + t216 * t221) + (Ifges(5,4) * t348 + pkin(9) * t262 + t373 * t344 + t578 * t347) * t348 + (-Ifges(5,4) * t345 + pkin(9) * t264 + (t233 / 0.2e1 + t235 / 0.2e1 + t449 * t345) * t347 + (-t229 / 0.2e1 - t231 / 0.2e1 - t448 * t345) * t344 + (m(6) * pkin(9) ^ 2 + Ifges(5,1) - Ifges(5,2) - t569) * t348) * t345;
t3 = (t509 - t483 / 0.2e1 + (-mrSges(5,1) / 0.2e1 - t404 / 0.2e1 + t511) * t345) * t173 + t579 + t581;
t378 = t3 * qJD(1) + t11 * qJD(2) + t17 * qJD(3);
t360 = t450 * t464 - t385 + t515;
t363 = t213 * t388 - t214 * t451;
t370 = -t450 * t466 + t430;
t376 = -t420 / 0.2e1 + t530 - t259 / 0.2e1;
t18 = t177 * t370 + t178 * t360 + t254 * t376 + t363;
t265 = t345 * t307;
t266 = t345 * t308;
t267 = t551 * t345;
t268 = t550 * t345;
t20 = t179 * t285 + t215 * t286 - t216 * t290 + t293 * t259 + (t425 - t289) * t180 + (pkin(9) * t260 + (t172 * mrSges(7,3) + t215 * mrSges(6,3) + t265 / 0.2e1 + t266 / 0.2e1 - t578) * t344 + (-t180 * mrSges(7,3) - t216 * mrSges(6,3) + t268 / 0.2e1 + t267 / 0.2e1 + (t261 + t504) * pkin(5) + t373) * t347) * t345;
t364 = t102 * t388 - t103 * t451;
t4 = t127 * t376 + t360 * t84 + t370 * t83 + t364;
t377 = -t4 * qJD(1) - t18 * qJD(2) + t20 * qJD(3);
t47 = (-t173 / 0.2e1 + t489 / 0.2e1 + t486 / 0.2e1) * t502;
t65 = (t344 * t285 - m(7) * (-t172 * t347 - t180 * t344) + t347 * t289) * t345;
t87 = (t410 + t474 / 0.2e1 + t473 / 0.2e1) * t502;
t375 = -qJD(1) * t47 - qJD(2) * t87 - qJD(3) * t65;
t350 = pkin(9) * t510 + (-t340 / 0.2e1 - t338 / 0.2e1) * pkin(10) * mrSges(6,3) + (mrSges(7,3) * t512 + t550 / 0.4e1 + t551 / 0.4e1 - t336 / 0.4e1 - t335 / 0.4e1 + (Ifges(6,2) / 0.4e1 + Ifges(7,2) / 0.4e1) * t344) * t344 + (t304 * t574 - t308 / 0.4e1 - t307 / 0.4e1 + (Ifges(6,1) / 0.4e1 + Ifges(7,1) / 0.4e1) * t347 + (-Ifges(6,4) / 0.4e1 - Ifges(7,4) / 0.4e1) * t344 + (t511 + t503 / 0.2e1) * pkin(5)) * t347;
t333 = Ifges(7,5) * t347;
t334 = Ifges(6,5) * t347;
t353 = (-t334 / 0.4e1 - t333 / 0.4e1) * t348 - t385 * t304 + pkin(4) * t530 + t293 * t573 + t285 * t512 + t327 * t531;
t356 = (t504 / 0.2e1 + t528) * pkin(5) - t228 / 0.4e1 - t230 / 0.4e1 + t267 / 0.4e1 + t268 / 0.4e1 + pkin(10) * t522;
t358 = t234 / 0.4e1 + t232 / 0.4e1 - t266 / 0.4e1 - t265 / 0.4e1 + pkin(10) * t516 + (t179 / 0.2e1 - t172 / 0.2e1) * mrSges(7,3);
t365 = mrSges(7,1) * t532 + t183 * t533 - t220 * mrSges(6,1) / 0.2e1 + t221 * t534;
t14 = (-t291 / 0.2e1 + m(7) * t532) * pkin(5) + ((0.3e1 / 0.4e1 * Ifges(6,6) + 0.3e1 / 0.4e1 * Ifges(7,6)) * t348 + t356) * t344 + (-t408 + t358) * t347 + (-Ifges(6,3) / 0.2e1 - Ifges(7,3) / 0.2e1 + t350) * t345 + t353 + t365;
t357 = -t305 / 0.2e1 + t409 + t452 * t344 - t547;
t28 = t357 * t254;
t49 = pkin(4) * t305 - t327 * t558 - t301 * t501 + (-t335 / 0.2e1 - t336 / 0.2e1 - t426) * t347 + (-pkin(5) * t503 + (Ifges(6,4) / 0.2e1 + Ifges(7,4) / 0.2e1) * t344 + (Ifges(7,2) / 0.2e1 - Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(7,1) / 0.2e1) * t347 - t427) * t344;
t6 = t357 * t127;
t369 = -t6 * qJD(1) - t28 * qJD(2) + t14 * qJD(3) - t49 * qJD(4);
t143 = -m(7) * t392 + mrSges(7,3) * t459;
t45 = (-t490 / 0.4e1 + t485 / 0.4e1 - t128 / 0.4e1) * t542;
t361 = m(7) * ((-t300 * t345 + t180) * t347 + (t304 * t345 - t172) * t344);
t59 = (t482 / 0.2e1 + t524) * t347 + (t484 / 0.2e1 + t517) * t344 + t454 - t361 / 0.2e1;
t85 = (t255 / 0.4e1 + t467 / 0.4e1 - t463 / 0.4e1) * t542;
t368 = qJD(1) * t45 - qJD(2) * t85 - qJD(3) * t59 + qJD(4) * t143;
t367 = qJD(4) * (-t572 * t459 + mrSges(5,2));
t299 = pkin(9) * t416;
t152 = t339 * pkin(9) * t173;
t88 = (-t473 - t474) * t453 + m(7) * t391;
t86 = (t255 + t393) * t537;
t60 = t285 * t506 + t361 / 0.2e1 + t344 * t518 + t454 + t547 * t348;
t48 = (-t486 - t489) * t453 + m(7) * t441;
t44 = (t128 + t398) * t537;
t29 = t546 * t254;
t19 = t549 * t177 + t544 * t178 + t545 * t254 + t363 + t572 * (t177 * t436 + t178 * t434);
t13 = t358 * t347 + t353 + t350 * t345 + pkin(5) * t514 + ((Ifges(6,6) / 0.4e1 + Ifges(7,6) / 0.4e1) * t348 + t356) * t344 + t175 * t458 - t365 + t569 * t507 - t570 * t465 / 0.2e1 + t571 * t462 / 0.2e1;
t10 = t254 * t352 + t349 + t576;
t7 = t546 * t127;
t5 = t545 * t127 + t544 * t84 + t549 * t83 + t364 + t572 * (t84 * t434 + t83 * t436);
t2 = mrSges(5,1) * t475 / 0.2e1 + t559 * t441 + (t483 / 0.2e1 + t509) * t173 + t579 - t581;
t1 = qJD(3) * t8 + qJD(4) * t9;
t16 = [qJD(3) * t12 + qJD(4) * t15, t1 (t560 * t102 + t561 * t103 + t567 * t174) * qJD(3) + t2 * qJD(4) + t5 * qJD(5) + t48 * qJD(6) + (mrSges(4,2) + t562 * t345 + (-t339 - t341) * mrSges(5,3)) * qJD(3) * t173 + 0.2e1 * ((t102 * t172 + t103 * t180 - t293 * t475) * t537 + (t102 * t215 + t103 * t216 - t152) * t539 + (-pkin(3) * t174 - t173 * t500 - t152) * t540) * qJD(3) + t397, t2 * qJD(3) + t7 * qJD(5) + t44 * qJD(6) + t128 * t407 + t127 * t367 + ((t127 * t392 + t128 * t327) * t537 + (-pkin(4) * t128 - t127 * t424) * t539) * t541 + t396, t5 * qJD(3) + t7 * qJD(4) + t389 * t84 + t423 * t83, qJD(3) * t48 + qJD(4) * t44; t1, qJD(3) * t39 + qJD(4) * t38 (-mrSges(4,2) * t446 + m(6) * t299 + m(5) * (t299 + (-pkin(3) * t346 + t500 * t505) * t343) + t567 * t468 + (-t562 + t504) * t415 + (m(6) * t216 + m(7) * t180 + t561) * t214 + (m(6) * t215 + m(7) * t172 + t560) * t213 + (t341 * t446 + t416) * mrSges(5,3)) * qJD(3) + t10 * qJD(4) + t19 * qJD(5) + t88 * qJD(6) + t394, t10 * qJD(3) + t29 * qJD(5) + t86 * qJD(6) + t255 * t407 + t254 * t367 + ((-pkin(4) * t255 - t254 * t424) * t539 + (t254 * t392 + t255 * t327) * t537) * t541 + t395, t19 * qJD(3) + t29 * qJD(4) + t177 * t423 + t178 * t389, qJD(3) * t88 + qJD(4) * t86; qJD(4) * t3 - qJD(5) * t4 - qJD(6) * t47 - t397, qJD(4) * t11 - qJD(5) * t18 - qJD(6) * t87 - t394, qJD(4) * t17 + qJD(5) * t20 - qJD(6) * t65, t13 * qJD(5) + t60 * qJD(6) + t378 + (-Ifges(5,6) * t345 + t327 * t263 + t300 * t291 + t294 * t301 - t304 * t287 - pkin(4) * t264 + m(7) * (t175 * t300 - t183 * t304 + t294 * t327) - t175 * t344 * mrSges(7,3) + pkin(9) * t487 + t183 * t495 + (m(6) * t552 + t347 * t288 - t344 * t292) * pkin(10) + (Ifges(5,5) + t426 * t347 + t427 * t344 + (-m(6) * pkin(4) + t481) * pkin(9)) * t348 + (t233 + t235) * t508 + (t571 * t344 + t582) * t507 + (t229 + t231) * t506 + t552 * mrSges(6,3)) * qJD(4), t13 * qJD(4) + t377 + (-mrSges(6,1) * t216 - mrSges(6,2) * t215 - mrSges(7,2) * t179 + (-t582 + (mrSges(7,3) * pkin(5) - t571) * t344) * t345 + t568 * t180) * qJD(5), qJD(4) * t60 + t375; -qJD(3) * t3 - qJD(5) * t6 + qJD(6) * t45 - t396, -qJD(3) * t11 - qJD(5) * t28 - qJD(6) * t85 - t395, qJD(5) * t14 - qJD(6) * t59 - t378, -qJD(5) * t49 + qJD(6) * t143, t369 + (-mrSges(7,2) * t300 - pkin(5) * t495 - pkin(10) * t496 + t333 + t334 + (mrSges(6,2) * pkin(10) - t570) * t344 - t568 * t304) * qJD(5), t368; qJD(3) * t4 + qJD(4) * t6, qJD(3) * t18 + qJD(4) * t28, -qJD(4) * t14 + qJD(6) * t205 - t377, qJD(6) * t269 - t369, 0, t387; qJD(3) * t47 - qJD(4) * t45, qJD(3) * t87 + qJD(4) * t85, qJD(4) * t59 - qJD(5) * t205 - t375, -qJD(5) * t269 - t368, -t387, 0;];
Cq  = t16;
