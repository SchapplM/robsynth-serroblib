% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:57:01
% EndTime: 2019-03-09 06:57:34
% DurationCPUTime: 17.95s
% Computational Cost: add. (15675->765), mult. (32343->1005), div. (0->0), fcn. (22619->18), ass. (0->364)
t572 = -mrSges(6,3) - mrSges(7,3);
t341 = cos(qJ(5));
t342 = cos(qJ(4));
t416 = qJD(4) * t342;
t404 = pkin(3) * t416;
t343 = cos(qJ(3));
t333 = sin(pkin(11));
t304 = pkin(1) * t333 + pkin(7);
t290 = t304 * qJD(1);
t386 = pkin(8) * qJD(1) + t290;
t338 = sin(qJ(3));
t420 = qJD(2) * t338;
t232 = t343 * t386 + t420;
t337 = sin(qJ(4));
t214 = t337 * t232;
t320 = t343 * qJD(2);
t231 = -t338 * t386 + t320;
t155 = t231 * t342 - t214;
t421 = qJD(1) * t343;
t422 = qJD(1) * t338;
t257 = -t337 * t422 + t342 * t421;
t271 = t337 * t343 + t338 * t342;
t258 = t271 * qJD(1);
t203 = pkin(4) * t258 - pkin(9) * t257;
t407 = pkin(3) * t422;
t183 = t203 + t407;
t336 = sin(qJ(5));
t88 = t341 * t155 + t336 * t183;
t575 = t341 * t404 - t88;
t87 = -t155 * t336 + t341 * t183;
t574 = -t336 * t404 - t87;
t415 = qJD(5) * t336;
t448 = t257 * t336;
t573 = t415 - t448;
t486 = pkin(3) * t337;
t310 = pkin(9) + t486;
t475 = -pkin(10) - t310;
t389 = qJD(5) * t475;
t411 = pkin(10) * t448;
t571 = t336 * t389 + t411 + t575;
t447 = t257 * t341;
t380 = t258 * pkin(5) - pkin(10) * t447;
t570 = t341 * t389 - t380 + t574;
t345 = -pkin(10) - pkin(9);
t397 = qJD(5) * t345;
t220 = qJD(3) * pkin(3) + t231;
t144 = t220 * t342 - t214;
t90 = t341 * t144 + t336 * t203;
t569 = t336 * t397 + t411 - t90;
t89 = -t144 * t336 + t341 * t203;
t568 = t341 * t397 - t380 - t89;
t331 = qJ(5) + qJ(6);
t321 = sin(t331);
t332 = qJ(3) + qJ(4);
t322 = sin(t332);
t472 = mrSges(6,2) * t336;
t567 = (-mrSges(7,2) * t321 - t472) * t322;
t340 = cos(qJ(6));
t335 = sin(qJ(6));
t329 = qJD(3) + qJD(4);
t221 = -t258 * t336 + t329 * t341;
t215 = t342 * t232;
t145 = t220 * t337 + t215;
t130 = pkin(9) * t329 + t145;
t334 = cos(pkin(11));
t305 = -pkin(1) * t334 - pkin(2);
t326 = t343 * pkin(3);
t286 = t305 - t326;
t259 = t286 * qJD(1);
t168 = -pkin(4) * t257 - pkin(9) * t258 + t259;
t85 = t130 * t341 + t168 * t336;
t66 = pkin(10) * t221 + t85;
t457 = t335 * t66;
t252 = qJD(5) - t257;
t222 = t258 * t341 + t329 * t336;
t84 = -t130 * t336 + t341 * t168;
t65 = -pkin(10) * t222 + t84;
t60 = pkin(5) * t252 + t65;
t19 = t340 * t60 - t457;
t455 = t340 * t66;
t20 = t335 * t60 + t455;
t549 = Ifges(5,6) * t329;
t566 = -t259 * mrSges(5,1) - t84 * mrSges(6,1) - t19 * mrSges(7,1) + t85 * mrSges(6,2) + t20 * mrSges(7,2) + t549 / 0.2e1;
t546 = t329 * Ifges(5,5);
t565 = -t259 * mrSges(5,2) + t144 * mrSges(5,3) - t546 / 0.2e1;
t413 = qJD(1) * qJD(3);
t275 = qJDD(1) * t343 - t338 * t413;
t276 = qJDD(1) * t338 + t343 * t413;
t170 = -qJD(4) * t258 + t275 * t342 - t337 * t276;
t167 = qJDD(5) - t170;
t164 = qJDD(6) + t167;
t502 = t164 / 0.2e1;
t269 = t337 * t338 - t342 * t343;
t355 = t269 * qJD(4);
t169 = -qJD(1) * t355 + t275 * t337 + t276 * t342;
t328 = qJDD(3) + qJDD(4);
t103 = qJD(5) * t221 + t169 * t341 + t328 * t336;
t104 = -qJD(5) * t222 - t169 * t336 + t328 * t341;
t139 = t221 * t335 + t222 * t340;
t37 = -qJD(6) * t139 - t103 * t335 + t104 * t340;
t514 = t37 / 0.2e1;
t384 = t340 * t221 - t222 * t335;
t36 = qJD(6) * t384 + t103 * t340 + t104 * t335;
t515 = t36 / 0.2e1;
t516 = Ifges(7,1) * t515 + Ifges(7,4) * t514 + Ifges(7,5) * t502;
t517 = Ifges(7,4) * t515 + Ifges(7,2) * t514 + Ifges(7,6) * t502;
t244 = qJD(6) + t252;
t547 = Ifges(6,3) * t252;
t548 = Ifges(6,6) * t221;
t564 = Ifges(6,5) * t222 + Ifges(7,5) * t139 + Ifges(7,6) * t384 + Ifges(7,3) * t244 + t547 + t548;
t458 = t258 * mrSges(5,3);
t540 = mrSges(5,1) * t329 + mrSges(6,1) * t221 - mrSges(6,2) * t222 - t458;
t323 = cos(t331);
t474 = mrSges(6,1) * t341;
t563 = (mrSges(7,1) * t323 + t474) * t322;
t562 = t573 * pkin(5);
t154 = t231 * t337 + t215;
t417 = qJD(4) * t337;
t561 = pkin(3) * t417 - t154;
t414 = qJD(5) * t341;
t243 = t290 * t343 + t420;
t288 = t304 * qJDD(1);
t188 = -qJD(3) * t243 + t343 * qJDD(2) - t288 * t338;
t153 = qJDD(3) * pkin(3) - pkin(8) * t276 + t188;
t419 = qJD(3) * t338;
t187 = qJD(3) * t320 + t338 * qJDD(2) + t343 * t288 - t290 * t419;
t159 = pkin(8) * t275 + t187;
t54 = t337 * t153 + t342 * t159 + t220 * t416 - t232 * t417;
t50 = pkin(9) * t328 + t54;
t289 = t305 * qJDD(1);
t233 = -pkin(3) * t275 + t289;
t81 = -pkin(4) * t170 - pkin(9) * t169 + t233;
t15 = -t130 * t415 + t168 * t414 + t336 * t81 + t341 * t50;
t16 = -qJD(5) * t85 - t336 * t50 + t341 * t81;
t560 = t15 * t341 - t16 * t336;
t324 = cos(t332);
t559 = t324 * mrSges(5,1) + (-mrSges(5,2) - t572) * t322;
t129 = -pkin(4) * t329 - t144;
t99 = -pkin(5) * t221 + t129;
t558 = -mrSges(7,1) * t99 + mrSges(7,3) * t20;
t557 = mrSges(7,2) * t99 - t19 * mrSges(7,3);
t535 = t324 * pkin(4) + t322 * pkin(9);
t311 = pkin(5) * t341 + pkin(4);
t537 = t324 * t311 - t322 * t345;
t555 = -m(6) * t535 - m(7) * t537;
t518 = m(7) * pkin(5);
t553 = -m(4) - m(3);
t508 = t103 / 0.2e1;
t507 = t104 / 0.2e1;
t501 = t167 / 0.2e1;
t552 = t275 / 0.2e1;
t263 = t475 * t336;
t325 = t341 * pkin(10);
t443 = t310 * t341;
t264 = t325 + t443;
t200 = t263 * t335 + t264 * t340;
t551 = -qJD(6) * t200 - t335 * t571 + t340 * t570;
t199 = t263 * t340 - t264 * t335;
t550 = qJD(6) * t199 + t335 * t570 + t340 * t571;
t545 = mrSges(6,1) + t518;
t295 = t345 * t336;
t296 = pkin(9) * t341 + t325;
t228 = t295 * t335 + t296 * t340;
t544 = -qJD(6) * t228 - t335 * t569 + t340 * t568;
t227 = t295 * t340 - t296 * t335;
t543 = qJD(6) * t227 + t335 * t568 + t340 * t569;
t149 = mrSges(5,1) * t328 - mrSges(5,3) * t169;
t52 = -mrSges(6,1) * t104 + mrSges(6,2) * t103;
t542 = t52 - t149;
t374 = mrSges(6,1) * t336 + mrSges(6,2) * t341;
t541 = t129 * t374;
t363 = t335 * t336 - t340 * t341;
t194 = t363 * t271;
t476 = pkin(8) + t304;
t261 = t476 * t338;
t262 = t476 * t343;
t198 = -t261 * t337 + t262 * t342;
t189 = t341 * t198;
t190 = pkin(4) * t269 - pkin(9) * t271 + t286;
t113 = t336 * t190 + t189;
t539 = -t342 * t261 - t262 * t337;
t538 = t561 + t562;
t364 = -t311 * t322 - t324 * t345;
t483 = pkin(4) * t322;
t485 = pkin(3) * t338;
t534 = -m(7) * (t364 - t485) - m(6) * (-t483 - t485) + t563;
t533 = -m(7) * t364 + t563;
t532 = -t145 + t562;
t330 = qJ(1) + pkin(11);
t317 = cos(t330);
t440 = t317 * t324;
t530 = t567 * t317 + t440 * t572;
t316 = sin(t330);
t442 = t316 * t324;
t529 = t567 * t316 + t442 * t572;
t528 = t187 * t343 - t188 * t338;
t527 = g(1) * t317 + g(2) * t316;
t526 = -m(7) - m(6) - m(5);
t525 = qJD(5) + qJD(6);
t434 = t324 * t341;
t435 = t324 * t336;
t436 = t323 * t324;
t438 = t321 * t324;
t524 = -mrSges(6,1) * t434 - mrSges(7,1) * t436 + mrSges(6,2) * t435 + mrSges(7,2) * t438 - t559;
t12 = pkin(10) * t104 + t15;
t7 = pkin(5) * t167 - pkin(10) * t103 + t16;
t3 = qJD(6) * t19 + t12 * t340 + t335 * t7;
t4 = -qJD(6) * t20 - t12 * t335 + t340 * t7;
t523 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t522 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t521 = -m(4) * pkin(7) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t294 = -mrSges(4,1) * t343 + mrSges(4,2) * t338;
t520 = m(4) * pkin(2) + mrSges(3,1) - t294 + t559;
t471 = mrSges(6,3) * t221;
t171 = -mrSges(6,2) * t252 + t471;
t470 = mrSges(6,3) * t222;
t172 = mrSges(6,1) * t252 - t470;
t350 = (-t336 * t85 - t341 * t84) * qJD(5) + t560;
t70 = mrSges(6,1) * t167 - mrSges(6,3) * t103;
t456 = t336 * t70;
t519 = m(6) * t350 - t171 * t415 - t172 * t414 - t456;
t513 = Ifges(6,1) * t508 + Ifges(6,4) * t507 + Ifges(6,5) * t501;
t463 = Ifges(7,4) * t139;
t68 = Ifges(7,2) * t384 + Ifges(7,6) * t244 + t463;
t512 = -t68 / 0.2e1;
t511 = t68 / 0.2e1;
t134 = Ifges(7,4) * t384;
t69 = Ifges(7,1) * t139 + Ifges(7,5) * t244 + t134;
t510 = -t69 / 0.2e1;
t509 = t69 / 0.2e1;
t506 = -t384 / 0.2e1;
t505 = t384 / 0.2e1;
t504 = -t139 / 0.2e1;
t503 = t139 / 0.2e1;
t499 = -t221 / 0.2e1;
t498 = -t222 / 0.2e1;
t497 = t222 / 0.2e1;
t496 = -t244 / 0.2e1;
t495 = t244 / 0.2e1;
t494 = -t252 / 0.2e1;
t492 = t257 / 0.2e1;
t490 = t258 / 0.2e1;
t339 = sin(qJ(1));
t487 = pkin(1) * t339;
t484 = pkin(3) * t342;
t482 = pkin(5) * t222;
t479 = g(3) * t322;
t344 = cos(qJ(1));
t327 = t344 * pkin(1);
t469 = Ifges(4,4) * t338;
t468 = Ifges(4,4) * t343;
t467 = Ifges(5,4) * t258;
t466 = Ifges(6,4) * t222;
t465 = Ifges(6,4) * t336;
t464 = Ifges(6,4) * t341;
t462 = pkin(5) * qJD(6);
t71 = -mrSges(6,2) * t167 + mrSges(6,3) * t104;
t454 = t341 * t71;
t209 = -qJD(3) * t269 - t355;
t450 = t209 * t341;
t445 = t271 * t336;
t444 = t271 * t341;
t441 = t316 * t336;
t439 = t317 * t336;
t123 = Ifges(6,2) * t221 + Ifges(6,6) * t252 + t466;
t433 = t336 * t123;
t213 = Ifges(6,4) * t221;
t124 = Ifges(6,1) * t222 + Ifges(6,5) * t252 + t213;
t430 = t341 * t124;
t223 = t316 * t438 + t317 * t323;
t224 = -t316 * t436 + t317 * t321;
t429 = -t223 * mrSges(7,1) + t224 * mrSges(7,2);
t225 = t316 * t323 - t317 * t438;
t226 = t316 * t321 + t317 * t436;
t428 = t225 * mrSges(7,1) - t226 * mrSges(7,2);
t418 = qJD(3) * t343;
t410 = Ifges(7,5) * t36 + Ifges(7,6) * t37 + Ifges(7,3) * t164;
t406 = pkin(3) * t419;
t400 = mrSges(4,3) * t422;
t399 = mrSges(4,3) * t421;
t398 = Ifges(6,5) * t103 + Ifges(6,6) * t104 + Ifges(6,3) * t167;
t396 = t271 * t414;
t393 = t430 / 0.2e1;
t391 = -t415 / 0.2e1;
t390 = qJD(3) * t476;
t248 = t338 * t390;
t249 = t343 * t390;
t116 = qJD(4) * t539 - t248 * t342 - t249 * t337;
t210 = t329 * t271;
t128 = pkin(4) * t210 - pkin(9) * t209 + t406;
t385 = -t116 * t336 + t341 * t128;
t112 = t341 * t190 - t198 * t336;
t377 = mrSges(4,1) * t338 + mrSges(4,2) * t343;
t375 = mrSges(5,1) * t322 + mrSges(5,2) * t324;
t373 = -mrSges(7,1) * t321 - mrSges(7,2) * t323;
t372 = Ifges(6,1) * t341 - t465;
t371 = t343 * Ifges(4,2) + t469;
t370 = -Ifges(6,2) * t336 + t464;
t369 = Ifges(4,5) * t343 - Ifges(4,6) * t338;
t368 = Ifges(6,5) * t341 - Ifges(6,6) * t336;
t86 = pkin(5) * t269 - pkin(10) * t444 + t112;
t93 = -pkin(10) * t445 + t113;
t44 = -t335 * t93 + t340 * t86;
t45 = t335 * t86 + t340 * t93;
t366 = -t336 * t84 + t341 * t85;
t270 = t335 * t341 + t336 * t340;
t55 = t153 * t342 - t337 * t159 - t220 * t417 - t232 * t416;
t360 = t410 + t523;
t238 = t316 * t341 - t317 * t435;
t236 = t316 * t435 + t317 * t341;
t359 = t209 * t336 + t396;
t358 = t271 * t415 - t450;
t357 = t305 * qJD(1) * t377;
t356 = t338 * (Ifges(4,1) * t343 - t469);
t38 = t341 * t116 + t336 * t128 + t190 * t414 - t198 * t415;
t51 = -pkin(4) * t328 - t55;
t117 = qJD(4) * t198 - t248 * t337 + t342 * t249;
t208 = t525 * t270;
t175 = t270 * t257;
t176 = t363 * t257;
t184 = Ifges(5,2) * t257 + t467 + t549;
t250 = Ifges(5,4) * t257;
t185 = t258 * Ifges(5,1) + t250 + t546;
t29 = -pkin(5) * t104 + t51;
t42 = Ifges(6,4) * t103 + Ifges(6,2) * t104 + Ifges(6,6) * t167;
t348 = (-Ifges(7,4) * t176 - Ifges(7,2) * t175) * t506 + (-(Ifges(7,1) * t503 + Ifges(7,4) * t505 + Ifges(7,5) * t495 + t509 + t557) * t525 + t29 * mrSges(7,1) - 0.2e1 * t517 - t3 * mrSges(7,3)) * t363 + (mrSges(7,2) * t29 - mrSges(7,3) * t4 + 0.2e1 * t516) * t270 + (-Ifges(7,1) * t176 - Ifges(7,4) * t175) * t504 + (t368 * t494 + t370 * t499 + t372 * t498 - t541 + t565) * t257 + (Ifges(6,5) * t498 + Ifges(7,5) * t504 + Ifges(6,6) * t499 + Ifges(7,6) * t506 + Ifges(6,3) * t494 + Ifges(7,3) * t496 + t566) * t258 + (t541 + t393) * qJD(5) + t51 * (t472 - t474) + t123 * t391 + (Ifges(6,2) * t341 + t465) * t507 + (Ifges(6,1) * t336 + t464) * t508 - t176 * t510 - t175 * t512 + t336 * t513 + t184 * t490 + t433 * t492 + (Ifges(6,5) * t336 + Ifges(6,6) * t341) * t501 + (-Ifges(7,5) * t176 - Ifges(7,6) * t175) * t496 - t99 * (mrSges(7,1) * t175 - mrSges(7,2) * t176) + t145 * t458 + t341 * t42 / 0.2e1 + Ifges(5,3) * t328 - t54 * mrSges(5,2) + t55 * mrSges(5,1) + (-t573 * t85 + (-t414 + t447) * t84 + t560) * mrSges(6,3) - (Ifges(7,4) * t503 + Ifges(7,2) * t505 + Ifges(7,6) * t495 + t511 + t558) * t208 + (t221 * t370 + t222 * t372 + t252 * t368) * qJD(5) / 0.2e1 - (-Ifges(5,2) * t258 + t185 + t250 + t430) * t257 / 0.2e1 - (Ifges(5,1) * t257 - t467 + t564) * t258 / 0.2e1 + Ifges(5,5) * t169 + Ifges(5,6) * t170 + (t175 * t20 - t176 * t19) * mrSges(7,3);
t346 = -pkin(8) - pkin(7);
t314 = Ifges(4,4) * t421;
t313 = t326 + pkin(2);
t312 = -pkin(4) - t484;
t293 = -qJD(3) * mrSges(4,2) + t399;
t291 = qJD(3) * mrSges(4,1) - t400;
t287 = -t311 - t484;
t283 = pkin(9) * t440;
t282 = pkin(9) * t442;
t256 = Ifges(4,1) * t422 + Ifges(4,5) * qJD(3) + t314;
t255 = Ifges(4,6) * qJD(3) + qJD(1) * t371;
t247 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t276;
t246 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t275;
t242 = -t290 * t338 + t320;
t239 = t317 * t434 + t441;
t237 = -t316 * t434 + t439;
t234 = -mrSges(5,2) * t329 + mrSges(5,3) * t257;
t202 = -mrSges(5,1) * t257 + mrSges(5,2) * t258;
t193 = t270 * t271;
t158 = pkin(5) * t445 - t539;
t150 = -mrSges(5,2) * t328 + mrSges(5,3) * t170;
t111 = mrSges(7,1) * t244 - mrSges(7,3) * t139;
t110 = -mrSges(7,2) * t244 + mrSges(7,3) * t384;
t80 = -mrSges(7,1) * t384 + mrSges(7,2) * t139;
t74 = pkin(5) * t359 + t117;
t64 = t194 * t525 - t270 * t209;
t63 = -t208 * t271 - t209 * t363;
t39 = -qJD(5) * t113 + t385;
t26 = -pkin(10) * t359 + t38;
t25 = -mrSges(7,2) * t164 + mrSges(7,3) * t37;
t24 = mrSges(7,1) * t164 - mrSges(7,3) * t36;
t22 = t340 * t65 - t457;
t21 = -t335 * t65 - t455;
t18 = -pkin(10) * t450 + pkin(5) * t210 + (-t189 + (pkin(10) * t271 - t190) * t336) * qJD(5) + t385;
t13 = -mrSges(7,1) * t37 + mrSges(7,2) * t36;
t6 = -qJD(6) * t45 + t18 * t340 - t26 * t335;
t5 = qJD(6) * t44 + t18 * t335 + t26 * t340;
t1 = [(-t433 / 0.2e1 + t393 + Ifges(5,1) * t490 + Ifges(5,4) * t492 + t185 / 0.2e1 - t565) * t209 + (Ifges(7,5) * t503 + Ifges(7,6) * t505 - Ifges(5,4) * t490 - Ifges(5,2) * t492 + Ifges(7,3) * t495 + Ifges(6,5) * t497 - t184 / 0.2e1 + t547 / 0.2e1 + t548 / 0.2e1 - mrSges(5,3) * t145 + t564 / 0.2e1 - t566) * t210 + t371 * t552 + t252 * (-Ifges(6,5) * t358 - Ifges(6,6) * t359) / 0.2e1 + (Ifges(7,1) * t63 + Ifges(7,4) * t64) * t503 + (Ifges(7,5) * t63 + Ifges(7,6) * t64) * t495 + t221 * (-Ifges(6,4) * t358 - Ifges(6,2) * t359) / 0.2e1 + (m(4) * ((-t242 * t343 - t243 * t338) * qJD(3) + t528) + t343 * t246 - t338 * t247 - t291 * t418 - t293 * t419) * t304 + (-t242 * t418 - t243 * t419 + t528) * mrSges(4,3) + (t356 + t343 * (-Ifges(4,2) * t338 + t468)) * t413 / 0.2e1 + (t410 + t398) * t269 / 0.2e1 + (t357 + t369 * qJD(3) / 0.2e1) * qJD(3) - (-m(5) * t55 + m(6) * t51 + t542) * t539 + t29 * (mrSges(7,1) * t193 - mrSges(7,2) * t194) + t276 * t468 / 0.2e1 + (t233 * mrSges(5,2) - t55 * mrSges(5,3) + Ifges(5,1) * t169 + Ifges(5,4) * t170 + Ifges(5,5) * t328 + t124 * t391 + t368 * t501 + t370 * t507 + t372 * t508 + t374 * t51) * t271 + (-m(5) * t144 + m(6) * t129 - t540) * t117 + (-t15 * t445 - t16 * t444 + t358 * t84 - t359 * t85) * mrSges(6,3) + (m(4) * t305 + t294) * t289 - t194 * t516 - t193 * t517 + t63 * t509 + t64 * t511 + t444 * t513 + (-Ifges(6,1) * t358 - Ifges(6,4) * t359) * t497 + (Ifges(7,4) * t63 + Ifges(7,2) * t64) * t505 - t255 * t419 / 0.2e1 + t343 * (Ifges(4,4) * t276 + Ifges(4,2) * t275) / 0.2e1 + m(7) * (t158 * t29 + t19 * t6 + t20 * t5 + t3 * t45 + t4 * t44 + t74 * t99) + (-Ifges(7,1) * t194 - Ifges(7,4) * t193) * t515 + (-Ifges(7,5) * t194 - Ifges(7,6) * t193) * t502 + (-t19 * t63 - t193 * t3 + t194 * t4 + t20 * t64) * mrSges(7,3) + (-Ifges(7,4) * t194 - Ifges(7,2) * t193) * t514 - t42 * t445 / 0.2e1 + (mrSges(5,1) * t233 - mrSges(5,3) * t54 - Ifges(5,4) * t169 + Ifges(6,5) * t508 + Ifges(7,5) * t515 - Ifges(5,2) * t170 - Ifges(5,6) * t328 + Ifges(6,6) * t507 + Ifges(7,6) * t514 + Ifges(6,3) * t501 + Ifges(7,3) * t502 + t522 + t523) * t269 + m(6) * (t112 * t16 + t113 * t15 + t38 * t85 + t39 * t84) + m(5) * (t116 * t145 + t198 * t54 + t233 * t286 + t259 * t406) + t256 * t418 / 0.2e1 - t123 * t396 / 0.2e1 + qJDD(3) * (Ifges(4,5) * t338 + Ifges(4,6) * t343) + t305 * (-mrSges(4,1) * t275 + mrSges(4,2) * t276) + t286 * (-mrSges(5,1) * t170 + mrSges(5,2) * t169) + t116 * t234 + t44 * t24 + t45 * t25 + t202 * t406 + t129 * (mrSges(6,1) * t359 - mrSges(6,2) * t358) + (-t441 * t518 - mrSges(2,1) * t344 - t239 * mrSges(6,1) - t226 * mrSges(7,1) + mrSges(2,2) * t339 - t238 * mrSges(6,2) - t225 * mrSges(7,2) + t526 * (t317 * t313 - t316 * t346 + t327) + t553 * t327 + t521 * t316 + (-t520 + t555) * t317) * g(2) + (Ifges(4,1) * t276 + Ifges(4,4) * t552) * t338 + t74 * t80 + t99 * (-mrSges(7,1) * t64 + mrSges(7,2) * t63) + t5 * t110 + t6 * t111 + t112 * t70 + t113 * t71 + (-t439 * t518 + mrSges(2,1) * t339 - t237 * mrSges(6,1) - t224 * mrSges(7,1) + mrSges(2,2) * t344 - t236 * mrSges(6,2) - t223 * mrSges(7,2) - t553 * t487 + t526 * (-t317 * t346 - t487) + t521 * t317 + (-m(7) * (-t313 - t537) - m(6) * (-t313 - t535) + m(5) * t313 + t520) * t316) * g(1) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t334 - 0.2e1 * mrSges(3,2) * t333 + m(3) * (t333 ^ 2 + t334 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + t158 * t13 + t38 * t171 + t39 * t172 + t198 * t150; m(3) * qJDD(2) + t63 * t110 + t64 * t111 - t193 * t24 - t194 * t25 + t338 * t246 + t343 * t247 + (-t291 * t338 + t293 * t343) * qJD(3) + (t13 + t542) * t269 + (t80 - t540) * t210 + (t171 * t341 - t172 * t336 + t234) * t209 + (-t456 + t454 + t150 + (-t171 * t336 - t172 * t341) * qJD(5)) * t271 + (t526 + t553) * g(3) + m(4) * (t187 * t338 + t188 * t343 + (-t242 * t338 + t243 * t343) * qJD(3)) + m(6) * (t129 * t210 + t209 * t366 + t269 * t51 + t271 * t350) + m(7) * (t19 * t64 - t193 * t4 - t194 * t3 + t20 * t63 + t210 * t99 + t269 * t29) + m(5) * (-t144 * t210 + t145 * t209 - t269 * t55 + t271 * t54); (t312 * t51 + (t129 * t337 + t342 * t366) * qJD(4) * pkin(3) - g(1) * t283 - g(2) * t282 - t129 * t154 - t84 * t87 - t85 * t88) * m(6) - (-Ifges(4,2) * t422 + t256 + t314) * t421 / 0.2e1 + t519 * t310 + t348 + (-t293 + t399) * t242 + (-t155 + t404) * t234 + (-t357 - t356 * qJD(1) / 0.2e1) * qJD(1) + (t316 * t534 + t529) * g(2) + (t317 * t534 + t530) * g(1) + t538 * t80 + (m(5) * t485 + t375 + t377) * t527 + ((t337 * t54 + t342 * t55 + (-t144 * t337 + t145 * t342) * qJD(4)) * pkin(3) + t144 * t154 - t145 * t155 - t259 * t407) * m(5) + t255 * t422 / 0.2e1 + (t291 + t400) * t243 - t369 * t413 / 0.2e1 - t202 * t407 + t149 * t484 + t150 * t486 + t71 * t443 + Ifges(4,3) * qJDD(3) + t312 * t52 + Ifges(4,5) * t276 + t287 * t13 + Ifges(4,6) * t275 + t550 * t110 + t551 * t111 + (t19 * t551 + t199 * t4 + t20 * t550 + t200 * t3 + t287 * t29 + t538 * t99) * m(7) + t574 * t172 + t575 * t171 - t540 * t561 + (-m(5) * t326 - m(6) * (t326 + t535) - m(7) * (t326 + t537) + t294 + t524) * g(3) - t187 * mrSges(4,2) + t188 * mrSges(4,1) + t199 * t24 + t200 * t25; -pkin(4) * t52 - t311 * t13 - t144 * t234 - t90 * t171 - t89 * t172 + t227 * t24 + t228 * t25 + t348 + t532 * t80 + t527 * t375 + t540 * t145 + t544 * t111 + t543 * t110 + (t316 * t533 + t529) * g(2) + (t533 * t317 + t530) * g(1) + (t19 * t544 + t20 * t543 + t227 * t4 + t228 * t3 - t29 * t311 + t532 * t99) * m(7) + (-t129 * t145 - t84 * t89 - t85 * t90 - pkin(4) * t51 - g(1) * (-t317 * t483 + t283) - g(2) * (-t316 * t483 + t282)) * m(6) + (t524 + t555) * g(3) + (t454 + t519) * pkin(9); (-t335 * t462 - t21) * t111 + (-mrSges(6,2) * t237 + t236 * t545 - t429) * g(2) + t360 - (Ifges(7,4) * t504 + Ifges(7,2) * t506 + Ifges(7,6) * t496 + t512 - t558) * t139 + (Ifges(7,1) * t504 + Ifges(7,4) * t506 + Ifges(7,5) * t496 + t510 - t557) * t384 + (-Ifges(6,2) * t222 + t124 + t213) * t499 + (mrSges(6,2) * t239 - t238 * t545 - t428) * g(1) + (t471 - t171) * t84 + (t470 + t172) * t85 + (t3 * t335 + t340 * t4 + (-t19 * t335 + t20 * t340) * qJD(6)) * t518 + (Ifges(6,5) * t221 - Ifges(6,6) * t222) * t494 + t123 * t497 + (Ifges(6,1) * t221 - t466) * t498 - t80 * t482 - m(7) * (t19 * t21 + t20 * t22 + t482 * t99) + (t340 * t462 - t22) * t110 + t398 + (t336 * t518 - t373 + t374) * t479 + t522 + (t24 * t340 + t25 * t335) * pkin(5) - t129 * (mrSges(6,1) * t222 + mrSges(6,2) * t221); -t99 * (mrSges(7,1) * t139 + mrSges(7,2) * t384) + (Ifges(7,1) * t384 - t463) * t504 + t68 * t503 + (Ifges(7,5) * t384 - Ifges(7,6) * t139) * t496 - t19 * t110 + t20 * t111 - g(1) * t428 - g(2) * t429 - t373 * t479 + (t139 * t20 + t19 * t384) * mrSges(7,3) + t360 + (-Ifges(7,2) * t139 + t134 + t69) * t506;];
tau  = t1;
