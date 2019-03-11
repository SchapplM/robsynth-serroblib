% Calculate vector of inverse dynamics joint torques for
% S6PRRRRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:37:35
% EndTime: 2019-03-09 00:38:18
% DurationCPUTime: 23.23s
% Computational Cost: add. (17152->762), mult. (39333->1051), div. (0->0), fcn. (30707->18), ass. (0->364)
t336 = sin(qJ(3));
t343 = -pkin(9) - pkin(8);
t297 = t343 * t336;
t341 = cos(qJ(3));
t298 = t343 * t341;
t335 = sin(qJ(4));
t340 = cos(qJ(4));
t224 = t335 * t297 - t340 * t298;
t277 = t335 * t341 + t336 * t340;
t412 = qJD(3) * t343;
t285 = t336 * t412;
t286 = t341 * t412;
t342 = cos(qJ(2));
t331 = sin(pkin(6));
t445 = qJD(1) * t331;
t410 = t342 * t445;
t543 = -qJD(4) * t224 + t277 * t410 - t285 * t335 + t340 * t286;
t276 = -t335 * t336 + t340 * t341;
t436 = qJD(4) * t340;
t437 = qJD(4) * t335;
t542 = -t276 * t410 + t340 * t285 + t335 * t286 + t297 * t436 + t298 * t437;
t358 = t277 * qJD(4);
t209 = -qJD(3) * t277 - t358;
t573 = -pkin(10) * t209 - t542;
t357 = t276 * qJD(4);
t208 = qJD(3) * t276 + t357;
t572 = -pkin(10) * t208 + t543;
t571 = mrSges(6,2) - mrSges(7,3);
t333 = sin(qJ(6));
t338 = cos(qJ(6));
t389 = -mrSges(7,1) * t338 + mrSges(7,2) * t333;
t570 = -mrSges(6,1) + t389;
t334 = sin(qJ(5));
t339 = cos(qJ(5));
t223 = t340 * t297 + t298 * t335;
t167 = -pkin(10) * t277 + t223;
t168 = pkin(10) * t276 + t224;
t378 = t339 * t167 - t168 * t334;
t569 = -qJD(5) * t378 - t572 * t334 + t339 * t573;
t439 = qJD(3) * t336;
t319 = pkin(3) * t439;
t176 = -pkin(4) * t209 + t319;
t337 = sin(qJ(2));
t411 = t337 * t445;
t375 = t339 * t276 - t277 * t334;
t91 = qJD(5) * t375 + t208 * t339 + t209 * t334;
t203 = t276 * t334 + t277 * t339;
t92 = qJD(5) * t203 + t208 * t334 - t339 * t209;
t568 = pkin(5) * t92 - pkin(11) * t91 + t176 - t411;
t269 = t276 * qJD(2);
t270 = t277 * qJD(2);
t393 = t339 * t269 - t270 * t334;
t166 = Ifges(6,4) * t393;
t376 = t269 * t334 + t339 * t270;
t327 = qJD(3) + qJD(4);
t322 = qJD(5) + t327;
t489 = Ifges(6,5) * t322;
t107 = Ifges(6,1) * t376 + t166 + t489;
t507 = pkin(3) * t341;
t316 = pkin(2) + t507;
t256 = -qJD(2) * t316 - t410;
t200 = -pkin(4) * t269 + t256;
t388 = mrSges(7,1) * t333 + mrSges(7,2) * t338;
t289 = qJD(2) * pkin(8) + t411;
t394 = pkin(9) * qJD(2) + t289;
t332 = cos(pkin(6));
t444 = qJD(1) * t332;
t409 = t336 * t444;
t216 = t341 * t394 + t409;
t204 = t335 * t216;
t307 = t341 * t444;
t215 = -t336 * t394 + t307;
t207 = qJD(3) * pkin(3) + t215;
t131 = t340 * t207 - t204;
t257 = t270 * pkin(10);
t115 = t131 - t257;
t102 = pkin(4) * t327 + t115;
t206 = t340 * t216;
t132 = t207 * t335 + t206;
t503 = pkin(10) * t269;
t116 = t132 + t503;
t472 = t116 * t334;
t58 = t102 * t339 - t472;
t53 = -pkin(5) * t322 - t58;
t364 = t53 * t388;
t151 = t322 * t338 - t333 * t376;
t150 = Ifges(7,4) * t151;
t152 = t322 * t333 + t338 * t376;
t169 = qJD(6) - t393;
t75 = t152 * Ifges(7,1) + t169 * Ifges(7,5) + t150;
t476 = t338 * t75;
t501 = t58 * mrSges(6,3);
t509 = t333 / 0.2e1;
t482 = t152 * Ifges(7,4);
t74 = t151 * Ifges(7,2) + t169 * Ifges(7,6) + t482;
t567 = -t476 / 0.2e1 + t74 * t509 + t501 - t107 / 0.2e1 + (-Ifges(6,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t376 - t166 / 0.2e1 - t364 - t200 * mrSges(6,2) - t489 / 0.2e1;
t487 = Ifges(6,6) * t322;
t492 = Ifges(6,4) * t376;
t106 = Ifges(6,2) * t393 + t487 + t492;
t456 = t339 * t116;
t59 = t102 * t334 + t456;
t500 = t59 * mrSges(6,3);
t54 = pkin(11) * t322 + t59;
t85 = -pkin(5) * t393 - pkin(11) * t376 + t200;
t23 = t333 * t85 + t338 * t54;
t554 = t23 * mrSges(7,2);
t22 = -t333 * t54 + t338 * t85;
t555 = t22 * mrSges(7,1);
t485 = Ifges(7,3) * t169;
t486 = Ifges(7,6) * t151;
t488 = Ifges(7,5) * t152;
t73 = t485 + t486 + t488;
t566 = t554 + t500 + t106 / 0.2e1 - t73 / 0.2e1 + t492 / 0.2e1 - t200 * mrSges(6,1) + t487 / 0.2e1 - t555;
t109 = t167 * t334 + t168 * t339;
t547 = -qJD(5) * t109 + t334 * t573 + t572 * t339;
t474 = mrSges(6,1) * t322 + mrSges(7,1) * t151 - mrSges(7,2) * t152 - mrSges(6,3) * t376;
t473 = cos(pkin(12));
t396 = t473 * t342;
t330 = sin(pkin(12));
t463 = t330 * t337;
t263 = -t332 * t463 + t396;
t329 = qJ(3) + qJ(4);
t323 = sin(t329);
t324 = cos(t329);
t464 = t330 * t331;
t563 = -t263 * t323 + t324 * t464;
t461 = t331 * t337;
t562 = -t323 * t461 + t324 * t332;
t528 = -m(4) * pkin(8) + m(5) * t343 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3) - mrSges(6,3) - t388;
t325 = qJ(5) + t329;
t311 = sin(t325);
t312 = cos(t325);
t391 = -mrSges(4,1) * t341 + mrSges(4,2) * t336;
t523 = m(7) * pkin(11);
t524 = m(7) * pkin(5);
t561 = -m(4) * pkin(2) - m(5) * t316 - t324 * mrSges(5,1) + t323 * mrSges(5,2) - mrSges(3,1) + t391 + (-t524 + t570) * t312 + (-t523 + t571) * t311;
t558 = -m(7) - m(6);
t431 = qJD(2) * qJD(3);
t287 = qJDD(2) * t341 - t336 * t431;
t557 = t287 / 0.2e1;
t288 = qJDD(2) * t336 + t341 * t431;
t556 = t288 / 0.2e1;
t326 = qJDD(3) + qJDD(4);
t321 = qJDD(5) + t326;
t155 = qJD(2) * t357 + t287 * t335 + t288 * t340;
t156 = -qJD(2) * t358 + t287 * t340 - t288 * t335;
t79 = qJD(5) * t393 + t155 * t339 + t156 * t334;
t48 = qJD(6) * t151 + t321 * t333 + t338 * t79;
t49 = -qJD(6) * t152 + t321 * t338 - t333 * t79;
t18 = -mrSges(7,1) * t49 + mrSges(7,2) * t48;
t68 = mrSges(6,1) * t321 - mrSges(6,3) * t79;
t553 = t18 - t68;
t552 = Ifges(4,2) * t341;
t525 = m(5) * pkin(3);
t551 = -mrSges(4,1) - t525;
t243 = -pkin(4) * t276 - t316;
t103 = -pkin(5) * t375 - pkin(11) * t203 + t243;
t57 = t103 * t333 + t109 * t338;
t549 = -qJD(6) * t57 + t333 * t569 + t338 * t568;
t56 = t103 * t338 - t109 * t333;
t548 = qJD(6) * t56 + t333 * t568 - t338 * t569;
t235 = -t311 * t461 + t312 * t332;
t236 = t311 * t332 + t312 * t461;
t541 = t570 * t235 + t236 * t571;
t190 = -t263 * t311 + t312 * t464;
t191 = t263 * t312 + t311 * t464;
t540 = t570 * t190 + t191 * t571;
t397 = t473 * t337;
t462 = t330 * t342;
t261 = t332 * t397 + t462;
t398 = t331 * t473;
t188 = -t261 * t311 - t312 * t398;
t189 = t261 * t312 - t311 * t398;
t539 = t570 * t188 + t189 * t571;
t432 = qJD(6) * t338;
t368 = t203 * t432 + t333 * t91;
t443 = qJD(2) * t331;
t404 = qJD(1) * t443;
t300 = t342 * t404;
t430 = qJDD(1) * t331;
t248 = t337 * t430 + t300;
t240 = qJDD(2) * pkin(8) + t248;
t429 = qJDD(1) * t332;
t141 = qJD(3) * t307 + t341 * t240 - t289 * t439 + t336 * t429;
t232 = t289 * t341 + t409;
t142 = -qJD(3) * t232 - t240 * t336 + t341 * t429;
t538 = t141 * t341 - t142 * t336;
t80 = -qJD(5) * t376 - t155 * t334 + t156 * t339;
t78 = qJDD(6) - t80;
t20 = mrSges(7,1) * t78 - mrSges(7,3) * t48;
t21 = -mrSges(7,2) * t78 + mrSges(7,3) * t49;
t537 = -t333 * t20 + t338 * t21;
t122 = -mrSges(6,1) * t393 + mrSges(6,2) * t376;
t197 = -mrSges(5,1) * t269 + mrSges(5,2) * t270;
t534 = t391 * qJD(2) + t122 + t197;
t157 = -mrSges(6,2) * t322 + mrSges(6,3) * t393;
t98 = -mrSges(7,2) * t169 + mrSges(7,3) * t151;
t99 = mrSges(7,1) * t169 - mrSges(7,3) * t152;
t533 = -t333 * t99 + t338 * t98 + t157;
t532 = -t562 * mrSges(5,1) - (-t323 * t332 - t324 * t461) * mrSges(5,2) + t541;
t531 = -t563 * mrSges(5,1) - (-t263 * t324 - t323 * t464) * mrSges(5,2) + t540;
t353 = -t261 * t323 - t324 * t398;
t530 = -t353 * mrSges(5,1) - (-t261 * t324 + t323 * t398) * mrSges(5,2) + t539;
t129 = qJDD(3) * pkin(3) - pkin(9) * t288 + t142;
t130 = pkin(9) * t287 + t141;
t42 = -qJD(4) * t132 + t340 * t129 - t130 * t335;
t32 = pkin(4) * t326 - pkin(10) * t155 + t42;
t41 = t335 * t129 + t340 * t130 + t207 * t436 - t216 * t437;
t36 = pkin(10) * t156 + t41;
t9 = -qJD(5) * t59 + t32 * t339 - t334 * t36;
t299 = t337 * t404;
t247 = t342 * t430 - t299;
t239 = -qJDD(2) * pkin(2) - t247;
t199 = -pkin(3) * t287 + t239;
t117 = -pkin(4) * t156 + t199;
t19 = -pkin(5) * t80 - pkin(11) * t79 + t117;
t434 = qJD(5) * t339;
t435 = qJD(5) * t334;
t8 = t102 * t434 - t116 * t435 + t334 * t32 + t339 * t36;
t5 = pkin(11) * t321 + t8;
t2 = qJD(6) * t22 + t19 * t333 + t338 * t5;
t3 = -qJD(6) * t23 + t19 * t338 - t333 * t5;
t529 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t123 = pkin(5) * t376 - pkin(11) * t393;
t382 = t22 * t338 + t23 * t333;
t350 = -qJD(6) * t382 - t3 * t333;
t349 = t2 * t338 + t350;
t433 = qJD(6) * t333;
t526 = m(7) * t349 - t99 * t432 - t98 * t433 + t537;
t344 = qJD(2) ^ 2;
t522 = t48 / 0.2e1;
t521 = t49 / 0.2e1;
t518 = t78 / 0.2e1;
t515 = -t151 / 0.2e1;
t514 = -t152 / 0.2e1;
t513 = t152 / 0.2e1;
t512 = -t169 / 0.2e1;
t510 = t270 / 0.2e1;
t508 = pkin(3) * t340;
t506 = pkin(4) * t270;
t505 = pkin(4) * t334;
t504 = pkin(4) * t339;
t496 = mrSges(7,3) * t333;
t495 = mrSges(7,3) * t338;
t494 = Ifges(4,4) * t336;
t493 = Ifges(4,4) * t341;
t491 = Ifges(7,4) * t333;
t490 = Ifges(7,4) * t338;
t484 = t131 * mrSges(5,3);
t483 = t132 * mrSges(5,3);
t481 = t270 * Ifges(5,4);
t469 = t203 * t333;
t468 = t203 * t338;
t460 = t331 * t341;
t459 = t331 * t342;
t458 = t334 * t335;
t457 = t335 * t339;
t134 = t340 * t215 - t204;
t291 = -pkin(3) * t336 - pkin(4) * t323;
t292 = pkin(4) * t324 + t507;
t449 = t263 * t291 + t292 * t464;
t446 = t291 * t461 + t332 * t292;
t315 = pkin(4) + t508;
t259 = pkin(3) * t457 + t334 * t315;
t442 = qJD(2) * t336;
t441 = qJD(2) * t337;
t440 = qJD(2) * t341;
t438 = qJD(3) * t341;
t428 = Ifges(7,5) * t48 + Ifges(7,6) * t49 + Ifges(7,3) * t78;
t318 = pkin(3) * t442;
t423 = mrSges(4,3) * t442;
t422 = mrSges(4,3) * t440;
t417 = t333 * t459;
t416 = t338 * t459;
t413 = t476 / 0.2e1;
t408 = t331 * t441;
t407 = t342 * t443;
t402 = -t433 / 0.2e1;
t183 = t188 * pkin(5);
t401 = t189 * pkin(11) + t183;
t400 = t190 * pkin(5) + pkin(11) * t191;
t229 = t235 * pkin(5);
t399 = pkin(11) * t236 + t229;
t395 = t431 / 0.2e1;
t133 = -t215 * t335 - t206;
t392 = t563 * pkin(4);
t387 = Ifges(7,1) * t338 - t491;
t386 = t494 + t552;
t385 = -Ifges(7,2) * t333 + t490;
t384 = Ifges(4,5) * t341 - Ifges(4,6) * t336;
t383 = Ifges(7,5) * t338 - Ifges(7,6) * t333;
t381 = -t22 * t333 + t23 * t338;
t264 = t332 * t341 - t336 * t461;
t265 = t332 * t336 + t337 * t460;
t170 = t264 * t340 - t265 * t335;
t171 = t264 * t335 + t265 * t340;
t377 = t339 * t170 - t171 * t334;
t112 = t170 * t334 + t171 * t339;
t374 = t562 * pkin(4);
t258 = -pkin(3) * t458 + t315 * t339;
t370 = t133 - t503;
t369 = t261 * t291 - t292 * t398;
t367 = t203 * t433 - t338 * t91;
t94 = -t112 * t333 - t416;
t366 = -t112 * t338 + t417;
t363 = t151 * t385;
t362 = t152 * t387;
t361 = t169 * t383;
t290 = -qJD(2) * pkin(2) - t410;
t360 = t290 * (mrSges(4,1) * t336 + mrSges(4,2) * t341);
t359 = t336 * (Ifges(4,1) * t341 - t494);
t96 = t123 + t506;
t351 = t353 * pkin(4);
t14 = t48 * Ifges(7,4) + t49 * Ifges(7,2) + t78 * Ifges(7,6);
t15 = t48 * Ifges(7,1) + t49 * Ifges(7,4) + t78 * Ifges(7,5);
t6 = -pkin(5) * t321 - t9;
t346 = t9 * mrSges(6,1) - t8 * mrSges(6,2) + t15 * t509 + t2 * t495 + t6 * t389 + t338 * t14 / 0.2e1 + Ifges(6,3) * t321 + (Ifges(7,1) * t333 + t490) * t522 + (Ifges(7,2) * t338 + t491) * t521 + (Ifges(7,5) * t333 + Ifges(7,6) * t338) * t518 + t74 * t402 + Ifges(6,6) * t80 + Ifges(6,5) * t79 + (t364 + t413) * qJD(6) + (t363 + t362 + t361) * qJD(6) / 0.2e1;
t162 = t269 * Ifges(5,2) + t327 * Ifges(5,6) + t481;
t255 = Ifges(5,4) * t269;
t163 = t270 * Ifges(5,1) + t327 * Ifges(5,5) + t255;
t345 = (-t22 * t432 - t23 * t433) * mrSges(7,3) + t346 - t256 * (mrSges(5,1) * t270 + mrSges(5,2) * t269) + Ifges(5,3) * t326 - t327 * (Ifges(5,5) * t269 - Ifges(5,6) * t270) / 0.2e1 - t3 * t496 + Ifges(5,5) * t155 + Ifges(5,6) * t156 - t270 * (Ifges(5,1) * t269 - t481) / 0.2e1 - t41 * mrSges(5,2) + t42 * mrSges(5,1) + t270 * t483 + t269 * t484 + t162 * t510 - (-Ifges(5,2) * t270 + t163 + t255) * t269 / 0.2e1 + (Ifges(7,5) * t514 + Ifges(7,6) * t515 + Ifges(7,3) * t512 + t566) * t376 + (t22 * t495 + t23 * t496 + t383 * t512 + t385 * t515 + t387 * t514 + t567) * t393;
t328 = -pkin(10) + t343;
t317 = Ifges(4,4) * t440;
t314 = -pkin(5) - t504;
t294 = -qJD(3) * mrSges(4,2) + t422;
t293 = qJD(3) * mrSges(4,1) - t423;
t284 = pkin(2) + t292;
t268 = Ifges(4,1) * t442 + Ifges(4,5) * qJD(3) + t317;
t267 = Ifges(4,6) * qJD(3) + qJD(2) * t386;
t262 = t332 * t462 + t397;
t260 = -t332 * t396 + t463;
t253 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t288;
t252 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t287;
t250 = -pkin(5) - t258;
t234 = mrSges(5,1) * t327 - mrSges(5,3) * t270;
t233 = -mrSges(5,2) * t327 + mrSges(5,3) * t269;
t231 = -t289 * t336 + t307;
t230 = t318 + t506;
t219 = -mrSges(4,1) * t287 + mrSges(4,2) * t288;
t212 = qJD(3) * t264 + t341 * t407;
t211 = -qJD(3) * t265 - t336 * t407;
t146 = -mrSges(5,2) * t326 + mrSges(5,3) * t156;
t145 = mrSges(5,1) * t326 - mrSges(5,3) * t155;
t120 = -t257 + t134;
t93 = -mrSges(5,1) * t156 + mrSges(5,2) * t155;
t90 = t318 + t96;
t88 = -qJD(4) * t171 + t211 * t340 - t212 * t335;
t87 = qJD(4) * t170 + t211 * t335 + t212 * t340;
t69 = -mrSges(6,2) * t321 + mrSges(6,3) * t80;
t65 = t339 * t120 + t334 * t370;
t63 = t115 * t339 - t472;
t62 = t115 * t334 + t456;
t38 = t123 * t333 + t338 * t58;
t37 = t123 * t338 - t333 * t58;
t31 = t333 * t90 + t338 * t65;
t30 = -t333 * t65 + t338 * t90;
t29 = t333 * t96 + t338 * t63;
t28 = -t333 * t63 + t338 * t96;
t26 = -mrSges(6,1) * t80 + mrSges(6,2) * t79;
t25 = qJD(5) * t112 + t334 * t87 - t339 * t88;
t24 = qJD(5) * t377 + t334 * t88 + t339 * t87;
t17 = qJD(6) * t366 - t24 * t333 + t338 * t408;
t16 = qJD(6) * t94 + t24 * t338 + t333 * t408;
t1 = [m(2) * qJDD(1) + t112 * t69 + t170 * t145 + t171 * t146 + t24 * t157 + t16 * t98 + t17 * t99 + t94 * t20 - t366 * t21 + t211 * t293 + t212 * t294 + t87 * t233 + t88 * t234 + t265 * t252 + t264 * t253 - t474 * t25 - t553 * t377 + (-m(2) - m(3) - m(4) - m(5) + t558) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t344 - t219 - t26 - t93) * t342 + (-mrSges(3,1) * t344 - mrSges(3,2) * qJDD(2) + qJD(2) * t534) * t337) * t331 + m(7) * (t16 * t23 + t17 * t22 - t2 * t366 + t25 * t53 + t3 * t94 - t377 * t6) + m(6) * (t377 * t9 + t112 * t8 + t24 * t59 - t25 * t58 + (-t117 * t342 + t200 * t441) * t331) + m(5) * (t131 * t88 + t132 * t87 + t170 * t42 + t171 * t41 + (-t199 * t342 + t256 * t441) * t331) + m(4) * (t141 * t265 + t142 * t264 + t211 * t231 + t212 * t232 + (-t239 * t342 + t290 * t441) * t331) + m(3) * (qJDD(1) * t332 ^ 2 + (t247 * t342 + t248 * t337) * t331); (t109 * t8 + t117 * t243 + t176 * t200 + t378 * t9 + t547 * t58 - t569 * t59) * m(6) - t569 * t157 + (t558 * (-t262 * t284 - t263 * t328) + t528 * t263 - t561 * t262) * g(1) + (t558 * (-t260 * t284 - t261 * t328) + t528 * t261 - t561 * t260) * g(2) + (t300 - t248) * mrSges(3,2) + (-pkin(2) * t239 - (t290 * t337 + (-t231 * t336 + t232 * t341) * t342) * t445) * m(4) - t92 * t554 + t393 * (Ifges(6,4) * t91 - Ifges(6,2) * t92) / 0.2e1 + t376 * (Ifges(6,1) * t91 - Ifges(6,4) * t92) / 0.2e1 + (-t293 * t438 - t294 * t439 + m(4) * ((-t231 * t341 - t232 * t336) * qJD(3) + t538) - t336 * t253 + t341 * t252) * pkin(8) + (-t231 * t438 - t232 * t439 + t538) * mrSges(4,3) - t368 * t74 / 0.2e1 + (t360 + t384 * qJD(3) / 0.2e1) * qJD(3) + t547 * t474 + (t558 * t284 * t459 + (t561 * t342 + (-t328 * t558 + t528) * t337) * t331) * g(3) + t548 * t98 + t549 * t99 + t542 * t233 + t543 * t234 + (t131 * t543 + t132 * t542 - t199 * t316 + t223 * t42 + t224 * t41 + t256 * t319) * m(5) + (t117 * mrSges(6,2) - t9 * mrSges(6,3) + Ifges(6,1) * t79 + Ifges(6,4) * t80 + Ifges(6,5) * t321 + t383 * t518 + t385 * t521 + t387 * t522 + t388 * t6 + t402 * t75) * t203 + (t299 + t247) * mrSges(3,1) + t197 * t319 + (t2 * t57 + t22 * t549 + t23 * t548 + t3 * t56 - t378 * t6 - t53 * t547) * m(7) - t553 * t378 - (-Ifges(6,6) * t321 - Ifges(6,4) * t79 - Ifges(6,2) * t80 + t117 * mrSges(6,1) + t428 / 0.2e1 - t8 * mrSges(6,3) + Ifges(7,3) * t518 + Ifges(7,6) * t521 + Ifges(7,5) * t522 + t529) * t375 + (Ifges(4,4) * t556 + Ifges(4,2) * t557 - t294 * t410 + t493 * t395) * t341 + (Ifges(4,1) * t288 + Ifges(4,4) * t557 + t293 * t410 - t395 * t552) * t336 + t53 * (mrSges(7,1) * t368 - mrSges(7,2) * t367) + t169 * (-Ifges(7,5) * t367 - Ifges(7,6) * t368 + Ifges(7,3) * t92) / 0.2e1 + t151 * (-Ifges(7,4) * t367 - Ifges(7,2) * t368 + Ifges(7,6) * t92) / 0.2e1 + (-t2 * t469 + t22 * t367 - t23 * t368 - t3 * t468) * mrSges(7,3) + (-m(5) * t256 - m(6) * t200 - t534) * t411 + qJDD(3) * (Ifges(4,5) * t336 + Ifges(4,6) * t341) + Ifges(3,3) * qJDD(2) + t327 * (Ifges(5,5) * t208 + Ifges(5,6) * t209) / 0.2e1 + t322 * (Ifges(6,5) * t91 - Ifges(6,6) * t92) / 0.2e1 - t316 * t93 + t269 * (Ifges(5,4) * t208 + Ifges(5,2) * t209) / 0.2e1 + t256 * (-mrSges(5,1) * t209 + mrSges(5,2) * t208) + t243 * t26 - pkin(2) * t219 + t223 * t145 + t224 * t146 + t208 * t163 / 0.2e1 + t209 * t162 / 0.2e1 + t200 * (mrSges(6,1) * t92 + mrSges(6,2) * t91) - t91 * t501 + t176 * t122 - t92 * t500 - t208 * t484 - t92 * t106 / 0.2e1 + t91 * t107 / 0.2e1 + t109 * t69 + t15 * t468 / 0.2e1 + t92 * t73 / 0.2e1 - t14 * t469 / 0.2e1 + t56 * t20 + t57 * t21 + t268 * t438 / 0.2e1 - t267 * t439 / 0.2e1 + t92 * t555 + t493 * t556 + t386 * t557 + t209 * t483 + t91 * t413 + t239 * t391 + t359 * t395 + (Ifges(5,1) * t208 + Ifges(5,4) * t209) * t510 + (-Ifges(7,1) * t367 - Ifges(7,4) * t368 + Ifges(7,5) * t92) * t513 + (mrSges(5,2) * t199 - mrSges(5,3) * t42 + Ifges(5,1) * t155 + Ifges(5,4) * t156 + Ifges(5,5) * t326) * t277 + (-mrSges(5,1) * t199 + mrSges(5,3) * t41 + Ifges(5,4) * t155 + Ifges(5,2) * t156 + Ifges(5,6) * t326) * t276; (t422 - t294) * t231 + (-m(6) * t369 - m(7) * (t369 + t401) - (-t261 * t341 + t336 * t398) * mrSges(4,2) + t551 * (-t261 * t336 - t341 * t398) + t530) * g(2) + (-(-t263 * t341 - t336 * t464) * mrSges(4,2) - m(6) * t449 - m(7) * (t400 + t449) + t551 * (-t263 * t336 + t330 * t460) + t531) * g(1) + (mrSges(4,2) * t265 - m(6) * t446 - m(7) * (t399 + t446) + t551 * t264 + t532) * g(3) + (-t200 * t230 + t258 * t9 + t259 * t8 - t59 * t65) * m(6) + (-m(6) * t58 + m(7) * t53 - t474) * (-t120 * t334 + t339 * t370 + t315 * t435 + (t335 * t434 + (t334 * t340 + t457) * qJD(4)) * pkin(3)) + (t335 * t146 + t233 * t436 - t234 * t437) * pkin(3) + (-t22 * t30 - t23 * t31 + t250 * t6) * m(7) + t526 * (pkin(11) + t259) + (t423 + t293) * t232 - (-Ifges(4,2) * t442 + t268 + t317) * t440 / 0.2e1 + (m(6) * t59 + m(7) * t381 + t533) * (t315 * t434 + (-t335 * t435 + (t339 * t340 - t458) * qJD(4)) * pkin(3)) + Ifges(4,3) * qJDD(3) + Ifges(4,5) * t288 + Ifges(4,6) * t287 + t258 * t68 + t259 * t69 + t250 * t18 - t134 * t233 - t133 * t234 - t230 * t122 - t344 * t359 / 0.2e1 - t65 * t157 - t141 * mrSges(4,2) + t142 * mrSges(4,1) - t31 * t98 - t30 * t99 + t267 * t442 / 0.2e1 - t384 * t431 / 0.2e1 - t197 * t318 - m(5) * (t131 * t133 + t132 * t134 + t256 * t318) + t145 * t508 + (t335 * t41 + t340 * t42 + (-t131 * t335 + t132 * t340) * qJD(4)) * t525 - qJD(2) * t360 + t345; -t122 * t506 - t131 * t233 + t132 * t234 - t63 * t157 + t314 * t18 - t28 * t99 - t29 * t98 + t68 * t504 + t69 * t505 + t345 + (t314 * t6 + (t334 * t53 + t339 * t381) * qJD(5) * pkin(4) - t22 * t28 - t23 * t29 - t53 * t62) * m(7) + ((t334 * t8 + t339 * t9 + (-t334 * t58 + t339 * t59) * qJD(5)) * pkin(4) - t200 * t506 + t58 * t62 - t59 * t63) * m(6) + t533 * pkin(4) * t434 + (-m(6) * t374 - m(7) * (t374 + t399) + t532) * g(3) + (-m(6) * t351 - m(7) * (t351 + t401) + t530) * g(2) + (-m(6) * t392 - m(7) * (t392 + t400) + t531) * g(1) + t526 * (pkin(11) + t505) + t474 * (-pkin(4) * t435 + t62); t346 + (-t361 / 0.2e1 - t363 / 0.2e1 - t362 / 0.2e1 + t382 * mrSges(7,3) + t567) * t393 + t350 * mrSges(7,3) + t349 * t523 - t6 * t524 + ((-t333 * t98 - t338 * t99) * qJD(6) + (-g(2) * t189 - g(3) * t236) * m(7) + t537) * pkin(11) + (-m(7) * t183 + t539) * g(2) + (-m(7) * t400 + t540) * g(1) + (-m(7) * t229 + t541) * g(3) - t58 * t157 - t38 * t98 - t37 * t99 - pkin(5) * t18 + (-t485 / 0.2e1 - t486 / 0.2e1 - t488 / 0.2e1 + t566) * t376 - m(7) * (t22 * t37 + t23 * t38 + t53 * t59) + t474 * t59; -t53 * (mrSges(7,1) * t152 + mrSges(7,2) * t151) + (Ifges(7,1) * t151 - t482) * t514 + t74 * t513 + (Ifges(7,5) * t151 - Ifges(7,6) * t152) * t512 - t22 * t98 + t23 * t99 - g(1) * ((-t191 * t333 + t262 * t338) * mrSges(7,1) + (-t191 * t338 - t262 * t333) * mrSges(7,2)) - g(2) * ((-t189 * t333 + t260 * t338) * mrSges(7,1) + (-t189 * t338 - t260 * t333) * mrSges(7,2)) - g(3) * ((-t236 * t333 - t416) * mrSges(7,1) + (-t236 * t338 + t417) * mrSges(7,2)) + (t151 * t22 + t152 * t23) * mrSges(7,3) + t428 + (-Ifges(7,2) * t152 + t150 + t75) * t515 + t529;];
tau  = t1;
