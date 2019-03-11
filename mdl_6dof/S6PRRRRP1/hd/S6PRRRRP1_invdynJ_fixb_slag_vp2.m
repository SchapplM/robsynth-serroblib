% Calculate vector of inverse dynamics joint torques for
% S6PRRRRP1
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
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:55:47
% EndTime: 2019-03-08 23:56:26
% DurationCPUTime: 22.05s
% Computational Cost: add. (10270->719), mult. (22924->934), div. (0->0), fcn. (17267->14), ass. (0->338)
t624 = Ifges(6,4) + Ifges(7,4);
t338 = sin(qJ(2));
t332 = sin(pkin(6));
t447 = qJD(1) * t332;
t416 = t338 * t447;
t337 = sin(qJ(3));
t441 = qJD(3) * t337;
t553 = pkin(3) * t441 - t416;
t341 = cos(qJ(3));
t435 = qJD(2) * qJD(3);
t293 = qJDD(2) * t341 - t337 * t435;
t294 = qJDD(2) * t337 + t341 * t435;
t336 = sin(qJ(4));
t340 = cos(qJ(4));
t285 = t336 * t337 - t340 * t341;
t352 = t285 * qJD(4);
t144 = -qJD(2) * t352 + t293 * t336 + t294 * t340;
t286 = t336 * t341 + t337 * t340;
t281 = t286 * qJD(2);
t330 = qJD(3) + qJD(4);
t335 = sin(qJ(5));
t339 = cos(qJ(5));
t232 = -t281 * t335 + t330 * t339;
t329 = qJDD(3) + qJDD(4);
t74 = qJD(5) * t232 + t144 * t339 + t329 * t335;
t537 = t74 / 0.2e1;
t233 = t281 * t339 + t330 * t335;
t75 = -qJD(5) * t233 - t144 * t335 + t329 * t339;
t536 = t75 / 0.2e1;
t353 = t286 * qJD(4);
t145 = -qJD(2) * t353 + t293 * t340 - t294 * t336;
t143 = qJDD(5) - t145;
t534 = t143 / 0.2e1;
t585 = Ifges(6,1) + Ifges(7,1);
t583 = Ifges(6,5) + Ifges(7,5);
t609 = Ifges(6,2) + Ifges(7,2);
t581 = Ifges(7,6) + Ifges(6,6);
t580 = Ifges(7,3) + Ifges(6,3);
t220 = -qJD(3) * t285 - t352;
t221 = qJD(3) * t286 + t353;
t623 = pkin(4) * t221 - pkin(10) * t220 + t553;
t343 = -pkin(9) - pkin(8);
t417 = qJD(3) * t343;
t291 = t337 * t417;
t292 = t341 * t417;
t342 = cos(qJ(2));
t415 = t342 * t447;
t307 = t343 * t337;
t308 = t343 * t341;
t554 = t340 * t307 + t308 * t336;
t556 = qJD(4) * t554 + t285 * t415 + t291 * t340 + t292 * t336;
t613 = mrSges(5,2) - mrSges(6,3) - mrSges(7,3);
t586 = -mrSges(6,2) - mrSges(7,2);
t620 = mrSges(6,1) + mrSges(7,1);
t622 = t335 * t586 + t339 * t620 + mrSges(5,1);
t621 = t583 * t534 + t624 * t536 + t537 * t585;
t619 = t624 * t232;
t436 = qJD(5) * t339;
t411 = t286 * t436;
t357 = t220 * t335 + t411;
t618 = t624 * t233;
t617 = t624 * t339;
t616 = t624 * t335;
t295 = qJD(2) * pkin(8) + t416;
t385 = pkin(9) * qJD(2) + t295;
t333 = cos(pkin(6));
t446 = qJD(1) * t333;
t414 = t337 * t446;
t226 = t341 * t385 + t414;
t213 = t336 * t226;
t316 = t341 * t446;
t225 = -t337 * t385 + t316;
t215 = qJD(3) * pkin(3) + t225;
t109 = t215 * t340 - t213;
t322 = pkin(3) * t341 + pkin(2);
t267 = -qJD(2) * t322 - t415;
t574 = t330 * Ifges(5,5);
t615 = t267 * mrSges(5,2) - t109 * mrSges(5,3) + t574 / 0.2e1;
t280 = t285 * qJD(2);
t271 = qJD(5) + t280;
t214 = t340 * t226;
t110 = t215 * t336 + t214;
t104 = pkin(10) * t330 + t110;
t140 = pkin(4) * t280 - pkin(10) * t281 + t267;
t45 = -t104 * t335 + t339 * t140;
t33 = -qJ(6) * t233 + t45;
t30 = pkin(5) * t271 + t33;
t46 = t104 * t339 + t140 * t335;
t34 = qJ(6) * t232 + t46;
t573 = t330 * Ifges(5,6);
t614 = -t267 * mrSges(5,1) - t45 * mrSges(6,1) - t30 * mrSges(7,1) + t46 * mrSges(6,2) + t34 * mrSges(7,2) + t573 / 0.2e1;
t545 = -m(4) * pkin(8) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t320 = pkin(5) * t339 + pkin(4);
t331 = qJ(3) + qJ(4);
t326 = sin(t331);
t327 = cos(t331);
t334 = -qJ(6) - pkin(10);
t375 = -mrSges(4,1) * t341 + mrSges(4,2) * t337;
t540 = m(4) * pkin(2) + mrSges(3,1) - t375 + (m(6) * pkin(4) + m(7) * t320 + mrSges(5,1)) * t327 + (m(6) * pkin(10) - m(7) * t334 - t613) * t326;
t578 = t232 * t609 + t271 * t581 + t618;
t610 = -t578 / 0.2e1;
t608 = t143 * t581 + t609 * t75 + t624 * t74;
t606 = t232 * t581 + t233 * t583 + t271 * t580;
t577 = t233 * t585 + t271 * t583 + t619;
t438 = qJD(4) * t340;
t427 = pkin(3) * t438;
t114 = t225 * t340 - t213;
t202 = pkin(4) * t281 + pkin(10) * t280;
t444 = qJD(2) * t337;
t430 = pkin(3) * t444;
t164 = t202 + t430;
t50 = t339 * t114 + t335 * t164;
t605 = t339 * t427 - t50;
t113 = t225 * t336 + t214;
t439 = qJD(4) * t336;
t604 = -pkin(3) * t439 + t113;
t475 = t280 * t335;
t603 = -qJ(6) * t475 + t339 * qJD(6);
t602 = -t335 * t581 + t339 * t583;
t601 = -t335 * t609 + t617;
t600 = t339 * t585 - t616;
t489 = sin(pkin(11));
t388 = t489 * t338;
t490 = cos(pkin(11));
t389 = t490 * t342;
t275 = -t333 * t388 + t389;
t391 = t332 * t489;
t599 = -t275 * t337 + t341 * t391;
t598 = -t335 * t556 + t339 * t623;
t204 = pkin(4) * t285 - pkin(10) * t286 - t322;
t597 = t204 * t436 + t335 * t623 + t339 * t556;
t467 = t332 * t338;
t276 = t333 * t341 - t337 * t467;
t387 = t489 * t342;
t390 = t490 * t338;
t273 = t333 * t390 + t387;
t392 = t332 * t490;
t205 = t273 * t326 + t327 * t392;
t206 = t273 * t327 - t326 * t392;
t596 = t205 * t622 + t613 * t206;
t207 = t275 * t326 - t327 * t391;
t208 = t275 * t327 + t326 * t391;
t595 = t207 * t622 + t613 * t208;
t256 = t326 * t467 - t333 * t327;
t257 = t326 * t333 + t327 * t467;
t594 = t256 * t622 + t613 * t257;
t103 = -pkin(4) * t330 - t109;
t372 = mrSges(7,1) * t335 + mrSges(7,2) * t339;
t373 = mrSges(6,1) * t335 + mrSges(6,2) * t339;
t56 = -pkin(5) * t232 + qJD(6) + t103;
t593 = t103 * t373 + t372 * t56;
t538 = m(7) * pkin(5);
t589 = t293 / 0.2e1;
t588 = t294 / 0.2e1;
t236 = t307 * t336 - t308 * t340;
t224 = t339 * t236;
t360 = -qJ(6) * t220 - qJD(6) * t286;
t579 = pkin(5) * t221 + t360 * t339 + (-t224 + (qJ(6) * t286 - t204) * t335) * qJD(5) + t598;
t576 = Ifges(5,4) * t280;
t575 = t280 * Ifges(5,2);
t572 = t341 * Ifges(4,2);
t570 = mrSges(7,1) + t538;
t112 = t335 * t204 + t224;
t569 = -qJD(5) * t112 + t598;
t568 = -qJ(6) * t411 + (-qJD(5) * t236 + t360) * t335 + t597;
t437 = qJD(5) * t335;
t567 = -t236 * t437 + t597;
t518 = pkin(3) * t336;
t319 = pkin(10) + t518;
t463 = -qJ(6) - t319;
t383 = qJD(5) * t463;
t566 = t335 * t383 + t603 + t605;
t328 = t339 * qJ(6);
t376 = t281 * pkin(5) + t280 * t328;
t49 = -t114 * t335 + t339 * t164;
t565 = (-qJD(6) - t427) * t335 + t339 * t383 - t376 - t49;
t393 = qJD(5) * t334;
t53 = t339 * t109 + t335 * t202;
t564 = t335 * t393 - t53 + t603;
t52 = -t109 * t335 + t339 * t202;
t563 = -qJD(6) * t335 + t339 * t393 - t376 - t52;
t130 = mrSges(5,1) * t329 - mrSges(5,3) * t144;
t29 = -mrSges(6,1) * t75 + mrSges(6,2) * t74;
t562 = t29 - t130;
t258 = pkin(5) * t475;
t426 = pkin(5) * t437;
t561 = t426 + t258 - t604;
t508 = mrSges(5,3) * t281;
t555 = mrSges(5,1) * t330 + mrSges(6,1) * t232 - mrSges(6,2) * t233 - t508;
t552 = t143 * t580 + t581 * t75 + t583 * t74;
t445 = qJD(2) * t332;
t403 = qJD(1) * t445;
t310 = t342 * t403;
t434 = qJDD(1) * t332;
t263 = t338 * t434 + t310;
t255 = qJDD(2) * pkin(8) + t263;
t433 = qJDD(1) * t333;
t120 = qJD(3) * t316 + t341 * t255 - t295 * t441 + t337 * t433;
t244 = t295 * t341 + t414;
t121 = -qJD(3) * t244 - t255 * t337 + t341 * t433;
t551 = t120 * t341 - t121 * t337;
t550 = -t46 * mrSges(6,3) - t34 * mrSges(7,3);
t549 = -mrSges(6,3) * t45 - mrSges(7,3) * t30;
t547 = m(7) + m(6) + m(5);
t546 = -mrSges(6,1) - t570;
t128 = -mrSges(7,1) * t232 + mrSges(7,2) * t233;
t541 = -m(7) * t56 - t128;
t100 = qJDD(3) * pkin(3) - pkin(9) * t294 + t121;
t105 = pkin(9) * t293 + t120;
t24 = t336 * t100 + t340 * t105 + t215 * t438 - t226 * t439;
t21 = pkin(10) * t329 + t24;
t309 = t338 * t403;
t262 = t342 * t434 - t309;
t254 = -qJDD(2) * pkin(2) - t262;
t203 = -pkin(3) * t293 + t254;
t43 = -pkin(4) * t145 - pkin(10) * t144 + t203;
t5 = -t104 * t437 + t140 * t436 + t339 * t21 + t335 * t43;
t3 = qJ(6) * t75 + qJD(6) * t232 + t5;
t6 = -qJD(5) * t46 - t21 * t335 + t339 * t43;
t539 = t6 * mrSges(6,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2);
t344 = qJD(2) ^ 2;
t1 = pkin(5) * t143 - qJ(6) * t74 - qJD(6) * t233 + t6;
t535 = t1 * mrSges(7,3);
t532 = -t232 / 0.2e1;
t531 = t232 / 0.2e1;
t530 = -t233 / 0.2e1;
t529 = t233 / 0.2e1;
t528 = -t271 / 0.2e1;
t527 = t271 / 0.2e1;
t526 = t280 / 0.2e1;
t524 = t281 / 0.2e1;
t517 = pkin(3) * t340;
t515 = pkin(10) * t339;
t513 = t335 * t6;
t512 = t339 * t5;
t507 = mrSges(6,3) * t232;
t506 = mrSges(6,3) * t233;
t505 = mrSges(7,3) * t232;
t504 = mrSges(7,3) * t233;
t503 = Ifges(4,4) * t337;
t502 = Ifges(4,4) * t341;
t496 = t110 * mrSges(5,3);
t493 = t281 * Ifges(5,4);
t481 = t220 * t339;
t478 = t273 * t335;
t477 = t275 * t335;
t474 = t280 * t339;
t473 = t286 * t335;
t472 = t286 * t339;
t470 = t319 * t339;
t469 = t327 * t335;
t468 = t327 * t339;
t465 = t335 * t342;
t464 = t339 * t342;
t461 = -t205 * t320 - t206 * t334;
t460 = -t207 * t320 - t208 * t334;
t453 = -t256 * t320 - t257 * t334;
t443 = qJD(2) * t338;
t442 = qJD(2) * t341;
t440 = qJD(3) * t341;
t424 = mrSges(4,3) * t444;
t423 = mrSges(4,3) * t442;
t421 = t332 * t465;
t419 = t332 * t464;
t413 = t332 * t443;
t412 = t342 * t445;
t28 = -t75 * mrSges(7,1) + t74 * mrSges(7,2);
t398 = -t437 / 0.2e1;
t396 = -t205 * pkin(4) + t206 * pkin(10);
t395 = -t207 * pkin(4) + pkin(10) * t208;
t394 = -t256 * pkin(4) + pkin(10) * t257;
t386 = t435 / 0.2e1;
t111 = t339 * t204 - t236 * t335;
t379 = t599 * pkin(3);
t369 = t503 + t572;
t366 = Ifges(4,5) * t341 - Ifges(4,6) * t337;
t277 = t333 * t337 + t341 * t467;
t363 = t340 * t276 - t277 * t336;
t172 = t276 * t336 + t277 * t340;
t359 = t276 * pkin(3);
t25 = t100 * t340 - t336 * t105 - t215 * t439 - t226 * t438;
t136 = -t172 * t335 - t419;
t358 = -t172 * t339 + t421;
t356 = t286 * t437 - t481;
t296 = -qJD(2) * pkin(2) - t415;
t355 = t296 * (mrSges(4,1) * t337 + mrSges(4,2) * t341);
t354 = t337 * (Ifges(4,1) * t341 - t503);
t350 = -t273 * t337 - t341 * t392;
t22 = -pkin(4) * t329 - t25;
t349 = t350 * pkin(3);
t119 = qJD(4) * t236 + t291 * t336 - t340 * t292;
t346 = -t513 + t512 + (-t335 * t46 - t339 * t45) * qJD(5);
t165 = t493 + t573 - t575;
t166 = t281 * Ifges(5,1) + t574 - t576;
t8 = -pkin(5) * t75 + qJDD(6) + t22;
t345 = t475 * t610 + (-t576 + t166) * t526 + (-t45 * t474 - t46 * t475 + t512) * mrSges(6,3) + (t3 * t339 - t30 * t474 - t34 * t475) * mrSges(7,3) + t22 * (-mrSges(6,1) * t339 + mrSges(6,2) * t335) + t8 * (-mrSges(7,1) * t339 + mrSges(7,2) * t335) + Ifges(5,3) * t329 + Ifges(5,5) * t144 + Ifges(5,6) * t145 + t593 * qJD(5) + t165 * t524 + t578 * t398 - (-Ifges(5,1) * t280 - t493 + t606) * t281 / 0.2e1 + t335 * t621 + t608 * t339 / 0.2e1 + (t339 * t609 + t616) * t536 + (t335 * t583 + t339 * t581) * t534 + (-Ifges(5,2) * t526 + t580 * t528 + t583 * t530 + t581 * t532 + t614) * t281 + (t335 * t585 + t617) * t537 - t24 * mrSges(5,2) + t25 * mrSges(5,1) + (-t602 * t528 - t600 * t530 - t601 * t532 + t593 + t615) * t280 + (t436 / 0.2e1 + t474 / 0.2e1) * t577 + (t232 * t601 + t233 * t600 + t271 * t602) * qJD(5) / 0.2e1;
t323 = Ifges(4,4) * t442;
t321 = -pkin(4) - t517;
t304 = t328 + t515;
t303 = t334 * t335;
t302 = -qJD(3) * mrSges(4,2) + t423;
t301 = qJD(3) * mrSges(4,1) - t424;
t300 = -t320 - t517;
t288 = t375 * qJD(2);
t284 = t328 + t470;
t283 = t463 * t335;
t279 = Ifges(4,1) * t444 + Ifges(4,5) * qJD(3) + t323;
t278 = Ifges(4,6) * qJD(3) + qJD(2) * t369;
t274 = t333 * t387 + t390;
t272 = -t333 * t389 + t388;
t265 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t294;
t264 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t293;
t245 = -mrSges(5,2) * t330 - mrSges(5,3) * t280;
t243 = -t295 * t337 + t316;
t229 = -mrSges(4,1) * t293 + mrSges(4,2) * t294;
t223 = qJD(3) * t276 + t341 * t412;
t222 = -qJD(3) * t277 - t337 * t412;
t199 = mrSges(5,1) * t280 + mrSges(5,2) * t281;
t156 = pkin(5) * t473 - t554;
t149 = mrSges(6,1) * t271 - t506;
t148 = mrSges(7,1) * t271 - t504;
t147 = -mrSges(6,2) * t271 + t507;
t146 = -mrSges(7,2) * t271 + t505;
t131 = -mrSges(5,2) * t329 + mrSges(5,3) * t145;
t84 = -qJ(6) * t473 + t112;
t76 = t110 - t258;
t62 = pkin(5) * t285 - t286 * t328 + t111;
t57 = -mrSges(5,1) * t145 + mrSges(5,2) * t144;
t55 = qJD(4) * t172 - t340 * t222 + t223 * t336;
t54 = qJD(4) * t363 + t222 * t336 + t223 * t340;
t51 = pkin(5) * t357 + t119;
t40 = -mrSges(6,2) * t143 + mrSges(6,3) * t75;
t39 = -mrSges(7,2) * t143 + mrSges(7,3) * t75;
t38 = mrSges(6,1) * t143 - mrSges(6,3) * t74;
t37 = mrSges(7,1) * t143 - mrSges(7,3) * t74;
t32 = qJD(5) * t358 - t335 * t54 + t339 * t413;
t31 = qJD(5) * t136 + t335 * t413 + t339 * t54;
t2 = [m(2) * qJDD(1) + t172 * t131 + t222 * t301 + t223 * t302 + t54 * t245 + t277 * t264 + t276 * t265 + (t149 + t148) * t32 + (t146 + t147) * t31 - (t39 + t40) * t358 + (t37 + t38) * t136 + (t128 - t555) * t55 - (t28 + t562) * t363 + (-m(2) - m(3) - m(4) - t547) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t344 - t229 - t57) * t342 + (-mrSges(3,1) * t344 - mrSges(3,2) * qJDD(2) + (t199 + t288) * qJD(2)) * t338) * t332 + m(5) * (-t109 * t55 + t110 * t54 + t363 * t25 + t172 * t24 + (-t203 * t342 + t267 * t443) * t332) + m(4) * (t120 * t277 + t121 * t276 + t222 * t243 + t223 * t244 + (-t254 * t342 + t296 * t443) * t332) + m(3) * (qJDD(1) * t333 ^ 2 + (t262 * t342 + t263 * t338) * t332) + m(7) * (t1 * t136 - t3 * t358 + t30 * t32 + t31 * t34 - t363 * t8 + t55 * t56) + m(6) * (t103 * t55 + t136 * t6 - t22 * t363 + t31 * t46 + t32 * t45 - t358 * t5); t502 * t588 + t369 * t589 + t472 * t621 + (-t478 * t538 - t547 * (-t272 * t322 - t273 * t343) - t620 * (-t272 * t468 + t478) + t586 * (t272 * t469 + t273 * t339) + t545 * t273 + t540 * t272) * g(2) + (-t477 * t538 - t547 * (-t274 * t322 - t275 * t343) - t620 * (-t274 * t468 + t477) + t586 * (t274 * t469 + t275 * t339) + t545 * t275 + t540 * t274) * g(1) + (t103 * t119 + t111 * t6 + t112 * t5 - t22 * t554 + t45 * t569 + t46 * t567) * m(6) + (-t109 * t119 + t110 * t556 - t203 * t322 + t236 * t24 + t25 * t554 + t267 * t553) * m(5) - t562 * t554 + (-t576 / 0.2e1 + t166 / 0.2e1 + Ifges(5,1) * t524 + t615) * t220 + (-t263 + t310) * mrSges(3,2) + t357 * t610 + (t262 + t309) * mrSges(3,1) + (-t1 * t472 - t3 * t473 + t30 * t356 - t34 * t357) * mrSges(7,3) + (-t356 * t624 - t357 * t609) * t531 + (-t356 * t585 - t357 * t624) * t529 + (t356 * t45 - t357 * t46 - t472 * t6 - t473 * t5) * mrSges(6,3) + ((m(5) * t109 - m(6) * t103 + t541 + t555) * t415 + t577 * t398 + t203 * mrSges(5,2) - t25 * mrSges(5,3) + Ifges(5,1) * t144 + Ifges(5,4) * t145 + Ifges(5,5) * t329 + t22 * t373 + t372 * t8 + t602 * t534 + t601 * t536 + t600 * t537) * t286 + ((-t464 * t620 - t465 * t586) * t327 + (-t322 * t547 - t540) * t342 + (t586 * t339 + t343 * t547 + (-t538 - t620) * t335 + t545) * t338) * g(3) * t332 + (-t496 + t575 / 0.2e1 - t165 / 0.2e1 - Ifges(5,4) * t524 + t581 * t531 + t583 * t529 + t580 * t527 + t606 / 0.2e1 - t614) * t221 + (t203 * mrSges(5,1) + t1 * mrSges(7,1) - t24 * mrSges(5,3) - Ifges(5,4) * t144 - Ifges(5,2) * t145 - Ifges(5,6) * t329 + t534 * t580 + t536 * t581 + t537 * t583 + t539 + t552 / 0.2e1) * t285 + (-t356 * t583 - t357 * t581) * t527 + (m(4) * ((-t243 * t341 - t244 * t337) * qJD(3) + t551) + t341 * t264 - t301 * t440 - t302 * t441 - t337 * t265) * pkin(8) + (-t243 * t440 - t244 * t441 + t551) * mrSges(4,3) + t553 * t199 + t279 * t440 / 0.2e1 - t278 * t441 / 0.2e1 + (t355 + t366 * qJD(3) / 0.2e1) * qJD(3) + t577 * t481 / 0.2e1 + t579 * t148 + (t1 * t62 + t156 * t8 + t3 * t84 + t30 * t579 + t34 * t568 + t51 * t56) * m(7) - t288 * t416 + t354 * t386 + (Ifges(4,4) * t588 + Ifges(4,2) * t589 - t302 * t415 + t502 * t386) * t341 + (Ifges(4,1) * t294 + Ifges(4,4) * t589 + t301 * t415 - t386 * t572) * t337 + t568 * t146 + t569 * t149 + t567 * t147 + Ifges(3,3) * qJDD(2) + qJDD(3) * (Ifges(4,5) * t337 + Ifges(4,6) * t341) + t254 * t375 - t322 * t57 + t103 * (mrSges(6,1) * t357 - mrSges(6,2) * t356) + t56 * (mrSges(7,1) * t357 - mrSges(7,2) * t356) + t236 * t131 - pkin(2) * t229 + t156 * t28 + (-pkin(2) * t254 - (t296 * t338 + (-t243 * t337 + t244 * t341) * t342) * t447) * m(4) + t51 * t128 + t111 * t38 + t112 * t40 - t608 * t473 / 0.2e1 + t62 * t37 - t555 * t119 + t556 * t245 + t84 * t39; (t22 * t321 + (t103 * t336 + (-t335 * t45 + t339 * t46) * t340) * qJD(4) * pkin(3) + t346 * t319 - t103 * t113 - t45 * t49 - t46 * t50) * m(6) + (t109 * t113 - t110 * t114 - t267 * t430 + (t24 * t336 + t25 * t340 + (-t109 * t336 + t110 * t340) * qJD(4)) * pkin(3)) * m(5) + (t427 - t114) * t245 + (t424 + t301) * t244 + (t423 - t302) * t243 + (-t319 * t436 - t335 * t427 - t49) * t149 - (-Ifges(4,2) * t444 + t279 + t323) * t442 / 0.2e1 + (-t319 * t38 - t535) * t335 + t549 * t436 + t550 * t437 + t278 * t444 / 0.2e1 - t366 * t435 / 0.2e1 + t281 * t496 + t40 * t470 - t344 * t354 / 0.2e1 - t199 * t430 - mrSges(6,3) * t513 + t345 + t561 * t128 + t565 * t148 + t566 * t146 + (t1 * t283 + t284 * t3 + t30 * t565 + t300 * t8 + t34 * t566 + t56 * t561) * m(7) + Ifges(4,3) * qJDD(3) + t321 * t29 + Ifges(4,5) * t294 + t300 * t28 + Ifges(4,6) * t293 + t283 * t37 + t284 * t39 - qJD(2) * t355 - t120 * mrSges(4,2) + t121 * mrSges(4,1) + t130 * t517 + t131 * t518 + (-m(5) * t359 - m(7) * (t359 + t453) - m(6) * (t359 + t394) - mrSges(4,1) * t276 + mrSges(4,2) * t277 + t594) * g(3) + (-m(5) * t349 - t350 * mrSges(4,1) - (-t273 * t341 + t337 * t392) * mrSges(4,2) - m(7) * (t349 + t461) - m(6) * (t349 + t396) + t596) * g(2) + (-m(5) * t379 - t599 * mrSges(4,1) - (-t275 * t341 - t337 * t391) * mrSges(4,2) - m(7) * (t379 + t460) - m(6) * (t379 + t395) + t595) * g(1) + t604 * t555 + (-t319 * t437 + t605) * t147; ((-pkin(10) * t149 + t549) * t339 + (pkin(5) * t128 - pkin(10) * t147 + t550) * t335) * qJD(5) + t564 * t146 + t345 + t596 * g(2) + t595 * g(1) + t594 * g(3) + (t555 + t508) * t110 + t563 * t148 + (-t6 * mrSges(6,3) - pkin(10) * t38 - t535) * t335 - t320 * t28 + t303 * t37 + t304 * t39 - t109 * t245 - t52 * t149 - t53 * t147 - t76 * t128 - pkin(4) * t29 + t40 * t515 + (-t460 * g(1) - t461 * g(2) - t453 * g(3) + t1 * t303 + t3 * t304 - t320 * t8 + (t426 - t76) * t56 + t564 * t34 + t563 * t30) * m(7) + (-pkin(4) * t22 + pkin(10) * t346 - g(1) * t395 - g(2) * t396 - g(3) * t394 - t103 * t110 - t45 * t52 - t46 * t53) * m(6); (t586 * (-t206 * t339 - t272 * t335) + t546 * (-t206 * t335 + t272 * t339)) * g(2) + (t586 * (-t208 * t339 - t274 * t335) + t546 * (-t208 * t335 + t274 * t339)) * g(1) + (t586 * (-t257 * t339 + t421) + t546 * (-t257 * t335 - t419)) * g(3) + (t232 * t585 - t618) * t530 + (t232 * t583 - t233 * t581) * t528 + t539 + (-t233 * t609 + t577 + t619) * t532 + t570 * t1 + t578 * t529 + (t149 + t506) * t46 + (-t147 + t507) * t45 - t103 * (mrSges(6,1) * t233 + mrSges(6,2) * t232) - t56 * (mrSges(7,1) * t233 + mrSges(7,2) * t232) + (-m(7) * (-t30 + t33) + t148 + t504) * t34 - t33 * t146 + t30 * t505 + t552 + (t233 * t541 + t37) * pkin(5); -t232 * t146 + t233 * t148 + (-g(1) * t207 - g(2) * t205 - g(3) * t256 - t34 * t232 + t30 * t233 + t8) * m(7) + t28;];
tau  = t2;
