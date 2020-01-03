% Calculate vector of inverse dynamics joint torques for
% S5RRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP11_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP11_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP11_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP11_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:14:27
% EndTime: 2019-12-31 22:15:26
% DurationCPUTime: 32.76s
% Computational Cost: add. (9505->787), mult. (23562->1044), div. (0->0), fcn. (18166->10), ass. (0->349)
t280 = cos(pkin(5));
t402 = qJD(1) * t280;
t267 = qJD(2) + t402;
t282 = sin(qJ(3));
t285 = cos(qJ(3));
t283 = sin(qJ(2));
t279 = sin(pkin(5));
t403 = qJD(1) * t279;
t370 = t283 * t403;
t185 = t267 * t282 + t285 * t370;
t286 = cos(qJ(2));
t391 = qJD(1) * qJD(2);
t222 = (qJDD(1) * t283 + t286 * t391) * t279;
t389 = qJDD(1) * t280;
t266 = qJDD(2) + t389;
t104 = -qJD(3) * t185 - t222 * t282 + t266 * t285;
t99 = qJDD(4) - t104;
t483 = t99 / 0.2e1;
t184 = t267 * t285 - t282 * t370;
t103 = qJD(3) * t184 + t222 * t285 + t266 * t282;
t281 = sin(qJ(4));
t284 = cos(qJ(4));
t366 = t286 * t403;
t311 = -qJD(3) + t366;
t143 = t284 * t185 - t281 * t311;
t221 = (-qJDD(1) * t286 + t283 * t391) * t279;
t206 = qJDD(3) + t221;
t47 = qJD(4) * t143 + t281 * t103 - t284 * t206;
t490 = t47 / 0.2e1;
t491 = -t47 / 0.2e1;
t142 = t281 * t185 + t284 * t311;
t397 = qJD(4) * t142;
t46 = t284 * t103 + t281 * t206 - t397;
t492 = t46 / 0.2e1;
t536 = Ifges(6,4) + Ifges(5,5);
t537 = Ifges(5,1) + Ifges(6,1);
t556 = Ifges(5,4) * t491 + Ifges(6,5) * t490 + t483 * t536 + t492 * t537;
t562 = Ifges(5,6) - Ifges(6,6);
t534 = -Ifges(5,3) - Ifges(6,2);
t481 = t103 / 0.2e1;
t567 = Ifges(4,4) * t481;
t467 = t206 / 0.2e1;
t480 = t104 / 0.2e1;
t179 = qJD(4) - t184;
t386 = pkin(1) * t402;
t390 = qJDD(1) * t279;
t551 = pkin(7) * t390 + qJD(2) * t386;
t552 = -pkin(7) * t279 * t391 + pkin(1) * t389;
t144 = t283 * t552 + t286 * t551;
t123 = pkin(8) * t266 + t144;
t131 = -pkin(1) * t390 + pkin(2) * t221 - pkin(8) * t222;
t275 = t280 * t283 * pkin(1);
t413 = t279 * t286;
t405 = pkin(7) * t413 + t275;
t218 = t405 * qJD(1);
t169 = t267 * pkin(8) + t218;
t177 = (-pkin(2) * t286 - pkin(8) * t283 - pkin(1)) * t403;
t398 = qJD(3) * t285;
t399 = qJD(3) * t282;
t30 = t285 * t123 + t282 * t131 - t169 * t399 + t177 * t398;
t24 = pkin(9) * t206 + t30;
t145 = -t283 * t551 + t286 * t552;
t124 = -t266 * pkin(2) - t145;
t29 = -t104 * pkin(3) - t103 * pkin(9) + t124;
t394 = qJD(4) * t284;
t396 = qJD(4) * t281;
t215 = -pkin(7) * t370 + t286 * t386;
t168 = -t267 * pkin(2) - t215;
t81 = -t184 * pkin(3) - t185 * pkin(9) + t168;
t101 = t285 * t169 + t282 * t177;
t83 = -pkin(9) * t311 + t101;
t3 = t284 * t24 + t281 * t29 + t81 * t394 - t396 * t83;
t1 = qJ(5) * t99 + qJD(5) * t179 + t3;
t34 = t281 * t81 + t284 * t83;
t4 = -qJD(4) * t34 - t24 * t281 + t284 * t29;
t2 = -pkin(4) * t99 + qJDD(5) - t4;
t543 = t4 * mrSges(5,1) - t2 * mrSges(6,1) - t3 * mrSges(5,2) + t1 * mrSges(6,3);
t566 = -0.2e1 * Ifges(4,2) * t480 - 0.2e1 * Ifges(4,6) * t467 + Ifges(5,6) * t491 + Ifges(6,6) * t490 - t534 * t483 + t536 * t492 + t543 - t567;
t299 = t279 * (pkin(2) * t283 - pkin(8) * t286);
t216 = qJD(1) * t299;
t134 = t285 * t215 + t282 * t216;
t565 = pkin(8) * t399 + pkin(9) * t370 + t134;
t336 = pkin(3) * t282 - pkin(9) * t285;
t564 = -pkin(8) * qJD(4) * t285 - (t275 + (pkin(7) + t336) * t413) * qJD(1) + t336 * qJD(3);
t341 = t282 * t366;
t561 = t341 - t399;
t313 = pkin(4) * t284 + qJ(5) * t281;
t331 = -t284 * mrSges(6,1) - t281 * mrSges(6,3);
t333 = mrSges(5,1) * t284 - mrSges(5,2) * t281;
t516 = m(6) * t313 - t331 + t333;
t541 = m(5) + m(6);
t560 = pkin(3) * t541 + mrSges(4,1) + t516;
t559 = -t46 * Ifges(6,5) / 0.2e1 - t99 * Ifges(6,6) / 0.2e1 + Ifges(5,4) * t492 + Ifges(5,6) * t483 + (Ifges(6,3) + Ifges(5,2)) * t491;
t531 = -t142 * t562 + t143 * t536 - t179 * t534;
t557 = t531 / 0.2e1;
t539 = t46 * t536 - t47 * t562 - t534 * t99;
t554 = -Ifges(5,4) + Ifges(6,5);
t140 = Ifges(5,4) * t142;
t434 = Ifges(6,5) * t142;
t530 = t143 * t537 + t179 * t536 - t140 + t434;
t553 = t144 * mrSges(3,2);
t458 = pkin(3) * t285;
t247 = -pkin(9) * t282 - pkin(2) - t458;
t521 = -t247 * t396 + t281 * t565 + t284 * t564;
t520 = t247 * t394 + t281 * t564 - t284 * t565;
t31 = -t282 * t123 + t131 * t285 - t169 * t398 - t177 * t399;
t550 = -t31 * mrSges(4,1) + t30 * mrSges(4,2);
t438 = Ifges(5,4) * t143;
t60 = -Ifges(5,2) * t142 + Ifges(5,6) * t179 + t438;
t487 = -t60 / 0.2e1;
t139 = Ifges(6,5) * t143;
t57 = Ifges(6,6) * t179 + Ifges(6,3) * t142 + t139;
t549 = t487 + t57 / 0.2e1;
t461 = sin(qJ(1));
t371 = t461 * t286;
t462 = cos(qJ(1));
t374 = t462 * t283;
t232 = t280 * t374 + t371;
t376 = t279 * t462;
t159 = t232 * t285 - t282 * t376;
t372 = t461 * t283;
t373 = t462 * t286;
t231 = -t280 * t373 + t372;
t112 = t159 * t281 - t231 * t284;
t547 = t159 * t284 + t231 * t281;
t100 = -t282 * t169 + t285 * t177;
t82 = pkin(3) * t311 - t100;
t32 = t142 * pkin(4) - t143 * qJ(5) + t82;
t546 = -t82 * mrSges(5,1) - t32 * mrSges(6,1);
t234 = -t280 * t372 + t373;
t375 = t279 * t461;
t163 = t234 * t285 + t282 * t375;
t545 = -g(1) * t163 - g(2) * t159;
t544 = Ifges(4,1) * t481 + Ifges(4,5) * t467;
t493 = Ifges(4,4) * t480 + t544;
t540 = -t311 / 0.2e1;
t446 = mrSges(4,3) - mrSges(3,2);
t538 = mrSges(5,3) + mrSges(6,2);
t532 = -qJ(5) * t561 - qJD(5) * t285 + t520;
t84 = -mrSges(6,2) * t142 + mrSges(6,3) * t179;
t445 = mrSges(5,3) * t142;
t85 = -mrSges(5,2) * t179 - t445;
t529 = -t85 - t84;
t444 = mrSges(5,3) * t143;
t86 = mrSges(5,1) * t179 - t444;
t87 = -mrSges(6,1) * t179 + mrSges(6,2) * t143;
t528 = t87 - t86;
t527 = t311 * Ifges(4,5);
t526 = t311 * Ifges(4,6);
t334 = mrSges(4,1) * t285 - mrSges(4,2) * t282;
t525 = t334 + mrSges(3,1);
t524 = pkin(4) * t561 - t521;
t312 = pkin(4) * t281 - qJ(5) * t284;
t523 = -qJD(5) * t281 + t179 * t312 - t101;
t133 = -t282 * t215 + t216 * t285;
t110 = -pkin(3) * t370 - t133;
t379 = t281 * t413;
t345 = t285 * t379;
t175 = qJD(1) * t345 - t284 * t370;
t407 = t285 * t286;
t187 = (t281 * t283 + t284 * t407) * t279;
t176 = qJD(1) * t187;
t302 = pkin(8) + t312;
t522 = -pkin(4) * t175 + qJ(5) * t176 - t110 + (qJD(4) * t313 - qJD(5) * t284) * t282 + t302 * t398;
t349 = mrSges(3,3) * t370;
t519 = -mrSges(3,1) * t267 - mrSges(4,1) * t184 + mrSges(4,2) * t185 + t349;
t367 = t281 * t398;
t289 = t282 * t394 + t367;
t518 = t175 - t289;
t395 = qJD(4) * t282;
t290 = -t281 * t395 + t284 * t398;
t517 = t176 - t290;
t414 = t279 * t283;
t268 = pkin(7) * t414;
t459 = pkin(1) * t286;
t237 = t280 * t459 - t268;
t330 = mrSges(6,1) * t281 - mrSges(6,3) * t284;
t332 = mrSges(5,1) * t281 + mrSges(5,2) * t284;
t515 = t32 * t330 + t82 * t332;
t514 = t281 * t536 + t284 * t562;
t513 = -t281 * t562 + t284 * t536;
t432 = Ifges(6,5) * t284;
t436 = Ifges(5,4) * t284;
t512 = t281 * t537 - t432 + t436;
t433 = Ifges(6,5) * t281;
t437 = Ifges(5,4) * t281;
t511 = t284 * t537 + t433 - t437;
t33 = -t281 * t83 + t284 * t81;
t510 = -t33 + qJD(5);
t420 = t184 * t284;
t509 = t394 - t420;
t421 = t184 * t281;
t508 = -t396 + t421;
t507 = -t281 * t4 + t284 * t3;
t506 = t1 * t284 + t2 * t281;
t505 = t60 / 0.2e1 - t57 / 0.2e1;
t443 = Ifges(3,4) * t283;
t504 = pkin(1) * (mrSges(3,1) * t283 + mrSges(3,2) * t286) - t283 * (Ifges(3,1) * t286 - t443) / 0.2e1;
t501 = pkin(8) * (-m(4) - t541) - t446;
t201 = t268 + (-pkin(2) - t459) * t280;
t229 = -t280 * t285 + t282 * t414;
t230 = t280 * t282 + t285 * t414;
t351 = -t229 * pkin(3) + t230 * pkin(9);
t106 = t201 - t351;
t202 = pkin(8) * t280 + t405;
t406 = pkin(2) * t413 + pkin(8) * t414;
t460 = pkin(1) * t279;
t203 = -t406 - t460;
t130 = t285 * t202 + t282 * t203;
t108 = -pkin(9) * t413 + t130;
t422 = t281 * t106 + t284 * t108;
t400 = qJD(2) * t279;
t369 = t283 * t400;
t217 = qJD(2) * t299;
t219 = t237 * qJD(2);
t69 = -t202 * t399 + t203 * t398 + t282 * t217 + t285 * t219;
t65 = pkin(9) * t369 + t69;
t368 = t286 * t400;
t154 = qJD(3) * t230 + t282 * t368;
t155 = -qJD(3) * t229 + t285 * t368;
t220 = t405 * qJD(2);
t77 = t154 * pkin(3) - t155 * pkin(9) + t220;
t15 = -qJD(4) * t422 - t281 * t65 + t284 * t77;
t350 = m(6) * pkin(4) + mrSges(5,1) + mrSges(6,1);
t342 = -m(6) * qJ(5) + mrSges(5,2) - mrSges(6,3);
t500 = -g(3) * t230 + t545;
t26 = -pkin(4) * t179 + t510;
t27 = qJ(5) * t179 + t34;
t498 = -t33 * mrSges(5,1) + t26 * mrSges(6,1) + t34 * mrSges(5,2) - t27 * mrSges(6,3);
t472 = -t179 / 0.2e1;
t477 = -t143 / 0.2e1;
t478 = t142 / 0.2e1;
t479 = -t142 / 0.2e1;
t497 = Ifges(5,6) * t478 + Ifges(6,6) * t479 - t472 * t534 + t477 * t536 + t498;
t496 = t279 ^ 2;
t441 = Ifges(4,4) * t185;
t93 = Ifges(4,2) * t184 + t441 - t526;
t484 = -t93 / 0.2e1;
t476 = t143 / 0.2e1;
t471 = t179 / 0.2e1;
t470 = -t184 / 0.2e1;
t469 = -t185 / 0.2e1;
t468 = t185 / 0.2e1;
t465 = t280 / 0.2e1;
t448 = qJD(3) / 0.2e1;
t442 = Ifges(3,4) * t286;
t440 = Ifges(4,4) * t282;
t439 = Ifges(4,4) * t285;
t435 = Ifges(4,5) * t185;
t431 = Ifges(3,6) * t267;
t430 = Ifges(4,6) * t184;
t429 = Ifges(4,3) * t283;
t428 = t184 * mrSges(4,3);
t427 = t185 * mrSges(4,3);
t25 = -pkin(3) * t206 - t31;
t426 = t25 * t282;
t425 = t267 * Ifges(3,5);
t424 = t282 * t31;
t423 = t285 * t30;
t128 = pkin(3) * t185 - pkin(9) * t184;
t53 = t284 * t100 + t281 * t128;
t418 = t231 * t282;
t233 = t280 * t371 + t374;
t416 = t233 * t282;
t415 = t247 * t284;
t412 = t281 * t282;
t411 = t281 * t285;
t410 = t282 * t284;
t409 = t282 * t286;
t408 = t284 * t285;
t196 = pkin(8) * t408 + t281 * t247;
t404 = t462 * pkin(1) + pkin(7) * t375;
t271 = pkin(3) * t413;
t387 = Ifges(4,5) * t103 + Ifges(4,6) * t104 + Ifges(4,3) * t206;
t385 = pkin(8) * t398;
t380 = t279 * t409;
t378 = Ifges(3,5) * t222 - Ifges(3,6) * t221 + Ifges(3,3) * t266;
t377 = t234 * pkin(2) + t404;
t21 = -t99 * mrSges(6,1) + t46 * mrSges(6,2);
t360 = t398 / 0.2e1;
t357 = -t395 / 0.2e1;
t356 = t394 / 0.2e1;
t355 = -t231 * pkin(2) + pkin(8) * t232;
t354 = -t233 * pkin(2) + pkin(8) * t234;
t129 = -t282 * t202 + t203 * t285;
t158 = -t232 * t282 - t285 * t376;
t348 = mrSges(3,3) * t366;
t337 = -pkin(1) * t461 + pkin(7) * t376;
t107 = -t129 + t271;
t335 = mrSges(4,1) * t229 + mrSges(4,2) * t230;
t329 = Ifges(4,1) * t285 - t440;
t324 = t286 * Ifges(3,2) + t443;
t323 = -Ifges(4,2) * t282 + t439;
t322 = -Ifges(5,2) * t281 + t436;
t321 = Ifges(5,2) * t284 + t437;
t318 = Ifges(4,5) * t285 - Ifges(4,6) * t282;
t315 = Ifges(6,3) * t281 + t432;
t314 = -Ifges(6,3) * t284 + t433;
t52 = -t100 * t281 + t128 * t284;
t49 = t106 * t284 - t108 * t281;
t303 = -t232 * pkin(2) + t337;
t70 = -t202 * t398 - t203 * t399 + t217 * t285 - t282 * t219;
t301 = -t541 * pkin(9) + mrSges(4,2) - t538;
t156 = t230 * t281 + t284 * t413;
t14 = t106 * t394 - t108 * t396 + t281 * t77 + t284 * t65;
t295 = t168 * (mrSges(4,1) * t282 + mrSges(4,2) * t285);
t66 = -pkin(3) * t369 - t70;
t264 = Ifges(3,4) * t366;
t243 = -pkin(3) - t313;
t235 = (-mrSges(3,1) * t286 + mrSges(3,2) * t283) * t279;
t214 = -t267 * mrSges(3,2) + t348;
t207 = t302 * t282;
t195 = -pkin(8) * t411 + t415;
t178 = Ifges(4,4) * t184;
t171 = -t415 + (pkin(8) * t281 + pkin(4)) * t285;
t170 = -qJ(5) * t285 + t196;
t167 = Ifges(3,1) * t370 + t264 + t425;
t166 = t324 * t403 + t431;
t162 = t234 * t282 - t285 * t375;
t157 = t230 * t284 - t379;
t147 = -mrSges(4,1) * t311 - t427;
t146 = mrSges(4,2) * t311 + t428;
t117 = t163 * t284 + t233 * t281;
t116 = t163 * t281 - t233 * t284;
t94 = Ifges(4,1) * t185 + t178 - t527;
t92 = -Ifges(4,3) * t311 + t430 + t435;
t80 = -qJD(4) * t156 + t155 * t284 + t281 * t369;
t79 = -qJD(4) * t379 + t155 * t281 + t230 * t394 - t284 * t369;
t76 = -mrSges(4,2) * t206 + mrSges(4,3) * t104;
t75 = mrSges(4,1) * t206 - mrSges(4,3) * t103;
t74 = mrSges(5,1) * t142 + mrSges(5,2) * t143;
t73 = mrSges(6,1) * t142 - mrSges(6,3) * t143;
t72 = pkin(4) * t143 + qJ(5) * t142;
t51 = pkin(4) * t156 - qJ(5) * t157 + t107;
t48 = -mrSges(4,1) * t104 + mrSges(4,2) * t103;
t40 = -pkin(4) * t229 - t49;
t39 = qJ(5) * t229 + t422;
t38 = -pkin(4) * t185 - t52;
t37 = qJ(5) * t185 + t53;
t22 = -mrSges(5,2) * t99 - mrSges(5,3) * t47;
t20 = mrSges(5,1) * t99 - mrSges(5,3) * t46;
t19 = -mrSges(6,2) * t47 + mrSges(6,3) * t99;
t18 = mrSges(5,1) * t47 + mrSges(5,2) * t46;
t17 = mrSges(6,1) * t47 - mrSges(6,3) * t46;
t16 = pkin(4) * t79 - qJ(5) * t80 - qJD(5) * t157 + t66;
t7 = -pkin(4) * t154 - t15;
t6 = qJ(5) * t154 + qJD(5) * t229 + t14;
t5 = pkin(4) * t47 - qJ(5) * t46 - qJD(5) * t143 + t25;
t8 = [t154 * t557 + (Ifges(4,5) * t155 - Ifges(4,6) * t154) * t540 + (t154 * t27 - t32 * t80) * mrSges(6,3) + (mrSges(5,2) * t25 + mrSges(6,2) * t2 - mrSges(5,3) * t4 - mrSges(6,3) * t5 + 0.2e1 * t556) * t157 + (-t154 * t34 + t80 * t82) * mrSges(5,2) + (-t31 * mrSges(4,3) + 0.2e1 * t493) * t230 + t519 * t220 + t124 * t335 + m(5) * (t107 * t25 + t14 * t34 + t15 * t33 + t3 * t422 + t4 * t49 + t66 * t82) + t422 * t22 + (t237 * mrSges(3,1) - t405 * mrSges(3,2) + Ifges(3,3) * t465 + (Ifges(3,5) * t283 + Ifges(3,6) * t286) * t279) * t266 - (t405 * mrSges(3,3) + Ifges(3,6) * t465 + (pkin(1) * mrSges(3,1) + t324) * t279) * t221 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t496 + t144 * t405 + t145 * t237 - t215 * t220 + t218 * t219) + ((-t215 * mrSges(3,3) + t167 / 0.2e1 + t425 / 0.2e1) * t286 + (-t166 / 0.2e1 + t92 / 0.2e1 - t218 * mrSges(3,3) - t101 * mrSges(4,2) + t100 * mrSges(4,1) + t430 / 0.2e1 + t435 / 0.2e1 + Ifges(4,3) * t448 - t431 / 0.2e1) * t283 + (t286 * (-Ifges(3,2) * t283 + t442) / 0.2e1 - t286 * t429 / 0.2e1 - t504) * t403) * t400 + t378 * t465 + (Ifges(4,1) * t155 - Ifges(4,4) * t154) * t468 + m(6) * (t1 * t39 + t16 * t32 + t2 * t40 + t26 * t7 + t27 * t6 + t5 * t51) + m(4) * (t100 * t70 + t101 * t69 + t124 * t201 + t129 * t31 + t130 * t30 + t168 * t220) + (Ifges(6,5) * t80 + Ifges(6,6) * t154) * t478 + (-t30 * mrSges(4,3) - t567 + t539 / 0.2e1 + t566) * t229 + (-t237 * mrSges(3,3) + Ifges(3,5) * t465 + (-pkin(1) * mrSges(3,2) + Ifges(3,1) * t283 + t442) * t279) * t222 + t154 * t484 + (-t235 * t460 + Ifges(2,3)) * qJDD(1) + t184 * (Ifges(4,4) * t155 - Ifges(4,2) * t154) / 0.2e1 + (-t100 * t155 - t101 * t154) * mrSges(4,3) + t145 * (mrSges(3,1) * t280 - mrSges(3,3) * t414) + (Ifges(5,4) * t80 + Ifges(5,6) * t154) * t479 + (-m(3) * t404 - t234 * mrSges(3,1) - mrSges(3,3) * t375 - t462 * mrSges(2,1) + t461 * mrSges(2,2) - m(4) * t377 - t163 * mrSges(4,1) + t501 * t233 - t350 * t117 + t342 * t116 + t301 * t162 - t541 * (t163 * pkin(3) + t377)) * g(2) + (-t154 * t534 + t536 * t80 - t562 * t79) * t471 + (mrSges(5,1) * t25 + mrSges(6,1) * t5 - mrSges(6,2) * t1 - mrSges(5,3) * t3 - Ifges(5,2) * t491 + Ifges(6,3) * t490 - t483 * t562 + t492 * t554 - t559) * t156 + (-m(4) * t303 + t159 * mrSges(4,1) - m(3) * t337 + t232 * mrSges(3,1) - mrSges(3,3) * t376 + t461 * mrSges(2,1) + t462 * mrSges(2,2) - t501 * t231 + t350 * t547 - t342 * t112 + t301 * t158 + t541 * (pkin(3) * t159 - t303)) * g(1) + (-mrSges(6,2) * t27 - mrSges(5,3) * t34 - Ifges(5,2) * t479 + Ifges(6,3) * t478 - t546 + t549) * t79 + (-t387 / 0.2e1 - Ifges(4,3) * t467 - Ifges(4,6) * t480 - Ifges(4,5) * t481 + t144 * mrSges(3,3) + t550) * t413 - t280 * t553 + t39 * t19 + t40 * t21 + t530 * t80 / 0.2e1 + (t154 * t536 + t537 * t80 + t554 * t79) * t476 + t49 * t20 + t51 * t17 + t16 * t73 + t66 * t74 + t6 * t84 + t14 * t85 + t15 * t86 + t7 * t87 + t107 * t18 + t129 * t75 + t130 * t76 + t69 * t146 + t70 * t147 + t26 * (-mrSges(6,1) * t154 + mrSges(6,2) * t80) + t33 * (mrSges(5,1) * t154 - mrSges(5,3) * t80) + t155 * t94 / 0.2e1 + t168 * (mrSges(4,1) * t154 + mrSges(4,2) * t155) + t201 * t48 + t219 * t214; t410 * t556 + (-t100 * t133 - t101 * t134 - pkin(2) * t124 + (-t424 + t423 + (-t100 * t285 - t101 * t282) * qJD(3)) * pkin(8)) * m(4) + (t184 * t323 + t185 * t329) * t448 + (t315 * t490 + t322 * t491 + t330 * t5 + t356 * t57 + t483 * t513 + t492 * t511 + t493 + t544) * t282 + (t281 * t57 + t94) * t360 + (-t385 - t133) * t147 + (-t3 * t412 + t33 * t517 + t34 * t518 - t4 * t410) * mrSges(5,3) + (-t1 * t412 + t2 * t410 - t26 * t517 + t27 * t518) * mrSges(6,2) + (-m(4) * t168 + t349 - t519) * t218 + mrSges(4,3) * t423 + t332 * t426 - t124 * t334 + (t385 - t110) * t74 - t553 + (t311 * (t286 * t318 + t429) + t283 * t166) * t403 / 0.2e1 + t520 * t85 + t521 * t86 + (-t110 * t82 + t195 * t4 + t196 * t3 + (t398 * t82 + t426) * pkin(8) + t520 * t34 + t521 * t33) * m(5) + t522 * t73 + t524 * t87 + (-t314 * t478 - t321 * t479 - t471 * t514 - t476 * t512) * t395 + t532 * t84 + (t1 * t170 + t171 * t2 + t207 * t5 + t26 * t524 + t27 * t532 + t32 * t522) * m(6) + (-t82 * mrSges(5,2) + t32 * mrSges(6,3) + Ifges(5,4) * t478 + Ifges(6,5) * t479 + t472 * t536 + t477 * t537) * t176 + t82 * (mrSges(5,1) * t289 + mrSges(5,2) * t290) + t32 * (mrSges(6,1) * t289 - mrSges(6,3) * t290) + t504 * qJD(1) ^ 2 * t496 + (-(mrSges(4,1) * t283 - mrSges(4,3) * t407) * t403 - mrSges(4,3) * t398) * t100 + t440 * t480 + (t348 - t214) * t215 + t367 * t487 + t378 - (t185 * (Ifges(4,5) * t283 + t286 * t329) + t184 * (Ifges(4,6) * t283 + t286 * t323) + t267 * (Ifges(3,5) * t286 - Ifges(3,6) * t283) + t283 * t92 + (-Ifges(3,2) * t370 + t282 * t531 + t285 * t94 + t167 + t264) * t286) * t403 / 0.2e1 - t101 * (-mrSges(4,2) * t283 - mrSges(4,3) * t409) * t403 + (-m(4) * t406 - (t283 * mrSges(4,3) + t286 * t334) * t279 + t235 - t538 * t380 - t541 * (pkin(9) * t380 + t285 * t271 + t406) - t350 * t187 + t342 * (-t284 * t414 + t345)) * g(3) + (-m(4) * t354 + t538 * t416 - t541 * (-pkin(9) * t416 - t233 * t458 + t354) - t446 * t234 + t525 * t233 - t350 * (-t233 * t408 + t234 * t281) + t342 * (-t233 * t411 - t234 * t284)) * g(1) + (-m(4) * t355 + t538 * t418 - t541 * (-pkin(9) * t418 - t231 * t458 + t355) - t446 * t232 + t525 * t231 - t350 * (-t231 * t408 + t232 * t281) + t342 * (-t231 * t411 - t232 * t284)) * g(2) - t539 * t285 / 0.2e1 + (t295 + t318 * t540 + (Ifges(6,6) * t282 + t285 * t315) * t478 + (Ifges(5,6) * t282 + t285 * t322) * t479 + (t282 * t536 + t285 * t511) * t476 + (-t282 * t534 + t285 * t513) * t471) * qJD(3) + t284 * t60 * t357 - mrSges(4,3) * t424 + (pkin(8) * t76 - t566) * t285 + (t93 / 0.2e1 + t497) * t341 + t439 * t481 + (-Ifges(5,2) * t478 + Ifges(6,3) * t479 - t472 * t562 + t477 * t554 + t505 + t546) * t175 - t559 * t412 - t295 * t366 + (-t75 + t18) * pkin(8) * t282 + t530 * (t281 * t357 + t284 * t360 - t176 / 0.2e1) - pkin(2) * t48 + (-t101 * mrSges(4,3) - pkin(8) * t146 + t484 - t498 + t557) * t399 + t145 * mrSges(3,1) - t134 * t146 + t170 * t19 + t171 * t21 + t195 * t20 + t196 * t22 + t207 * t17; t281 * t556 + (-t322 / 0.2e1 + t315 / 0.2e1) * t397 + t387 + t5 * t331 - t25 * t333 + (((t26 * t284 - t27 * t281) * qJD(4) + t506) * m(6) + ((-t281 * t34 - t284 * t33) * qJD(4) + t507) * m(5) + t528 * t394 + t529 * t396 + (t22 + t19) * t284 + (-t20 + t21) * t281 + t545 * t541) * pkin(9) + (mrSges(4,2) * t163 + t162 * t560) * g(1) + (mrSges(4,2) * t159 - t158 * t560) * g(2) + t523 * t73 + (t143 * t511 + t179 * t513) * qJD(4) / 0.2e1 + t512 * t492 + t514 * t483 + t515 * qJD(4) + (t26 * t509 + t27 * t508 + t500 + t506) * mrSges(6,2) + (-t33 * t509 + t34 * t508 + t500 + t507) * mrSges(5,3) + t505 * t421 + (-t441 + t531) * t469 + t93 * t468 + (t356 - t420 / 0.2e1) * t530 - t550 + t314 * t490 + t321 * t491 + (-m(5) * t82 + t147 + t427 - t74) * t101 + (-Ifges(4,2) * t470 - t526 / 0.2e1 - t168 * mrSges(4,1) + t497) * t185 + (Ifges(4,1) * t469 + t322 * t478 + t315 * t479 + t527 / 0.2e1 - t168 * mrSges(4,2) + t511 * t477 + t513 * t472 - t515) * t184 + (t428 - t146) * t100 + (t94 + t178) * t470 + (t243 * t5 - t26 * t38 - t27 * t37 + t523 * t32) * m(6) + (-pkin(3) * t25 - t33 * t52 - t34 * t53) * m(5) + (t516 * t229 - t351 * t541 + t335) * g(3) + t559 * t284 + t549 * t396 - pkin(3) * t18 - t37 * t84 - t53 * t85 - t52 * t86 - t38 * t87 + t243 * t17; (t142 * t26 + t143 * t27) * mrSges(6,2) + t543 + (-t445 + t529) * t33 + (t350 * t112 + t342 * t547) * g(2) + (t444 - t528) * t34 + (-pkin(4) * t2 + qJ(5) * t1 - t26 * t34 + t27 * t510 - t32 * t72) * m(6) + t60 * t476 + (Ifges(6,3) * t143 - t434) * t479 + (-t142 * t537 + t139 - t438 + t57) * t477 + (t156 * t350 + t157 * t342) * g(3) + (t116 * t350 + t117 * t342) * g(1) + (-t142 * t536 - t143 * t562) * t472 + (-Ifges(5,2) * t143 - t140 + t530) * t478 + qJ(5) * t19 - pkin(4) * t21 - t72 * t73 + qJD(5) * t84 - t32 * (mrSges(6,1) * t143 + mrSges(6,3) * t142) - t82 * (mrSges(5,1) * t143 - mrSges(5,2) * t142) + t539; t143 * t73 - t179 * t84 + (-g(1) * t116 - g(2) * t112 - g(3) * t156 + t32 * t143 - t27 * t179 + t2) * m(6) + t21;];
tau = t8;
