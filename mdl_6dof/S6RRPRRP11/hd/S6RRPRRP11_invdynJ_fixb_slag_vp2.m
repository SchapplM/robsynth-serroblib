% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP11
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
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP11_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP11_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:45:27
% EndTime: 2019-03-09 12:46:19
% DurationCPUTime: 35.37s
% Computational Cost: add. (10107->790), mult. (20764->1021), div. (0->0), fcn. (13142->10), ass. (0->348)
t292 = sin(qJ(4));
t296 = cos(qJ(4));
t297 = cos(qJ(2));
t394 = qJD(1) * t297;
t214 = -qJD(2) * t292 - t296 * t394;
t291 = sin(qJ(5));
t295 = cos(qJ(5));
t392 = qJD(2) * t296;
t316 = t292 * t394 - t392;
t134 = t214 * t291 - t295 * t316;
t269 = pkin(7) * t394;
t224 = pkin(3) * t394 + t269;
t289 = qJD(2) * qJ(3);
t190 = t289 + t224;
t140 = -pkin(4) * t214 + t190;
t293 = sin(qJ(2));
t384 = qJD(1) * qJD(2);
t228 = qJDD(1) * t293 + t297 * t384;
t213 = qJDD(4) + t228;
t201 = qJDD(5) + t213;
t227 = -t297 * qJDD(1) + t293 * t384;
t124 = qJD(4) * t214 + qJDD(2) * t296 + t227 * t292;
t125 = qJD(4) * t316 - qJDD(2) * t292 + t227 * t296;
t352 = t295 * t214 + t291 * t316;
t45 = qJD(5) * t352 + t124 * t295 + t125 * t291;
t210 = t228 * pkin(7);
t358 = qJDD(3) + t210;
t494 = pkin(2) + pkin(8);
t127 = pkin(3) * t228 - qJDD(2) * t494 + t358;
t274 = t293 * qJD(1);
t426 = qJDD(1) * pkin(1);
t309 = -qJ(3) * t228 - qJD(3) * t274 - t426;
t94 = t227 * t494 + t309;
t277 = t293 * qJ(3);
t357 = -pkin(1) - t277;
t151 = (-t297 * t494 + t357) * qJD(1);
t267 = pkin(7) * t274;
t223 = -pkin(3) * t274 - t267;
t515 = qJD(3) - t223;
t165 = -qJD(2) * t494 + t515;
t98 = t151 * t296 + t165 * t292;
t29 = -qJD(4) * t98 + t296 * t127 - t292 * t94;
t19 = pkin(4) * t213 - pkin(9) * t124 + t29;
t389 = qJD(4) * t296;
t390 = qJD(4) * t292;
t28 = t292 * t127 - t151 * t390 + t165 * t389 + t296 * t94;
t23 = pkin(9) * t125 + t28;
t255 = t274 + qJD(4);
t97 = -t151 * t292 + t296 * t165;
t84 = pkin(9) * t316 + t97;
t75 = pkin(4) * t255 + t84;
t85 = pkin(9) * t214 + t98;
t83 = t295 * t85;
t31 = t291 * t75 + t83;
t6 = -qJD(5) * t31 + t295 * t19 - t23 * t291;
t2 = pkin(5) * t201 - qJ(6) * t45 - qJD(6) * t134 + t6;
t46 = -qJD(5) * t134 - t124 * t291 + t125 * t295;
t385 = qJD(5) * t295;
t386 = qJD(5) * t291;
t5 = t291 * t19 + t295 * t23 + t75 * t385 - t386 * t85;
t3 = qJ(6) * t46 + qJD(6) * t352 + t5;
t379 = qJD(4) + qJD(5);
t246 = t274 + t379;
t472 = -t246 / 0.2e1;
t486 = -t134 / 0.2e1;
t489 = -t352 / 0.2e1;
t545 = Ifges(6,3) + Ifges(7,3);
t546 = Ifges(6,6) + Ifges(7,6);
t548 = Ifges(6,5) + Ifges(7,5);
t523 = t201 * t545 + t45 * t548 + t46 * t546;
t552 = Ifges(6,1) + Ifges(7,1);
t550 = Ifges(6,4) + Ifges(7,4);
t562 = t352 * t550;
t536 = t134 * t552 + t246 * t548 + t562;
t547 = Ifges(6,2) + Ifges(7,2);
t561 = t134 * t550;
t78 = -pkin(5) * t352 + qJD(6) + t140;
t575 = t6 * mrSges(6,1) + t2 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2) + t523 + (-t134 * t546 + t548 * t352) * t472 - t140 * (mrSges(6,1) * t134 + mrSges(6,2) * t352) - t78 * (mrSges(7,1) * t134 + mrSges(7,2) * t352) + (t552 * t352 - t561) * t486 + (-t134 * t547 + t536 + t562) * t489;
t554 = -mrSges(7,1) - mrSges(6,1);
t553 = -mrSges(7,2) - mrSges(6,2);
t574 = -mrSges(6,3) - mrSges(7,3);
t268 = pkin(2) * t274;
t427 = qJ(3) * t297;
t329 = pkin(8) * t293 - t427;
t171 = qJD(1) * t329 + t268;
t121 = -t171 * t292 + t296 * t224;
t414 = t292 * t293;
t323 = pkin(4) * t297 - pkin(9) * t414;
t455 = pkin(9) + t494;
t573 = -qJD(1) * t323 + t455 * t390 - t121;
t122 = t296 * t171 + t292 * t224;
t233 = t455 * t296;
t367 = t296 * t274;
t572 = pkin(9) * t367 + qJD(4) * t233 + t122;
t409 = t295 * t296;
t137 = -t291 * t390 - t292 * t386 + t379 * t409;
t417 = t291 * t292;
t160 = -t274 * t417 + t295 * t367;
t398 = t160 + t137;
t325 = t291 * t296 + t295 * t292;
t138 = t379 * t325;
t315 = t325 * t293;
t161 = qJD(1) * t315;
t397 = -t161 - t138;
t551 = -Ifges(4,4) + Ifges(3,5);
t549 = Ifges(4,5) - Ifges(3,6);
t341 = t297 * mrSges(4,2) - t293 * mrSges(4,3);
t345 = mrSges(3,1) * t297 - mrSges(3,2) * t293;
t571 = t341 - t345;
t570 = qJD(3) + t267;
t290 = qJ(4) + qJ(5);
t275 = sin(t290);
t467 = pkin(4) * t292;
t229 = pkin(5) * t275 + t467;
t276 = cos(t290);
t342 = mrSges(5,1) * t292 + mrSges(5,2) * t296;
t569 = -m(6) * t467 - m(7) * t229 + t275 * t554 + t553 * t276 - t342;
t299 = -pkin(9) - pkin(8);
t285 = -qJ(6) + t299;
t568 = -m(6) * (-pkin(2) + t299) + m(5) * t494 + mrSges(5,3) - m(7) * (-pkin(2) + t285) - t574;
t81 = t291 * t85;
t30 = t295 * t75 - t81;
t563 = qJ(6) * t134;
t24 = t30 - t563;
t21 = pkin(5) * t246 + t24;
t532 = qJ(6) * t352;
t25 = t31 + t532;
t443 = Ifges(4,6) * t297;
t331 = -t293 * Ifges(4,2) - t443;
t567 = t25 * mrSges(7,2) + t31 * mrSges(6,2) + Ifges(4,4) * qJD(2) / 0.2e1 + qJD(1) * t331 / 0.2e1 - t21 * mrSges(7,1) - t30 * mrSges(6,1);
t537 = t246 * t546 + t352 * t547 + t561;
t566 = t537 / 0.2e1;
t565 = -t550 * t46 / 0.2e1 - t552 * t45 / 0.2e1 - t548 * t201 / 0.2e1;
t294 = sin(qJ(1));
t564 = g(2) * t294;
t232 = t455 * t292;
t142 = -t295 * t232 - t291 * t233;
t539 = -qJD(5) * t142 + t291 * t572 + t295 * t573;
t538 = t232 * t386 - t233 * t385 + t291 * t573 - t295 * t572;
t279 = t296 * pkin(4);
t263 = t279 + pkin(3);
t526 = pkin(4) * t389 + t263 * t274 + t570;
t324 = -t409 + t417;
t560 = -t30 * t397 - t31 * t398 + t324 * t6 - t325 * t5;
t559 = t2 * t324 - t21 * t397 - t25 * t398 - t3 * t325;
t558 = t297 * t574 + t571;
t454 = mrSges(5,3) * t214;
t146 = -mrSges(5,2) * t255 + t454;
t453 = mrSges(5,3) * t316;
t147 = mrSges(5,1) * t255 + t453;
t327 = t296 * t146 - t292 * t147;
t91 = mrSges(5,1) * t213 - mrSges(5,3) * t124;
t92 = -mrSges(5,2) * t213 + mrSges(5,3) * t125;
t556 = t327 * qJD(4) + t292 * t92 + t296 * t91;
t544 = t201 * t546 + t45 * t550 + t46 * t547;
t34 = -mrSges(7,2) * t201 + mrSges(7,3) * t46;
t35 = -mrSges(6,2) * t201 + mrSges(6,3) * t46;
t542 = t34 + t35;
t541 = -pkin(5) * t394 - qJ(6) * t397 + qJD(6) * t324 + t539;
t540 = -qJ(6) * t398 - qJD(6) * t325 + t538;
t501 = m(6) * pkin(4);
t533 = -t501 - mrSges(5,1);
t531 = pkin(5) * t398 + t526;
t493 = pkin(3) + pkin(7);
t242 = t493 * t293;
t219 = t296 * t242;
t282 = t297 * pkin(2);
t396 = t282 + t277;
t348 = pkin(8) * t297 + t396;
t203 = -pkin(1) - t348;
t356 = pkin(9) * t297 - t203;
t110 = pkin(4) * t293 + t292 * t356 + t219;
t218 = t292 * t242;
t136 = t296 * t203 + t218;
t406 = t296 * t297;
t118 = -pkin(9) * t406 + t136;
t62 = t291 * t110 + t295 * t118;
t139 = -mrSges(5,1) * t214 - mrSges(5,2) * t316;
t370 = mrSges(4,1) * t394;
t239 = -qJD(2) * mrSges(4,3) - t370;
t530 = -t239 + t139;
t529 = qJD(2) * mrSges(3,2) - mrSges(3,3) * t394 + t239;
t371 = mrSges(4,1) * t274;
t528 = mrSges(3,3) * t274 + t371 + (-mrSges(3,1) + mrSges(4,2)) * qJD(2);
t527 = t553 * t275 * t297;
t444 = Ifges(4,6) * t293;
t525 = t293 * (-Ifges(4,2) * t297 + t444) + t297 * (Ifges(4,3) * t293 - t443);
t524 = t293 * t549 + t297 * t551;
t298 = cos(qJ(1));
t412 = t293 * t294;
t168 = t275 * t298 + t276 * t412;
t169 = -t275 * t412 + t276 * t298;
t522 = t168 * t554 + t553 * t169;
t411 = t293 * t298;
t166 = -t275 * t294 + t276 * t411;
t167 = t275 * t411 + t276 * t294;
t521 = t166 * t554 - t553 * t167;
t209 = t227 * pkin(7);
t520 = -t209 * t297 + t210 * t293;
t162 = -qJDD(2) * qJ(3) - qJD(2) * qJD(3) + t209;
t170 = -qJDD(2) * pkin(2) + t358;
t519 = -t162 * t297 + t170 * t293;
t517 = t28 * t292 + t29 * t296;
t516 = g(1) * t298 + t564;
t372 = m(4) + m(5) + m(6) + m(7);
t436 = t316 * Ifges(5,4);
t114 = t214 * Ifges(5,2) + t255 * Ifges(5,6) - t436;
t204 = Ifges(5,4) * t214;
t115 = -Ifges(5,1) * t316 + t255 * Ifges(5,5) + t204;
t330 = -t297 * Ifges(4,3) - t444;
t513 = Ifges(4,5) * qJD(2) + qJD(1) * t330 + t296 * t114 + t292 * t115;
t266 = Ifges(3,4) * t394;
t510 = Ifges(3,1) * t274 + Ifges(3,5) * qJD(2) - Ifges(5,5) * t316 + t214 * Ifges(5,6) + t255 * Ifges(5,3) + t134 * t548 + t246 * t545 + t352 * t546 + t266;
t509 = t297 * t379;
t431 = t297 * mrSges(4,3);
t508 = -t431 + t569 * t297 + (m(4) * pkin(2) - mrSges(4,2) + t568) * t293;
t230 = pkin(5) * t276 + t279;
t503 = -m(5) * pkin(3) - m(6) * t263 - m(7) * (pkin(3) + t230) - mrSges(4,1) + mrSges(2,2) - mrSges(3,3);
t500 = m(7) * pkin(5);
t499 = t45 / 0.2e1;
t498 = t46 / 0.2e1;
t497 = -t124 * Ifges(5,4) / 0.2e1 - t125 * Ifges(5,2) / 0.2e1 - t213 * Ifges(5,6) / 0.2e1;
t492 = t124 / 0.2e1;
t491 = t125 / 0.2e1;
t488 = t352 / 0.2e1;
t485 = t134 / 0.2e1;
t478 = t201 / 0.2e1;
t477 = t213 / 0.2e1;
t475 = -t316 / 0.2e1;
t471 = t246 / 0.2e1;
t469 = pkin(4) * t316;
t466 = pkin(4) * t295;
t464 = pkin(7) * t293;
t280 = t297 * pkin(7);
t38 = t295 * t84 - t81;
t452 = mrSges(6,3) * t352;
t451 = mrSges(6,3) * t134;
t450 = mrSges(7,3) * t352;
t449 = mrSges(7,3) * t134;
t448 = Ifges(3,4) * t293;
t447 = Ifges(3,4) * t297;
t446 = Ifges(5,4) * t292;
t445 = Ifges(5,4) * t296;
t418 = t285 * t297;
t413 = t292 * t297;
t410 = t294 * t296;
t405 = t296 * t298;
t404 = t297 * t298;
t403 = t297 * t299;
t243 = t297 * pkin(3) + t280;
t395 = t298 * pkin(1) + t294 * pkin(7);
t393 = qJD(2) * t293;
t391 = qJD(2) * t297;
t388 = qJD(4) * t297;
t376 = pkin(4) * t386;
t375 = pkin(4) * t385;
t369 = Ifges(5,5) * t124 + Ifges(5,6) * t125 + Ifges(5,3) * t213;
t184 = pkin(4) * t406 + t243;
t366 = t292 * t388;
t15 = -t46 * mrSges(7,1) + t45 * mrSges(7,2);
t360 = -t389 / 0.2e1;
t258 = qJ(3) + t467;
t37 = -t291 * t84 - t83;
t355 = -t384 / 0.2e1;
t176 = t228 * mrSges(4,1) + qJDD(2) * mrSges(4,2);
t61 = t295 * t110 - t118 * t291;
t351 = pkin(2) * t393 - qJD(3) * t293;
t148 = qJD(2) * t329 + t351;
t226 = t493 * t391;
t353 = -t148 * t292 + t296 * t226;
t141 = t232 * t291 - t295 * t233;
t344 = mrSges(3,1) * t293 + mrSges(3,2) * t297;
t343 = mrSges(5,1) * t296 - mrSges(5,2) * t292;
t340 = Ifges(5,1) * t296 - t446;
t339 = Ifges(5,1) * t292 + t445;
t338 = t297 * Ifges(3,2) + t448;
t336 = -Ifges(5,2) * t292 + t445;
t335 = Ifges(5,2) * t296 + t446;
t333 = Ifges(5,5) * t296 - Ifges(5,6) * t292;
t332 = Ifges(5,5) * t292 + Ifges(5,6) * t296;
t328 = t292 * t97 - t296 * t98;
t231 = -qJD(2) * pkin(2) + t570;
t237 = -t269 - t289;
t326 = t231 * t297 + t237 * t293;
t322 = t357 - t282;
t321 = pkin(1) * t344;
t192 = -t292 * t294 + t293 * t405;
t194 = t292 * t298 + t293 * t410;
t55 = t323 * qJD(2) + (t296 * t356 - t218) * qJD(4) + t353;
t312 = t293 * t392 + t366;
t70 = t296 * t148 - t203 * t390 + t292 * t226 + t242 * t389;
t57 = pkin(9) * t312 + t70;
t9 = t110 * t385 - t118 * t386 + t291 * t55 + t295 * t57;
t191 = t322 * qJD(1);
t320 = t191 * (-mrSges(4,2) * t293 - t431);
t319 = t293 * (Ifges(3,1) * t297 - t448);
t311 = t292 * t393 - t296 * t388;
t308 = Ifges(5,5) * t297 + t293 * t339;
t307 = Ifges(5,6) * t297 + t293 * t335;
t306 = Ifges(5,3) * t297 + t293 * t332;
t130 = -pkin(3) * t227 - t162;
t10 = -qJD(5) * t62 - t291 * t57 + t295 * t55;
t143 = -pkin(4) * t366 + (-pkin(7) - t263) * t393;
t304 = -qJD(4) * t328 + t517;
t74 = -pkin(4) * t125 + t130;
t262 = pkin(5) + t466;
t234 = -pkin(1) - t396;
t225 = t493 * t393;
t221 = -qJ(3) * t394 + t268;
t220 = t341 * qJD(1);
t202 = t343 * t297;
t195 = -t292 * t412 + t405;
t193 = t292 * t411 + t410;
t185 = Ifges(3,6) * qJD(2) + qJD(1) * t338;
t177 = -qJ(3) * t391 + t351;
t175 = mrSges(4,1) * t227 - qJDD(2) * mrSges(4,3);
t173 = t325 * t297;
t172 = t324 * t297;
t163 = pkin(5) * t325 + t258;
t135 = -t203 * t292 + t219;
t126 = -pkin(5) * t172 + t184;
t123 = pkin(2) * t227 + t309;
t108 = -qJ(6) * t325 + t142;
t107 = qJ(6) * t324 + t141;
t104 = mrSges(6,1) * t246 - t451;
t103 = mrSges(7,1) * t246 - t449;
t102 = -mrSges(6,2) * t246 + t452;
t101 = -mrSges(7,2) * t246 + t450;
t95 = pkin(5) * t134 - t469;
t87 = -t324 * t393 + t325 * t509;
t86 = qJD(2) * t315 + t324 * t509;
t77 = -mrSges(6,1) * t352 + mrSges(6,2) * t134;
t76 = -mrSges(7,1) * t352 + mrSges(7,2) * t134;
t71 = -qJD(4) * t136 + t353;
t69 = -mrSges(5,1) * t125 + mrSges(5,2) * t124;
t60 = t124 * Ifges(5,1) + t125 * Ifges(5,4) + t213 * Ifges(5,5);
t58 = -pkin(5) * t87 + t143;
t50 = qJ(6) * t172 + t62;
t49 = pkin(5) * t293 + qJ(6) * t173 + t61;
t33 = mrSges(6,1) * t201 - mrSges(6,3) * t45;
t32 = mrSges(7,1) * t201 - mrSges(7,3) * t45;
t27 = t38 - t563;
t26 = t37 - t532;
t20 = -pkin(5) * t46 + qJDD(6) + t74;
t16 = -mrSges(6,1) * t46 + mrSges(6,2) * t45;
t8 = qJ(6) * t87 + qJD(6) * t172 + t9;
t7 = pkin(5) * t391 - qJ(6) * t86 + qJD(6) * t173 + t10;
t1 = [t544 * t172 / 0.2e1 + t114 * t366 / 0.2e1 + (t545 * t471 + t546 * t488 - t98 * mrSges(5,2) + t97 * mrSges(5,1) + t548 * t485 + t510 / 0.2e1 + t528 * pkin(7) + t231 * mrSges(4,1) - t567) * t391 + (t172 * t546 - t173 * t548 + t293 * t545) * t478 + (t172 * t547 - t173 * t550 + t293 * t546) * t498 + (t293 * t551 - t297 * t549) * qJDD(2) / 0.2e1 + t29 * (mrSges(5,1) * t293 + mrSges(5,3) * t413) + t173 * t565 + m(5) * (t130 * t243 + t135 * t29 + t136 * t28 - t190 * t225 + t70 * t98 + t71 * t97) + t228 * t447 / 0.2e1 + (-Ifges(3,4) * t227 + Ifges(3,5) * qJDD(2) + t369 + t523) * t293 / 0.2e1 + t6 * (mrSges(6,1) * t293 + mrSges(6,3) * t173) + t2 * (mrSges(7,1) * t293 + mrSges(7,3) * t173) + t74 * (-mrSges(6,1) * t172 - mrSges(6,2) * t173) + t20 * (-mrSges(7,1) * t172 - mrSges(7,2) * t173) + (t172 * t550 - t173 * t552 + t293 * t548) * t499 + t524 * qJD(2) ^ 2 / 0.2e1 + t525 * t355 + t98 * mrSges(5,3) * t312 + (qJD(2) * t308 - t340 * t388) * t475 + (Ifges(5,3) * t293 - t297 * t332) * t477 + (-qJDD(2) * mrSges(3,1) + t176) * t464 + (t297 * (-Ifges(3,2) * t293 + t447) + t319) * t384 / 0.2e1 + (t548 * t471 + t550 * t488 + t536 / 0.2e1 - t21 * mrSges(7,3) - t30 * mrSges(6,3) + t552 * t485 + mrSges(6,2) * t140 + mrSges(7,2) * t78) * t86 + t28 * (-mrSges(5,2) * t293 - mrSges(5,3) * t406) + t345 * t426 + t519 * mrSges(4,1) - t321 * t384 + t406 * t497 + qJD(2) * t320 + (-t195 * mrSges(5,1) + t194 * mrSges(5,2) + t554 * t169 - t553 * t168 + (-m(6) * (-t258 * t293 - pkin(1)) - m(5) * t357 - m(7) * (-pkin(1) + (-qJ(3) - t229) * t293) + m(3) * pkin(1) - m(4) * t322 + mrSges(2,1) + t568 * t297 - t571) * t294 + ((-m(3) - t372) * pkin(7) + t503) * t298) * g(1) - t60 * t413 / 0.2e1 + (-t185 / 0.2e1 + t513 / 0.2e1 + t529 * pkin(7) + t237 * mrSges(4,1)) * t393 + t214 * (qJD(2) * t307 - t336 * t388) / 0.2e1 + t255 * (qJD(2) * t306 - t333 * t388) / 0.2e1 + (Ifges(5,6) * t293 - t297 * t335) * t491 + (Ifges(5,5) * t293 - t297 * t339) * t492 + (-m(3) * t395 - t193 * mrSges(5,1) - t192 * mrSges(5,2) + (-m(5) * pkin(8) - mrSges(5,3)) * t404 + t554 * t167 + t553 * t166 - t372 * (pkin(2) * t404 + qJ(3) * t411 + t395) + t503 * t294 + (-m(7) * (t229 * t293 - t418) - m(6) * (pkin(4) * t414 - t403) - mrSges(2,1) + t558) * t298) * g(2) + (-qJDD(2) * mrSges(3,2) - t175) * t280 + m(7) * (t126 * t20 + t2 * t49 + t21 * t7 + t25 * t8 + t3 * t50 + t58 * t78) + m(6) * (t10 * t30 + t140 * t143 + t184 * t74 + t31 * t9 + t5 * t62 + t6 * t61) + (-t140 * mrSges(6,1) - t78 * mrSges(7,1) + t31 * mrSges(6,3) + t25 * mrSges(7,3) + t546 * t471 + t550 * t485 + t547 * t488 + t566) * t87 + t228 * t293 * Ifges(3,1) - t227 * t338 / 0.2e1 + t123 * t341 + t227 * t330 / 0.2e1 - t228 * t331 / 0.2e1 - t97 * mrSges(5,3) * t311 + m(4) * (t123 * t234 + t177 * t191 + (qJD(2) * t326 + t519) * pkin(7)) + (-t227 * t280 + t228 * t464 + t520) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t520) + Ifges(2,3) * qJDD(1) + t190 * (-mrSges(5,1) * t312 + mrSges(5,2) * t311) + t297 * (Ifges(3,4) * t228 - Ifges(3,2) * t227 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t297 * (Ifges(4,5) * qJDD(2) - Ifges(4,6) * t228 + Ifges(4,3) * t227) / 0.2e1 - t293 * (Ifges(4,4) * qJDD(2) - Ifges(4,2) * t228 + Ifges(4,6) * t227) / 0.2e1 + t5 * (-mrSges(6,2) * t293 + mrSges(6,3) * t172) + t3 * (-mrSges(7,2) * t293 + mrSges(7,3) * t172) + t243 * t69 + t234 * (-mrSges(4,2) * t227 - mrSges(4,3) * t228) - pkin(1) * (mrSges(3,1) * t227 + mrSges(3,2) * t228) + t177 * t220 - t225 * t139 + t130 * t202 + t184 * t16 + t143 * t77 + t70 * t146 + t71 * t147 + t297 * t115 * t360 + t49 * t32 + t50 * t34 + t61 * t33 + t62 * t35 + t58 * t76 + t8 * t101 + t9 * t102 + t7 * t103 + t10 * t104 + t126 * t15 + t135 * t91 + t136 * t92; (-t175 + t69) * qJ(3) + (-t160 / 0.2e1 - t137 / 0.2e1) * t537 + (t545 * t472 + t548 * t486 + t546 * t489 + t567) * t394 + (t160 * t546 + t161 * t548) * t472 + t549 * t227 + (t160 * t547 + t161 * t550) * t489 + t551 * t228 + (t130 * qJ(3) - t121 * t97 - t122 * t98 + t190 * t515) * m(5) + t538 * t102 + (t526 * t140 + t141 * t6 + t142 * t5 + t258 * t74 + t539 * t30 + t538 * t31) * m(6) + t539 * t104 + t540 * t101 + t541 * t103 + (t107 * t2 + t108 * t3 + t163 * t20 + t541 * t21 + t540 * t25 + t531 * t78) * m(7) + (-t372 * t427 + t508) * t564 + t324 * t565 - t115 * t390 / 0.2e1 + (-t319 / 0.2e1 + t321 + t525 / 0.2e1) * qJD(1) ^ 2 + (t160 * t550 + t161 * t552) * t486 + t524 * t355 + t526 * t77 - t528 * t269 - t529 * t267 + t530 * qJD(3) + t333 * t477 - t231 * t370 - t237 * t371 + (mrSges(6,1) * t398 + mrSges(6,2) * t397) * t140 + (mrSges(7,1) * t398 + mrSges(7,2) * t397) * t78 + t531 * t76 - (-Ifges(3,2) * t274 + t266 + t510) * t394 / 0.2e1 - (t214 * t335 + t255 * t332 - t316 * t339) * qJD(4) / 0.2e1 - (t214 * t307 + t255 * t306 - t308 * t316) * qJD(1) / 0.2e1 - t544 * t325 / 0.2e1 + (-t324 * t548 - t325 * t546) * t478 + (-t324 * t550 - t325 * t547) * t498 + t20 * (mrSges(7,1) * t325 - mrSges(7,2) * t324) + t74 * (mrSges(6,1) * t325 - mrSges(6,2) * t324) + (-t324 * t552 - t325 * t550) * t499 + (-t137 * t546 - t138 * t548) * t471 + (-t137 * t547 - t138 * t550) * t488 + (-t161 / 0.2e1 - t138 / 0.2e1) * t536 + (-t137 * t550 - t138 * t552) * t485 + (-t98 * (mrSges(5,3) * t293 * t296 - mrSges(5,2) * t297) - t320 - m(4) * t326 * pkin(7) - t97 * (mrSges(5,1) * t297 - mrSges(5,3) * t414)) * qJD(1) + t292 * t497 - t513 * t274 / 0.2e1 + (-pkin(2) * t170 - qJ(3) * t162 - qJD(3) * t237 - t191 * t221) * m(4) + (-m(7) * (t396 - t418) - m(6) * (t396 - t403) - m(4) * t396 - m(5) * t348 - t297 * mrSges(5,3) + t569 * t293 + t558) * g(3) + t336 * t491 + t340 * t492 + t559 * mrSges(7,3) + t560 * mrSges(6,3) + (Ifges(3,3) + Ifges(4,1)) * qJDD(2) + (-qJ(3) * t372 * t404 + t298 * t508) * g(1) + (-m(5) * t304 - t556) * t494 + t114 * t360 + t185 * t274 / 0.2e1 + t255 * t343 * t190 + t130 * t342 + t516 * t344 + (-t389 * t98 + t390 * t97 - t517) * mrSges(5,3) + t296 * t60 / 0.2e1 + t258 * t16 - t221 * t220 - t223 * t139 + t209 * mrSges(3,2) - t210 * mrSges(3,1) - pkin(2) * t176 + t170 * mrSges(4,2) - t162 * mrSges(4,3) + t163 * t15 + t141 * t33 + t142 * t35 - t122 * t146 - t121 * t147 + t107 * t32 + t108 * t34; (-t76 - t77 - t530) * qJD(2) + t176 + ((t220 + t327) * qJD(1) - t516 * t372) * t293 - (t32 + t33) * t324 + t542 * t325 + t372 * t297 * g(3) + (t102 + t101) * t398 + (t104 + t103) * t397 + (-qJD(2) * t78 - t559) * m(7) + (-qJD(2) * t140 - t560) * m(6) + (-qJD(2) * t190 - t274 * t328 + t304) * m(5) + (qJD(2) * t237 + t191 * t274 + t170) * m(4) + t556; -m(6) * (-t140 * t469 + t30 * t37 + t31 * t38) + t134 * t566 + (-t453 + t147) * t98 + (t375 - t27) * t101 + (t134 * t25 + t21 * t352) * mrSges(7,3) + (t134 * t31 + t30 * t352) * mrSges(6,3) + (t202 + (m(6) * t279 + m(7) * t230 - t276 * t554) * t297 + t527) * g(3) + (-t21 * t26 - t25 * t27 - t78 * t95 + t2 * t262 + (t291 * t3 + (-t21 * t291 + t25 * t295) * qJD(5)) * pkin(4)) * m(7) + t575 + (-t37 - t376) * t104 + t77 * t469 + t114 * t475 + t33 * t466 + (t454 - t146) * t97 + t316 * (Ifges(5,1) * t214 + t436) / 0.2e1 - (Ifges(5,2) * t316 + t115 + t204) * t214 / 0.2e1 - t255 * (Ifges(5,5) * t214 + Ifges(5,6) * t316) / 0.2e1 - t190 * (-mrSges(5,1) * t316 + mrSges(5,2) * t214) + (-t26 - t376) * t103 + (t291 * t5 + t295 * t6 + (-t291 * t30 + t295 * t31) * qJD(5)) * t501 + t542 * pkin(4) * t291 + (-m(7) * (t229 * t298 + t230 * t412) - mrSges(5,2) * t195 + t533 * t194 + t522) * g(2) + (-m(7) * (-t229 * t294 + t230 * t411) + mrSges(5,2) * t193 + t533 * t192 + t521) * g(1) + (t375 - t38) * t102 + t369 + t262 * t32 - t28 * mrSges(5,2) + t29 * mrSges(5,1) - t95 * t76; t2 * t500 + t21 * t450 - t24 * t101 + t537 * t485 + (t451 + t104) * t31 + (t452 - t102) * t30 + (-m(7) * (-t21 + t24) + t449 + t103) * t25 + ((t500 - t554) * t276 * t297 + t527) * g(3) + (-t168 * t500 + t522) * g(2) + (-t166 * t500 + t521) * g(1) + (t32 + (-m(7) * t78 - t76) * t134) * pkin(5) + t575; -t352 * t101 + t134 * t103 + (-g(3) * t293 + t21 * t134 - t25 * t352 - t297 * t516 + t20) * m(7) + t15;];
tau  = t1;
