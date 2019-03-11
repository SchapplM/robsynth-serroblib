% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR8_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR8_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR8_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR8_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR8_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:22:41
% EndTime: 2019-03-09 05:23:27
% DurationCPUTime: 31.73s
% Computational Cost: add. (12085->818), mult. (24435->1094), div. (0->0), fcn. (16507->14), ass. (0->368)
t537 = -Ifges(4,2) / 0.2e1;
t305 = sin(qJ(3));
t309 = cos(qJ(3));
t352 = pkin(3) * t309 + pkin(8) * t305;
t258 = t352 * qJD(1);
t311 = -pkin(1) - pkin(7);
t276 = qJD(1) * t311 + qJD(2);
t308 = cos(qJ(4));
t304 = sin(qJ(4));
t407 = t304 * t309;
t168 = t308 * t258 - t276 * t407;
t302 = -qJ(5) - pkin(8);
t357 = qJD(4) * t302;
t404 = t305 * t308;
t536 = -(pkin(4) * t309 + qJ(5) * t404) * qJD(1) - t168 - qJD(5) * t304 + t308 * t357;
t398 = t308 * t309;
t169 = t304 * t258 + t276 * t398;
t293 = t305 * qJD(1);
t371 = t304 * t293;
t384 = qJD(5) * t308;
t535 = qJ(5) * t371 - t304 * t357 + t169 - t384;
t279 = t293 + qJD(4);
t274 = qJD(6) + t279;
t441 = -t274 / 0.2e1;
t390 = qJD(3) * t308;
t392 = qJD(1) * t309;
t252 = -t304 * t392 + t390;
t253 = qJD(3) * t304 + t308 * t392;
t300 = sin(pkin(10));
t301 = cos(pkin(10));
t162 = t252 * t300 + t253 * t301;
t447 = -t162 / 0.2e1;
t353 = t301 * t252 - t253 * t300;
t449 = -t353 / 0.2e1;
t303 = sin(qJ(6));
t307 = cos(qJ(6));
t95 = t162 * t307 + t303 * t353;
t455 = -t95 / 0.2e1;
t507 = -t162 * t303 + t307 * t353;
t457 = -t507 / 0.2e1;
t534 = Ifges(6,5) * t447 + Ifges(7,5) * t455 + Ifges(6,6) * t449 + Ifges(7,6) * t457 + Ifges(7,3) * t441;
t533 = -t252 / 0.2e1;
t532 = -t253 / 0.2e1;
t439 = -t279 / 0.2e1;
t420 = Ifges(4,4) * t309;
t531 = t305 * t537 + t420 / 0.2e1;
t247 = t300 * t308 + t301 * t304;
t222 = t247 * qJD(1);
t191 = t305 * t222;
t221 = t247 * qJD(4);
t519 = t191 + t221;
t332 = t300 * t304 - t301 * t308;
t491 = t305 * t332;
t192 = qJD(1) * t491;
t480 = qJD(4) * t332;
t530 = t192 + t480;
t529 = Ifges(5,3) + Ifges(6,3);
t499 = t300 * t535 + t301 * t536;
t498 = t300 * t536 - t301 * t535;
t381 = qJD(1) * qJD(3);
t260 = qJDD(1) * t309 - t305 * t381;
t153 = qJD(4) * t252 + qJDD(3) * t304 + t260 * t308;
t261 = -t305 * qJDD(1) - t309 * t381;
t244 = qJDD(4) - t261;
t351 = pkin(3) * t305 - pkin(8) * t309;
t265 = qJ(2) + t351;
t226 = t265 * qJD(1);
t264 = t305 * t276;
t233 = qJD(3) * pkin(8) + t264;
t148 = t226 * t304 + t233 * t308;
t382 = qJD(1) * qJD(2);
t277 = qJDD(1) * qJ(2) + t382;
t155 = -pkin(3) * t261 - pkin(8) * t260 + t277;
t272 = qJDD(1) * t311 + qJDD(2);
t389 = qJD(3) * t309;
t176 = t305 * t272 + t276 * t389;
t172 = qJDD(3) * pkin(8) + t176;
t62 = -qJD(4) * t148 + t308 * t155 - t172 * t304;
t36 = pkin(4) * t244 - qJ(5) * t153 - qJD(5) * t253 + t62;
t154 = -qJD(4) * t253 + qJDD(3) * t308 - t260 * t304;
t386 = qJD(4) * t308;
t387 = qJD(4) * t304;
t61 = t304 * t155 + t308 * t172 + t226 * t386 - t233 * t387;
t45 = qJ(5) * t154 + qJD(5) * t252 + t61;
t12 = t300 * t36 + t301 * t45;
t85 = -t153 * t300 + t154 * t301;
t10 = pkin(9) * t85 + t12;
t516 = pkin(9) * t162;
t147 = t308 * t226 - t233 * t304;
t113 = -qJ(5) * t253 + t147;
t104 = pkin(4) * t279 + t113;
t114 = qJ(5) * t252 + t148;
t107 = t300 * t114;
t55 = t301 * t104 - t107;
t40 = pkin(5) * t279 - t516 + t55;
t502 = pkin(9) * t353;
t409 = t301 * t114;
t56 = t300 * t104 + t409;
t46 = t56 + t502;
t13 = -t303 * t46 + t307 * t40;
t11 = -t300 * t45 + t301 * t36;
t86 = t153 * t301 + t154 * t300;
t9 = pkin(5) * t244 - pkin(9) * t86 + t11;
t2 = qJD(6) * t13 + t10 * t307 + t303 * t9;
t14 = t303 * t40 + t307 * t46;
t3 = -qJD(6) * t14 - t10 * t303 + t307 * t9;
t528 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t299 = qJ(4) + pkin(10);
t289 = cos(t299);
t428 = pkin(4) * t308;
t263 = pkin(5) * t289 + t428;
t286 = pkin(3) + t428;
t527 = -m(6) * t286 - m(7) * (pkin(3) + t263);
t526 = -m(6) * t302 - m(7) * (-pkin(9) + t302) + mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t525 = t62 * mrSges(5,1) + t11 * mrSges(6,1) - t61 * mrSges(5,2) - t12 * mrSges(6,2);
t524 = t14 * mrSges(7,2) + t56 * mrSges(6,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t531 + Ifges(5,5) * t532 + Ifges(5,6) * t533 + t529 * t439 - t13 * mrSges(7,1) - t55 * mrSges(6,1) + t534;
t523 = t309 / 0.2e1;
t452 = -m(6) - m(7);
t522 = t452 - m(5);
t521 = -pkin(5) * t392 + pkin(9) * t530 + t499;
t520 = -pkin(9) * t519 + t498;
t206 = t247 * t309;
t485 = t332 * qJD(1) - qJD(3) * t206 + t305 * t480;
t208 = t332 * t309;
t484 = -qJD(3) * t208 - t221 * t305 - t222;
t306 = sin(qJ(1));
t310 = cos(qJ(1));
t518 = g(1) * t306 - g(2) * t310;
t473 = m(6) * pkin(4);
t429 = pkin(4) * t304;
t517 = m(6) * t429;
t91 = -mrSges(5,1) * t154 + mrSges(5,2) * t153;
t515 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t260 - t91;
t486 = -t264 + (t371 + t387) * pkin(4);
t240 = Ifges(5,4) * t252;
t144 = t253 * Ifges(5,1) + t279 * Ifges(5,5) + t240;
t421 = Ifges(4,4) * t305;
t344 = t309 * Ifges(4,1) - t421;
t514 = Ifges(4,5) * qJD(3) + qJD(1) * t344 + t308 * t144;
t513 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t252 + mrSges(5,2) * t253 + mrSges(4,3) * t392;
t385 = qJD(4) * t309;
t367 = t308 * t385;
t391 = qJD(3) * t305;
t370 = t304 * t391;
t512 = t367 - t370;
t511 = Ifges(5,5) * t153 + Ifges(6,5) * t86 + Ifges(5,6) * t154 + Ifges(6,6) * t85 + t244 * t529;
t175 = t272 * t309 - t276 * t391;
t334 = t175 * t309 + t176 * t305;
t349 = mrSges(4,1) * t309 - mrSges(4,2) * t305;
t510 = (-Ifges(4,1) * t305 - t420) * t523 + qJ(2) * t349;
t508 = -m(4) + t522;
t288 = sin(t299);
t262 = pkin(5) * t288 + t429;
t506 = m(7) * t262 + mrSges(2,1) - mrSges(3,2) + mrSges(4,3);
t348 = mrSges(4,1) * t305 + mrSges(4,2) * t309;
t505 = -m(5) * t351 + t305 * t527 + t309 * t526 + mrSges(2,2) - mrSges(3,3) - t348;
t22 = qJD(6) * t507 + t303 * t85 + t307 * t86;
t472 = t22 / 0.2e1;
t23 = -qJD(6) * t95 - t303 * t86 + t307 * t85;
t471 = t23 / 0.2e1;
t459 = t85 / 0.2e1;
t458 = t86 / 0.2e1;
t451 = t153 / 0.2e1;
t450 = t154 / 0.2e1;
t238 = qJDD(6) + t244;
t445 = t238 / 0.2e1;
t444 = t244 / 0.2e1;
t44 = -t85 * mrSges(6,1) + t86 * mrSges(6,2);
t8 = -t23 * mrSges(7,1) + t22 * mrSges(7,2);
t503 = -t44 - t8;
t269 = t302 * t304;
t271 = t302 * t308;
t173 = t301 * t269 + t271 * t300;
t137 = -pkin(9) * t247 + t173;
t174 = t300 * t269 - t301 * t271;
t138 = -pkin(9) * t332 + t174;
t68 = t137 * t303 + t138 * t307;
t501 = -qJD(6) * t68 - t303 * t520 + t307 * t521;
t67 = t137 * t307 - t138 * t303;
t500 = qJD(6) * t67 + t303 * t521 + t307 * t520;
t205 = t247 * t305;
t126 = -t205 * t303 - t307 * t491;
t497 = -qJD(6) * t126 - t303 * t484 + t307 * t485;
t124 = -t205 * t307 + t303 * t491;
t496 = qJD(6) * t124 + t303 * t485 + t307 * t484;
t430 = pkin(4) * t301;
t284 = pkin(5) + t430;
t431 = pkin(4) * t300;
t212 = t284 * t307 - t303 * t431;
t59 = -t113 * t300 - t409;
t49 = t59 - t502;
t60 = t301 * t113 - t107;
t50 = t60 - t516;
t495 = t212 * qJD(6) - t303 * t49 - t307 * t50;
t213 = t284 * t303 + t307 * t431;
t494 = -t213 * qJD(6) + t303 * t50 - t307 * t49;
t493 = -t473 - mrSges(5,1);
t489 = t309 * t518;
t487 = pkin(5) * t519 + t486;
t402 = t305 * t311;
t186 = t304 * t265 + t308 * t402;
t111 = mrSges(5,1) * t244 - mrSges(5,3) * t153;
t112 = -mrSges(5,2) * t244 + mrSges(5,3) * t154;
t483 = -t304 * t111 + t308 * t112;
t482 = -t304 * t62 + t308 * t61;
t291 = qJ(6) + t299;
t282 = sin(t291);
t283 = cos(t291);
t347 = -mrSges(5,1) * t308 + mrSges(5,2) * t304;
t479 = m(5) * pkin(3) + mrSges(6,1) * t289 + mrSges(7,1) * t283 - mrSges(6,2) * t288 - mrSges(7,2) * t282 - t347 - t527;
t478 = -m(5) * pkin(8) - t526;
t476 = qJD(1) ^ 2;
t475 = Ifges(7,4) * t472 + Ifges(7,2) * t471 + Ifges(7,6) * t445;
t474 = Ifges(7,1) * t472 + Ifges(7,4) * t471 + Ifges(7,5) * t445;
t470 = Ifges(6,4) * t458 + Ifges(6,2) * t459 + Ifges(6,6) * t444;
t469 = Ifges(6,1) * t458 + Ifges(6,4) * t459 + Ifges(6,5) * t444;
t433 = Ifges(7,4) * t95;
t42 = Ifges(7,2) * t507 + t274 * Ifges(7,6) + t433;
t468 = -t42 / 0.2e1;
t467 = t42 / 0.2e1;
t88 = Ifges(7,4) * t507;
t43 = t95 * Ifges(7,1) + t274 * Ifges(7,5) + t88;
t466 = -t43 / 0.2e1;
t465 = t43 / 0.2e1;
t464 = Ifges(5,1) * t451 + Ifges(5,4) * t450 + Ifges(5,5) * t444;
t83 = t162 * Ifges(6,4) + Ifges(6,2) * t353 + t279 * Ifges(6,6);
t463 = -t83 / 0.2e1;
t462 = t83 / 0.2e1;
t84 = t162 * Ifges(6,1) + Ifges(6,4) * t353 + t279 * Ifges(6,5);
t461 = -t84 / 0.2e1;
t460 = t84 / 0.2e1;
t456 = t507 / 0.2e1;
t454 = t95 / 0.2e1;
t453 = -m(3) - m(4);
t448 = t353 / 0.2e1;
t446 = t162 / 0.2e1;
t442 = t253 / 0.2e1;
t440 = t274 / 0.2e1;
t438 = t279 / 0.2e1;
t436 = mrSges(6,3) * t55;
t435 = mrSges(7,3) * t13;
t434 = mrSges(7,3) * t14;
t432 = pkin(4) * t253;
t425 = g(3) * t309;
t424 = t56 * mrSges(6,3);
t245 = qJD(3) * t352 + qJD(2);
t319 = -qJD(4) * t186 + t308 * t245;
t358 = -t304 * t311 + pkin(4);
t369 = t305 * t390;
t77 = qJ(5) * t369 + (qJ(5) * t387 + qJD(3) * t358 - t384) * t309 + t319;
t388 = qJD(3) * t311;
t368 = t309 * t388;
t373 = t304 * t245 + t265 * t386 + t308 * t368;
t87 = -qJ(5) * t367 + (-qJD(5) * t309 + (qJ(5) * qJD(3) - qJD(4) * t311) * t305) * t304 + t373;
t38 = t300 * t77 + t301 * t87;
t419 = Ifges(5,4) * t304;
t418 = Ifges(5,4) * t308;
t417 = t147 * mrSges(5,3);
t416 = t148 * mrSges(5,3);
t415 = t253 * Ifges(5,4);
t410 = t276 * t309;
t406 = t304 * t310;
t405 = t305 * t306;
t403 = t305 * t310;
t401 = t306 * t308;
t397 = t308 * t310;
t396 = t309 * t311;
t243 = t308 * t265;
t156 = -qJ(5) * t398 + t305 * t358 + t243;
t167 = -qJ(5) * t407 + t186;
t97 = t300 * t156 + t301 * t167;
t187 = -t282 * t405 + t283 * t310;
t188 = t282 * t310 + t283 * t405;
t395 = t187 * mrSges(7,1) - t188 * mrSges(7,2);
t189 = t282 * t403 + t283 * t306;
t190 = -t282 * t306 + t283 * t403;
t394 = t189 * mrSges(7,1) + t190 * mrSges(7,2);
t393 = t310 * pkin(1) + t306 * qJ(2);
t383 = qJDD(1) * mrSges(3,2);
t379 = Ifges(7,5) * t22 + Ifges(7,6) * t23 + Ifges(7,3) * t238;
t375 = t304 * t402;
t278 = t305 * t388;
t360 = -t385 / 0.2e1;
t37 = -t300 * t87 + t301 * t77;
t356 = -t381 / 0.2e1;
t354 = (t277 + t382) * qJ(2);
t96 = t301 * t156 - t167 * t300;
t248 = pkin(4) * t407 - t396;
t346 = mrSges(5,1) * t304 + mrSges(5,2) * t308;
t345 = -mrSges(7,1) * t282 - mrSges(7,2) * t283;
t343 = Ifges(5,1) * t308 - t419;
t342 = Ifges(5,1) * t304 + t418;
t340 = -Ifges(5,2) * t304 + t418;
t339 = Ifges(5,2) * t308 + t419;
t338 = -Ifges(4,5) * t305 - Ifges(4,6) * t309;
t337 = Ifges(5,5) * t308 - Ifges(5,6) * t304;
t336 = Ifges(5,5) * t304 + Ifges(5,6) * t308;
t65 = pkin(5) * t305 + pkin(9) * t208 + t96;
t66 = -pkin(9) * t206 + t97;
t30 = -t303 * t66 + t307 * t65;
t31 = t303 * t65 + t307 * t66;
t234 = -qJD(3) * pkin(3) - t410;
t335 = t147 * t308 + t148 * t304;
t183 = -mrSges(5,2) * t279 + mrSges(5,3) * t252;
t184 = mrSges(5,1) * t279 - mrSges(5,3) * t253;
t333 = -t304 * t183 - t308 * t184;
t125 = -t206 * t307 + t208 * t303;
t127 = -t206 * t303 - t208 * t307;
t157 = -t247 * t303 - t307 * t332;
t158 = t247 * t307 - t303 * t332;
t331 = t379 + t528;
t229 = t304 * t403 + t401;
t227 = -t304 * t405 + t397;
t171 = -qJDD(3) * pkin(3) - t175;
t329 = t305 * (-Ifges(4,2) * t309 - t421);
t177 = pkin(4) * t512 + t278;
t322 = t304 * t385 + t369;
t166 = -pkin(4) * t252 + qJD(5) + t234;
t316 = Ifges(5,5) * t309 - t305 * t343;
t315 = Ifges(5,6) * t309 - t305 * t340;
t314 = Ifges(5,3) * t309 - t305 * t337;
t100 = -pkin(4) * t154 + qJDD(5) + t171;
t313 = -qJD(4) * t335 + t482;
t287 = -pkin(1) * qJDD(1) + qJDD(2);
t267 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t293;
t257 = t348 * qJD(1);
t239 = t346 * t309;
t230 = -t304 * t306 + t305 * t397;
t228 = t305 * t401 + t406;
t211 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t261;
t204 = -t288 * t306 + t289 * t403;
t203 = t288 * t403 + t289 * t306;
t202 = t288 * t310 + t289 * t405;
t201 = -t288 * t405 + t289 * t310;
t193 = pkin(5) * t332 - t286;
t185 = t243 - t375;
t164 = pkin(5) * t206 + t248;
t143 = t252 * Ifges(5,2) + t279 * Ifges(5,6) + t415;
t134 = qJD(3) * t491 - t221 * t309;
t132 = t247 * t391 + t309 * t480;
t123 = mrSges(6,1) * t279 - mrSges(6,3) * t162;
t122 = -mrSges(6,2) * t279 + mrSges(6,3) * t353;
t119 = t191 * t303 + t192 * t307;
t118 = t191 * t307 - t192 * t303;
t117 = pkin(5) * t162 + t432;
t116 = -t304 * t368 + t319;
t115 = -qJD(4) * t375 + t373;
t101 = -pkin(5) * t353 + t166;
t99 = -pkin(5) * t132 + t177;
t98 = -mrSges(6,1) * t353 + mrSges(6,2) * t162;
t90 = -qJD(6) * t158 - t221 * t307 + t303 * t480;
t89 = qJD(6) * t157 - t221 * t303 - t307 * t480;
t74 = mrSges(7,1) * t274 - mrSges(7,3) * t95;
t73 = -mrSges(7,2) * t274 + mrSges(7,3) * t507;
t71 = t153 * Ifges(5,4) + t154 * Ifges(5,2) + t244 * Ifges(5,6);
t64 = mrSges(6,1) * t244 - mrSges(6,3) * t86;
t63 = -mrSges(6,2) * t244 + mrSges(6,3) * t85;
t54 = -qJD(6) * t127 + t132 * t307 - t134 * t303;
t52 = qJD(6) * t125 + t132 * t303 + t134 * t307;
t48 = -pkin(5) * t85 + t100;
t47 = -mrSges(7,1) * t507 + mrSges(7,2) * t95;
t29 = pkin(9) * t132 + t38;
t28 = pkin(5) * t389 - pkin(9) * t134 + t37;
t18 = -mrSges(7,2) * t238 + mrSges(7,3) * t23;
t17 = mrSges(7,1) * t238 - mrSges(7,3) * t22;
t5 = -qJD(6) * t31 + t28 * t307 - t29 * t303;
t4 = qJD(6) * t30 + t28 * t303 + t29 * t307;
t1 = [m(5) * (t234 * t311 * t391 + t115 * t148 + t116 * t147 + t185 * t62 + t186 * t61) + (-Ifges(6,5) * t208 - Ifges(6,6) * t206) * t444 + (-Ifges(6,4) * t208 - Ifges(6,2) * t206) * t459 + (Ifges(6,4) * t134 + Ifges(6,2) * t132) * t448 + t261 * t531 + (Ifges(7,4) * t52 + Ifges(7,2) * t54) * t456 + (Ifges(7,4) * t127 + Ifges(7,2) * t125) * t471 + (Ifges(6,5) * t134 + Ifges(6,6) * t132 + qJD(3) * t314 - t336 * t385) * t438 + (Ifges(2,3) + Ifges(3,1)) * qJDD(1) + (t370 / 0.2e1 + t308 * t360) * t143 + (Ifges(4,1) * t260 + Ifges(4,4) * t261) * t523 + (Ifges(7,1) * t52 + Ifges(7,4) * t54) * t454 + (Ifges(7,1) * t127 + Ifges(7,4) * t125) * t472 + (t147 * t322 - t148 * t512 - t398 * t62 - t407 * t61) * mrSges(5,3) + t234 * (mrSges(5,1) * t512 - mrSges(5,2) * t322) + (-Ifges(6,1) * t208 - Ifges(6,4) * t206) * t458 + (Ifges(6,1) * t134 + Ifges(6,4) * t132) * t446 + t100 * (mrSges(6,1) * t206 - mrSges(6,2) * t208) + (Ifges(6,5) * t458 + Ifges(6,6) * t459 + t529 * t444 + t379 / 0.2e1 + t511 / 0.2e1 + Ifges(7,3) * t445 + Ifges(5,6) * t450 + Ifges(5,5) * t451 - Ifges(4,4) * t260 / 0.2e1 + t261 * t537 - Ifges(4,6) * qJDD(3) + Ifges(7,6) * t471 + Ifges(7,5) * t472 + t525 + t528) * t305 + m(4) * (t311 * t334 + t354) + m(3) * (-pkin(1) * t287 + t354) + (t125 * t2 - t127 * t3 - t13 * t52 + t14 * t54) * mrSges(7,3) - t334 * mrSges(4,3) + t513 * t278 - t514 * t391 / 0.2e1 + t510 * t381 + (t11 * t208 - t12 * t206 + t132 * t56 - t134 * t55) * mrSges(6,3) + qJD(3) ^ 2 * t338 / 0.2e1 + m(6) * (t100 * t248 + t11 * t96 + t12 * t97 + t166 * t177 + t37 * t55 + t38 * t56) + m(7) * (t101 * t99 + t13 * t5 + t14 * t4 + t164 * t48 + t2 * t31 + t3 * t30) + (qJD(3) * t316 - t342 * t385) * t442 + (t348 + 0.2e1 * mrSges(3,3)) * t277 + t260 * t344 / 0.2e1 - t71 * t407 / 0.2e1 + (-m(5) * t171 * t311 + Ifges(4,5) * qJDD(3) + t337 * t444 + t340 * t450 + t343 * t451) * t309 + t329 * t356 + (Ifges(7,5) * t52 + Ifges(7,6) * t54) * t440 + (Ifges(7,5) * t127 + Ifges(7,6) * t125) * t445 + (t147 * mrSges(5,1) - t148 * mrSges(5,2) + Ifges(6,5) * t446 + Ifges(7,5) * t454 + Ifges(6,6) * t448 + Ifges(7,6) * t456 + Ifges(6,3) * t438 + Ifges(7,3) * t440 - t524) * t389 - pkin(1) * t383 + t252 * (qJD(3) * t315 - t339 * t385) / 0.2e1 + t287 * mrSges(3,2) + qJ(2) * (-mrSges(4,1) * t261 + mrSges(4,2) * t260) + qJD(2) * t257 + t248 * t44 + t171 * t239 + t185 * t111 + t186 * t112 + t177 * t98 + t115 * t183 + t116 * t184 + t166 * (-mrSges(6,1) * t132 + mrSges(6,2) * t134) + t164 * t8 + t37 * t123 + t48 * (-mrSges(7,1) * t125 + mrSges(7,2) * t127) + t38 * t122 + t101 * (-mrSges(7,1) * t54 + mrSges(7,2) * t52) + t96 * t64 + t97 * t63 + t99 * t47 + t4 * t73 + t5 * t74 + t134 * t460 + t132 * t462 + t515 * t396 + t304 * t144 * t360 + t267 * t368 + (-t406 * t473 - m(3) * t393 - t228 * mrSges(5,1) - t202 * mrSges(6,1) - t188 * mrSges(7,1) - t227 * mrSges(5,2) - t201 * mrSges(6,2) - t187 * mrSges(7,2) + t508 * (t310 * pkin(7) + t393) - t506 * t310 + t505 * t306) * g(2) + t398 * t464 + t52 * t465 + t54 * t467 - t208 * t469 - t206 * t470 + t127 * t474 + t125 * t475 + t211 * t402 + (-t230 * mrSges(5,1) - t204 * mrSges(6,1) - t190 * mrSges(7,1) + t229 * mrSges(5,2) + t203 * mrSges(6,2) + t189 * mrSges(7,2) + (m(3) * pkin(1) + t311 * t508 + t506 + t517) * t306 + ((-m(3) + t508) * qJ(2) + t505) * t310) * g(1) + t30 * t17 + t31 * t18; t383 + t124 * t17 + t126 * t18 - t205 * t64 - t491 * t63 + t497 * t74 + t496 * t73 + t485 * t123 + t484 * t122 + (qJ(2) * t453 - mrSges(3,3)) * t476 + (-t257 + t333) * qJD(1) + ((t183 * t308 - t184 * t304 + t267) * qJD(3) + t503 + t515) * t309 + (t211 + t333 * qJD(4) + (t47 + t98 + t513) * qJD(3) + t483) * t305 + m(3) * t287 + m(4) * t334 - t518 * (-t453 - t522) + (t101 * t391 + t124 * t3 + t126 * t2 + t13 * t497 + t14 * t496 - t309 * t48) * m(7) + (-t100 * t309 - t11 * t205 - t12 * t491 + t166 * t391 + t484 * t56 + t485 * t55) * m(6) + ((-t171 + (-t147 * t304 + t148 * t308) * qJD(3)) * t309 + (qJD(3) * t234 + t313) * t305 - t335 * qJD(1)) * m(5); (t519 * mrSges(6,1) - mrSges(6,2) * t530) * t166 - t513 * t264 + t518 * (t305 * t478 - t309 * t479 - t349) + (-t118 * t14 + t119 * t13 + t157 * t2 - t158 * t3) * mrSges(7,3) + (Ifges(7,4) * t119 + Ifges(7,2) * t118) * t457 + (Ifges(6,1) * t192 + Ifges(6,4) * t191) * t447 + t514 * t293 / 0.2e1 + (t329 / 0.2e1 - t510) * t476 + t279 * (t234 * t346 - t304 * t143 / 0.2e1) - t387 * t416 + (Ifges(6,3) * t439 + t524 + t534) * t392 + (Ifges(7,1) * t119 + Ifges(7,4) * t118) * t455 + (Ifges(7,5) * t158 + Ifges(7,6) * t157) * t445 + t339 * t450 + t342 * t451 + (-t148 * (mrSges(5,3) * t304 * t305 - mrSges(5,2) * t309) - t147 * (mrSges(5,1) * t309 + mrSges(5,3) * t404)) * qJD(1) + t171 * t347 + ((-t119 + t89) * mrSges(7,2) + (t118 - t90) * mrSges(7,1)) * t101 + (Ifges(6,4) * t192 + Ifges(6,2) * t191) * t449 + (-t417 + t144 / 0.2e1) * t386 + (t252 * t340 + t253 * t343 + t279 * t337) * qJD(4) / 0.2e1 - (t252 * t315 + t253 * t316 + t279 * t314) * qJD(1) / 0.2e1 + t338 * t356 - t267 * t410 - (Ifges(6,1) * t446 + Ifges(6,4) * t448 + Ifges(6,5) * t438 - t436 + t460) * t480 + (-t11 * t247 - t12 * t332 - t191 * t56 + t192 * t55) * mrSges(6,3) + (Ifges(6,5) * t247 - Ifges(6,6) * t332 + t336) * t444 + t100 * (mrSges(6,1) * t332 + mrSges(6,2) * t247) + (Ifges(6,1) * t247 - Ifges(6,4) * t332) * t458 + (Ifges(6,4) * t247 - Ifges(6,2) * t332) * t459 - t332 * t470 + (t305 * t479 + t309 * t478 + t348) * g(3) + (Ifges(7,4) * t454 + Ifges(7,2) * t456 + Ifges(7,6) * t440 + t434 + t467) * t90 + (Ifges(7,1) * t454 + Ifges(7,4) * t456 + Ifges(7,5) * t440 - t435 + t465) * t89 + (Ifges(7,5) * t119 + Ifges(7,6) * t118) * t441 + (Ifges(6,5) * t192 + Ifges(6,6) * t191) * t439 + (-pkin(3) * t171 - t147 * t168 - t148 * t169 - t234 * t264) * m(5) + t308 * t71 / 0.2e1 - t286 * t44 + Ifges(4,5) * t260 + Ifges(4,6) * t261 + t193 * t8 - t176 * mrSges(4,2) - t169 * t183 - t168 * t184 + t173 * t64 + t174 * t63 + t175 * mrSges(4,1) + t48 * (-mrSges(7,1) * t157 + mrSges(7,2) * t158) + Ifges(4,3) * qJDD(3) - pkin(3) * t91 + t67 * t17 + t68 * t18 + t192 * t461 + t191 * t463 + t482 * mrSges(5,3) + (m(5) * t313 - t183 * t387 - t184 * t386 + t483) * pkin(8) + t486 * t98 + t487 * t47 + t498 * t122 + t499 * t123 + t304 * t464 + t119 * t466 + t118 * t468 + t247 * t469 + (Ifges(7,4) * t158 + Ifges(7,2) * t157) * t471 + (Ifges(7,1) * t158 + Ifges(7,4) * t157) * t472 + t158 * t474 + t157 * t475 + (-t100 * t286 + t11 * t173 + t12 * t174 + t166 * t486 + t498 * t56 + t499 * t55) * m(6) + t500 * t73 + (t101 * t487 + t13 * t501 + t14 * t500 + t193 * t48 + t2 * t68 + t3 * t67) * m(7) + t501 * t74 - (Ifges(6,4) * t446 + Ifges(6,2) * t448 + Ifges(6,6) * t438 + t424 + t462) * t221; (Ifges(5,1) * t252 - t415) * t532 + (-Ifges(5,2) * t253 + t144 + t240) * t533 + (-mrSges(7,2) * t101 + Ifges(7,1) * t455 + Ifges(7,4) * t457 + Ifges(7,5) * t441 + t435 + t466) * t507 + t525 - (mrSges(7,1) * t101 + Ifges(7,4) * t455 + Ifges(7,2) * t457 + Ifges(7,6) * t441 - t434 + t468) * t95 + t331 + t253 * t416 + t252 * t417 + t64 * t430 + t63 * t431 + t511 + t143 * t442 + (Ifges(5,5) * t252 - Ifges(5,6) * t253) * t439 + (-mrSges(6,2) * t166 + Ifges(6,1) * t447 + Ifges(6,4) * t449 + Ifges(6,5) * t439 + t436 + t461) * t353 - (mrSges(6,1) * t166 + Ifges(6,4) * t447 + Ifges(6,2) * t449 + Ifges(6,6) * t439 - t424 + t463) * t162 - t98 * t432 - m(6) * (t166 * t432 + t55 * t59 + t56 * t60) - t234 * (mrSges(5,1) * t253 + mrSges(5,2) * t252) + g(3) * t239 + t212 * t17 + t213 * t18 - t147 * t183 + t148 * t184 - t60 * t122 - t59 * t123 - t117 * t47 + (mrSges(6,1) * t288 + mrSges(6,2) * t289 - t345 + t517) * t425 + (-t203 * mrSges(6,1) - t204 * mrSges(6,2) - mrSges(5,2) * t230 - m(7) * (t262 * t403 + t263 * t306) - t394 + t493 * t229) * g(2) + (-m(7) * (-t262 * t405 + t263 * t310) - t395 - t201 * mrSges(6,1) + t202 * mrSges(6,2) + mrSges(5,2) * t228 + t493 * t227) * g(1) + t494 * t74 + t495 * t73 + (-t101 * t117 + t13 * t494 + t14 * t495 + t2 * t213 + t212 * t3 + t262 * t425) * m(7) + (t11 * t301 + t12 * t300) * t473; t452 * t305 * g(3) - t353 * t122 + t162 * t123 - t507 * t73 + t95 * t74 + (t13 * t95 - t14 * t507 + t48 + t489) * m(7) + (t162 * t55 - t353 * t56 + t100 + t489) * m(6) - t503; -t101 * (mrSges(7,1) * t95 + mrSges(7,2) * t507) + (Ifges(7,1) * t507 - t433) * t455 + t42 * t454 + (Ifges(7,5) * t507 - Ifges(7,6) * t95) * t441 - t13 * t73 + t14 * t74 - g(1) * t395 - g(2) * t394 - t345 * t425 + (t13 * t507 + t14 * t95) * mrSges(7,3) + t331 + (-Ifges(7,2) * t95 + t43 + t88) * t457;];
tau  = t1;
