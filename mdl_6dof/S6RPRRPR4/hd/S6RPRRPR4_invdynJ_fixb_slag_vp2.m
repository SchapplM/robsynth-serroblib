% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:09:04
% EndTime: 2019-03-09 05:09:38
% DurationCPUTime: 23.94s
% Computational Cost: add. (19493->744), mult. (48143->964), div. (0->0), fcn. (38514->18), ass. (0->344)
t307 = sin(pkin(11));
t309 = cos(pkin(11));
t313 = sin(qJ(6));
t317 = cos(qJ(6));
t340 = t307 * t313 - t309 * t317;
t237 = t340 * qJD(6);
t308 = sin(pkin(10));
t315 = sin(qJ(3));
t310 = cos(pkin(10));
t318 = cos(qJ(3));
t392 = t310 * t318;
t339 = t308 * t315 - t392;
t334 = t339 * qJD(1);
t444 = cos(qJ(4));
t227 = t444 * t334;
t255 = t308 * t318 + t310 * t315;
t236 = t255 * qJD(1);
t314 = sin(qJ(4));
t199 = t236 * t314 + t227;
t515 = t340 * t199;
t534 = t515 + t237;
t254 = t307 * t317 + t309 * t313;
t238 = t254 * qJD(6);
t516 = t254 * t199;
t533 = t516 + t238;
t381 = qJD(1) * qJD(2);
t279 = qJ(2) * qJDD(1) + t381;
t518 = t199 * t307;
t532 = pkin(5) * t518;
t531 = pkin(9) * t518;
t530 = -mrSges(6,3) - mrSges(7,3);
t302 = t308 ^ 2;
t303 = t310 ^ 2;
t385 = t302 + t303;
t305 = pkin(10) + qJ(3);
t298 = qJ(4) + t305;
t286 = sin(t298);
t304 = pkin(11) + qJ(6);
t294 = sin(t304);
t426 = mrSges(7,2) * t294;
t427 = mrSges(6,2) * t307;
t529 = (-t426 - t427) * t286;
t431 = pkin(7) + qJ(2);
t272 = t431 * t308;
t256 = qJD(1) * t272;
t274 = t431 * t310;
t257 = qJD(1) * t274;
t211 = -t315 * t256 + t318 * t257;
t183 = -pkin(8) * t334 + t211;
t175 = t314 * t183;
t401 = t257 * t315;
t210 = -t318 * t256 - t401;
t182 = -pkin(8) * t236 + t210;
t177 = qJD(3) * pkin(3) + t182;
t128 = t177 * t444 - t175;
t306 = qJD(3) + qJD(4);
t116 = -t306 * pkin(4) + qJD(5) - t128;
t328 = t314 * t334;
t323 = t236 * t444 - t328;
t180 = t306 * t307 + t309 * t323;
t360 = t309 * t306 - t307 * t323;
t490 = mrSges(5,1) * t306 + mrSges(6,1) * t360 - mrSges(6,2) * t180 - mrSges(5,3) * t323;
t528 = m(6) * t116 - t490;
t517 = t199 * t309;
t527 = pkin(5) * t323 + pkin(9) * t517;
t239 = t339 * qJD(3);
t206 = -qJD(1) * t239 + qJDD(1) * t255;
t240 = t255 * qJD(3);
t207 = -qJD(1) * t240 - qJDD(1) * t339;
t120 = t444 * t206 + (-qJD(4) * t236 + t207) * t314 - qJD(4) * t227;
t365 = qJD(4) * t444;
t121 = -qJD(4) * t328 + t314 * t206 - t444 * t207 + t236 * t365;
t176 = t444 * t183;
t129 = t314 * t177 + t176;
t122 = t306 * qJ(5) + t129;
t290 = t310 * pkin(2) + pkin(1);
t226 = pkin(3) * t339 - t290;
t212 = qJD(1) * t226 + qJD(2);
t130 = t199 * pkin(4) - qJ(5) * t323 + t212;
t72 = -t122 * t307 + t309 * t130;
t53 = pkin(5) * t199 - pkin(9) * t180 + t72;
t73 = t309 * t122 + t307 * t130;
t57 = pkin(9) * t360 + t73;
t17 = -t313 * t57 + t317 * t53;
t18 = t313 * t53 + t317 * t57;
t300 = qJDD(3) + qJDD(4);
t105 = -t120 * t307 + t300 * t309;
t361 = pkin(7) * qJDD(1) + t279;
t230 = t361 * t308;
t231 = t361 * t310;
t164 = -qJD(3) * t211 - t318 * t230 - t231 * t315;
t119 = qJDD(3) * pkin(3) - pkin(8) * t206 + t164;
t384 = qJD(3) * t318;
t163 = -qJD(3) * t401 - t315 * t230 + t318 * t231 - t256 * t384;
t124 = pkin(8) * t207 + t163;
t383 = qJD(4) * t314;
t49 = t314 * t119 + t444 * t124 + t177 * t365 - t183 * t383;
t46 = qJ(5) * t300 + qJD(5) * t306 + t49;
t261 = -qJDD(1) * t290 + qJDD(2);
t184 = -pkin(3) * t207 + t261;
t54 = pkin(4) * t121 - qJ(5) * t120 - qJD(5) * t323 + t184;
t15 = t307 * t54 + t309 * t46;
t11 = pkin(9) * t105 + t15;
t106 = t120 * t309 + t300 * t307;
t14 = -t307 * t46 + t309 * t54;
t6 = pkin(5) * t121 - pkin(9) * t106 + t14;
t2 = qJD(6) * t17 + t11 * t317 + t313 * t6;
t3 = -qJD(6) * t18 - t11 * t313 + t317 * t6;
t50 = t119 * t444 - t314 * t124 - t177 * t383 - t183 * t365;
t47 = -t300 * pkin(4) + qJDD(5) - t50;
t30 = -t105 * pkin(5) + t47;
t344 = -t14 * t307 + t15 * t309;
t430 = mrSges(6,1) * t309;
t349 = t427 - t430;
t37 = t106 * Ifges(6,4) + t105 * Ifges(6,2) + t121 * Ifges(6,6);
t422 = Ifges(6,4) * t309;
t423 = Ifges(6,4) * t307;
t445 = t309 / 0.2e1;
t195 = qJD(6) + t199;
t454 = t195 / 0.2e1;
t455 = -t195 / 0.2e1;
t134 = t180 * t317 + t313 * t360;
t459 = t134 / 0.2e1;
t460 = -t134 / 0.2e1;
t509 = -t180 * t313 + t317 * t360;
t461 = t509 / 0.2e1;
t462 = -t509 / 0.2e1;
t463 = t121 / 0.2e1;
t115 = qJDD(6) + t121;
t464 = t115 / 0.2e1;
t465 = t106 / 0.2e1;
t466 = t105 / 0.2e1;
t131 = Ifges(7,4) * t509;
t64 = Ifges(7,1) * t134 + Ifges(7,5) * t195 + t131;
t467 = t64 / 0.2e1;
t421 = Ifges(7,4) * t134;
t63 = Ifges(7,2) * t509 + Ifges(7,6) * t195 + t421;
t469 = t63 / 0.2e1;
t43 = -qJD(6) * t134 + t105 * t317 - t106 * t313;
t471 = t43 / 0.2e1;
t42 = qJD(6) * t509 + t105 * t313 + t106 * t317;
t472 = t42 / 0.2e1;
t473 = Ifges(6,1) * t465 + Ifges(6,4) * t466 + Ifges(6,5) * t463;
t474 = Ifges(7,1) * t472 + Ifges(7,4) * t471 + Ifges(7,5) * t464;
t475 = Ifges(7,4) * t472 + Ifges(7,2) * t471 + Ifges(7,6) * t464;
t89 = -pkin(5) * t360 + t116;
t526 = (Ifges(7,5) * t254 - Ifges(7,6) * t340) * t464 + (Ifges(7,4) * t254 - Ifges(7,2) * t340) * t471 + (Ifges(7,1) * t254 - Ifges(7,4) * t340) * t472 + t30 * (mrSges(7,1) * t340 + mrSges(7,2) * t254) - t340 * t475 + (t534 * t17 - t533 * t18 - t2 * t340 - t254 * t3) * mrSges(7,3) + (t533 * mrSges(7,1) - t534 * mrSges(7,2)) * t89 - t515 * t64 / 0.2e1 - t516 * t63 / 0.2e1 + (Ifges(6,5) * t307 + Ifges(6,6) * t309) * t463 + (Ifges(6,1) * t307 + t422) * t465 + (-Ifges(7,5) * t237 - Ifges(7,6) * t238) * t454 + (-Ifges(7,1) * t237 - Ifges(7,4) * t238) * t459 + (-Ifges(7,4) * t237 - Ifges(7,2) * t238) * t461 + t37 * t445 + (Ifges(7,5) * t515 + Ifges(7,6) * t516) * t455 + (Ifges(7,4) * t515 + Ifges(7,2) * t516) * t462 + (Ifges(7,1) * t515 + Ifges(7,4) * t516) * t460 + Ifges(5,5) * t120 - Ifges(5,6) * t121 - t49 * mrSges(5,2) + t50 * mrSges(5,1) + (Ifges(6,2) * t309 + t423) * t466 - t237 * t467 - t238 * t469 + t307 * t473 + t254 * t474 + Ifges(5,3) * t300 + (-t517 * t72 - t518 * t73 + t344) * mrSges(6,3) + t47 * t349;
t456 = -t180 / 0.2e1;
t452 = -t199 / 0.2e1;
t457 = -t360 / 0.2e1;
t525 = m(5) * t128;
t287 = cos(t298);
t432 = t309 * pkin(5);
t289 = pkin(4) + t432;
t311 = -pkin(9) - qJ(5);
t488 = -t286 * t311 + t287 * t289;
t523 = m(7) * t488;
t521 = t180 * Ifges(6,5);
t522 = Ifges(6,6) * t360;
t500 = t134 * Ifges(7,5) + Ifges(7,6) * t509 + t199 * Ifges(6,3) + t195 * Ifges(7,3) + t521 + t522;
t520 = t306 * Ifges(5,5);
t519 = t306 * Ifges(5,6);
t514 = -t287 * mrSges(5,1) + (mrSges(5,2) - mrSges(7,3)) * t286;
t296 = cos(t304);
t429 = mrSges(7,1) * t296;
t512 = (t429 + t430) * t286;
t348 = mrSges(6,1) * t307 + mrSges(6,2) * t309;
t335 = t116 * t348;
t504 = Ifges(6,4) * t456 + Ifges(6,2) * t457 + Ifges(6,6) * t452;
t99 = t180 * Ifges(6,1) + Ifges(6,4) * t360 + t199 * Ifges(6,5);
t511 = t307 * t504 + t99 * t445 + t335;
t510 = t212 * mrSges(5,2) - t128 * mrSges(5,3);
t158 = pkin(4) * t323 + qJ(5) * t199;
t508 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t301 = -pkin(8) - t431;
t443 = m(3) * qJ(2);
t507 = -m(4) * t431 + m(6) * t301 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3) - t348 - t443;
t295 = sin(t305);
t297 = cos(t305);
t353 = mrSges(4,1) * t297 - mrSges(4,2) * t295;
t354 = -mrSges(3,1) * t310 + mrSges(3,2) * t308;
t506 = -m(3) * pkin(1) - m(4) * t290 - mrSges(2,1) - t353 + t354 - (m(6) * pkin(4) - t349) * t287 + t514;
t505 = t212 * mrSges(5,1) + t72 * mrSges(6,1) + t17 * mrSges(7,1) - t73 * mrSges(6,2) - t18 * mrSges(7,2) - t129 * mrSges(5,3);
t503 = m(6) + m(7);
t450 = -t323 / 0.2e1;
t439 = pkin(3) * t314;
t288 = qJ(5) + t439;
t244 = (-pkin(9) - t288) * t307;
t299 = t309 * pkin(9);
t397 = t288 * t309;
t245 = t299 + t397;
t203 = t244 * t317 - t245 * t313;
t357 = pkin(3) * t365;
t280 = t357 + qJD(5);
t137 = t182 * t444 - t175;
t442 = pkin(3) * t236;
t141 = t158 + t442;
t79 = -t137 * t307 + t309 * t141;
t56 = t527 + t79;
t80 = t309 * t137 + t307 * t141;
t66 = t80 + t531;
t499 = qJD(6) * t203 - t280 * t340 - t313 * t56 - t317 * t66;
t204 = t244 * t313 + t245 * t317;
t498 = -qJD(6) * t204 - t254 * t280 + t313 * t66 - t317 * t56;
t271 = t311 * t307;
t412 = qJ(5) * t309;
t273 = t299 + t412;
t213 = t271 * t317 - t273 * t313;
t85 = -t128 * t307 + t309 * t158;
t58 = t527 + t85;
t86 = t309 * t128 + t307 * t158;
t71 = t86 + t531;
t497 = -qJD(5) * t340 + qJD(6) * t213 - t313 * t58 - t317 * t71;
t215 = t271 * t313 + t273 * t317;
t496 = -qJD(5) * t254 - qJD(6) * t215 + t313 * t71 - t317 * t58;
t214 = -t318 * t272 - t274 * t315;
t192 = -pkin(8) * t255 + t214;
t216 = -t315 * t272 + t318 * t274;
t193 = -pkin(8) * t339 + t216;
t489 = t444 * t192 - t314 * t193;
t341 = -t286 * t289 - t287 * t311;
t438 = pkin(4) * t286;
t440 = pkin(3) * t295;
t486 = -m(7) * (t341 - t440) - m(6) * (-t438 - t440) + t512;
t485 = -m(7) * t341 + t512;
t319 = cos(qJ(1));
t398 = t287 * t319;
t484 = t319 * t529 + t398 * t530;
t316 = sin(qJ(1));
t399 = t287 * t316;
t483 = t316 * t529 + t399 * t530;
t145 = -mrSges(6,2) * t199 + mrSges(6,3) * t360;
t146 = mrSges(6,1) * t199 - mrSges(6,3) * t180;
t482 = t145 * t309 - t146 * t307;
t481 = g(1) * t319 + g(2) * t316;
t277 = t286 * mrSges(6,3);
t480 = -t277 + t514 + (t426 - t429 + t349) * t287;
t78 = -mrSges(7,1) * t509 + mrSges(7,2) * t134;
t478 = m(7) * t89 + t528 + t78;
t345 = Ifges(6,5) * t309 - Ifges(6,6) * t307;
t346 = -Ifges(6,2) * t307 + t422;
t347 = Ifges(6,1) * t309 - t423;
t446 = -t306 / 0.2e1;
t477 = Ifges(5,1) * t450 + Ifges(5,5) * t446 + t345 * t452 + t346 * t457 + t347 * t456 - t510;
t451 = t199 / 0.2e1;
t476 = Ifges(6,5) * t456 + Ifges(7,5) * t460 - Ifges(5,2) * t451 - Ifges(5,6) * t446 + Ifges(6,6) * t457 + Ifges(7,6) * t462 + Ifges(6,3) * t452 + Ifges(7,3) * t455 - t505;
t449 = t323 / 0.2e1;
t447 = t236 / 0.2e1;
t441 = pkin(3) * t240;
t284 = pkin(3) * t297;
t437 = pkin(5) * t307;
t189 = -t272 * t384 + qJD(2) * t392 + (-qJD(2) * t308 - qJD(3) * t274) * t315;
t170 = -pkin(8) * t240 + t189;
t190 = -t255 * qJD(2) - qJD(3) * t216;
t171 = pkin(8) * t239 + t190;
t76 = qJD(4) * t489 + t444 * t170 + t314 * t171;
t336 = -t314 * t255 - t339 * t444;
t165 = qJD(4) * t336 - t239 * t444 - t314 * t240;
t209 = t255 * t444 - t314 * t339;
t166 = qJD(4) * t209 - t314 * t239 + t240 * t444;
t83 = pkin(4) * t166 - qJ(5) * t165 - qJD(5) * t209 + t441;
t32 = t307 * t83 + t309 * t76;
t425 = Ifges(4,4) * t236;
t424 = Ifges(5,4) * t323;
t416 = t211 * mrSges(4,3);
t415 = t302 * mrSges(3,3);
t414 = t303 * mrSges(3,3);
t68 = mrSges(6,1) * t121 - mrSges(6,3) * t106;
t413 = t307 * t68;
t409 = t165 * t307;
t408 = t165 * t309;
t403 = t209 * t307;
t402 = t209 * t309;
t275 = t286 * qJ(5);
t396 = t294 * t316;
t395 = t294 * t319;
t394 = t296 * t316;
t393 = t296 * t319;
t147 = -pkin(4) * t336 - qJ(5) * t209 + t226;
t153 = t314 * t192 + t193 * t444;
t88 = t307 * t147 + t309 * t153;
t386 = t287 * pkin(4) + t275;
t380 = qJDD(1) * t308;
t379 = qJDD(1) * t310;
t378 = Ifges(7,5) * t42 + Ifges(7,6) * t43 + Ifges(7,3) * t115;
t377 = t444 * pkin(3);
t12 = -t43 * mrSges(7,1) + t42 * mrSges(7,2);
t31 = -t307 * t76 + t309 * t83;
t364 = -t207 * mrSges(4,1) + t206 * mrSges(4,2);
t363 = t121 * mrSges(5,1) + t120 * mrSges(5,2);
t65 = -t105 * mrSges(6,1) + t106 * mrSges(6,2);
t87 = t309 * t147 - t153 * t307;
t136 = t182 * t314 + t176;
t291 = -t377 - pkin(4);
t355 = -mrSges(3,1) * t379 + mrSges(3,2) * t380;
t350 = mrSges(5,1) * t286 + mrSges(5,2) * t287;
t343 = -t307 * t72 + t309 * t73;
t61 = -pkin(5) * t336 - pkin(9) * t402 + t87;
t75 = -pkin(9) * t403 + t88;
t26 = -t313 * t75 + t317 * t61;
t27 = t313 * t61 + t317 * t75;
t329 = mrSges(4,3) * t334;
t77 = qJD(4) * t153 + t314 * t170 - t444 * t171;
t293 = -qJDD(1) * pkin(1) + qJDD(2);
t270 = t291 - t432;
t265 = qJ(5) * t398;
t264 = qJ(5) * t399;
t262 = -qJD(1) * t290 + qJD(2);
t260 = t284 + t290;
t241 = t319 * t260;
t232 = Ifges(4,4) * t334;
t222 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t236;
t221 = -qJD(3) * mrSges(4,2) - t329;
t220 = t287 * t393 + t396;
t219 = -t287 * t395 + t394;
t218 = -t287 * t394 + t395;
t217 = t287 * t396 + t393;
t197 = t236 * Ifges(4,1) + Ifges(4,5) * qJD(3) - t232;
t196 = -Ifges(4,2) * t334 + Ifges(4,6) * qJD(3) + t425;
t194 = Ifges(5,4) * t199;
t185 = -mrSges(5,2) * t306 - mrSges(5,3) * t199;
t161 = t340 * t209;
t160 = t254 * t209;
t159 = mrSges(5,1) * t199 + mrSges(5,2) * t323;
t155 = Ifges(5,1) * t323 - t194 + t520;
t154 = -t199 * Ifges(5,2) + t424 + t519;
t108 = -mrSges(5,2) * t300 - mrSges(5,3) * t121;
t107 = mrSges(5,1) * t300 - mrSges(5,3) * t120;
t103 = pkin(5) * t403 - t489;
t93 = t136 - t532;
t92 = mrSges(7,1) * t195 - mrSges(7,3) * t134;
t91 = -mrSges(7,2) * t195 + mrSges(7,3) * t509;
t90 = t129 - t532;
t70 = -t165 * t254 + t209 * t237;
t69 = -t165 * t340 - t209 * t238;
t67 = -mrSges(6,2) * t121 + mrSges(6,3) * t105;
t55 = pkin(5) * t409 + t77;
t29 = -mrSges(7,2) * t115 + mrSges(7,3) * t43;
t28 = mrSges(7,1) * t115 - mrSges(7,3) * t42;
t24 = -pkin(9) * t409 + t32;
t19 = pkin(5) * t166 - pkin(9) * t408 + t31;
t5 = -qJD(6) * t27 + t19 * t317 - t24 * t313;
t4 = qJD(6) * t26 + t19 * t313 + t24 * t317;
t1 = [-(Ifges(6,5) * t106 + Ifges(6,6) * t105 + Ifges(6,3) * t121 + t378) * t336 / 0.2e1 + (mrSges(4,1) * t261 - mrSges(4,3) * t163 - Ifges(4,4) * t206 - Ifges(4,2) * t207 - Ifges(4,6) * qJDD(3)) * t339 + (mrSges(4,2) * t261 - mrSges(4,3) * t164 + Ifges(4,1) * t206 + Ifges(4,4) * t207 + Ifges(4,5) * qJDD(3)) * t255 + (t184 * mrSges(5,2) - t50 * mrSges(5,3) + Ifges(5,1) * t120 - Ifges(5,4) * t121 + Ifges(5,5) * t300 + t345 * t463 + t346 * t466 + t347 * t465 + t348 * t47) * t209 + (-m(6) * t241 - t220 * mrSges(7,1) - t219 * mrSges(7,2) + (-m(7) - m(5)) * (-t301 * t316 + t241) + (-m(7) * t437 + t507) * t316 + (-t523 - (m(6) * qJ(5) + mrSges(6,3)) * t286 + t506) * t319) * g(2) + (Ifges(5,4) * t452 + Ifges(5,1) * t449 + t345 * t451 + t360 * t346 / 0.2e1 + t155 / 0.2e1 + t520 / 0.2e1 + t180 * t347 / 0.2e1 + t510 + t511) * t165 + t385 * t279 * mrSges(3,3) + (-Ifges(7,5) * t161 - Ifges(7,6) * t160) * t464 + (-Ifges(7,4) * t161 - Ifges(7,2) * t160) * t471 + (-Ifges(7,1) * t161 - Ifges(7,4) * t160) * t472 + (t500 / 0.2e1 - Ifges(5,2) * t452 + Ifges(7,3) * t454 + Ifges(7,5) * t459 + Ifges(7,6) * t461 - Ifges(5,4) * t449 + Ifges(6,3) * t451 + t522 / 0.2e1 - t154 / 0.2e1 - t519 / 0.2e1 + t521 / 0.2e1 + t505) * t166 + (-t184 * mrSges(5,1) - t14 * mrSges(6,1) + t15 * mrSges(6,2) + t49 * mrSges(5,3) + Ifges(5,4) * t120 - Ifges(6,5) * t465 - Ifges(7,5) * t472 - Ifges(5,2) * t121 + Ifges(5,6) * t300 - Ifges(6,6) * t466 - Ifges(7,6) * t471 - Ifges(6,3) * t463 - Ifges(7,3) * t464 - t508) * t336 + (-t218 * mrSges(7,1) - t217 * mrSges(7,2) + (-m(7) * (-t301 + t437) + m(5) * t301 + t507) * t319 + (-m(7) * (-t260 - t488) - m(6) * (-t260 - t275) + t277 + m(5) * t260 - t506) * t316) * g(1) + (-t160 * t2 + t161 * t3 - t17 * t69 + t18 * t70) * mrSges(7,3) + t30 * (mrSges(7,1) * t160 - mrSges(7,2) * t161) + m(6) * (t14 * t87 + t15 * t88 + t31 * t72 + t32 * t73) + (-t14 * t402 - t15 * t403 - t408 * t72 - t409 * t73) * mrSges(6,3) + (Ifges(7,4) * t69 + Ifges(7,2) * t70) * t461 + m(4) * (t163 * t216 + t164 * t214 + t189 * t211 + t190 * t210 - t261 * t290) + m(7) * (t103 * t30 + t17 * t5 + t18 * t4 + t2 * t27 + t26 * t3 + t55 * t89) - (-Ifges(4,4) * t239 - Ifges(4,2) * t240) * t334 / 0.2e1 + (Ifges(7,1) * t69 + Ifges(7,4) * t70) * t459 + t214 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t206) + t216 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t207) + m(5) * (t129 * t76 + t153 * t49 + t184 * t226 + t212 * t441) + (-t525 + t528) * t77 + m(3) * (-pkin(1) * t293 + (t279 + t381) * qJ(2) * t385) + Ifges(2,3) * qJDD(1) + t159 * t441 + t76 * t185 + t153 * t108 + t32 * t145 + t31 * t146 - t37 * t403 / 0.2e1 - (-m(5) * t50 + m(6) * t47 - t107 + t65) * t489 + t103 * t12 - t240 * t416 + t279 * t414 + t279 * t415 + t210 * t239 * mrSges(4,3) + t4 * t91 + t5 * t92 + t87 * t68 + t88 * t67 + t89 * (-mrSges(7,1) * t70 + mrSges(7,2) * t69) + t55 * t78 + t226 * t363 - t290 * t364 + t27 * t29 + t26 * t28 + t69 * t467 + t70 * t469 + t402 * t473 - t161 * t474 - t160 * t475 + t189 * t221 + t190 * t222 - t239 * t197 / 0.2e1 - t240 * t196 / 0.2e1 + (Ifges(7,5) * t69 + Ifges(7,6) * t70) * t454 + (Ifges(3,4) * t308 + Ifges(3,2) * t310) * t379 + (Ifges(3,1) * t308 + Ifges(3,4) * t310) * t380 + (-Ifges(4,1) * t239 - Ifges(4,4) * t240) * t447 + qJD(3) * (-Ifges(4,5) * t239 - Ifges(4,6) * t240) / 0.2e1 + t262 * (mrSges(4,1) * t240 - mrSges(4,2) * t239) + t293 * t354 - pkin(1) * t355; t355 + m(6) * (t14 * t309 + t15 * t307) + m(3) * t293 + m(5) * t184 + t363 + t364 + t145 * t517 - t146 * t518 + t236 * t222 - t340 * t28 + t254 * t29 + t307 * t67 + t221 * t334 + t309 * t68 - t533 * t92 - t534 * t91 + (-t17 * t533 - t18 * t534 + t2 * t254 - t3 * t340) * m(7) + (t210 * t236 + t211 * t334 + t261) * m(4) + (-t385 * t443 - t414 - t415) * qJD(1) ^ 2 - (-m(5) * t129 - m(6) * t343 - t185) * t199 + (-t478 + t525) * t323 + (-g(1) * t316 + g(2) * t319) * (m(5) + m(4) + m(3) + t503); (-Ifges(4,2) * t236 + t197 - t232) * t334 / 0.2e1 + (-t329 - t221) * t210 - (Ifges(5,4) * t451 - t155 / 0.2e1 - t335 + t477) * t199 + (-m(5) * t284 - m(6) * (t284 + t386) - m(7) * (t284 + t488) - t353 + t480) * g(3) + Ifges(4,5) * t206 + Ifges(4,6) * t207 + (t128 * t136 - t129 * t137 - t212 * t442 + (t444 * t50 + t314 * t49 + (-t128 * t314 + t129 * t444) * qJD(4)) * pkin(3)) * m(5) + (t476 + t154 / 0.2e1 - Ifges(5,4) * t450) * t323 + t526 + (m(5) * t440 + mrSges(4,1) * t295 + mrSges(4,2) * t297 + t350) * t481 + t108 * t439 + t196 * t447 + t478 * pkin(3) * t383 + (t357 - t137) * t185 - t236 * (-Ifges(4,1) * t334 - t425) / 0.2e1 - t262 * (t236 * mrSges(4,1) - mrSges(4,2) * t334) - qJD(3) * (-Ifges(4,5) * t334 - Ifges(4,6) * t236) / 0.2e1 - t159 * t442 + t203 * t28 + t204 * t29 - t288 * t413 - t163 * mrSges(4,2) + t164 * mrSges(4,1) - t80 * t145 - t79 * t146 + t99 * t517 / 0.2e1 + t518 * t504 + t236 * t416 + t67 * t397 - t93 * t78 + Ifges(4,3) * qJDD(3) + t211 * t222 + t482 * t280 + (t316 * t486 + t483) * g(2) + (t319 * t486 + t484) * g(1) + t490 * t136 + t270 * t12 + t291 * t65 + t498 * t92 + t499 * t91 + (t17 * t498 + t18 * t499 + t2 * t204 + t203 * t3 + t270 * t30 - t89 * t93) * m(7) + t500 * t450 + t107 * t377 + (-g(1) * t265 - g(2) * t264 - t116 * t136 + t280 * t343 + t288 * t344 + t291 * t47 - t72 * t79 - t73 * t80) * m(6); (-m(6) * t386 + t480 - t523) * g(3) + (-t477 + t511) * t199 + (-t116 * t129 - t72 * t85 - t73 * t86 - g(1) * (-t319 * t438 + t265) - g(2) * (-t316 * t438 + t264) - pkin(4) * t47 + qJ(5) * t344 + qJD(5) * t343) * m(6) + t476 * t323 + t526 + (-t194 + t155) * t451 + t213 * t28 + t215 * t29 + t154 * t449 - t128 * t185 - qJ(5) * t413 - t86 * t145 - t85 * t146 + t67 * t412 - t90 * t78 - pkin(4) * t65 + t481 * t350 + t482 * qJD(5) + (t316 * t485 + t483) * g(2) + (t319 * t485 + t484) * g(1) + t490 * t129 - t289 * t12 + t496 * t92 + t497 * t91 + (t17 * t496 + t18 * t497 + t2 * t215 + t213 * t3 - t289 * t30 - t89 * t90) * m(7) + (-t424 + t500) * t450; t134 * t92 - t509 * t91 - t360 * t145 + t180 * t146 + t12 + t65 + (g(3) * t287 - t286 * t481) * t503 + (t134 * t17 - t18 * t509 + t30) * m(7) + (t180 * t72 - t360 * t73 + t47) * m(6); -t89 * (mrSges(7,1) * t134 + mrSges(7,2) * t509) + (Ifges(7,1) * t509 - t421) * t460 + t63 * t459 + (Ifges(7,5) * t509 - Ifges(7,6) * t134) * t455 - t17 * t91 + t18 * t92 - g(1) * (mrSges(7,1) * t219 - mrSges(7,2) * t220) - g(2) * (-mrSges(7,1) * t217 + mrSges(7,2) * t218) - g(3) * (-mrSges(7,1) * t294 - mrSges(7,2) * t296) * t286 + (t134 * t18 + t17 * t509) * mrSges(7,3) + t378 + (-Ifges(7,2) * t134 + t131 + t64) * t462 + t508;];
tau  = t1;
