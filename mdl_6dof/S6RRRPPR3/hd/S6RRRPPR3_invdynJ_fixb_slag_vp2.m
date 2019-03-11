% Calculate vector of inverse dynamics joint torques for
% S6RRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:28:35
% EndTime: 2019-03-09 15:28:59
% DurationCPUTime: 17.08s
% Computational Cost: add. (7998->706), mult. (17004->837), div. (0->0), fcn. (11302->10), ass. (0->315)
t508 = mrSges(4,1) + mrSges(5,1);
t507 = -mrSges(6,1) - mrSges(5,3);
t488 = Ifges(5,1) + Ifges(4,1);
t487 = Ifges(6,1) + Ifges(5,3);
t484 = -Ifges(4,5) - Ifges(5,4);
t506 = -Ifges(6,5) + Ifges(5,6);
t493 = -m(7) - m(6);
t505 = -m(7) * pkin(5) + t507;
t504 = -Ifges(4,6) + Ifges(5,6);
t282 = sin(qJ(3));
t283 = sin(qJ(2));
t432 = cos(qJ(3));
t433 = cos(qJ(2));
t213 = t282 * t433 + t283 * t432;
t193 = t213 * qJD(1);
t183 = Ifges(5,5) * t193;
t334 = t432 * t433;
t377 = qJD(1) * t283;
t192 = -qJD(1) * t334 + t282 * t377;
t275 = qJD(2) + qJD(3);
t412 = Ifges(6,4) * t193;
t503 = t192 * t487 + t275 * t506 + t183 - t412;
t212 = t282 * t283 - t334;
t371 = qJD(1) * qJD(2);
t502 = t433 * qJDD(1) - t283 * t371;
t277 = qJ(2) + qJ(3);
t272 = cos(t277);
t501 = (mrSges(6,2) - mrSges(7,3)) * t272;
t281 = sin(qJ(6));
t420 = mrSges(7,2) * t281;
t285 = cos(qJ(6));
t423 = mrSges(7,1) * t285;
t500 = -t420 + t423;
t162 = t192 * t285 - t275 * t281;
t185 = Ifges(4,4) * t192;
t408 = t192 * Ifges(5,5);
t188 = qJD(6) + t193;
t481 = t188 * Ifges(7,3);
t161 = -t192 * t281 - t275 * t285;
t482 = t161 * Ifges(7,6);
t499 = t162 * Ifges(7,5) + t193 * t488 - t275 * t484 - t185 + t408 + t481 + t482;
t271 = sin(t277);
t498 = -t508 * t272 + (mrSges(4,2) + t507) * t271;
t349 = qJD(1) * t433;
t223 = -qJD(1) * pkin(1) - pkin(2) * t349;
t117 = t192 * pkin(3) - t193 * qJ(4) + t223;
t316 = qJD(5) - t117;
t68 = -pkin(4) * t192 + t316;
t497 = -t223 * mrSges(4,1) - t117 * mrSges(5,1) - t68 * mrSges(6,2);
t320 = Ifges(7,5) * t285 - Ifges(7,6) * t281;
t410 = Ifges(7,4) * t285;
t321 = -Ifges(7,2) * t281 + t410;
t411 = Ifges(7,4) * t281;
t322 = Ifges(7,1) * t285 - t411;
t323 = mrSges(7,1) * t281 + mrSges(7,2) * t285;
t443 = -t188 / 0.2e1;
t445 = -t162 / 0.2e1;
t446 = -t161 / 0.2e1;
t409 = t162 * Ifges(7,4);
t62 = t161 * Ifges(7,2) + t188 * Ifges(7,6) + t409;
t270 = t275 * qJ(4);
t287 = -pkin(8) - pkin(7);
t224 = t287 * t283;
t215 = qJD(1) * t224;
t203 = qJD(2) * pkin(2) + t215;
t367 = t433 * pkin(7);
t225 = pkin(8) * t433 + t367;
t216 = t225 * qJD(1);
t354 = t432 * t216;
t146 = t282 * t203 + t354;
t402 = qJ(5) * t192;
t98 = t146 + t402;
t80 = -t270 - t98;
t69 = pkin(5) * t275 - t80;
t496 = t320 * t443 + t321 * t446 + t322 * t445 + t281 * t62 / 0.2e1 - t323 * t69;
t452 = -pkin(4) - pkin(9);
t44 = pkin(5) * t193 + t192 * t452 + t316;
t276 = -pkin(3) + t452;
t194 = t282 * t216;
t382 = -t432 * t203 + t194;
t397 = t193 * qJ(5);
t97 = t382 - t397;
t311 = qJD(4) + t97;
t66 = t275 * t276 + t311;
t22 = -t281 * t66 + t285 * t44;
t23 = t281 * t44 + t285 * t66;
t495 = t68 * mrSges(6,1) + t22 * mrSges(7,1) + t223 * mrSges(4,2) - t23 * mrSges(7,2) - t117 * mrSges(5,3);
t494 = -m(5) - m(6);
t348 = qJD(2) * t433;
t370 = t283 * qJDD(1);
t300 = qJD(1) * t348 + t370;
t113 = qJD(3) * t192 - t282 * t502 - t432 * t300;
t492 = t113 / 0.2e1;
t294 = t213 * qJD(3);
t114 = t282 * t370 - t432 * t502 + (t282 * t348 + t294) * qJD(1);
t491 = -t114 / 0.2e1;
t274 = qJDD(2) + qJDD(3);
t490 = -t274 / 0.2e1;
t489 = t502 / 0.2e1;
t486 = -Ifges(4,4) + Ifges(5,5);
t485 = Ifges(6,4) - Ifges(5,5);
t86 = -mrSges(6,1) * t274 - mrSges(6,3) * t114;
t88 = -mrSges(5,2) * t114 + mrSges(5,3) * t274;
t483 = t88 - t86;
t417 = mrSges(6,3) * t192;
t403 = -mrSges(6,1) * t275 + mrSges(7,1) * t161 - mrSges(7,2) * t162 - t417;
t418 = mrSges(4,3) * t193;
t421 = mrSges(5,2) * t193;
t479 = t275 * t508 - t418 - t421;
t416 = mrSges(6,3) * t193;
t169 = mrSges(6,2) * t275 - t416;
t115 = -mrSges(7,2) * t188 + mrSges(7,3) * t161;
t116 = mrSges(7,1) * t188 - mrSges(7,3) * t162;
t315 = t115 * t285 - t116 * t281;
t478 = t169 + t315;
t419 = mrSges(4,3) * t192;
t166 = -mrSges(4,2) * t275 - t419;
t422 = mrSges(5,2) * t192;
t170 = mrSges(5,3) * t275 - t422;
t477 = -t170 - t166;
t286 = cos(qJ(1));
t389 = t272 * t286;
t391 = t271 * t286;
t476 = pkin(3) * t389 + qJ(4) * t391;
t473 = t271 * pkin(5) + t272 * pkin(9);
t428 = pkin(7) * t283;
t472 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t377) * t367 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t349) * t428;
t209 = t502 * pkin(7);
t210 = t300 * pkin(7);
t471 = t433 * t209 + t210 * t283;
t108 = qJDD(6) - t113;
t51 = qJD(6) * t161 + t114 * t285 - t274 * t281;
t28 = mrSges(7,1) * t108 - mrSges(7,3) * t51;
t52 = -qJD(6) * t162 - t114 * t281 - t274 * t285;
t29 = -mrSges(7,2) * t108 + mrSges(7,3) * t52;
t470 = -t281 * t28 + t285 * t29;
t284 = sin(qJ(1));
t469 = g(1) * t286 + g(2) * t284;
t310 = mrSges(3,1) * t283 + mrSges(3,2) * t433;
t414 = Ifges(3,4) * t283;
t468 = pkin(1) * t310 - t283 * (Ifges(3,1) * t433 - t414) / 0.2e1;
t390 = t272 * t284;
t226 = qJ(4) * t390;
t365 = t272 * t423;
t392 = t271 * t284;
t467 = -m(7) * t226 - mrSges(6,2) * t392 - t284 * t365 + t390 * t505;
t228 = qJ(4) * t389;
t466 = -m(7) * t228 - mrSges(6,2) * t391 - t286 * t365 + t389 * t505;
t327 = m(7) * t276 - mrSges(7,3);
t293 = t271 * t327 - t272 * t420;
t288 = -pkin(3) - pkin(4);
t359 = t288 * t271;
t430 = pkin(2) * t283;
t465 = m(7) * t430 - t293 - m(5) * (-pkin(3) * t271 - t430) + t271 * mrSges(5,1) - m(6) * (t359 - t430);
t134 = t270 + t146;
t464 = m(5) * t134 + t170;
t463 = -t271 * t500 + t498 + t501;
t158 = qJDD(2) * pkin(2) - pkin(8) * t300 - t210;
t160 = pkin(8) * t502 + t209;
t347 = qJD(3) * t432;
t375 = qJD(3) * t282;
t41 = t158 * t432 - t282 * t160 - t203 * t375 - t216 * t347;
t297 = qJDD(4) - t41;
t291 = t113 * qJ(5) - t193 * qJD(5) + t297;
t13 = t274 * t276 + t291;
t278 = qJDD(1) * pkin(1);
t182 = -pkin(2) * t502 - t278;
t30 = t114 * pkin(3) + t113 * qJ(4) - t193 * qJD(4) + t182;
t307 = qJDD(5) - t30;
t6 = -pkin(5) * t113 + t114 * t452 + t307;
t2 = qJD(6) * t22 + t13 * t285 + t281 * t6;
t3 = -qJD(6) * t23 - t13 * t281 + t285 * t6;
t462 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t118 = -pkin(3) * t275 + qJD(4) + t382;
t461 = m(5) * t118 - t479;
t222 = -mrSges(3,1) * t433 + t283 * mrSges(3,2);
t460 = -m(3) * pkin(1) - mrSges(2,1) + t222 + t498;
t459 = m(6) * t80 - m(7) * t69 + t403;
t319 = t22 * t285 + t23 * t281;
t424 = t2 * t285;
t292 = -qJD(6) * t319 - t281 * t3 + t424;
t373 = qJD(6) * t285;
t374 = qJD(6) * t281;
t458 = m(7) * t292 - t115 * t374 - t116 * t373 + t470;
t457 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(5,2) - mrSges(3,3) - mrSges(4,3) + mrSges(6,3) + t493 * (-qJ(5) - t287);
t456 = -t459 + t464;
t454 = t51 / 0.2e1;
t453 = t52 / 0.2e1;
t451 = t108 / 0.2e1;
t450 = -t113 / 0.2e1;
t449 = t114 / 0.2e1;
t444 = t162 / 0.2e1;
t442 = -t192 / 0.2e1;
t441 = t192 / 0.2e1;
t440 = -t193 / 0.2e1;
t439 = t193 / 0.2e1;
t436 = t274 / 0.2e1;
t435 = -t275 / 0.2e1;
t434 = t275 / 0.2e1;
t431 = pkin(2) * t282;
t429 = pkin(4) * t193;
t425 = g(3) * t272;
t259 = t272 * pkin(3);
t280 = qJ(4) + pkin(5);
t273 = t433 * pkin(2);
t262 = t273 + pkin(1);
t415 = mrSges(7,3) * t281;
t413 = Ifges(4,4) * t193;
t155 = qJD(2) * t213 + t294;
t399 = t155 * t281;
t398 = t155 * t285;
t396 = t193 * t285;
t394 = t212 * t281;
t393 = t212 * t285;
t249 = t271 * qJ(4);
t388 = t281 * t284;
t387 = t281 * t286;
t385 = t284 * t285;
t384 = t285 * t286;
t138 = t193 * pkin(3) + t192 * qJ(4);
t378 = t259 + t249;
t376 = qJD(2) * t283;
t369 = Ifges(7,5) * t51 + Ifges(7,6) * t52 + Ifges(7,3) * t108;
t366 = t432 * pkin(2);
t364 = pkin(2) * t377;
t363 = pkin(2) * t376;
t360 = Ifges(3,4) * t433;
t258 = t272 * pkin(4);
t357 = t258 + t378;
t356 = t273 + t378;
t345 = -t373 / 0.2e1;
t255 = qJ(4) + t431;
t343 = -t113 * mrSges(6,1) + t114 * mrSges(6,2);
t340 = -t262 - t249;
t151 = t215 * t282 + t354;
t163 = -t432 * t224 + t225 * t282;
t229 = t286 * t262;
t339 = -t284 * t287 + t229;
t338 = -t212 * pkin(3) + t213 * qJ(4) + t262;
t336 = pkin(2) * t347;
t261 = -t366 - pkin(3);
t331 = g(1) * t284 - g(2) * t286;
t248 = -pkin(4) + t261;
t325 = mrSges(4,1) * t271 + mrSges(4,2) * t272;
t318 = t22 * t281 - t23 * t285;
t317 = t357 + t473;
t132 = -qJ(5) * t213 + t163;
t59 = pkin(5) * t213 + t212 * t452 + t338;
t32 = t132 * t285 + t281 * t59;
t31 = -t132 * t281 + t285 * t59;
t314 = -t115 * t281 - t116 * t285;
t119 = t364 + t138;
t309 = t433 * Ifges(3,2) + t414;
t308 = Ifges(3,5) * t433 - Ifges(3,6) * t283;
t152 = t215 * t432 - t194;
t164 = t282 * t224 + t225 * t432;
t306 = t212 * t373 + t399;
t305 = t212 * t374 - t398;
t154 = t212 * t275;
t56 = t155 * pkin(3) + t154 * qJ(4) - t213 * qJD(4) + t363;
t40 = t282 * t158 + t432 * t160 + t203 * t347 - t216 * t375;
t217 = qJD(2) * t224;
t296 = qJD(2) * t225;
t71 = t432 * t217 + t224 * t347 - t225 * t375 - t282 * t296;
t53 = -pkin(5) * t192 + t193 * t452 - t138;
t35 = t274 * qJ(4) + t275 * qJD(4) + t40;
t72 = qJD(3) * t164 + t282 * t217 + t432 * t296;
t19 = -qJ(5) * t114 - qJD(5) * t192 - t35;
t10 = t51 * Ifges(7,1) + t52 * Ifges(7,4) + t108 * Ifges(7,5);
t184 = Ifges(6,4) * t192;
t126 = -t193 * Ifges(6,2) - t275 * Ifges(6,6) + t184;
t127 = -Ifges(4,2) * t192 + Ifges(4,6) * t275 + t413;
t14 = pkin(5) * t274 - t19;
t16 = t274 * t288 + t291;
t36 = -t274 * pkin(3) + t297;
t159 = Ifges(7,4) * t161;
t63 = Ifges(7,1) * t162 + Ifges(7,5) * t188 + t159;
t67 = t275 * t288 + t311;
t9 = t51 * Ifges(7,4) + t52 * Ifges(7,2) + t108 * Ifges(7,6);
t289 = (-Ifges(6,5) + t504) * t114 + (Ifges(6,5) * t434 + t23 * t415 + t435 * t504 + t496 + t497) * t193 + (-Ifges(6,6) + t484) * t113 + (-Ifges(7,5) * t445 + Ifges(6,2) * t439 + Ifges(6,6) * t434 - Ifges(7,6) * t446 - Ifges(7,3) * t443 + t435 * t484 + t495) * t192 + (t193 * t487 + t126 + t184 - t408) * t442 + (-t192 * t488 + t183 - t413 + t503) * t440 + (-Ifges(4,2) * t193 - t185 + t499) * t441 + t14 * t500 + t496 * qJD(6) + (t23 * t374 - t424 + (t373 + t396) * t22) * mrSges(7,3) + (t412 + t127) * t439 + (-t396 / 0.2e1 + t345) * t63 - t40 * mrSges(4,2) + t41 * mrSges(4,1) - t36 * mrSges(5,1) + (-Ifges(7,5) * t281 - Ifges(7,6) * t285) * t451 + (-Ifges(7,2) * t285 - t411) * t453 + (-Ifges(7,1) * t281 - t410) * t454 + t35 * mrSges(5,3) - t19 * mrSges(6,1) + t16 * mrSges(6,2) - t285 * t9 / 0.2e1 + t382 * t419 + (Ifges(5,2) + Ifges(4,3) + Ifges(6,3)) * t274 + t3 * t415 + t80 * t416 + t146 * t418 + t134 * t421 + t118 * t422 - t281 * t10 / 0.2e1 - t67 * t417;
t264 = Ifges(3,4) * t349;
t247 = pkin(5) + t255;
t191 = Ifges(3,1) * t377 + Ifges(3,5) * qJD(2) + t264;
t190 = Ifges(3,6) * qJD(2) + qJD(1) * t309;
t181 = t271 * t384 - t388;
t180 = -t271 * t387 - t385;
t179 = -t271 * t385 - t387;
t178 = t271 * t388 - t384;
t143 = mrSges(4,1) * t192 + mrSges(4,2) * t193;
t141 = mrSges(5,1) * t192 - mrSges(5,3) * t193;
t137 = mrSges(6,1) * t193 + mrSges(6,2) * t192;
t112 = t152 + t397;
t111 = t151 + t402;
t110 = -pkin(4) * t212 + t338;
t102 = t113 * mrSges(5,2);
t101 = t113 * mrSges(6,3);
t95 = -t138 - t429;
t87 = -mrSges(4,2) * t274 - mrSges(4,3) * t114;
t85 = -t274 * mrSges(5,1) - t102;
t84 = mrSges(4,1) * t274 + mrSges(4,3) * t113;
t83 = t274 * mrSges(6,2) + t101;
t79 = -t119 - t429;
t48 = t53 - t364;
t38 = t154 * qJ(5) - t213 * qJD(5) + t72;
t33 = -pkin(4) * t155 - t56;
t27 = t281 * t53 + t285 * t98;
t26 = -t281 * t98 + t285 * t53;
t25 = t111 * t285 + t281 * t48;
t24 = -t111 * t281 + t285 * t48;
t18 = -pkin(5) * t154 + t155 * t452 - t56;
t17 = -mrSges(7,1) * t52 + mrSges(7,2) * t51;
t12 = -pkin(4) * t114 + t307;
t5 = -qJD(6) * t32 + t18 * t285 - t281 * t38;
t4 = qJD(6) * t31 + t18 * t281 + t285 * t38;
t1 = [m(4) * (-t182 * t262 + t223 * t363) + (Ifges(3,1) * t300 + Ifges(3,4) * t489) * t283 + (Ifges(6,6) * t274 - Ifges(6,2) * t113 + t12 * mrSges(6,1) + Ifges(4,4) * t491 - t30 * mrSges(5,3) + t182 * mrSges(4,2) + Ifges(5,5) * t449 + Ifges(7,3) * t451 + Ifges(7,6) * t453 + Ifges(7,5) * t454 + t36 * mrSges(5,2) - t41 * mrSges(4,3) - t16 * mrSges(6,3) + t488 * t450 - t484 * t436 + (t491 - t449) * Ifges(6,4) + t462) * t213 + (Ifges(4,2) * t114 + Ifges(6,5) * t490 + Ifges(6,4) * t492 + t12 * mrSges(6,2) + t30 * mrSges(5,1) + t182 * mrSges(4,1) + Ifges(5,5) * t450 + t320 * t451 + t321 * t453 + t322 * t454 + t14 * t323 + t62 * t345 + Ifges(5,6) * t436 - t19 * mrSges(6,3) - t35 * mrSges(5,2) - t40 * mrSges(4,3) + t487 * t449 + (t490 - t436) * Ifges(4,6) + (t492 - t450) * Ifges(4,4)) * t212 + (-t488 * t113 + t486 * t114 - t484 * t274 + t369) * t213 / 0.2e1 + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t471) + (t308 * qJD(2) / 0.2e1 - t472) * qJD(2) - t468 * t371 + (m(4) * t146 + t166 + t464) * t71 + (-t179 * mrSges(7,1) - t178 * mrSges(7,2) + ((m(4) + m(5)) * t287 + t457) * t286 + (-m(6) * t340 - m(5) * (t340 - t259) - m(7) * (-t271 * t280 - t262) + m(4) * t262 + (-m(6) * t288 - mrSges(6,2) - t327) * t272 - t460) * t284) * g(1) + m(7) * (t2 * t32 + t22 * t5 + t23 * t4 + t3 * t31) + m(6) * (t110 * t12 + t132 * t16 + t33 * t68 + t38 * t67) + t161 * (-Ifges(7,4) * t305 - Ifges(7,2) * t306) / 0.2e1 + (-t181 * mrSges(7,1) - t180 * mrSges(7,2) - m(5) * (t339 + t476) - m(4) * t339 + t493 * (pkin(4) * t389 + t229 + t476) + t457 * t284 + (-m(7) * t473 + t460 + t501) * t286) * g(2) + (-t459 * qJ(5) - t127 / 0.2e1 - Ifges(4,2) * t442 + Ifges(6,5) * t435 + Ifges(6,4) * t440 - mrSges(6,3) * t80 - mrSges(5,2) * t134 - mrSges(4,3) * t146 + t487 * t441 + t486 * t439 + t504 * t434 - t497 + t503 / 0.2e1) * t155 + (-Ifges(7,1) * t305 - Ifges(7,4) * t306) * t444 + (-m(6) * t19 + m(7) * t14 + t17 - t86) * (t212 * qJ(5) + t164) + (qJD(1) * (-Ifges(3,2) * t283 + t360) + t191) * t348 / 0.2e1 - (qJD(6) * t63 + t9) * t394 / 0.2e1 + (-t499 / 0.2e1 + t126 / 0.2e1 - Ifges(4,4) * t442 - Ifges(7,5) * t444 - t482 / 0.2e1 - t481 / 0.2e1 + Ifges(6,6) * t435 + Ifges(6,2) * t440 + t67 * mrSges(6,3) - mrSges(5,2) * t118 - mrSges(4,3) * t382 + t485 * t441 - t488 * t439 + t484 * t434 - t495) * t154 + t300 * t360 / 0.2e1 + (-t2 * t394 + t22 * t305 - t23 * t306 - t3 * t393) * mrSges(7,3) + (-m(4) * t41 + m(5) * t36 - t84 + t85) * t163 + t459 * (-t212 * qJD(5) - t71) + (t300 * t428 + t367 * t502 + t471) * mrSges(3,3) + t433 * (Ifges(3,4) * t300 + Ifges(3,2) * t502) / 0.2e1 - pkin(1) * (-mrSges(3,1) * t502 + mrSges(3,2) * t300) + t309 * t489 + (m(4) * t40 + m(5) * t35 + t87 + t88) * t164 - t262 * (mrSges(4,1) * t114 - mrSges(4,2) * t113) + t188 * (-Ifges(7,5) * t305 - Ifges(7,6) * t306) / 0.2e1 - t222 * t278 + t63 * t398 / 0.2e1 - t62 * t399 / 0.2e1 + t10 * t393 / 0.2e1 - t190 * t376 / 0.2e1 + Ifges(2,3) * qJDD(1) + t38 * t169 + t33 * t137 + t56 * t141 + t132 * t83 + t4 * t115 + t5 * t116 + t110 * t343 + t32 * t29 + t31 * t28 + (-mrSges(3,1) * t428 - mrSges(3,2) * t367 + Ifges(3,5) * t283 + Ifges(3,6) * t433) * qJDD(2) + (t485 * t113 + t487 * t114 + t506 * t274) * t212 / 0.2e1 + t143 * t363 + m(5) * (t117 * t56 - t30 * t338) - t338 * (mrSges(5,1) * t114 + mrSges(5,3) * t113) + (m(4) * t382 + t461) * t72 + t69 * (mrSges(7,1) * t306 - mrSges(7,2) * t305); (t228 * t494 + t465 * t286 + t466) * g(1) + (t226 * t494 + t465 * t284 + t467) * g(2) + t483 * t255 + t403 * t112 + t477 * t152 + t479 * t151 + (t222 - m(5) * t356 - m(4) * t273 - m(6) * (t258 + t356) - m(7) * (t273 + t317) + t463) * g(3) + (m(6) * t67 - m(7) * t318 + t461 + t478) * pkin(2) * t375 + t458 * (-pkin(9) + t248) + t456 * (t336 + qJD(4)) + (m(4) * t430 + t310 + t325) * t469 - (-Ifges(3,2) * t377 + t191 + t264) * t349 / 0.2e1 + (qJD(1) * t468 + t472) * qJD(1) + (-t112 * t69 + t14 * t247 - t22 * t24 - t23 * t25) * m(7) + (-t111 * t67 + t112 * t80 + t16 * t248 - t19 * t255 - t68 * t79) * m(6) + (-t117 * t119 - t118 * t151 - t134 * t152 + t255 * t35 + t261 * t36) * m(5) + Ifges(3,5) * t300 + Ifges(3,6) * t502 + t166 * t336 + t261 * t85 + t247 * t17 + t248 * t83 - t209 * mrSges(3,2) - t210 * mrSges(3,1) + t190 * t377 / 0.2e1 - t308 * t371 / 0.2e1 - t111 * t169 - t143 * t364 - t79 * t137 - t119 * t141 - t25 * t115 - t24 * t116 + Ifges(3,3) * qJDD(2) + (-t382 * t151 - t146 * t152 - t223 * t364 + (t432 * t41 + t282 * t40 + (t146 * t432 + t282 * t382) * qJD(3)) * pkin(2)) * m(4) + t84 * t366 + t87 * t431 + t289; -pkin(3) * t85 - t27 * t115 - t26 * t116 - t95 * t137 - t138 * t141 - t98 * t169 + t280 * t17 + t288 * t83 + t289 - t403 * t97 + t469 * t325 + t479 * t146 - t477 * t382 + t483 * qJ(4) + (t14 * t280 - t22 * t26 - t23 * t27 + t69 * t97) * m(7) + (-t19 * qJ(4) + t16 * t288 - t67 * t98 - t68 * t95 - t80 * t97) * m(6) + (-pkin(3) * t36 + qJ(4) * t35 - t117 * t138 - t118 * t146 + t134 * t382) * m(5) + t456 * qJD(4) + (-m(5) * (-pkin(3) * t392 + t226) + mrSges(5,1) * t392 - m(6) * (t284 * t359 + t226) - t293 * t284 + t467) * g(2) + (-m(5) * (-pkin(3) * t391 + t228) + mrSges(5,1) * t391 - m(6) * (t286 * t359 + t228) - t293 * t286 + t466) * g(1) + t458 * t276 + (-m(5) * t378 - m(6) * t357 - m(7) * t317 + t463) * g(3); t101 - t102 + (-mrSges(5,1) + mrSges(6,2)) * t274 + t314 * qJD(6) + (-t170 + t403) * t275 + (-t137 + t141 + t314) * t193 + (-t193 * t319 - t275 * t69 + t292) * m(7) + (-t193 * t68 + t275 * t80 + t16) * m(6) + (t117 * t193 - t134 * t275 + t36) * m(5) + (-t271 * t469 + t425) * (m(5) - t493) + t470; t285 * t28 + t281 * t29 + t403 * t192 + t315 * qJD(6) + t478 * t193 + t343 + (-t188 * t318 - t192 * t69 + t2 * t281 + t285 * t3 + t331) * m(7) + (t192 * t80 + t193 * t67 + t12 + t331) * m(6); -t69 * (mrSges(7,1) * t162 + mrSges(7,2) * t161) + (Ifges(7,1) * t161 - t409) * t445 + t62 * t444 + (Ifges(7,5) * t161 - Ifges(7,6) * t162) * t443 - t22 * t115 + t23 * t116 - g(1) * (mrSges(7,1) * t180 - mrSges(7,2) * t181) - g(2) * (-mrSges(7,1) * t178 + mrSges(7,2) * t179) - t323 * t425 + (t161 * t22 + t162 * t23) * mrSges(7,3) + t369 + (-Ifges(7,2) * t162 + t159 + t63) * t446 + t462;];
tau  = t1;
