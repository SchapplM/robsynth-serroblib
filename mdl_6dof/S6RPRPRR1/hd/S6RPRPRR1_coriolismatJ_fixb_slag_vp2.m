% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:33:57
% EndTime: 2019-03-09 03:34:07
% DurationCPUTime: 8.09s
% Computational Cost: add. (21929->449), mult. (40965->585), div. (0->0), fcn. (47497->10), ass. (0->246)
t269 = sin(qJ(6));
t272 = cos(qJ(6));
t251 = -mrSges(7,1) * t272 + mrSges(7,2) * t269;
t266 = sin(pkin(11));
t267 = cos(pkin(11));
t271 = sin(qJ(3));
t326 = sin(pkin(10)) * pkin(1) + pkin(7);
t312 = qJ(4) + t326;
t300 = t271 * t312;
t273 = cos(qJ(3));
t452 = t273 * t312;
t208 = -t266 * t300 + t267 * t452;
t247 = -t266 * t271 + t267 * t273;
t170 = t247 * pkin(8) + t208;
t270 = sin(qJ(5));
t411 = cos(qJ(5));
t248 = -t266 * t273 - t267 * t271;
t462 = -t266 * t452 - t267 * t300;
t477 = t248 * pkin(8) + t462;
t488 = t170 * t411 + t270 * t477;
t491 = t488 * t251;
t493 = t488 * mrSges(6,1);
t117 = t170 * t270 - t411 * t477;
t499 = t117 * mrSges(6,2);
t500 = t491 - t493 + t499;
t259 = pkin(3) * t267 + pkin(4);
t409 = pkin(3) * t266;
t238 = t270 * t259 + t409 * t411;
t498 = t117 * t238;
t497 = t117 * t269;
t496 = t117 * t272;
t303 = t270 * t247 - t248 * t411;
t371 = t117 * t303;
t374 = t488 * t117;
t494 = t491 / 0.2e1 - t493 / 0.2e1;
t237 = t259 * t411 - t270 * t409;
t235 = -pkin(5) - t237;
t492 = t235 * t488;
t224 = t411 * t247 + t248 * t270;
t373 = t488 * t224;
t380 = t272 * mrSges(7,2);
t385 = t269 * mrSges(7,1);
t252 = t380 + t385;
t151 = t252 * t303;
t379 = t272 * mrSges(7,3);
t464 = mrSges(7,1) * t303 - t224 * t379;
t467 = t252 * t224;
t339 = -cos(pkin(10)) * pkin(1) - pkin(2);
t250 = -pkin(3) * t273 + t339;
t227 = -t247 * pkin(4) + t250;
t131 = -pkin(5) * t224 - pkin(9) * t303 + t227;
t50 = t131 * t272 - t269 * t488;
t490 = t117 * t467 + t488 * t151 + t50 * t464;
t430 = -t303 / 0.2e1;
t489 = 0.2e1 * t430;
t487 = t248 * mrSges(5,1) - t247 * mrSges(5,2);
t264 = t269 ^ 2;
t265 = t272 ^ 2;
t350 = t264 + t265;
t468 = t224 * t350;
t445 = mrSges(7,3) * t468;
t451 = t303 * t251;
t484 = t445 + t451;
t412 = t272 / 0.2e1;
t415 = -t269 / 0.2e1;
t429 = t303 / 0.2e1;
t260 = Ifges(7,5) * t272;
t397 = Ifges(7,6) * t269;
t450 = t260 - t397;
t261 = Ifges(7,4) * t272;
t449 = -t269 * Ifges(7,2) + t261;
t463 = Ifges(7,6) * t303 + t224 * t449;
t399 = Ifges(7,4) * t269;
t257 = Ifges(7,1) * t272 - t399;
t465 = Ifges(7,5) * t303 + t224 * t257;
t472 = -t224 / 0.2e1;
t483 = Ifges(6,4) * t489 + t465 * t412 + t463 * t415 + t227 * mrSges(6,1) + t450 * t429 + (Ifges(6,2) + Ifges(7,3)) * t472;
t307 = t380 / 0.2e1 + t385 / 0.2e1;
t296 = t307 * t224;
t479 = -t467 / 0.2e1 - t296;
t482 = t479 * qJD(6);
t481 = t463 / 0.4e1;
t480 = t465 / 0.4e1;
t432 = t467 / 0.2e1;
t340 = t451 / 0.2e1 + t445 / 0.2e1;
t438 = m(7) / 0.2e1;
t459 = pkin(5) * t303;
t475 = pkin(9) * t468 - t459;
t454 = t303 * mrSges(6,1);
t469 = t224 * mrSges(6,2);
t476 = -t454 / 0.2e1 - t469 / 0.2e1;
t308 = t438 * t475 + t340 + t476;
t478 = -t454 - t469;
t426 = -t224 / 0.4e1;
t471 = t224 / 0.2e1;
t423 = t224 / 0.4e1;
t367 = t224 * t303;
t253 = Ifges(7,5) * t269 + Ifges(7,6) * t272;
t466 = t253 * t303;
t160 = -pkin(9) * t224 + t459;
t457 = mrSges(7,2) * t303;
t325 = (t265 / 0.2e1 + t264 / 0.2e1) * mrSges(7,3);
t419 = t253 / 0.2e1;
t453 = t419 - Ifges(6,6);
t328 = t350 * t237;
t401 = Ifges(7,1) * t269;
t256 = t261 + t401;
t359 = t269 * t224;
t153 = -mrSges(7,3) * t359 - t457;
t448 = t153 * t412 + t415 * t464;
t114 = -Ifges(7,5) * t224 + t257 * t303;
t357 = t272 * t114;
t111 = -Ifges(7,6) * t224 + t303 * t449;
t362 = t269 * t111;
t447 = t357 / 0.4e1 - t362 / 0.4e1;
t62 = t160 * t272 + t497;
t63 = t160 * t269 - t496;
t319 = -t269 * t62 + t272 * t63;
t446 = Ifges(6,1) * t429 + Ifges(6,4) * t224 - t362 / 0.2e1 + t357 / 0.2e1;
t444 = (-t260 / 0.2e1 + t397 / 0.2e1) * t224;
t442 = 0.2e1 * m(7);
t441 = 2 * qJD(3);
t440 = m(5) / 0.2e1;
t439 = m(6) / 0.2e1;
t437 = pkin(5) / 0.2e1;
t436 = Ifges(7,3) / 0.2e1;
t435 = -t62 / 0.2e1;
t434 = t63 / 0.2e1;
t433 = t117 / 0.2e1;
t431 = t151 / 0.2e1;
t422 = -t235 / 0.2e1;
t236 = pkin(9) + t238;
t421 = -t236 / 0.2e1;
t420 = -t237 / 0.2e1;
t417 = -t256 / 0.4e1;
t416 = t260 / 0.4e1;
t414 = t269 / 0.2e1;
t413 = -t272 / 0.2e1;
t408 = pkin(5) * t467;
t405 = pkin(5) * t252;
t263 = t271 * pkin(3);
t398 = Ifges(7,2) * t272;
t389 = t236 * t62;
t388 = t236 * t63;
t387 = t237 * t50;
t51 = t131 * t269 + t272 * t488;
t386 = t237 * t51;
t384 = t269 * mrSges(7,3);
t231 = -pkin(4) * t248 + t263;
t132 = t160 + t231;
t56 = t132 * t272 + t497;
t382 = t269 * t56;
t378 = t272 * t51;
t57 = t132 * t269 - t496;
t377 = t272 * t57;
t349 = t303 * t384;
t154 = mrSges(7,2) * t224 - t349;
t157 = -mrSges(7,1) * t224 - t303 * t379;
t304 = t154 * t415 + t157 * t413;
t23 = t451 * t472 + (t303 * t325 - t304) * t303;
t375 = qJD(1) * t23;
t121 = t224 * t238 - t237 * t303;
t282 = (t269 * t57 + t272 * t56) * t438 + t153 * t414 + t464 * t412 + t263 * t440;
t311 = (t247 * t266 + t248 * t267) * pkin(3) * t440 + (t235 * t303 + t236 * t468) * t438;
t289 = t311 + t340;
t301 = t478 + t487;
t17 = 0.2e1 * (t121 / 0.4e1 - t231 / 0.4e1) * m(6) - t282 + t289 + t301;
t369 = t17 * qJD(1);
t152 = -t224 * t384 - t457;
t279 = -m(7) * (t269 * t63 + t272 * t62) / 0.2e1 + t152 * t415 + t464 * t413 + t476;
t21 = t279 + t308;
t368 = t21 * qJD(1);
t364 = t235 * t252;
t355 = t272 * t154;
t360 = t269 * t157;
t305 = t355 / 0.2e1 - t360 / 0.2e1;
t25 = -t296 - t305;
t363 = t25 * qJD(1);
t361 = t269 * t464;
t254 = t398 + t399;
t358 = t269 * t254;
t356 = t272 * t152;
t354 = t272 * t256;
t287 = m(7) * (t303 * t468 - t367);
t37 = t287 / 0.2e1;
t353 = t37 * qJD(1);
t348 = t121 * t439;
t347 = t436 + Ifges(6,2) / 0.2e1;
t343 = t260 / 0.2e1;
t341 = mrSges(6,1) / 0.2e1 - t251 / 0.2e1;
t338 = t252 * t433;
t243 = -t358 / 0.2e1;
t333 = t423 + t426;
t327 = t449 * t412 + t257 * t414 + t243 + t354 / 0.2e1;
t2 = t51 * t153 + t57 * t154 + t56 * t157 + (mrSges(4,1) * t339 - Ifges(4,4) * t271) * t271 + m(7) * (t50 * t56 + t51 * t57 + t374) + (-mrSges(5,2) * t263 - Ifges(5,4) * t248) * t248 + (t339 * mrSges(4,2) + Ifges(4,4) * t273 + (Ifges(4,1) - Ifges(4,2)) * t271) * t273 + (-mrSges(5,1) * t263 + Ifges(5,4) * t247 + (-Ifges(5,1) + Ifges(5,2)) * t248) * t247 + (t231 * mrSges(6,2) + Ifges(6,1) * t471 + t483) * t303 + (-t231 * mrSges(6,1) - t347 * t303 + t472 * t450 + t446) * t224 + (m(5) * t263 - t487) * t250 + t490 + (m(6) * t231 + t469) * t227;
t321 = -t269 * t50 + t378;
t286 = t224 * t321 + t371;
t320 = t377 - t382;
t7 = t151 * t429 + t467 * t472 + t448 * t303 + t305 * t224 + (t303 * t320 + t286 - t373) * t438;
t322 = t2 * qJD(1) + t7 * qJD(2);
t10 = (t356 / 0.2e1 - t361 / 0.2e1 + (t117 + t319) * t438 + t431) * t303 - ((t488 - t321) * t438 + t432 - t305) * t224;
t5 = t62 * t157 + t63 * t154 + t51 * t152 + m(7) * (t50 * t62 + t51 * t63 + t374) + t483 * t303 - (-t227 * mrSges(6,2) - t444 + (-Ifges(6,1) / 0.2e1 + t347) * t303 - t446) * t224 + t490;
t318 = t5 * qJD(1) + t10 * qJD(2);
t6 = t51 * t157 + t117 * t451 + (t114 * t414 + t111 * t412 + mrSges(7,3) * t378 - t224 * t419 + (-t256 * t413 + t243) * t303) * t303 + (-t349 - t154) * t50;
t317 = -t6 * qJD(1) - t23 * qJD(2);
t316 = t7 * qJD(1) + qJD(2) * t287;
t39 = m(7) * (0.1e1 - t350) * t367;
t315 = t10 * qJD(1) - t39 * qJD(2);
t13 = (mrSges(6,3) * t303 + t151) * t303 + (t247 ^ 2 + t248 ^ 2) * mrSges(5,3) + (t224 * mrSges(6,3) + t355 - t360) * t224 + m(7) * t286 + m(6) * (t371 + t373) + m(5) * (t208 * t247 + t248 * t462);
t314 = -qJD(1) * t13 - qJD(2) * t37;
t309 = mrSges(7,1) * t435 + mrSges(7,2) * t434;
t306 = t450 * t426;
t297 = -t354 / 0.2e1 + t358 / 0.2e1 - Ifges(6,5);
t275 = ((t235 + t328) * t438 - t341) * t303 - ((-t236 * t350 + t238) * t438 + mrSges(6,2) / 0.2e1 - t325) * t224;
t20 = t275 - t308;
t278 = (t492 + t498) * t438 + t235 * t432 + t238 * t431 + t494;
t292 = t157 * t420 + t464 * t421 + t480;
t293 = t481 + t236 * t152 / 0.2e1 + t237 * t154 / 0.2e1;
t3 = t408 / 0.2e1 + (t430 + t429) * Ifges(6,6) + (t471 + t472) * Ifges(6,5) + (-t117 / 0.2e1 + t433) * mrSges(6,2) + (m(7) * t437 + t341) * t488 + (-pkin(9) * t153 / 0.2e1 - t463 / 0.4e1 - t333 * t256 + (t434 - t57 / 0.2e1) * mrSges(7,3) + (-pkin(9) * t57 / 0.4e1 + t388 / 0.4e1 + t386 / 0.4e1) * t442 + t293) * t272 + (pkin(9) * t464 / 0.2e1 - t465 / 0.4e1 + t333 * t254 + (t435 + t56 / 0.2e1) * mrSges(7,3) + (pkin(9) * t56 / 0.4e1 - t389 / 0.4e1 - t387 / 0.4e1) * t442 + t292) * t269 + t278;
t283 = -t237 * mrSges(6,2) + mrSges(7,3) * t328 + (t251 - mrSges(6,1)) * t238;
t41 = m(7) * (t235 * t238 + t236 * t328) + t283;
t295 = t3 * qJD(1) + t20 * qJD(2) + t41 * qJD(3);
t133 = t327 + t364;
t73 = t432 - t296;
t280 = (t257 / 0.4e1 - t254 / 0.4e1 - t398 / 0.4e1) * t272 + (t417 - t449 / 0.4e1 - t261 / 0.2e1 - t401 / 0.4e1) * t269;
t274 = (-t236 * t325 + t280) * t303 + t338 + t451 * t422;
t288 = Ifges(7,3) * t430 - t56 * mrSges(7,1) / 0.2e1 + t57 * mrSges(7,2) / 0.2e1;
t9 = -t224 * t416 + (t154 * t421 - t111 / 0.4e1 + (t423 + t471) * Ifges(7,6)) * t269 + (t157 * t421 + t114 / 0.4e1 + Ifges(7,5) * t472) * t272 + t274 + t288;
t294 = t9 * qJD(1) - t73 * qJD(2) + t133 * qJD(3);
t276 = -pkin(9) * t325 + t280;
t277 = pkin(9) * t304 + t437 * t451 + t338 + t447;
t12 = -(-0.3e1 / 0.4e1 * t397 + t416 + t343) * t224 + (-Ifges(7,3) / 0.2e1 + t276) * t303 + t277 + t309;
t161 = t327 - t405;
t70 = (t422 + t437) * t252 + (mrSges(7,2) * t420 - t256 / 0.2e1 - t449 / 0.2e1) * t272 + (mrSges(7,1) * t420 - t257 / 0.2e1 + t254 / 0.2e1) * t269;
t72 = (-t252 / 0.2e1 + t307) * t224;
t285 = t12 * qJD(1) + t72 * qJD(2) - t70 * qJD(3) + t161 * qJD(5);
t71 = t364 / 0.2e1 - t405 / 0.2e1 - t307 * t237 + t327;
t26 = -t296 + t305;
t22 = -t279 + t308;
t19 = t275 + t308;
t18 = t231 * t439 + t282 + t289 + t348;
t11 = t277 + t306 - t309 + (t276 + t436) * t303 - t444;
t8 = t306 + t224 * t343 - Ifges(7,6) * t359 / 0.2e1 + t304 * t236 + t274 - t288 + t447;
t4 = t466 / 0.2e1 + t278 + t354 * t423 + t358 * t426 - t408 / 0.2e1 + t499 - pkin(5) * t488 * t438 + (t481 + (t386 + t388) * t438 - t224 * t417 + t293) * t272 + (t480 + (-t387 - t389) * t438 + t254 * t426 + t292) * t269 + (t320 * t438 + t448) * pkin(9) + 0.2e1 * t471 * Ifges(6,5) + (t377 / 0.2e1 - t382 / 0.2e1 + t435 * t269 + t434 * t272) * mrSges(7,3) + Ifges(6,6) * t489 + t494;
t1 = qJD(3) * t7 + qJD(4) * t37 + qJD(5) * t10 - qJD(6) * t23;
t14 = [qJD(3) * t2 + qJD(4) * t13 + qJD(5) * t5 - qJD(6) * t6, t1, t18 * qJD(4) + t4 * qJD(5) + t8 * qJD(6) + ((t236 * t320 + t492) * t438 + (-t237 * t488 - t498) * t439) * t441 + t322 + (-t208 * mrSges(5,1) - t462 * mrSges(5,2) + Ifges(5,5) * t247 + Ifges(5,6) * t248 + t235 * t467 + (t57 * mrSges(7,3) + t236 * t153 + t463 / 0.2e1) * t272 + (-t236 * t464 - t56 * mrSges(7,3) + t465 / 0.2e1) * t269 + (-t238 * mrSges(6,3) + t453) * t303 + (-t237 * mrSges(6,3) - t297) * t224 + (m(5) * (-t208 * t267 + t266 * t462) + (-t247 * t267 + t248 * t266) * mrSges(5,3)) * pkin(3) + (-mrSges(4,1) * t326 + Ifges(4,5)) * t273 + (mrSges(4,2) * t326 - Ifges(4,6)) * t271 + t500) * qJD(3), qJD(3) * t18 + qJD(5) * t22 + qJD(6) * t26 - t314, t4 * qJD(3) + t22 * qJD(4) + t11 * qJD(6) + t318 + (t463 * t412 + t465 * t414 - t297 * t224 + t453 * t303 + (-m(7) * t488 - t467) * pkin(5) + (m(7) * t319 + t356 - t361) * pkin(9) + t319 * mrSges(7,3) + t500) * qJD(5), t8 * qJD(3) + t26 * qJD(4) + t11 * qJD(5) + (-mrSges(7,1) * t51 - mrSges(7,2) * t50 - t466) * qJD(6) + t317; t1, qJD(3) * t287 - qJD(5) * t39 (-t271 * mrSges(4,1) - t273 * mrSges(4,2) + t301 + t484) * qJD(3) + t19 * qJD(5) + t482 + (t348 + t311) * t441 + t316, t353, t19 * qJD(3) + (m(7) * t475 + t478 + t484) * qJD(5) + t482 + t315, qJD(6) * t451 - t375 + (qJD(3) + qJD(5)) * t479; qJD(4) * t17 + qJD(5) * t3 + qJD(6) * t9 - t322, qJD(5) * t20 - qJD(6) * t73 - t316, qJD(5) * t41 + qJD(6) * t133, t369 (m(7) * (-pkin(5) * t238 + pkin(9) * t328) + t283) * qJD(5) + t71 * qJD(6) + t295, t71 * qJD(5) + (t236 * t251 + t450) * qJD(6) + t294; -qJD(3) * t17 - qJD(5) * t21 - qJD(6) * t25 + t314, -t353, -t369, 0, -t368, -qJD(6) * t252 - t363; -qJD(3) * t3 + qJD(4) * t21 + qJD(6) * t12 - t318, -qJD(3) * t20 + qJD(6) * t72 - t315, -qJD(6) * t70 - t295, t368, t161 * qJD(6) (pkin(9) * t251 + t450) * qJD(6) + t285; -qJD(3) * t9 + qJD(4) * t25 - qJD(5) * t12 - t317, qJD(3) * t73 - qJD(5) * t72 + t375, qJD(5) * t70 - t294, t363, -t285, 0;];
Cq  = t14;
