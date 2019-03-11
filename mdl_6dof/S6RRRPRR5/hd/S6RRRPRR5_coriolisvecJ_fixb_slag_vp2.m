% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:20:39
% EndTime: 2019-03-09 18:21:10
% DurationCPUTime: 15.85s
% Computational Cost: add. (13785->670), mult. (32420->872), div. (0->0), fcn. (22764->8), ass. (0->326)
t315 = sin(qJ(6));
t316 = sin(qJ(5));
t319 = cos(qJ(6));
t320 = cos(qJ(5));
t342 = t315 * t320 + t319 * t316;
t468 = qJD(5) + qJD(6);
t225 = t468 * t342;
t317 = sin(qJ(3));
t318 = sin(qJ(2));
t431 = cos(qJ(3));
t432 = cos(qJ(2));
t275 = t317 * t432 + t318 * t431;
t258 = t275 * qJD(1);
t332 = t342 * t258;
t521 = t225 + t332;
t384 = qJD(3) * t317;
t375 = pkin(2) * t384;
t354 = t431 * t432;
t340 = qJD(1) * t354;
t386 = qJD(1) * t318;
t257 = t317 * t386 - t340;
t214 = pkin(3) * t258 + qJ(4) * t257;
t307 = pkin(2) * t386;
t186 = t214 + t307;
t253 = t258 * pkin(9);
t153 = t186 + t253;
t286 = (-pkin(8) - pkin(7)) * t318;
t277 = qJD(1) * t286;
t377 = t432 * pkin(7);
t287 = pkin(8) * t432 + t377;
t278 = t287 * qJD(1);
t370 = t431 * t278;
t221 = t277 * t317 + t370;
t429 = pkin(4) * t257;
t172 = t221 - t429;
t84 = t320 * t153 + t316 * t172;
t520 = t316 * t375 - t84;
t381 = qJD(6) * t315;
t383 = qJD(5) * t316;
t389 = t319 * t320;
t224 = -t315 * t383 - t316 * t381 + t389 * t468;
t341 = t315 * t316 - t389;
t479 = t341 * t258;
t507 = t224 - t479;
t519 = t257 * pkin(5) + pkin(10) * t383;
t314 = qJD(2) + qJD(3);
t232 = t257 * t320 - t314 * t316;
t455 = pkin(3) + pkin(9);
t264 = qJD(2) * pkin(2) + t277;
t246 = t431 * t264;
t391 = t317 * t278;
t219 = -t246 + t391;
t427 = t258 * pkin(4);
t338 = t219 + t427;
t500 = qJD(4) + t338;
t138 = -t314 * t455 + t500;
t305 = -pkin(2) * t432 - pkin(1);
t285 = qJD(1) * t305;
t325 = -t258 * qJ(4) + t285;
t141 = t257 * t455 + t325;
t70 = t138 * t316 + t141 * t320;
t54 = pkin(10) * t232 + t70;
t406 = t315 * t54;
t254 = qJD(5) + t258;
t233 = t257 * t316 + t314 * t320;
t69 = t320 * t138 - t141 * t316;
t53 = -pkin(10) * t233 + t69;
t47 = pkin(5) * t254 + t53;
t19 = t319 * t47 - t406;
t404 = t319 * t54;
t20 = t315 * t47 + t404;
t227 = t314 * t275;
t212 = t227 * qJD(1);
t129 = -qJD(5) * t233 + t212 * t320;
t382 = qJD(5) * t320;
t390 = t317 * t318;
t349 = t314 * t390;
t211 = qJD(1) * t349 - t314 * t340;
t380 = qJD(1) * qJD(2);
t363 = t318 * t380;
t297 = pkin(2) * t363;
t334 = qJ(4) * t211 - qJD(4) * t258 + t297;
t66 = t212 * t455 + t334;
t220 = t317 * t264 + t370;
t326 = qJD(2) * t287;
t263 = t431 * t326;
t279 = qJD(2) * t286;
t267 = qJD(1) * t279;
t123 = qJD(1) * t263 + qJD(3) * t220 + t317 * t267;
t88 = -t211 * pkin(4) + t123;
t15 = t138 * t382 - t141 * t383 + t316 * t88 + t320 * t66;
t10 = pkin(10) * t129 + t15;
t128 = qJD(5) * t232 + t212 * t316;
t402 = qJD(5) * t70;
t16 = -t316 * t66 + t320 * t88 - t402;
t7 = -pkin(5) * t211 - pkin(10) * t128 + t16;
t3 = qJD(6) * t19 + t10 * t319 + t315 * t7;
t481 = qJD(6) * t20;
t4 = -t10 * t315 + t319 * t7 - t481;
t518 = t19 * t521 - t20 * t507 - t342 * t3 + t341 * t4;
t517 = mrSges(5,2) - mrSges(4,1);
t167 = t320 * t172;
t428 = pkin(10) * t258;
t376 = t431 * pkin(2);
t304 = -t376 - pkin(3);
t298 = -pkin(9) + t304;
t473 = -t298 * t383 + t320 * t375;
t516 = t473 - t167 - (-t153 - t428) * t316 + t519;
t395 = t258 * t320;
t242 = pkin(10) * t395;
t423 = -pkin(10) + t298;
t266 = t423 * t320;
t515 = qJD(5) * t266 - t242 + t520;
t170 = t220 - t429;
t159 = t320 * t170;
t164 = t214 + t253;
t514 = t383 * t455 - t159 - (-t164 - t428) * t316 + t519;
t422 = -pkin(10) - t455;
t281 = t422 * t320;
t87 = t320 * t164 + t316 * t170;
t513 = qJD(5) * t281 - t242 - t87;
t312 = t314 * qJ(4);
t155 = t170 + t312;
t108 = -pkin(5) * t232 + t155;
t145 = t232 * t315 + t233 * t319;
t359 = t319 * t232 - t233 * t315;
t38 = qJD(6) * t359 + t128 * t319 + t129 * t315;
t39 = -qJD(6) * t145 - t128 * t315 + t129 * t319;
t378 = Ifges(7,5) * t38 + Ifges(7,6) * t39 - Ifges(7,3) * t211;
t416 = Ifges(7,4) * t145;
t244 = qJD(6) + t254;
t442 = -t244 / 0.2e1;
t450 = -t145 / 0.2e1;
t466 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t512 = t378 + t466 + (Ifges(7,5) * t359 - Ifges(7,6) * t145) * t442 + (t145 * t20 + t19 * t359) * mrSges(7,3) - t108 * (mrSges(7,1) * t145 + mrSges(7,2) * t359) + (Ifges(7,1) * t359 - t416) * t450;
t511 = -Ifges(4,5) + Ifges(5,4);
t510 = -Ifges(4,6) + Ifges(5,5);
t509 = t316 * t69;
t231 = Ifges(6,4) * t232;
t119 = Ifges(6,1) * t233 + Ifges(6,5) * t254 + t231;
t407 = t258 * Ifges(5,6);
t508 = t314 * Ifges(5,5) + t257 * Ifges(5,3) + t316 * t119 - t407;
t371 = t431 * t277;
t222 = t371 - t391;
t364 = qJD(3) * t431;
t355 = pkin(2) * t364;
t290 = t355 + qJD(4);
t505 = -t290 + t222;
t372 = -pkin(5) * t320 - pkin(4);
t504 = pkin(5) * t382 - t258 * t372 + t391;
t259 = t317 * t326;
t122 = -qJD(1) * t259 + t264 * t364 + t431 * t267 - t278 * t384;
t115 = -t314 * qJD(4) - t122;
t469 = -qJD(4) - t219;
t185 = -pkin(3) * t314 - t469;
t250 = Ifges(5,6) * t257;
t192 = t314 * Ifges(5,4) - t258 * Ifges(5,2) + t250;
t408 = t258 * Ifges(4,4);
t193 = -t257 * Ifges(4,2) + t314 * Ifges(4,6) + t408;
t251 = Ifges(4,4) * t257;
t345 = Ifges(6,5) * t316 + Ifges(6,6) * t320;
t418 = Ifges(6,4) * t316;
t346 = Ifges(6,2) * t320 + t418;
t417 = Ifges(6,4) * t320;
t347 = Ifges(6,1) * t316 + t417;
t361 = -t383 / 0.2e1;
t80 = -pkin(4) * t212 - t115;
t42 = -pkin(5) * t129 + t80;
t421 = mrSges(5,1) * t257;
t434 = -t314 / 0.2e1;
t436 = t258 / 0.2e1;
t437 = -t258 / 0.2e1;
t438 = t257 / 0.2e1;
t439 = -t257 / 0.2e1;
t440 = -t254 / 0.2e1;
t441 = t244 / 0.2e1;
t444 = -t233 / 0.2e1;
t445 = -t232 / 0.2e1;
t448 = -t211 / 0.2e1;
t449 = t145 / 0.2e1;
t45 = t128 * Ifges(6,4) + t129 * Ifges(6,2) - t211 * Ifges(6,6);
t451 = t359 / 0.2e1;
t452 = -t359 / 0.2e1;
t453 = t129 / 0.2e1;
t454 = t128 / 0.2e1;
t139 = Ifges(7,4) * t359;
t62 = Ifges(7,1) * t145 + Ifges(7,5) * t244 + t139;
t456 = t62 / 0.2e1;
t457 = -t62 / 0.2e1;
t61 = Ifges(7,2) * t359 + Ifges(7,6) * t244 + t416;
t458 = -t61 / 0.2e1;
t459 = Ifges(6,1) * t454 + Ifges(6,4) * t453 + Ifges(6,5) * t448;
t460 = t39 / 0.2e1;
t461 = t38 / 0.2e1;
t462 = Ifges(7,1) * t461 + Ifges(7,4) * t460 + Ifges(7,5) * t448;
t463 = Ifges(7,4) * t461 + Ifges(7,2) * t460 + Ifges(7,6) * t448;
t178 = t257 * pkin(3) + t325;
t495 = t69 * mrSges(6,1) + t19 * mrSges(7,1) + t285 * mrSges(4,2) - t70 * mrSges(6,2) - t20 * mrSges(7,2) + t219 * mrSges(4,3) - t178 * mrSges(5,3);
t497 = t285 * mrSges(4,1) - t178 * mrSges(5,2);
t485 = t254 * Ifges(6,3);
t486 = t232 * Ifges(6,6);
t499 = t258 * Ifges(4,1) + t314 * Ifges(4,5) + t233 * Ifges(6,5) + t145 * Ifges(7,5) + Ifges(7,6) * t359 + t244 * Ifges(7,3) - t251 + t485 + t486;
t348 = mrSges(6,1) * t320 - mrSges(6,2) * t316;
t419 = Ifges(6,4) * t233;
t118 = Ifges(6,2) * t232 + Ifges(6,6) * t254 + t419;
t388 = t320 * t118;
t502 = -t388 / 0.2e1 + t155 * t348;
t503 = (mrSges(6,3) * t509 + Ifges(5,3) * t439 + t345 * t440 + t346 * t445 + t347 * t444 + t434 * t510 - t497 + t502) * t258 + t510 * t212 + (-Ifges(6,5) * t444 - Ifges(7,5) * t450 + Ifges(5,2) * t436 - Ifges(6,6) * t445 - Ifges(7,6) * t452 - Ifges(6,3) * t440 - Ifges(7,3) * t442 + t434 * t511 + t495) * t257 + t511 * t211 + t517 * t123 + t507 * t458 + (-Ifges(4,1) * t257 - t408 + t508) * t437 + t502 * qJD(5) + (-Ifges(4,2) * t258 - t251 + t499) * t438 + t518 * mrSges(7,3) - (t232 * t346 + t233 * t347 + t254 * t345) * qJD(5) / 0.2e1 + (t250 + t192) * t439 + (Ifges(7,5) * t332 - Ifges(7,6) * t479) * t442 + (Ifges(7,4) * t332 - Ifges(7,2) * t479) * t452 + (Ifges(7,1) * t332 - Ifges(7,4) * t479) * t450 - t342 * t463 + (Ifges(6,5) * t320 - Ifges(7,5) * t341 - Ifges(6,6) * t316 - Ifges(7,6) * t342) * t448 + t42 * (mrSges(7,1) * t342 - mrSges(7,2) * t341) + (-Ifges(7,4) * t341 - Ifges(7,2) * t342) * t460 + (-Ifges(7,1) * t341 - Ifges(7,4) * t342) * t461 - t341 * t462 + (-Ifges(7,5) * t225 - Ifges(7,6) * t224) * t441 + (-Ifges(7,1) * t225 - Ifges(7,4) * t224) * t449 + (-Ifges(7,4) * t225 - Ifges(7,2) * t224) * t451 + (t407 + t193) * t436 + (t383 * t69 - t395 * t70) * mrSges(6,3) + t332 * t457 + t119 * t361 + t80 * (mrSges(6,1) * t316 + mrSges(6,2) * t320) - t316 * t45 / 0.2e1 + (mrSges(7,1) * t507 - mrSges(7,2) * t521) * t108 - t115 * mrSges(5,3) - t122 * mrSges(4,2) + t185 * t421 + (-Ifges(6,2) * t316 + t417) * t453 + (Ifges(6,1) * t320 - t418) * t454 - t225 * t456 + t320 * t459;
t501 = t375 - t221;
t498 = -Ifges(7,2) * t145 + t139;
t265 = t423 * t316;
t210 = t265 * t319 + t266 * t315;
t493 = -qJD(6) * t210 - t315 * t515 + t319 * t516;
t209 = -t265 * t315 + t266 * t319;
t492 = qJD(6) * t209 + t315 * t516 + t319 * t515;
t280 = t422 * t316;
t230 = t280 * t319 + t281 * t315;
t483 = -qJD(6) * t230 - t315 * t513 + t319 * t514;
t229 = -t280 * t315 + t281 * t319;
t482 = qJD(6) * t229 + t315 * t514 + t319 * t513;
t273 = -t354 + t390;
t199 = t341 * t273;
t478 = -t371 + t290 + t504;
t331 = -t275 * qJ(4) + t305;
t171 = t273 * t455 + t331;
t236 = -t431 * t286 + t287 * t317;
t196 = pkin(4) * t275 + t236;
t183 = t316 * t196;
t95 = t320 * t171 + t183;
t238 = -mrSges(4,2) * t314 - mrSges(4,3) * t257;
t240 = -mrSges(5,3) * t314 + t421;
t477 = t238 - t240;
t409 = t258 * mrSges(4,3);
t476 = t258 * mrSges(5,1) + t314 * t517 + t409;
t475 = t427 - t505;
t474 = qJD(4) - t246 + t504;
t366 = qJD(1) * t432;
t472 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t386) * t377 + t318 * pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t366);
t89 = -mrSges(6,1) * t211 - mrSges(6,3) * t128;
t90 = mrSges(6,2) * t211 + mrSges(6,3) * t129;
t471 = t316 * t90 + t320 * t89;
t470 = t15 * t316 + t16 * t320;
t237 = t317 * t286 + t287 * t431;
t467 = t123 * t275 - t211 * t236 - t212 * t237;
t465 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t443 = t233 / 0.2e1;
t433 = t314 / 0.2e1;
t430 = pkin(2) * t317;
t420 = Ifges(3,4) * t318;
t411 = t211 * mrSges(5,1);
t401 = t123 * t236;
t198 = -t312 - t220;
t399 = t198 * t258;
t394 = t273 * t316;
t393 = t273 * t320;
t154 = -mrSges(6,1) * t232 + mrSges(6,2) * t233;
t387 = t154 - t240;
t385 = qJD(2) * t318;
t310 = pkin(2) * t385;
t374 = Ifges(3,4) * t432;
t373 = Ifges(6,5) * t128 + Ifges(6,6) * t129 - Ifges(6,3) * t211;
t365 = qJD(2) * t432;
t299 = qJ(4) + t430;
t360 = -pkin(10) * t273 - t171;
t353 = qJD(1) * t365;
t184 = t320 * t196;
t73 = pkin(5) * t275 + t316 * t360 + t184;
t85 = pkin(10) * t393 + t95;
t29 = -t315 * t85 + t319 * t73;
t30 = t315 * t73 + t319 * t85;
t344 = -t320 * t70 + t509;
t174 = -mrSges(6,2) * t254 + mrSges(6,3) * t232;
t175 = mrSges(6,1) * t254 - mrSges(6,3) * t233;
t343 = t174 * t320 - t175 * t316;
t337 = Ifges(3,5) * t432 - Ifges(3,6) * t318;
t336 = t227 * t316 + t273 * t382;
t335 = -t227 * t320 + t273 * t383;
t226 = -t314 * t354 + t349;
t333 = qJ(4) * t226 - qJD(4) * t275 + t310;
t147 = qJD(3) * t237 + t317 * t279 + t263;
t103 = -t226 * pkin(4) + t147;
t79 = t227 * t455 + t333;
t23 = t316 * t103 - t171 * t383 + t196 * t382 + t320 * t79;
t146 = -t431 * t279 - t286 * t364 + t287 * t384 + t259;
t330 = pkin(1) * (mrSges(3,1) * t318 + mrSges(3,2) * t432);
t329 = t318 * (Ifges(3,1) * t432 - t420);
t328 = (Ifges(3,2) * t432 + t420) * qJD(1);
t102 = -pkin(4) * t227 - t146;
t324 = -qJD(5) * t344 + t470;
t323 = m(6) * t324;
t313 = t316 * pkin(5);
t306 = Ifges(3,4) * t366;
t300 = qJ(4) + t313;
t282 = t299 + t313;
t256 = Ifges(3,1) * t386 + Ifges(3,5) * qJD(2) + t306;
t255 = Ifges(3,6) * qJD(2) + t328;
t218 = t273 * pkin(3) + t331;
t217 = -mrSges(5,2) * t257 - mrSges(5,3) * t258;
t216 = mrSges(4,1) * t257 + mrSges(4,2) * t258;
t200 = t342 * t273;
t197 = -t273 * pkin(4) + t237;
t156 = t273 * t372 + t237;
t110 = mrSges(7,1) * t244 - mrSges(7,3) * t145;
t109 = -mrSges(7,2) * t244 + mrSges(7,3) * t359;
t106 = pkin(3) * t227 + t333;
t99 = t320 * t103;
t94 = -t171 * t316 + t184;
t91 = pkin(3) * t212 + t334;
t86 = -t164 * t316 + t159;
t83 = -t153 * t316 + t167;
t77 = -mrSges(7,1) * t359 + mrSges(7,2) * t145;
t65 = -mrSges(6,1) * t129 + mrSges(6,2) * t128;
t63 = pkin(5) * t335 + t102;
t56 = -t225 * t273 - t227 * t341;
t55 = -t199 * t468 + t342 * t227;
t32 = mrSges(7,2) * t211 + mrSges(7,3) * t39;
t31 = -mrSges(7,1) * t211 - mrSges(7,3) * t38;
t24 = -qJD(5) * t95 - t316 * t79 + t99;
t22 = t319 * t53 - t406;
t21 = -t315 * t53 - t404;
t18 = -pkin(10) * t335 + t23;
t14 = -pkin(5) * t226 + t99 + (-pkin(10) * t227 - t79) * t316 + (t320 * t360 - t183) * qJD(5);
t13 = -mrSges(7,1) * t39 + mrSges(7,2) * t38;
t6 = -qJD(6) * t30 + t14 * t319 - t18 * t315;
t5 = qJD(6) * t29 + t14 * t315 + t18 * t319;
t1 = [(t198 * t227 + t467) * mrSges(5,1) + (Ifges(5,3) * t438 - t193 / 0.2e1 - Ifges(4,4) * t436 - Ifges(4,2) * t439 + t510 * t433 + t497) * t227 + (t337 * qJD(2) / 0.2e1 - t472) * qJD(2) + (m(4) * t219 + m(5) * t185 + t476) * t147 + (-m(4) * t220 + m(5) * t198 - t477) * t146 + (-t220 * t227 + t467) * mrSges(4,3) + (-t495 + t511 * t433 - t485 / 0.2e1 + t192 / 0.2e1 - Ifges(4,1) * t436 + Ifges(5,2) * t437 - Ifges(4,4) * t439 - Ifges(7,3) * t441 - Ifges(6,5) * t443 - Ifges(7,5) * t449 - Ifges(7,6) * t451 - t486 / 0.2e1 - t185 * mrSges(5,1) - t499 / 0.2e1 + t438 * Ifges(5,6)) * t226 + (t211 * t273 + t227 * t437) * Ifges(5,6) + (t388 + t508) * t227 / 0.2e1 + (Ifges(6,1) * t336 - Ifges(6,4) * t335) * t443 + t216 * t310 + t232 * (Ifges(6,4) * t336 - Ifges(6,2) * t335) / 0.2e1 + t155 * (mrSges(6,1) * t335 + mrSges(6,2) * t336) + (mrSges(4,1) * t297 + t115 * mrSges(5,1) - t91 * mrSges(5,2) - t122 * mrSges(4,3) + t118 * t361 + t345 * t448 + t346 * t453 + t347 * t454 - t80 * t348 + (Ifges(5,3) + Ifges(4,2)) * t212 + (t211 / 0.2e1 - t448) * Ifges(4,4)) * t273 + (-0.2e1 * t330 + t329) * t380 + (qJD(1) * (Ifges(3,1) * t318 + t374) + t256) * t365 / 0.2e1 - (t328 + t255) * t385 / 0.2e1 + (qJD(5) * t119 + t45) * t393 / 0.2e1 + t254 * (Ifges(6,5) * t336 - Ifges(6,6) * t335) / 0.2e1 + (t15 * t393 - t16 * t394 - t335 * t70 - t336 * t69) * mrSges(6,3) + (-Ifges(3,2) * t318 + t374) * t353 + m(4) * (t122 * t237 + 0.2e1 * t285 * t310 + t401) + m(5) * (t106 * t178 - t115 * t237 + t218 * t91 + t401) + (-t19 * t55 - t199 * t3 + t20 * t56 - t200 * t4) * mrSges(7,3) + (Ifges(7,1) * t200 - Ifges(7,4) * t199) * t461 + (Ifges(7,4) * t200 - Ifges(7,2) * t199) * t460 + (Ifges(7,5) * t200 - Ifges(7,6) * t199) * t448 + t42 * (mrSges(7,1) * t199 + mrSges(7,2) * t200) + (Ifges(7,1) * t55 + Ifges(7,4) * t56) * t449 + m(7) * (t108 * t63 + t156 * t42 + t19 * t6 + t20 * t5 + t29 * t4 + t3 * t30) + m(6) * (t102 * t155 + t15 * t95 + t16 * t94 + t197 * t80 + t23 * t70 + t24 * t69) + (Ifges(7,4) * t55 + Ifges(7,2) * t56) * t451 + t305 * (mrSges(4,1) * t212 - mrSges(4,2) * t211) + (Ifges(7,5) * t55 + Ifges(7,6) * t56) * t441 + (mrSges(4,2) * t297 - t91 * mrSges(5,3) + Ifges(6,6) * t453 + Ifges(6,5) * t454 + Ifges(7,6) * t460 + Ifges(7,5) * t461 + (Ifges(7,3) + Ifges(6,3) + Ifges(4,1)) * t448 + t465 + t466 + t373 / 0.2e1 + t378 / 0.2e1 + (-Ifges(5,2) - Ifges(4,1) / 0.2e1) * t211 + (-Ifges(5,6) - Ifges(4,4)) * t212) * t275 + t29 * t31 + t30 * t32 + t56 * t61 / 0.2e1 + t63 * t77 + t94 * t89 + t95 * t90 + t108 * (-mrSges(7,1) * t56 + mrSges(7,2) * t55) + t5 * t109 + t6 * t110 + t102 * t154 + t156 * t13 + t23 * t174 + t24 * t175 + t197 * t65 + t106 * t217 + t218 * (-mrSges(5,2) * t212 + mrSges(5,3) * t211) + t55 * t456 + t394 * t459 + t200 * t462 - t199 * t463; t478 * t77 + (-t382 * t70 - t470) * mrSges(6,3) + t471 * t298 + (-t83 + t473) * t175 + t475 * t154 + (t299 * t80 + (t316 * t70 + t320 * t69) * t375 + t324 * t298 - t69 * t83 - t70 * t84 + t475 * t155) * m(6) - t477 * t222 + t501 * t476 + (-t115 * t299 + t123 * t304 - t178 * t186 + t185 * t501 + t198 * t505) * m(5) + (t298 * t382 + t520) * t174 + t238 * t355 + Ifges(3,5) * t353 - t304 * t411 + t255 * t386 / 0.2e1 - (-Ifges(3,2) * t386 + t256 + t306) * t366 / 0.2e1 - t337 * t380 / 0.2e1 + ((-t431 * t123 + t122 * t317 + (t219 * t317 + t220 * t431) * qJD(3)) * pkin(2) - t219 * t221 - t220 * t222 - t285 * t307) * m(4) - Ifges(3,6) * t363 + (t472 + (-t329 / 0.2e1 + t330) * qJD(1)) * qJD(1) - t216 * t307 + t503 - t290 * t240 + t299 * t65 + (-t212 * t299 - t399) * mrSges(5,1) + t492 * t109 + t493 * t110 + (t108 * t478 + t19 * t493 + t20 * t492 + t209 * t4 + t210 * t3 + t282 * t42) * m(7) + t220 * t409 + (t211 * t376 - t212 * t430) * mrSges(4,3) + (-mrSges(3,1) * t353 + mrSges(3,2) * t363) * pkin(7) + t209 * t31 + t210 * t32 - t186 * t217 + t282 * t13; t474 * t77 + (t409 - t476) * t220 + t477 * t219 + (-pkin(3) * t123 - qJ(4) * t115 - t178 * t214 - t185 * t220 + t198 * t469) * m(5) + (t80 * qJ(4) + t155 * t500 - t69 * t86 - t70 * t87) * m(6) + t338 * t154 + (-(qJD(5) * t174 + t89) * t455 + (-t16 - t402) * mrSges(6,3)) * t320 - t455 * t323 + (-t15 * mrSges(6,3) - (-qJD(5) * t175 + t90) * t455) * t316 + (pkin(3) * t211 - qJ(4) * t212 - t399) * mrSges(5,1) + t387 * qJD(4) + t503 + t300 * t13 + t482 * t109 + t483 * t110 + (t108 * t474 + t19 * t483 + t20 * t482 + t229 * t4 + t230 * t3 + t300 * t42) * m(7) + qJ(4) * t65 - t87 * t174 - t86 * t175 - t214 * t217 + t229 * t31 + t230 * t32; -t411 + t342 * t32 - t341 * t31 - t521 * t110 + t507 * t109 + t343 * qJD(5) + (-t77 - t387) * t314 + (t217 + t343) * t258 + t323 - m(6) * (t155 * t314 + t258 * t344) + (-t108 * t314 - t518) * m(7) + (t178 * t258 + t198 * t314 + t123) * m(5) + t471; t465 - t145 * t458 + t373 + (t232 * t69 + t233 * t70) * mrSges(6,3) + t359 * t457 - m(7) * (t19 * t21 + t20 * t22) + (-Ifges(6,2) * t233 + t119 + t231) * t445 + (-t233 * t77 + t319 * t31 + t315 * t32 + (t109 * t319 - t110 * t315) * qJD(6) + (-t108 * t233 - t19 * t381 + t3 * t315 + (t4 + t481) * t319) * m(7)) * pkin(5) - t22 * t109 - t21 * t110 - t69 * t174 + t70 * t175 - t155 * (mrSges(6,1) * t233 + mrSges(6,2) * t232) + (Ifges(6,5) * t232 - Ifges(6,6) * t233) * t440 + t118 * t443 + (Ifges(6,1) * t232 - t419) * t444 + t498 * t452 + t512; t61 * t449 - t19 * t109 + t20 * t110 + (t498 + t62) * t452 + t512;];
tauc  = t1(:);
