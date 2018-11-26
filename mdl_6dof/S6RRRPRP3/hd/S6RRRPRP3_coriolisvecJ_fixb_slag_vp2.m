% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:42
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:42:13
% EndTime: 2018-11-23 17:42:24
% DurationCPUTime: 10.95s
% Computational Cost: add. (13495->625), mult. (33475->822), div. (0->0), fcn. (24277->8), ass. (0->288)
t328 = sin(qJ(3));
t329 = sin(qJ(2));
t331 = cos(qJ(3));
t332 = cos(qJ(2));
t300 = t328 * t332 + t331 * t329;
t287 = t300 * qJD(1);
t324 = qJD(2) + qJD(3);
t325 = sin(pkin(10));
t326 = cos(pkin(10));
t263 = t287 * t326 + t324 * t325;
t327 = sin(qJ(5));
t330 = cos(qJ(5));
t339 = t287 * t325 - t324 * t326;
t182 = t327 * t263 + t330 * t339;
t299 = t328 * t329 - t331 * t332;
t257 = t324 * t299;
t243 = t257 * qJD(1);
t297 = t325 * t327 - t330 * t326;
t86 = -qJD(5) * t182 + t243 * t297;
t440 = t86 / 0.2e1;
t298 = t325 * t330 + t326 * t327;
t471 = t330 * t263 - t327 * t339;
t87 = qJD(5) * t471 - t243 * t298;
t439 = -t87 / 0.2e1;
t258 = t324 * t300;
t244 = t258 * qJD(1);
t421 = t244 / 0.2e1;
t460 = Ifges(6,1) + Ifges(7,1);
t480 = Ifges(6,4) - Ifges(7,5);
t458 = Ifges(7,4) + Ifges(6,5);
t363 = qJD(1) * t332;
t364 = qJD(1) * t329;
t286 = t328 * t364 - t331 * t363;
t281 = qJD(5) + t286;
t419 = t281 / 0.2e1;
t424 = t471 / 0.2e1;
t426 = t182 / 0.2e1;
t427 = -t182 / 0.2e1;
t478 = Ifges(6,6) - Ifges(7,6);
t435 = -pkin(8) - pkin(7);
t313 = t435 * t332;
t303 = qJD(1) * t313;
t288 = t328 * t303;
t312 = t435 * t329;
t302 = qJD(1) * t312;
t294 = qJD(2) * pkin(2) + t302;
t251 = t331 * t294 + t288;
t217 = -t324 * pkin(3) + qJD(4) - t251;
t168 = pkin(4) * t339 + t217;
t321 = -pkin(2) * t332 - pkin(1);
t311 = qJD(1) * t321;
t213 = t286 * pkin(3) - t287 * qJ(4) + t311;
t368 = t331 * t303;
t252 = t328 * t294 - t368;
t226 = qJ(4) * t324 + t252;
t139 = t325 * t213 + t326 * t226;
t108 = -pkin(9) * t339 + t139;
t138 = t326 * t213 - t226 * t325;
t98 = pkin(4) * t286 - pkin(9) * t263 + t138;
t31 = t108 * t330 + t327 * t98;
t26 = qJ(6) * t281 + t31;
t47 = t182 * pkin(5) - qJ(6) * t471 + t168;
t179 = Ifges(7,5) * t471;
t89 = Ifges(7,6) * t281 + Ifges(7,3) * t182 + t179;
t395 = Ifges(6,4) * t471;
t92 = -Ifges(6,2) * t182 + Ifges(6,6) * t281 + t395;
t483 = -mrSges(6,3) * t31 + mrSges(6,1) * t168 + mrSges(7,1) * t47 - mrSges(7,2) * t26 - t92 / 0.2e1 + t89 / 0.2e1;
t488 = -Ifges(6,2) * t427 + Ifges(7,3) * t426 - t419 * t478 - t424 * t480 + t483;
t379 = t243 * t326;
t362 = qJD(2) * t329;
t354 = pkin(2) * t362;
t472 = qJD(1) * t354;
t124 = pkin(3) * t244 + qJ(4) * t243 - qJD(4) * t287 + t472;
t353 = qJD(2) * t435;
t348 = qJD(1) * t353;
t295 = t329 * t348;
t336 = t332 * t348;
t360 = qJD(3) * t331;
t361 = qJD(3) * t328;
t166 = t294 * t360 + t331 * t295 + t303 * t361 + t328 * t336;
t156 = qJD(4) * t324 + t166;
t60 = t326 * t124 - t156 * t325;
t29 = pkin(4) * t244 + pkin(9) * t379 + t60;
t358 = qJD(5) * t330;
t359 = qJD(5) * t327;
t380 = t243 * t325;
t61 = t325 * t124 + t326 * t156;
t38 = pkin(9) * t380 + t61;
t5 = -t108 * t359 + t327 * t29 + t330 * t38 + t98 * t358;
t1 = qJ(6) * t244 + qJD(6) * t281 + t5;
t167 = qJD(3) * t252 + t295 * t328 - t331 * t336;
t116 = -pkin(4) * t380 + t167;
t12 = pkin(5) * t87 - qJ(6) * t86 - qJD(6) * t471 + t116;
t438 = t87 / 0.2e1;
t465 = -Ifges(6,6) / 0.2e1;
t487 = mrSges(6,1) * t116 + mrSges(7,1) * t12 - mrSges(7,2) * t1 - Ifges(6,4) * t86 / 0.2e1 + t244 * t465 + 0.2e1 * Ifges(7,3) * t438 - mrSges(6,3) * t5 + (-t439 + t438) * Ifges(6,2) + (-t480 + Ifges(7,5)) * t440 + (-t478 + Ifges(7,6)) * t421;
t479 = Ifges(7,2) + Ifges(6,3);
t317 = pkin(2) * t328 + qJ(4);
t292 = (-pkin(9) - t317) * t325;
t323 = t326 * pkin(9);
t370 = t317 * t326;
t293 = t323 + t370;
t237 = t292 * t327 + t293 * t330;
t314 = pkin(2) * t360 + qJD(4);
t172 = qJD(5) * t237 + t298 * t314;
t245 = pkin(3) * t287 + qJ(4) * t286;
t218 = pkin(2) * t364 + t245;
t254 = t302 * t331 + t288;
t154 = t326 * t218 - t254 * t325;
t375 = t286 * t326;
t347 = t287 * pkin(4) + pkin(9) * t375;
t112 = t154 + t347;
t155 = t325 * t218 + t326 * t254;
t376 = t286 * t325;
t357 = pkin(9) * t376;
t125 = t357 + t155;
t39 = t112 * t330 - t125 * t327;
t486 = -t39 - t172;
t338 = t330 * t292 - t293 * t327;
t171 = qJD(5) * t338 - t297 * t314;
t40 = t327 * t112 + t330 * t125;
t485 = -t40 + t171;
t211 = t298 * t286;
t212 = t297 * t286;
t270 = pkin(4) * t376;
t282 = t297 * qJD(5);
t283 = t298 * qJD(5);
t484 = -qJD(6) * t298 + t270 + (t212 + t282) * qJ(6) + (t211 + t283) * pkin(5);
t6 = -qJD(5) * t31 + t29 * t330 - t327 * t38;
t3 = -pkin(5) * t244 - t6;
t482 = mrSges(6,2) * t116 + mrSges(7,2) * t3 - mrSges(6,3) * t6 - mrSges(7,3) * t12 + Ifges(7,5) * t438 + 0.2e1 * t421 * t458 + 0.2e1 * t440 * t460 + (Ifges(6,4) + t480) * t439;
t30 = -t108 * t327 + t330 * t98;
t448 = qJD(6) - t30;
t25 = -pkin(5) * t281 + t448;
t180 = Ifges(6,4) * t182;
t393 = Ifges(7,5) * t182;
t456 = t281 * t458 + t460 * t471 - t180 + t393;
t481 = t419 * t458 + t424 * t460 + t456 / 0.2e1 + mrSges(6,2) * t168 + mrSges(7,2) * t25 - mrSges(6,3) * t30 - mrSges(7,3) * t47 + Ifges(6,4) * t427 + Ifges(7,5) * t426;
t277 = t287 * qJ(6);
t475 = -t277 + t485;
t407 = pkin(5) * t287;
t474 = t407 - t486;
t159 = t326 * t245 - t251 * t325;
t115 = t159 + t347;
t160 = t325 * t245 + t326 * t251;
t131 = t357 + t160;
t307 = (-pkin(9) - qJ(4)) * t325;
t308 = qJ(4) * t326 + t323;
t261 = t307 * t327 + t308 * t330;
t452 = -qJD(4) * t298 - qJD(5) * t261 - t115 * t330 + t131 * t327;
t337 = t330 * t307 - t308 * t327;
t451 = -qJD(4) * t297 + qJD(5) * t337 - t327 * t115 - t330 * t131;
t253 = t302 * t328 - t368;
t473 = pkin(2) * t361 - t253 + t484;
t455 = mrSges(4,2) * t287;
t454 = -t277 + t451;
t453 = t407 - t452;
t450 = -t252 + t484;
t346 = mrSges(5,1) * t325 + mrSges(5,2) * t326;
t449 = t217 * t346;
t250 = t299 * pkin(3) - t300 * qJ(4) + t321;
t266 = t312 * t328 - t313 * t331;
t173 = t326 * t250 - t266 * t325;
t373 = t300 * t326;
t130 = pkin(4) * t299 - pkin(9) * t373 + t173;
t174 = t325 * t250 + t326 * t266;
t374 = t300 * t325;
t145 = -pkin(9) * t374 + t174;
t447 = t327 * t130 + t330 * t145;
t446 = t331 * t312 + t313 * t328;
t445 = t479 * t244 + t458 * t86 - t478 * t87;
t146 = pkin(3) * t258 + qJ(4) * t257 - qJD(4) * t300 + t354;
t304 = t329 * t353;
t305 = t332 * t353;
t186 = qJD(3) * t446 + t304 * t331 + t305 * t328;
t77 = t326 * t146 - t186 * t325;
t44 = pkin(4) * t258 + t257 * t323 + t77;
t377 = t257 * t325;
t78 = t325 * t146 + t326 * t186;
t62 = pkin(9) * t377 + t78;
t10 = -qJD(5) * t447 - t327 * t62 + t330 * t44;
t444 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t1 * mrSges(7,3);
t443 = m(4) / 0.2e1;
t434 = pkin(1) * mrSges(3,1);
t433 = pkin(1) * mrSges(3,2);
t425 = -t471 / 0.2e1;
t420 = -t281 / 0.2e1;
t417 = -t286 / 0.2e1;
t416 = t286 / 0.2e1;
t415 = -t287 / 0.2e1;
t412 = -t325 / 0.2e1;
t411 = t326 / 0.2e1;
t409 = m(4) * t311;
t408 = pkin(2) * t331;
t406 = t60 * mrSges(5,3);
t57 = mrSges(6,1) * t244 - mrSges(6,3) * t86;
t58 = -t244 * mrSges(7,1) + t86 * mrSges(7,2);
t405 = t58 - t57;
t56 = -mrSges(7,2) * t87 + mrSges(7,3) * t244;
t59 = -mrSges(6,2) * t244 - mrSges(6,3) * t87;
t404 = t59 + t56;
t403 = mrSges(4,3) * t286;
t402 = mrSges(4,3) * t287;
t401 = mrSges(6,3) * t182;
t400 = mrSges(6,3) * t471;
t399 = Ifges(3,4) * t329;
t398 = Ifges(4,4) * t287;
t397 = Ifges(5,4) * t325;
t396 = Ifges(5,4) * t326;
t394 = Ifges(5,5) * t263;
t392 = Ifges(5,2) * t325;
t391 = t287 * Ifges(4,1);
t390 = t324 * Ifges(4,5);
t389 = t324 * Ifges(4,6);
t388 = t325 * t60;
t387 = t326 * t61;
t386 = Ifges(3,5) * qJD(2);
t385 = Ifges(3,6) * qJD(2);
t384 = qJD(2) * mrSges(3,1);
t383 = qJD(2) * mrSges(3,2);
t382 = t138 * t325;
t381 = t167 * t446;
t371 = t314 * t326;
t157 = -mrSges(5,2) * t244 + mrSges(5,3) * t380;
t369 = t326 * t157;
t148 = -mrSges(7,2) * t182 + mrSges(7,3) * t281;
t149 = -mrSges(6,2) * t281 - t401;
t367 = -t148 - t149;
t150 = mrSges(6,1) * t281 - t400;
t151 = -mrSges(7,1) * t281 + mrSges(7,2) * t471;
t366 = t150 - t151;
t365 = -mrSges(4,1) * t324 + mrSges(5,1) * t339 + t263 * mrSges(5,2) + t402;
t153 = -mrSges(5,1) * t380 - mrSges(5,2) * t379;
t318 = -pkin(4) * t326 - pkin(3);
t352 = t386 / 0.2e1;
t351 = -t385 / 0.2e1;
t350 = (Ifges(5,4) * t263 - Ifges(5,2) * t339 + Ifges(5,6) * t286) * t412;
t349 = (Ifges(5,1) * t263 - Ifges(5,4) * t339 + Ifges(5,5) * t286) * t411;
t24 = t87 * mrSges(6,1) + t86 * mrSges(6,2);
t23 = t87 * mrSges(7,1) - t86 * mrSges(7,3);
t215 = pkin(4) * t374 - t446;
t345 = Ifges(5,1) * t326 - t397;
t344 = -t392 + t396;
t343 = Ifges(5,5) * t326 - Ifges(5,6) * t325;
t48 = t130 * t330 - t145 * t327;
t9 = t130 * t358 - t145 * t359 + t327 * t44 + t330 * t62;
t334 = Ifges(5,6) * t339;
t242 = pkin(5) * t297 - qJ(6) * t298 + t318;
t187 = qJD(3) * t266 + t304 * t328 - t331 * t305;
t134 = -pkin(4) * t377 + t187;
t109 = t244 * Ifges(5,6) - t243 * t344;
t110 = t244 * Ifges(5,5) - t243 * t345;
t162 = Ifges(5,3) * t286 - t334 + t394;
t223 = -t286 * Ifges(4,2) + t389 + t398;
t279 = Ifges(4,4) * t286;
t224 = -t279 + t390 + t391;
t90 = Ifges(6,5) * t471 - t182 * Ifges(6,6) + t281 * Ifges(6,3);
t91 = Ifges(7,4) * t471 + t281 * Ifges(7,2) + t182 * Ifges(7,6);
t333 = (Ifges(6,2) * t426 - Ifges(7,3) * t427 + t420 * t478 + t425 * t480 + t483) * t211 - t481 * t282 + t482 * t298 - t30 * (mrSges(6,1) * t287 - mrSges(6,3) * t212) - t25 * (-mrSges(7,1) * t287 + mrSges(7,2) * t212) - t139 * (-mrSges(5,2) * t287 + mrSges(5,3) * t376) - t138 * (mrSges(5,1) * t287 + mrSges(5,3) * t375) - t311 * (mrSges(4,1) * t287 - mrSges(4,2) * t286) - t456 * t212 / 0.2e1 + (-t168 * t212 + t287 * t31) * mrSges(6,2) + (Ifges(6,4) * t212 + Ifges(6,6) * t287) * t426 + (-Ifges(4,2) * t287 + t224 - t279) * t416 + (t212 * t47 - t26 * t287) * mrSges(7,3) + (-mrSges(5,1) * t326 + mrSges(5,2) * t325 - mrSges(4,1)) * t167 + (-Ifges(4,1) * t286 + t162 - t398 + t90 + t91) * t415 + t487 * t297 + t488 * t283 + (Ifges(7,5) * t212 + Ifges(7,6) * t287) * t427 + (Ifges(5,2) * t326 + t397) * t380 / 0.2e1 + mrSges(5,3) * t387 + (t212 * t458 + t287 * t479) * t420 + (t212 * t460 + t287 * t458) * t425 + t286 * t349 + t286 * t350 - t251 * t403 + t109 * t411 + (Ifges(5,3) * t287 - t286 * t343) * t417 + (Ifges(5,5) * t325 + Ifges(5,6) * t326) * t421 - t166 * mrSges(4,2) + t339 * (Ifges(5,6) * t287 - t286 * t344) / 0.2e1 - t263 * (Ifges(5,5) * t287 - t286 * t345) / 0.2e1 - Ifges(4,5) * t243 - Ifges(4,6) * t244 + t286 * t449 + t287 * t223 / 0.2e1 - t324 * (-Ifges(4,5) * t286 - Ifges(4,6) * t287) / 0.2e1 + t325 * t110 / 0.2e1 - (Ifges(5,1) * t325 + t396) * t379 / 0.2e1;
t322 = Ifges(3,4) * t363;
t320 = -pkin(3) - t408;
t310 = mrSges(3,3) * t363 - t383;
t309 = -mrSges(3,3) * t364 + t384;
t306 = t318 - t408;
t285 = Ifges(3,1) * t364 + t322 + t386;
t284 = t385 + (t332 * Ifges(3,2) + t399) * qJD(1);
t267 = -mrSges(4,2) * t324 - t403;
t247 = mrSges(4,1) * t286 + t455;
t230 = t297 * t300;
t229 = t298 * t300;
t220 = t242 - t408;
t206 = mrSges(5,1) * t286 - mrSges(5,3) * t263;
t205 = -t286 * mrSges(5,2) - mrSges(5,3) * t339;
t199 = t253 - t270;
t194 = -t270 + t252;
t158 = mrSges(5,1) * t244 + mrSges(5,3) * t379;
t114 = -t257 * t298 + t358 * t373 - t359 * t374;
t113 = t257 * t297 - t283 * t300;
t106 = pkin(5) * t229 + qJ(6) * t230 + t215;
t103 = mrSges(6,1) * t182 + mrSges(6,2) * t471;
t102 = mrSges(7,1) * t182 - mrSges(7,3) * t471;
t101 = pkin(5) * t471 + qJ(6) * t182;
t46 = -pkin(5) * t299 - t48;
t45 = qJ(6) * t299 + t447;
t13 = pkin(5) * t114 - qJ(6) * t113 + qJD(6) * t230 + t134;
t8 = -pkin(5) * t258 - t10;
t7 = qJ(6) * t258 + qJD(6) * t299 + t9;
t2 = [t481 * t113 - t482 * t230 - (t321 * mrSges(4,2) + (Ifges(4,1) + Ifges(5,1) * t326 ^ 2 / 0.2e1 + (-t396 + t392 / 0.2e1) * t325) * t300) * t243 + (t167 * t300 + t243 * t446 - t244 * t266 + t251 * t257 - t252 * t258) * mrSges(4,3) + (mrSges(4,1) * t472 + mrSges(5,1) * t60 - mrSges(5,2) * t61 - mrSges(4,3) * t166 + Ifges(6,6) * t439 + Ifges(7,6) * t438 + t458 * t440 + t479 * t421 - (-Ifges(4,4) + t343) * t243 + t444 + t445 / 0.2e1 + (Ifges(4,2) + Ifges(5,3)) * t244) * t299 + m(7) * (t1 * t45 + t106 * t12 + t13 * t47 + t25 * t8 + t26 * t7 + t3 * t46) + m(6) * (t10 * t30 + t116 * t215 + t134 * t168 + t31 * t9 + t447 * t5 + t48 * t6) + t447 * t59 + (-t300 * t244 - t257 * t417 + t258 * t415) * Ifges(4,4) + (t90 / 0.2e1 + t91 / 0.2e1 + t30 * mrSges(6,1) - t25 * mrSges(7,1) + t26 * mrSges(7,3) - t31 * mrSges(6,2) + t162 / 0.2e1 - t223 / 0.2e1 + t311 * mrSges(4,1) - t389 / 0.2e1 - t139 * mrSges(5,2) + t138 * mrSges(5,1) + t394 / 0.2e1 - t334 / 0.2e1 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t286 + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t281 + (Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1) * t471 + (Ifges(7,6) / 0.2e1 + t465) * t182) * t258 - t446 * t153 + t321 * t244 * mrSges(4,1) + t487 * t229 + t488 * t114 + m(5) * (t138 * t77 + t139 * t78 + t173 * t60 + t174 * t61 + t187 * t217 - t381) + m(4) * (t166 * t266 + t186 * t252 - t187 * t251 - t381) + (-t284 / 0.2e1 - pkin(7) * t310 + t351 + (-0.2e1 * t434 - 0.3e1 / 0.2e1 * t399 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t332) * qJD(1) + (0.2e1 * t409 + t247 + t455) * pkin(2)) * t362 + t45 * t56 + t48 * t57 + t46 * t58 + (t110 * t411 + t109 * t412 + t167 * t346 + t343 * t421 + (-t325 * t61 - t326 * t60) * mrSges(5,3)) * t300 + t13 * t102 + t106 * t23 - (t391 / 0.2e1 + t349 + t350 + t224 / 0.2e1 + t311 * mrSges(4,2) + t390 / 0.2e1 + t449 + t263 * t345 / 0.2e1 + t343 * t416 - t339 * t344 / 0.2e1 + (-t138 * t326 - t139 * t325) * mrSges(5,3)) * t257 + t134 * t103 + t7 * t148 + t9 * t149 + t10 * t150 + t8 * t151 + t173 * t158 + t174 * t157 + t78 * t205 + t77 * t206 + (t285 / 0.2e1 - pkin(7) * t309 + t352 + (-0.2e1 * t433 + 0.3e1 / 0.2e1 * Ifges(3,4) * t332) * qJD(1)) * t332 * qJD(2) + t215 * t24 + t186 * t267 + t365 * t187; t485 * t149 + t486 * t150 + 0.2e1 * ((t166 * t328 - t167 * t331) * t443 + ((-t251 * t328 + t252 * t331) * t443 + (m(5) * t217 + m(6) * t168) * t328 / 0.2e1) * qJD(3)) * pkin(2) - m(4) * (-t251 * t253 + t252 * t254) + t317 * t369 + t473 * t102 + t474 * t151 + (t1 * t237 + t12 * t220 + t25 * t474 + t26 * t475 - t338 * t3 + t47 * t473) * m(7) + t475 * t148 - t405 * t338 + m(6) * (t116 * t306 + t171 * t31 - t172 * t30 + t237 * t5 + t338 * t6) + ((t243 * t331 - t244 * t328) * mrSges(4,3) + (t267 * t331 + (t103 + t365) * t328) * qJD(3)) * pkin(2) + t333 - m(6) * (t168 * t199 + t30 * t39 + t31 * t40) - m(5) * (t138 * t154 + t139 * t155 + t217 * t253) + ((-t322 / 0.2e1 - t285 / 0.2e1 + t352 + qJD(1) * t433 + (t309 - t384) * pkin(7)) * t332 + (t284 / 0.2e1 + t351 + (t434 + t399 / 0.2e1 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t332) * qJD(1) + (t310 + t383) * pkin(7) + (-t247 - t409) * pkin(2)) * t329) * qJD(1) + t404 * t237 + (-t317 * t158 - t314 * t206 - t406) * t325 + t252 * t402 - t199 * t103 - t154 * t206 + t220 * t23 - t254 * t267 + t306 * t24 + t320 * t153 - t365 * t253 + (-t155 + t371) * t205 + m(5) * (t139 * t371 + t167 * t320 - t314 * t382 - t317 * t388 + t370 * t61); t452 * t150 + (-qJ(4) * t158 - qJD(4) * t206 - t406) * t325 + t451 * t149 + t454 * t148 - t405 * t337 + t450 * t102 + t453 * t151 + (-t365 + t402) * t252 + t333 + (qJD(4) * t326 - t160) * t205 + qJ(4) * t369 + t404 * t261 - pkin(3) * t153 - t194 * t103 - t159 * t206 + t242 * t23 - t251 * t267 + t318 * t24 + (t1 * t261 + t12 * t242 + t25 * t453 + t26 * t454 - t3 * t337 + t450 * t47) * m(7) + (t116 * t318 - t168 * t194 + t261 * t5 + t30 * t452 + t31 * t451 + t337 * t6) * m(6) + (-t138 * t159 - t139 * t160 - t217 * t252 - pkin(3) * t167 + (t139 * t326 - t382) * qJD(4) + (t387 - t388) * qJ(4)) * m(5); t366 * t471 - t367 * t182 + t339 * t205 + t263 * t206 + t153 + t23 + t24 + (t182 * t26 - t25 * t471 + t12) * m(7) + (t182 * t31 + t30 * t471 + t116) * m(6) + (t138 * t263 + t139 * t339 + t167) * m(5); (-t182 * t460 + t179 - t395 + t89) * t425 + t444 + (-t182 * t458 - t471 * t478) * t420 + qJ(6) * t56 - pkin(5) * t58 + (t182 * t25 + t26 * t471) * mrSges(7,2) + (t366 + t400) * t31 + (t367 - t401) * t30 - t101 * t102 + (-Ifges(6,2) * t471 - t180 + t456) * t426 + qJD(6) * t148 + t92 * t424 + (Ifges(7,3) * t471 - t393) * t427 - t47 * (mrSges(7,1) * t471 + mrSges(7,3) * t182) - t168 * (mrSges(6,1) * t471 - mrSges(6,2) * t182) + (-pkin(5) * t3 + qJ(6) * t1 - t101 * t47 - t25 * t31 + t26 * t448) * m(7) + t445; t471 * t102 - t281 * t148 + 0.2e1 * (t3 / 0.2e1 + t47 * t424 + t26 * t420) * m(7) + t58;];
tauc  = t2(:);
