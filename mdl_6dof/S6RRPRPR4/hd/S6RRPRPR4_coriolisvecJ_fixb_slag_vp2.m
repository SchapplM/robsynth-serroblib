% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2018-11-23 17:02
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:01:53
% EndTime: 2018-11-23 17:02:14
% DurationCPUTime: 22.06s
% Computational Cost: add. (19977->790), mult. (60146->1118), div. (0->0), fcn. (48738->12), ass. (0->356)
t285 = sin(qJ(2));
t282 = cos(pkin(6));
t397 = pkin(1) * t282;
t273 = t285 * t397;
t281 = sin(pkin(6));
t288 = cos(qJ(2));
t360 = t281 * t288;
t389 = pkin(8) + qJ(3);
t453 = t360 * t389 + t273;
t222 = t453 * qJD(1);
t368 = sin(pkin(11));
t213 = t368 * t222;
t274 = t288 * t397;
t269 = qJD(1) * t274;
t336 = t389 * t285;
t320 = t281 * t336;
t221 = -qJD(1) * t320 + t269;
t370 = cos(pkin(11));
t170 = t221 * t370 - t213;
t294 = -t285 * t368 + t288 * t370;
t292 = t281 * t294;
t238 = qJD(1) * t292;
t354 = qJD(1) * t281;
t461 = t285 * t370 + t288 * t368;
t239 = t461 * t354;
t340 = t285 * t354;
t324 = pkin(2) * t340;
t187 = pkin(3) * t239 - pkin(9) * t238 + t324;
t284 = sin(qJ(4));
t287 = cos(qJ(4));
t110 = t170 * t287 + t187 * t284;
t341 = t368 * pkin(2);
t277 = t341 + pkin(9);
t358 = qJ(5) + t277;
t325 = qJD(4) * t358;
t365 = t238 * t284;
t484 = qJ(5) * t365 + qJD(5) * t287 - t284 * t325 - t110;
t109 = -t170 * t284 + t187 * t287;
t364 = t238 * t287;
t483 = -pkin(4) * t239 + qJ(5) * t364 - qJD(5) * t284 - t287 * t325 - t109;
t280 = sin(pkin(12));
t369 = cos(pkin(12));
t297 = -t280 * t284 + t287 * t369;
t182 = t297 * t238;
t255 = t297 * qJD(4);
t482 = t182 - t255;
t260 = t280 * t287 + t284 * t369;
t181 = t260 * t238;
t254 = t260 * qJD(4);
t462 = t254 - t181;
t466 = t280 * t483 + t369 * t484;
t332 = t370 * t222;
t169 = t221 * t368 + t332;
t351 = qJD(4) * t284;
t459 = -t169 + (t351 - t365) * pkin(4);
t235 = -t294 * t354 + qJD(4);
t408 = -t235 / 0.2e1;
t271 = qJD(1) * t282 + qJD(2);
t204 = -t239 * t284 + t271 * t287;
t305 = -t239 * t287 - t271 * t284;
t298 = t204 * t280 - t305 * t369;
t420 = -t298 / 0.2e1;
t326 = t204 * t369 + t280 * t305;
t422 = -t326 / 0.2e1;
t481 = -Ifges(6,4) * t420 - Ifges(6,2) * t422 - Ifges(6,6) * t408;
t137 = qJD(6) - t326;
t423 = t137 / 0.2e1;
t283 = sin(qJ(6));
t286 = cos(qJ(6));
t118 = t235 * t283 + t286 * t298;
t427 = t118 / 0.2e1;
t117 = t235 * t286 - t283 * t298;
t429 = t117 / 0.2e1;
t480 = Ifges(7,5) * t427 + Ifges(7,6) * t429 + Ifges(7,3) * t423;
t473 = Ifges(5,3) + Ifges(6,3);
t479 = -pkin(10) * t239 + t466;
t478 = pkin(5) * t462 + pkin(10) * t482 + t459;
t206 = pkin(2) * t271 + t221;
t149 = t206 * t370 - t213;
t143 = -pkin(3) * t271 - t149;
t115 = -pkin(4) * t204 + qJD(5) + t143;
t150 = t206 * t368 + t332;
t144 = pkin(9) * t271 + t150;
t261 = (-pkin(2) * t288 - pkin(1)) * t281;
t253 = qJD(1) * t261 + qJD(3);
t165 = -pkin(3) * t238 - pkin(9) * t239 + t253;
t99 = -t144 * t284 + t165 * t287;
t81 = qJ(5) * t305 + t99;
t71 = pkin(4) * t235 + t81;
t100 = t144 * t287 + t165 * t284;
t82 = qJ(5) * t204 + t100;
t76 = t369 * t82;
t35 = t280 * t71 + t76;
t31 = pkin(10) * t235 + t35;
t63 = -pkin(5) * t326 - pkin(10) * t298 + t115;
t16 = -t283 * t31 + t286 * t63;
t17 = t283 * t63 + t286 * t31;
t477 = mrSges(6,1) * t115 + mrSges(7,1) * t16 - mrSges(7,2) * t17 - mrSges(6,3) * t35 + t480 - t481;
t467 = -t280 * t484 + t369 * t483;
t407 = t235 / 0.2e1;
t419 = t298 / 0.2e1;
t421 = t326 / 0.2e1;
t476 = -Ifges(6,4) * t419 - Ifges(6,2) * t421 - Ifges(6,6) * t407 + t477 + t480;
t424 = -t137 / 0.2e1;
t428 = -t118 / 0.2e1;
t430 = -t117 / 0.2e1;
t475 = Ifges(7,5) * t428 + Ifges(7,6) * t430 + Ifges(7,3) * t424 - t477 + t481;
t241 = qJD(2) * t292;
t233 = qJD(1) * t241;
t154 = qJD(4) * t204 + t233 * t287;
t155 = qJD(4) * t305 - t233 * t284;
t106 = t154 * t369 + t155 * t280;
t245 = t461 * t281;
t240 = qJD(2) * t245;
t232 = qJD(1) * t240;
t61 = qJD(6) * t117 + t106 * t286 + t232 * t283;
t439 = t61 / 0.2e1;
t62 = -qJD(6) * t118 - t106 * t283 + t232 * t286;
t438 = t62 / 0.2e1;
t105 = t154 * t280 - t155 * t369;
t432 = t105 / 0.2e1;
t22 = -mrSges(7,1) * t62 + mrSges(7,2) * t61;
t90 = mrSges(6,1) * t232 - mrSges(6,3) * t106;
t472 = t22 - t90;
t343 = t370 * pkin(2);
t279 = -t343 - pkin(3);
t264 = -pkin(4) * t287 + t279;
t197 = -pkin(5) * t297 - pkin(10) * t260 + t264;
t258 = t358 * t287;
t327 = t358 * t284;
t203 = t258 * t369 - t280 * t327;
t136 = t197 * t283 + t203 * t286;
t471 = -qJD(6) * t136 - t283 * t479 + t286 * t478;
t135 = t197 * t286 - t203 * t283;
t470 = qJD(6) * t135 + t283 * t478 + t286 * t479;
t469 = -Ifges(5,5) * t305 + Ifges(6,5) * t298 + t204 * Ifges(5,6) + Ifges(6,6) * t326 + t235 * t473;
t468 = pkin(5) * t239 - t467;
t120 = mrSges(6,1) * t235 - mrSges(6,3) * t298;
t69 = -mrSges(7,1) * t117 + mrSges(7,2) * t118;
t465 = t69 - t120;
t133 = -t182 * t283 + t239 * t286;
t348 = qJD(6) * t286;
t302 = t255 * t283 + t260 * t348;
t464 = t133 + t302;
t134 = t182 * t286 + t239 * t283;
t349 = qJD(6) * t283;
t301 = -t255 * t286 + t260 * t349;
t463 = t134 + t301;
t387 = mrSges(4,3) * t239;
t357 = -mrSges(4,1) * t271 - mrSges(5,1) * t204 - mrSges(5,2) * t305 + t387;
t218 = pkin(2) * t282 + t274 - t320;
t355 = pkin(8) * t360 + t273;
t236 = qJ(3) * t360 + t355;
t180 = t218 * t368 + t236 * t370;
t168 = pkin(9) * t282 + t180;
t189 = -pkin(3) * t292 - pkin(9) * t245 + t261;
t112 = t168 * t287 + t189 * t284;
t460 = Ifges(5,5) * t154 + Ifges(6,5) * t106 + Ifges(5,6) * t155 - Ifges(6,6) * t105 + t232 * t473;
t266 = qJD(2) * t269;
t293 = (-qJD(2) * t336 + qJD(3) * t288) * t281;
t198 = qJD(1) * t293 + t266;
t361 = t281 * t285;
t208 = -qJD(2) * t453 - qJD(3) * t361;
t289 = qJD(1) * t208;
t131 = t198 * t370 + t289 * t368;
t352 = qJD(2) * t281;
t335 = qJD(1) * t352;
t319 = t285 * t335;
t308 = pkin(2) * t319;
t166 = pkin(3) * t232 - pkin(9) * t233 + t308;
t350 = qJD(4) * t287;
t56 = t131 * t287 - t144 * t351 + t165 * t350 + t166 * t284;
t57 = -qJD(4) * t100 - t131 * t284 + t166 * t287;
t458 = -t284 * t57 + t287 * t56;
t28 = mrSges(7,1) * t105 - mrSges(7,3) * t61;
t29 = -mrSges(7,2) * t105 + mrSges(7,3) * t62;
t457 = -t28 * t283 + t286 * t29;
t130 = t198 * t368 - t289 * t370;
t95 = -pkin(4) * t155 + t130;
t26 = pkin(5) * t105 - pkin(10) * t106 + t95;
t27 = pkin(4) * t232 - qJ(5) * t154 + qJD(5) * t305 + t57;
t33 = qJ(5) * t155 + qJD(5) * t204 + t56;
t8 = t27 * t280 + t33 * t369;
t6 = pkin(10) * t232 + t8;
t1 = qJD(6) * t16 + t26 * t283 + t286 * t6;
t2 = -qJD(6) * t17 + t26 * t286 - t283 * t6;
t456 = t1 * t286 - t2 * t283;
t455 = qJD(4) - t238;
t375 = t280 * t82;
t34 = t369 * t71 - t375;
t454 = t115 * mrSges(6,2) - t34 * mrSges(6,3);
t452 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t13 = Ifges(7,5) * t61 + Ifges(7,6) * t62 + Ifges(7,3) * t105;
t431 = t106 / 0.2e1;
t433 = -t105 / 0.2e1;
t450 = t95 * mrSges(6,1) - t8 * mrSges(6,3) - Ifges(6,4) * t431 + Ifges(7,5) * t439 - Ifges(6,2) * t433 + Ifges(7,6) * t438 + Ifges(7,3) * t432 + t13 / 0.2e1 + t452;
t447 = -0.2e1 * pkin(1);
t445 = Ifges(7,1) * t439 + Ifges(7,4) * t438 + Ifges(7,5) * t432;
t382 = Ifges(7,4) * t118;
t54 = Ifges(7,2) * t117 + Ifges(7,6) * t137 + t382;
t441 = t54 / 0.2e1;
t116 = Ifges(7,4) * t117;
t55 = Ifges(7,1) * t118 + Ifges(7,5) * t137 + t116;
t440 = -t55 / 0.2e1;
t80 = Ifges(6,1) * t298 + Ifges(6,4) * t326 + t235 * Ifges(6,5);
t435 = -t80 / 0.2e1;
t434 = t80 / 0.2e1;
t377 = t305 * Ifges(5,4);
t122 = t204 * Ifges(5,2) + Ifges(5,6) * t235 - t377;
t426 = t122 / 0.2e1;
t200 = Ifges(5,4) * t204;
t123 = -Ifges(5,1) * t305 + Ifges(5,5) * t235 + t200;
t425 = t123 / 0.2e1;
t418 = t154 / 0.2e1;
t417 = t155 / 0.2e1;
t216 = -t245 * t284 + t282 * t287;
t217 = t245 * t287 + t282 * t284;
t156 = -t216 * t369 + t217 * t280;
t416 = -t156 / 0.2e1;
t157 = t216 * t280 + t217 * t369;
t415 = t157 / 0.2e1;
t414 = -t204 / 0.2e1;
t413 = t305 / 0.2e1;
t412 = -t305 / 0.2e1;
t411 = t216 / 0.2e1;
t410 = t217 / 0.2e1;
t409 = t232 / 0.2e1;
t405 = t239 / 0.2e1;
t401 = t282 / 0.2e1;
t400 = -t285 / 0.2e1;
t399 = t286 / 0.2e1;
t396 = pkin(4) * t305;
t395 = pkin(4) * t280;
t177 = qJD(4) * t216 + t241 * t287;
t270 = qJD(2) * t274;
t207 = t270 + t293;
t147 = t207 * t370 + t208 * t368;
t338 = t285 * t352;
t323 = pkin(2) * t338;
t188 = pkin(3) * t240 - pkin(9) * t241 + t323;
t65 = -qJD(4) * t112 - t147 * t284 + t188 * t287;
t40 = pkin(4) * t240 - qJ(5) * t177 - qJD(5) * t217 + t65;
t178 = -qJD(4) * t217 - t241 * t284;
t64 = t147 * t287 - t168 * t351 + t188 * t284 + t189 * t350;
t46 = qJ(5) * t178 + qJD(5) * t216 + t64;
t12 = t280 * t40 + t369 * t46;
t111 = -t168 * t284 + t189 * t287;
t87 = -pkin(4) * t292 - qJ(5) * t217 + t111;
t97 = qJ(5) * t216 + t112;
t48 = t280 * t87 + t369 * t97;
t388 = mrSges(4,3) * t238;
t386 = mrSges(5,3) * t204;
t385 = Ifges(3,4) * t285;
t384 = Ifges(5,4) * t284;
t383 = Ifges(5,4) * t287;
t381 = Ifges(7,4) * t283;
t380 = Ifges(7,4) * t286;
t379 = Ifges(3,5) * t288;
t378 = t100 * mrSges(5,3);
t376 = t239 * Ifges(4,4);
t367 = t326 * t283;
t366 = t326 * t286;
t363 = t260 * t283;
t362 = t260 * t286;
t346 = t55 * t399;
t342 = t369 * pkin(4);
t339 = t288 * t354;
t333 = -t349 / 0.2e1;
t60 = t105 * mrSges(6,1) + mrSges(6,2) * t106;
t322 = mrSges(3,3) * t340;
t321 = mrSges(3,3) * t339;
t317 = mrSges(7,1) * t283 + mrSges(7,2) * t286;
t316 = Ifges(5,1) * t287 - t384;
t315 = Ifges(7,1) * t286 - t381;
t314 = -Ifges(5,2) * t284 + t383;
t313 = -Ifges(7,2) * t283 + t380;
t312 = Ifges(5,5) * t287 - Ifges(5,6) * t284;
t311 = Ifges(7,5) * t286 - Ifges(7,6) * t283;
t310 = -t16 * t283 + t17 * t286;
t43 = -pkin(10) * t292 + t48;
t179 = t218 * t370 - t236 * t368;
t167 = -pkin(3) * t282 - t179;
t128 = -pkin(4) * t216 + t167;
t68 = pkin(5) * t156 - pkin(10) * t157 + t128;
t19 = t283 * t68 + t286 * t43;
t18 = -t283 * t43 + t286 * t68;
t74 = -mrSges(7,2) * t137 + mrSges(7,3) * t117;
t75 = mrSges(7,1) * t137 - mrSges(7,3) * t118;
t309 = -t283 * t75 + t286 * t74;
t146 = t207 * t368 - t208 * t370;
t127 = t157 * t286 - t283 * t292;
t126 = -t157 * t283 - t286 * t292;
t163 = -mrSges(5,2) * t235 + t386;
t164 = mrSges(5,1) * t235 + mrSges(5,3) * t305;
t306 = t163 * t287 - t164 * t284;
t119 = -mrSges(6,2) * t235 + mrSges(6,3) * t326;
t304 = -t119 - t309;
t30 = -pkin(5) * t235 - t34;
t303 = t30 * t317;
t7 = t27 * t369 - t280 * t33;
t11 = -t280 * t46 + t369 * t40;
t47 = -t280 * t97 + t369 * t87;
t299 = t271 * (-Ifges(3,6) * t285 + t379);
t108 = -pkin(4) * t178 + t146;
t252 = t355 * qJD(2);
t242 = -pkin(8) * t319 + t266;
t243 = qJD(1) * t252;
t291 = -t243 * mrSges(3,1) - t130 * mrSges(4,1) - t242 * mrSges(3,2) - t131 * mrSges(4,2);
t290 = (-t16 * t286 - t17 * t283) * qJD(6) + t456;
t278 = -t342 - pkin(5);
t268 = Ifges(3,4) * t339;
t265 = t335 * t379;
t256 = -pkin(8) * t361 + t274;
t251 = -pkin(8) * t338 + t270;
t250 = t355 * qJD(1);
t249 = -pkin(8) * t340 + t269;
t247 = -mrSges(3,2) * t271 + t321;
t246 = mrSges(3,1) * t271 - t322;
t234 = Ifges(4,4) * t238;
t227 = Ifges(4,5) * t233;
t226 = Ifges(4,6) * t232;
t223 = t233 * mrSges(4,2);
t220 = Ifges(3,1) * t340 + Ifges(3,5) * t271 + t268;
t219 = Ifges(3,6) * t271 + (Ifges(3,2) * t288 + t385) * t354;
t209 = -mrSges(4,2) * t271 + t388;
t202 = t258 * t280 + t327 * t369;
t192 = -mrSges(4,1) * t238 + mrSges(4,2) * t239;
t186 = t239 * Ifges(4,1) + t271 * Ifges(4,5) + t234;
t185 = t238 * Ifges(4,2) + t271 * Ifges(4,6) + t376;
t125 = -mrSges(5,2) * t232 + mrSges(5,3) * t155;
t124 = mrSges(5,1) * t232 - mrSges(5,3) * t154;
t114 = t177 * t369 + t178 * t280;
t113 = t177 * t280 - t178 * t369;
t107 = -mrSges(5,1) * t155 + mrSges(5,2) * t154;
t93 = Ifges(5,1) * t154 + Ifges(5,4) * t155 + Ifges(5,5) * t232;
t92 = Ifges(5,4) * t154 + Ifges(5,2) * t155 + Ifges(5,6) * t232;
t91 = -mrSges(6,1) * t326 + mrSges(6,2) * t298;
t89 = -mrSges(6,2) * t232 - mrSges(6,3) * t105;
t73 = pkin(5) * t298 - pkin(10) * t326 - t396;
t67 = -qJD(6) * t127 - t114 * t283 + t240 * t286;
t66 = qJD(6) * t126 + t114 * t286 + t240 * t283;
t52 = t106 * Ifges(6,1) - t105 * Ifges(6,4) + t232 * Ifges(6,5);
t51 = t106 * Ifges(6,4) - t105 * Ifges(6,2) + t232 * Ifges(6,6);
t42 = pkin(5) * t292 - t47;
t38 = t369 * t81 - t375;
t37 = t280 * t81 + t76;
t36 = pkin(5) * t113 - pkin(10) * t114 + t108;
t21 = t283 * t73 + t286 * t38;
t20 = -t283 * t38 + t286 * t73;
t14 = Ifges(7,4) * t61 + Ifges(7,2) * t62 + Ifges(7,6) * t105;
t10 = pkin(10) * t240 + t12;
t9 = -pkin(5) * t240 - t11;
t5 = -pkin(5) * t232 - t7;
t4 = -qJD(6) * t19 - t10 * t283 + t286 * t36;
t3 = qJD(6) * t18 + t10 * t286 + t283 * t36;
t15 = [(Ifges(7,5) * t66 + Ifges(7,6) * t67) * t423 + (Ifges(7,5) * t127 + Ifges(7,6) * t126) * t432 + t261 * t223 + (t288 * t220 / 0.2e1 + t219 * t400 + t299 / 0.2e1 + t285 * pkin(2) * t192 + (-t249 * t288 - t250 * t285) * mrSges(3,3)) * t352 + (Ifges(7,1) * t66 + Ifges(7,4) * t67) * t427 + (Ifges(7,1) * t127 + Ifges(7,4) * t126) * t439 + (t1 * t126 - t127 * t2 - t16 * t66 + t17 * t67) * mrSges(7,3) + t357 * t146 + (t242 * t288 + t243 * t285) * t281 * mrSges(3,3) + (-mrSges(4,3) * t179 + Ifges(4,1) * t245 + Ifges(4,4) * t292 + Ifges(4,5) * t401) * t233 + (-Ifges(4,4) * t245 - Ifges(4,6) * t282 / 0.2e1 + t261 * mrSges(4,1) + Ifges(6,5) * t415 + Ifges(6,6) * t416 + Ifges(5,5) * t410 + Ifges(5,6) * t411 - t180 * mrSges(4,3) - (Ifges(4,2) + Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t292) * t232 + ((-t256 * mrSges(3,3) + Ifges(3,5) * t401 + (mrSges(3,2) * t447 + 0.3e1 / 0.2e1 * Ifges(3,4) * t288) * t281) * t288 + (-t355 * mrSges(3,3) - Ifges(3,6) * t282 + (mrSges(3,1) * t447 - 0.3e1 / 0.2e1 * t385 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t288) * t281 + (m(4) * t261 - mrSges(4,1) * t292 + mrSges(4,2) * t245) * pkin(2)) * t285) * t335 + (Ifges(6,1) * t157 - Ifges(6,5) * t292) * t431 + (t114 * t115 + t157 * t95 - t240 * t35 + t292 * t8) * mrSges(6,2) + (Ifges(6,4) * t157 - Ifges(6,6) * t292) * t433 + (t130 * t245 + t131 * t292 - t149 * t241 - t150 * t240) * mrSges(4,3) + t7 * (-mrSges(6,1) * t292 - mrSges(6,3) * t157) + t56 * (mrSges(5,2) * t292 + mrSges(5,3) * t216) + t57 * (-mrSges(5,1) * t292 - mrSges(5,3) * t217) + (Ifges(5,4) * t217 + Ifges(5,2) * t216 - Ifges(5,6) * t292) * t417 + (Ifges(5,1) * t217 + Ifges(5,4) * t216 - Ifges(5,5) * t292) * t418 - t460 * t292 / 0.2e1 + t450 * t156 + (Ifges(7,4) * t66 + Ifges(7,2) * t67) * t429 + (Ifges(7,4) * t127 + Ifges(7,2) * t126) * t438 + (Ifges(6,1) * t114 + Ifges(6,5) * t240) * t419 + t476 * t113 + m(7) * (t1 * t19 + t16 * t4 + t17 * t3 + t18 * t2 + t30 * t9 + t42 * t5) + m(6) * (t108 * t115 + t11 * t34 + t12 * t35 + t128 * t95 + t47 * t7 + t48 * t8) + m(5) * (t100 * t64 + t111 * t57 + t112 * t56 + t130 * t167 + t143 * t146 + t65 * t99) + (Ifges(6,4) * t114 + Ifges(6,6) * t240) * t421 + m(3) * (t242 * t355 - t243 * t256 - t249 * t252 + t250 * t251) + m(4) * (-t130 * t179 + t131 * t180 - t146 * t149 + t147 * t150 + t253 * t323) + t271 * (Ifges(4,5) * t241 - Ifges(4,6) * t240) / 0.2e1 + t251 * t247 - t252 * t246 + t253 * (mrSges(4,1) * t240 + mrSges(4,2) * t241) - t240 * t185 / 0.2e1 + t238 * (Ifges(4,4) * t241 - Ifges(4,2) * t240) / 0.2e1 + t241 * t186 / 0.2e1 + t34 * (mrSges(6,1) * t240 - mrSges(6,3) * t114) + t99 * (mrSges(5,1) * t240 - mrSges(5,3) * t177) + t100 * (-mrSges(5,2) * t240 + mrSges(5,3) * t178) + t204 * (Ifges(5,4) * t177 + Ifges(5,2) * t178 + Ifges(5,6) * t240) / 0.2e1 + (Ifges(4,1) * t241 - Ifges(4,4) * t240) * t405 + t93 * t410 + t92 * t411 + (Ifges(5,1) * t177 + Ifges(5,4) * t178 + Ifges(5,5) * t240) * t412 + t52 * t415 + t51 * t416 + t177 * t425 + t178 * t426 + t114 * t434 + t67 * t441 + t127 * t445 + t469 * t240 / 0.2e1 + t18 * t28 + t19 * t29 + (Ifges(5,5) * t177 + Ifges(6,5) * t114 + Ifges(5,6) * t178 + t240 * t473) * t407 + t42 * t22 + t66 * t55 / 0.2e1 + t30 * (-mrSges(7,1) * t67 + mrSges(7,2) * t66) + t9 * t69 + t3 * t74 + t4 * t75 + (t265 / 0.2e1 + t227 / 0.2e1 - t226 / 0.2e1 + t291) * t282 + t48 * t89 + t47 * t90 + t108 * t91 + t12 * t119 + t11 * t120 + t111 * t124 + t112 * t125 + t126 * t14 / 0.2e1 + t5 * (-mrSges(7,1) * t126 + mrSges(7,2) * t127) + t128 * t60 + t64 * t163 + t65 * t164 + t167 * t107 + t143 * (-mrSges(5,1) * t178 + mrSges(5,2) * t177) + t147 * t209 + t130 * (-mrSges(5,1) * t216 + mrSges(5,2) * t217); (-t122 / 0.2e1 - t378) * t351 + (pkin(1) * (mrSges(3,1) * t285 + mrSges(3,2) * t288) + (Ifges(3,1) * t288 - t385) * t400) * qJD(1) ^ 2 * t281 ^ 2 + t265 + (t149 * t169 - t150 * t170 - t253 * t324 + (-t130 * t370 + t131 * t368) * pkin(2)) * m(4) + t150 * t387 - t226 + t227 - t14 * t363 / 0.2e1 + t455 * t143 * (mrSges(5,1) * t284 + mrSges(5,2) * t287) + (t321 - t247) * t249 + (t55 * t333 + t52 / 0.2e1 + t95 * mrSges(6,2) + t5 * t317 + Ifges(6,5) * t409 + Ifges(6,1) * t431 + t311 * t432 + Ifges(6,4) * t433 + t313 * t438 + t315 * t439 - t7 * mrSges(6,3)) * t260 + t291 + (-Ifges(7,5) * t301 - Ifges(7,6) * t302) * t423 + (Ifges(7,5) * t134 + Ifges(7,6) * t133) * t424 - t99 * mrSges(5,1) * t239 + t100 * mrSges(5,2) * t239 + t219 * t340 / 0.2e1 + t475 * t181 + t476 * t254 + (-t100 * t110 - t109 * t99 + t130 * t279 - t143 * t169) * m(5) + (-Ifges(7,4) * t301 - Ifges(7,2) * t302) * t429 + (Ifges(7,4) * t134 + Ifges(7,2) * t133) * t430 - (-Ifges(4,2) * t239 + t186 + t234) * t238 / 0.2e1 + (-t115 * t182 + t239 * t35) * mrSges(6,2) + (Ifges(6,4) * t182 + Ifges(6,6) * t239) * t422 - t271 * (Ifges(4,5) * t238 - Ifges(4,6) * t239) / 0.2e1 - t253 * (mrSges(4,1) * t239 + mrSges(4,2) * t238) - t34 * (mrSges(6,1) * t239 - mrSges(6,3) * t182) + (Ifges(5,5) * t239 + t238 * t316) * t413 + (Ifges(5,6) * t239 + t238 * t314) * t414 + (Ifges(6,5) * t182 + t238 * t312 + t239 * t473) * t408 + (Ifges(6,1) * t182 + Ifges(6,5) * t239) * t420 - (Ifges(4,1) * t238 - t376 + t469) * t239 / 0.2e1 - t123 * t364 / 0.2e1 - t192 * t324 + (-Ifges(7,1) * t301 - Ifges(7,4) * t302) * t427 + (Ifges(7,1) * t134 + Ifges(7,4) * t133) * t428 - Ifges(3,6) * t319 + (t322 + t246) * t250 + t130 * (-mrSges(5,1) * t287 + mrSges(5,2) * t284) + t287 * t92 / 0.2e1 + t284 * t93 / 0.2e1 + t279 * t107 + t264 * t60 + t149 * t388 + (-t232 * t341 - t233 * t343) * mrSges(4,3) + (Ifges(6,1) * t419 + Ifges(6,4) * t421 + Ifges(6,5) * t407 + t346 + t434 + t454) * t255 + (m(5) * ((-t100 * t284 - t287 * t99) * qJD(4) + t458) - t164 * t350 - t163 * t351 - t284 * t124 + t287 * t125) * t277 + (t100 * t365 + (-t350 + t364) * t99 + t458) * mrSges(5,3) + t459 * t91 - t357 * t169 - t464 * t54 / 0.2e1 + (mrSges(7,1) * t464 - mrSges(7,2) * t463) * t30 + (-t1 * t363 + t16 * t463 - t17 * t464 - t2 * t362) * mrSges(7,3) + t185 * t405 + (Ifges(5,5) * t284 + Ifges(5,6) * t287) * t409 + (Ifges(5,2) * t287 + t384) * t417 + (Ifges(5,1) * t284 + t383) * t418 + t350 * t425 + t365 * t426 + t182 * t435 + t134 * t440 + t362 * t445 + t466 * t119 + t467 * t120 + (t115 * t459 - t202 * t7 + t203 * t8 + t264 * t95 + t34 * t467 + t35 * t466) * m(6) + t468 * t69 - ((-Ifges(3,2) * t340 + t220 + t268) * t288 + t299) * t354 / 0.2e1 + (t204 * t314 + t235 * t312 - t305 * t316) * qJD(4) / 0.2e1 - (-t51 / 0.2e1 - Ifges(6,6) * t409 + t450) * t297 + t470 * t74 + t471 * t75 + (t1 * t136 + t135 * t2 + t16 * t471 + t17 * t470 + t202 * t5 + t30 * t468) * m(7) + t472 * t202 + t135 * t28 + t136 * t29 - t110 * t163 - t109 * t164 + t203 * t89 - t170 * t209; t232 * mrSges(4,1) - t182 * t119 + t287 * t124 + t284 * t125 - t133 * t75 - t134 * t74 + t223 - t472 * t297 + t306 * qJD(4) - t304 * t255 - (t91 + t357) * t239 + (-t209 - t306) * t238 + (t89 + (-t283 * t74 - t286 * t75) * qJD(6) + t457) * t260 + t465 * t462 + (-t133 * t16 - t134 * t17 + t255 * t310 + t260 * t290 - t297 * t5 + t30 * t462) * m(7) + (-t115 * t239 + t260 * t8 + t297 * t7 - t462 * t34 - t35 * t482) * m(6) + (-t143 * t239 + t284 * t56 + t287 * t57 + t455 * (t100 * t287 - t284 * t99)) * m(5) + (t149 * t239 - t150 * t238 + t308) * m(4); t91 * t396 + t460 + (t117 * t313 + t118 * t315 + t137 * t311) * qJD(6) / 0.2e1 + (t346 + t303) * qJD(6) + (-t16 * t20 - t17 * t21 + t278 * t5 - t30 * t37) * m(7) + t475 * t298 + (t386 - t163) * t99 + t54 * t333 + t90 * t342 + t5 * (-mrSges(7,1) * t286 + mrSges(7,2) * t283) + t278 * t22 + (Ifges(5,1) * t204 + t377) * t413 + (t115 * t396 + t34 * t37 - t35 * t38 + (t280 * t8 + t369 * t7) * pkin(4)) * m(6) + (Ifges(6,1) * t420 + Ifges(6,4) * t422 + Ifges(6,5) * t408 + t311 * t424 + t313 * t430 + t315 * t428 - t303 + t435 - t454) * t326 + t7 * mrSges(6,1) - t8 * mrSges(6,2) + ((-t349 + t367) * t17 + (-t348 + t366) * t16 + t456) * mrSges(7,3) + (m(7) * t290 - t348 * t75 - t349 * t74 + t457) * (pkin(10) + t395) + t89 * t395 + t14 * t399 + t122 * t412 + (Ifges(7,5) * t283 + Ifges(7,6) * t286) * t432 + (Ifges(7,2) * t286 + t381) * t438 + (Ifges(7,1) * t283 + t380) * t439 + t366 * t440 + t367 * t441 + t283 * t445 - t465 * t37 + (Ifges(5,2) * t305 + t123 + t200) * t414 + (Ifges(5,5) * t204 + Ifges(5,6) * t305) * t408 - t143 * (-mrSges(5,1) * t305 + mrSges(5,2) * t204) - t305 * t378 - t56 * mrSges(5,2) + t57 * mrSges(5,1) - t21 * t74 - t20 * t75 - t38 * t119 + t100 * t164; t286 * t28 + t283 * t29 - t465 * t298 + t309 * qJD(6) + t304 * t326 + t60 + (t1 * t283 + t137 * t310 + t2 * t286 - t298 * t30) * m(7) + (t298 * t34 - t326 * t35 + t95) * m(6); -t30 * (mrSges(7,1) * t118 + mrSges(7,2) * t117) + (Ifges(7,1) * t117 - t382) * t428 + t54 * t427 + (Ifges(7,5) * t117 - Ifges(7,6) * t118) * t424 - t16 * t74 + t17 * t75 + (t117 * t16 + t118 * t17) * mrSges(7,3) + t13 + (-Ifges(7,2) * t118 + t116 + t55) * t430 + t452;];
tauc  = t15(:);
