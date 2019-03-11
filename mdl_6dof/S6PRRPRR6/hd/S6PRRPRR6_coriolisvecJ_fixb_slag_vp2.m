% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRPRR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:22:54
% EndTime: 2019-03-08 22:23:39
% DurationCPUTime: 22.40s
% Computational Cost: add. (13191->703), mult. (36072->1025), div. (0->0), fcn. (29809->14), ass. (0->353)
t264 = cos(pkin(7));
t272 = cos(qJ(3));
t273 = cos(qJ(2));
t353 = t272 * t273;
t268 = sin(qJ(3));
t269 = sin(qJ(2));
t356 = t268 * t269;
t287 = -t264 * t356 + t353;
t262 = sin(pkin(6));
t349 = qJD(1) * t262;
t197 = t287 * t349;
t261 = sin(pkin(7));
t344 = qJD(3) * t268;
t332 = t261 * t344;
t343 = qJD(3) * t272;
t333 = t264 * t343;
t219 = pkin(2) * t333 - pkin(9) * t332;
t467 = -t197 + t219;
t304 = pkin(3) * t268 - qJ(4) * t272;
t278 = qJD(3) * t304 - qJD(4) * t268;
t337 = t269 * t349;
t328 = t261 * t337;
t466 = t278 * t261 - t328;
t465 = -qJD(4) * t264 - t467;
t260 = sin(pkin(13));
t263 = cos(pkin(13));
t435 = t465 * t260 + t263 * t466;
t434 = t260 * t466 - t465 * t263;
t358 = t263 * t272;
t285 = (pkin(4) * t268 - pkin(10) * t358) * t261;
t281 = qJD(3) * t285;
t464 = t281 + t435;
t334 = t261 * t343;
t325 = t260 * t334;
t463 = -pkin(10) * t325 + t434;
t267 = sin(qJ(5));
t271 = cos(qJ(5));
t234 = t260 * t271 + t263 * t267;
t359 = t261 * t272;
t284 = t234 * t359;
t187 = qJD(2) * t284;
t228 = t234 * qJD(5);
t352 = t187 - t228;
t233 = t260 * t267 - t271 * t263;
t283 = t233 * t359;
t188 = qJD(2) * t283;
t227 = t233 * qJD(5);
t351 = -t188 + t227;
t254 = qJD(2) * t264 + qJD(3);
t347 = qJD(2) * t261;
t336 = t268 * t347;
t205 = t254 * t263 - t260 * t336;
t206 = t254 * t260 + t263 * t336;
t329 = t271 * t205 - t206 * t267;
t455 = t329 / 0.2e1;
t335 = t272 * t347;
t249 = qJD(5) - t335;
t457 = t249 / 0.2e1;
t448 = Ifges(6,4) * t455 + Ifges(6,5) * t457;
t296 = t205 * t267 + t271 * t206;
t456 = t296 / 0.2e1;
t451 = Ifges(6,1) * t456;
t420 = t451 + t448;
t265 = cos(pkin(6));
t348 = qJD(1) * t265;
t338 = t261 * t348;
t231 = pkin(9) * t347 + t337;
t243 = qJD(2) * pkin(2) + t273 * t349;
t357 = t264 * t268;
t433 = t272 * t231 + t243 * t357;
t157 = t268 * t338 + t433;
t141 = qJ(4) * t254 + t157;
t253 = t264 * t348;
t305 = -pkin(3) * t272 - qJ(4) * t268;
t170 = t253 + (qJD(2) * t305 - t243) * t261;
t85 = -t141 * t260 + t263 * t170;
t68 = -pkin(4) * t335 - pkin(10) * t206 + t85;
t86 = t263 * t141 + t260 * t170;
t71 = pkin(10) * t205 + t86;
t32 = t267 * t68 + t271 * t71;
t26 = pkin(11) * t249 + t32;
t266 = sin(qJ(6));
t270 = cos(qJ(6));
t156 = -t268 * t231 + t272 * (t243 * t264 + t338);
t138 = -pkin(3) * t254 + qJD(4) - t156;
t104 = -pkin(4) * t205 + t138;
t41 = -pkin(5) * t329 - pkin(11) * t296 + t104;
t11 = -t26 * t266 + t270 * t41;
t279 = qJD(3) * t283;
t102 = -qJD(2) * t279 + qJD(5) * t329;
t280 = qJD(3) * t284;
t103 = qJD(2) * t280 + qJD(5) * t296;
t354 = t269 * t272;
t355 = t268 * t273;
t289 = t264 * t354 + t355;
t282 = t289 * qJD(2);
t324 = t265 * t332;
t109 = t433 * qJD(3) + (t262 * t282 + t324) * qJD(1);
t345 = qJD(3) * t261;
t330 = qJD(2) * t345;
t321 = t272 * t330;
t293 = t260 * t321;
t93 = pkin(4) * t293 + t109;
t33 = pkin(5) * t103 - pkin(11) * t102 + t93;
t322 = t268 * t330;
t341 = qJD(5) * t271;
t342 = qJD(5) * t267;
t167 = (t278 + t337) * t347;
t323 = t265 * t334;
t346 = qJD(2) * t262;
t331 = qJD(1) * t346;
t108 = t243 * t333 + qJD(1) * t323 + t331 * t353 + (-t264 * t269 * t331 - qJD(3) * t231) * t268;
t96 = qJD(4) * t254 + t108;
t69 = t263 * t167 - t260 * t96;
t55 = qJD(2) * t281 + t69;
t70 = t260 * t167 + t263 * t96;
t61 = -pkin(10) * t293 + t70;
t7 = t267 * t55 + t271 * t61 + t68 * t341 - t342 * t71;
t5 = pkin(11) * t322 + t7;
t1 = qJD(6) * t11 + t266 * t33 + t270 * t5;
t462 = t1 * mrSges(7,2);
t12 = t26 * t270 + t266 * t41;
t2 = -qJD(6) * t12 - t266 * t5 + t270 * t33;
t461 = t2 * mrSges(7,1);
t218 = t304 * t347;
t116 = t263 * t156 + t260 * t218;
t327 = t260 * t335;
t101 = -pkin(10) * t327 + t116;
t396 = pkin(10) + qJ(4);
t245 = t396 * t260;
t246 = t396 * t263;
t294 = -t271 * t245 - t246 * t267;
t115 = -t156 * t260 + t263 * t218;
t90 = qJD(2) * t285 + t115;
t445 = -qJD(4) * t233 + qJD(5) * t294 - t271 * t101 - t267 * t90;
t147 = qJD(6) - t329;
t378 = t147 * Ifges(7,3);
t125 = t249 * t266 + t270 * t296;
t379 = t125 * Ifges(7,5);
t124 = t249 * t270 - t266 * t296;
t380 = t124 * Ifges(7,6);
t36 = t378 + t379 + t380;
t370 = t249 * Ifges(6,6);
t377 = t329 * Ifges(6,2);
t389 = Ifges(6,4) * t296;
t78 = t370 + t377 + t389;
t460 = -t78 / 0.2e1 + t36 / 0.2e1;
t31 = -t267 * t71 + t271 * t68;
t459 = -t104 * mrSges(6,2) + t31 * mrSges(6,3);
t458 = -t104 * mrSges(6,1) - t11 * mrSges(7,1) + t12 * mrSges(7,2) + t32 * mrSges(6,3);
t230 = pkin(2) * t357 + pkin(9) * t359;
t210 = qJ(4) * t264 + t230;
t211 = (-pkin(2) + t305) * t261;
t154 = -t210 * t260 + t263 * t211;
t360 = t261 * t268;
t226 = t260 * t264 + t263 * t360;
t114 = -pkin(4) * t359 - pkin(10) * t226 + t154;
t155 = t263 * t210 + t260 * t211;
t224 = -t260 * t360 + t263 * t264;
t127 = pkin(10) * t224 + t155;
t438 = t114 * t341 - t127 * t342 + t267 * t464 + t463 * t271;
t454 = -pkin(11) * t336 + t445;
t134 = pkin(4) * t327 + t157;
t453 = -pkin(5) * t352 + pkin(11) * t351 - t134;
t220 = t230 * qJD(3);
t184 = pkin(4) * t325 + t220;
t196 = t289 * t349;
t452 = t184 - t196;
t450 = pkin(11) * t332 + t438;
t186 = -t245 * t267 + t246 * t271;
t444 = -qJD(4) * t234 - qJD(5) * t186 + t101 * t267 - t271 * t90;
t295 = t271 * t224 - t226 * t267;
t121 = qJD(5) * t295 - t279;
t172 = t224 * t267 + t226 * t271;
t122 = qJD(5) * t172 + t280;
t449 = pkin(5) * t122 - pkin(11) * t121 + t452;
t259 = -pkin(4) * t263 - pkin(3);
t177 = pkin(5) * t233 - pkin(11) * t234 + t259;
t118 = t177 * t270 - t186 * t266;
t447 = qJD(6) * t118 + t266 * t453 + t270 * t454;
t119 = t177 * t266 + t186 * t270;
t446 = -qJD(6) * t119 - t266 * t454 + t270 * t453;
t436 = t267 * t114 + t271 * t127;
t439 = -qJD(5) * t436 - t463 * t267 + t271 * t464;
t443 = pkin(5) * t336 - t444;
t59 = -pkin(11) * t359 + t436;
t255 = pkin(9) * t360;
t401 = pkin(2) * t272;
t213 = t255 + (-pkin(3) - t401) * t264;
t173 = -pkin(4) * t224 + t213;
t76 = -pkin(5) * t295 - pkin(11) * t172 + t173;
t28 = t266 * t76 + t270 * t59;
t442 = -qJD(6) * t28 - t266 * t450 + t270 * t449;
t27 = -t266 * t59 + t270 * t76;
t441 = qJD(6) * t27 + t266 * t449 + t270 * t450;
t440 = -pkin(5) * t332 - t439;
t303 = t11 * t270 + t12 * t266;
t437 = t303 * mrSges(7,3);
t432 = t264 * t353 - t356;
t431 = t1 * t270 - t2 * t266;
t25 = -pkin(5) * t249 - t31;
t307 = Ifges(7,5) * t270 - Ifges(7,6) * t266;
t386 = Ifges(7,4) * t270;
t309 = -Ifges(7,2) * t266 + t386;
t387 = Ifges(7,4) * t266;
t312 = Ifges(7,1) * t270 - t387;
t314 = mrSges(7,1) * t266 + mrSges(7,2) * t270;
t388 = Ifges(7,4) * t125;
t37 = Ifges(7,2) * t124 + Ifges(7,6) * t147 + t388;
t120 = Ifges(7,4) * t124;
t38 = Ifges(7,1) * t125 + Ifges(7,5) * t147 + t120;
t403 = t270 / 0.2e1;
t406 = -t266 / 0.2e1;
t413 = t147 / 0.2e1;
t415 = t125 / 0.2e1;
t417 = t124 / 0.2e1;
t276 = t25 * t314 + t307 * t413 + t309 * t417 + t312 * t415 + t37 * t406 + t38 * t403;
t430 = t420 + t276 + t448;
t8 = -qJD(5) * t32 - t267 * t61 + t271 * t55;
t429 = -t8 * mrSges(6,1) + t7 * mrSges(6,2) - Ifges(6,5) * t102 + Ifges(6,6) * t103;
t369 = t249 * Ifges(6,3);
t374 = t296 * Ifges(6,5);
t376 = t329 * Ifges(6,6);
t77 = t369 + t374 + t376;
t428 = Ifges(6,2) / 0.2e1;
t51 = -qJD(6) * t125 - t102 * t266 + t270 * t322;
t47 = Ifges(7,6) * t51;
t50 = qJD(6) * t124 + t102 * t270 + t266 * t322;
t48 = Ifges(7,5) * t50;
t13 = Ifges(7,3) * t103 + t47 + t48;
t427 = t13 / 0.2e1;
t15 = t50 * Ifges(7,1) + t51 * Ifges(7,4) + t103 * Ifges(7,5);
t426 = t15 / 0.2e1;
t424 = -t37 / 0.2e1;
t423 = t50 / 0.2e1;
t422 = t51 / 0.2e1;
t419 = t103 / 0.2e1;
t418 = -t124 / 0.2e1;
t416 = -t125 / 0.2e1;
t414 = -t147 / 0.2e1;
t412 = t295 / 0.2e1;
t411 = t172 / 0.2e1;
t410 = t224 / 0.2e1;
t409 = t226 / 0.2e1;
t408 = -t260 / 0.2e1;
t407 = t263 / 0.2e1;
t404 = t268 / 0.2e1;
t398 = t31 * mrSges(6,1);
t397 = t32 * mrSges(6,2);
t18 = -mrSges(7,1) * t51 + mrSges(7,2) * t50;
t87 = mrSges(6,1) * t322 - mrSges(6,3) * t102;
t395 = -t87 + t18;
t394 = mrSges(4,2) * t268;
t393 = mrSges(5,2) * t263;
t392 = Ifges(4,4) * t268;
t391 = Ifges(5,4) * t260;
t390 = Ifges(5,4) * t263;
t385 = Ifges(5,5) * t268;
t384 = Ifges(5,6) * t268;
t383 = t102 * Ifges(6,1);
t382 = t102 * Ifges(6,4);
t381 = t103 * Ifges(6,4);
t373 = t205 * Ifges(5,6);
t372 = t206 * Ifges(5,5);
t368 = t254 * Ifges(4,5);
t367 = t254 * Ifges(4,6);
t366 = t272 * Ifges(4,2);
t129 = mrSges(6,1) * t249 - mrSges(6,3) * t296;
t65 = -mrSges(7,1) * t124 + mrSges(7,2) * t125;
t365 = t129 - t65;
t191 = mrSges(5,1) * t293 + t321 * t393;
t49 = t103 * mrSges(6,1) + t102 * mrSges(6,2);
t364 = t191 + t49;
t180 = -t262 * t432 - t265 * t359;
t363 = t109 * t180;
t361 = t260 * t272;
t350 = mrSges(4,1) * t254 + mrSges(5,1) * t205 - mrSges(5,2) * t206 - mrSges(4,3) * t336;
t83 = -mrSges(6,1) * t329 + mrSges(6,2) * t296;
t340 = -t83 + t350;
t326 = t261 * t269 * t346;
t319 = t461 - t462;
t318 = -t1 * t266 - t2 * t270;
t317 = -t109 * mrSges(4,1) - t108 * mrSges(4,2);
t316 = mrSges(4,1) * t268 + mrSges(4,2) * t272;
t315 = mrSges(7,1) * t270 - mrSges(7,2) * t266;
t313 = Ifges(5,1) * t263 - t391;
t311 = Ifges(7,1) * t266 + t386;
t310 = -Ifges(5,2) * t260 + t390;
t308 = Ifges(7,2) * t270 + t387;
t306 = Ifges(7,5) * t266 + Ifges(7,6) * t270;
t302 = t11 * t266 - t12 * t270;
t72 = -mrSges(7,2) * t147 + mrSges(7,3) * t124;
t73 = mrSges(7,1) * t147 - mrSges(7,3) * t125;
t301 = -t266 * t73 + t270 * t72;
t288 = t264 * t355 + t354;
t181 = t262 * t288 + t265 * t360;
t225 = -t261 * t262 * t273 + t264 * t265;
t142 = -t181 * t260 + t225 * t263;
t143 = t181 * t263 + t225 * t260;
t81 = t142 * t267 + t143 * t271;
t56 = t180 * t270 - t266 * t81;
t57 = t180 * t266 + t270 * t81;
t62 = t114 * t271 - t127 * t267;
t297 = t271 * t142 - t143 * t267;
t292 = mrSges(5,1) * t268 - mrSges(5,3) * t358;
t291 = -mrSges(5,2) * t268 - mrSges(5,3) * t361;
t144 = -t172 * t266 - t270 * t359;
t290 = -t172 * t270 + t266 * t359;
t277 = -t389 / 0.2e1 + t380 / 0.2e1 + t379 / 0.2e1 + t378 / 0.2e1 - t370 / 0.2e1 + t460;
t251 = Ifges(4,4) * t335;
t248 = Ifges(4,5) * t321;
t247 = Ifges(6,3) * t322;
t229 = t264 * t401 - t255;
t217 = (-mrSges(4,1) * t272 + t394) * t347;
t216 = -mrSges(4,2) * t254 + mrSges(4,3) * t335;
t209 = t316 * t330;
t200 = -t243 * t261 + t253;
t199 = t292 * t330;
t198 = t291 * t330;
t190 = Ifges(4,1) * t336 + t251 + t368;
t189 = t367 + (t366 + t392) * t347;
t179 = -mrSges(5,1) * t335 - mrSges(5,3) * t206;
t178 = mrSges(5,2) * t335 + mrSges(5,3) * t205;
t175 = (t272 * t313 + t385) * t330;
t174 = (t272 * t310 + t384) * t330;
t169 = -t188 * t270 + t266 * t336;
t168 = t188 * t266 + t270 * t336;
t137 = Ifges(5,1) * t206 + Ifges(5,4) * t205 - Ifges(5,5) * t335;
t136 = Ifges(5,4) * t206 + Ifges(5,2) * t205 - Ifges(5,6) * t335;
t135 = -Ifges(5,3) * t335 + t372 + t373;
t133 = t323 + (t287 * qJD(2) + qJD(3) * t432) * t262;
t132 = t324 + (qJD(3) * t288 + t282) * t262;
t128 = -mrSges(6,2) * t249 + mrSges(6,3) * t329;
t112 = t133 * t263 + t260 * t326;
t111 = -t133 * t260 + t263 * t326;
t88 = -mrSges(6,2) * t322 - mrSges(6,3) * t103;
t84 = pkin(5) * t296 - pkin(11) * t329;
t67 = qJD(6) * t290 - t121 * t266 + t270 * t332;
t66 = qJD(6) * t144 + t121 * t270 + t266 * t332;
t58 = pkin(5) * t359 - t62;
t45 = Ifges(6,5) * t322 - t381 + t383;
t44 = -t103 * Ifges(6,2) + Ifges(6,6) * t322 + t382;
t35 = -mrSges(7,2) * t103 + mrSges(7,3) * t51;
t34 = mrSges(7,1) * t103 - mrSges(7,3) * t50;
t30 = qJD(5) * t81 - t271 * t111 + t112 * t267;
t29 = qJD(5) * t297 + t111 * t267 + t112 * t271;
t17 = t266 * t84 + t270 * t31;
t16 = -t266 * t31 + t270 * t84;
t14 = t50 * Ifges(7,4) + t51 * Ifges(7,2) + t103 * Ifges(7,6);
t10 = -qJD(6) * t57 + t132 * t270 - t266 * t29;
t9 = qJD(6) * t56 + t132 * t266 + t270 * t29;
t6 = -pkin(5) * t322 - t8;
t3 = [t10 * t73 + t111 * t179 + t112 * t178 + t29 * t128 + t133 * t216 + t142 * t199 + t143 * t198 + t225 * t209 + t56 * t34 + t57 * t35 + t9 * t72 + t81 * t88 - t395 * t297 - t365 * t30 + (-mrSges(3,1) * t269 - mrSges(3,2) * t273) * qJD(2) ^ 2 * t262 + t364 * t180 - t340 * t132 + (t217 * t262 * t269 + (t180 * t272 - t181 * t268) * qJD(3) * mrSges(4,3)) * t347 + m(4) * (t108 * t181 + t363 - t132 * t156 + t133 * t157 + (qJD(1) * t225 + t200) * t326) + m(6) * (t104 * t132 + t180 * t93 + t29 * t32 + t297 * t8 - t30 * t31 + t7 * t81) + m(5) * (t111 * t85 + t112 * t86 + t132 * t138 + t142 * t69 + t143 * t70 + t363) + m(7) * (t1 * t57 + t10 * t11 + t12 * t9 + t2 * t56 + t25 * t30 - t297 * t6); (t248 / 0.2e1 + t317) * t264 + m(4) * (t108 * t230 - t109 * t229 - t156 * t220 + t157 * t219) + (t1 * t144 - t11 * t66 + t12 * t67 + t2 * t290) * mrSges(7,3) + (-Ifges(6,4) * t456 + Ifges(7,5) * t415 - Ifges(6,2) * t455 - Ifges(6,6) * t457 + Ifges(7,6) * t417 + Ifges(7,3) * t413 - t458 + t460) * t122 - m(4) * (-t156 * t196 + t157 * t197 + t200 * t328) + t467 * t216 + t436 * t88 + t6 * (-mrSges(7,1) * t144 - mrSges(7,2) * t290) - t290 * t426 + t102 * (Ifges(6,1) * t172 + Ifges(6,4) * t295) / 0.2e1 + t93 * (-mrSges(6,1) * t295 + mrSges(6,2) * t172) - t103 * (Ifges(6,4) * t172 + Ifges(6,2) * t295) / 0.2e1 + (-Ifges(7,5) * t290 + Ifges(7,6) * t144 - Ifges(7,3) * t295) * t419 + (-Ifges(7,4) * t290 + Ifges(7,2) * t144 - Ifges(7,6) * t295) * t422 + (-Ifges(7,1) * t290 + Ifges(7,4) * t144 - Ifges(7,5) * t295) * t423 - t295 * t427 + (Ifges(7,5) * t66 + Ifges(7,6) * t67) * t413 + (-t172 * t8 + t295 * t7) * mrSges(6,3) + (Ifges(7,4) * t66 + Ifges(7,2) * t67) * t417 + (t224 * t70 - t226 * t69) * mrSges(5,3) - t295 * t461 + t109 * (-mrSges(5,1) * t224 + mrSges(5,2) * t226) + t213 * t191 + t155 * t198 + t154 * t199 + (t104 * t452 + t173 * t93 + t439 * t31 + t438 * t32 + t436 * t7 + t62 * t8) * m(6) + t184 * t83 + t175 * t409 + t174 * t410 + t45 * t411 + t44 * t412 + t295 * t462 + t173 * t49 + t144 * t14 / 0.2e1 + t62 * t87 + t25 * (-mrSges(7,1) * t67 + mrSges(7,2) * t66) + t67 * t37 / 0.2e1 + t66 * t38 / 0.2e1 + t58 * t18 + t28 * t35 + t27 * t34 + (-t217 * t337 + t109 * mrSges(4,3) * t268 - pkin(2) * t209 + (t108 * mrSges(4,3) + t70 * mrSges(5,2) - t69 * mrSges(5,1) - t247 / 0.2e1 + t429) * t272) * t261 + t434 * t178 + t435 * t179 + (t109 * t213 + t154 * t69 + t155 * t70 + t434 * t86 + t435 * t85 + (-t196 + t220) * t138) * m(5) + (Ifges(7,1) * t66 + Ifges(7,4) * t67) * t415 + t438 * t128 + t439 * t129 + t440 * t65 + t441 * t72 + t442 * t73 + (t1 * t28 + t11 * t442 + t441 * t12 + t2 * t27 + t440 * t25 + t58 * t6) * m(7) + t340 * t196 - t350 * t220 + (0.2e1 * t420 - t459) * t121 + ((t205 * t310 / 0.2e1 + t200 * mrSges(4,2) + t368 / 0.2e1 + t206 * t313 / 0.2e1 + t138 * (mrSges(5,1) * t260 + t393) + t137 * t407 + t136 * t408 - t156 * mrSges(4,3) + t190 / 0.2e1 + (-t260 * t86 - t263 * t85) * mrSges(5,3)) * t272 + (t373 / 0.2e1 + t200 * mrSges(4,1) - t86 * mrSges(5,2) - t367 / 0.2e1 + t85 * mrSges(5,1) + t372 / 0.2e1 + t374 / 0.2e1 + t376 / 0.2e1 - t397 + t398 + t369 / 0.2e1 + t77 / 0.2e1 + t135 / 0.2e1 - t189 / 0.2e1 - t157 * mrSges(4,3)) * t268) * t345 + ((-m(4) * pkin(2) + t394) * t328 + (-Ifges(4,6) * t264 + Ifges(6,5) * t411 + Ifges(6,6) * t412 + Ifges(5,5) * t409 + Ifges(5,6) * t410 - t230 * mrSges(4,3) - 0.3e1 / 0.2e1 * Ifges(4,4) * t360) * t344 + (-mrSges(4,1) * t328 + (Ifges(4,5) * t264 / 0.2e1 - t229 * mrSges(4,3) + (Ifges(5,4) * t226 + Ifges(5,2) * t224) * t408 + (Ifges(5,1) * t226 + Ifges(5,4) * t224) * t407 + (-0.3e1 / 0.2e1 * t263 * Ifges(5,5) + 0.3e1 / 0.2e1 * t260 * Ifges(5,6) + 0.3e1 / 0.2e1 * Ifges(4,4)) * t359 + (-Ifges(6,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,3) + 0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t360) * qJD(3)) * t272) * t347; t248 + ((t266 * t227 - t168) * mrSges(7,3) + t352 * mrSges(7,2)) * t12 + ((t270 * t227 + t169) * mrSges(7,3) - t352 * mrSges(7,1)) * t11 + (t109 * mrSges(5,2) + t175 / 0.2e1 - qJ(4) * t199 - qJD(4) * t179 - t69 * mrSges(5,3)) * t260 - t395 * t294 - t249 * (-Ifges(6,5) * t188 - Ifges(6,6) * t187) / 0.2e1 - t329 * (-Ifges(6,4) * t188 - Ifges(6,2) * t187) / 0.2e1 - t296 * (-Ifges(6,1) * t188 - Ifges(6,4) * t187) / 0.2e1 + t188 * t420 + t259 * t49 - t156 * t216 - pkin(3) * t191 - t187 * t36 / 0.2e1 + t187 * t78 / 0.2e1 - t115 * t179 + t186 * t88 + (Ifges(7,1) * t169 + Ifges(7,4) * t168 + Ifges(7,5) * t187) * t416 + (Ifges(7,4) * t169 + Ifges(7,2) * t168 + Ifges(7,6) * t187) * t418 + t168 * t424 + (Ifges(7,5) * t169 + Ifges(7,6) * t168 + Ifges(7,3) * t187) * t414 - t116 * t178 - t25 * (-mrSges(7,1) * t168 + mrSges(7,2) * t169) - t169 * t38 / 0.2e1 - t134 * t83 + t118 * t34 + t119 * t35 + (-t109 * mrSges(5,1) + t174 / 0.2e1 + qJ(4) * t198 + qJD(4) * t178 + t70 * mrSges(5,3)) * t263 + t317 + (t268 * t397 - t268 * t398 - t137 * t358 / 0.2e1 + t136 * t361 / 0.2e1 + (t272 * (Ifges(5,5) * t358 - Ifges(5,6) * t361 + Ifges(5,3) * t268) / 0.2e1 + t366 * t404) * t347 - t205 * (Ifges(5,4) * t358 - Ifges(5,2) * t361 + t384) / 0.2e1 - t200 * t316 - t86 * t291 - t254 * (Ifges(4,5) * t272 - Ifges(4,6) * t268) / 0.2e1 - t85 * t292 - t206 * (Ifges(5,1) * t358 - Ifges(5,4) * t361 + t385) / 0.2e1 - t138 * (mrSges(5,1) * t361 + mrSges(5,2) * t358) + t189 * t404 + (t156 * t272 + t157 * t268) * mrSges(4,3) + ((Ifges(5,5) * t260 / 0.2e1 + Ifges(5,6) * t407 + Ifges(6,5) * t234 / 0.2e1 - Ifges(6,6) * t233 / 0.2e1 - Ifges(4,6)) * t268 + ((Ifges(5,2) * t263 + t391) * t408 + (Ifges(5,1) * t260 + t390) * t407) * t272) * qJD(3) - (t190 + t251) * t272 / 0.2e1 - ((Ifges(4,1) * t272 - t392) * t347 + 0.2e1 * t77 + t135) * t268 / 0.2e1) * t347 - (t451 + t430) * t227 + t443 * t65 + t444 * t129 + t445 * t128 + (-t104 * t134 + t186 * t7 + t259 * t93 + t294 * t8 + t31 * t444 + t32 * t445) * m(6) + t446 * t73 + t447 * t72 + (t1 * t119 + t11 * t446 + t118 * t2 + t12 * t447 + t25 * t443 - t294 * t6) * m(7) + t350 * t157 + (t31 * t351 + t32 * t352) * mrSges(6,3) + (-mrSges(6,1) * t352 - mrSges(6,2) * t351) * t104 + (-t377 / 0.2e1 + t277) * t228 + (-pkin(3) * t109 + (-t260 * t85 + t263 * t86) * qJD(4) + (-t260 * t69 + t263 * t70) * qJ(4) - t115 * t85 - t116 * t86 - t138 * t157) * m(5) + (t6 * t314 + t307 * t419 + t312 * t423 + t309 * t422 + t93 * mrSges(6,2) + t45 / 0.2e1 + t383 / 0.2e1 - t381 / 0.2e1 + t15 * t403 + t14 * t406 - t8 * mrSges(6,3) + t318 * mrSges(7,3) + (mrSges(7,3) * t302 + t25 * t315 + t270 * t424 + t306 * t414 + t308 * t418 + t311 * t416 + t38 * t406) * qJD(6)) * t234 + (t48 / 0.2e1 + t47 / 0.2e1 + t93 * mrSges(6,1) + t427 - t44 / 0.2e1 - t382 / 0.2e1 - t7 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + t428) * t103 + t319) * t233; -t205 * t178 + t206 * t179 + t266 * t35 + t270 * t34 + t365 * t296 + t301 * qJD(6) + (-t128 - t301) * t329 + t364 + (-t147 * t302 - t296 * t25 - t318) * m(7) + (t296 * t31 - t32 * t329 + t93) * m(6) + (-t205 * t86 + t206 * t85 + t109) * m(5); t247 + (-pkin(5) * t6 - t11 * t16 - t12 * t17 - t25 * t32) * m(7) - t6 * t315 + ((-Ifges(6,1) / 0.2e1 + t428) * t296 + t437 - t430 + t459) * t329 + (-t277 + t458) * t296 - t429 + t306 * t419 + t308 * t422 + t311 * t423 + t266 * t426 + t14 * t403 - t31 * t128 - t17 * t72 - t16 * t73 - pkin(5) * t18 + ((-m(7) * t303 - t266 * t72 - t270 * t73) * qJD(6) + m(7) * t431 - t266 * t34 + t270 * t35) * pkin(11) + (t276 - t437) * qJD(6) + t365 * t32 + t431 * mrSges(7,3); -t25 * (mrSges(7,1) * t125 + mrSges(7,2) * t124) + (Ifges(7,1) * t124 - t388) * t416 + t37 * t415 + (Ifges(7,5) * t124 - Ifges(7,6) * t125) * t414 - t11 * t72 + t12 * t73 + (t11 * t124 + t12 * t125) * mrSges(7,3) + t319 + t13 + (-Ifges(7,2) * t125 + t120 + t38) * t418;];
tauc  = t3(:);
