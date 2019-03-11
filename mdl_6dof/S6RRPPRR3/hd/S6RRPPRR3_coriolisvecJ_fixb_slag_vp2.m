% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:54:45
% EndTime: 2019-03-09 08:55:23
% DurationCPUTime: 20.21s
% Computational Cost: add. (18917->744), mult. (57235->1047), div. (0->0), fcn. (46953->12), ass. (0->354)
t276 = sin(pkin(11));
t277 = sin(pkin(6));
t285 = cos(qJ(2));
t360 = cos(pkin(11));
t326 = t360 * t285;
t316 = t277 * t326;
t282 = sin(qJ(2));
t345 = qJD(1) * t277;
t336 = t282 * t345;
t234 = -qJD(1) * t316 + t276 * t336;
t275 = sin(pkin(12));
t278 = cos(pkin(12));
t281 = sin(qJ(5));
t284 = cos(qJ(5));
t257 = t275 * t284 + t278 * t281;
t173 = t257 * t234;
t251 = t257 * qJD(5);
t445 = t173 + t251;
t256 = t275 * t281 - t284 * t278;
t174 = t256 * t234;
t250 = t256 * qJD(5);
t470 = t174 + t250;
t395 = pkin(2) * t276;
t272 = qJ(4) + t395;
t389 = pkin(9) + t272;
t252 = t389 * t275;
t253 = t389 * t278;
t295 = -t284 * t252 - t253 * t281;
t279 = cos(pkin(6));
t396 = pkin(1) * t279;
t270 = t282 * t396;
t349 = t277 * t285;
t390 = pkin(8) + qJ(3);
t439 = t349 * t390 + t270;
t220 = t439 * qJD(1);
t208 = t276 * t220;
t271 = t285 * t396;
t266 = qJD(1) * t271;
t331 = t390 * t282;
t320 = t277 * t331;
t219 = -qJD(1) * t320 + t266;
t165 = t219 * t360 - t208;
t241 = (t276 * t285 + t282 * t360) * t277;
t237 = qJD(1) * t241;
t324 = pkin(2) * t336;
t177 = pkin(3) * t237 + qJ(4) * t234 + t324;
t102 = -t165 * t275 + t278 * t177;
t354 = t234 * t278;
t81 = pkin(4) * t237 + pkin(9) * t354 + t102;
t103 = t278 * t165 + t275 * t177;
t355 = t234 * t275;
t91 = pkin(9) * t355 + t103;
t451 = -qJD(4) * t256 + qJD(5) * t295 - t281 * t81 - t284 * t91;
t469 = -pkin(10) * t237 + t451;
t327 = t360 * t220;
t164 = t219 * t276 + t327;
t124 = -pkin(4) * t355 + t164;
t468 = pkin(5) * t445 + t470 * pkin(10) - t124;
t268 = qJD(1) * t279 + qJD(2);
t201 = pkin(2) * t268 + t219;
t151 = t201 * t360 - t208;
t145 = -t268 * pkin(3) + qJD(4) - t151;
t325 = -t237 * t275 + t278 * t268;
t106 = -pkin(4) * t325 + t145;
t198 = t237 * t278 + t268 * t275;
t152 = t276 * t201 + t327;
t146 = qJ(4) * t268 + t152;
t258 = (-pkin(2) * t285 - pkin(1)) * t277;
t249 = qJD(1) * t258 + qJD(3);
t159 = t234 * pkin(3) - t237 * qJ(4) + t249;
t92 = -t146 * t275 + t278 * t159;
t63 = pkin(4) * t234 - pkin(9) * t198 + t92;
t93 = t278 * t146 + t275 * t159;
t74 = pkin(9) * t325 + t93;
t31 = t281 * t63 + t284 * t74;
t232 = qJD(5) + t234;
t28 = pkin(10) * t232 + t31;
t280 = sin(qJ(6));
t283 = cos(qJ(6));
t462 = t198 * t284 + t281 * t325;
t463 = -t198 * t281 + t284 * t325;
t51 = -pkin(5) * t463 - pkin(10) * t462 + t106;
t17 = t28 * t283 + t280 * t51;
t460 = t17 * mrSges(7,2);
t16 = -t28 * t280 + t283 * t51;
t461 = t16 * mrSges(7,1);
t467 = -t106 * mrSges(6,1) + t31 * mrSges(6,3) + t460 - t461;
t30 = -t281 * t74 + t284 * t63;
t466 = t106 * mrSges(6,2) - t30 * mrSges(6,3);
t136 = qJD(6) - t463;
t192 = -t252 * t281 + t253 * t284;
t452 = -qJD(4) * t257 - qJD(5) * t192 + t281 * t91 - t284 * t81;
t404 = t232 / 0.2e1;
t410 = t462 / 0.2e1;
t411 = t463 / 0.2e1;
t412 = t136 / 0.2e1;
t110 = t232 * t280 + t283 * t462;
t414 = t110 / 0.2e1;
t109 = t232 * t283 - t280 * t462;
t416 = t109 / 0.2e1;
t377 = t136 * Ifges(7,3);
t378 = t110 * Ifges(7,5);
t379 = t109 * Ifges(7,6);
t43 = t377 + t378 + t379;
t369 = t232 * Ifges(6,6);
t384 = Ifges(6,4) * t462;
t72 = Ifges(6,2) * t463 + t369 + t384;
t465 = -Ifges(6,4) * t410 + Ifges(7,5) * t414 - Ifges(6,2) * t411 - Ifges(6,6) * t404 + Ifges(7,6) * t416 + Ifges(7,3) * t412 - t72 / 0.2e1 + t43 / 0.2e1 - t467;
t135 = Ifges(6,4) * t463;
t370 = t232 * Ifges(6,5);
t73 = Ifges(6,1) * t462 + t135 + t370;
t464 = Ifges(6,1) * t410 + Ifges(6,4) * t411 + Ifges(6,5) * t404 + t73 / 0.2e1 + t466;
t235 = qJD(2) * t241;
t229 = qJD(1) * t235;
t343 = qJD(2) * t277;
t236 = (-t276 * t282 + t326) * t343;
t230 = qJD(1) * t236;
t98 = qJD(5) * t463 - t230 * t256;
t49 = qJD(6) * t109 + t229 * t280 + t283 * t98;
t426 = t49 / 0.2e1;
t50 = -qJD(6) * t110 + t229 * t283 - t280 * t98;
t425 = t50 / 0.2e1;
t99 = qJD(5) * t462 + t230 * t257;
t418 = t99 / 0.2e1;
t459 = t92 * mrSges(5,1);
t458 = t93 * mrSges(5,2);
t20 = -mrSges(7,1) * t50 + mrSges(7,2) * t49;
t82 = mrSges(6,1) * t229 - mrSges(6,3) * t98;
t457 = t20 - t82;
t337 = t360 * pkin(2);
t274 = -t337 - pkin(3);
t260 = -t278 * pkin(4) + t274;
t188 = t256 * pkin(5) - t257 * pkin(10) + t260;
t132 = t188 * t283 - t192 * t280;
t456 = qJD(6) * t132 + t468 * t280 + t283 * t469;
t133 = t188 * t280 + t192 * t283;
t455 = -qJD(6) * t133 - t280 * t469 + t468 * t283;
t112 = mrSges(6,1) * t232 - mrSges(6,3) * t462;
t57 = -mrSges(7,1) * t109 + mrSges(7,2) * t110;
t361 = t112 - t57;
t453 = pkin(5) * t237 - t452;
t450 = t145 * (mrSges(5,1) * t275 + mrSges(5,2) * t278);
t128 = -t174 * t280 + t237 * t283;
t340 = qJD(6) * t257;
t293 = -t280 * t250 + t283 * t340;
t447 = t128 + t293;
t129 = t174 * t283 + t237 * t280;
t292 = t283 * t250 + t280 * t340;
t446 = t129 + t292;
t348 = -mrSges(4,1) * t268 - mrSges(5,1) * t325 + mrSges(5,2) * t198 + mrSges(4,3) * t237;
t397 = t283 / 0.2e1;
t383 = Ifges(7,4) * t110;
t44 = Ifges(7,2) * t109 + Ifges(7,6) * t136 + t383;
t427 = -t44 / 0.2e1;
t108 = Ifges(7,4) * t109;
t45 = Ifges(7,1) * t110 + Ifges(7,5) * t136 + t108;
t444 = t280 * t427 + t45 * t397;
t160 = -mrSges(5,2) * t234 + mrSges(5,3) * t325;
t161 = mrSges(5,1) * t234 - mrSges(5,3) * t198;
t443 = t160 * t278 - t161 * t275;
t263 = qJD(2) * t266;
t288 = (-qJD(2) * t331 + qJD(3) * t285) * t277;
t189 = qJD(1) * t288 + t263;
t350 = t277 * t282;
t203 = -qJD(2) * t439 - qJD(3) * t350;
t287 = qJD(1) * t203;
t126 = t360 * t189 + t276 * t287;
t118 = qJD(4) * t268 + t126;
t330 = qJD(1) * t343;
t319 = t282 * t330;
t297 = pkin(2) * t319;
t127 = pkin(3) * t229 - qJ(4) * t230 - qJD(4) * t237 + t297;
t69 = -t118 * t275 + t278 * t127;
t70 = t278 * t118 + t275 * t127;
t442 = -t275 * t69 + t278 * t70;
t125 = t189 * t276 - t360 * t287;
t357 = t230 * t275;
t107 = pkin(4) * t357 + t125;
t29 = pkin(5) * t99 - pkin(10) * t98 + t107;
t341 = qJD(5) * t284;
t342 = qJD(5) * t281;
t356 = t230 * t278;
t55 = pkin(4) * t229 - pkin(9) * t356 + t69;
t59 = -pkin(9) * t357 + t70;
t7 = t281 * t55 + t284 * t59 + t63 * t341 - t342 * t74;
t5 = pkin(10) * t229 + t7;
t1 = qJD(6) * t16 + t280 * t29 + t283 * t5;
t2 = -qJD(6) * t17 - t280 * t5 + t283 * t29;
t315 = t1 * t283 - t2 * t280;
t8 = -qJD(5) * t31 - t281 * t59 + t284 * t55;
t441 = t8 * mrSges(6,1) - t7 * mrSges(6,2) + Ifges(6,5) * t98 - Ifges(6,6) * t99;
t27 = -pkin(5) * t232 - t30;
t306 = t16 * t283 + t17 * t280;
t307 = Ifges(7,5) * t283 - Ifges(7,6) * t280;
t381 = Ifges(7,4) * t283;
t309 = -Ifges(7,2) * t280 + t381;
t382 = Ifges(7,4) * t280;
t311 = Ifges(7,1) * t283 - t382;
t313 = mrSges(7,1) * t280 + mrSges(7,2) * t283;
t440 = -t306 * mrSges(7,3) + t27 * t313 + t307 * t412 + t309 * t416 + t311 * t414 + t444;
t216 = pkin(2) * t279 + t271 - t320;
t346 = pkin(8) * t349 + t270;
t233 = qJ(3) * t349 + t346;
t172 = t276 * t216 + t360 * t233;
t162 = qJ(4) * t279 + t172;
t240 = t276 * t350 - t316;
t182 = t240 * pkin(3) - t241 * qJ(4) + t258;
t104 = -t162 * t275 + t278 * t182;
t212 = t241 * t278 + t275 * t279;
t80 = pkin(4) * t240 - pkin(9) * t212 + t104;
t105 = t278 * t162 + t275 * t182;
t211 = -t241 * t275 + t278 * t279;
t89 = pkin(9) * t211 + t105;
t388 = t281 * t80 + t284 * t89;
t267 = qJD(2) * t271;
t202 = t267 + t288;
t150 = t360 * t202 + t276 * t203;
t138 = qJD(4) * t279 + t150;
t334 = t282 * t343;
t323 = pkin(2) * t334;
t148 = pkin(3) * t235 - qJ(4) * t236 - qJD(4) * t241 + t323;
t84 = -t138 * t275 + t278 * t148;
t62 = -pkin(9) * t236 * t278 + pkin(4) * t235 + t84;
t353 = t236 * t275;
t85 = t278 * t138 + t275 * t148;
t68 = -pkin(9) * t353 + t85;
t12 = -qJD(5) * t388 - t281 * t68 + t284 * t62;
t438 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t419 = -t99 / 0.2e1;
t420 = t98 / 0.2e1;
t436 = t107 * mrSges(6,2) - t8 * mrSges(6,3) + Ifges(6,1) * t420 + Ifges(6,4) * t419;
t13 = Ifges(7,5) * t49 + Ifges(7,6) * t50 + Ifges(7,3) * t99;
t434 = t107 * mrSges(6,1) - t7 * mrSges(6,3) - Ifges(6,4) * t420 + Ifges(7,5) * t426 - Ifges(6,2) * t419 + Ifges(7,6) * t425 + Ifges(7,3) * t418 + t13 / 0.2e1 + t438;
t432 = -0.2e1 * pkin(1);
t430 = Ifges(7,1) * t426 + Ifges(7,4) * t425 + Ifges(7,5) * t418;
t429 = -t43 / 0.2e1;
t423 = t72 / 0.2e1;
t422 = -t73 / 0.2e1;
t417 = -t109 / 0.2e1;
t415 = -t110 / 0.2e1;
t413 = -t136 / 0.2e1;
t296 = t284 * t211 - t212 * t281;
t409 = t296 / 0.2e1;
t154 = t211 * t281 + t212 * t284;
t408 = t154 / 0.2e1;
t407 = t211 / 0.2e1;
t406 = t212 / 0.2e1;
t405 = t229 / 0.2e1;
t403 = t234 / 0.2e1;
t401 = -t275 / 0.2e1;
t400 = t278 / 0.2e1;
t399 = t279 / 0.2e1;
t398 = -t282 / 0.2e1;
t387 = Ifges(3,4) * t282;
t231 = Ifges(4,4) * t234;
t386 = Ifges(5,4) * t275;
t385 = Ifges(5,4) * t278;
t380 = Ifges(3,5) * t285;
t376 = t463 * Ifges(6,6);
t375 = t462 * Ifges(6,5);
t374 = t151 * mrSges(4,3);
t373 = t152 * mrSges(4,3);
t372 = t325 * Ifges(5,6);
t371 = t198 * Ifges(5,5);
t368 = t232 * Ifges(6,3);
t367 = t237 * Ifges(4,1);
t366 = t237 * Ifges(4,4);
t365 = t268 * Ifges(4,5);
t364 = t268 * Ifges(4,6);
t352 = t257 * t280;
t351 = t257 * t283;
t168 = mrSges(5,1) * t357 + mrSges(5,2) * t356;
t335 = t285 * t345;
t333 = (t198 * Ifges(5,4) + Ifges(5,2) * t325 + Ifges(5,6) * t234) * t401;
t332 = (t198 * Ifges(5,1) + Ifges(5,4) * t325 + Ifges(5,5) * t234) * t400;
t46 = t99 * mrSges(6,1) + t98 * mrSges(6,2);
t149 = t202 * t276 - t360 * t203;
t322 = mrSges(3,3) * t336;
t321 = mrSges(3,3) * t335;
t113 = pkin(4) * t353 + t149;
t312 = t278 * Ifges(5,1) - t386;
t310 = -t275 * Ifges(5,2) + t385;
t308 = Ifges(5,5) * t278 - Ifges(5,6) * t275;
t305 = -t16 * t280 + t17 * t283;
t25 = mrSges(7,1) * t99 - mrSges(7,3) * t49;
t26 = -mrSges(7,2) * t99 + mrSges(7,3) * t50;
t304 = -t280 * t25 + t283 * t26;
t303 = t275 * t92 - t278 * t93;
t33 = pkin(10) * t240 + t388;
t171 = t216 * t360 - t276 * t233;
t163 = -t279 * pkin(3) - t171;
t121 = -t211 * pkin(4) + t163;
t56 = -pkin(5) * t296 - t154 * pkin(10) + t121;
t19 = t280 * t56 + t283 * t33;
t18 = -t280 * t33 + t283 * t56;
t65 = -mrSges(7,2) * t136 + mrSges(7,3) * t109;
t66 = mrSges(7,1) * t136 - mrSges(7,3) * t110;
t302 = -t280 * t66 + t283 * t65;
t301 = -t280 * t65 - t283 * t66;
t37 = -t281 * t89 + t284 * t80;
t120 = t154 * t283 + t240 * t280;
t119 = -t154 * t280 + t240 * t283;
t111 = -mrSges(6,2) * t232 + mrSges(6,3) * t463;
t294 = -t111 - t302;
t11 = t281 * t62 + t284 * t68 + t80 * t341 - t342 * t89;
t291 = t268 * (-Ifges(3,6) * t282 + t380);
t238 = -pkin(8) * t319 + t263;
t248 = t346 * qJD(2);
t239 = qJD(1) * t248;
t289 = -t239 * mrSges(3,1) - t238 * mrSges(3,2) - t126 * mrSges(4,2);
t265 = Ifges(3,4) * t335;
t262 = t330 * t380;
t254 = -pkin(8) * t350 + t271;
t247 = -pkin(8) * t334 + t267;
t246 = t346 * qJD(1);
t245 = -pkin(8) * t336 + t266;
t244 = -t268 * mrSges(3,2) + t321;
t243 = mrSges(3,1) * t268 - t322;
t226 = Ifges(4,5) * t230;
t225 = Ifges(4,6) * t229;
t224 = Ifges(6,3) * t229;
t223 = t230 * mrSges(4,2);
t218 = Ifges(3,1) * t336 + t268 * Ifges(3,5) + t265;
t217 = Ifges(3,6) * t268 + (Ifges(3,2) * t285 + t387) * t345;
t205 = -mrSges(4,2) * t268 - mrSges(4,3) * t234;
t183 = mrSges(4,1) * t234 + mrSges(4,2) * t237;
t179 = mrSges(5,1) * t229 - mrSges(5,3) * t356;
t178 = -mrSges(5,2) * t229 - mrSges(5,3) * t357;
t176 = -t231 + t365 + t367;
t175 = -t234 * Ifges(4,2) + t364 + t366;
t131 = t229 * Ifges(5,5) + t230 * t312;
t130 = t229 * Ifges(5,6) + t230 * t310;
t114 = t234 * Ifges(5,3) + t371 + t372;
t101 = qJD(5) * t154 + t236 * t257;
t100 = qJD(5) * t296 - t236 * t256;
t87 = pkin(5) * t462 - pkin(10) * t463;
t86 = -mrSges(6,1) * t463 + mrSges(6,2) * t462;
t83 = -mrSges(6,2) * t229 - mrSges(6,3) * t99;
t71 = t368 + t375 + t376;
t53 = -qJD(6) * t120 - t100 * t280 + t235 * t283;
t52 = qJD(6) * t119 + t100 * t283 + t235 * t280;
t42 = t98 * Ifges(6,1) - t99 * Ifges(6,4) + t229 * Ifges(6,5);
t41 = t98 * Ifges(6,4) - t99 * Ifges(6,2) + t229 * Ifges(6,6);
t34 = pkin(5) * t101 - pkin(10) * t100 + t113;
t32 = -pkin(5) * t240 - t37;
t22 = t280 * t87 + t283 * t30;
t21 = -t280 * t30 + t283 * t87;
t14 = t49 * Ifges(7,4) + t50 * Ifges(7,2) + t99 * Ifges(7,6);
t10 = -pkin(5) * t235 - t12;
t9 = pkin(10) * t235 + t11;
t6 = -pkin(5) * t229 - t8;
t4 = -qJD(6) * t19 - t280 * t9 + t283 * t34;
t3 = qJD(6) * t18 + t280 * t34 + t283 * t9;
t15 = [(Ifges(7,1) * t52 + Ifges(7,4) * t53) * t414 + (Ifges(7,1) * t120 + Ifges(7,4) * t119) * t426 + t37 * t82 + (Ifges(7,5) * t52 + Ifges(7,6) * t53) * t412 + (Ifges(7,5) * t120 + Ifges(7,6) * t119) * t418 + (Ifges(7,4) * t52 + Ifges(7,2) * t53) * t416 + (Ifges(7,4) * t120 + Ifges(7,2) * t119) * t425 + t464 * t100 + t465 * t101 + m(6) * (t106 * t113 + t107 * t121 + t11 * t31 + t12 * t30 + t37 * t8 + t388 * t7) + t388 * t83 + ((Ifges(3,5) * t399 - t254 * mrSges(3,3) + (mrSges(3,2) * t432 + 0.3e1 / 0.2e1 * Ifges(3,4) * t285) * t277) * t285 + (-t346 * mrSges(3,3) - Ifges(3,6) * t279 + (mrSges(3,1) * t432 - 0.3e1 / 0.2e1 * t387 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t285) * t277 + (m(4) * t258 + mrSges(4,1) * t240 + mrSges(4,2) * t241) * pkin(2)) * t282) * t330 - t434 * t296 + m(3) * (t238 * t346 - t239 * t254 - t245 * t248 + t246 * t247) + (t262 / 0.2e1 + t226 / 0.2e1 - t225 / 0.2e1 + t289) * t279 + m(4) * (-t125 * t171 + t126 * t172 - t149 * t151 + t150 * t152 + t249 * t323) + (t238 * t285 + t239 * t282) * t277 * mrSges(3,3) + (t1 * t119 - t120 * t2 - t16 * t52 + t17 * t53) * mrSges(7,3) + (-mrSges(4,1) * t279 - mrSges(5,1) * t211 + mrSges(5,2) * t212 + mrSges(4,3) * t241) * t125 + t348 * t149 + t258 * t223 + t436 * t154 + (t282 * pkin(2) * t183 + t285 * t218 / 0.2e1 + t217 * t398 + t291 / 0.2e1 + (-t245 * t285 - t246 * t282) * mrSges(3,3)) * t343 + ((Ifges(5,4) * t212 + Ifges(5,2) * t211) * t401 - t171 * mrSges(4,3) + Ifges(4,1) * t241 + Ifges(4,5) * t399 + (Ifges(5,1) * t212 + Ifges(5,4) * t211) * t400) * t230 + (-Ifges(4,4) * t241 - Ifges(4,6) * t279 / 0.2e1 + Ifges(6,5) * t408 + Ifges(6,6) * t409 + Ifges(5,5) * t406 + Ifges(5,6) * t407 - t172 * mrSges(4,3) + t258 * mrSges(4,1)) * t229 + m(5) * (t104 * t69 + t105 * t70 + t125 * t163 + t145 * t149 + t84 * t92 + t85 * t93) + m(7) * (t1 * t19 + t10 * t27 + t16 * t4 + t17 * t3 + t18 * t2 + t32 * t6) + t4 * t66 + t3 * t65 + t27 * (-mrSges(7,1) * t53 + mrSges(7,2) * t52) + t53 * t44 / 0.2e1 + t10 * t57 + t52 * t45 / 0.2e1 + t32 * t20 + t18 * t25 + t19 * t26 + (t211 * t70 - t212 * t69) * mrSges(5,3) + (-t70 * mrSges(5,2) + t69 * mrSges(5,1) + t224 / 0.2e1 - t126 * mrSges(4,3) + (-Ifges(4,4) + t308) * t230 + (Ifges(6,3) / 0.2e1 + Ifges(5,3) + Ifges(4,2)) * t229 + t441) * t240 + t11 * t111 + t12 * t112 + t113 * t86 + t119 * t14 / 0.2e1 + t6 * (-mrSges(7,1) * t119 + mrSges(7,2) * t120) + t121 * t46 + (t332 + t333 - t374 + t176 / 0.2e1 - t231 / 0.2e1 + t367 / 0.2e1 + t249 * mrSges(4,2) + t365 / 0.2e1 + t450 + t325 * t310 / 0.2e1 + t198 * t312 / 0.2e1 + t308 * t403 + (-t275 * t93 - t278 * t92) * mrSges(5,3)) * t236 + t131 * t406 + t130 * t407 + t42 * t408 + t41 * t409 + t120 * t430 + t85 * t160 + t84 * t161 + t163 * t168 + (-t373 + t71 / 0.2e1 + t30 * mrSges(6,1) - t31 * mrSges(6,2) + t114 / 0.2e1 + t376 / 0.2e1 + t375 / 0.2e1 + t368 / 0.2e1 - t175 / 0.2e1 - t366 / 0.2e1 + t249 * mrSges(4,1) - t364 / 0.2e1 - t458 + t459 + t372 / 0.2e1 + t371 / 0.2e1 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t234) * t235 + t105 * t178 + t104 * t179 + t150 * t205 + t247 * t244 - t248 * t243; (t278 * t178 - t275 * t179) * t272 + t234 * t450 + t237 * t458 + t217 * t336 / 0.2e1 - t225 + t226 + (-Ifges(7,5) * t292 - Ifges(7,6) * t293) * t412 - t234 * (Ifges(5,3) * t237 - t234 * t308) / 0.2e1 - (t444 + t464) * t250 + t465 * t251 - t234 * t374 - t31 * (-mrSges(6,2) * t237 + mrSges(6,3) * t173) - t232 * (Ifges(6,5) * t174 + Ifges(6,6) * t173 + Ifges(6,3) * t237) / 0.2e1 - (t280 * t45 + t283 * t44) * t340 / 0.2e1 - (-Ifges(4,1) * t234 + t114 - t366 + t71) * t237 / 0.2e1 - t325 * (Ifges(5,6) * t237 - t234 * t310) / 0.2e1 + (Ifges(7,5) * t129 + Ifges(7,6) * t128 - Ifges(7,3) * t173) * t413 + (Ifges(7,1) * t129 + Ifges(7,4) * t128 - Ifges(7,5) * t173) * t415 + (Ifges(7,4) * t129 + Ifges(7,2) * t128 - Ifges(7,6) * t173) * t417 - t106 * (-mrSges(6,1) * t173 + mrSges(6,2) * t174) + t262 - ((-Ifges(3,2) * t336 + t218 + t265) * t285 + t291) * t345 / 0.2e1 + (-t106 * t124 + t107 * t260 + t192 * t7 + t295 * t8 + t30 * t452 + t31 * t451) * m(6) + (t1 * t133 + t132 * t2 + t16 * t455 + t17 * t456 + t27 * t453 - t295 * t6) * m(7) - t457 * t295 - Ifges(3,6) * t319 - t463 * (Ifges(6,4) * t174 + Ifges(6,2) * t173 + Ifges(6,6) * t237) / 0.2e1 - t462 * (Ifges(6,1) * t174 + Ifges(6,4) * t173 + Ifges(6,5) * t237) / 0.2e1 + t234 * t332 + t234 * t333 - t14 * t352 / 0.2e1 - t198 * (Ifges(5,5) * t237 - t234 * t312) / 0.2e1 + (-t229 * t395 - t230 * t337) * mrSges(4,3) + (-Ifges(6,6) * t405 - t41 / 0.2e1 + t434) * t256 + (-Ifges(7,4) * t292 - Ifges(7,2) * t293) * t416 + t289 + (t6 * t313 + Ifges(6,5) * t405 + t307 * t418 + t309 * t425 + t311 * t426 + t42 / 0.2e1 + t436) * t257 + (pkin(1) * (mrSges(3,1) * t282 + mrSges(3,2) * t285) + (Ifges(3,1) * t285 - t387) * t398) * qJD(1) ^ 2 * t277 ^ 2 + (-Ifges(7,1) * t292 - Ifges(7,4) * t293) * t414 - t183 * t324 + (Ifges(5,1) * t275 + t385) * t356 / 0.2e1 - (Ifges(5,2) * t278 + t386) * t357 / 0.2e1 + ((-t125 * t360 + t126 * t276) * pkin(2) + t151 * t164 - t152 * t165 - t249 * t324) * m(4) + t237 * t373 + t173 * t461 + (-mrSges(5,1) * t278 + mrSges(5,2) * t275 - mrSges(4,1)) * t125 + (-t354 * t92 - t355 * t93 + t442) * mrSges(5,3) + (-t303 * qJD(4) - t102 * t92 - t103 * t93 + t125 * t274 - t145 * t164 + t272 * t442) * m(5) + t443 * qJD(4) - t348 * t164 + (mrSges(7,1) * t447 - mrSges(7,2) * t446) * t27 + (-t1 * t352 + t16 * t446 - t17 * t447 - t2 * t351) * mrSges(7,3) - t124 * t86 + t130 * t400 + (Ifges(5,5) * t275 + Ifges(5,6) * t278) * t405 + t174 * t422 - t173 * t423 + t128 * t427 - t173 * t429 + t351 * t430 - t129 * t45 / 0.2e1 + t451 * t111 + t132 * t25 + t452 * t112 + t133 * t26 + t453 * t57 + (-Ifges(4,2) * t237 + t176 - t231) * t403 - t103 * t160 - t102 * t161 + t455 * t66 + t456 * t65 - t237 * t459 - t173 * t460 + t192 * t83 - t165 * t205 - t30 * (mrSges(6,1) * t237 - mrSges(6,3) * t174) + t237 * t175 / 0.2e1 - t249 * (mrSges(4,1) * t237 - mrSges(4,2) * t234) + (t321 - t244) * t245 + t260 * t46 - t268 * (-Ifges(4,5) * t234 - Ifges(4,6) * t237) / 0.2e1 + t274 * t168 + t275 * t131 / 0.2e1 + (t322 + t243) * t246; t229 * mrSges(4,1) - t174 * t111 - t128 * t66 - t129 * t65 + t275 * t178 + t278 * t179 + t223 + t457 * t256 + t294 * t250 + (-t86 - t348) * t237 + (t205 + t443) * t234 + (qJD(6) * t301 + t304 + t83) * t257 - t361 * t445 + (t256 * t6 - t305 * t250 + (-qJD(6) * t306 + t315) * t257 - t128 * t16 - t129 * t17 + t445 * t27) * m(7) + (-t106 * t237 - t256 * t8 + t257 * t7 - t445 * t30 - t31 * t470) * m(6) + (-t145 * t237 - t234 * t303 + t275 * t70 + t278 * t69) * m(5) + (t151 * t237 + t152 * t234 + t297) * m(4); t302 * qJD(6) + t361 * t462 + t294 * t463 - t325 * t160 + t198 * t161 + t283 * t25 + t280 * t26 + t168 + t46 + (t1 * t280 + t136 * t305 + t2 * t283 - t462 * t27) * m(7) + (t30 * t462 - t31 * t463 + t107) * m(6) + (t198 * t92 - t325 * t93 + t125) * m(5); t315 * mrSges(7,3) + t361 * t31 + t224 + (t304 + (-m(7) * t306 + t301) * qJD(6) + m(7) * t315) * pkin(10) + (-pkin(5) * t6 - t16 * t21 - t17 * t22 - t27 * t31) * m(7) + t441 - t21 * t66 - t22 * t65 - pkin(5) * t20 + (t422 - t135 / 0.2e1 - t370 / 0.2e1 + (Ifges(6,2) / 0.2e1 - Ifges(6,1) / 0.2e1) * t462 - t440 - t466) * t463 + t440 * qJD(6) - t30 * t111 + t14 * t397 + (Ifges(7,5) * t280 + Ifges(7,6) * t283) * t418 + (Ifges(7,2) * t283 + t382) * t425 + (Ifges(7,1) * t280 + t381) * t426 + t280 * t430 + (t369 / 0.2e1 - t379 / 0.2e1 - t378 / 0.2e1 - t377 / 0.2e1 + t429 + t423 + t384 / 0.2e1 + t467) * t462 + t6 * (-mrSges(7,1) * t283 + mrSges(7,2) * t280); -t27 * (mrSges(7,1) * t110 + mrSges(7,2) * t109) + (Ifges(7,1) * t109 - t383) * t415 + t44 * t414 + (Ifges(7,5) * t109 - Ifges(7,6) * t110) * t413 - t16 * t65 + t17 * t66 + (t109 * t16 + t110 * t17) * mrSges(7,3) + t13 + (-Ifges(7,2) * t110 + t108 + t45) * t417 + t438;];
tauc  = t15(:);
