% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:15:52
% EndTime: 2019-03-09 10:16:30
% DurationCPUTime: 19.66s
% Computational Cost: add. (16855->727), mult. (43244->997), div. (0->0), fcn. (32783->10), ass. (0->338)
t286 = sin(qJ(2));
t289 = cos(qJ(2));
t353 = sin(pkin(10));
t314 = qJD(1) * t353;
t354 = cos(pkin(10));
t315 = qJD(1) * t354;
t247 = -t286 * t314 + t289 * t315;
t248 = -t286 * t315 - t289 * t314;
t338 = qJD(1) * t286;
t331 = pkin(2) * t338;
t182 = -pkin(3) * t248 - pkin(8) * t247 + t331;
t371 = -qJ(3) - pkin(7);
t272 = t371 * t289;
t266 = qJD(1) * t272;
t252 = t353 * t266;
t271 = t371 * t286;
t265 = qJD(1) * t271;
t203 = t265 * t354 + t252;
t285 = sin(qJ(4));
t288 = cos(qJ(4));
t130 = t288 * t182 - t203 * t285;
t325 = t353 * pkin(2);
t277 = t325 + pkin(8);
t341 = qJ(5) + t277;
t311 = qJD(4) * t341;
t346 = t247 * t288;
t474 = pkin(4) * t248 + qJ(5) * t346 - qJD(5) * t285 - t288 * t311 - t130;
t131 = t285 * t182 + t288 * t203;
t347 = t247 * t285;
t473 = -qJ(5) * t347 - qJD(5) * t288 + t285 * t311 + t131;
t282 = sin(pkin(11));
t283 = cos(pkin(11));
t262 = t282 * t288 + t283 * t285;
t172 = t262 * t247;
t246 = t262 * qJD(4);
t436 = t246 - t172;
t297 = t282 * t285 - t283 * t288;
t173 = t297 * t247;
t250 = t297 * qJD(4);
t435 = t250 - t173;
t284 = sin(qJ(6));
t287 = cos(qJ(6));
t240 = qJD(4) - t247;
t213 = qJD(2) * t288 + t248 * t285;
t214 = qJD(2) * t285 - t248 * t288;
t150 = t213 * t282 + t214 * t283;
t460 = pkin(9) * t150;
t280 = -pkin(2) * t289 - pkin(1);
t339 = qJD(1) * t280;
t267 = qJD(3) + t339;
t171 = -t247 * pkin(3) + t248 * pkin(8) + t267;
t257 = qJD(2) * pkin(2) + t265;
t317 = t354 * t266;
t199 = t353 * t257 - t317;
t192 = qJD(2) * pkin(8) + t199;
t120 = t288 * t171 - t192 * t285;
t101 = -qJ(5) * t214 + t120;
t82 = pkin(4) * t240 + t101;
t121 = t171 * t285 + t192 * t288;
t102 = qJ(5) * t213 + t121;
t95 = t282 * t102;
t48 = t283 * t82 - t95;
t32 = pkin(5) * t240 - t460 + t48;
t312 = t283 * t213 - t214 * t282;
t454 = pkin(9) * t312;
t344 = t283 * t102;
t49 = t282 * t82 + t344;
t33 = t49 + t454;
t11 = -t284 * t33 + t287 * t32;
t12 = t284 * t32 + t287 * t33;
t263 = t286 * t354 + t289 * t353;
t249 = t263 * qJD(2);
t232 = qJD(1) * t249;
t261 = t286 * t353 - t289 * t354;
t251 = t261 * qJD(2);
t233 = qJD(1) * t251;
t165 = qJD(4) * t213 - t233 * t288;
t166 = -qJD(4) * t214 + t233 * t285;
t108 = -t165 * t282 + t166 * t283;
t109 = t165 * t283 + t166 * t282;
t457 = -t150 * t284 + t287 * t312;
t30 = qJD(6) * t457 + t108 * t284 + t109 * t287;
t94 = t150 * t287 + t284 * t312;
t31 = -qJD(6) * t94 + t108 * t287 - t109 * t284;
t332 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t232;
t382 = Ifges(7,4) * t94;
t231 = qJD(6) + t240;
t392 = -t231 / 0.2e1;
t407 = -t94 / 0.2e1;
t318 = qJD(2) * t371;
t244 = qJD(3) * t289 + t286 * t318;
t221 = t244 * qJD(1);
t245 = -t286 * qJD(3) + t289 * t318;
t290 = qJD(1) * t245;
t164 = t221 * t354 + t290 * t353;
t333 = qJD(1) * qJD(2);
t321 = t286 * t333;
t310 = pkin(2) * t321;
t168 = pkin(3) * t232 + pkin(8) * t233 + t310;
t69 = -qJD(4) * t121 - t164 * t285 + t288 * t168;
t38 = pkin(4) * t232 - qJ(5) * t165 - qJD(5) * t214 + t69;
t334 = qJD(4) * t288;
t335 = qJD(4) * t285;
t68 = t288 * t164 + t285 * t168 + t171 * t334 - t192 * t335;
t41 = qJ(5) * t166 + qJD(5) * t213 + t68;
t16 = t282 * t38 + t283 * t41;
t10 = pkin(9) * t108 + t16;
t15 = -t282 * t41 + t283 * t38;
t8 = pkin(5) * t232 - pkin(9) * t109 + t15;
t2 = qJD(6) * t11 + t10 * t287 + t284 * t8;
t3 = -qJD(6) * t12 - t10 * t284 + t287 * t8;
t453 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t198 = t257 * t354 + t252;
t191 = -qJD(2) * pkin(3) - t198;
t146 = -t213 * pkin(4) + qJD(5) + t191;
t85 = -pkin(5) * t312 + t146;
t472 = t332 + t453 + (Ifges(7,5) * t457 - Ifges(7,6) * t94) * t392 + (t11 * t457 + t12 * t94) * mrSges(7,3) - t85 * (mrSges(7,1) * t94 + mrSges(7,2) * t457) + (Ifges(7,1) * t457 - t382) * t407;
t448 = t473 * t282 + t283 * t474;
t447 = t282 * t474 - t473 * t283;
t471 = t335 - t347;
t88 = Ifges(7,4) * t457;
t470 = -Ifges(7,2) * t94 + t88;
t384 = -t285 / 0.2e1;
t468 = -Ifges(6,3) - Ifges(5,3);
t387 = -t247 / 0.2e1;
t464 = Ifges(4,2) * t387;
t462 = pkin(5) * t248 + pkin(9) * t435 + t448;
t461 = -pkin(9) * t436 + t447;
t459 = t150 * Ifges(6,4);
t202 = t265 * t353 - t317;
t439 = pkin(4) * t471 - t202;
t362 = t214 * Ifges(5,4);
t138 = t213 * Ifges(5,2) + t240 * Ifges(5,6) + t362;
t308 = mrSges(5,1) * t285 + mrSges(5,2) * t288;
t458 = t138 * t384 + t191 * t308;
t456 = -Ifges(4,6) / 0.2e1;
t421 = t30 / 0.2e1;
t420 = t31 / 0.2e1;
t404 = t108 / 0.2e1;
t403 = t109 / 0.2e1;
t390 = t232 / 0.2e1;
t336 = qJD(2) * t286;
t455 = pkin(2) * t336;
t258 = t341 * t285;
t259 = t341 * t288;
t193 = -t283 * t258 - t259 * t282;
t156 = -pkin(9) * t262 + t193;
t194 = -t282 * t258 + t283 * t259;
t157 = -pkin(9) * t297 + t194;
t103 = t156 * t287 - t157 * t284;
t452 = qJD(6) * t103 + t284 * t462 + t287 * t461;
t104 = t156 * t284 + t157 * t287;
t451 = -qJD(6) * t104 - t284 * t461 + t287 * t462;
t450 = Ifges(6,4) * t312;
t449 = t267 * mrSges(4,2);
t278 = pkin(4) * t283 + pkin(5);
t380 = pkin(4) * t282;
t242 = t278 * t287 - t284 * t380;
t54 = -t101 * t282 - t344;
t34 = t54 - t454;
t55 = t283 * t101 - t95;
t35 = t55 - t460;
t446 = t242 * qJD(6) - t284 * t34 - t287 * t35;
t243 = t278 * t284 + t287 * t380;
t445 = -t243 * qJD(6) + t284 * t35 - t287 * t34;
t444 = Ifges(4,5) * qJD(2);
t119 = -t172 * t284 - t173 * t287;
t200 = -t262 * t284 - t287 * t297;
t135 = qJD(6) * t200 - t246 * t284 - t250 * t287;
t441 = t135 - t119;
t118 = -t172 * t287 + t173 * t284;
t201 = t262 * t287 - t284 * t297;
t136 = -qJD(6) * t201 - t246 * t287 + t250 * t284;
t440 = t136 - t118;
t438 = Ifges(5,5) * t165 + Ifges(5,6) * t166;
t197 = t261 * pkin(3) - t263 * pkin(8) + t280;
t207 = t271 * t353 - t272 * t354;
t204 = t288 * t207;
t145 = t285 * t197 + t204;
t437 = pkin(5) * t436 + t439;
t369 = mrSges(4,3) * t248;
t340 = -qJD(2) * mrSges(4,1) - mrSges(5,1) * t213 + mrSges(5,2) * t214 - t369;
t434 = Ifges(6,5) * t109 + Ifges(6,6) * t108 - t232 * t468 + t438;
t433 = -t285 * t69 + t288 * t68;
t337 = qJD(1) * t289;
t351 = Ifges(3,6) * qJD(2);
t368 = Ifges(3,4) * t286;
t432 = t351 / 0.2e1 + (t289 * Ifges(3,2) + t368) * qJD(1) / 0.2e1 + pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t337);
t431 = qJD(2) * t456;
t357 = t248 * Ifges(4,4);
t430 = t357 / 0.2e1 + t464;
t429 = t69 * mrSges(5,1) + t15 * mrSges(6,1) - t68 * mrSges(5,2) - t16 * mrSges(6,2);
t303 = Ifges(5,5) * t288 - Ifges(5,6) * t285;
t366 = Ifges(5,4) * t288;
t305 = -Ifges(5,2) * t285 + t366;
t367 = Ifges(5,4) * t285;
t307 = Ifges(5,1) * t288 - t367;
t208 = Ifges(5,4) * t213;
t139 = t214 * Ifges(5,1) + t240 * Ifges(5,5) + t208;
t342 = t288 * t139;
t388 = t240 / 0.2e1;
t393 = t214 / 0.2e1;
t428 = t303 * t388 + t307 * t393 + t213 * t305 / 0.2e1 + t342 / 0.2e1 + t458;
t427 = t267 * mrSges(4,1) + t120 * mrSges(5,1) + t48 * mrSges(6,1) + t11 * mrSges(7,1) - t121 * mrSges(5,2) - t49 * mrSges(6,2) - t12 * mrSges(7,2) + t430 + t431;
t426 = -0.2e1 * pkin(1);
t424 = Ifges(7,4) * t421 + Ifges(7,2) * t420 + Ifges(7,6) * t390;
t423 = Ifges(7,1) * t421 + Ifges(7,4) * t420 + Ifges(7,5) * t390;
t422 = Ifges(5,3) / 0.2e1;
t46 = Ifges(7,2) * t457 + Ifges(7,6) * t231 + t382;
t419 = -t46 / 0.2e1;
t418 = t46 / 0.2e1;
t47 = Ifges(7,1) * t94 + Ifges(7,5) * t231 + t88;
t417 = -t47 / 0.2e1;
t416 = t47 / 0.2e1;
t415 = Ifges(6,4) * t403 + Ifges(6,2) * t404 + Ifges(6,6) * t390;
t414 = Ifges(6,1) * t403 + Ifges(6,4) * t404 + Ifges(6,5) * t390;
t80 = Ifges(6,2) * t312 + t240 * Ifges(6,6) + t459;
t413 = -t80 / 0.2e1;
t412 = t80 / 0.2e1;
t81 = t150 * Ifges(6,1) + t240 * Ifges(6,5) + t450;
t411 = -t81 / 0.2e1;
t410 = t81 / 0.2e1;
t409 = -t457 / 0.2e1;
t408 = t457 / 0.2e1;
t406 = t94 / 0.2e1;
t402 = -t312 / 0.2e1;
t401 = t312 / 0.2e1;
t400 = -t150 / 0.2e1;
t399 = t150 / 0.2e1;
t398 = t165 / 0.2e1;
t397 = t166 / 0.2e1;
t395 = -t213 / 0.2e1;
t394 = -t214 / 0.2e1;
t391 = t231 / 0.2e1;
t389 = -t240 / 0.2e1;
t386 = t248 / 0.2e1;
t383 = t288 / 0.2e1;
t381 = pkin(4) * t214;
t379 = pkin(7) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t338);
t373 = t457 * Ifges(7,6);
t372 = t94 * Ifges(7,5);
t296 = qJ(5) * t251 - qJD(5) * t263;
t181 = t244 * t354 + t245 * t353;
t183 = pkin(3) * t249 + pkin(8) * t251 + t455;
t313 = -t181 * t285 + t288 * t183;
t350 = qJ(5) * t263;
t61 = pkin(4) * t249 + t296 * t288 + (-t204 + (-t197 + t350) * t285) * qJD(4) + t313;
t324 = t263 * t334;
t327 = t288 * t181 + t285 * t183 + t197 * t334;
t65 = -qJ(5) * t324 + (-qJD(4) * t207 + t296) * t285 + t327;
t24 = t282 * t61 + t283 * t65;
t370 = mrSges(4,3) * t247;
t239 = Ifges(4,4) * t247;
t365 = t312 * Ifges(6,6);
t364 = t150 * Ifges(6,5);
t363 = t213 * Ifges(5,6);
t361 = t214 * Ifges(5,5);
t360 = t231 * Ifges(7,3);
t358 = t248 * Ifges(4,1);
t352 = Ifges(3,5) * qJD(2);
t163 = t221 * t353 - t354 * t290;
t206 = -t354 * t271 - t272 * t353;
t349 = t163 * t206;
t345 = t263 * t285;
t144 = t288 * t197 - t207 * t285;
t114 = pkin(4) * t261 - t288 * t350 + t144;
t123 = -qJ(5) * t345 + t145;
t71 = t282 * t114 + t283 * t123;
t326 = t354 * pkin(2);
t9 = -t31 * mrSges(7,1) + t30 * mrSges(7,2);
t320 = t289 * t333;
t23 = -t282 * t65 + t283 * t61;
t316 = t232 * mrSges(4,1) - t233 * mrSges(4,2);
t66 = -t108 * mrSges(6,1) + t109 * mrSges(6,2);
t70 = t283 * t114 - t123 * t282;
t279 = -t326 - pkin(3);
t309 = mrSges(5,1) * t288 - mrSges(5,2) * t285;
t306 = Ifges(5,1) * t285 + t366;
t304 = Ifges(5,2) * t288 + t367;
t302 = Ifges(5,5) * t285 + Ifges(5,6) * t288;
t186 = t297 * t263;
t52 = pkin(5) * t261 + pkin(9) * t186 + t70;
t185 = t262 * t263;
t60 = -pkin(9) * t185 + t71;
t21 = -t284 * t60 + t287 * t52;
t22 = t284 * t52 + t287 * t60;
t301 = -t285 * t68 - t288 * t69;
t180 = t244 * t353 - t354 * t245;
t300 = -t120 * t288 - t121 * t285;
t299 = t120 * t285 - t121 * t288;
t169 = -mrSges(5,2) * t240 + mrSges(5,3) * t213;
t170 = mrSges(5,1) * t240 - mrSges(5,3) * t214;
t298 = t169 * t288 - t170 * t285;
t126 = -t185 * t287 + t186 * t284;
t127 = -t185 * t284 - t186 * t287;
t179 = pkin(4) * t345 + t206;
t268 = -t288 * pkin(4) + t279;
t132 = t180 + (-t251 * t285 + t324) * pkin(4);
t116 = -t166 * pkin(4) + t163;
t281 = Ifges(3,4) * t337;
t256 = Ifges(3,1) * t338 + t281 + t352;
t219 = -qJD(2) * mrSges(4,2) + t370;
t212 = pkin(5) * t297 + t268;
t195 = -mrSges(4,1) * t247 - mrSges(4,2) * t248;
t188 = t239 - t358 + t444;
t141 = -mrSges(5,2) * t232 + mrSges(5,3) * t166;
t140 = mrSges(5,1) * t232 - mrSges(5,3) * t165;
t137 = t240 * Ifges(5,3) + t361 + t363;
t134 = mrSges(6,1) * t240 - mrSges(6,3) * t150;
t133 = -mrSges(6,2) * t240 + mrSges(6,3) * t312;
t129 = t185 * pkin(5) + t179;
t128 = pkin(5) * t150 + t381;
t125 = -t246 * t263 + t251 * t297;
t124 = t250 * t263 + t251 * t262;
t113 = -mrSges(5,1) * t166 + mrSges(5,2) * t165;
t100 = -mrSges(6,1) * t312 + mrSges(6,2) * t150;
t87 = t165 * Ifges(5,1) + t166 * Ifges(5,4) + t232 * Ifges(5,5);
t86 = t165 * Ifges(5,4) + t166 * Ifges(5,2) + t232 * Ifges(5,6);
t84 = mrSges(6,1) * t232 - mrSges(6,3) * t109;
t83 = -mrSges(6,2) * t232 + mrSges(6,3) * t108;
t79 = t240 * Ifges(6,3) + t364 + t365;
t76 = mrSges(7,1) * t231 - mrSges(7,3) * t94;
t75 = -mrSges(7,2) * t231 + mrSges(7,3) * t457;
t74 = -qJD(4) * t145 + t313;
t73 = -t207 * t335 + t327;
t72 = -t124 * pkin(5) + t132;
t67 = -t108 * pkin(5) + t116;
t53 = -mrSges(7,1) * t457 + mrSges(7,2) * t94;
t45 = t360 + t372 + t373;
t44 = -qJD(6) * t127 + t124 * t287 - t125 * t284;
t43 = qJD(6) * t126 + t124 * t284 + t125 * t287;
t26 = -mrSges(7,2) * t232 + mrSges(7,3) * t31;
t25 = mrSges(7,1) * t232 - mrSges(7,3) * t30;
t20 = pkin(9) * t124 + t24;
t19 = pkin(5) * t249 - pkin(9) * t125 + t23;
t5 = -qJD(6) * t22 + t19 * t287 - t20 * t284;
t4 = qJD(6) * t21 + t19 * t284 + t20 * t287;
t1 = [t280 * t316 + (-t11 * t43 + t12 * t44 + t126 * t2 - t127 * t3) * mrSges(7,3) + (-Ifges(4,4) * t232 - Ifges(4,1) * t233 + t305 * t397 + t307 * t398 + t86 * t384 + t87 * t383 + (mrSges(4,3) + t308) * t163 + t301 * mrSges(5,3) + (t302 * t389 + t306 * t394 + t304 * t395 + t191 * t309 + t139 * t384 - t288 * t138 / 0.2e1 + t299 * mrSges(5,3)) * qJD(4)) * t263 + (-mrSges(4,3) * t233 + t113) * t206 + (mrSges(4,1) * qJD(1) * t455 - mrSges(4,3) * t164 + Ifges(4,4) * t233 + Ifges(6,5) * t403 + Ifges(7,5) * t421 + Ifges(6,6) * t404 + Ifges(7,6) * t420 + (Ifges(6,3) + Ifges(7,3)) * t390 + t429 + t453) * t261 - t186 * t414 - t185 * t415 + t43 * t416 + t44 * t418 + t127 * t423 + t126 * t424 + (t422 + Ifges(4,2)) * t232 * t261 + t340 * t180 + (t137 / 0.2e1 + t360 / 0.2e1 + t45 / 0.2e1 + t79 / 0.2e1 + t372 / 0.2e1 + t373 / 0.2e1 + t365 / 0.2e1 + t364 / 0.2e1 + t427 + t363 / 0.2e1 + t361 / 0.2e1 + (t422 + Ifges(6,3) / 0.2e1) * t240 + t430) * t249 + (Ifges(7,1) * t127 + Ifges(7,4) * t126) * t421 + (Ifges(7,4) * t127 + Ifges(7,2) * t126) * t420 + t181 * t219 + t179 * t66 + t73 * t169 + t74 * t170 + t145 * t141 + t146 * (-mrSges(6,1) * t124 + mrSges(6,2) * t125) + t144 * t140 + t132 * t100 + t24 * t133 + t23 * t134 + t67 * (-mrSges(7,1) * t126 + mrSges(7,2) * t127) + t129 * t9 + t85 * (-mrSges(7,1) * t44 + mrSges(7,2) * t43) + t71 * t83 + t70 * t84 + t4 * t75 + t5 * t76 + t72 * t53 + (t198 * t251 - t199 * t249 - t207 * t232) * mrSges(4,3) + (-Ifges(6,4) * t186 - Ifges(6,2) * t185) * t404 + (-Ifges(6,1) * t186 - Ifges(6,4) * t185) * t403 + (-Ifges(6,5) * t186 + Ifges(7,5) * t127 - Ifges(6,6) * t185 + Ifges(7,6) * t126 + t263 * t303) * t390 + (t124 * t49 - t125 * t48 + t15 * t186 - t16 * t185) * mrSges(6,3) + t116 * (mrSges(6,1) * t185 - mrSges(6,2) * t186) + t21 * t25 + t22 * t26 + (Ifges(6,1) * t125 + Ifges(6,4) * t124) * t399 + (Ifges(6,4) * t125 + Ifges(6,2) * t124) * t401 + (Ifges(7,1) * t43 + Ifges(7,4) * t44) * t406 + (Ifges(7,4) * t43 + Ifges(7,2) * t44) * t408 + t125 * t410 + t124 * t412 + (Ifges(6,5) * t125 + Ifges(6,6) * t124) * t388 + (Ifges(7,5) * t43 + Ifges(7,6) * t44) * t391 + (-t351 / 0.2e1 + (mrSges(3,1) * t426 - 0.3e1 / 0.2e1 * t368 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t289) * qJD(1) + (qJD(1) * mrSges(4,2) * t263 + t195 + m(4) * (t267 + t339)) * pkin(2) - t432) * t336 + (t332 + t434 + t438) * t261 / 0.2e1 + m(4) * (t164 * t207 - t180 * t198 + t181 * t199 + t349) + m(5) * (t120 * t74 + t121 * t73 + t144 * t69 + t145 * t68 + t180 * t191 + t349) - (t449 - t358 / 0.2e1 + t239 / 0.2e1 + t188 / 0.2e1 + t300 * mrSges(5,3) + t428) * t251 + m(7) * (t11 * t5 + t12 * t4 + t129 * t67 + t2 * t22 + t21 * t3 + t72 * t85) + m(6) * (t116 * t179 + t132 * t146 + t15 * t70 + t16 * t71 + t23 * t48 + t24 * t49) + (-Ifges(4,5) * t251 / 0.2e1 + t249 * t456 + (-t379 + t256 / 0.2e1 + t352 / 0.2e1 + (mrSges(3,2) * t426 + 0.3e1 / 0.2e1 * Ifges(3,4) * t289) * qJD(1)) * t289) * qJD(2); (-t471 * t121 + (-t334 + t346) * t120 + t433) * mrSges(5,3) - t195 * t331 - Ifges(3,6) * t321 - t297 * t415 + (Ifges(6,5) * t262 + Ifges(7,5) * t201 - Ifges(6,6) * t297 + Ifges(7,6) * t200 + t302) * t390 - (Ifges(3,5) * t289 - Ifges(3,6) * t286) * t333 / 0.2e1 + (m(5) * t277 * t300 + t428) * qJD(4) + t116 * (mrSges(6,1) * t297 + mrSges(6,2) * t262) + (Ifges(6,1) * t262 - Ifges(6,4) * t297) * t403 + (Ifges(6,4) * t262 - Ifges(6,2) * t297) * t404 + (-t15 * t262 - t16 * t297 + t435 * t48 - t436 * t49) * mrSges(6,3) + (-t232 * t325 + t233 * t326) * mrSges(4,3) + (Ifges(7,5) * t119 + Ifges(7,6) * t118) * t392 + (Ifges(7,1) * t119 + Ifges(7,4) * t118) * t407 - t172 * t413 + t262 * t414 + t135 * t416 + t119 * t417 + t136 * t418 + t118 * t419 + (Ifges(7,4) * t201 + Ifges(7,2) * t200) * t420 + (Ifges(7,1) * t201 + Ifges(7,4) * t200) * t421 + t201 * t423 + t200 * t424 + (Ifges(7,4) * t119 + Ifges(7,2) * t118) * t409 + (-t449 - t444 / 0.2e1 + t307 * t394 + t305 * t395 + Ifges(4,1) * t386 + t303 * t389 - t458) * t247 + t279 * t113 + t268 * t66 - m(5) * (t120 * t130 + t121 * t131) + t285 * t87 / 0.2e1 - Ifges(4,6) * t232 - Ifges(4,5) * t233 + (t357 + t137 + t79 + t45) * t386 - t203 * t219 + t212 * t9 + t67 * (-mrSges(7,1) * t200 + mrSges(7,2) * t201) + t193 * t84 + t194 * t83 - t131 * t169 - t130 * t170 - t164 * mrSges(4,2) + Ifges(3,5) * t320 + t103 * t25 + t104 * t26 + (-Ifges(6,1) * t250 - Ifges(6,4) * t246) * t399 + (-Ifges(6,5) * t250 - Ifges(6,6) * t246) * t388 + (-Ifges(6,4) * t250 - Ifges(6,2) * t246) * t401 + (-Ifges(6,5) * t173 - Ifges(6,6) * t172) * t389 + (-Ifges(6,1) * t173 - Ifges(6,4) * t172) * t400 + (-Ifges(6,4) * t173 - Ifges(6,2) * t172) * t402 - (-Ifges(3,2) * t338 + t256 + t281) * t337 / 0.2e1 + (-t286 * (Ifges(3,1) * t289 - t368) / 0.2e1 + pkin(1) * (mrSges(3,1) * t286 + mrSges(3,2) * t289)) * qJD(1) ^ 2 + (t239 + t342 + t188) * t387 + (t198 * t202 - t199 * t203 - t267 * t331 + (-t163 * t354 + t164 * t353) * pkin(2)) * m(4) + (Ifges(7,1) * t135 + Ifges(7,4) * t136) * t406 + (Ifges(7,4) * t135 + Ifges(7,2) * t136) * t408 - t250 * t410 - t173 * t411 - t246 * t412 + t304 * t397 + t306 * t398 + (Ifges(7,5) * t135 + Ifges(7,6) * t136) * t391 + t198 * t370 + t337 * t379 + t86 * t383 + (m(5) * t279 - mrSges(4,1) - t309) * t163 + (-mrSges(3,1) * t320 + mrSges(3,2) * t321) * pkin(7) + (-Ifges(5,5) * t394 - Ifges(6,5) * t400 - Ifges(7,5) * t407 - Ifges(5,6) * t395 - Ifges(6,6) * t402 - Ifges(7,6) * t409 - Ifges(7,3) * t392 + t389 * t468 + t427 + t431 + t464) * t248 + t432 * t338 + (m(5) * t433 - t285 * t140 + t288 * t141 - t169 * t335 - t170 * t334) * t277 + (mrSges(6,1) * t436 - mrSges(6,2) * t435) * t146 + (-m(5) * t191 - t340) * t202 + t437 * t53 + t439 * t100 + (-mrSges(7,1) * t440 + mrSges(7,2) * t441) * t85 + (-t11 * t441 + t12 * t440 + t2 * t200 - t201 * t3) * mrSges(7,3) + t447 * t133 + t448 * t134 + (t116 * t268 + t146 * t439 + t15 * t193 + t16 * t194 + t447 * t49 + t448 * t48) * m(6) + t451 * t76 + (t103 * t3 + t104 * t2 + t11 * t451 + t12 * t452 + t212 * t67 + t437 * t85) * m(7) + t452 * t75 - t199 * t369; (-t219 - t298) * t247 + t298 * qJD(4) + t440 * t76 + t441 * t75 + (t53 + t100 + t340) * t248 + t288 * t140 + t285 * t141 + t262 * t83 - t297 * t84 + t200 * t25 + t201 * t26 - t436 * t134 - t435 * t133 + t316 + (t11 * t440 + t12 * t441 + t2 * t201 + t200 * t3 + t248 * t85) * m(7) + (t146 * t248 - t15 * t297 + t16 * t262 - t435 * t49 - t436 * t48) * m(6) + (t191 * t248 - t240 * t299 - t301) * m(5) + (-t198 * t248 - t199 * t247 + t310) * m(4); -t150 * t413 + (Ifges(5,5) * t213 + Ifges(6,5) * t312 - Ifges(5,6) * t214 - Ifges(6,6) * t150) * t389 - t146 * (mrSges(6,1) * t150 + mrSges(6,2) * t312) + (t150 * t49 + t312 * t48) * mrSges(6,3) + (-Ifges(6,2) * t150 + t450) * t402 + t434 + (Ifges(6,1) * t312 - t459) * t400 + (t120 * t213 + t121 * t214) * mrSges(5,3) - t94 * t419 + t472 + ((t15 * t283 + t16 * t282) * pkin(4) - t146 * t381 - t48 * t54 - t49 * t55) * m(6) + (-t100 * t214 + t282 * t83 + t283 * t84) * pkin(4) + t457 * t417 + t470 * t409 + t242 * t25 + t243 * t26 - t191 * (mrSges(5,1) * t214 + mrSges(5,2) * t213) - t120 * t169 + t121 * t170 - t55 * t133 - t54 * t134 - t128 * t53 + t312 * t411 + (-Ifges(5,2) * t214 + t139 + t208) * t395 + t138 * t393 + (Ifges(5,1) * t213 - t362) * t394 + t429 + t445 * t76 + t446 * t75 + (t11 * t445 + t12 * t446 - t128 * t85 + t2 * t243 + t242 * t3) * m(7); -t312 * t133 + t150 * t134 - t457 * t75 + t94 * t76 + t66 + t9 + (t11 * t94 - t12 * t457 + t67) * m(7) + (t150 * t48 - t312 * t49 + t116) * m(6); t46 * t406 - t11 * t75 + t12 * t76 + (t47 + t470) * t409 + t472;];
tauc  = t1(:);
