% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:38:12
% EndTime: 2019-03-09 10:38:59
% DurationCPUTime: 25.63s
% Computational Cost: add. (12382->725), mult. (37316->954), div. (0->0), fcn. (29275->10), ass. (0->342)
t245 = sin(qJ(2));
t248 = cos(qJ(2));
t241 = sin(pkin(6));
t341 = qJD(1) * t241;
t356 = sin(pkin(11));
t297 = t356 * t341;
t357 = cos(pkin(11));
t298 = t357 * t341;
t199 = -t245 * t298 - t248 * t297;
t244 = sin(qJ(4));
t247 = cos(qJ(4));
t242 = cos(pkin(6));
t307 = qJD(1) * t242 + qJD(2);
t164 = -t199 * t244 - t247 * t307;
t406 = -t164 / 0.2e1;
t165 = -t247 * t199 + t244 * t307;
t163 = qJD(6) + t165;
t407 = t163 / 0.2e1;
t198 = -t245 * t297 + t248 * t298;
t195 = qJD(4) - t198;
t243 = sin(qJ(6));
t246 = cos(qJ(6));
t119 = t164 * t243 + t195 * t246;
t410 = t119 / 0.2e1;
t118 = t164 * t246 - t195 * t243;
t412 = t118 / 0.2e1;
t488 = Ifges(5,4) * t406 + Ifges(7,5) * t410 + Ifges(7,6) * t412 + Ifges(7,3) * t407;
t491 = Ifges(6,6) * t406 + t488;
t403 = t165 / 0.2e1;
t490 = Ifges(5,1) * t403;
t399 = t195 / 0.2e1;
t400 = -t195 / 0.2e1;
t391 = pkin(1) * t242;
t236 = t248 * t391;
t232 = qJD(1) * t236;
t382 = pkin(8) + qJ(3);
t319 = t382 * t245;
t304 = t241 * t319;
t266 = t242 * pkin(2) - t304;
t166 = qJD(2) * pkin(2) + qJD(1) * t266 + t232;
t235 = t245 * t391;
t349 = t241 * t248;
t443 = t349 * t382 + t235;
t182 = t443 * qJD(1);
t314 = t357 * t182;
t107 = t356 * t166 + t314;
t100 = pkin(9) * t307 + t107;
t222 = (-pkin(2) * t248 - pkin(1)) * t241;
t342 = qJD(1) * t222;
t216 = qJD(3) + t342;
t255 = -pkin(3) * t198 + pkin(9) * t199 + t216;
t59 = t100 * t244 - t247 * t255;
t457 = -qJD(5) - t59;
t45 = -pkin(4) * t195 - t457;
t418 = pkin(4) + pkin(10);
t274 = pkin(5) * t165 + t59;
t477 = qJD(5) + t274;
t34 = -t195 * t418 + t477;
t172 = t356 * t182;
t106 = t166 * t357 - t172;
t99 = -pkin(3) * t307 - t106;
t252 = -t165 * qJ(5) + t99;
t39 = t164 * t418 + t252;
t5 = -t243 * t39 + t246 * t34;
t6 = t243 * t34 + t246 * t39;
t485 = t45 * mrSges(6,1) + t5 * mrSges(7,1) - t6 * mrSges(7,2) + t59 * mrSges(5,3) + Ifges(6,4) * t400 + Ifges(5,5) * t399 + Ifges(6,2) * t403 + t490 + t491;
t58 = t164 * pkin(4) + t252;
t489 = t99 * mrSges(5,2) - t58 * mrSges(6,3) + t485;
t470 = Ifges(5,6) - Ifges(6,5);
t487 = Ifges(6,6) + Ifges(5,4);
t181 = -qJD(1) * t304 + t232;
t134 = t181 * t356 + t314;
t337 = qJD(4) * t244;
t354 = t198 * t244;
t486 = -qJD(5) * t244 - t134 + (t337 - t354) * pkin(4);
t340 = qJD(2) * t241;
t450 = -t356 * t245 + t357 * t248;
t201 = t450 * t340;
t256 = qJD(1) * t201;
t453 = qJD(4) * t307 + t256;
t116 = -t199 * t337 - t247 * t453;
t417 = -t116 / 0.2e1;
t336 = qJD(4) * t247;
t117 = -t199 * t336 + t244 * t453;
t415 = -t117 / 0.2e1;
t205 = (t245 * t357 + t248 * t356) * t241;
t200 = qJD(2) * t205;
t192 = qJD(1) * t200;
t401 = t192 / 0.2e1;
t483 = Ifges(6,4) - Ifges(5,5);
t469 = Ifges(5,3) + Ifges(6,1);
t482 = t486 + t195 * (pkin(10) * t244 - qJ(5) * t247);
t135 = t181 * t357 - t172;
t125 = t244 * t135;
t328 = t245 * t341;
t309 = pkin(2) * t328;
t147 = -pkin(3) * t199 - pkin(9) * t198 + t309;
t329 = t356 * pkin(2);
t238 = t329 + pkin(9);
t383 = pkin(5) + t238;
t221 = t383 * t247;
t481 = qJD(4) * t221 - t125 - (pkin(5) * t198 - t147) * t247 - t418 * t199;
t60 = t247 * t100 + t244 * t255;
t48 = -t195 * qJ(5) - t60;
t445 = t48 * mrSges(6,1) - t60 * mrSges(5,3);
t480 = -t99 * mrSges(5,1) + t58 * mrSges(6,2) - t445;
t416 = t116 / 0.2e1;
t414 = t117 / 0.2e1;
t229 = qJD(2) * t232;
t262 = (-qJD(2) * t319 + qJD(3) * t248) * t241;
t158 = qJD(1) * t262 + t229;
t350 = t241 * t245;
t168 = -qJD(2) * t443 - qJD(3) * t350;
t250 = qJD(1) * t168;
t95 = t158 * t357 + t250 * t356;
t479 = qJD(4) * t255 + t95;
t94 = t158 * t356 - t357 * t250;
t260 = t116 * qJ(5) - t165 * qJD(5) + t94;
t31 = t117 * pkin(4) + t260;
t52 = -qJD(6) * t119 + t117 * t246 - t192 * t243;
t424 = t52 / 0.2e1;
t51 = qJD(6) * t118 + t117 * t243 + t192 * t246;
t425 = t51 / 0.2e1;
t21 = t117 * t418 + t260;
t326 = t245 * t340;
t303 = qJD(1) * t326;
t276 = pkin(2) * t303;
t131 = t192 * pkin(3) - pkin(9) * t256 + t276;
t24 = -t100 * t336 + t131 * t247 - t244 * t479;
t7 = -pkin(5) * t116 - t192 * t418 - t24;
t1 = qJD(6) * t5 + t21 * t246 + t243 * t7;
t2 = -qJD(6) * t6 - t21 * t243 + t246 * t7;
t442 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t473 = -Ifges(6,4) / 0.2e1;
t9 = Ifges(7,5) * t51 + Ifges(7,6) * t52 - Ifges(7,3) * t116;
t476 = mrSges(5,2) * t94 - mrSges(6,3) * t31 + Ifges(7,5) * t425 + Ifges(7,6) * t424 + t192 * t473 + t442 - Ifges(6,2) * t416 + t9 / 0.2e1 + (0.2e1 * Ifges(5,1) + Ifges(7,3) + Ifges(6,2)) * t417 + t487 * t415 + (-t483 + Ifges(5,5)) * t401;
t474 = -0.2e1 * Ifges(6,3) * t414 - t94 * mrSges(5,1) + t31 * mrSges(6,2) - t487 * t416 - (t414 - t415) * Ifges(5,2) + 0.2e1 * t401 * t470;
t330 = t357 * pkin(2);
t239 = -t330 - pkin(3);
t268 = -t244 * qJ(5) + t239;
t206 = -t247 * t418 + t268;
t220 = t383 * t244;
t157 = t206 * t246 + t220 * t243;
t468 = -qJD(6) * t157 - t243 * t482 + t246 * t481;
t156 = -t206 * t243 + t220 * t246;
t467 = qJD(6) * t156 + t243 * t481 + t246 * t482;
t466 = -t164 * t470 - t165 * t483 + t195 * t469;
t194 = Ifges(4,4) * t198;
t464 = Ifges(4,5) * t256;
t463 = Ifges(4,2) * t198;
t462 = Ifges(4,6) * t192;
t461 = t307 * Ifges(4,5);
t460 = t307 * Ifges(4,6);
t72 = t247 * t135 + t244 * t147;
t61 = qJ(5) * t199 - t72;
t459 = pkin(5) * t354 - t383 * t337 + t61;
t353 = t198 * t247;
t458 = (-t336 + t353) * qJ(5) + t486;
t364 = t199 * mrSges(4,3);
t456 = -mrSges(4,1) * t307 + mrSges(5,1) * t164 + mrSges(5,2) * t165 - t364;
t152 = t199 * t243 + t246 * t354;
t335 = qJD(6) * t247;
t322 = t243 * t335;
t455 = -t246 * t337 + t152 - t322;
t153 = -t199 * t246 + t243 * t354;
t454 = -t243 * t337 + t246 * t335 + t153;
t452 = t99 * (mrSges(5,1) * t244 + mrSges(5,2) * t247) + t58 * (-mrSges(6,2) * t244 - mrSges(6,3) * t247);
t451 = -t244 * t470 - t247 * t483;
t449 = t116 * t483 - t117 * t470 + t192 * t469;
t334 = qJD(1) * qJD(2);
t318 = t248 * t334;
t302 = t241 * t318;
t448 = Ifges(3,5) * t302 - Ifges(3,6) * t303 - t462 + t464;
t23 = -t100 * t337 + t244 * t131 + t247 * t479;
t447 = t23 * t247 - t24 * t244;
t14 = -qJ(5) * t192 - qJD(5) * t195 - t23;
t20 = -pkin(4) * t192 - t24;
t446 = -t14 * t247 + t20 * t244;
t296 = t1 * t243 + t2 * t246;
t80 = t195 * Ifges(6,5) - t165 * Ifges(6,6) + t164 * Ifges(6,3);
t83 = t165 * Ifges(5,4) - t164 * Ifges(5,2) + t195 * Ifges(5,6);
t444 = t83 / 0.2e1 - t80 / 0.2e1;
t78 = -mrSges(7,2) * t163 + mrSges(7,3) * t118;
t79 = mrSges(7,1) * t163 - mrSges(7,3) * t119;
t279 = t243 * t79 - t246 * t78;
t32 = -mrSges(7,1) * t116 - mrSges(7,3) * t51;
t33 = mrSges(7,2) * t116 + mrSges(7,3) * t52;
t280 = t243 * t33 + t246 * t32;
t89 = -t116 * mrSges(6,1) + t192 * mrSges(6,2);
t439 = t279 * qJD(6) - t280 - t89;
t437 = t216 * mrSges(4,2) + t461 / 0.2e1;
t436 = t24 * mrSges(5,1) - t23 * mrSges(5,2) + t20 * mrSges(6,2) - t14 * mrSges(6,3);
t202 = -pkin(8) * t303 + t229;
t343 = pkin(8) * t349 + t235;
t212 = t343 * qJD(2);
t203 = qJD(1) * t212;
t435 = -t203 * mrSges(3,1) - t94 * mrSges(4,1) - t202 * mrSges(3,2) - t95 * mrSges(4,2);
t284 = Ifges(7,5) * t243 + Ifges(7,6) * t246;
t373 = Ifges(7,4) * t243;
t287 = Ifges(7,2) * t246 + t373;
t372 = Ifges(7,4) * t246;
t291 = Ifges(7,1) * t243 + t372;
t294 = mrSges(7,1) * t246 - mrSges(7,2) * t243;
t295 = t243 * t5 - t246 * t6;
t374 = Ifges(7,4) * t119;
t43 = Ifges(7,2) * t118 + Ifges(7,6) * t163 + t374;
t359 = t246 * t43;
t389 = t164 * pkin(5);
t36 = -t48 - t389;
t113 = Ifges(7,4) * t118;
t44 = Ifges(7,1) * t119 + Ifges(7,5) * t163 + t113;
t360 = t243 * t44;
t408 = -t163 / 0.2e1;
t411 = -t119 / 0.2e1;
t413 = -t118 / 0.2e1;
t433 = t284 * t408 + t287 * t413 + t291 * t411 + t36 * t294 + t295 * mrSges(7,3) - t360 / 0.2e1 - t359 / 0.2e1;
t431 = t216 * mrSges(4,1) + t45 * mrSges(6,2) - t460 / 0.2e1 - t48 * mrSges(6,3) - t59 * mrSges(5,1) - t60 * mrSges(5,2);
t429 = t241 ^ 2;
t11 = Ifges(7,1) * t51 + Ifges(7,4) * t52 - Ifges(7,5) * t116;
t427 = t11 / 0.2e1;
t426 = t43 / 0.2e1;
t420 = -t83 / 0.2e1;
t405 = t164 / 0.2e1;
t404 = -t165 / 0.2e1;
t398 = -t198 / 0.2e1;
t397 = -t199 / 0.2e1;
t19 = -mrSges(7,1) * t52 + mrSges(7,2) * t51;
t88 = mrSges(6,1) * t117 - mrSges(6,3) * t192;
t381 = t19 - t88;
t380 = mrSges(4,3) * t192;
t379 = Ifges(3,4) * t245;
t378 = Ifges(3,4) * t248;
t377 = Ifges(4,4) * t199;
t376 = Ifges(5,4) * t244;
t375 = Ifges(5,4) * t247;
t371 = Ifges(6,6) * t244;
t370 = Ifges(6,6) * t247;
t365 = t198 * mrSges(4,3);
t127 = mrSges(6,1) * t164 - mrSges(6,3) * t195;
t70 = -mrSges(7,1) * t118 + mrSges(7,2) * t119;
t358 = t70 - t127;
t355 = qJ(5) * t164;
t348 = t243 * t247;
t347 = t246 * t247;
t178 = t236 + t266;
t197 = qJ(3) * t349 + t343;
t143 = t356 * t178 + t357 * t197;
t133 = pkin(9) * t242 + t143;
t204 = t450 * t241;
t151 = -pkin(3) * t204 - pkin(9) * t205 + t222;
t75 = t247 * t133 + t244 * t151;
t129 = -mrSges(5,2) * t195 - mrSges(5,3) * t164;
t346 = -t127 + t129;
t128 = mrSges(6,1) * t165 + mrSges(6,2) * t195;
t130 = mrSges(5,1) * t195 - mrSges(5,3) * t165;
t345 = -t128 + t130;
t327 = t248 * t341;
t325 = t248 * t340;
t324 = t238 * t337;
t323 = t238 * t336;
t311 = t192 * mrSges(4,1) + mrSges(4,2) * t256;
t74 = -t244 * t133 + t151 * t247;
t71 = t147 * t247 - t125;
t308 = pkin(2) * t326;
t306 = mrSges(3,3) * t328;
t305 = mrSges(3,3) * t327;
t299 = t325 / 0.2e1;
t65 = qJ(5) * t204 - t75;
t293 = Ifges(5,1) * t247 - t376;
t292 = Ifges(7,1) * t246 - t373;
t290 = -Ifges(5,2) * t244 + t375;
t288 = -Ifges(7,2) * t243 + t372;
t285 = Ifges(7,5) * t246 - Ifges(7,6) * t243;
t283 = -Ifges(6,2) * t247 + t371;
t282 = Ifges(6,3) * t244 - t370;
t176 = t205 * t247 + t242 * t244;
t38 = pkin(5) * t176 + t204 * t418 - t74;
t175 = t205 * t244 - t242 * t247;
t142 = t178 * t357 - t356 * t197;
t132 = -t242 * pkin(3) - t142;
t261 = -t176 * qJ(5) + t132;
t57 = t175 * t418 + t261;
t15 = -t243 * t57 + t246 * t38;
t16 = t243 * t38 + t246 * t57;
t233 = qJD(2) * t236;
t167 = t233 + t262;
t104 = t167 * t356 - t357 * t168;
t149 = t175 * t246 + t204 * t243;
t150 = t175 * t243 - t204 * t246;
t105 = t167 * t357 + t168 * t356;
t148 = pkin(3) * t200 - pkin(9) * t201 + t308;
t30 = -t244 * t105 - t133 * t336 + t148 * t247 - t151 * t337;
t29 = t247 * t105 - t133 * t337 + t244 * t148 + t151 * t336;
t267 = pkin(1) * t429 * (mrSges(3,1) * t245 + mrSges(3,2) * t248);
t265 = t245 * t429 * (Ifges(3,1) * t248 - t379);
t141 = -t205 * t337 + (qJD(4) * t242 + t201) * t247;
t259 = -t141 * qJ(5) - t176 * qJD(5) + t104;
t258 = t241 * t307 * (Ifges(3,5) * t248 - Ifges(3,6) * t245);
t257 = (Ifges(3,6) * t242 + (Ifges(3,2) * t248 + t379) * t241) * qJD(1);
t22 = -qJ(5) * t200 + qJD(5) * t204 - t29;
t253 = mrSges(4,3) * t256;
t231 = Ifges(3,4) * t327;
t219 = -t247 * pkin(4) + t268;
t217 = -pkin(8) * t350 + t236;
t211 = -pkin(8) * t326 + t233;
t210 = t343 * qJD(1);
t209 = -pkin(8) * t328 + t232;
t208 = -mrSges(3,2) * t307 + t305;
t207 = mrSges(3,1) * t307 - t306;
t180 = Ifges(3,1) * t328 + Ifges(3,5) * t307 + t231;
t179 = Ifges(3,6) * qJD(2) + t257;
t169 = -mrSges(4,2) * t307 + t365;
t154 = -mrSges(4,1) * t198 - mrSges(4,2) * t199;
t146 = -Ifges(4,1) * t199 + t194 + t461;
t145 = -t377 + t460 + t463;
t140 = qJD(4) * t176 + t201 * t244;
t103 = -mrSges(6,2) * t164 - mrSges(6,3) * t165;
t101 = pkin(4) * t165 + t355;
t87 = -mrSges(5,2) * t192 - mrSges(5,3) * t117;
t86 = mrSges(5,1) * t192 + mrSges(5,3) * t116;
t76 = t165 * t418 + t355;
t73 = t175 * pkin(4) + t261;
t69 = mrSges(5,1) * t117 - mrSges(5,2) * t116;
t68 = -mrSges(6,2) * t117 + mrSges(6,3) * t116;
t66 = pkin(4) * t204 - t74;
t64 = qJD(6) * t149 + t140 * t243 + t200 * t246;
t63 = -qJD(6) * t150 + t140 * t246 - t200 * t243;
t62 = pkin(4) * t199 - t71;
t46 = -pkin(5) * t175 - t65;
t41 = t60 - t389;
t35 = t140 * pkin(4) + t259;
t28 = t140 * t418 + t259;
t27 = -pkin(4) * t200 - t30;
t26 = t243 * t41 + t246 * t76;
t25 = -t243 * t76 + t246 * t41;
t13 = -pkin(5) * t140 - t22;
t12 = pkin(5) * t141 - t200 * t418 - t30;
t10 = Ifges(7,4) * t51 + Ifges(7,2) * t52 - Ifges(7,6) * t116;
t8 = -pkin(5) * t117 - t14;
t4 = -qJD(6) * t16 + t12 * t246 - t243 * t28;
t3 = qJD(6) * t15 + t12 * t243 + t246 * t28;
t17 = [(-t107 * mrSges(4,3) + t466 / 0.2e1 - Ifges(4,4) * t397 + Ifges(5,5) * t403 + Ifges(6,4) * t404 + Ifges(6,5) * t405 + Ifges(5,6) * t406 - t145 / 0.2e1 - t463 / 0.2e1 + t469 * t399 + t431) * t200 - (t257 + t179) * t326 / 0.2e1 + (t448 / 0.2e1 + qJD(1) * Ifges(3,5) * t299 + t464 / 0.2e1 - t462 / 0.2e1 + t435) * t242 + (-m(4) * t106 + m(5) * t99 + t456) * t104 + (t14 * mrSges(6,1) - t23 * mrSges(5,3) - Ifges(5,4) * t417 + Ifges(6,6) * t416 - t474) * t175 + (-Ifges(5,4) * t403 + Ifges(6,6) * t404 + Ifges(6,3) * t405 - Ifges(5,2) * t406 + t420 + t80 / 0.2e1 - t470 * t399 - t480) * t140 + (mrSges(6,1) * t20 - mrSges(5,3) * t24 + Ifges(5,4) * t415 - Ifges(6,6) * t414 + t476) * t176 + (t180 + (Ifges(3,1) * t245 + t378) * t341) * t299 - (t449 / 0.2e1 + mrSges(4,1) * t276 - mrSges(4,3) * t95 - Ifges(4,4) * t256 + Ifges(6,4) * t416 + Ifges(5,5) * t417 + Ifges(6,5) * t414 + Ifges(4,2) * t192 + Ifges(5,6) * t415 + t469 * t401 + t436) * t204 + (t202 * t349 + t203 * t350 - t209 * t325 - t210 * t326 - t217 * t302 - t303 * t343) * mrSges(3,3) + m(3) * (t202 * t343 - t203 * t217 - t209 * t212 + t210 * t211) + m(5) * (t132 * t94 + t23 * t75 + t24 * t74 + t29 * t60 - t30 * t59) + m(4) * (t105 * t107 - t142 * t94 + t143 * t95 + (t216 + t342) * t308) + m(7) * (t1 * t16 + t13 * t36 + t15 * t2 + t3 * t6 + t4 * t5 + t46 * t8) + m(6) * (t14 * t65 + t20 * t66 + t22 * t48 + t27 * t45 + t31 * t73 + t35 * t58) + (Ifges(7,5) * t64 + Ifges(7,6) * t63) * t407 + (Ifges(7,5) * t150 + Ifges(7,6) * t149) * t417 + (mrSges(4,2) * t276 + mrSges(4,3) * t94 + Ifges(4,1) * t256 - Ifges(4,4) * t192) * t205 - t143 * t380 - t142 * t253 + t154 * t308 + (Ifges(7,1) * t64 + Ifges(7,4) * t63) * t410 + (Ifges(7,1) * t150 + Ifges(7,4) * t149) * t425 + (Ifges(7,4) * t64 + Ifges(7,2) * t63) * t412 + (Ifges(7,4) * t150 + Ifges(7,2) * t149) * t424 + t429 * (-Ifges(3,2) * t245 + t378) * t318 + qJD(2) * t258 / 0.2e1 + (t1 * t149 - t150 * t2 - t5 * t64 + t6 * t63) * mrSges(7,3) + t222 * t311 + (-0.2e1 * t267 + t265) * t334 + t211 * t208 - t212 * t207 + (-t106 * mrSges(4,3) + Ifges(4,1) * t397 + t146 / 0.2e1 + t194 / 0.2e1 + t437) * t201 + (-Ifges(6,2) * t404 - Ifges(6,6) * t405 - t399 * t483 + t488 + t489 + t490) * t141 + t15 * t32 + t16 * t33 + t46 * t19 + t64 * t44 / 0.2e1 + t36 * (-mrSges(7,1) * t63 + mrSges(7,2) * t64) + t13 * t70 + t73 * t68 + t3 * t78 + t4 * t79 + t74 * t86 + t75 * t87 + t65 * t88 + t66 * t89 + t35 * t103 + t63 * t426 + t150 * t427 + t22 * t127 + t27 * t128 + t29 * t129 + t30 * t130 + t132 * t69 + t149 * t10 / 0.2e1 + t8 * (-mrSges(7,1) * t149 + mrSges(7,2) * t150) + t105 * t169; (t306 + t207) * t210 + t459 * t70 + (t305 - t208) * t209 + (-t323 - t71) * t130 + (-t324 - t72) * t129 + (t323 - t62) * t128 + (t324 - t61) * t127 + (-t258 / 0.2e1 + (t267 - t265 / 0.2e1) * qJD(1)) * qJD(1) + (Ifges(7,5) * t153 + Ifges(7,6) * t152) * t408 - (-Ifges(3,2) * t328 + t180 + t231) * t327 / 0.2e1 + (t80 + t359 + t360) * t337 / 0.2e1 - (qJD(6) * t44 + t10) * t347 / 0.2e1 + (Ifges(7,4) * t153 + Ifges(7,2) * t152) * t413 + (Ifges(7,1) * t153 + Ifges(7,4) * t152) * t411 + (t451 * t399 + (t293 / 0.2e1 - t283 / 0.2e1) * t165 + (Ifges(7,3) * t247 + t244 * t284) * t407 + (Ifges(7,5) * t247 + t244 * t291) * t410 + (Ifges(7,6) * t247 + t244 * t287) * t412 + t452 + (-t290 / 0.2e1 + t282 / 0.2e1) * t164 + (m(5) * (-t244 * t60 + t247 * t59) + m(6) * (t244 * t48 + t247 * t45)) * t238) * qJD(4) + ((t356 * t95 - t357 * t94) * pkin(2) + t106 * t134 - t107 * t135 - t216 * t309) * m(4) + (mrSges(7,1) * t455 - mrSges(7,2) * t454) * t36 + (-t1 * t347 + t2 * t348 + t454 * t5 - t455 * t6) * mrSges(7,3) - t456 * t134 + (t420 + t445) * t337 + (t282 * t406 + t283 * t403 + t290 * t405 + t293 * t404 + t400 * t451 - t437 - t452) * t198 + ((-t88 + t87) * t238 - t284 * t417 - t287 * t424 - t291 * t425 + t8 * t294 + t474) * t247 + ((-t86 + t89) * t238 + t476) * t244 + t444 * t354 + (t194 + t146) * t398 + t435 + (-t285 * t407 - t288 * t412 - t292 * t410) * t335 - t371 * t414 + t376 * t415 - t370 * t416 + t375 * t417 - t329 * t380 + (-t354 * t48 + t446) * mrSges(6,1) + (t354 * t60 + t447) * mrSges(5,3) - t11 * t348 / 0.2e1 + (-t134 * t99 + t238 * t447 + t239 * t94 + t59 * t71 - t60 * t72) * m(5) + (t219 * t31 + t238 * t446 - t45 * t62 + t458 * t58 - t48 * t61) * m(6) - t154 * t309 + t458 * t103 - t107 * t364 + t106 * t365 + t179 * t328 / 0.2e1 - t253 * t330 + (Ifges(7,5) * t411 + Ifges(7,6) * t413 + Ifges(7,3) * t408 - t485) * t353 + t485 * t336 + t448 + t239 * t69 + t219 * t68 + t221 * t19 + t145 * t397 + (Ifges(4,1) * t198 + t377 + t466) * t199 / 0.2e1 + t467 * t78 + t468 * t79 + (t1 * t157 + t156 * t2 + t221 * t8 + t36 * t459 + t467 * t6 + t468 * t5) * m(7) + (-Ifges(6,4) * t403 - Ifges(5,5) * t404 - Ifges(6,5) * t406 + Ifges(4,2) * t398 - Ifges(5,6) * t405 - t400 * t469 + t431) * t199 + t322 * t426 - t152 * t43 / 0.2e1 - t153 * t44 / 0.2e1 + t156 * t32 + t157 * t33 - t135 * t169; -t152 * t79 - t153 * t78 - t198 * t169 + (t103 + t456) * t199 + (t87 + t345 * t198 + (t243 * t78 + t246 * t79 - t345) * qJD(4) + t381) * t244 + (t86 - t195 * (-t70 - t346) + t439) * t247 + t311 + ((t8 + (t243 * t6 + t246 * t5) * qJD(4)) * t244 + (qJD(4) * t36 + qJD(6) * t295 - t296) * t247 - t152 * t5 - t153 * t6 - t36 * t353) * m(7) + (-t14 * t244 + t199 * t58 - t20 * t247 + t195 * (t244 * t45 - t247 * t48)) * m(6) + (t199 * t99 + t23 * t244 + t24 * t247 + t195 * (t244 * t59 + t247 * t60)) * m(5) + (-t106 * t199 - t107 * t198 + t276) * m(4); (t8 * qJ(5) - t25 * t5 - t26 * t6 + t477 * t36) * m(7) - (m(7) * t296 + t280 + (-m(7) * t295 - t279) * qJD(6)) * t418 + ((Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1) * t165 + (-Ifges(6,5) / 0.2e1 + Ifges(5,6) / 0.2e1) * t195 + (-Ifges(6,3) / 0.2e1 - Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1) * t164 + t433 + t444 + t480) * t165 + t381 * qJ(5) + ((t473 + Ifges(5,5) / 0.2e1) * t195 + t489 + t491) * t164 + t433 * qJD(6) + (-pkin(4) * t20 - qJ(5) * t14 - t101 * t58 - t45 * t60 + t457 * t48) * m(6) + t358 * qJD(5) + t345 * t60 + t346 * t59 + t436 + t449 + t8 * (mrSges(7,1) * t243 + mrSges(7,2) * t246) - t243 * t10 / 0.2e1 - t296 * mrSges(7,3) + t274 * t70 - t26 * t78 - t25 * t79 - pkin(4) * t89 - t101 * t103 + t285 * t417 + t288 * t424 + t292 * t425 + t246 * t427; -t358 * t195 + (t103 - t279) * t165 + (-t163 * t295 - t195 * t36 + t296) * m(7) + (t165 * t58 + t195 * t48 + t20) * m(6) - t439; -t36 * (mrSges(7,1) * t119 + mrSges(7,2) * t118) + (Ifges(7,1) * t118 - t374) * t411 + t43 * t410 + (Ifges(7,5) * t118 - Ifges(7,6) * t119) * t408 - t5 * t78 + t6 * t79 + (t118 * t5 + t119 * t6) * mrSges(7,3) + t9 + (-Ifges(7,2) * t119 + t113 + t44) * t413 + t442;];
tauc  = t17(:);
