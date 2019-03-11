% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:47:58
% EndTime: 2019-03-09 00:48:36
% DurationCPUTime: 18.98s
% Computational Cost: add. (14211->700), mult. (34924->989), div. (0->0), fcn. (25855->12), ass. (0->324)
t283 = sin(qJ(3));
t288 = cos(qJ(3));
t313 = pkin(3) * t283 - pkin(9) * t288;
t250 = t313 * qJD(3);
t254 = -pkin(3) * t288 - pkin(9) * t283 - pkin(2);
t282 = sin(qJ(4));
t287 = cos(qJ(4));
t336 = qJD(4) * t287;
t337 = qJD(4) * t282;
t339 = qJD(3) * t287;
t143 = t282 * t250 + t254 * t336 + (-t283 * t339 - t288 * t337) * pkin(8);
t284 = sin(qJ(2));
t278 = sin(pkin(6));
t345 = qJD(1) * t278;
t289 = cos(qJ(2));
t350 = t288 * t289;
t202 = (t282 * t284 + t287 * t350) * t345;
t492 = t143 - t202;
t201 = (-t282 * t350 + t284 * t287) * t345;
t351 = t287 * t288;
t271 = pkin(8) * t351;
t302 = pkin(4) * t283 - pkin(10) * t351;
t340 = qJD(3) * t283;
t392 = pkin(8) * t282;
t346 = t287 * t250 + t340 * t392;
t390 = pkin(10) * t283;
t491 = t302 * qJD(3) + (-t271 + (-t254 + t390) * t282) * qJD(4) + t346 - t201;
t338 = qJD(3) * t288;
t454 = t282 * t338 + t283 * t336;
t490 = -pkin(10) * t454 + t492;
t327 = t284 * t345;
t252 = qJD(2) * pkin(8) + t327;
t279 = cos(pkin(6));
t344 = qJD(1) * t279;
t206 = -t283 * t252 + t288 * t344;
t247 = t313 * qJD(2);
t149 = -t206 * t282 + t287 * t247;
t417 = -pkin(10) - pkin(9);
t328 = qJD(4) * t417;
t489 = -qJD(2) * t302 + t287 * t328 - t149;
t150 = t287 * t206 + t282 * t247;
t341 = qJD(2) * t288;
t323 = t282 * t341;
t488 = -pkin(10) * t323 - t282 * t328 + t150;
t343 = qJD(2) * t283;
t239 = -t282 * t343 + t339;
t240 = qJD(3) * t282 + t287 * t343;
t281 = sin(qJ(5));
t286 = cos(qJ(5));
t173 = t239 * t281 + t240 * t286;
t280 = sin(qJ(6));
t285 = cos(qJ(6));
t314 = t286 * t239 - t240 * t281;
t109 = t173 * t285 + t280 * t314;
t445 = pkin(11) * t314;
t207 = t288 * t252 + t283 * t344;
t193 = qJD(3) * pkin(9) + t207;
t326 = t289 * t345;
t209 = qJD(2) * t254 - t326;
t133 = t193 * t287 + t209 * t282;
t113 = pkin(10) * t239 + t133;
t103 = t286 * t113;
t132 = -t193 * t282 + t287 * t209;
t112 = -pkin(10) * t240 + t132;
t269 = qJD(4) - t341;
t96 = pkin(4) * t269 + t112;
t57 = t281 * t96 + t103;
t42 = t57 + t445;
t367 = t280 * t42;
t260 = qJD(5) + t269;
t462 = pkin(11) * t173;
t101 = t281 * t113;
t56 = t286 * t96 - t101;
t41 = t56 - t462;
t38 = pkin(5) * t260 + t41;
t16 = t285 * t38 - t367;
t365 = t285 * t42;
t17 = t280 * t38 + t365;
t331 = qJD(2) * qJD(3);
t318 = t283 * t331;
t266 = Ifges(7,3) * t318;
t380 = Ifges(7,4) * t109;
t255 = qJD(6) + t260;
t401 = -t255 / 0.2e1;
t411 = t109 / 0.2e1;
t412 = -t109 / 0.2e1;
t451 = -t173 * t280 + t285 * t314;
t414 = -t451 / 0.2e1;
t330 = qJD(3) * qJD(4);
t199 = t287 * t330 + (-t283 * t337 + t287 * t338) * qJD(2);
t355 = t278 * t289;
t324 = qJD(2) * t355;
t295 = qJD(1) * (qJD(3) * t279 + t324);
t159 = -t252 * t340 + t288 * t295;
t203 = (t250 + t327) * qJD(2);
t69 = -qJD(4) * t133 - t159 * t282 + t287 * t203;
t55 = pkin(4) * t318 - pkin(10) * t199 + t69;
t200 = -qJD(2) * t454 - t282 * t330;
t68 = t287 * t159 - t193 * t337 + t282 * t203 + t209 * t336;
t60 = pkin(10) * t200 + t68;
t15 = -qJD(5) * t57 - t281 * t60 + t286 * t55;
t80 = qJD(5) * t314 + t199 * t286 + t200 * t281;
t8 = pkin(5) * t318 - pkin(11) * t80 + t15;
t334 = qJD(5) * t286;
t335 = qJD(5) * t281;
t14 = -t113 * t335 + t281 * t55 + t286 * t60 + t96 * t334;
t81 = -qJD(5) * t173 - t199 * t281 + t200 * t286;
t9 = pkin(11) * t81 + t14;
t2 = qJD(6) * t16 + t280 * t8 + t285 * t9;
t3 = -qJD(6) * t17 - t280 * t9 + t285 * t8;
t32 = qJD(6) * t451 + t280 * t81 + t285 * t80;
t33 = -qJD(6) * t109 - t280 * t80 + t285 * t81;
t431 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t32 + Ifges(7,6) * t33;
t53 = Ifges(7,2) * t451 + Ifges(7,6) * t255 + t380;
t98 = Ifges(7,4) * t451;
t54 = Ifges(7,1) * t109 + Ifges(7,5) * t255 + t98;
t192 = -qJD(3) * pkin(3) - t206;
t152 = -pkin(4) * t239 + t192;
t97 = -pkin(5) * t314 + t152;
t487 = (Ifges(7,1) * t451 - t380) * t412 + (Ifges(7,5) * t451 - Ifges(7,6) * t109) * t401 + (t109 * t17 + t16 * t451) * mrSges(7,3) - t97 * (mrSges(7,1) * t109 + mrSges(7,2) * t451) + t53 * t411 + t266 + t431 + (-Ifges(7,2) * t109 + t54 + t98) * t414;
t238 = t287 * t254;
t176 = -t287 * t390 + t238 + (-pkin(4) - t392) * t288;
t211 = t282 * t254 + t271;
t352 = t282 * t283;
t187 = -pkin(10) * t352 + t211;
t468 = t176 * t334 - t187 * t335 + t281 * t491 + t490 * t286;
t117 = t281 * t176 + t286 * t187;
t467 = -qJD(5) * t117 - t490 * t281 + t286 * t491;
t303 = t281 * t282 - t286 * t287;
t432 = qJD(4) + qJD(5);
t180 = t432 * t303;
t297 = t288 * t303;
t213 = qJD(2) * t297;
t486 = t180 - t213;
t242 = t281 * t287 + t282 * t286;
t181 = t432 * t242;
t298 = t288 * t242;
t212 = qJD(2) * t298;
t477 = t181 - t212;
t167 = Ifges(6,4) * t314;
t267 = Ifges(6,3) * t318;
t381 = Ifges(6,4) * t173;
t399 = -t260 / 0.2e1;
t408 = -t173 / 0.2e1;
t410 = -t314 / 0.2e1;
t429 = t15 * mrSges(6,1) - t14 * mrSges(6,2) + Ifges(6,5) * t80 + Ifges(6,6) * t81;
t95 = t173 * Ifges(6,1) + t260 * Ifges(6,5) + t167;
t485 = t267 + t429 + (Ifges(6,5) * t314 - Ifges(6,6) * t173) * t399 + (t173 * t57 + t314 * t56) * mrSges(6,3) + (-Ifges(6,2) * t173 + t167 + t95) * t410 - t152 * (mrSges(6,1) * t173 + mrSges(6,2) * t314) + (Ifges(6,1) * t314 - t381) * t408 + t487;
t123 = -qJD(3) * t297 - t181 * t283;
t483 = pkin(5) * t340 - pkin(11) * t123 + t467;
t216 = t303 * t283;
t124 = -qJD(3) * t298 + t216 * t432;
t482 = -pkin(11) * t124 - t468;
t258 = t417 * t282;
t259 = t417 * t287;
t198 = t281 * t258 - t286 * t259;
t457 = -qJD(5) * t198 + t488 * t281 + t286 * t489;
t456 = t258 * t334 + t259 * t335 + t281 * t489 - t488 * t286;
t480 = -pkin(5) * t343 + pkin(11) * t486 + t457;
t479 = pkin(11) * t477 - t456;
t474 = -t341 / 0.2e1;
t116 = t286 * t176 - t187 * t281;
t84 = -pkin(5) * t288 + pkin(11) * t216 + t116;
t215 = t242 * t283;
t89 = -pkin(11) * t215 + t117;
t44 = -t280 * t89 + t285 * t84;
t473 = qJD(6) * t44 + t280 * t483 - t285 * t482;
t45 = t280 * t84 + t285 * t89;
t472 = -qJD(6) * t45 + t280 * t482 + t285 * t483;
t464 = -Ifges(4,1) / 0.2e1;
t397 = -t269 / 0.2e1;
t463 = Ifges(4,4) * t474;
t197 = t286 * t258 + t259 * t281;
t153 = -pkin(11) * t242 + t197;
t154 = -pkin(11) * t303 + t198;
t88 = t153 * t280 + t154 * t285;
t461 = -qJD(6) * t88 + t479 * t280 + t480 * t285;
t87 = t153 * t285 - t154 * t280;
t460 = qJD(6) * t87 + t480 * t280 - t479 * t285;
t273 = pkin(4) * t286 + pkin(5);
t332 = qJD(6) * t285;
t333 = qJD(6) * t280;
t353 = t281 * t285;
t62 = -t112 * t281 - t103;
t46 = t62 - t445;
t63 = t286 * t112 - t101;
t47 = t63 - t462;
t459 = t280 * t47 - t285 * t46 - t273 * t333 + (-t281 * t332 + (-t280 * t286 - t353) * qJD(5)) * pkin(4);
t354 = t280 * t281;
t458 = -t280 * t46 - t285 * t47 + t273 * t332 + (-t281 * t333 + (t285 * t286 - t354) * qJD(5)) * pkin(4);
t177 = pkin(4) * t323 + t207;
t455 = pkin(4) * t337 + t477 * pkin(5) - t177;
t253 = -qJD(2) * pkin(2) - t326;
t373 = t240 * Ifges(5,4);
t156 = t239 * Ifges(5,2) + t269 * Ifges(5,6) + t373;
t235 = Ifges(5,4) * t239;
t157 = t240 * Ifges(5,1) + t269 * Ifges(5,5) + t235;
t304 = t132 * t287 + t133 * t282;
t306 = Ifges(5,5) * t287 - Ifges(5,6) * t282;
t382 = Ifges(5,4) * t287;
t308 = -Ifges(5,2) * t282 + t382;
t383 = Ifges(5,4) * t282;
t310 = Ifges(5,1) * t287 - t383;
t311 = mrSges(5,1) * t282 + mrSges(5,2) * t287;
t394 = t287 / 0.2e1;
t395 = -t282 / 0.2e1;
t396 = t269 / 0.2e1;
t402 = t240 / 0.2e1;
t292 = -t304 * mrSges(5,3) + t156 * t395 + t157 * t394 + t192 * t311 + t239 * t308 / 0.2e1 + t310 * t402 + t306 * t396;
t361 = Ifges(4,5) * qJD(3);
t453 = t206 * mrSges(4,3) + t343 * t464 + t463 - t361 / 0.2e1 - t253 * mrSges(4,2) - t292;
t425 = t32 / 0.2e1;
t424 = t33 / 0.2e1;
t419 = t80 / 0.2e1;
t418 = t81 / 0.2e1;
t94 = Ifges(6,2) * t314 + t260 * Ifges(6,6) + t381;
t448 = t94 / 0.2e1;
t447 = t318 / 0.2e1;
t446 = Ifges(5,3) * t397;
t444 = qJD(2) / 0.2e1;
t61 = -mrSges(7,1) * t451 + mrSges(7,2) * t109;
t435 = -m(7) * t97 - t61;
t114 = -mrSges(6,1) * t314 + mrSges(6,2) * t173;
t316 = m(6) * t152 + t114;
t433 = -t282 * t69 + t287 * t68;
t428 = -t69 * mrSges(5,1) + t68 * mrSges(5,2) - Ifges(5,5) * t199 - Ifges(5,6) * t200;
t427 = Ifges(7,4) * t425 + Ifges(7,2) * t424 + Ifges(7,6) * t447;
t426 = Ifges(7,1) * t425 + Ifges(7,4) * t424 + Ifges(7,5) * t447;
t423 = Ifges(6,4) * t419 + Ifges(6,2) * t418 + Ifges(6,6) * t447;
t422 = Ifges(6,1) * t419 + Ifges(6,4) * t418 + Ifges(6,5) * t447;
t413 = t451 / 0.2e1;
t409 = t314 / 0.2e1;
t407 = t173 / 0.2e1;
t406 = t199 / 0.2e1;
t405 = t200 / 0.2e1;
t404 = -t239 / 0.2e1;
t403 = -t240 / 0.2e1;
t400 = t255 / 0.2e1;
t398 = t260 / 0.2e1;
t391 = pkin(8) * t288;
t384 = Ifges(4,4) * t283;
t141 = -t212 * t285 + t213 * t280;
t175 = t242 * t285 - t280 * t303;
t75 = -qJD(6) * t175 + t180 * t280 - t181 * t285;
t363 = t141 - t75;
t142 = -t212 * t280 - t213 * t285;
t174 = -t242 * t280 - t285 * t303;
t74 = qJD(6) * t174 - t180 * t285 - t181 * t280;
t362 = t142 - t74;
t360 = Ifges(4,6) * qJD(3);
t160 = t252 * t338 + t283 * t295;
t356 = t278 * t284;
t220 = -t279 * t288 + t283 * t356;
t359 = t160 * t220;
t347 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t239 + mrSges(5,2) * t240 + mrSges(4,3) * t343;
t251 = pkin(4) * t352 + t283 * pkin(8);
t342 = qJD(2) * t284;
t329 = mrSges(4,3) * t341;
t208 = pkin(4) * t454 + pkin(8) * t338;
t274 = -pkin(4) * t287 - pkin(3);
t325 = t278 * t342;
t320 = t361 / 0.2e1;
t319 = -t360 / 0.2e1;
t317 = 0.3e1 / 0.2e1 * t341;
t312 = mrSges(5,1) * t287 - mrSges(5,2) * t282;
t309 = Ifges(5,1) * t282 + t382;
t307 = Ifges(5,2) * t287 + t383;
t305 = Ifges(5,5) * t282 + Ifges(5,6) * t287;
t221 = t279 * t283 + t288 * t356;
t185 = -t221 * t282 - t287 * t355;
t300 = -t221 * t287 + t282 * t355;
t119 = t185 * t286 + t281 * t300;
t120 = t185 * t281 - t286 * t300;
t65 = t119 * t285 - t120 * t280;
t66 = t119 * t280 + t120 * t285;
t147 = -t215 * t285 + t216 * t280;
t148 = -t215 * t280 - t216 * t285;
t299 = -m(4) * t206 + m(5) * t192 + t347;
t122 = -pkin(4) * t200 + t160;
t291 = t133 * mrSges(5,2) + t17 * mrSges(7,2) + t207 * mrSges(4,3) + t57 * mrSges(6,2) + t446 - t240 * Ifges(5,5) - t239 * Ifges(5,6) + t360 / 0.2e1 + (t288 * Ifges(4,2) + t384) * t444 - t255 * Ifges(7,3) - t109 * Ifges(7,5) - t451 * Ifges(7,6) - t260 * Ifges(6,3) - t173 * Ifges(6,5) - t314 * Ifges(6,6) - t132 * mrSges(5,1) - t16 * mrSges(7,1) - t253 * mrSges(4,1) - t56 * mrSges(6,1);
t290 = qJD(2) ^ 2;
t268 = Ifges(5,3) * t318;
t257 = -qJD(3) * mrSges(4,2) + t329;
t246 = (-mrSges(4,1) * t288 + mrSges(4,2) * t283) * qJD(2);
t231 = (mrSges(4,1) * t283 + mrSges(4,2) * t288) * t331;
t218 = pkin(4) * t353 + t273 * t280;
t217 = -pkin(4) * t354 + t273 * t285;
t214 = pkin(5) * t303 + t274;
t210 = -t282 * t391 + t238;
t205 = mrSges(5,1) * t269 - mrSges(5,3) * t240;
t204 = -mrSges(5,2) * t269 + mrSges(5,3) * t239;
t184 = -qJD(3) * t220 + t288 * t324;
t183 = qJD(3) * t221 + t283 * t324;
t178 = pkin(5) * t215 + t251;
t169 = -mrSges(5,2) * t318 + mrSges(5,3) * t200;
t168 = mrSges(5,1) * t318 - mrSges(5,3) * t199;
t146 = mrSges(6,1) * t260 - mrSges(6,3) * t173;
t145 = -mrSges(6,2) * t260 + mrSges(6,3) * t314;
t144 = -qJD(4) * t211 + t346;
t140 = pkin(4) * t240 + pkin(5) * t173;
t138 = -mrSges(5,1) * t200 + mrSges(5,2) * t199;
t126 = t199 * Ifges(5,1) + t200 * Ifges(5,4) + Ifges(5,5) * t318;
t125 = t199 * Ifges(5,4) + t200 * Ifges(5,2) + Ifges(5,6) * t318;
t100 = qJD(4) * t185 + t184 * t287 + t282 * t325;
t99 = qJD(4) * t300 - t184 * t282 + t287 * t325;
t90 = -pkin(5) * t124 + t208;
t86 = mrSges(7,1) * t255 - mrSges(7,3) * t109;
t85 = -mrSges(7,2) * t255 + mrSges(7,3) * t451;
t77 = -mrSges(6,2) * t318 + mrSges(6,3) * t81;
t76 = mrSges(6,1) * t318 - mrSges(6,3) * t80;
t58 = -pkin(5) * t81 + t122;
t49 = -qJD(6) * t148 - t123 * t280 + t124 * t285;
t48 = qJD(6) * t147 + t123 * t285 + t124 * t280;
t43 = -mrSges(6,1) * t81 + mrSges(6,2) * t80;
t35 = -qJD(5) * t120 - t100 * t281 + t286 * t99;
t34 = qJD(5) * t119 + t100 * t286 + t281 * t99;
t28 = -mrSges(7,2) * t318 + mrSges(7,3) * t33;
t27 = mrSges(7,1) * t318 - mrSges(7,3) * t32;
t19 = t285 * t41 - t367;
t18 = -t280 * t41 - t365;
t12 = -mrSges(7,1) * t33 + mrSges(7,2) * t32;
t7 = -qJD(6) * t66 - t280 * t34 + t285 * t35;
t6 = qJD(6) * t65 + t280 * t35 + t285 * t34;
t1 = [-t221 * mrSges(4,3) * t318 + t65 * t27 + t66 * t28 + t6 * t85 + t7 * t86 + t119 * t76 + t120 * t77 + t34 * t145 + t35 * t146 + t185 * t168 - t300 * t169 + t100 * t204 + t99 * t205 + t184 * t257 + ((-mrSges(3,2) * t290 - t231) * t289 + (-mrSges(3,1) * t290 + qJD(2) * t246) * t284) * t278 + (qJD(3) * t329 + t12 + t138 + t43) * t220 + (t61 + t114 + t347) * t183 + m(7) * (t16 * t7 + t17 * t6 + t183 * t97 + t2 * t66 + t220 * t58 + t3 * t65) + m(6) * (t119 * t15 + t120 * t14 + t122 * t220 + t152 * t183 + t34 * t57 + t35 * t56) + m(5) * (t100 * t133 + t132 * t99 + t183 * t192 + t185 * t69 - t300 * t68 + t359) + m(4) * (t159 * t221 + t359 - t183 * t206 + t184 * t207 + (t253 - t326) * t325); t472 * t86 + (t16 * t472 + t17 * t473 + t178 * t58 + t2 * t45 + t3 * t44 + t90 * t97) * m(7) + t473 * t85 + (t310 * t406 + t308 * t405 + t125 * t395 + t126 * t394 + (-t282 * t68 - t287 * t69) * mrSges(5,3) + (mrSges(4,3) + t311) * t160 + (mrSges(4,2) * t342 + (-t299 - t316 + t435) * t289) * t345 + (t157 * t395 - t287 * t156 / 0.2e1 + t192 * t312 + t307 * t404 + t309 * t403 + t305 * t397 + (t132 * t282 - t133 * t287) * mrSges(5,3)) * qJD(4) + (Ifges(4,1) * t317 - 0.3e1 / 0.2e1 * Ifges(4,2) * t341 + t319 - 0.3e1 / 0.2e1 * Ifges(4,4) * t343 - t291 + (t474 + t396) * Ifges(5,3) + (-Ifges(6,5) * t216 + Ifges(7,5) * t148 - Ifges(6,6) * t215 + Ifges(7,6) * t147 + t306 * t283 + (-Ifges(6,3) - Ifges(7,3)) * t288) * t444) * qJD(3) + (t138 + (m(4) + m(5)) * t160 + (-m(4) * t207 - t257) * qJD(3)) * pkin(8)) * t283 + t467 * t146 + (t116 * t15 + t117 * t14 + t122 * t251 + t152 * t208 + t467 * t56 + t468 * t57) * m(6) + t468 * t145 + m(5) * (t132 * t144 + t133 * t143 + t210 * t69 + t211 * t68) + (t147 * t2 - t148 * t3 - t16 * t48 + t17 * t49) * mrSges(7,3) + (t144 - t201) * t205 + t492 * t204 - m(5) * (t132 * t201 + t133 * t202) + t124 * t448 + (Ifges(7,1) * t148 + Ifges(7,4) * t147) * t425 + t148 * t426 + t147 * t427 + (-t123 * t56 + t124 * t57 - t14 * t215 + t15 * t216) * mrSges(6,3) + (-Ifges(6,4) * t216 - Ifges(6,2) * t215) * t418 + (-Ifges(6,1) * t216 - Ifges(6,4) * t215) * t419 + t122 * (mrSges(6,1) * t215 - mrSges(6,2) * t216) + (t159 * mrSges(4,3) - t266 / 0.2e1 - t267 / 0.2e1 - t268 / 0.2e1 + t428 - t429 - t431) * t288 + (Ifges(6,5) * t123 + Ifges(6,6) * t124) * t398 + (Ifges(7,5) * t48 + Ifges(7,6) * t49) * t400 + (Ifges(6,1) * t123 + Ifges(6,4) * t124) * t407 + (Ifges(6,4) * t123 + Ifges(6,2) * t124) * t409 + (Ifges(7,1) * t48 + Ifges(7,4) * t49) * t411 + (Ifges(7,4) * t48 + Ifges(7,2) * t49) * t413 - t216 * t422 - t215 * t423 + (Ifges(7,4) * t148 + Ifges(7,2) * t147) * t424 + t44 * t27 + t45 * t28 + t49 * t53 / 0.2e1 + t48 * t54 / 0.2e1 + t90 * t61 + t97 * (-mrSges(7,1) * t49 + mrSges(7,2) * t48) + t116 * t76 + t117 * t77 + t123 * t95 / 0.2e1 + t58 * (-mrSges(7,1) * t147 + mrSges(7,2) * t148) + t152 * (-mrSges(6,1) * t124 + mrSges(6,2) * t123) - m(4) * (t207 * t288 * t326 + t253 * t327) + t178 * t12 + t208 * t114 + t210 * t168 + t211 * t169 - pkin(2) * t231 + (-t257 * t350 + (-mrSges(4,1) * t341 - t246) * t284) * t345 + t251 * t43 + m(4) * (-pkin(2) * qJD(1) * t325 + t159 * t391) + (Ifges(4,4) * t317 + t299 * pkin(8) + t320 - t453) * t338; (-Ifges(6,1) * t213 - Ifges(6,4) * t212) * t408 + (-Ifges(6,4) * t213 - Ifges(6,2) * t212) * t410 + (-t180 / 0.2e1 + t213 / 0.2e1) * t95 + (-Ifges(6,5) * t213 - Ifges(6,6) * t212) * t399 + t122 * (mrSges(6,1) * t303 + mrSges(6,2) * t242) - t303 * t423 + (-Ifges(6,5) * t180 - Ifges(6,6) * t181) * t398 + (-Ifges(6,1) * t180 - Ifges(6,4) * t181) * t407 + (-Ifges(6,4) * t180 - Ifges(6,2) * t181) * t409 + (-t181 / 0.2e1 + t212 / 0.2e1) * t94 + (Ifges(6,4) * t242 - Ifges(6,2) * t303) * t418 + (Ifges(6,1) * t242 - Ifges(6,4) * t303) * t419 + (t75 / 0.2e1 - t141 / 0.2e1) * t53 + (-t168 * t282 + t169 * t287) * pkin(9) - m(5) * (t132 * t149 + t133 * t150 + t192 * t207) + (-t14 * t303 - t15 * t242 - t477 * t57 + t486 * t56) * mrSges(6,3) + (mrSges(6,1) * t477 - mrSges(6,2) * t486) * t152 + ((t320 + t463 + t453) * t288 + ((t384 / 0.2e1 + (t464 + Ifges(4,2) / 0.2e1) * t288) * qJD(2) + t319 + t446 + t291) * t283 + (Ifges(6,5) * t242 + Ifges(7,5) * t175 - Ifges(6,6) * t303 + Ifges(7,6) * t174 + t305) * t340 / 0.2e1) * qJD(2) + t460 * t85 + t461 * t86 + (t16 * t461 + t17 * t460 + t2 * t88 + t214 * t58 + t3 * t87 + t455 * t97) * m(7) + t456 * t145 + t457 * t146 + (t122 * t274 + t14 * t198 + t15 * t197 - t152 * t177 + t456 * t57 + t457 * t56) * m(6) + t455 * t61 + t175 * t426 + t174 * t427 + (-mrSges(4,1) - t312) * t160 + t125 * t394 + (Ifges(7,5) * t74 + Ifges(7,6) * t75) * t400 + (Ifges(7,5) * t142 + Ifges(7,6) * t141) * t401 + t307 * t405 + t309 * t406 + (Ifges(7,1) * t74 + Ifges(7,4) * t75) * t411 + (Ifges(7,1) * t142 + Ifges(7,4) * t141) * t412 + (Ifges(7,4) * t74 + Ifges(7,2) * t75) * t413 + (Ifges(7,4) * t142 + Ifges(7,2) * t141) * t414 + t242 * t422 + (Ifges(7,4) * t175 + Ifges(7,2) * t174) * t424 + (Ifges(7,1) * t175 + Ifges(7,4) * t174) * t425 + m(5) * (-pkin(3) * t160 + pkin(9) * t433) + t433 * mrSges(5,3) + t87 * t27 + t88 * t28 + (t74 / 0.2e1 - t142 / 0.2e1) * t54 + (t316 * t282 * pkin(4) + (-m(5) * t304 - t282 * t204 - t287 * t205) * pkin(9) + t292) * qJD(4) - pkin(3) * t138 - t159 * mrSges(4,2) + t58 * (-mrSges(7,1) * t174 + mrSges(7,2) * t175) - t177 * t114 + t197 * t76 + t198 * t77 - t150 * t204 - t149 * t205 + t214 * t12 - t347 * t207 - t206 * t257 + t274 * t43 + (mrSges(7,1) * t363 - mrSges(7,2) * t362) * t97 + (t16 * t362 - t17 * t363 + t174 * t2 - t175 * t3) * mrSges(7,3) + t282 * t126 / 0.2e1; t173 * t448 + (-Ifges(5,2) * t240 + t157 + t235) * t404 + t268 - m(6) * (t56 * t62 + t57 * t63) + (t132 * t239 + t133 * t240) * mrSges(5,3) + t485 + t458 * t85 + (-t140 * t97 + t16 * t459 + t17 * t458 + t2 * t218 + t217 * t3) * m(7) + t459 * t86 - t428 + (Ifges(5,5) * t239 - Ifges(5,6) * t240) * t397 + t156 * t402 + (Ifges(5,1) * t239 - t373) * t403 + (t281 * t77 + t286 * t76 + (t145 * t286 - t146 * t281) * qJD(5) + m(6) * (t14 * t281 + t15 * t286 + t334 * t57 - t335 * t56) - t316 * t240) * pkin(4) - t140 * t61 - t63 * t145 - t62 * t146 - t132 * t204 + t133 * t205 + t217 * t27 + t218 * t28 - t192 * (mrSges(5,1) * t240 + mrSges(5,2) * t239); -m(7) * (t16 * t18 + t17 * t19) + t94 * t407 + (t285 * t27 + t280 * t28 + (-t280 * t86 + t285 * t85) * qJD(6) + m(7) * (-t16 * t333 + t17 * t332 + t2 * t280 + t285 * t3) + t435 * t173) * pkin(5) - t19 * t85 - t18 * t86 - t56 * t145 + t57 * t146 + t485; -t16 * t85 + t17 * t86 + t487;];
tauc  = t1(:);
