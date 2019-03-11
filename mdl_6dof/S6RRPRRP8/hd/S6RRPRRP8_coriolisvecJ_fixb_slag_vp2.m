% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:22:00
% EndTime: 2019-03-09 12:22:37
% DurationCPUTime: 19.12s
% Computational Cost: add. (13717->705), mult. (34486->934), div. (0->0), fcn. (25311->8), ass. (0->307)
t282 = cos(pkin(10));
t281 = sin(pkin(10));
t285 = sin(qJ(2));
t345 = qJD(1) * t285;
t324 = t281 * t345;
t238 = t282 * qJD(2) - t324;
t323 = t282 * t345;
t337 = t281 * qJD(2);
t239 = t323 + t337;
t284 = sin(qJ(4));
t286 = cos(qJ(4));
t179 = t238 * t286 - t239 * t284;
t180 = t238 * t284 + t239 * t286;
t283 = sin(qJ(5));
t386 = cos(qJ(5));
t112 = -t386 * t179 + t283 * t180;
t287 = cos(qJ(2));
t309 = t281 * t284 - t282 * t286;
t296 = t309 * t287;
t293 = qJD(2) * t296;
t132 = -qJD(1) * t293 + qJD(4) * t179;
t244 = t281 * t286 + t282 * t284;
t297 = t244 * t287;
t294 = qJD(2) * t297;
t133 = -qJD(1) * t294 - qJD(4) * t180;
t53 = -qJD(5) * t112 + t132 * t386 + t283 * t133;
t414 = t53 / 0.2e1;
t430 = Ifges(6,1) + Ifges(7,1);
t451 = t414 * t430;
t208 = qJD(1) * t297;
t226 = t244 * qJD(4);
t347 = -t208 + t226;
t437 = t179 * t283 + t180 * t386;
t54 = qJD(5) * t437 + t283 * t132 - t133 * t386;
t412 = t54 / 0.2e1;
t429 = Ifges(7,4) + Ifges(6,5);
t428 = Ifges(7,5) - Ifges(6,4);
t225 = t309 * qJD(4);
t301 = -t283 * t244 - t309 * t386;
t116 = qJD(5) * t301 - t225 * t386 - t283 * t226;
t209 = qJD(1) * t296;
t142 = -t283 * t208 - t209 * t386;
t349 = t116 - t142;
t182 = t244 * t386 - t283 * t309;
t117 = qJD(5) * t182 - t283 * t225 + t226 * t386;
t141 = t208 * t386 - t209 * t283;
t348 = t117 - t141;
t336 = qJD(1) * qJD(2);
t318 = t285 * t336;
t450 = t428 * t412 + t451 + t429 * t318 / 0.2e1;
t427 = -Ifges(6,6) + Ifges(7,6);
t449 = Ifges(6,3) + Ifges(7,2);
t379 = pkin(8) + qJ(3);
t251 = t379 * t281;
t252 = t379 * t282;
t339 = qJD(4) * t286;
t341 = qJD(3) * t282;
t342 = qJD(3) * t281;
t146 = -t251 * t339 + t286 * t341 + (-qJD(4) * t252 - t342) * t284;
t125 = -pkin(9) * t226 + t146;
t196 = -t284 * t251 + t286 * t252;
t147 = -t244 * qJD(3) - qJD(4) * t196;
t290 = pkin(9) * t225 + t147;
t311 = pkin(2) * t285 - qJ(3) * t287;
t246 = t311 * qJD(1);
t198 = pkin(7) * t324 + t282 * t246;
t351 = t282 * t287;
t308 = pkin(3) * t285 - pkin(8) * t351;
t164 = qJD(1) * t308 + t198;
t227 = t281 * t246;
t352 = t282 * t285;
t353 = t281 * t287;
t298 = -pkin(7) * t352 - pkin(8) * t353;
t183 = qJD(1) * t298 + t227;
t105 = t286 * t164 - t183 * t284;
t85 = pkin(4) * t345 + pkin(9) * t209 + t105;
t106 = t284 * t164 + t286 * t183;
t89 = -pkin(9) * t208 + t106;
t195 = -t286 * t251 - t252 * t284;
t157 = -pkin(9) * t244 + t195;
t158 = -pkin(9) * t309 + t196;
t96 = t283 * t157 + t158 * t386;
t445 = -qJD(5) * t96 + (t290 - t85) * t386 + (-t125 + t89) * t283;
t344 = qJD(1) * t287;
t278 = pkin(7) * t344;
t332 = pkin(3) * t344;
t232 = t281 * t332 + t278;
t439 = pkin(4) * t347 - t232;
t448 = Ifges(3,5) / 0.2e1;
t447 = Ifges(4,3) / 0.2e1;
t444 = pkin(5) * t345 - t445;
t443 = pkin(5) * t348 - qJ(6) * t349 - qJD(6) * t182 + t439;
t109 = Ifges(6,4) * t112;
t272 = qJD(4) - t344;
t263 = qJD(5) + t272;
t364 = Ifges(7,5) * t112;
t426 = t263 * t429 + t430 * t437 - t109 + t364;
t442 = Ifges(7,3) * t412;
t373 = mrSges(6,3) * t437;
t100 = mrSges(6,1) * t263 - t373;
t101 = -mrSges(7,1) * t263 + mrSges(7,2) * t437;
t350 = t100 - t101;
t441 = -t105 + t147;
t440 = -t106 + t146;
t277 = pkin(7) * t345;
t249 = -qJD(2) * pkin(2) + qJD(3) + t277;
t197 = -t238 * pkin(3) + t249;
t129 = -pkin(4) * t179 + t197;
t250 = -pkin(2) * t287 - t285 * qJ(3) - pkin(1);
t231 = t250 * qJD(1);
t257 = qJD(2) * qJ(3) + t278;
t185 = t282 * t231 - t281 * t257;
t139 = -t239 * pkin(8) + t185 - t332;
t186 = t281 * t231 + t282 * t257;
t143 = pkin(8) * t238 + t186;
t87 = t284 * t139 + t286 * t143;
t76 = pkin(9) * t179 + t87;
t359 = t283 * t76;
t86 = t286 * t139 - t143 * t284;
t75 = -pkin(9) * t180 + t86;
t73 = pkin(4) * t272 + t75;
t21 = t386 * t73 - t359;
t424 = qJD(6) - t21;
t19 = -t263 * pkin(5) + t424;
t37 = t112 * pkin(5) - qJ(6) * t437 + t129;
t438 = t129 * mrSges(6,2) + t19 * mrSges(7,2) - t21 * mrSges(6,3) - t37 * mrSges(7,3) + t426 / 0.2e1;
t70 = pkin(5) * t437 + qJ(6) * t112;
t436 = -Ifges(6,6) / 0.2e1;
t435 = Ifges(7,6) / 0.2e1;
t406 = -t112 / 0.2e1;
t405 = t112 / 0.2e1;
t434 = -t179 / 0.2e1;
t433 = t179 / 0.2e1;
t392 = t263 / 0.2e1;
t402 = t437 / 0.2e1;
t322 = -Ifges(3,6) * qJD(2) / 0.2e1;
t432 = qJD(2) * t448;
t98 = -mrSges(7,2) * t112 + mrSges(7,3) * t263;
t374 = mrSges(6,3) * t112;
t99 = -mrSges(6,2) * t263 - t374;
t377 = t98 + t99;
t237 = t282 * t250;
t184 = -pkin(8) * t352 + t237 + (-pkin(7) * t281 - pkin(3)) * t287;
t204 = pkin(7) * t351 + t281 * t250;
t354 = t281 * t285;
t192 = -pkin(8) * t354 + t204;
t122 = t284 * t184 + t286 * t192;
t423 = t318 * t449 + t427 * t54 + t429 * t53;
t223 = qJD(2) * t311 - t285 * qJD(3);
t210 = t223 * qJD(1);
t248 = (qJD(3) - t277) * qJD(2);
t162 = t282 * t210 - t281 * t248;
t295 = t308 * qJD(2);
t136 = qJD(1) * t295 + t162;
t163 = t281 * t210 + t282 * t248;
t317 = t287 * t336;
t315 = t281 * t317;
t140 = -pkin(8) * t315 + t163;
t340 = qJD(4) * t284;
t31 = t284 * t136 + t139 * t339 + t286 * t140 - t143 * t340;
t32 = -qJD(4) * t87 + t286 * t136 - t140 * t284;
t422 = -t32 * mrSges(5,1) + t31 * mrSges(5,2);
t343 = qJD(2) * t285;
t26 = pkin(4) * t318 - pkin(9) * t132 + t32;
t28 = pkin(9) * t133 + t31;
t319 = qJD(5) * t386;
t338 = qJD(5) * t283;
t5 = t283 * t26 + t386 * t28 + t73 * t319 - t338 * t76;
t2 = qJ(6) * t318 + qJD(6) * t263 + t5;
t326 = t386 * t76;
t22 = t283 * t73 + t326;
t6 = -qJD(5) * t22 + t26 * t386 - t283 * t28;
t3 = -pkin(5) * t318 - t6;
t421 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3);
t276 = Ifges(3,4) * t344;
t363 = Ifges(4,2) * t281;
t370 = Ifges(4,4) * t282;
t312 = -t363 + t370;
t371 = Ifges(4,4) * t281;
t313 = Ifges(4,1) * t282 - t371;
t376 = mrSges(4,2) * t282;
t388 = t282 / 0.2e1;
t389 = -t281 / 0.2e1;
t420 = -(t185 * t282 + t186 * t281) * mrSges(4,3) + t249 * (mrSges(4,1) * t281 + t376) + Ifges(3,1) * t345 / 0.2e1 + t276 / 0.2e1 + t432 + t238 * t312 / 0.2e1 + t239 * t313 / 0.2e1 + (Ifges(4,4) * t239 + Ifges(4,2) * t238 - Ifges(4,6) * t344) * t389 + (Ifges(4,1) * t239 + Ifges(4,4) * t238 - Ifges(4,5) * t344) * t388;
t121 = t286 * t184 - t284 * t192;
t219 = t309 * t285;
t92 = -pkin(4) * t287 + t219 * pkin(9) + t121;
t218 = t244 * t285;
t97 = -pkin(9) * t218 + t122;
t378 = t283 * t92 + t386 * t97;
t159 = -t226 * t285 - t293;
t330 = pkin(7) * t343;
t190 = t282 * t223 + t281 * t330;
t154 = t295 + t190;
t213 = t281 * t223;
t165 = qJD(2) * t298 + t213;
t68 = -qJD(4) * t122 + t286 * t154 - t165 * t284;
t43 = pkin(4) * t343 - pkin(9) * t159 + t68;
t160 = t225 * t285 - t294;
t67 = t284 * t154 + t286 * t165 + t184 * t339 - t192 * t340;
t55 = pkin(9) * t160 + t67;
t10 = -qJD(5) * t378 - t283 * t55 + t386 * t43;
t20 = t263 * qJ(6) + t22;
t327 = t435 + t436;
t328 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t329 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t372 = Ifges(3,4) * t285;
t419 = -t327 * t112 - t328 * t263 - t329 * t437 - t185 * mrSges(4,1) - t20 * mrSges(7,3) - t21 * mrSges(6,1) - t86 * mrSges(5,1) - Ifges(5,6) * t179 - Ifges(5,3) * t272 - Ifges(5,5) * t180 + t344 * t447 - Ifges(4,6) * t238 - Ifges(4,5) * t239 - t322 + (t287 * Ifges(3,2) + t372) * qJD(1) / 0.2e1 - Ifges(6,6) * t406 - Ifges(7,6) * t405 + t186 * mrSges(4,2) + t19 * mrSges(7,1) + t22 * mrSges(6,2) + t87 * mrSges(5,2) - t429 * t402 - t449 * t392;
t108 = Ifges(7,5) * t437;
t61 = Ifges(7,6) * t263 + Ifges(7,3) * t112 + t108;
t368 = Ifges(6,4) * t437;
t64 = -Ifges(6,2) * t112 + Ifges(6,6) * t263 + t368;
t417 = t129 * mrSges(6,1) + t37 * mrSges(7,1) + t61 / 0.2e1 - t64 / 0.2e1 - t20 * mrSges(7,2) - t22 * mrSges(6,3);
t416 = Ifges(7,5) * t414 + t318 * t435 + t442;
t415 = -t53 * Ifges(6,4) / 0.2e1 + Ifges(6,2) * t412 + t318 * t436;
t413 = -t54 / 0.2e1;
t410 = pkin(1) * mrSges(3,1);
t409 = pkin(1) * mrSges(3,2);
t403 = -t437 / 0.2e1;
t401 = t132 / 0.2e1;
t400 = t133 / 0.2e1;
t398 = -t180 / 0.2e1;
t397 = t180 / 0.2e1;
t395 = -t218 / 0.2e1;
t394 = -t219 / 0.2e1;
t393 = -t263 / 0.2e1;
t391 = -t272 / 0.2e1;
t390 = t272 / 0.2e1;
t385 = pkin(4) * t180;
t383 = pkin(4) * t283;
t36 = t283 * t85 + t386 * t89;
t375 = mrSges(5,3) * t180;
t369 = Ifges(5,4) * t180;
t366 = Ifges(4,5) * t282;
t361 = Ifges(4,6) * t281;
t356 = qJD(2) * mrSges(3,2);
t346 = -t209 + t225;
t212 = mrSges(4,1) * t315 + t317 * t376;
t271 = pkin(7) * t317;
t222 = pkin(3) * t315 + t271;
t233 = (pkin(3) * t337 + pkin(7) * qJD(2)) * t287;
t247 = pkin(3) * t354 + t285 * pkin(7);
t333 = t386 * pkin(4);
t325 = Ifges(5,5) * t132 + Ifges(5,6) * t133 + Ifges(5,3) * t318;
t274 = -pkin(3) * t282 - pkin(2);
t17 = t54 * mrSges(6,1) + t53 * mrSges(6,2);
t16 = t54 * mrSges(7,1) - t53 * mrSges(7,3);
t82 = -t133 * mrSges(5,1) + t132 * mrSges(5,2);
t316 = pkin(4) * t319;
t314 = m(4) * t249 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t238 + mrSges(4,2) * t239 + mrSges(3,3) * t345;
t107 = -pkin(4) * t133 + t222;
t134 = -pkin(4) * t160 + t233;
t188 = pkin(4) * t218 + t247;
t211 = pkin(4) * t309 + t274;
t40 = -mrSges(7,1) * t318 + t53 * mrSges(7,2);
t307 = t366 / 0.2e1 - t361 / 0.2e1;
t59 = -t283 * t97 + t386 * t92;
t9 = t283 * t43 + t92 * t319 - t338 * t97 + t386 * t55;
t303 = t157 * t386 - t283 * t158;
t302 = -t218 * t386 + t283 * t219;
t153 = -t283 * t218 - t219 * t386;
t299 = t179 * mrSges(5,3);
t291 = t421 + t423;
t275 = -t333 - pkin(5);
t273 = qJ(6) + t383;
t264 = t316 + qJD(6);
t258 = mrSges(3,3) * t344 - t356;
t221 = (mrSges(4,1) * t285 - mrSges(4,3) * t351) * t336;
t220 = (-mrSges(4,2) * t285 - mrSges(4,3) * t353) * t336;
t207 = -mrSges(4,1) * t344 - t239 * mrSges(4,3);
t206 = mrSges(4,2) * t344 + t238 * mrSges(4,3);
t203 = -pkin(7) * t353 + t237;
t199 = -pkin(7) * t323 + t227;
t194 = (Ifges(4,5) * t285 + t287 * t313) * t336;
t193 = (Ifges(4,6) * t285 + t287 * t312) * t336;
t191 = -t282 * t330 + t213;
t177 = Ifges(5,4) * t179;
t150 = mrSges(5,1) * t272 - t375;
t149 = -t272 * mrSges(5,2) + t299;
t120 = -mrSges(5,2) * t318 + mrSges(5,3) * t133;
t119 = mrSges(5,1) * t318 - mrSges(5,3) * t132;
t118 = -mrSges(5,1) * t179 + t180 * mrSges(5,2);
t104 = t180 * Ifges(5,1) + t272 * Ifges(5,5) + t177;
t103 = Ifges(5,2) * t179 + Ifges(5,6) * t272 + t369;
t94 = -pkin(5) * t301 - qJ(6) * t182 + t211;
t83 = -pkin(5) * t302 - qJ(6) * t153 + t188;
t81 = t132 * Ifges(5,1) + t133 * Ifges(5,4) + Ifges(5,5) * t318;
t80 = t132 * Ifges(5,4) + t133 * Ifges(5,2) + Ifges(5,6) * t318;
t79 = qJD(5) * t153 + t283 * t159 - t160 * t386;
t78 = qJD(5) * t302 + t159 * t386 + t283 * t160;
t72 = mrSges(6,1) * t112 + mrSges(6,2) * t437;
t71 = mrSges(7,1) * t112 - mrSges(7,3) * t437;
t58 = t385 + t70;
t57 = t287 * pkin(5) - t59;
t56 = -qJ(6) * t287 + t378;
t41 = -mrSges(6,2) * t318 - mrSges(6,3) * t54;
t39 = mrSges(6,1) * t318 - mrSges(6,3) * t53;
t38 = -mrSges(7,2) * t54 + mrSges(7,3) * t318;
t33 = qJ(6) * t345 + t36;
t29 = qJD(5) * t303 + t125 * t386 + t283 * t290;
t25 = t386 * t75 - t359;
t24 = t283 * t75 + t326;
t18 = pkin(5) * t79 - qJ(6) * t78 - qJD(6) * t153 + t134;
t11 = pkin(5) * t54 - qJ(6) * t53 - qJD(6) * t437 + t107;
t8 = -pkin(5) * t343 - t10;
t7 = qJ(6) * t343 - qJD(6) * t287 + t9;
t1 = [((-pkin(7) * t258 + t322 - t419) * qJD(2) + t193 * t389 + t194 * t388 + pkin(7) * t212 + (-t162 * t282 - t163 * t281) * mrSges(4,3) + (-0.2e1 * t410 + Ifges(5,5) * t394 + Ifges(5,6) * t395 + (-0.3e1 / 0.2e1 * Ifges(3,4) + t307) * t285 + t329 * t153 - t327 * t302) * t336) * t285 - (mrSges(6,1) * t107 + mrSges(7,1) * t11 - mrSges(7,2) * t2 - mrSges(6,3) * t5 - Ifges(6,2) * t413 + t414 * t428 + t415 + t416 + t442) * t302 + (-t159 * t86 + t160 * t87 - t218 * t31 + t219 * t32) * mrSges(5,3) + (-Ifges(5,4) * t219 - Ifges(5,2) * t218) * t400 + (-Ifges(5,1) * t219 - Ifges(5,4) * t218) * t401 + t222 * (mrSges(5,1) * t218 - mrSges(5,2) * t219) + t378 * t41 + m(6) * (t10 * t21 + t107 * t188 + t129 * t134 + t22 * t9 + t378 * t5 + t59 * t6) + (-Ifges(6,2) * t406 + Ifges(7,3) * t405 + t392 * t427 + t402 * t428 + t417) * t79 + t153 * t450 + (t163 * mrSges(4,2) - t162 * mrSges(4,1) + (t314 * pkin(7) + t420 + t432) * qJD(2) - Ifges(5,6) * t400 - Ifges(5,5) * t401 - Ifges(7,6) * t412 - Ifges(6,6) * t413 - t429 * t414 + (-0.2e1 * t409 + (0.3e1 / 0.2e1 * Ifges(3,4) - 0.3e1 / 0.2e1 * t366 + 0.3e1 / 0.2e1 * t361) * t287 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + Ifges(4,1) * t282 ^ 2 / 0.2e1 - Ifges(5,3) / 0.2e1 + (m(4) * pkin(7) + t376) * pkin(7) + (pkin(7) * mrSges(4,1) - t370 + t363 / 0.2e1) * t281 - t328) * t285) * t336 - t421 + t422) * t287 - (t325 + t423) * t287 / 0.2e1 + (Ifges(6,4) * t406 + Ifges(7,5) * t405 + t429 * t392 + t430 * t402 + t438) * t78 + m(4) * (t162 * t203 + t163 * t204 + t185 * t190 + t186 * t191) + (mrSges(6,2) * t107 + mrSges(7,2) * t3 - mrSges(6,3) * t6 - mrSges(7,3) * t11 + Ifges(6,4) * t413 + Ifges(7,5) * t412 + t451) * t153 + (Ifges(5,5) * t159 + Ifges(5,6) * t160) * t390 + t81 * t394 + t80 * t395 + (Ifges(5,1) * t159 + Ifges(5,4) * t160) * t397 + (Ifges(5,4) * t159 + Ifges(5,2) * t160) * t433 + t56 * t38 + t57 * t40 + t59 * t39 + t18 * t71 + t83 * t16 + t7 * t98 + t9 * t99 + t10 * t100 + t8 * t101 + t121 * t119 + t122 * t120 + t134 * t72 + t67 * t149 + t68 * t150 + t159 * t104 / 0.2e1 + t160 * t103 / 0.2e1 + t188 * t17 + t197 * (-mrSges(5,1) * t160 + mrSges(5,2) * t159) + t191 * t206 + t190 * t207 + t204 * t220 + t203 * t221 + t233 * t118 + t247 * t82 + m(5) * (t121 * t32 + t122 * t31 + t197 * t233 + t222 * t247 + t67 * t87 + t68 * t86) + m(7) * (t11 * t83 + t18 * t37 + t19 * t8 + t2 * t56 + t20 * t7 + t3 * t57); t439 * t72 + ((-t276 / 0.2e1 + (t287 * t307 + t409) * qJD(1) + (t448 + (Ifges(4,1) * t281 + t370) * t388 + (Ifges(4,2) * t282 + t371) * t389) * qJD(2) + ((-m(4) * pkin(2) - mrSges(4,1) * t282 + mrSges(4,2) * t281 - mrSges(3,1)) * qJD(2) - t314) * pkin(7) - t420) * t287 + ((t410 + t372 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + t447 + Ifges(3,2) / 0.2e1) * t287) * qJD(1) + t322 + (t258 + t356) * pkin(7) + t419) * t285 + (Ifges(4,5) * t281 + Ifges(5,5) * t244 + Ifges(4,6) * t282 - Ifges(5,6) * t309 + t182 * t429 - t301 * t427) * t343 / 0.2e1) * qJD(1) + (t208 / 0.2e1 - t226 / 0.2e1) * t103 + (-Ifges(5,5) * t209 - Ifges(5,6) * t208) * t391 + (-Ifges(5,1) * t209 - Ifges(5,4) * t208) * t398 + (-Ifges(5,4) * t209 - Ifges(5,2) * t208) * t434 + (t209 / 0.2e1 - t225 / 0.2e1) * t104 + (t61 - t64) * (t117 / 0.2e1 - t141 / 0.2e1) - (t40 - t39) * t303 + (t182 * t3 + t19 * t349 + t2 * t301 - t20 * t348) * mrSges(7,2) + (-t182 * t6 - t21 * t349 - t22 * t348 + t301 * t5) * mrSges(6,3) + (Ifges(7,5) * t182 - Ifges(7,3) * t301) * t412 + (Ifges(6,4) * t182 + Ifges(6,2) * t301) * t413 + t11 * (-mrSges(7,1) * t301 - mrSges(7,3) * t182) + t107 * (-mrSges(6,1) * t301 + mrSges(6,2) * t182) - t301 * t415 - t301 * t416 + (-t244 * t32 - t309 * t31 + t346 * t86 - t347 * t87) * mrSges(5,3) + (Ifges(5,4) * t244 - Ifges(5,2) * t309) * t400 + (Ifges(5,1) * t244 - Ifges(5,4) * t309) * t401 + t222 * (mrSges(5,1) * t309 + mrSges(5,2) * t244) - t309 * t80 / 0.2e1 + (t116 * t429 + t117 * t427) * t392 + (t141 * t427 + t142 * t429) * t393 + (t141 * t428 + t142 * t430) * t403 + (t116 * t430 + t117 * t428) * t402 + (t182 * t430 - t301 * t428) * t414 + (-Ifges(5,5) * t225 - Ifges(5,6) * t226) * t390 + (-Ifges(5,1) * t225 - Ifges(5,4) * t226) * t397 + m(4) * (-t185 * t342 + t186 * t341 + (-t162 * t281 + t163 * t282) * qJ(3)) + t426 * (t116 / 0.2e1 - t142 / 0.2e1) + t443 * t71 + (t11 * t94 + t2 * t96 - t3 * t303 + t443 * t37 + (t29 - t33) * t20 + t444 * t19) * m(7) + t444 * t101 + (t107 * t211 + t5 * t96 + t6 * t303 + (t29 - t36) * t22 + t445 * t21 + t439 * t129) * m(6) + t445 * t100 + t182 * t450 + t440 * t149 + (t195 * t32 + t196 * t31 - t197 * t232 + t222 * t274 + t440 * t87 + t441 * t86) * m(5) + t441 * t150 + (t38 + t41) * t96 - m(4) * (t185 * t198 + t186 * t199) + (Ifges(6,4) * t116 + Ifges(7,5) * t142 - Ifges(6,2) * t117 + Ifges(7,3) * t141) * t406 + (Ifges(6,4) * t142 + Ifges(7,5) * t116 - Ifges(6,2) * t141 + Ifges(7,3) * t117) * t405 + t377 * t29 + (mrSges(6,1) * t348 + mrSges(6,2) * t349) * t129 + (mrSges(7,1) * t348 - mrSges(7,3) * t349) * t37 + (mrSges(5,1) * t347 - mrSges(5,2) * t346) * t197 + (-Ifges(5,4) * t225 - Ifges(5,2) * t226) * t433 + (-qJ(3) * t221 - qJD(3) * t207 - t162 * mrSges(4,3) + t194 / 0.2e1) * t281 + t94 * t16 - t33 * t98 - t36 * t99 + (qJ(3) * t220 + qJD(3) * t206 + t163 * mrSges(4,3) + t193 / 0.2e1) * t282 + t195 * t119 + t196 * t120 - t199 * t206 - t198 * t207 + t211 * t17 - pkin(2) * t212 - t232 * t118 + t244 * t81 / 0.2e1 + t274 * t82; t350 * t437 + t377 * t112 - t179 * t149 + t180 * t150 - t238 * t206 + t239 * t207 + t16 + t17 + t212 + t82 + (t112 * t20 - t19 * t437 + t11) * m(7) + (t112 * t22 + t21 * t437 + t107) * m(6) + (-t179 * t87 + t180 * t86 + t222) * m(5) + (t185 * t239 - t186 * t238 + t271) * m(4); (t375 + t150) * t87 + (t299 - t149) * t86 + (t2 * t273 + t275 * t3 - t37 * t58 + (-t25 + t264) * t20) * m(7) - t377 * t25 + (-Ifges(6,2) * t405 + Ifges(7,3) * t406 + t393 * t427 + t403 * t428 - t417) * t437 - t422 + t325 + (-t129 * t385 + t21 * t24 - t22 * t25 + (t386 * t6 + t283 * t5 + (-t21 * t283 + t22 * t386) * qJD(5)) * pkin(4)) * m(6) + t291 + t39 * t333 + t99 * t316 - t72 * t385 + t41 * t383 + (Ifges(5,5) * t179 - Ifges(5,6) * t180) * t391 + t103 * t397 + (Ifges(5,1) * t179 - t369) * t398 - t58 * t71 - t197 * (t180 * mrSges(5,1) + mrSges(5,2) * t179) + t264 * t98 + t273 * t38 + t275 * t40 + (-Ifges(5,2) * t180 + t104 + t177) * t434 + (-m(7) * t19 + t350) * (-pkin(4) * t338 + t24) + (-Ifges(6,4) * t405 - Ifges(7,5) * t406 - t429 * t393 - t430 * t403 + t438) * t112; (t112 * t19 + t20 * t437) * mrSges(7,2) + t291 + (-t374 - t377) * t21 + (t350 + t373) * t22 + qJ(6) * t38 - pkin(5) * t40 - t70 * t71 + qJD(6) * t98 - t37 * (mrSges(7,1) * t437 + mrSges(7,3) * t112) + (Ifges(7,3) * t437 - t364) * t406 + t64 * t402 - t129 * (mrSges(6,1) * t437 - mrSges(6,2) * t112) + (-t112 * t429 + t427 * t437) * t393 + (-pkin(5) * t3 + qJ(6) * t2 - t19 * t22 + t20 * t424 - t37 * t70) * m(7) + (-Ifges(6,2) * t437 - t109 + t426) * t405 + (-t112 * t430 + t108 - t368 + t61) * t403; t437 * t71 - t263 * t98 + 0.2e1 * (t3 / 0.2e1 + t37 * t402 + t20 * t393) * m(7) + t40;];
tauc  = t1(:);
