% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRRPP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:58:33
% EndTime: 2019-03-08 22:58:46
% DurationCPUTime: 5.75s
% Computational Cost: add. (4495->539), mult. (10161->675), div. (0->0), fcn. (7789->10), ass. (0->255)
t197 = cos(qJ(4));
t194 = sin(qJ(4));
t269 = -t194 * qJ(5) - pkin(3);
t377 = pkin(4) * t197 - t269;
t198 = cos(qJ(3));
t380 = -pkin(3) * t198 - pkin(2);
t192 = sin(pkin(6));
t196 = sin(qJ(2));
t199 = cos(qJ(2));
t318 = t198 * t199;
t101 = (t194 * t196 + t197 * t318) * t192;
t195 = sin(qJ(3));
t250 = pkin(3) * t195 - pkin(9) * t198;
t145 = t250 * qJD(3);
t150 = -pkin(9) * t195 + t380;
t299 = qJD(4) * t197;
t379 = -qJD(1) * t101 + t194 * t145 + t150 * t299;
t185 = t198 * qJDD(2);
t296 = qJD(2) * qJD(3);
t136 = t195 * t296 + qJDD(4) - t185;
t355 = pkin(4) + qJ(6);
t274 = t355 * t136;
t304 = qJD(3) * t194;
t307 = qJD(2) * t195;
t140 = t197 * t307 + t304;
t309 = qJD(1) * t192;
t146 = qJD(2) * pkin(8) + t196 * t309;
t339 = cos(pkin(6));
t260 = qJD(1) * t339;
t169 = t195 * t260;
t97 = t198 * t146 + t169;
t88 = qJD(3) * pkin(9) + t97;
t275 = t199 * t309;
t99 = qJD(2) * t150 - t275;
t29 = t194 * t88 - t197 * t99;
t237 = t140 * pkin(5) + t29;
t316 = qJD(5) + t237;
t300 = qJD(4) * t195;
t378 = qJD(2) * t300 - qJDD(3);
t295 = t195 * qJDD(2);
t306 = qJD(2) * t198;
t59 = t194 * (qJD(3) * (qJD(4) + t306) + t295) + t378 * t197;
t319 = t197 * t198;
t178 = pkin(8) * t319;
t301 = qJD(4) * t194;
t370 = t194 * t318 - t196 * t197;
t376 = qJD(4) * t178 - t197 * t145 + t150 * t301 - t370 * t309;
t375 = t198 * qJD(5) - t379;
t374 = -t195 * t146 + t198 * t260;
t135 = t140 ^ 2;
t174 = -qJD(4) + t306;
t171 = t174 ^ 2;
t373 = -t135 - t171;
t129 = t136 * qJ(5);
t162 = qJD(5) * t174;
t372 = t162 - t129;
t322 = t194 * t198;
t292 = pkin(4) * t322;
t371 = pkin(4) * t301 - qJD(2) * t292 - t194 * qJD(5) - t97;
t369 = -pkin(5) * t59 + qJDD(6);
t326 = t192 * t198;
t126 = t195 * t339 + t196 * t326;
t308 = qJD(2) * t192;
t281 = t196 * t308;
t325 = t192 * t199;
t286 = t197 * t325;
t327 = t192 * t196;
t125 = t195 * t327 - t198 * t339;
t280 = t199 * t308;
t77 = -qJD(3) * t125 + t198 * t280;
t12 = -qJD(4) * t286 - t126 * t301 + t194 * t281 + t197 * t77;
t272 = t198 * t296;
t298 = t197 * qJD(3);
t58 = -qJD(4) * t298 + (-t272 - t295) * t197 + t378 * t194;
t78 = qJD(3) * t126 + t195 * t280;
t80 = t126 * t197 - t194 * t325;
t226 = -t12 * t174 + t125 * t58 + t136 * t80 - t78 * t140;
t11 = qJD(4) * t80 + t77 * t194 - t197 * t281;
t138 = t194 * t307 - t298;
t79 = t126 * t194 + t286;
t227 = t11 * t174 + t125 * t59 - t136 * t79 + t78 * t138;
t200 = qJD(3) ^ 2;
t297 = qJD(1) * qJD(2);
t273 = t196 * t297;
t243 = -qJDD(1) * t325 + t192 * t273;
t191 = sin(pkin(10));
t338 = cos(pkin(10));
t244 = t339 * t338;
t120 = t191 * t196 - t199 * t244;
t265 = t191 * t339;
t122 = t196 * t338 + t199 * t265;
t249 = g(1) * t122 + g(2) * t120;
t368 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t200 + t192 * (-g(3) * t199 + t273) - t243 + t249;
t366 = t138 ^ 2;
t365 = 0.2e1 * t129;
t364 = pkin(5) + pkin(9);
t362 = t58 * pkin(5);
t359 = pkin(5) * t138;
t358 = pkin(9) * t136;
t357 = t136 * pkin(4);
t7 = t174 * t355 + t316;
t356 = t174 * t7;
t276 = t194 * t300;
t282 = -pkin(8) * t194 - pkin(4);
t290 = pkin(5) * t319;
t354 = -pkin(5) * t276 + t198 * qJD(6) + (t290 + (-qJ(6) + t282) * t195) * qJD(3) + t376;
t177 = pkin(8) * t322;
t291 = pkin(5) * t322;
t321 = t195 * t197;
t353 = (-pkin(5) * t321 - t177) * qJD(4) + (-t291 + (-pkin(8) * t197 + qJ(5)) * t195) * qJD(3) - t375;
t303 = qJD(3) * t195;
t352 = -qJ(5) * t303 + (t195 * t298 + t198 * t301) * pkin(8) + t375;
t351 = t282 * t303 + t376;
t30 = t194 * t99 + t197 * t88;
t337 = qJ(5) * t197;
t239 = qJ(6) * t194 - t337;
t225 = t239 * t198;
t350 = -qJD(2) * t225 + qJD(4) * t239 - t197 * qJD(6) + t371;
t349 = pkin(9) * qJD(4);
t348 = qJ(5) * t59;
t347 = qJD(2) * pkin(2);
t20 = qJ(5) * t174 - t30;
t346 = t174 * t20;
t345 = t174 * t30;
t344 = t58 * t194;
t343 = -qJ(5) * t299 + t306 * t337 + t371;
t142 = t250 * qJD(2);
t340 = t194 * t142 + t197 * t374;
t342 = -t364 * t301 - (qJ(5) * t195 - t291) * qJD(2) - t340;
t157 = t364 * t197;
t84 = t194 * t374;
t261 = -t197 * t142 + t84;
t341 = qJD(4) * t157 - (-t195 * t355 + t290) * qJD(2) - t261;
t336 = qJ(6) * t197;
t334 = t120 * t195;
t333 = t122 * t195;
t332 = t125 * t197;
t331 = t138 * qJ(5);
t330 = t138 * t174;
t329 = t140 * t138;
t328 = t174 * t140;
t323 = t194 * t195;
t317 = -qJD(5) - t29;
t16 = t30 - t359;
t315 = -qJD(6) - t16;
t314 = qJDD(1) - g(3);
t312 = pkin(4) * t323 + t195 * pkin(8);
t311 = t194 * t150 + t178;
t189 = t195 ^ 2;
t310 = -t198 ^ 2 + t189;
t305 = qJD(3) * t140;
t302 = qJD(3) * t198;
t121 = t191 * t199 + t196 * t244;
t264 = t192 * t338;
t72 = -t121 * t195 - t198 * t264;
t294 = t377 * t72;
t123 = -t196 * t265 + t199 * t338;
t74 = -t123 * t195 + t191 * t326;
t293 = t377 * t74;
t87 = -qJD(3) * pkin(3) - t374;
t209 = -t140 * qJ(5) + t87;
t17 = t138 * t355 + t209;
t289 = t17 * t301;
t288 = t17 * t299;
t287 = t195 * t325;
t284 = -pkin(4) * t332 + t125 * t269;
t279 = t174 * t304;
t278 = t174 * t298;
t277 = t174 * t301;
t271 = t199 * t296;
t73 = t121 * t198 - t195 * t264;
t36 = -t120 * t197 + t194 * t73;
t37 = t120 * t194 + t197 * t73;
t268 = -t36 * pkin(4) + qJ(5) * t37;
t75 = t191 * t192 * t195 + t123 * t198;
t38 = -t122 * t197 + t194 * t75;
t39 = t122 * t194 + t197 * t75;
t267 = -t38 * pkin(4) + qJ(5) * t39;
t266 = -t79 * pkin(4) + qJ(5) * t80;
t103 = qJDD(2) * pkin(8) + (qJDD(1) * t196 + t199 * t297) * t192;
t258 = qJDD(1) * t339;
t247 = t195 * t258;
t24 = qJDD(3) * pkin(9) + qJD(3) * t374 + t198 * t103 + t247;
t47 = qJD(2) * t145 + qJDD(2) * t150 + t243;
t263 = t194 * t24 - t197 * t47 + t88 * t299 + t99 * t301;
t262 = -t194 * t47 - t197 * t24 - t99 * t299 + t88 * t301;
t259 = t150 * t197 - t177;
t255 = t195 * pkin(4) * t299 + pkin(8) * t302 + qJ(5) * t276 + qJD(3) * t292;
t254 = t195 * t275;
t253 = t138 * t275;
t252 = t140 * t275;
t91 = qJ(5) * t198 - t311;
t246 = -qJDD(5) - t263;
t19 = pkin(4) * t174 - t317;
t240 = t19 * t197 + t194 * t20;
t4 = t262 + t372;
t201 = qJD(2) ^ 2;
t238 = qJDD(2) * t199 - t196 * t201;
t234 = t194 * t136 - t174 * t299;
t233 = t197 * t136 + t277;
t232 = t11 * t140 - t12 * t138 - t79 * t58 - t59 * t80;
t100 = t370 * t192;
t51 = -t120 * t322 - t121 * t197;
t53 = -t122 * t322 - t123 * t197;
t231 = g(1) * t53 + g(2) * t51 + g(3) * t100;
t52 = -t120 * t319 + t121 * t194;
t54 = -t122 * t319 + t123 * t194;
t230 = -g(1) * t54 - g(2) * t52 - g(3) * t101;
t229 = -g(1) * t74 - g(2) * t72 + g(3) * t125;
t228 = g(1) * t75 + g(2) * t73 + g(3) * t126;
t218 = -qJD(3) * t169 - t195 * t103 - t146 * t302 + t198 * t258;
t25 = -qJDD(3) * pkin(3) - t218;
t203 = t58 * qJ(5) - t140 * qJD(5) + t25;
t3 = t138 * qJD(6) + t355 * t59 + t203;
t224 = t229 - t3;
t223 = t192 * pkin(3) * t318 + pkin(2) * t325 + t101 * pkin(4) + pkin(8) * t327 + pkin(9) * t287 + qJ(5) * t100;
t221 = -g(3) * t325 + t249;
t220 = -t174 * t87 - t358;
t31 = t138 * pkin(4) + t209;
t219 = t174 * t31 + t358;
t217 = t52 * pkin(4) + pkin(8) * t121 - pkin(9) * t334 + qJ(5) * t51 + t120 * t380;
t216 = t54 * pkin(4) + pkin(8) * t123 - pkin(9) * t333 + qJ(5) * t53 + t122 * t380;
t215 = t136 - t329;
t214 = g(1) * t38 + g(2) * t36 + g(3) * t79 - t263;
t213 = g(1) * t39 + g(2) * t37 + g(3) * t80 + t262;
t212 = -t174 * t349 - t229;
t6 = t59 * pkin(4) + t203;
t210 = -t212 - t6;
t208 = -qJDD(5) + t214;
t27 = -t58 - t330;
t147 = -t275 - t347;
t207 = -pkin(8) * qJDD(3) + (t147 + t275 - t347) * qJD(3);
t206 = t140 * t31 - t208;
t205 = t17 * t140 - t208 - t362;
t204 = -t138 * t17 - t213 + t369;
t202 = -t328 - t59;
t188 = t198 * pkin(4);
t156 = t364 * t194;
t130 = -t197 * t355 + t269;
t115 = -qJ(5) * t321 + t312;
t92 = t188 - t259;
t86 = t195 * t239 + t312;
t71 = pkin(4) * t140 + t331;
t65 = -pkin(5) * t323 - t91;
t60 = t198 * qJ(6) + t177 + t188 + (pkin(5) * t195 - t150) * t197;
t48 = t140 * t355 + t331;
t46 = (-qJ(5) * t302 - qJD(5) * t195) * t197 + t255;
t45 = -pkin(4) * t307 + t261;
t44 = -qJ(5) * t307 - t340;
t14 = qJD(3) * t225 + (qJD(6) * t194 + (qJ(6) * qJD(4) - qJD(5)) * t197) * t195 + t255;
t9 = qJD(6) - t20 - t359;
t5 = -t246 - t357;
t2 = -t4 + t369;
t1 = qJD(6) * t174 - t246 - t274 - t362;
t8 = [t314, 0, t238 * t192 (-qJDD(2) * t196 - t199 * t201) * t192, 0, 0, 0, 0, 0, -t78 * qJD(3) - t125 * qJDD(3) + (-t195 * t271 + t198 * t238) * t192, -t77 * qJD(3) - t126 * qJDD(3) + (-t195 * t238 - t198 * t271) * t192, 0, 0, 0, 0, 0, t227, -t226, t232, -t227, t226, t11 * t19 - t12 * t20 + t125 * t6 + t31 * t78 - t4 * t80 + t5 * t79 - g(3), t232, t226, t227, t1 * t79 + t11 * t7 + t12 * t9 + t125 * t3 + t17 * t78 + t2 * t80 - g(3); 0, qJDD(2), t314 * t325 + t249, g(1) * t123 + g(2) * t121 - t314 * t327, qJDD(2) * t189 + 0.2e1 * t195 * t272, 0.2e1 * t185 * t195 - 0.2e1 * t296 * t310, qJDD(3) * t195 + t198 * t200, qJDD(3) * t198 - t195 * t200, 0, t207 * t195 + t198 * t368, -t195 * t368 + t207 * t198, -t58 * t321 + (t198 * t298 - t276) * t140 (-t138 * t197 - t140 * t194) * t302 + (t344 - t197 * t59 + (t138 * t194 - t140 * t197) * qJD(4)) * t195 (t58 - t278) * t198 + (t233 + t305) * t195 (t59 + t279) * t198 + (-qJD(3) * t138 - t234) * t195, -t136 * t198 - t174 * t303, t259 * t136 + t376 * t174 + ((pkin(8) * t138 + t194 * t87) * qJD(3) + t263) * t198 + (-t253 + t87 * t299 - t29 * qJD(3) + t25 * t194 + (t59 - t279) * pkin(8)) * t195 + t230, -t311 * t136 + t379 * t174 + (t87 * t298 + (-t277 + t305) * pkin(8) - t262) * t198 + (-t252 - t87 * t301 - t30 * qJD(3) + t25 * t197 + (-t58 - t278) * pkin(8)) * t195 + t231, -t58 * t92 + t59 * t91 + t351 * t140 + t352 * t138 + t240 * t302 + (t194 * t4 + t197 * t5 + (-t19 * t194 + t197 * t20) * qJD(4) + t221) * t195, -t115 * t59 + t92 * t136 - t46 * t138 + (-t304 * t31 - t5) * t198 - t351 * t174 + (qJD(3) * t19 - t6 * t194 - t299 * t31 + t253) * t195 - t230, t115 * t58 - t91 * t136 - t46 * t140 + (-t298 * t31 + t4) * t198 + t352 * t174 + (-qJD(3) * t20 - t6 * t197 + t301 * t31 + t252) * t195 - t231, t6 * t115 + t4 * t91 + t5 * t92 - g(1) * t216 - g(2) * t217 - g(3) * t223 + (t46 - t254) * t31 + t352 * t20 + t351 * t19, -t60 * t58 - t65 * t59 + t354 * t140 - t353 * t138 + (-t194 * t9 + t197 * t7) * t302 + (t1 * t197 - t194 * t2 + (-t194 * t7 - t197 * t9) * qJD(4) + t221) * t195, t65 * t136 - t14 * t140 + t86 * t58 + (-t17 * t298 - t2) * t198 - t353 * t174 + (qJD(3) * t9 - t3 * t197 + t252 + t289) * t195 - t231, -t60 * t136 + t14 * t138 + t86 * t59 + (t17 * t304 + t1) * t198 + t354 * t174 + (-qJD(3) * t7 + t3 * t194 - t253 + t288) * t195 + t230, t3 * t86 + t1 * t60 + t2 * t65 - g(1) * (-pkin(5) * t333 + qJ(6) * t54 + t216) - g(2) * (-pkin(5) * t334 + qJ(6) * t52 + t217) - g(3) * (pkin(5) * t287 + qJ(6) * t101 + t223) + t353 * t9 + t354 * t7 + (t14 - t254) * t17; 0, 0, 0, 0, -t195 * t201 * t198, t310 * t201, t295, t185, qJDD(3), t97 * qJD(3) - t147 * t307 + t218 + t229, -t247 + (-qJD(2) * t147 - t103) * t198 + t228, -t197 * t328 - t344 (-t58 + t330) * t197 + (-t59 + t328) * t194 (-t140 * t195 + t174 * t319) * qJD(2) + t234 (t138 * t195 - t174 * t322) * qJD(2) + t233, t174 * t307, t29 * t307 - pkin(3) * t59 - t97 * t138 - t84 * t174 + t220 * t194 + (-t25 + (t142 + t349) * t174 + t229) * t197, pkin(3) * t58 - t340 * t174 + t30 * t307 - t97 * t140 + t220 * t197 + (t212 + t25) * t194, -t44 * t138 - t45 * t140 + (-t4 - t174 * t19 + (qJD(4) * t140 - t59) * pkin(9)) * t197 + (t5 - t346 + (qJD(4) * t138 - t58) * pkin(9)) * t194 - t228, -t138 * t343 + t45 * t174 - t19 * t307 + t194 * t219 - t197 * t210 + t377 * t59, -t140 * t343 - t44 * t174 + t194 * t210 + t197 * t219 + t20 * t307 - t377 * t58, -t6 * t377 - t20 * t44 - t19 * t45 - g(1) * t293 - g(2) * t294 - g(3) * t284 + t343 * t31 + (qJD(4) * t240 + t5 * t194 - t4 * t197 - t228) * pkin(9), -t156 * t58 - t157 * t59 + t341 * t140 - t342 * t138 + (t2 - t356) * t197 + (t9 * t174 + t1) * t194 - t228, -t288 + t130 * t58 + t157 * t136 - t342 * t174 - t350 * t140 + (t17 * t319 - t195 * t9) * qJD(2) + t224 * t194, t289 + t130 * t59 - t156 * t136 + t341 * t174 + t350 * t138 + (-t17 * t322 + t195 * t7) * qJD(2) + t224 * t197, t3 * t130 + t1 * t156 + t2 * t157 - g(1) * (t336 * t74 + t364 * t75 + t293) - g(2) * (t336 * t72 + t364 * t73 + t294) - g(3) * (-qJ(6) * t332 + t126 * t364 + t284) + t342 * t9 + t341 * t7 + t350 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t329, t135 - t366, t27, t202, t136, -t140 * t87 + t214 - t345, t138 * t87 + t174 * t29 + t213, pkin(4) * t58 - t348 + (-t20 - t30) * t140 + (t19 + t317) * t138, t138 * t71 + t206 + t345 - 0.2e1 * t357, -t31 * t138 + t71 * t140 + t174 * t317 - t162 - t213 + t365, -t5 * pkin(4) - g(1) * t267 - g(2) * t268 - g(3) * t266 - t4 * qJ(5) - t19 * t30 + t20 * t317 - t31 * t71, -t348 + t355 * t58 + (t9 + t315) * t140 + (t7 - t316) * t138, t140 * t48 - t174 * t237 - 0.2e1 * t162 + t204 + t365, -t138 * t48 + (-0.2e1 * qJD(6) - t16) * t174 + 0.2e1 * t274 - t205, -t1 * t355 + t2 * qJ(5) - t17 * t48 - g(1) * (-qJ(6) * t38 + t267) - g(2) * (-qJ(6) * t36 + t268) - g(3) * (-qJ(6) * t79 + t266) + t316 * t9 + t315 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t215, t373, t206 - t346 - t357, t27, t373, -t215 (qJD(6) + t9) * t174 - t274 + t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, t136 + t329, -t171 - t366, t204 - t356 - t372;];
tau_reg  = t8;
