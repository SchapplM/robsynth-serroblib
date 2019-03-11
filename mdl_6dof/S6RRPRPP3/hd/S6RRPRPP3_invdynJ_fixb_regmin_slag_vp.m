% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRPP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRPP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:57:13
% EndTime: 2019-03-09 09:57:26
% DurationCPUTime: 5.80s
% Computational Cost: add. (7020->565), mult. (15555->670), div. (0->0), fcn. (11050->10), ass. (0->264)
t219 = sin(qJ(4));
t216 = cos(pkin(9));
t222 = cos(qJ(2));
t312 = qJD(1) * qJD(2);
t295 = t222 * t312;
t220 = sin(qJ(2));
t310 = t220 * qJDD(1);
t254 = t295 + t310;
t215 = sin(pkin(9));
t311 = t215 * qJDD(2);
t234 = t216 * t254 + t311;
t319 = qJD(1) * t220;
t301 = t215 * t319;
t316 = qJD(2) * t216;
t147 = -t301 + t316;
t304 = t216 * t319;
t317 = qJD(2) * t215;
t148 = t304 + t317;
t368 = cos(qJ(4));
t258 = -t219 * t147 - t148 * t368;
t323 = t254 * t215;
t270 = qJDD(2) * t216 - t323;
t229 = qJD(4) * t258 - t219 * t234 + t368 * t270;
t318 = qJD(1) * t222;
t194 = -qJD(4) + t318;
t346 = t258 * t194;
t387 = t229 - t346;
t306 = t368 * t216;
t286 = t222 * t306;
t305 = t215 * t318;
t300 = qJD(4) * t368;
t313 = qJD(4) * t219;
t377 = -t215 * t313 + t216 * t300;
t326 = -qJD(1) * t286 + t219 * t305 + t377;
t89 = -t368 * t147 + t148 * t219;
t373 = t89 ^ 2;
t87 = t258 ^ 2;
t386 = t194 * t89;
t356 = pkin(8) + qJ(3);
t166 = t356 * t215;
t167 = t356 * t216;
t271 = pkin(2) * t220 - qJ(3) * t222;
t154 = t271 * qJD(1);
t110 = pkin(7) * t301 + t216 * t154;
t335 = t216 * t222;
t264 = pkin(3) * t220 - pkin(8) * t335;
t77 = qJD(1) * t264 + t110;
t138 = t215 * t154;
t336 = t216 * t220;
t337 = t215 * t222;
t256 = -pkin(7) * t336 - pkin(8) * t337;
t94 = qJD(1) * t256 + t138;
t385 = qJD(3) * t306 - t166 * t300 - t368 * t94 + (-qJD(3) * t215 - qJD(4) * t167 - t77) * t219;
t205 = t222 * qJDD(1);
t253 = t220 * t312 - t205;
t150 = qJDD(4) + t253;
t358 = pkin(4) + qJ(6);
t299 = t358 * t150;
t152 = t215 * t368 + t219 * t216;
t137 = t152 * qJD(4);
t246 = t222 * t152;
t325 = -qJD(1) * t246 + t137;
t223 = cos(qJ(1));
t333 = t220 * t223;
t221 = sin(qJ(1));
t334 = t220 * t221;
t384 = g(1) * t333 + g(2) * t334;
t272 = pkin(2) * t222 + qJ(3) * t220;
t162 = -pkin(1) - t272;
t140 = t162 * qJD(1);
t201 = pkin(7) * t318;
t168 = qJD(2) * qJ(3) + t201;
t97 = t216 * t140 - t168 * t215;
t55 = -pkin(3) * t318 - pkin(8) * t148 + t97;
t98 = t215 * t140 + t216 * t168;
t62 = pkin(8) * t147 + t98;
t24 = t219 * t62 - t368 * t55;
t265 = -pkin(5) * t258 + t24;
t328 = qJD(5) + t265;
t383 = -qJ(6) * t229 + t89 * qJD(6);
t351 = qJ(5) * t319 - t385;
t105 = -t219 * t166 + t167 * t368;
t382 = qJD(3) * t152 + qJD(4) * t105 - t219 * t94 + t368 * t77;
t188 = t194 ^ 2;
t381 = -t87 - t188;
t363 = g(2) * t221;
t278 = g(1) * t223 + t363;
t380 = t222 * t278;
t261 = t278 * t220;
t141 = t150 * qJ(5);
t176 = qJD(5) * t194;
t379 = t176 - t141;
t142 = pkin(3) * t305 + t201;
t378 = -t326 * qJ(5) - t152 * qJD(5) - t142;
t376 = pkin(7) * t295 + qJDD(3);
t375 = pkin(5) * t229 + qJDD(6);
t374 = -0.2e1 * pkin(1);
t372 = 0.2e1 * t141;
t370 = pkin(5) * t89;
t35 = -t147 * t300 + t148 * t313 - t219 * t270 - t368 * t234;
t369 = t35 * pkin(5);
t367 = pkin(3) * t215;
t366 = pkin(7) * t147;
t365 = g(1) * t221;
t362 = g(2) * t223;
t361 = g(3) * t220;
t211 = g(3) * t222;
t360 = t150 * pkin(4);
t359 = t258 * t89;
t357 = pkin(5) + t356;
t151 = t215 * t219 - t306;
t355 = t151 * qJD(6) + t325 * t358 + t378;
t354 = -pkin(5) * t325 - t351;
t298 = t358 * t220;
t353 = pkin(5) * t326 + qJD(1) * t298 + t382;
t352 = pkin(4) * t325 + t378;
t25 = t219 * t55 + t368 * t62;
t350 = pkin(4) * t319 + t382;
t348 = qJ(5) * t229;
t347 = t194 * t25;
t345 = t89 * qJ(5);
t114 = pkin(7) * t335 + t215 * t162;
t338 = t215 * t220;
t103 = -pkin(8) * t338 + t114;
t146 = t216 * t162;
t96 = -pkin(8) * t336 + t146 + (-pkin(7) * t215 - pkin(3)) * t222;
t344 = t368 * t103 + t219 * t96;
t123 = -pkin(7) * t253 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t134 = qJD(2) * t271 - t220 * qJD(3);
t85 = qJD(1) * t134 + qJDD(1) * t162;
t47 = t216 * t123 + t215 * t85;
t343 = qJDD(2) * pkin(2);
t212 = pkin(9) + qJ(4);
t203 = sin(t212);
t342 = t203 * t220;
t341 = t203 * t222;
t204 = cos(t212);
t340 = t204 * t220;
t339 = t204 * t222;
t332 = t221 * t222;
t198 = t216 * pkin(3) + pkin(2);
t171 = t222 * t198;
t331 = t222 * t223;
t330 = t223 * t203;
t329 = -qJD(5) - t24;
t15 = t25 - t370;
t327 = -qJD(6) - t15;
t315 = qJD(2) * t220;
t307 = pkin(7) * t315;
t101 = t216 * t134 + t215 * t307;
t199 = pkin(7) * t310;
t322 = -t199 - t211;
t314 = qJD(2) * t222;
t303 = t215 * t314;
t143 = pkin(3) * t303 + pkin(7) * t314;
t155 = pkin(3) * t338 + t220 * pkin(7);
t321 = t223 * pkin(1) + t221 * pkin(7);
t213 = t220 ^ 2;
t320 = -t222 ^ 2 + t213;
t297 = t357 * t223;
t294 = -pkin(1) - t171;
t132 = t199 - t343 + t376;
t293 = -t132 - t211;
t46 = -t215 * t123 + t216 * t85;
t28 = pkin(3) * t253 - pkin(8) * t234 + t46;
t37 = pkin(8) * t270 + t47;
t292 = -t219 * t28 - t55 * t300 + t62 * t313 - t368 * t37;
t291 = t219 * t37 - t368 * t28 + t62 * t300 + t55 * t313;
t290 = -qJD(2) * pkin(2) + qJD(3);
t124 = t203 * t332 + t204 * t223;
t125 = t204 * t332 - t330;
t289 = -t124 * pkin(4) + qJ(5) * t125;
t126 = -t221 * t204 + t222 * t330;
t127 = t203 * t221 + t204 * t331;
t288 = -t126 * pkin(4) + qJ(5) * t127;
t287 = -qJ(5) * t203 - t198;
t285 = g(3) * (pkin(4) * t339 + qJ(5) * t341 + t171);
t196 = g(1) * t334;
t283 = -g(2) * t333 + t196;
t20 = qJ(5) * t194 - t25;
t281 = -t219 * t103 + t368 * t96;
t42 = qJ(5) * t222 - t344;
t280 = g(1) * t124 - g(2) * t126;
t279 = g(1) * t125 - g(2) * t127;
t277 = -t362 + t365;
t276 = -t87 - t373;
t161 = pkin(7) * t319 + t290;
t131 = -t219 * t338 + t220 * t306;
t274 = -t131 * qJ(5) + t155;
t273 = -qJDD(5) - t291;
t268 = -g(3) * t341 + t203 * t384;
t267 = -g(3) * t339 + t204 * t384;
t43 = t222 * pkin(4) - t281;
t4 = t292 + t379;
t263 = -t152 * qJ(5) - t198;
t66 = qJD(2) * t264 + t101;
t121 = t215 * t134;
t78 = qJD(2) * t256 + t121;
t260 = t103 * t300 + t219 * t78 + t96 * t313 - t368 * t66;
t259 = -pkin(7) * qJDD(2) + t312 * t374;
t104 = t166 * t368 + t219 * t167;
t257 = -t103 * t313 + t219 * t66 + t96 * t300 + t368 * t78;
t209 = t223 * pkin(7);
t255 = -t125 * pkin(4) - t124 * qJ(5) + t223 * t367 - t334 * t356 + t209;
t252 = t216 * t310 + t311;
t225 = qJD(1) ^ 2;
t251 = pkin(1) * t225 + t278;
t109 = -t147 * pkin(3) + t161;
t224 = qJD(2) ^ 2;
t250 = pkin(7) * t224 + qJDD(1) * t374 + t362;
t71 = -qJD(2) * t286 + t137 * t220 + t219 * t303;
t249 = t71 * qJ(5) - t131 * qJD(5) + t143;
t247 = t127 * pkin(4) + t126 * qJ(5) + t198 * t331 + t221 * t367 + t321;
t245 = -t104 * t150 + t267;
t244 = t105 * t150 + t268;
t243 = t150 + t359;
t242 = -t361 - t380;
t241 = g(1) * t126 + g(2) * t124 + g(3) * t342 - t291;
t17 = -t35 - t386;
t240 = t35 - t386;
t239 = -t261 - t343;
t238 = qJ(5) * t258 + t109;
t237 = t293 + t261;
t236 = -qJDD(5) + t241;
t235 = g(1) * t127 + g(2) * t125 + g(3) * t340 + t292;
t67 = -pkin(3) * t270 + t132;
t30 = pkin(4) * t89 + t238;
t233 = -t258 * t30 - t236;
t10 = -qJ(5) * t315 + qJD(5) * t222 - t257;
t16 = t358 * t89 + t238;
t232 = -t16 * t258 - t236 - t369;
t230 = -t16 * t89 - t235 + t375;
t6 = -pkin(4) * t229 + t35 * qJ(5) + qJD(5) * t258 + t67;
t228 = t229 + t346;
t226 = t211 - t261 + t6;
t169 = qJ(5) * t340;
t130 = t152 * t220;
t113 = -pkin(7) * t337 + t146;
t111 = -pkin(7) * t304 + t138;
t102 = -t216 * t307 + t121;
t86 = pkin(4) * t151 + t263;
t72 = qJD(2) * t246 + t377 * t220;
t70 = -t151 * pkin(5) + t105;
t69 = t152 * pkin(5) + t104;
t60 = t151 * t358 + t263;
t56 = pkin(4) * t130 + t274;
t45 = -pkin(4) * t258 + t345;
t44 = t130 * t358 + t274;
t31 = -pkin(5) * t130 - t42;
t29 = t131 * pkin(5) + t222 * qJ(6) + t43;
t23 = -t258 * t358 + t345;
t19 = pkin(4) * t194 - t329;
t18 = pkin(4) * t72 + t249;
t13 = qJD(6) - t20 - t370;
t12 = t194 * t358 + t328;
t11 = -pkin(4) * t315 + t260;
t9 = t130 * qJD(6) + t358 * t72 + t249;
t8 = -pkin(5) * t72 - t10;
t7 = -t71 * pkin(5) - qJD(2) * t298 + t222 * qJD(6) + t260;
t5 = -t273 - t360;
t3 = t6 + t383;
t2 = -t4 + t375;
t1 = qJD(6) * t194 - t273 - t299 - t369;
t14 = [qJDD(1), t277, t278, qJDD(1) * t213 + 0.2e1 * t220 * t295, 0.2e1 * t205 * t220 - 0.2e1 * t312 * t320, qJDD(2) * t220 + t222 * t224, qJDD(2) * t222 - t220 * t224, 0, t259 * t220 + (-t250 + t365) * t222, t220 * t250 + t222 * t259 - t196, -t278 * t215 + (-pkin(7) * t270 + t132 * t215 + (qJD(1) * t113 + t97) * qJD(2)) * t220 + (-t101 * qJD(1) - t113 * qJDD(1) - t46 + t277 * t216 + (t161 * t215 - t366) * qJD(2)) * t222, -t278 * t216 + (t132 * t216 + (-qJD(1) * t114 - t98) * qJD(2) + t252 * pkin(7)) * t220 + (t102 * qJD(1) + t114 * qJDD(1) + t47 - t277 * t215 + (t161 * t216 + (t148 + t304) * pkin(7)) * qJD(2)) * t222, t102 * t147 - t114 * t323 - t101 * t148 + (-qJDD(2) * t113 - t220 * t47 - t314 * t98) * t215 + (t114 * qJDD(2) - t113 * t254 - t46 * t220 - t314 * t97) * t216 + t283, t47 * t114 + t98 * t102 + t46 * t113 + t97 * t101 - g(1) * t209 - g(2) * (t223 * t272 + t321) - t162 * t365 + (t132 * t220 + t161 * t314) * pkin(7), -t131 * t35 + t258 * t71, t130 * t35 + t131 * t229 + t258 * t72 + t71 * t89, t131 * t150 + t194 * t71 + t222 * t35 - t258 * t315, -t130 * t150 + t194 * t72 - t222 * t229 - t315 * t89, -t150 * t222 - t194 * t315, t109 * t72 + t67 * t130 + t143 * t89 + t150 * t281 - t155 * t229 + t194 * t260 + t222 * t291 - t24 * t315 + t279, -t109 * t71 + t67 * t131 - t143 * t258 - t150 * t344 - t155 * t35 + t194 * t257 - t222 * t292 - t25 * t315 - t280, t10 * t89 - t11 * t258 + t130 * t4 + t131 * t5 - t19 * t71 + t20 * t72 - t229 * t42 - t35 * t43 + t283, -t11 * t194 - t130 * t6 + t150 * t43 - t18 * t89 + t19 * t315 - t222 * t5 + t229 * t56 - t30 * t72 - t279, t10 * t194 - t131 * t6 - t150 * t42 + t18 * t258 - t20 * t315 + t222 * t4 + t30 * t71 + t35 * t56 + t280, t6 * t56 + t30 * t18 + t4 * t42 + t20 * t10 + t5 * t43 + t19 * t11 - g(1) * (t221 * t294 + t255) - g(2) * (t333 * t356 + t247) t1 * t131 - t12 * t71 - t13 * t72 - t130 * t2 + t229 * t31 - t258 * t7 - t29 * t35 - t8 * t89 + t283, t13 * t315 - t131 * t3 + t150 * t31 + t16 * t71 - t194 * t8 - t2 * t222 + t258 * t9 + t35 * t44 + t280, t1 * t222 - t12 * t315 + t130 * t3 - t150 * t29 + t16 * t72 + t194 * t7 - t229 * t44 + t89 * t9 + t279, t3 * t44 + t16 * t9 + t1 * t29 + t12 * t7 + t2 * t31 + t13 * t8 - g(1) * (-t125 * qJ(6) + t255) - g(2) * (t127 * qJ(6) + t220 * t297 + t247) - (-pkin(5) * t220 + t294) * t365; 0, 0, 0, -t220 * t225 * t222, t320 * t225, t310, t205, qJDD(2), t220 * t251 + t322, t361 + (-pkin(7) * qJDD(1) + t251) * t222, t215 * qJ(3) * t205 - pkin(2) * t323 + (t237 + t343) * t216 + ((-qJ(3) * t317 - t97) * t220 + (t366 + t110 + (qJD(3) - t161) * t215) * t222) * qJD(1), -t271 * t216 * qJDD(1) + (t239 - t293) * t215 + ((-qJ(3) * t316 + t98) * t220 + (-pkin(7) * t148 - t111 + (-t161 + t290) * t216) * t222) * qJD(1), t110 * t148 - t111 * t147 + (qJ(3) * t270 + qJD(3) * t147 + t318 * t97 + t47) * t216 + (qJ(3) * t234 + qJD(3) * t148 + t318 * t98 - t46) * t215 + t242, -t161 * t201 - t97 * t110 - t98 * t111 + (-t215 * t97 + t216 * t98) * qJD(3) + t237 * pkin(2) + (-t46 * t215 + t47 * t216 + t242) * qJ(3), -t35 * t152 - t258 * t326, t35 * t151 + t152 * t229 + t258 * t325 - t326 * t89, t152 * t150 - t194 * t326 + t258 * t319, -t151 * t150 + t194 * t325 + t319 * t89, t194 * t319, t325 * t109 - t142 * t89 + t67 * t151 + t382 * t194 + t198 * t229 + t24 * t319 + t245, t326 * t109 + t142 * t258 + t67 * t152 + t194 * t385 + t198 * t35 + t25 * t319 - t244, -t104 * t35 + t105 * t229 + t4 * t151 + t5 * t152 + t19 * t326 + t20 * t325 - t258 * t350 + t351 * t89 + t242, -t6 * t151 - t19 * t319 - t194 * t350 + t229 * t86 - t30 * t325 - t352 * t89 - t245, -t6 * t152 + t194 * t351 + t20 * t319 + t258 * t352 - t30 * t326 + t86 * t35 + t244, t6 * t86 - t4 * t105 + t5 * t104 - t285 + t352 * t30 - t356 * t380 + t351 * t20 + t350 * t19 + (-g(3) * t356 + t278 * (pkin(4) * t204 - t287)) * t220, t1 * t152 + t12 * t326 - t13 * t325 - t2 * t151 + t229 * t70 - t258 * t353 - t69 * t35 - t354 * t89 + t242, -t13 * t319 + t70 * t150 - t3 * t152 - t16 * t326 - t194 * t354 + t258 * t355 + t60 * t35 + t268, t12 * t319 - t69 * t150 + t3 * t151 + t16 * t325 + t194 * t353 - t229 * t60 + t355 * t89 + t267, t3 * t60 + t1 * t69 + t2 * t70 - t285 + t355 * t16 + t354 * t13 + t353 * t12 + (-g(3) * qJ(6) * t204 - g(1) * t297 - t357 * t363) * t222 + (-g(3) * t357 + t278 * (t204 * t358 - t287)) * t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148 * t318 - t270 (-t147 + t316) * t318 + t252, -t147 ^ 2 - t148 ^ 2, -t98 * t147 + t97 * t148 + t239 - t322 + t376, 0, 0, 0, 0, 0, -t387, -t240, t276, t387, t240, t19 * t258 - t20 * t89 + t226, t276, t240, -t387, t12 * t258 + t13 * t89 + t226 + t383; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t359, t87 - t373, t17, t228, t150, t109 * t258 + t241 - t347, t109 * t89 + t194 * t24 + t235, pkin(4) * t35 + t348 - (-t20 - t25) * t258 + (t19 + t329) * t89, t45 * t89 + t233 + t347 - 0.2e1 * t360, t194 * t329 - t258 * t45 - t30 * t89 - t176 - t235 + t372, -t4 * qJ(5) - t5 * pkin(4) - t30 * t45 - t19 * t25 - g(1) * t288 - g(2) * t289 - g(3) * (-pkin(4) * t342 + t169) + t329 * t20, t348 + t358 * t35 - (t13 + t327) * t258 + (t12 - t328) * t89, -t194 * t265 - t23 * t258 - 0.2e1 * t176 + t230 + t372, -t23 * t89 + (-0.2e1 * qJD(6) - t15) * t194 + 0.2e1 * t299 - t232, -t1 * t358 + t2 * qJ(5) - t16 * t23 - g(1) * (-qJ(6) * t126 + t288) - g(2) * (-qJ(6) * t124 + t289) - g(3) * (-t203 * t298 + t169) + t328 * t13 + t327 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t243, t381, -t194 * t20 + t233 - t360, t17, t381, -t243 (qJD(6) + t13) * t194 - t299 + t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, t150 - t359, -t188 - t373, -t12 * t194 + t230 - t379;];
tau_reg  = t14;
