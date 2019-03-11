% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRP7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:21:04
% EndTime: 2019-03-09 06:21:17
% DurationCPUTime: 5.46s
% Computational Cost: add. (10729->519), mult. (25011->635), div. (0->0), fcn. (19416->14), ass. (0->259)
t216 = sin(qJ(4));
t221 = -pkin(9) - pkin(8);
t293 = qJD(4) * t221;
t212 = sin(pkin(10));
t217 = sin(qJ(3));
t213 = cos(pkin(10));
t362 = cos(qJ(3));
t289 = t362 * t213;
t248 = -t217 * t212 + t289;
t395 = t248 * qJD(1);
t331 = t395 * t216;
t170 = t362 * t212 + t217 * t213;
t374 = t170 * qJD(1);
t111 = pkin(3) * t374 - pkin(8) * t395;
t219 = cos(qJ(4));
t351 = pkin(7) + qJ(2);
t182 = t351 * t212;
t171 = qJD(1) * t182;
t183 = t351 * t213;
t172 = qJD(1) * t183;
t380 = -t362 * t171 - t217 * t172;
t341 = t216 * t111 + t219 * t380;
t400 = -pkin(9) * t331 - t216 * t293 + t341;
t330 = t395 * t219;
t97 = t219 * t111;
t399 = -pkin(4) * t374 + pkin(9) * t330 + t216 * t380 + t219 * t293 - t97;
t163 = t170 * qJD(3);
t300 = t212 * qJDD(1);
t261 = -qJDD(1) * t289 + t217 * t300;
t114 = qJD(1) * t163 + t261;
t108 = qJDD(4) + t114;
t105 = qJDD(5) + t108;
t151 = qJD(4) - t395;
t147 = qJD(5) + t151;
t215 = sin(qJ(5));
t361 = cos(qJ(5));
t288 = t361 * t216;
t174 = t215 * t219 + t288;
t287 = t361 * t219;
t321 = t215 * t216;
t245 = t287 - t321;
t371 = qJD(4) + qJD(5);
t280 = t361 * qJD(5);
t378 = t361 * qJD(4) + t280;
t340 = -t219 * t378 + t245 * t395 + t321 * t371;
t269 = t174 * t105 - t340 * t147;
t128 = t219 * qJD(3) - t216 * t374;
t129 = qJD(3) * t216 + t219 * t374;
t247 = t215 * t128 + t361 * t129;
t344 = t374 * t247;
t398 = t269 + t344;
t120 = t371 * t174;
t339 = -t174 * t395 + t120;
t270 = t245 * t105 - t339 * t147;
t73 = -t361 * t128 + t129 * t215;
t345 = t374 * t73;
t397 = t270 - t345;
t305 = qJD(4) * t216;
t396 = t305 - t331;
t109 = -qJD(3) * pkin(3) - t380;
t69 = -t128 * pkin(4) + t109;
t31 = t73 * pkin(5) - qJ(6) * t247 + t69;
t394 = t31 * t73;
t393 = t69 * t73;
t352 = t247 * t73;
t337 = qJDD(1) * pkin(1);
t202 = qJDD(2) - t337;
t218 = sin(qJ(1));
t220 = cos(qJ(1));
t379 = g(1) * t218 - g(2) * t220;
t260 = -t202 + t379;
t162 = t248 * qJD(3);
t237 = t170 * qJDD(1);
t226 = qJD(1) * t162 + t237;
t392 = qJD(3) * qJD(4) + t226;
t304 = qJD(4) * t219;
t319 = t216 * t162;
t391 = t170 * t304 + t319;
t210 = pkin(10) + qJ(3);
t203 = sin(t210);
t324 = t203 * t220;
t325 = t203 * t218;
t390 = -g(1) * t324 - g(2) * t325;
t184 = t221 * t216;
t185 = t221 * t219;
t127 = t215 * t184 - t361 * t185;
t211 = qJ(4) + qJ(5);
t205 = sin(t211);
t264 = g(1) * t220 + g(2) * t218;
t204 = cos(t210);
t354 = g(3) * t204;
t233 = t264 * t203 - t354;
t389 = t127 * t105 + t233 * t205;
t388 = qJD(3) * t374;
t363 = t247 ^ 2;
t387 = -t73 ^ 2 + t363;
t225 = (qJD(3) * (qJD(4) + t395) + t237) * t219 + (-qJD(4) * t374 + qJDD(3)) * t216;
t294 = t392 * t216 + t374 * t304;
t257 = t219 * qJDD(3) - t294;
t303 = qJD(5) * t215;
t27 = -t128 * t280 + t129 * t303 - t215 * t257 - t361 * t225;
t13 = t147 * t73 - t27;
t43 = pkin(5) * t247 + qJ(6) * t73;
t138 = t147 * qJD(6);
t93 = t105 * qJ(6);
t385 = t138 + t93;
t94 = t219 * t108;
t384 = t151 * t305 - t94;
t246 = t361 * t184 + t215 * t185;
t383 = -t246 * qJD(5) - t399 * t215 + t400 * t361;
t382 = -t127 * qJD(5) + t400 * t215 + t399 * t361;
t117 = -t217 * t171 + t362 * t172;
t268 = t396 * pkin(4) - t117;
t320 = t216 * t108;
t285 = t170 * t305;
t312 = t219 * t162;
t377 = -t285 + t312;
t376 = -t362 * t182 - t217 * t183;
t311 = t219 * t220;
t318 = t216 * t218;
t152 = t204 * t318 + t311;
t313 = t218 * t219;
t317 = t216 * t220;
t154 = -t204 * t317 + t313;
t375 = -g(1) * t154 + g(2) * t152;
t99 = t105 * pkin(5);
t373 = t99 - qJDD(6);
t372 = qJ(2) * qJDD(1);
t28 = t128 * t303 + t129 * t280 + t215 * t225 - t361 * t257;
t370 = t147 * t247 - t28;
t206 = cos(t211);
t322 = t206 * t220;
t323 = t205 * t218;
t134 = t204 * t323 + t322;
t310 = t220 * t205;
t314 = t218 * t206;
t136 = t204 * t310 - t314;
t110 = qJD(3) * pkin(8) + t117;
t302 = qJD(1) * qJD(2);
t364 = t351 * qJDD(1) + t302;
t144 = t364 * t212;
t145 = t364 * t213;
t250 = -t217 * t144 + t362 * t145;
t61 = qJDD(3) * pkin(8) + qJD(3) * t380 + t250;
t299 = t213 * qJDD(1);
t65 = -pkin(2) * t299 + t114 * pkin(3) - pkin(8) * t226 + t202;
t197 = t213 * pkin(2) + pkin(1);
t179 = -t197 * qJD(1) + qJD(2);
t86 = -pkin(3) * t395 - pkin(8) * t374 + t179;
t243 = -t110 * t305 + t216 * t65 + t219 * t61 + t86 * t304;
t12 = t257 * pkin(9) + t243;
t54 = -t110 * t216 + t219 * t86;
t45 = -pkin(9) * t129 + t54;
t37 = pkin(4) * t151 + t45;
t55 = t219 * t110 + t216 * t86;
t46 = pkin(9) * t128 + t55;
t262 = -t110 * t304 + t219 * t65;
t9 = t108 * pkin(4) - pkin(9) * t225 - t216 * t61 - t86 * t305 + t262;
t278 = t215 * t12 + t46 * t280 + t37 * t303 - t361 * t9;
t327 = t203 * t205;
t235 = g(1) * t136 + g(2) * t134 + g(3) * t327 - t278;
t228 = t247 * t31 - t235 - t373;
t369 = -t69 * t247 + t235;
t368 = -t151 ^ 2 * t219 - t320;
t355 = g(3) * t203;
t367 = t264 * t204 + t355;
t366 = -t245 * t27 - t247 * t339;
t113 = -pkin(3) * t248 - pkin(8) * t170 - t197;
t125 = -t217 * t182 + t362 * t183;
t118 = t219 * t125;
t87 = t248 * qJD(2) + qJD(3) * t376;
t112 = pkin(3) * t163 - pkin(8) * t162;
t98 = t219 * t112;
t24 = -pkin(9) * t312 + pkin(4) * t163 - t216 * t87 + t98 + (-t118 + (pkin(9) * t170 - t113) * t216) * qJD(4);
t103 = t219 * t113;
t328 = t170 * t219;
t50 = -pkin(4) * t248 - pkin(9) * t328 - t125 * t216 + t103;
t308 = t216 * t113 + t118;
t329 = t170 * t216;
t57 = -pkin(9) * t329 + t308;
t253 = t215 * t50 + t361 * t57;
t242 = t216 * t112 + t113 * t304 - t125 * t305 + t219 * t87;
t30 = -pkin(9) * t391 + t242;
t365 = -qJD(5) * t253 - t215 * t30 + t361 * t24;
t353 = g(3) * t216;
t349 = t339 * pkin(5) + t340 * qJ(6) - qJD(6) * t174 + t268;
t348 = -qJ(6) * t374 - t383;
t347 = pkin(5) * t374 - t382;
t295 = t361 * t46;
t17 = t215 * t37 + t295;
t346 = t147 * t17;
t342 = t215 * t46;
t19 = t361 * t45 - t342;
t338 = pkin(4) * t280 + qJD(6) - t19;
t336 = qJDD(3) * pkin(3);
t334 = t128 * t374;
t333 = t129 * t151;
t332 = t129 * t374;
t326 = t203 * t206;
t16 = t361 * t37 - t342;
t309 = qJD(6) - t16;
t307 = t212 ^ 2 + t213 ^ 2;
t306 = qJD(3) * t217;
t201 = t219 * pkin(4) + pkin(3);
t283 = t109 * t304;
t282 = qJD(3) * t362;
t279 = t361 * t12 + t215 * t9 + t37 * t280 - t46 * t303;
t277 = pkin(4) * t216 + t351;
t275 = -qJD(4) * t86 - t61;
t274 = t307 * qJD(1) ^ 2;
t88 = t170 * qJD(2) - t182 * t306 + t183 * t282;
t272 = -t174 * t28 + t340 * t73;
t271 = 0.2e1 * t307;
t18 = t215 * t45 + t295;
t267 = pkin(4) * t303 - t18;
t266 = -g(1) * t134 + g(2) * t136;
t135 = t204 * t314 - t310;
t137 = t204 * t322 + t323;
t265 = g(1) * t135 - g(2) * t137;
t259 = t151 * t331 - t384;
t89 = pkin(4) * t329 - t376;
t258 = pkin(5) * t206 + qJ(6) * t205 + t201;
t255 = -t215 * t57 + t361 * t50;
t63 = t391 * pkin(4) + t88;
t252 = t201 * t204 - t203 * t221 + t197;
t251 = t215 * t24 + t50 * t280 + t361 * t30 - t57 * t303;
t241 = -t362 * t144 - t217 * t145 + t171 * t306 - t172 * t282;
t240 = t246 * t105 + (-t354 - t390) * t206;
t236 = t260 + t337;
t234 = g(1) * t137 + g(2) * t135 + g(3) * t326 - t279;
t62 = -t241 - t336;
t231 = t147 * t16 + t234;
t229 = t271 * t302 - t264;
t227 = -g(1) * (-t136 * pkin(5) + qJ(6) * t137) - g(2) * (-t134 * pkin(5) + qJ(6) * t135) - g(3) * (-pkin(5) * t327 + qJ(6) * t326);
t32 = -t257 * pkin(4) + t62;
t224 = t225 * t219;
t200 = -t361 * pkin(4) - pkin(5);
t196 = pkin(4) * t215 + qJ(6);
t177 = -t197 * qJDD(1) + qJDD(2);
t155 = t204 * t311 + t318;
t153 = -t204 * t313 + t317;
t115 = -pkin(5) * t245 - qJ(6) * t174 - t201;
t101 = t245 * t170;
t100 = t174 * t170;
t40 = t100 * pkin(5) - t101 * qJ(6) + t89;
t39 = t162 * t288 - t215 * t285 - t303 * t329 + (t162 * t215 + t378 * t170) * t219;
t38 = t120 * t170 - t162 * t287 + t215 * t319;
t34 = pkin(4) * t129 + t43;
t26 = pkin(5) * t248 - t255;
t25 = -qJ(6) * t248 + t253;
t15 = t147 * qJ(6) + t17;
t14 = -t147 * pkin(5) + t309;
t6 = pkin(5) * t39 + qJ(6) * t38 - qJD(6) * t101 + t63;
t5 = t28 * pkin(5) + t27 * qJ(6) - qJD(6) * t247 + t32;
t4 = -t163 * pkin(5) - t365;
t3 = qJ(6) * t163 - qJD(6) * t248 + t251;
t2 = t278 - t373;
t1 = t279 + t385;
t7 = [qJDD(1), t379, t264, t236 * t213, -t236 * t212, t271 * t372 + t229, t260 * pkin(1) + (t307 * t372 + t229) * qJ(2), t162 * t374 + t170 * t226, -t170 * t114 + t162 * t395 - t163 * t374 + t226 * t248, qJD(3) * t162 + t170 * qJDD(3), -qJD(3) * t163 + qJDD(3) * t248, 0, -qJD(3) * t88 + qJDD(3) * t376 - t114 * t197 + t163 * t179 - t177 * t248 + t204 * t379, -g(1) * t325 + g(2) * t324 - t87 * qJD(3) - t125 * qJDD(3) + t179 * t162 + t177 * t170 - t197 * t226, t129 * t377 + t170 * t224, t128 * t377 - t129 * t391 - t225 * t329 + t257 * t328, t129 * t163 + t151 * t377 + t170 * t94 - t225 * t248, t128 * t163 - t151 * t391 - t170 * t320 - t248 * t257, -t108 * t248 + t151 * t163 (-t125 * t304 + t98) * t151 + t103 * t108 - t262 * t248 + t54 * t163 - t88 * t128 + t376 * t257 + t170 * t283 - g(1) * t153 - g(2) * t155 + ((-qJD(4) * t113 - t87) * t151 - t125 * t108 - t275 * t248 + t62 * t170 + t109 * t162) * t216, -g(1) * t152 - g(2) * t154 - t308 * t108 + t109 * t377 + t88 * t129 - t242 * t151 - t55 * t163 - t225 * t376 + t243 * t248 + t62 * t328, -t101 * t27 - t247 * t38, t100 * t27 - t101 * t28 - t247 * t39 + t38 * t73, t101 * t105 - t147 * t38 + t163 * t247 + t248 * t27, -t100 * t105 - t147 * t39 - t163 * t73 + t248 * t28, -t105 * t248 + t147 * t163, t32 * t100 + t255 * t105 + t147 * t365 + t16 * t163 + t248 * t278 + t89 * t28 + t69 * t39 + t63 * t73 + t265, t32 * t101 - t253 * t105 - t251 * t147 - t17 * t163 + t247 * t63 + t248 * t279 - t89 * t27 - t69 * t38 + t266, t100 * t5 - t105 * t26 - t14 * t163 - t147 * t4 + t2 * t248 + t28 * t40 + t31 * t39 + t6 * t73 + t265, -t1 * t100 + t101 * t2 - t14 * t38 - t15 * t39 + t203 * t379 + t247 * t4 - t25 * t28 - t26 * t27 - t3 * t73, -t1 * t248 - t101 * t5 + t105 * t25 + t147 * t3 + t15 * t163 - t247 * t6 + t27 * t40 + t31 * t38 - t266, t1 * t25 + t15 * t3 + t5 * t40 + t31 * t6 + t2 * t26 + t14 * t4 - g(1) * (-pkin(5) * t135 - qJ(6) * t134) - g(2) * (pkin(5) * t137 + qJ(6) * t136) + (-g(1) * t277 - g(2) * t252) * t220 + (g(1) * t252 - g(2) * t277) * t218; 0, 0, 0, -t299, t300, -t274, -qJ(2) * t274 - t260, 0, 0, 0, 0, 0, t261 + 0.2e1 * t388, 0.2e1 * t395 * qJD(3) + t237, 0, 0, 0, 0, 0, t259 + t334, -t332 + t368, 0, 0, 0, 0, 0, t397, -t398, t397, t272 - t366, t398, t1 * t174 + t339 * t14 - t340 * t15 - t2 * t245 - t31 * t374 - t379; 0, 0, 0, 0, 0, 0, 0, -t374 * t395, t374 ^ 2 - t395 ^ 2, t237, -t261, qJDD(3), t117 * qJD(3) - t179 * t374 + t233 + t241, -t179 * t395 - t250 + t367, t216 * t225 + t333 * t219, t216 * t257 + t224 - t396 * t129 + (t304 - t330) * t128, -t332 - t368, t259 - t334, -t151 * t374, -pkin(3) * t294 - t54 * t374 + t117 * t128 - pkin(8) * t320 + (t233 + t336 - t62) * t219 + (-t97 + (t109 + t380) * t216 - pkin(8) * t304) * t151, -pkin(3) * t225 - t109 * t330 - t117 * t129 + t341 * t151 + t55 * t374 + t204 * t353 + t283 + t384 * pkin(8) + (t62 + t390) * t216, -t174 * t27 - t247 * t340, t272 + t366, t269 - t344, t270 + t345, -t147 * t374, t382 * t147 - t16 * t374 - t201 * t28 - t245 * t32 + t268 * t73 + t339 * t69 + t240, t383 * t147 + t17 * t374 + t32 * t174 + t201 * t27 + t268 * t247 - t340 * t69 - t389, t115 * t28 + t14 * t374 - t347 * t147 - t245 * t5 + t339 * t31 + t349 * t73 + t240, t1 * t245 - t127 * t28 - t340 * t14 - t339 * t15 + t174 * t2 + t246 * t27 + t247 * t347 - t348 * t73 - t367, t115 * t27 + t348 * t147 - t15 * t374 - t174 * t5 - t247 * t349 + t340 * t31 + t389, t1 * t127 + t5 * t115 - t2 * t246 + t349 * t31 + t348 * t15 + t347 * t14 + (-g(3) * t258 + t221 * t264) * t204 + (g(3) * t221 + t258 * t264) * t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129 * t128, -t128 ^ 2 + t129 ^ 2, t216 * qJDD(3) - t128 * t151 + t392 * t219 - t374 * t305, t257 + t333, t108, -t109 * t129 + t151 * t55 + (t275 + t355) * t216 + t262 + t375, g(1) * t155 - g(2) * t153 - t109 * t128 + t151 * t54 + t219 * t355 - t243, t352, t387, t13, t370, t105, t18 * t147 + (t361 * t105 - t129 * t73 - t147 * t303) * pkin(4) + t369, t19 * t147 + t393 + (-t215 * t105 - t129 * t247 - t147 * t280) * pkin(4) + t234, -t105 * t200 - t147 * t267 - t34 * t73 - t228, -t196 * t28 - t200 * t27 + (t15 + t267) * t247 + (t14 - t338) * t73, t105 * t196 + t338 * t147 + t247 * t34 - t234 + t385 - t394, t1 * t196 + t2 * t200 - t31 * t34 - t14 * t18 + t338 * t15 + (t14 * t303 + t203 * t353 + t375) * pkin(4) + t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t352, t387, t13, t370, t105, t346 + t369, t231 + t393, -t43 * t73 - t228 + t346 + t99, pkin(5) * t27 - qJ(6) * t28 + (t15 - t17) * t247 + (t14 - t309) * t73, t247 * t43 + 0.2e1 * t138 - t231 - t394 + 0.2e1 * t93, -t2 * pkin(5) + t1 * qJ(6) - t14 * t17 + t15 * t309 - t31 * t43 + t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - qJDD(5) - t261 + t352 - t388, t13, -t147 ^ 2 - t363, -t147 * t15 + t228;];
tau_reg  = t7;
