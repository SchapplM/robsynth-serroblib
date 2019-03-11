% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:26:54
% EndTime: 2019-03-09 02:27:06
% DurationCPUTime: 6.08s
% Computational Cost: add. (8706->597), mult. (15842->762), div. (0->0), fcn. (9943->12), ass. (0->287)
t198 = sin(qJ(5));
t201 = cos(qJ(5));
t322 = qJD(4) * t201;
t199 = sin(qJ(4));
t326 = qJD(1) * t199;
t130 = t198 * t326 + t322;
t313 = t198 * qJD(4);
t131 = t201 * t326 - t313;
t197 = sin(qJ(6));
t200 = cos(qJ(6));
t248 = t130 * t197 - t200 * t131;
t79 = t200 * t130 + t131 * t197;
t371 = t248 * t79;
t202 = cos(qJ(4));
t325 = qJD(1) * t202;
t292 = t198 * t325;
t380 = pkin(9) + pkin(8);
t293 = qJD(5) * t380;
t259 = -pkin(4) * t199 + pkin(8) * t202;
t137 = t259 * qJD(1);
t204 = -pkin(1) - pkin(2);
t156 = qJD(1) * t204 + qJD(2);
t195 = sin(pkin(10));
t354 = cos(pkin(10));
t275 = t354 * qJD(1);
t112 = qJ(2) * t275 + t195 * t156;
t102 = -qJD(1) * pkin(7) + t112;
t97 = t199 * t102;
t83 = qJD(3) * t202 - t97;
t50 = t198 * t137 + t201 * t83;
t397 = pkin(9) * t292 + t198 * t293 + t50;
t333 = t201 * t202;
t244 = -pkin(5) * t199 + pkin(9) * t333;
t49 = t201 * t137 - t198 * t83;
t396 = -qJD(1) * t244 - t201 * t293 - t49;
t274 = t354 * qJD(2);
t164 = qJD(1) * t274;
t153 = qJDD(1) * t204 + qJDD(2);
t272 = t354 * qJDD(1);
t332 = qJ(2) * t272 + t195 * t153;
t93 = t164 + t332;
t395 = -qJDD(1) * pkin(7) + qJD(3) * qJD(4) + t93;
t336 = t198 * t202;
t120 = -t195 * t336 - t201 * t354;
t277 = t202 * t354;
t288 = t199 * t322;
t356 = qJD(5) * t120 - t195 * t288 - (t195 * t198 + t201 * t277) * qJD(1);
t121 = t195 * t333 - t198 * t354;
t290 = t199 * t313;
t355 = -qJD(5) * t121 + t195 * t290 - (t195 * t201 - t198 * t277) * qJD(1);
t309 = qJD(1) * qJD(4);
t281 = t202 * t309;
t303 = t199 * qJDD(1);
t234 = t281 + t303;
t223 = t234 * t201;
t160 = qJD(5) + t325;
t319 = qJD(5) * t201;
t298 = -t199 * qJDD(3) - t202 * t395;
t323 = qJD(4) * t199;
t36 = -t102 * t323 - t298;
t32 = qJDD(4) * pkin(8) + t36;
t180 = t202 * qJDD(1);
t386 = -t199 * t309 + t180;
t310 = qJD(1) * qJD(2);
t162 = t195 * t310;
t304 = t195 * qJDD(1);
t252 = -qJ(2) * t304 + t153 * t354;
t92 = -t162 + t252;
t87 = qJDD(1) * pkin(3) - t92;
t46 = pkin(4) * t386 + t234 * pkin(8) + t87;
t176 = t195 * qJ(2);
t111 = -qJD(1) * t176 + t156 * t354;
t101 = qJD(1) * pkin(3) - t111;
t260 = t202 * pkin(4) + t199 * pkin(8);
t75 = qJD(1) * t260 + t101;
t301 = t198 * t46 + t201 * t32 + t75 * t319;
t320 = qJD(5) * t198;
t312 = t199 * qJD(3);
t84 = t102 * t202 + t312;
t74 = qJD(4) * pkin(8) + t84;
t240 = t320 * t74 - t301;
t34 = -t198 * t74 + t201 * t75;
t394 = t34 * t160 + t240;
t393 = t248 ^ 2 - t79 ^ 2;
t29 = pkin(9) * t131 + t34;
t27 = pkin(5) * t160 + t29;
t35 = t198 * t75 + t201 * t74;
t30 = pkin(9) * t130 + t35;
t318 = qJD(6) * t197;
t127 = -qJDD(5) - t386;
t315 = t130 * qJD(5);
t68 = t198 * qJDD(4) - t223 + t315;
t45 = t201 * t46;
t9 = -qJD(5) * t35 - t198 * t32 + t45;
t6 = -t127 * pkin(5) - t68 * pkin(9) + t9;
t170 = qJD(5) * t313;
t268 = t201 * qJDD(4) - t170;
t308 = qJD(1) * qJD(5);
t280 = t201 * t308;
t352 = qJD(5) * t74;
t7 = (t199 * t280 + t268) * pkin(9) + (pkin(9) * t234 - t352) * t198 + t301;
t1 = (qJD(6) * t27 + t7) * t200 + t197 * t6 - t30 * t318;
t193 = qJ(5) + qJ(6);
t183 = cos(t193);
t188 = g(3) * t199;
t73 = -qJD(4) * pkin(4) - t83;
t51 = -pkin(5) * t130 + t73;
t377 = sin(qJ(1));
t378 = cos(qJ(1));
t128 = -t195 * t377 - t354 * t378;
t129 = t195 * t378 - t354 * t377;
t182 = sin(t193);
t340 = t183 * t202;
t58 = t128 * t182 + t129 * t340;
t60 = -t128 * t340 + t129 * t182;
t392 = g(1) * t60 - g(2) * t58 - t183 * t188 - t51 * t79 - t1;
t155 = qJD(6) + t160;
t228 = t198 * t303 + t268;
t285 = t199 * t319;
t289 = t202 * t313;
t230 = t285 + t289;
t208 = qJD(1) * t230 + t228;
t317 = qJD(6) * t200;
t236 = t130 * t317 + t131 * t318 + t197 * t208 + t200 * t68;
t391 = -t155 * t79 + t236;
t321 = qJD(4) * t202;
t271 = -t202 * qJDD(3) + t102 * t321 + t199 * t395;
t351 = qJDD(4) * pkin(4);
t33 = t271 - t351;
t389 = qJD(5) * pkin(8) * t160 + t33;
t388 = t198 * t234;
t373 = g(2) * t129;
t257 = g(1) * t128 + t373;
t229 = t257 * t202 + t188;
t191 = t199 ^ 2;
t192 = t202 ^ 2;
t387 = t191 + t192;
t141 = t354 * qJ(2) + t195 * t204;
t133 = -pkin(7) + t141;
t385 = t133 * t321 + t199 * t274;
t287 = t201 * t321;
t384 = -t199 * t320 + t287;
t302 = qJD(5) + qJD(6);
t364 = t200 * t30;
t11 = t197 * t27 + t364;
t2 = -qJD(6) * t11 - t197 * t7 + t200 * t6;
t341 = t182 * t202;
t57 = -t128 * t183 + t129 * t341;
t59 = t128 * t341 + t129 * t183;
t383 = -g(1) * t59 - g(2) * t57 - t182 * t188 - t51 * t248 + t2;
t18 = qJD(6) * t248 + t197 * t68 - t200 * t208;
t382 = t155 * t248 - t18;
t222 = t201 * (t199 * t308 + qJDD(4)) - t170;
t372 = g(3) * t202;
t219 = t199 * t257 - t372;
t233 = qJDD(1) * t198 + t280;
t381 = -t130 * qJD(4) + t202 * t233;
t212 = -qJD(4) * (t199 * t84 + t202 * t83) + t271 * t199 + t36 * t202;
t347 = t129 * t201;
t69 = t128 * t336 + t347;
t379 = g(1) * t69;
t376 = pkin(5) * t198;
t374 = g(2) * t128;
t64 = t120 * t200 - t121 * t197;
t370 = qJD(6) * t64 + t197 * t355 + t356 * t200;
t65 = t120 * t197 + t121 * t200;
t369 = -qJD(6) * t65 - t197 * t356 + t355 * t200;
t150 = t380 * t198;
t151 = t380 * t201;
t95 = -t150 * t200 - t151 * t197;
t368 = qJD(6) * t95 + t197 * t396 - t397 * t200;
t96 = -t150 * t197 + t151 * t200;
t367 = -qJD(6) * t96 + t197 * t397 + t396 * t200;
t366 = t197 * t30;
t365 = t199 * t83;
t363 = t201 * t34;
t362 = t33 * t199;
t360 = t34 * t198;
t359 = t35 * t160;
t135 = t197 * t201 + t198 * t200;
t86 = t302 * t135;
t358 = t135 * t325 + t86;
t334 = t200 * t201;
t338 = t197 * t198;
t357 = -t197 * t292 + t200 * t319 + t201 * t317 - t302 * t338 + t325 * t334;
t140 = t204 * t354 - t176;
t132 = pkin(3) - t140;
t103 = t132 + t260;
t113 = t133 * t333;
t55 = t198 * t103 + t113;
t353 = pkin(1) * qJDD(1);
t350 = t128 * t198;
t349 = t128 * t201;
t348 = t129 * t198;
t346 = t130 * t131;
t345 = t130 * t160;
t344 = t130 * t201;
t343 = t131 * t160;
t342 = t131 * t198;
t339 = t195 * t199;
t337 = t198 * t199;
t335 = t199 * t201;
t331 = t378 * pkin(1) + t377 * qJ(2);
t330 = g(1) * t377 - g(2) * t378;
t329 = t191 - t192;
t205 = qJD(4) ^ 2;
t206 = qJD(1) ^ 2;
t328 = t205 + t206;
t327 = qJD(1) * t101;
t324 = qJD(4) * t131;
t314 = t131 * qJD(5);
t311 = qJ(2) * qJDD(1);
t306 = qJDD(4) * t199;
t305 = qJDD(4) * t202;
t299 = 0.2e1 * t310;
t297 = t197 * t337;
t296 = t199 * t206 * t202;
t295 = t378 * pkin(2) + t331;
t294 = pkin(7) + t376;
t291 = t160 * t322;
t276 = qJD(1) * t55 + t35;
t273 = t68 - t315;
t270 = qJDD(1) * t387;
t269 = qJDD(2) - t353;
t267 = -t128 * pkin(3) + t295;
t266 = t199 * t281;
t265 = pkin(5) * t320 - t312 - (-qJD(1) * t376 + t102) * t202;
t264 = -pkin(1) * t377 + t378 * qJ(2);
t263 = t202 * t274;
t261 = t199 * t275;
t258 = -g(1) * t129 + t374;
t114 = t135 * t199;
t40 = -qJD(6) * t297 + (t302 * t335 + t289) * t200 + t384 * t197;
t256 = t114 * t236 + t248 * t40;
t115 = t199 * t334 - t297;
t39 = t197 * t289 + t199 * t86 - t200 * t287;
t255 = -t115 * t18 - t39 * t79;
t99 = t201 * t103;
t47 = pkin(9) * t335 + t99 + (-t133 * t198 + pkin(5)) * t202;
t48 = pkin(9) * t337 + t55;
t20 = -t197 * t48 + t200 * t47;
t21 = t197 * t47 + t200 * t48;
t124 = -qJDD(6) + t127;
t251 = -t114 * t124 + t155 * t40;
t250 = t115 * t124 + t155 * t39;
t175 = pkin(5) * t201 + pkin(4);
t247 = t175 * t202 + t199 * t380;
t245 = 0.2e1 * qJD(4) * t275;
t243 = t202 * t236 - t248 * t323;
t242 = -t18 * t202 - t323 * t79;
t239 = -t198 * t127 + t160 * t319;
t238 = t201 * t127 + t160 * t320;
t237 = t129 * pkin(7) + t267;
t235 = -t111 * t195 + t112 * t354;
t232 = g(1) * t378 + g(2) * t377;
t231 = -t257 + t327;
t227 = -pkin(2) * t377 + t264;
t226 = t387 * t354;
t225 = pkin(8) * t127 + t160 * t73;
t117 = t195 * qJD(2) + qJD(4) * t259;
t224 = t201 * t117 + t133 * t290 - t198 * t263;
t221 = -t195 * t321 + t261;
t220 = t129 * pkin(3) + t227;
t217 = -t277 * t84 + t354 * t365;
t215 = t222 - t314;
t214 = t128 * pkin(7) + t220;
t213 = -t240 * t201 + (-t198 * t35 - t363) * qJD(5) - t9 * t198;
t25 = t198 * t117 + t201 * t263 + t103 * t319 + (-t202 * t320 - t288) * t133;
t211 = -qJDD(4) * t133 + (-qJD(1) * t132 - t101 - t274) * qJD(4);
t210 = qJDD(1) * t132 - t133 * t205 + t162 + t258 + t87;
t144 = -t199 * t205 + t305;
t143 = -t202 * t205 - t306;
t134 = -t334 + t338;
t100 = (t133 - t376) * t199;
t70 = -t128 * t333 + t348;
t56 = -pkin(5) * t230 + t385;
t54 = -t133 * t336 + t99;
t26 = -qJD(5) * t55 + t224;
t22 = -pkin(5) * t208 + t33;
t19 = pkin(9) * t230 + t25;
t16 = t244 * qJD(4) + (-t113 + (-pkin(9) * t199 - t103) * t198) * qJD(5) + t224;
t13 = t200 * t29 - t366;
t12 = -t197 * t29 - t364;
t10 = t200 * t27 - t366;
t4 = -qJD(6) * t21 + t200 * t16 - t197 * t19;
t3 = qJD(6) * t20 + t197 * t16 + t200 * t19;
t5 = [0, 0, 0, 0, 0, qJDD(1), t330, t232, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -qJDD(2) + t330 + 0.2e1 * t353, 0, -t232 + t299 + 0.2e1 * t311, -t269 * pkin(1) - g(1) * t264 - g(2) * t331 + (t299 + t311) * qJ(2), 0, 0, 0, 0, 0, qJDD(1), -t140 * qJDD(1) + 0.2e1 * t162 - t252 + t258, qJDD(1) * t141 + 0.2e1 * t164 + t257 + t332, 0, -g(1) * t227 - g(2) * t295 + qJD(2) * t235 + t92 * t140 + t93 * t141, qJDD(1) * t191 + 0.2e1 * t266, 0.2e1 * t180 * t199 - 0.2e1 * t309 * t329, t143, qJDD(1) * t192 - 0.2e1 * t266, -t144, 0, t199 * t211 + t202 * t210, -t199 * t210 + t202 * t211, -t133 * t270 - t226 * t310 - t212 - t257, t87 * t132 - g(1) * t214 - g(2) * t237 + (t101 * t195 - t217) * qJD(2) + t212 * t133, t131 * t384 - t68 * t335 (-t342 - t344) * t321 + ((t68 + t315) * t198 + (-t222 - t314 - t388) * t201) * t199 (t68 - t291) * t202 + (t238 + t324) * t199, t130 * t230 + t208 * t337 ((t160 + t325) * t313 + t268) * t202 + (t239 + t381) * t199, -t127 * t202 - t160 * t323, -g(1) * t350 - g(2) * t70 - t54 * t127 + t26 * t160 + (-g(1) * t347 + t9 + (-t130 * t133 - t198 * t73) * qJD(4)) * t202 + (-t130 * t274 - t73 * t319 - t34 * qJD(4) - t33 * t198 + (-t198 * t281 - t199 * t233 - t268) * t133) * t199, -g(1) * t349 - g(2) * t69 + t55 * t127 - t25 * t160 + (g(1) * t348 + t240 + (-t131 * t133 - t201 * t73) * qJD(4)) * t202 + (t35 * qJD(4) - t131 * t274 + t133 * t68 - t33 * t201 + t320 * t73) * t199, t25 * t130 + t55 * t268 + t26 * t131 - t54 * t68 + (t198 * t276 + t363) * t321 + (t9 * t201 + (qJDD(1) * t55 - t240) * t198 + (t201 * t276 - t360) * qJD(5) + t258) * t199, -t240 * t55 + t35 * t25 + t9 * t54 + t34 * t26 + t133 * t362 - g(1) * (t129 * t260 + t214) - g(2) * (-t128 * t260 + t237) + t385 * t73, -t115 * t236 + t248 * t39, -t255 + t256, t243 + t250, -t114 * t18 + t40 * t79, t242 + t251, -t124 * t202 - t155 * t323, -g(1) * t58 - g(2) * t60 - t10 * t323 + t100 * t18 - t22 * t114 - t20 * t124 + t4 * t155 + t2 * t202 - t51 * t40 - t56 * t79, g(1) * t57 - g(2) * t59 - t1 * t202 + t100 * t236 + t11 * t323 - t22 * t115 + t21 * t124 - t3 * t155 + t248 * t56 + t51 * t39, t1 * t114 - t10 * t39 + t11 * t40 + t2 * t115 - t21 * t18 + t199 * t258 - t20 * t236 - t248 * t4 + t3 * t79, t1 * t21 + t11 * t3 + t2 * t20 + t10 * t4 + t22 * t100 + t51 * t56 - g(1) * t220 - g(2) * t267 + (-g(1) * t247 - g(2) * t294) * t129 + (-g(1) * t294 + g(2) * t247) * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t206, -qJ(2) * t206 + t269 - t330, 0, 0, 0, 0, 0, 0, -t195 * t206 - t272, -t206 * t354 + t304, 0, -qJD(1) * t235 + t93 * t195 + t354 * t92 - t330, 0, 0, 0, 0, 0, 0, t199 * t245 - t202 * t272 + (-t202 * t328 - t306) * t195, t202 * t245 + t199 * t272 + (t199 * t328 - t305) * t195, -t195 * t270 + t206 * t226, -t87 * t354 + t217 * qJD(1) + (t212 - t327) * t195 - t330, 0, 0, 0, 0, 0, 0, -t120 * t127 + t221 * t130 + t355 * t160 - t208 * t339, t121 * t127 + t131 * t221 - t160 * t356 + t339 * t68, -t120 * t68 + t121 * t208 + t130 * t356 + t131 * t355, -t73 * t261 + t9 * t120 - t240 * t121 + t356 * t35 + t355 * t34 + (t321 * t73 + t362) * t195 - t330, 0, 0, 0, 0, 0, 0, t79 * t261 - t64 * t124 + (t18 * t199 - t321 * t79) * t195 + t369 * t155, -t248 * t261 + t65 * t124 + (t199 * t236 + t248 * t321) * t195 - t370 * t155, -t18 * t65 - t236 * t64 - t248 * t369 + t370 * t79, -t51 * t261 + t1 * t65 + t2 * t64 + (t199 * t22 + t321 * t51) * t195 + t370 * t11 + t369 * t10 - t330; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) + g(3), 0, 0, 0, 0, 0, 0, t144, t143, 0, t36 * t199 - t271 * t202 + g(3) + (t202 * t84 - t365) * qJD(4), 0, 0, 0, 0, 0, 0 ((-t160 + t325) * t313 + t268) * t202 + (-t239 + t381) * t199 (-t68 - t291) * t202 + (t238 - t324) * t199 (-t342 + t344) * t321 + (t273 * t198 + (t215 + t388) * t201) * t199, g(3) + (-t33 + (t201 * t35 - t360) * qJD(4)) * t202 + (qJD(4) * t73 + t213) * t199, 0, 0, 0, 0, 0, 0, t242 - t251, -t243 + t250, t255 + t256, t1 * t115 - t10 * t40 - t11 * t39 - t114 * t2 - t202 * t22 + t323 * t51 + g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t296, t329 * t206, -t303, t296, -t180, qJDD(4), t84 * qJD(4) + t199 * t231 - t271 + t372, -t188 + (t83 + t97) * qJD(4) + t231 * t202 + t298, 0, 0, t68 * t198 - t201 * t343 (t68 + t345) * t201 + (t314 + (t285 + (t131 + t313) * t202) * qJD(1) + t228) * t198 (-t131 * t199 + t160 * t333) * qJD(1) + t239, t222 * t201 + (t223 - t345) * t198 (t130 * t199 - t160 * t336) * qJD(1) - t238, t160 * t326, t34 * t326 - pkin(4) * t170 + t84 * t130 - t49 * t160 + (pkin(4) * t234 + t225) * t198 + (t351 + t372 + (pkin(4) * t308 - t257) * t199 - t389) * t201, -t35 * t326 - pkin(4) * t68 + t131 * t84 + t160 * t50 + t225 * t201 + (t219 + t389) * t198, -t50 * t130 - t49 * t131 + (pkin(8) * t215 - t394) * t201 + (-t9 - t359 + (t223 + t273) * pkin(8)) * t198 + t229, -t34 * t49 - t35 * t50 - t73 * t84 + (-t33 - t219) * pkin(4) + (t213 + t229) * pkin(8), t135 * t236 + t248 * t357, -t134 * t236 - t135 * t18 - t248 * t358 + t357 * t79, -t135 * t124 + t155 * t357 + t248 * t326, t18 * t134 - t358 * t79, t134 * t124 - t155 * t358 + t326 * t79, t155 * t326, t10 * t326 - t124 * t95 + t134 * t22 + t155 * t367 - t175 * t18 - t183 * t219 - t265 * t79 + t358 * t51, -t11 * t326 + t124 * t96 + t135 * t22 - t155 * t368 - t175 * t236 + t182 * t219 + t248 * t265 + t357 * t51, -t1 * t134 - t10 * t357 - t11 * t358 - t135 * t2 - t18 * t96 - t236 * t95 - t248 * t367 + t368 * t79 + t229, g(3) * t247 + t1 * t96 + t10 * t367 + t11 * t368 - t22 * t175 + t2 * t95 + t265 * t51 - t257 * (t175 * t199 - t202 * t380); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t346, -t130 ^ 2 + t131 ^ 2, t68 - t345, -t346, t208 - t343, -t127, -t379 + t73 * t131 + t359 + t45 + (-t352 + t374) * t201 + (-qJD(5) * t75 - t202 * t373 - t188 - t32) * t198, -t73 * t130 + g(1) * t70 - g(2) * (t129 * t333 + t350) - g(3) * t335 + t394, 0, 0, -t371, t393, t391, t371, t382, -t124, -t12 * t155 + (-t124 * t200 - t131 * t79 - t155 * t318) * pkin(5) + t383, t13 * t155 + (t124 * t197 + t131 * t248 - t155 * t317) * pkin(5) + t392, t10 * t79 + t11 * t248 + t12 * t248 - t13 * t79 + (-t236 * t200 - t18 * t197 + (t197 * t248 + t200 * t79) * qJD(6)) * pkin(5), -t10 * t12 - t11 * t13 + (t1 * t197 + t2 * t200 + t51 * t131 - t379 - g(2) * (t129 * t336 - t349) - g(3) * t337 + (-t10 * t197 + t11 * t200) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t371, t393, t391, t371, t382, -t124, t11 * t155 + t383, t10 * t155 + t392, 0, 0;];
tau_reg  = t5;
