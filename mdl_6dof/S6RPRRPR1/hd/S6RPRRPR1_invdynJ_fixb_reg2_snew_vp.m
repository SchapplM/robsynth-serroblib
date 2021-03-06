% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRPR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:56:47
% EndTime: 2019-05-05 21:57:00
% DurationCPUTime: 6.55s
% Computational Cost: add. (40264->498), mult. (87174->731), div. (0->0), fcn. (62601->12), ass. (0->301)
t374 = 2 * qJD(5);
t288 = sin(pkin(11));
t294 = sin(qJ(4));
t297 = cos(qJ(4));
t298 = cos(qJ(3));
t295 = sin(qJ(3));
t332 = qJD(1) * t295;
t254 = t297 * t298 * qJD(1) - t294 * t332;
t255 = (t298 * t294 + t295 * t297) * qJD(1);
t290 = cos(pkin(11));
t226 = -t290 * t254 + t288 * t255;
t228 = t288 * t254 + t290 * t255;
t191 = t228 * t226;
t283 = qJDD(3) + qJDD(4);
t366 = -t191 + t283;
t373 = t288 * t366;
t372 = t290 * t366;
t293 = sin(qJ(6));
t327 = qJD(1) * qJD(3);
t320 = t298 * t327;
t326 = t295 * qJDD(1);
t260 = t320 + t326;
t279 = t298 * qJDD(1);
t321 = t295 * t327;
t309 = t279 - t321;
t312 = t294 * t260 - t297 * t309;
t211 = -t255 * qJD(4) - t312;
t212 = t254 * qJD(4) + t297 * t260 + t294 * t309;
t314 = -t290 * t211 + t288 * t212;
t169 = qJDD(6) + t314;
t284 = qJD(3) + qJD(4);
t296 = cos(qJ(6));
t206 = t293 * t228 - t296 * t284;
t208 = t296 * t228 + t293 * t284;
t174 = t208 * t206;
t367 = t169 - t174;
t371 = t293 * t367;
t236 = t254 * t255;
t365 = t236 + t283;
t370 = t294 * t365;
t369 = t296 * t367;
t368 = t297 * t365;
t353 = t284 * t228;
t150 = t314 + t353;
t185 = t226 * pkin(5) - t228 * pkin(9);
t363 = t284 ^ 2;
t300 = qJD(1) ^ 2;
t360 = sin(qJ(1));
t361 = cos(qJ(1));
t305 = t361 * g(1) + t360 * g(2);
t258 = -t300 * pkin(1) - t305;
t289 = sin(pkin(10));
t291 = cos(pkin(10));
t304 = t360 * g(1) - t361 * g(2);
t302 = qJDD(1) * pkin(1) + t304;
t333 = t291 * t258 + t289 * t302;
t230 = -t300 * pkin(2) + qJDD(1) * pkin(7) + t333;
t286 = -g(3) + qJDD(2);
t209 = t295 * t230 - t298 * t286;
t335 = t298 * t300;
t271 = t295 * t335;
t264 = qJDD(3) + t271;
t186 = (-t260 + t320) * pkin(8) + t264 * pkin(3) - t209;
t266 = qJD(3) * pkin(3) - pkin(8) * t332;
t340 = t295 * t286;
t188 = t340 + t309 * pkin(8) - qJD(3) * t266 + (-pkin(3) * t335 + t230) * t298;
t142 = t294 * t186 + t297 * t188;
t240 = t284 * pkin(4) - t255 * qJ(5);
t249 = t254 ^ 2;
t120 = -t249 * pkin(4) + t211 * qJ(5) - t284 * t240 + t142;
t141 = -t297 * t186 + t294 * t188;
t244 = t284 * t254;
t198 = -t244 + t212;
t301 = pkin(4) * t365 - qJ(5) * t198 - t141;
t68 = -0.2e1 * qJD(5) * t226 + t290 * t120 + t288 * t301;
t50 = -pkin(5) * t363 + t283 * pkin(9) - t226 * t185 + t68;
t313 = -t289 * t258 + t291 * t302;
t229 = -qJDD(1) * pkin(2) - t300 * pkin(7) - t313;
t362 = t298 ^ 2;
t281 = t362 * t300;
t192 = -t309 * pkin(3) - pkin(8) * t281 + t266 * t332 + t229;
t134 = -t211 * pkin(4) - t249 * qJ(5) + t255 * t240 + qJDD(5) + t192;
t172 = t288 * t211 + t290 * t212;
t218 = t284 * t226;
t154 = t172 - t218;
t81 = t150 * pkin(5) - pkin(9) * t154 + t134;
t39 = t293 * t50 - t296 * t81;
t40 = t293 * t81 + t296 * t50;
t22 = t293 * t39 + t296 * t40;
t364 = t244 + t212;
t317 = t288 * t120 - t290 * t301;
t67 = t228 * t374 + t317;
t223 = qJD(6) + t226;
t315 = t293 * t172 - t296 * t283;
t123 = (qJD(6) - t223) * t208 + t315;
t203 = t206 ^ 2;
t204 = t208 ^ 2;
t222 = t223 ^ 2;
t224 = t226 ^ 2;
t225 = t228 ^ 2;
t250 = t255 ^ 2;
t359 = pkin(5) * t288;
t49 = -t283 * pkin(5) - t363 * pkin(9) + (t374 + t185) * t228 + t317;
t46 = t293 * t49;
t42 = t288 * t68 - t290 * t67;
t358 = t294 * t42;
t96 = -t297 * t141 + t294 * t142;
t357 = t295 * t96;
t47 = t296 * t49;
t356 = t297 * t42;
t355 = t223 * t293;
t354 = t223 * t296;
t352 = t284 * t288;
t351 = t284 * t290;
t350 = t284 * t294;
t349 = t284 * t297;
t348 = t288 * t134;
t183 = t191 + t283;
t347 = t288 * t183;
t346 = t290 * t134;
t345 = t290 * t183;
t130 = t169 + t174;
t344 = t293 * t130;
t343 = t294 * t192;
t233 = -t236 + t283;
t342 = t294 * t233;
t341 = t295 * t264;
t339 = t296 * t130;
t338 = t297 * t192;
t337 = t297 * t233;
t265 = qJDD(3) - t271;
t336 = t298 * t265;
t328 = qJD(6) + t223;
t14 = t288 * t22 - t290 * t49;
t325 = pkin(4) * t14 - pkin(5) * t49 + pkin(9) * t22;
t324 = t288 * t174;
t323 = t290 * t174;
t322 = -pkin(5) * t290 - pkin(4);
t43 = t288 * t67 + t290 * t68;
t308 = -t296 * t172 - t293 * t283;
t128 = t328 * t206 + t308;
t167 = -t204 - t222;
t95 = -t293 * t167 - t339;
t64 = t290 * t128 + t288 * t95;
t319 = pkin(4) * t64 + pkin(5) * t128 + pkin(9) * t95 + t46;
t124 = -t328 * t208 - t315;
t159 = -t222 - t203;
t90 = t296 * t159 - t371;
t59 = t290 * t124 + t288 * t90;
t318 = pkin(4) * t59 + pkin(5) * t124 + pkin(9) * t90 - t47;
t97 = t294 * t141 + t297 * t142;
t210 = t298 * t230 + t340;
t170 = t295 * t209 + t298 * t210;
t214 = -t225 - t363;
t161 = t290 * t214 - t347;
t311 = pkin(4) * t161 - t68;
t149 = t203 + t204;
t140 = -t206 * qJD(6) - t308;
t179 = t223 * t206;
t127 = t140 + t179;
t78 = -t123 * t296 + t293 * t127;
t54 = t290 * t149 + t288 * t78;
t310 = pkin(4) * t54 + pkin(5) * t149 + pkin(9) * t78 + t22;
t261 = t279 - 0.2e1 * t321;
t21 = t293 * t40 - t296 * t39;
t181 = -t363 - t224;
t136 = t288 * t181 + t372;
t307 = pkin(4) * t136 - t67;
t306 = -t314 + t353;
t303 = (-qJD(4) + t284) * t255 - t312;
t299 = qJD(3) ^ 2;
t285 = t295 ^ 2;
t280 = t285 * t300;
t269 = -t281 - t299;
t268 = -t280 - t299;
t263 = t280 + t281;
t262 = (t285 + t362) * qJDD(1);
t259 = 0.2e1 * t320 + t326;
t242 = -t250 + t363;
t241 = t249 - t363;
t239 = -t250 - t363;
t238 = -t295 * t268 - t336;
t237 = t298 * t269 - t341;
t235 = t250 - t249;
t231 = -t363 - t249;
t216 = -t225 + t363;
t215 = t224 - t363;
t213 = -t249 - t250;
t200 = -t294 * t239 - t337;
t199 = t297 * t239 - t342;
t193 = (qJD(4) + t284) * t255 + t312;
t190 = t297 * t231 - t370;
t189 = t294 * t231 + t368;
t187 = t225 - t224;
t178 = -t204 + t222;
t177 = t203 - t222;
t176 = (-t226 * t290 + t228 * t288) * t284;
t175 = (-t226 * t288 - t228 * t290) * t284;
t173 = t204 - t203;
t168 = -t224 - t225;
t166 = t290 * t215 - t347;
t165 = -t288 * t216 + t372;
t164 = t288 * t215 + t345;
t163 = t290 * t216 + t373;
t162 = -t288 * t214 - t345;
t158 = -t295 * t199 + t298 * t200;
t157 = t294 * t198 + t297 * t303;
t156 = -t297 * t198 + t294 * t303;
t155 = t172 + t218;
t147 = t290 * t172 - t228 * t352;
t146 = t288 * t172 + t228 * t351;
t145 = t226 * t351 + t288 * t314;
t144 = t226 * t352 - t290 * t314;
t143 = -t295 * t189 + t298 * t190;
t139 = -t208 * qJD(6) - t315;
t137 = t290 * t181 - t373;
t133 = (-t206 * t296 + t208 * t293) * t223;
t132 = (-t206 * t293 - t208 * t296) * t223;
t126 = t140 - t179;
t119 = t296 * t140 - t208 * t355;
t118 = t293 * t140 + t208 * t354;
t117 = -t293 * t139 + t206 * t354;
t116 = t296 * t139 + t206 * t355;
t112 = -t294 * t161 + t297 * t162;
t111 = t297 * t161 + t294 * t162;
t110 = -t295 * t156 + t298 * t157;
t109 = t290 * t133 + t288 * t169;
t108 = t288 * t133 - t290 * t169;
t107 = t288 * t155 + t290 * t306;
t106 = -t290 * t150 - t288 * t154;
t105 = -t290 * t155 + t288 * t306;
t104 = -t288 * t150 + t290 * t154;
t103 = pkin(4) * t105;
t102 = t296 * t177 - t344;
t101 = -t293 * t178 + t369;
t100 = t293 * t177 + t339;
t99 = t296 * t178 + t371;
t98 = -qJ(5) * t161 + t346;
t94 = t296 * t167 - t344;
t92 = -t294 * t136 + t297 * t137;
t91 = t297 * t136 + t294 * t137;
t89 = t293 * t159 + t369;
t87 = -qJ(5) * t136 + t348;
t86 = t290 * t119 + t324;
t85 = t290 * t117 - t324;
t84 = t288 * t119 - t323;
t83 = t288 * t117 + t323;
t82 = -pkin(4) * t154 + qJ(5) * t162 + t348;
t79 = -pkin(4) * t150 + qJ(5) * t137 - t346;
t77 = t296 * t124 - t293 * t126;
t76 = -t123 * t293 - t296 * t127;
t75 = t293 * t124 + t296 * t126;
t73 = t290 * t102 - t288 * t123;
t72 = t290 * t101 + t288 * t127;
t71 = t288 * t102 + t290 * t123;
t70 = t288 * t101 - t290 * t127;
t69 = -t295 * t111 + t298 * t112;
t65 = -t288 * t128 + t290 * t95;
t62 = -t294 * t105 + t297 * t107;
t61 = t297 * t105 + t294 * t107;
t60 = -t288 * t124 + t290 * t90;
t57 = t288 * t173 + t290 * t77;
t56 = -t290 * t173 + t288 * t77;
t55 = -t288 * t149 + t290 * t78;
t52 = t298 * t97 - t357;
t51 = -t295 * t91 + t298 * t92;
t45 = -pkin(9) * t94 + t47;
t44 = -pkin(9) * t89 + t46;
t41 = pkin(4) * t42;
t36 = -t294 * t64 + t297 * t65;
t35 = t294 * t65 + t297 * t64;
t34 = -t295 * t61 + t298 * t62;
t33 = -t294 * t59 + t297 * t60;
t32 = t294 * t60 + t297 * t59;
t31 = -pkin(4) * t134 + qJ(5) * t43;
t30 = -t294 * t54 + t297 * t55;
t29 = t294 * t55 + t297 * t54;
t28 = -qJ(5) * t105 - t42;
t27 = -pkin(5) * t94 + t40;
t26 = -pkin(5) * t89 + t39;
t25 = -pkin(4) * t168 + qJ(5) * t107 + t43;
t24 = t297 * t43 - t358;
t23 = t294 * t43 + t356;
t19 = -t295 * t35 + t298 * t36;
t18 = -t295 * t32 + t298 * t33;
t17 = -t295 * t29 + t298 * t30;
t16 = -pkin(9) * t76 - t21;
t15 = t290 * t22 + t288 * t49;
t12 = -qJ(5) * t64 - t288 * t27 + t290 * t45;
t11 = -qJ(5) * t59 - t288 * t26 + t290 * t44;
t10 = -pkin(4) * t94 + qJ(5) * t65 + t290 * t27 + t288 * t45;
t9 = -pkin(4) * t89 + qJ(5) * t60 + t290 * t26 + t288 * t44;
t8 = -qJ(5) * t54 + t290 * t16 + t76 * t359;
t7 = qJ(5) * t55 + t288 * t16 + t322 * t76;
t6 = -t295 * t23 + t298 * t24;
t5 = -t294 * t14 + t297 * t15;
t4 = t297 * t14 + t294 * t15;
t3 = -qJ(5) * t14 + (-pkin(9) * t290 + t359) * t21;
t2 = qJ(5) * t15 + (-pkin(9) * t288 + t322) * t21;
t1 = -t295 * t4 + t298 * t5;
t13 = [0, 0, 0, 0, 0, qJDD(1), t304, t305, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t291 * qJDD(1) - t289 * t300) + t313, pkin(1) * (-t289 * qJDD(1) - t291 * t300) - t333, 0, pkin(1) * (t289 * t333 + t291 * t313), (t260 + t320) * t295, t298 * t259 + t295 * t261, t341 + t298 * (-t280 + t299), t261 * t298, t295 * (t281 - t299) + t336, 0, -t298 * t229 + pkin(2) * t261 + pkin(7) * t237 + pkin(1) * (t289 * t237 + t291 * t261), t295 * t229 - pkin(2) * t259 + pkin(7) * t238 + pkin(1) * (t289 * t238 - t291 * t259), pkin(2) * t263 + pkin(7) * t262 + pkin(1) * (t289 * t262 + t291 * t263) + t170, -pkin(2) * t229 + pkin(7) * t170 + pkin(1) * (t289 * t170 - t291 * t229), t295 * (t297 * t212 - t255 * t350) + t298 * (t294 * t212 + t255 * t349), t295 * (-t297 * t193 - t294 * t364) + t298 * (-t294 * t193 + t297 * t364), t295 * (-t294 * t242 + t368) + t298 * (t297 * t242 + t370), t295 * (-t294 * t211 - t254 * t349) + t298 * (t297 * t211 - t254 * t350), t295 * (t297 * t241 - t342) + t298 * (t294 * t241 + t337), (t295 * (t254 * t297 + t255 * t294) + t298 * (t254 * t294 - t255 * t297)) * t284, t295 * (-pkin(8) * t189 + t343) + t298 * (-pkin(3) * t193 + pkin(8) * t190 - t338) - pkin(2) * t193 + pkin(7) * t143 + pkin(1) * (t289 * t143 - t291 * t193), t295 * (-pkin(8) * t199 + t338) + t298 * (-pkin(3) * t364 + pkin(8) * t200 + t343) - pkin(2) * t364 + pkin(7) * t158 + pkin(1) * (t289 * t158 - t291 * t364), t295 * (-pkin(8) * t156 - t96) + t298 * (-pkin(3) * t213 + pkin(8) * t157 + t97) - pkin(2) * t213 + pkin(7) * t110 + pkin(1) * (t289 * t110 - t291 * t213), -pkin(8) * t357 + t298 * (-pkin(3) * t192 + pkin(8) * t97) - pkin(2) * t192 + pkin(7) * t52 + pkin(1) * (-t291 * t192 + t289 * t52), t295 * (-t294 * t146 + t297 * t147) + t298 * (t297 * t146 + t294 * t147), t295 * (-t294 * t104 + t297 * t106) + t298 * (t297 * t104 + t294 * t106), t295 * (-t294 * t163 + t297 * t165) + t298 * (t297 * t163 + t294 * t165), t295 * (-t294 * t144 + t297 * t145) + t298 * (t297 * t144 + t294 * t145), t295 * (-t294 * t164 + t297 * t166) + t298 * (t297 * t164 + t294 * t166), t295 * (-t294 * t175 + t297 * t176) + t298 * (t297 * t175 + t294 * t176), t295 * (-pkin(8) * t91 - t294 * t79 + t297 * t87) + t298 * (-pkin(3) * t150 + pkin(8) * t92 + t294 * t87 + t297 * t79) - pkin(2) * t150 + pkin(7) * t51 + pkin(1) * (-t291 * t150 + t289 * t51), t295 * (-pkin(8) * t111 - t294 * t82 + t297 * t98) + t298 * (-pkin(3) * t154 + pkin(8) * t112 + t294 * t98 + t297 * t82) - pkin(2) * t154 + pkin(7) * t69 + pkin(1) * (-t154 * t291 + t289 * t69), t295 * (-pkin(8) * t61 - t294 * t25 + t297 * t28) + t298 * (-pkin(3) * t168 + pkin(8) * t62 + t297 * t25 + t294 * t28) - pkin(2) * t168 + pkin(7) * t34 + pkin(1) * (-t291 * t168 + t289 * t34), t295 * (-pkin(8) * t23 - qJ(5) * t356 - t294 * t31) + t298 * (-pkin(3) * t134 + pkin(8) * t24 - qJ(5) * t358 + t297 * t31) - pkin(2) * t134 + pkin(7) * t6 + pkin(1) * (-t291 * t134 + t289 * t6), t295 * (-t294 * t84 + t297 * t86) + t298 * (t294 * t86 + t297 * t84), t295 * (-t294 * t56 + t297 * t57) + t298 * (t294 * t57 + t297 * t56), t295 * (-t294 * t70 + t297 * t72) + t298 * (t294 * t72 + t297 * t70), t295 * (-t294 * t83 + t297 * t85) + t298 * (t294 * t85 + t297 * t83), t295 * (-t294 * t71 + t297 * t73) + t298 * (t294 * t73 + t297 * t71), t295 * (-t294 * t108 + t297 * t109) + t298 * (t297 * t108 + t294 * t109), t295 * (-pkin(8) * t32 + t297 * t11 - t294 * t9) + t298 * (-pkin(3) * t89 + pkin(8) * t33 + t294 * t11 + t297 * t9) - pkin(2) * t89 + pkin(7) * t18 + pkin(1) * (t289 * t18 - t291 * t89), t295 * (-pkin(8) * t35 - t294 * t10 + t297 * t12) + t298 * (-pkin(3) * t94 + pkin(8) * t36 + t297 * t10 + t294 * t12) - pkin(2) * t94 + pkin(7) * t19 + pkin(1) * (t289 * t19 - t291 * t94), t295 * (-pkin(8) * t29 - t294 * t7 + t297 * t8) + t298 * (-pkin(3) * t76 + pkin(8) * t30 + t294 * t8 + t297 * t7) - pkin(2) * t76 + pkin(7) * t17 + pkin(1) * (t289 * t17 - t291 * t76), t295 * (-pkin(8) * t4 - t294 * t2 + t297 * t3) + t298 * (-pkin(3) * t21 + pkin(8) * t5 + t297 * t2 + t294 * t3) - pkin(2) * t21 + pkin(7) * t1 + pkin(1) * (t289 * t1 - t291 * t21); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t286, 0, 0, 0, 0, 0, 0, t298 * t264 + t295 * t269, -t295 * t265 + t298 * t268, 0, -t298 * t209 + t295 * t210, 0, 0, 0, 0, 0, 0, t298 * t189 + t295 * t190, t298 * t199 + t295 * t200, t298 * t156 + t295 * t157, t295 * t97 + t298 * t96, 0, 0, 0, 0, 0, 0, t295 * t92 + t298 * t91, t298 * t111 + t295 * t112, t295 * t62 + t298 * t61, t298 * t23 + t295 * t24, 0, 0, 0, 0, 0, 0, t295 * t33 + t298 * t32, t295 * t36 + t298 * t35, t298 * t29 + t295 * t30, t295 * t5 + t298 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t271, t280 - t281, t326, t271, t279, qJDD(3), -t209, -t210, 0, 0, -t236, t235, t198, t236, t303, t283, pkin(3) * t189 - t141, pkin(3) * t199 - t142, pkin(3) * t156, pkin(3) * t96, t191, t187, t155, -t191, t306, t283, pkin(3) * t91 + t307, pkin(3) * t111 + t311, pkin(3) * t61 + t103, pkin(3) * t23 + t41, t118, t75, t99, t116, t100, t132, pkin(3) * t32 + t318, pkin(3) * t35 + t319, pkin(3) * t29 + t310, pkin(3) * t4 + t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t236, t235, t198, t236, t303, t283, -t141, -t142, 0, 0, t191, t187, t155, -t191, t306, t283, t307, t311, t103, t41, t118, t75, t99, t116, t100, t132, t318, t319, t310, t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, t154, t168, t134, 0, 0, 0, 0, 0, 0, t89, t94, t76, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, t173, t127, -t174, -t123, t169, -t39, -t40, 0, 0;];
tauJ_reg  = t13;
