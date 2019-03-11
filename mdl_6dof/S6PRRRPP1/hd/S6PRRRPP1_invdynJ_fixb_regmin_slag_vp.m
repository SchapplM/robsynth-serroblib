% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% 
% Output:
% tau_reg [6x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRRPP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:46:56
% EndTime: 2019-03-08 22:47:10
% DurationCPUTime: 6.09s
% Computational Cost: add. (6335->527), mult. (14543->715), div. (0->0), fcn. (11086->14), ass. (0->267)
t233 = sin(qJ(3));
t236 = cos(qJ(3));
t274 = pkin(3) * t233 - pkin(9) * t236;
t181 = t274 * qJD(3);
t232 = sin(qJ(4));
t234 = sin(qJ(2));
t235 = cos(qJ(4));
t317 = qJD(3) * t233;
t228 = sin(pkin(6));
t324 = qJD(1) * t228;
t237 = cos(qJ(2));
t336 = t236 * t237;
t374 = pkin(8) * t232;
t396 = (-t232 * t336 + t234 * t235) * t324 - t235 * t181 - t317 * t374;
t186 = -pkin(3) * t236 - pkin(9) * t233 - pkin(2);
t312 = qJD(4) * t235;
t339 = t232 * t234;
t395 = -(t235 * t336 + t339) * t324 + t232 * t181 + t186 * t312;
t320 = qJD(2) * t236;
t295 = t232 * t320;
t314 = qJD(4) * t232;
t394 = -t295 + t314;
t337 = t235 * t236;
t209 = pkin(8) * t337;
t376 = pkin(4) * t233;
t266 = -qJ(5) * t337 + t376;
t310 = t235 * qJD(5);
t393 = -t233 * t310 + t266 * qJD(3) + (-t209 + (qJ(5) * t233 - t186) * t232) * qJD(4) - t396;
t338 = t233 * t235;
t392 = -(-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t338 - (-qJD(5) * t233 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t236) * t232 - t395;
t372 = qJ(5) + pkin(9);
t391 = t372 * t233 + pkin(2);
t390 = pkin(4) * t232 + pkin(8);
t206 = -qJD(4) + t320;
t311 = t235 * qJD(3);
t321 = qJD(2) * t233;
t177 = t232 * t321 - t311;
t318 = qJD(3) * t232;
t179 = t235 * t321 + t318;
t227 = sin(pkin(11));
t229 = cos(pkin(11));
t95 = t229 * t177 + t179 * t227;
t389 = t206 * t95;
t173 = t227 * t235 + t229 * t232;
t157 = t173 * qJD(4);
t331 = t173 * t320 - t157;
t342 = t229 * t235;
t330 = t394 * t227 - t229 * t312 + t320 * t342;
t313 = qJD(4) * t233;
t388 = -qJD(2) * t313 + qJDD(3);
t215 = pkin(4) * t235 + pkin(3);
t387 = -t215 * t236 - t391;
t306 = t233 * qJDD(2);
t89 = (qJD(3) * (qJD(4) + t320) + t306) * t232 - t388 * t235;
t268 = -t177 * t227 + t229 * t179;
t386 = t268 ^ 2;
t369 = t392 * t227 + t229 * t393;
t368 = t227 * t393 - t392 * t229;
t283 = qJD(4) * t372;
t150 = -t232 * t283 + t310;
t254 = -t232 * qJD(5) - t235 * t283;
t180 = t274 * qJD(2);
t164 = t235 * t180;
t182 = qJD(2) * pkin(8) + t234 * t324;
t230 = cos(pkin(6));
t323 = qJD(1) * t230;
t385 = -t233 * t182 + t236 * t323;
t67 = qJD(2) * t266 - t232 * t385 + t164;
t333 = t232 * t180 + t235 * t385;
t74 = -qJ(5) * t295 + t333;
t366 = (t254 - t67) * t229 + (-t150 + t74) * t227;
t202 = t233 * t323;
t131 = t236 * t182 + t202;
t384 = t394 * pkin(4) - t131;
t383 = pkin(4) * t89 + qJDD(5);
t121 = qJD(3) * pkin(9) + t131;
t291 = t237 * t324;
t132 = qJD(2) * t186 - t291;
t309 = qJD(1) * qJD(2);
t139 = qJDD(2) * pkin(8) + (qJDD(1) * t234 + t237 * t309) * t228;
t307 = qJDD(1) * t230;
t285 = t233 * t307;
t59 = qJDD(3) * pkin(9) + qJD(3) * t385 + t236 * t139 + t285;
t289 = t234 * t309;
t343 = t228 * t237;
t269 = -qJDD(1) * t343 + t228 * t289;
t81 = qJD(2) * t181 + qJDD(2) * t186 + t269;
t304 = t132 * t312 + t232 * t81 + t235 * t59;
t258 = t121 * t314 - t304;
t11 = -qJ(5) * t89 - qJD(5) * t177 - t258;
t221 = t236 * qJDD(2);
t308 = qJD(2) * qJD(3);
t171 = t233 * t308 + qJDD(4) - t221;
t71 = t121 * t235 + t132 * t232;
t76 = t235 * t81;
t244 = -qJD(4) * t71 - t232 * t59 + t76;
t288 = t236 * t308;
t88 = qJD(4) * t311 + (t288 + t306) * t235 + t388 * t232;
t8 = t171 * pkin(4) - t88 * qJ(5) - t179 * qJD(5) + t244;
t3 = -t227 * t11 + t229 * t8;
t303 = -qJDD(6) + t3;
t120 = -qJD(3) * pkin(3) - t385;
t86 = pkin(4) * t177 + qJD(5) + t120;
t31 = pkin(5) * t95 - qJ(6) * t268 + t86;
t357 = sin(pkin(10));
t277 = t357 * t237;
t358 = cos(pkin(10));
t280 = t358 * t234;
t154 = t230 * t280 + t277;
t282 = t228 * t358;
t107 = t154 * t236 - t233 * t282;
t278 = t357 * t234;
t279 = t358 * t237;
t153 = -t230 * t279 + t278;
t224 = qJ(4) + pkin(11);
t217 = sin(t224);
t218 = cos(t224);
t62 = t107 * t217 - t153 * t218;
t156 = -t230 * t278 + t279;
t281 = t228 * t357;
t109 = t156 * t236 + t233 * t281;
t155 = t230 * t277 + t280;
t64 = t109 * t217 - t155 * t218;
t344 = t228 * t234;
t161 = t230 * t233 + t236 * t344;
t99 = t161 * t217 + t218 * t343;
t382 = g(1) * t64 + g(2) * t62 + g(3) * t99 - t31 * t268 + t303;
t238 = qJD(3) ^ 2;
t273 = g(1) * t155 + g(2) * t153;
t381 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t238 + t228 * (-g(3) * t237 + t289) - t269 + t273;
t253 = g(3) * t343 - t273;
t379 = qJD(4) * (pkin(8) * t206 + t121) - t253;
t4 = t229 * t11 + t227 * t8;
t375 = pkin(5) * t171;
t53 = -qJ(5) * t177 + t71;
t47 = t229 * t53;
t69 = -t121 * t232 + t235 * t132;
t52 = -qJ(5) * t179 + t69;
t23 = t227 * t52 + t47;
t373 = t23 * t268;
t371 = qJ(6) * t317 - qJD(6) * t236 + t368;
t370 = -pkin(5) * t317 - t369;
t39 = -pkin(4) * t206 + t52;
t18 = t227 * t39 + t47;
t367 = -pkin(5) * t331 + qJ(6) * t330 - qJD(6) * t173 + t384;
t30 = t227 * t67 + t229 * t74;
t365 = pkin(5) * t321 - t366;
t27 = qJ(6) * t321 + t30;
t85 = t229 * t150 + t227 * t254;
t364 = t85 - t27;
t363 = qJD(2) * pkin(2);
t362 = t227 * t53;
t361 = t88 * t232;
t175 = t235 * t186;
t104 = -qJ(5) * t338 + t175 + (-pkin(4) - t374) * t236;
t327 = t232 * t186 + t209;
t340 = t232 * t233;
t114 = -qJ(5) * t340 + t327;
t55 = t227 * t104 + t229 * t114;
t106 = t154 * t233 + t236 * t282;
t360 = -t106 * t215 + t107 * t372;
t108 = t156 * t233 - t236 * t281;
t359 = -t108 * t215 + t109 * t372;
t355 = qJDD(3) * pkin(3);
t354 = t107 * t232;
t353 = t109 * t232;
t352 = t153 * t235;
t351 = t155 * t235;
t350 = t177 * t206;
t349 = t179 * t206;
t348 = t179 * t235;
t346 = t217 * t236;
t345 = t218 * t236;
t24 = t229 * t52 - t362;
t335 = qJD(6) - t24;
t334 = qJDD(1) - g(3);
t160 = -t230 * t236 + t233 * t344;
t332 = -t160 * t215 + t161 * t372;
t326 = pkin(4) * t340 + t233 * pkin(8);
t225 = t233 ^ 2;
t325 = -t236 ^ 2 + t225;
t322 = qJD(2) * t228;
t319 = qJD(3) * t177;
t316 = qJD(3) * t236;
t315 = qJD(4) * t206;
t305 = t171 * qJ(6) + t4;
t301 = t235 * t343;
t293 = t232 * t316;
t299 = pkin(4) * t293 + pkin(8) * t316 + t312 * t376;
t297 = t234 * t322;
t296 = t237 * t322;
t294 = t206 * t311;
t292 = t236 * t311;
t290 = t232 * t372;
t287 = t237 * t308;
t36 = t227 * t88 + t229 * t89;
t275 = t233 * t291;
t272 = g(1) * t156 + g(2) * t154;
t271 = -t95 ^ 2 - t386;
t270 = -pkin(5) * t218 - qJ(6) * t217;
t17 = t229 * t39 - t362;
t37 = -t227 * t89 + t229 * t88;
t54 = t104 * t229 - t114 * t227;
t172 = t227 * t232 - t342;
t239 = qJD(2) ^ 2;
t267 = qJDD(2) * t237 - t234 * t239;
t112 = -t161 * t232 - t301;
t263 = -t161 * t235 + t232 * t343;
t261 = t232 * t171 - t206 * t312;
t260 = t235 * t171 + t206 * t314;
t110 = -qJD(3) * t160 + t236 * t296;
t45 = qJD(4) * t263 - t110 * t232 + t235 * t297;
t46 = qJD(4) * t112 + t110 * t235 + t232 * t297;
t20 = t227 * t46 - t229 * t45;
t22 = t227 * t45 + t229 * t46;
t56 = -t229 * t112 - t227 * t263;
t57 = t112 * t227 - t229 * t263;
t259 = t20 * t268 - t22 * t95 - t57 * t36 + t37 * t56;
t257 = g(1) * t108 + g(2) * t106 + g(3) * t160;
t256 = g(1) * t109 + g(2) * t107 + g(3) * t161;
t255 = -qJD(3) * t202 - t233 * t139 - t182 * t316 + t236 * t307;
t252 = -g(3) * t344 - t272;
t251 = -pkin(9) * t171 - t120 * t206;
t250 = pkin(8) * t344 + t391 * t343 + (pkin(4) * t339 + t215 * t336) * t228;
t249 = t153 * t387 + t154 * t390;
t248 = t155 * t387 + t156 * t390;
t246 = t253 * t233;
t60 = -t255 - t355;
t188 = t372 * t235;
t117 = t188 * t227 + t229 * t290;
t118 = t229 * t188 - t227 * t290;
t243 = t117 * t37 - t118 * t36 - t85 * t95 - t256;
t242 = pkin(9) * t315 + t257 - t60;
t183 = -t291 - t363;
t241 = -pkin(8) * qJDD(3) + (t183 + t291 - t363) * qJD(3);
t33 = t60 + t383;
t240 = t255 + t257;
t5 = pkin(5) * t36 - qJ(6) * t37 - qJD(6) * t268 + t33;
t213 = -pkin(4) * t229 - pkin(5);
t210 = pkin(4) * t227 + qJ(6);
t147 = -t227 * t340 + t229 * t338;
t146 = t173 * t233;
t143 = pkin(4) * t351;
t141 = pkin(4) * t352;
t124 = (t217 * t234 + t218 * t336) * t228;
t123 = (t217 * t336 - t218 * t234) * t228;
t111 = qJD(3) * t161 + t233 * t296;
t100 = t161 * t218 - t217 * t343;
t92 = pkin(5) * t172 - qJ(6) * t173 - t215;
t83 = t157 * t233 + t227 * t293 - t229 * t292;
t82 = t172 * t313 - t173 * t316;
t80 = -t155 * t345 + t156 * t217;
t79 = -t155 * t346 - t156 * t218;
t78 = -t153 * t345 + t154 * t217;
t77 = -t153 * t346 - t154 * t218;
t73 = pkin(5) * t146 - qJ(6) * t147 + t326;
t65 = t109 * t218 + t155 * t217;
t63 = t107 * t218 + t153 * t217;
t50 = pkin(5) * t236 - t54;
t49 = -qJ(6) * t236 + t55;
t34 = pkin(4) * t179 + pkin(5) * t268 + qJ(6) * t95;
t26 = -pkin(5) * t82 + qJ(6) * t83 - qJD(6) * t147 + t299;
t15 = -qJ(6) * t206 + t18;
t13 = pkin(5) * t206 + qJD(6) - t17;
t2 = -t303 - t375;
t1 = -qJD(6) * t206 + t305;
t6 = [t334, 0, t267 * t228 (-qJDD(2) * t234 - t237 * t239) * t228, 0, 0, 0, 0, 0, -t111 * qJD(3) - t160 * qJDD(3) + (-t233 * t287 + t236 * t267) * t228, -t110 * qJD(3) - t161 * qJDD(3) + (-t233 * t267 - t236 * t287) * t228, 0, 0, 0, 0, 0, t111 * t177 + t112 * t171 + t160 * t89 - t206 * t45, t111 * t179 + t160 * t88 + t171 * t263 + t206 * t46, t259, t111 * t86 + t160 * t33 - t17 * t20 + t18 * t22 - t3 * t56 + t4 * t57 - g(3), t111 * t95 + t160 * t36 - t171 * t56 + t20 * t206, t259, -t111 * t268 - t160 * t37 + t171 * t57 - t206 * t22, t1 * t57 + t111 * t31 + t13 * t20 + t15 * t22 + t160 * t5 + t2 * t56 - g(3); 0, qJDD(2), t334 * t343 + t273, -t334 * t344 + t272, qJDD(2) * t225 + 0.2e1 * t233 * t288, 0.2e1 * t221 * t233 - 0.2e1 * t308 * t325, qJDD(3) * t233 + t236 * t238, qJDD(3) * t236 - t233 * t238, 0, t241 * t233 + t381 * t236, -t381 * t233 + t241 * t236, t88 * t338 + (-t232 * t313 + t292) * t179 (-t177 * t235 - t179 * t232) * t316 + (-t361 - t235 * t89 + (t177 * t232 - t348) * qJD(4)) * t233 (-t88 - t294) * t236 + (qJD(3) * t179 + t260) * t233 (t206 * t318 + t89) * t236 + (-t261 - t319) * t233, -t171 * t236 - t206 * t317, t175 * t171 + t396 * t206 + (t186 * t315 + t252) * t232 + (pkin(8) * t319 - t76 + (-pkin(8) * t171 + qJD(3) * t120 + qJD(4) * t132 + t59) * t232 + t379 * t235) * t236 + (pkin(8) * t89 + qJD(3) * t69 + t120 * t312 - t177 * t291 + t60 * t232) * t233, -t327 * t171 + t395 * t206 + t252 * t235 + ((pkin(8) * t179 + t120 * t235) * qJD(3) - t379 * t232 + t304) * t236 + (-t179 * t291 - t120 * t314 - t71 * qJD(3) + t60 * t235 + (t88 - t294) * pkin(8)) * t233, -t146 * t4 - t147 * t3 + t17 * t83 + t18 * t82 - t268 * t369 - t36 * t55 - t368 * t95 - t37 * t54 - t246, t4 * t55 + t3 * t54 + t33 * t326 - g(1) * t248 - g(2) * t249 - g(3) * t250 + (-t275 + t299) * t86 + t368 * t18 + t369 * t17, -g(1) * t80 - g(2) * t78 - g(3) * t124 + t146 * t5 - t171 * t50 + t2 * t236 + t26 * t95 - t31 * t82 + t36 * t73 + (-qJD(3) * t13 - t291 * t95) * t233 + t370 * t206, -t1 * t146 - t13 * t83 + t147 * t2 + t15 * t82 + t268 * t370 - t36 * t49 + t37 * t50 - t371 * t95 - t246, -g(1) * t79 - g(2) * t77 - g(3) * t123 - t1 * t236 - t147 * t5 + t171 * t49 - t26 * t268 + t31 * t83 - t37 * t73 + (qJD(3) * t15 + t268 * t291) * t233 - t371 * t206, t1 * t49 + t5 * t73 + t2 * t50 - g(1) * (pkin(5) * t80 + qJ(6) * t79 + t248) - g(2) * (pkin(5) * t78 + qJ(6) * t77 + t249) - g(3) * (pkin(5) * t124 + qJ(6) * t123 + t250) + (t26 - t275) * t31 + t371 * t15 + t370 * t13; 0, 0, 0, 0, -t233 * t239 * t236, t325 * t239, t306, t221, qJDD(3), qJD(3) * t131 - t183 * t321 + t240, -t285 + (-qJD(2) * t183 - t139) * t236 + t256, -t206 * t348 + t361 (t88 + t350) * t235 + (-t89 + t349) * t232 (-t179 * t233 + t206 * t337) * qJD(2) + t261 (-t206 * t232 * t236 + t177 * t233) * qJD(2) + t260, t206 * t321, -t69 * t321 - pkin(3) * t89 - t131 * t177 + t164 * t206 + (-t206 * t385 + t251) * t232 + t242 * t235, -pkin(3) * t88 - t131 * t179 - t206 * t333 - t232 * t242 + t235 * t251 + t321 * t71, t17 * t330 - t172 * t4 - t173 * t3 + t18 * t331 - t268 * t366 + t30 * t95 + t243, t4 * t118 - t3 * t117 - t33 * t215 - g(1) * t359 - g(2) * t360 - g(3) * t332 + t384 * t86 + (t85 - t30) * t18 + t366 * t17, -t117 * t171 + t13 * t321 + t5 * t172 + t206 * t365 + t218 * t257 - t31 * t331 + t92 * t36 + t367 * t95, -t1 * t172 - t13 * t330 + t15 * t331 + t173 * t2 + t268 * t365 + t27 * t95 + t243, t118 * t171 - t15 * t321 - t173 * t5 - t206 * t364 + t217 * t257 - t268 * t367 + t31 * t330 - t37 * t92, t1 * t118 + t5 * t92 + t2 * t117 - g(1) * (t108 * t270 + t359) - g(2) * (t106 * t270 + t360) - g(3) * (t160 * t270 + t332) + t367 * t31 + t364 * t15 + t365 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179 * t177, -t177 ^ 2 + t179 ^ 2, t88 - t350, -t349 - t89, t171, -t71 * t206 - t120 * t179 - g(1) * (t351 - t353) - g(2) * (t352 - t354) - g(3) * t112 + t244, -t69 * t206 + t120 * t177 - g(1) * (-t109 * t235 - t155 * t232) - g(2) * (-t107 * t235 - t153 * t232) - g(3) * t263 + t258, t18 * t268 - t373 + (-t227 * t36 - t229 * t37) * pkin(4) + (-t17 + t24) * t95, -g(1) * t143 - g(2) * t141 + t17 * t23 - t18 * t24 + (g(3) * t301 - t86 * t179 + t4 * t227 + t3 * t229 + t232 * t256) * pkin(4), -t23 * t206 - t34 * t95 + (pkin(5) - t213) * t171 + t382, t15 * t268 - t210 * t36 + t213 * t37 - t373 + (t13 - t335) * t95, -g(1) * t65 - g(2) * t63 - g(3) * t100 + t210 * t171 - t31 * t95 + t34 * t268 + (-0.2e1 * qJD(6) + t24) * t206 + t305, t1 * t210 + t2 * t213 - t31 * t34 - t13 * t23 - g(1) * (-pkin(4) * t353 - pkin(5) * t64 + qJ(6) * t65 + t143) - g(2) * (-pkin(4) * t354 - pkin(5) * t62 + qJ(6) * t63 + t141) - g(3) * (pkin(4) * t112 - t99 * pkin(5) + t100 * qJ(6)) + t335 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t271, t17 * t268 + t18 * t95 - t240 - t355 + t383, -t206 * t268 + t36, t271, -t37 - t389, -t13 * t268 + t15 * t95 - t257 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t268 * t95 - t171, t37 - t389, -t206 ^ 2 - t386, t15 * t206 - t375 - t382;];
tau_reg  = t6;
