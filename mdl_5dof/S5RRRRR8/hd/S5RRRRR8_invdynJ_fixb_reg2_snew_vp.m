% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:25:50
% EndTime: 2019-12-31 22:26:04
% DurationCPUTime: 5.79s
% Computational Cost: add. (33632->455), mult. (67966->639), div. (0->0), fcn. (49108->10), ass. (0->288)
t262 = sin(qJ(5));
t264 = sin(qJ(3));
t265 = sin(qJ(2));
t269 = cos(qJ(3));
t270 = cos(qJ(2));
t238 = (t270 * t264 + t265 * t269) * qJD(1);
t259 = qJD(2) + qJD(3);
t263 = sin(qJ(4));
t268 = cos(qJ(4));
t221 = t238 * t263 - t268 * t259;
t223 = t238 * t268 + t259 * t263;
t267 = cos(qJ(5));
t193 = t267 * t221 + t223 * t262;
t195 = -t221 * t262 + t223 * t267;
t160 = t195 * t193;
t254 = t265 * qJDD(1);
t310 = qJD(1) * qJD(2);
t299 = t270 * t310;
t243 = t254 + t299;
t255 = t270 * qJDD(1);
t300 = t265 * t310;
t244 = t255 - t300;
t291 = t243 * t264 - t269 * t244;
t202 = -qJD(3) * t238 - t291;
t201 = qJDD(4) - t202;
t200 = qJDD(5) + t201;
t352 = -t160 + t200;
t358 = t262 * t352;
t199 = t223 * t221;
t350 = -t199 + t201;
t357 = t263 * t350;
t315 = qJD(1) * t265;
t236 = -t269 * t270 * qJD(1) + t264 * t315;
t216 = t238 * t236;
t309 = qJDD(2) + qJDD(3);
t349 = -t216 + t309;
t356 = t264 * t349;
t355 = t267 * t352;
t354 = t268 * t350;
t353 = t269 * t349;
t287 = t243 * t269 + t244 * t264;
t203 = -qJD(3) * t236 + t287;
t229 = t259 * t236;
t186 = t203 - t229;
t261 = t270 ^ 2;
t272 = qJD(1) ^ 2;
t266 = sin(qJ(1));
t343 = cos(qJ(1));
t298 = g(1) * t266 - t343 * g(2);
t282 = qJDD(1) * pkin(1) + t298;
t283 = qJD(2) * pkin(2) - pkin(7) * t315;
t205 = pkin(2) * t244 + (pkin(7) * t261 + pkin(6)) * t272 - t283 * t315 + t282;
t120 = -t186 * pkin(8) + (t238 * t259 - t202) * pkin(3) - t205;
t284 = t343 * g(1) + t266 * g(2);
t335 = qJDD(1) * pkin(6);
t275 = -t272 * pkin(1) - t284 + t335;
t226 = -t265 * g(3) + t270 * t275;
t257 = t261 * t272;
t198 = -pkin(2) * t257 + t244 * pkin(7) - qJD(2) * t283 + t226;
t317 = t269 * t198;
t274 = t265 * t275;
t318 = t265 * t272;
t341 = t243 * pkin(7);
t348 = qJDD(2) * pkin(2) - t274 + (pkin(2) * t318 + pkin(7) * t310 - g(3)) * t270 - t341;
t163 = t348 * t264 + t317;
t214 = pkin(3) * t236 - pkin(8) * t238;
t346 = t259 ^ 2;
t131 = -t346 * pkin(3) + t309 * pkin(8) - t236 * t214 + t163;
t75 = -t268 * t120 + t131 * t263;
t76 = t263 * t120 + t268 * t131;
t45 = t263 * t75 + t268 * t76;
t279 = -t268 * t203 - t263 * t309;
t172 = -t221 * qJD(4) - t279;
t293 = t203 * t263 - t268 * t309;
t281 = qJD(4) * t223 + t293;
t115 = -t193 * qJD(5) + t267 * t172 - t262 * t281;
t233 = qJD(4) + t236;
t231 = qJD(5) + t233;
t179 = t231 * t193;
t351 = -t179 + t115;
t209 = t233 * t221;
t157 = t209 + t172;
t294 = t172 * t262 + t267 * t281;
t92 = (qJD(5) - t231) * t195 + t294;
t153 = (qJD(4) - t233) * t223 + t293;
t191 = t193 ^ 2;
t192 = t195 ^ 2;
t347 = t221 ^ 2;
t219 = t223 ^ 2;
t230 = t231 ^ 2;
t232 = t233 ^ 2;
t234 = t236 ^ 2;
t235 = t238 ^ 2;
t62 = pkin(4) * t350 - t157 * pkin(9) - t75;
t206 = pkin(4) * t233 - pkin(9) * t223;
t67 = -t347 * pkin(4) - pkin(9) * t281 - t233 * t206 + t76;
t30 = t262 * t67 - t267 * t62;
t31 = t262 * t62 + t267 * t67;
t16 = t262 * t31 - t267 * t30;
t345 = pkin(4) * t16;
t95 = t179 + t115;
t54 = -t262 * t92 - t267 * t95;
t344 = pkin(4) * t54;
t342 = pkin(3) * t264;
t340 = t16 * t263;
t339 = t16 * t268;
t162 = t198 * t264 - t269 * t348;
t130 = -t309 * pkin(3) - t346 * pkin(8) + t214 * t238 + t162;
t78 = pkin(4) * t281 - t347 * pkin(9) + t206 * t223 + t130;
t338 = t262 * t78;
t337 = t267 * t78;
t336 = -pkin(3) * t130 + pkin(8) * t45;
t116 = -t162 * t269 + t163 * t264;
t334 = t116 * t265;
t136 = t160 + t200;
t333 = t136 * t262;
t332 = t136 * t267;
t165 = t199 + t201;
t331 = t165 * t263;
t330 = t165 * t268;
t329 = t205 * t264;
t328 = t205 * t269;
t212 = t216 + t309;
t327 = t212 * t264;
t326 = t212 * t269;
t325 = t231 * t262;
t324 = t231 * t267;
t323 = t233 * t263;
t322 = t233 * t268;
t321 = t259 * t264;
t320 = t259 * t269;
t126 = t263 * t130;
t249 = t270 * t318;
t246 = qJDD(2) + t249;
t319 = t265 * t246;
t127 = t268 * t130;
t316 = t270 * (qJDD(2) - t249);
t314 = qJD(3) + t259;
t312 = qJD(4) + t233;
t308 = t264 * t160;
t307 = t269 * t160;
t306 = t264 * t199;
t305 = t269 * t199;
t190 = -t219 - t232;
t134 = -t190 * t263 - t330;
t158 = t312 * t221 + t279;
t304 = pkin(3) * t158 + pkin(8) * t134 + t126;
t181 = -t232 - t347;
t124 = t181 * t268 - t357;
t154 = -t312 * t223 - t293;
t303 = pkin(3) * t154 + pkin(8) * t124 - t127;
t302 = -pkin(3) * t269 - pkin(2);
t17 = t262 * t30 + t267 * t31;
t125 = -t191 - t192;
t56 = t262 * t95 - t267 * t92;
t10 = -pkin(4) * t125 + pkin(9) * t56 + t17;
t12 = -pkin(9) * t54 - t16;
t28 = -t263 * t54 + t268 * t56;
t297 = -pkin(3) * t125 + pkin(8) * t28 + t268 * t10 + t263 * t12;
t146 = -t230 - t191;
t100 = t146 * t267 - t358;
t91 = (qJD(5) + t231) * t195 + t294;
t34 = -pkin(4) * t91 + pkin(9) * t100 - t337;
t99 = t146 * t262 + t355;
t47 = -pkin(9) * t99 + t338;
t60 = t100 * t268 - t263 * t99;
t296 = -pkin(3) * t91 + pkin(8) * t60 + t263 * t47 + t268 * t34;
t167 = -t192 - t230;
t104 = -t167 * t262 - t332;
t36 = -pkin(4) * t351 + pkin(9) * t104 + t338;
t103 = t167 * t267 - t333;
t57 = -pkin(9) * t103 + t337;
t65 = -t103 * t263 + t104 * t268;
t295 = -pkin(3) * t351 + pkin(8) * t65 + t263 * t57 + t268 * t36;
t117 = t162 * t264 + t269 * t163;
t225 = t270 * g(3) + t274;
t292 = t265 * t225 + t270 * t226;
t109 = -t153 * t268 + t157 * t263;
t174 = t219 + t347;
t290 = pkin(3) * t174 + pkin(8) * t109 + t45;
t289 = t263 * t76 - t268 * t75;
t286 = pkin(4) * t99 - t30;
t14 = -pkin(4) * t78 + pkin(9) * t17;
t8 = t17 * t268 - t340;
t285 = -pkin(3) * t78 + pkin(8) * t8 - pkin(9) * t340 + t268 * t14;
t280 = pkin(4) * t103 - t31;
t278 = (-qJD(3) + t259) * t238 - t291;
t271 = qJD(2) ^ 2;
t260 = t265 ^ 2;
t256 = t260 * t272;
t245 = t255 - 0.2e1 * t300;
t242 = t254 + 0.2e1 * t299;
t239 = pkin(6) * t272 + t282;
t228 = -t235 + t346;
t227 = t234 - t346;
t224 = -t235 - t346;
t215 = t235 - t234;
t210 = -t346 - t234;
t208 = -t219 + t232;
t207 = -t232 + t347;
t204 = -t234 - t235;
t197 = t219 - t347;
t189 = -t224 * t264 - t326;
t188 = t224 * t269 - t327;
t187 = t203 + t229;
t185 = -t314 * t236 + t287;
t182 = t314 * t238 + t291;
t178 = t210 * t269 - t356;
t177 = t210 * t264 + t353;
t176 = -t192 + t230;
t175 = t191 - t230;
t169 = (-t221 * t268 + t223 * t263) * t233;
t168 = (-t221 * t263 - t223 * t268) * t233;
t159 = t192 - t191;
t156 = -t209 + t172;
t150 = t172 * t268 - t223 * t323;
t149 = t172 * t263 + t223 * t322;
t148 = t221 * t322 + t263 * t281;
t147 = t221 * t323 - t268 * t281;
t145 = t187 * t264 + t269 * t278;
t144 = -t187 * t269 + t264 * t278;
t143 = t207 * t268 - t331;
t142 = -t208 * t263 + t354;
t141 = t207 * t263 + t330;
t140 = t208 * t268 + t357;
t139 = (-t193 * t267 + t195 * t262) * t231;
t138 = (-t193 * t262 - t195 * t267) * t231;
t133 = t190 * t268 - t331;
t123 = t181 * t263 + t354;
t114 = -qJD(5) * t195 - t294;
t113 = t175 * t267 - t333;
t112 = -t176 * t262 + t355;
t111 = t175 * t262 + t332;
t110 = t176 * t267 + t358;
t108 = t154 * t268 - t156 * t263;
t107 = -t153 * t263 - t157 * t268;
t106 = t154 * t263 + t156 * t268;
t102 = t134 * t269 - t158 * t264;
t101 = t134 * t264 + t158 * t269;
t98 = t124 * t269 - t154 * t264;
t97 = t124 * t264 + t154 * t269;
t88 = t115 * t267 - t195 * t325;
t87 = t115 * t262 + t195 * t324;
t86 = -t114 * t262 + t193 * t324;
t85 = t114 * t267 + t193 * t325;
t84 = -t138 * t263 + t139 * t268;
t83 = t138 * t268 + t139 * t263;
t82 = -pkin(8) * t133 + t127;
t81 = t109 * t269 - t174 * t264;
t80 = t109 * t264 + t174 * t269;
t79 = -pkin(8) * t123 + t126;
t72 = -t111 * t263 + t113 * t268;
t71 = -t110 * t263 + t112 * t268;
t70 = t111 * t268 + t113 * t263;
t69 = t110 * t268 + t112 * t263;
t68 = -pkin(3) * t133 + t76;
t66 = -pkin(3) * t123 + t75;
t64 = t103 * t268 + t104 * t263;
t59 = t100 * t263 + t268 * t99;
t55 = -t262 * t351 - t267 * t91;
t53 = -t262 * t91 + t267 * t351;
t51 = -t263 * t87 + t268 * t88;
t50 = -t263 * t85 + t268 * t86;
t49 = t263 * t88 + t268 * t87;
t48 = t263 * t86 + t268 * t85;
t42 = t264 * t351 + t269 * t65;
t41 = t264 * t65 - t269 * t351;
t39 = -t130 * t269 + t264 * t45;
t38 = t264 * t91 + t269 * t60;
t37 = t264 * t60 - t269 * t91;
t32 = -pkin(8) * t107 - t289;
t27 = -t263 * t53 + t268 * t55;
t26 = t263 * t56 + t268 * t54;
t25 = t263 * t55 + t268 * t53;
t23 = t125 * t264 + t269 * t28;
t22 = -t125 * t269 + t264 * t28;
t21 = -pkin(3) * t26 - t344;
t20 = -pkin(3) * t64 - t280;
t19 = -pkin(3) * t59 - t286;
t18 = -pkin(8) * t64 - t263 * t36 + t268 * t57;
t15 = -pkin(8) * t59 - t263 * t34 + t268 * t47;
t7 = t17 * t263 + t339;
t5 = t264 * t78 + t269 * t8;
t4 = t264 * t8 - t269 * t78;
t3 = -pkin(3) * t7 - t345;
t2 = -pkin(8) * t26 - t10 * t263 + t12 * t268;
t1 = -pkin(8) * t7 - pkin(9) * t339 - t14 * t263;
t6 = [0, 0, 0, 0, 0, qJDD(1), t298, t284, 0, 0, (t243 + t299) * t265, t242 * t270 + t245 * t265, t319 + t270 * (-t256 + t271), (t244 - t300) * t270, t265 * (t257 - t271) + t316, 0, t270 * t239 + pkin(1) * t245 + pkin(6) * (t270 * (-t257 - t271) - t319), -t265 * t239 - pkin(1) * t242 + pkin(6) * (-t316 - t265 * (-t256 - t271)), pkin(1) * (t256 + t257) + (t260 + t261) * t335 + t292, pkin(1) * t239 + pkin(6) * t292, t265 * (t203 * t269 - t238 * t321) + t270 * (t203 * t264 + t238 * t320), t265 * (-t182 * t269 - t186 * t264) + t270 * (-t182 * t264 + t186 * t269), t265 * (-t228 * t264 + t353) + t270 * (t228 * t269 + t356), t265 * (-t202 * t264 + t236 * t320) + t270 * (t202 * t269 + t236 * t321), t265 * (t227 * t269 - t327) + t270 * (t227 * t264 + t326), (t265 * (-t236 * t269 + t238 * t264) + t270 * (-t236 * t264 - t238 * t269)) * t259, t265 * (-pkin(7) * t177 - t329) + t270 * (-pkin(2) * t182 + pkin(7) * t178 + t328) - pkin(1) * t182 + pkin(6) * (-t177 * t265 + t178 * t270), t265 * (-pkin(7) * t188 - t328) + t270 * (-pkin(2) * t185 + pkin(7) * t189 - t329) - pkin(1) * t185 + pkin(6) * (-t188 * t265 + t189 * t270), t265 * (-pkin(7) * t144 - t116) + t270 * (-pkin(2) * t204 + pkin(7) * t145 + t117) - pkin(1) * t204 + pkin(6) * (-t144 * t265 + t145 * t270), -pkin(7) * t334 + t270 * (pkin(2) * t205 + pkin(7) * t117) + pkin(1) * t205 + pkin(6) * (t117 * t270 - t334), t265 * (t150 * t269 + t306) + t270 * (t150 * t264 - t305), t265 * (t108 * t269 + t197 * t264) + t270 * (t108 * t264 - t197 * t269), t265 * (t142 * t269 + t157 * t264) + t270 * (t142 * t264 - t157 * t269), t265 * (t148 * t269 - t306) + t270 * (t148 * t264 + t305), t265 * (t143 * t269 - t153 * t264) + t270 * (t143 * t264 + t153 * t269), t265 * (t169 * t269 + t201 * t264) + t270 * (t169 * t264 - t201 * t269), t265 * (-pkin(7) * t97 - t264 * t66 + t269 * t79) + t270 * (-pkin(2) * t123 + pkin(7) * t98 + t264 * t79 + t269 * t66) - pkin(1) * t123 + pkin(6) * (-t265 * t97 + t270 * t98), t265 * (-pkin(7) * t101 - t264 * t68 + t269 * t82) + t270 * (-pkin(2) * t133 + pkin(7) * t102 + t264 * t82 + t269 * t68) - pkin(1) * t133 + pkin(6) * (-t101 * t265 + t102 * t270), t265 * (-pkin(7) * t80 + t269 * t32) + t270 * (pkin(7) * t81 + t264 * t32) + pkin(6) * (-t265 * t80 + t270 * t81) + (t265 * t342 + t270 * t302 - pkin(1)) * t107, (t265 * (-pkin(8) * t269 + t342) + t270 * (-pkin(8) * t264 + t302) - pkin(1)) * t289 + (pkin(6) + pkin(7)) * (-t265 * t39 + t270 * (t130 * t264 + t269 * t45)), t265 * (t269 * t51 + t308) + t270 * (t264 * t51 - t307), t265 * (t159 * t264 + t269 * t27) + t270 * (-t159 * t269 + t264 * t27), t265 * (t264 * t95 + t269 * t71) + t270 * (t264 * t71 - t269 * t95), t265 * (t269 * t50 - t308) + t270 * (t264 * t50 + t307), t265 * (-t264 * t92 + t269 * t72) + t270 * (t264 * t72 + t269 * t92), t265 * (t200 * t264 + t269 * t84) + t270 * (-t200 * t269 + t264 * t84), t265 * (-pkin(7) * t37 + t15 * t269 - t19 * t264) + t270 * (-pkin(2) * t59 + pkin(7) * t38 + t15 * t264 + t19 * t269) - pkin(1) * t59 + pkin(6) * (-t265 * t37 + t270 * t38), t265 * (-pkin(7) * t41 + t18 * t269 - t20 * t264) + t270 * (-pkin(2) * t64 + pkin(7) * t42 + t18 * t264 + t20 * t269) - pkin(1) * t64 + pkin(6) * (-t265 * t41 + t270 * t42), t265 * (-pkin(7) * t22 + t2 * t269 - t21 * t264) + t270 * (-pkin(2) * t26 + pkin(7) * t23 + t2 * t264 + t21 * t269) - pkin(1) * t26 + pkin(6) * (-t22 * t265 + t23 * t270), t265 * (-pkin(7) * t4 + t1 * t269 - t264 * t3) + t270 * (-pkin(2) * t7 + pkin(7) * t5 + t1 * t264 + t269 * t3) - pkin(1) * t7 + pkin(6) * (-t265 * t4 + t270 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t249, t256 - t257, t254, t249, t255, qJDD(2), -t225, -t226, 0, 0, t216, t215, t187, -t216, t278, t309, pkin(2) * t177 - t162, -t317 - t264 * (pkin(7) * t299 - t225 - t341) + (-t264 * t246 + t188) * pkin(2), pkin(2) * t144, pkin(2) * t116, t149, t106, t140, t147, t141, t168, pkin(2) * t97 + t303, pkin(2) * t101 + t304, pkin(2) * t80 + t290, pkin(2) * t39 + t336, t49, t25, t69, t48, t70, t83, pkin(2) * t37 + t296, pkin(2) * t41 + t295, pkin(2) * t22 + t297, pkin(2) * t4 + t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t216, t215, t187, -t216, t278, t309, -t162, -t163, 0, 0, t149, t106, t140, t147, t141, t168, t303, t304, t290, t336, t49, t25, t69, t48, t70, t83, t296, t295, t297, t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199, t197, t157, -t199, -t153, t201, -t75, -t76, 0, 0, t160, t159, t95, -t160, -t92, t200, t286, t280, t344, t345; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, t159, t95, -t160, -t92, t200, -t30, -t31, 0, 0;];
tauJ_reg = t6;
