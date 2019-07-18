% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRRPR2
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
% Datum: 2019-05-05 22:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRRPR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:06:35
% EndTime: 2019-05-05 22:06:55
% DurationCPUTime: 8.13s
% Computational Cost: add. (49761->498), mult. (101047->721), div. (0->0), fcn. (70599->12), ass. (0->307)
t282 = cos(qJ(3));
t257 = qJD(1) * t282 - qJD(4);
t357 = t257 ^ 2;
t277 = sin(qJ(4));
t281 = cos(qJ(4));
t278 = sin(qJ(3));
t322 = qJD(1) * t278;
t237 = -t281 * qJD(3) + t277 * t322;
t238 = qJD(3) * t277 + t281 * t322;
t270 = sin(pkin(11));
t272 = cos(pkin(11));
t216 = t272 * t237 + t238 * t270;
t359 = t216 ^ 2;
t174 = -t357 - t359;
t318 = qJD(1) * qJD(3);
t261 = t278 * t318;
t316 = t282 * qJDD(1);
t244 = -t261 + t316;
t236 = -qJDD(4) + t244;
t218 = -t237 * t270 + t238 * t272;
t336 = t218 * t216;
t290 = -t236 - t336;
t365 = t272 * t290;
t121 = t174 * t270 + t365;
t284 = qJD(1) ^ 2;
t279 = sin(qJ(1));
t283 = cos(qJ(1));
t310 = t279 * g(1) - g(2) * t283;
t239 = qJDD(1) * pkin(1) + t310;
t300 = g(1) * t283 + g(2) * t279;
t240 = -pkin(1) * t284 - t300;
t271 = sin(pkin(10));
t273 = cos(pkin(10));
t306 = t273 * t239 - t271 * t240;
t206 = -qJDD(1) * pkin(2) - t284 * pkin(7) - t306;
t297 = -t244 + t261;
t311 = t282 * t318;
t317 = t278 * qJDD(1);
t243 = t311 + t317;
t298 = t243 + t311;
t171 = pkin(3) * t297 - pkin(8) * t298 + t206;
t323 = t271 * t239 + t273 * t240;
t207 = -pkin(2) * t284 + qJDD(1) * pkin(7) + t323;
t301 = -t282 * pkin(3) - t278 * pkin(8);
t305 = t284 * t301 + t207;
t325 = -g(3) + qJDD(2);
t309 = t278 * t325;
t356 = qJD(3) ^ 2;
t184 = -t356 * pkin(3) + qJDD(3) * pkin(8) + t305 * t282 + t309;
t129 = -t281 * t171 + t184 * t277;
t294 = -qJDD(3) * t277 - t243 * t281;
t212 = -qJD(4) * t237 - t294;
t226 = t237 * t257;
t190 = t212 - t226;
t334 = t238 * t237;
t289 = -t236 - t334;
t104 = t289 * pkin(4) - qJ(5) * t190 - t129;
t130 = t277 * t171 + t281 * t184;
t223 = -pkin(4) * t257 - qJ(5) * t238;
t295 = qJDD(3) * t281 - t243 * t277;
t287 = -qJD(4) * t238 + t295;
t358 = t237 ^ 2;
t116 = -t358 * pkin(4) + qJ(5) * t287 + t257 * t223 + t130;
t299 = -0.2e1 * qJD(5) * t218 + t272 * t104 - t116 * t270;
t369 = pkin(4) * t121 + t299;
t276 = sin(qJ(6));
t280 = cos(qJ(6));
t180 = t280 * t216 + t218 * t276;
t182 = -t216 * t276 + t218 * t280;
t132 = t182 * t180;
t232 = -qJDD(6) + t236;
t361 = -t132 - t232;
t368 = t276 * t361;
t367 = t280 * t361;
t366 = t270 * t290;
t364 = t277 * t289;
t363 = t281 * t289;
t172 = -t212 * t270 + t272 * t287;
t173 = t272 * t212 + t270 * t287;
t115 = -qJD(6) * t180 + t172 * t276 + t173 * t280;
t251 = -qJD(6) + t257;
t160 = t180 * t251;
t362 = t115 + t160;
t200 = t216 * t257;
t150 = -t200 + t173;
t360 = t200 + t173;
t177 = t180 ^ 2;
t178 = t182 ^ 2;
t215 = t218 ^ 2;
t235 = t238 ^ 2;
t250 = t251 ^ 2;
t321 = qJD(5) * t216;
t210 = -0.2e1 * t321;
t324 = t270 * t104 + t272 * t116;
t63 = t210 + t324;
t39 = t270 * t63 + t272 * t299;
t355 = pkin(4) * t39;
t335 = t218 * t257;
t293 = t172 - t335;
t110 = -t150 * t272 + t270 * t293;
t354 = pkin(4) * t110;
t46 = t290 * pkin(5) - pkin(9) * t150 + t299;
t196 = -pkin(5) * t257 - pkin(9) * t218;
t47 = -pkin(5) * t359 + pkin(9) * t172 + t196 * t257 + t63;
t27 = t276 * t47 - t280 * t46;
t28 = t276 * t46 + t280 * t47;
t13 = -t27 * t280 + t276 * t28;
t352 = t13 * t270;
t351 = t13 * t272;
t260 = t282 * t325;
t183 = -qJDD(3) * pkin(3) - t356 * pkin(8) + t278 * t305 - t260;
t126 = -t287 * pkin(4) - t358 * qJ(5) + t223 * t238 + qJDD(5) + t183;
t77 = -t172 * pkin(5) - t359 * pkin(9) + t196 * t218 + t126;
t350 = t276 * t77;
t349 = t277 * t39;
t348 = t280 * t77;
t347 = t281 * t39;
t123 = -t132 + t232;
t346 = t123 * t276;
t345 = t123 * t280;
t344 = t126 * t270;
t343 = t126 * t272;
t167 = t236 - t336;
t342 = t167 * t270;
t341 = t167 * t272;
t340 = t183 * t277;
t339 = t183 * t281;
t202 = t236 - t334;
t338 = t202 * t277;
t337 = t202 * t281;
t333 = t251 * t276;
t332 = t251 * t280;
t331 = t257 * t270;
t330 = t257 * t272;
t329 = t257 * t277;
t328 = t257 * t281;
t256 = t282 * t284 * t278;
t249 = qJDD(3) + t256;
t327 = t278 * t249;
t248 = -t256 + qJDD(3);
t326 = t282 * t248;
t319 = qJD(4) - t257;
t315 = t282 * t132;
t314 = t282 * t336;
t313 = t282 * t334;
t312 = pkin(1) * t271 + pkin(7);
t14 = t27 * t276 + t280 * t28;
t40 = -t270 * t299 + t272 * t63;
t96 = t129 * t277 + t281 * t130;
t307 = -t280 * t172 + t173 * t276;
t194 = t207 * t278 - t260;
t195 = t282 * t207 + t309;
t153 = t278 * t194 + t282 * t195;
t8 = t14 * t270 + t351;
t304 = pkin(4) * t8 + pkin(5) * t13;
t288 = (-qJD(6) - t251) * t182 - t307;
t93 = t115 - t160;
t56 = t276 * t288 - t280 * t93;
t58 = t276 * t93 + t280 * t288;
t35 = t270 * t58 + t272 * t56;
t303 = pkin(4) * t35 + pkin(5) * t56;
t193 = -t215 - t357;
t134 = t193 * t272 + t342;
t302 = pkin(4) * t134 - t324;
t296 = t129 * t281 - t130 * t277;
t292 = -pkin(1) * t273 - pkin(2) + t301;
t128 = -t250 - t177;
t87 = t128 * t276 + t367;
t88 = t128 * t280 - t368;
t53 = t270 * t88 + t272 * t87;
t291 = pkin(4) * t53 + pkin(5) * t87 - t27;
t155 = -t178 - t250;
t100 = t155 * t280 + t346;
t101 = -t155 * t276 + t345;
t64 = t100 * t272 + t101 * t270;
t286 = pkin(4) * t64 + pkin(5) * t100 - t28;
t285 = (-qJD(4) - t257) * t238 + t295;
t267 = t282 ^ 2;
t266 = t278 ^ 2;
t264 = t267 * t284;
t262 = t266 * t284;
t254 = -t264 - t356;
t253 = -t262 - t356;
t247 = t262 + t264;
t246 = (t266 + t267) * qJDD(1);
t245 = -0.2e1 * t261 + t316;
t242 = 0.2e1 * t311 + t317;
t227 = t282 * t236;
t225 = -t235 + t357;
t224 = -t357 + t358;
t222 = -t253 * t278 - t326;
t221 = t254 * t282 - t327;
t220 = t235 - t358;
t219 = -t235 - t357;
t214 = -t357 - t358;
t201 = t235 + t358;
t198 = -t215 + t357;
t197 = -t357 + t359;
t191 = t319 * t237 + t294;
t189 = t212 + t226;
t187 = -t319 * t238 + t295;
t185 = t215 - t359;
t176 = -t219 * t277 + t337;
t175 = t219 * t281 + t338;
t166 = t214 * t281 - t364;
t165 = t214 * t277 + t363;
t159 = -t178 + t250;
t158 = t177 - t250;
t157 = (t216 * t272 - t218 * t270) * t257;
t156 = (t216 * t270 + t218 * t272) * t257;
t154 = -t215 - t359;
t152 = t190 * t277 + t281 * t285;
t145 = -t172 - t335;
t144 = t173 * t272 + t218 * t331;
t143 = t173 * t270 - t218 * t330;
t142 = -t172 * t270 - t216 * t330;
t141 = t172 * t272 - t216 * t331;
t140 = t197 * t272 + t342;
t139 = -t198 * t270 + t365;
t138 = t197 * t270 - t341;
t137 = t198 * t272 + t366;
t136 = t176 * t282 - t191 * t278;
t135 = -t193 * t270 + t341;
t133 = t166 * t282 - t187 * t278;
t131 = t178 - t177;
t122 = t174 * t272 - t366;
t119 = (t180 * t280 - t182 * t276) * t251;
t118 = (t180 * t276 + t182 * t280) * t251;
t117 = -t177 - t178;
t114 = -qJD(6) * t182 - t307;
t112 = t150 * t270 + t272 * t293;
t111 = -t145 * t272 - t270 * t360;
t109 = -t145 * t270 + t272 * t360;
t108 = t158 * t280 + t346;
t107 = -t159 * t276 + t367;
t106 = t158 * t276 - t345;
t105 = t159 * t280 + t368;
t99 = -t134 * t277 + t135 * t281;
t98 = t134 * t281 + t135 * t277;
t97 = -qJ(5) * t134 + t343;
t89 = (qJD(6) - t251) * t182 + t307;
t86 = t115 * t280 + t182 * t333;
t85 = t115 * t276 - t182 * t332;
t84 = -t114 * t276 - t180 * t332;
t83 = t114 * t280 - t180 * t333;
t82 = -qJ(5) * t121 + t344;
t81 = -t121 * t277 + t122 * t281;
t80 = t121 * t281 + t122 * t277;
t79 = -t118 * t270 + t119 * t272;
t78 = t118 * t272 + t119 * t270;
t75 = t278 * t360 + t282 * t99;
t74 = -pkin(4) * t360 + qJ(5) * t135 + t344;
t73 = -pkin(4) * t145 + qJ(5) * t122 - t343;
t72 = t145 * t278 + t282 * t81;
t71 = -t110 * t277 + t112 * t281;
t70 = t110 * t281 + t112 * t277;
t69 = -t106 * t270 + t108 * t272;
t68 = -t105 * t270 + t107 * t272;
t67 = t106 * t272 + t108 * t270;
t66 = t105 * t272 + t107 * t270;
t65 = -t100 * t270 + t101 * t272;
t60 = t154 * t278 + t282 * t71;
t59 = -pkin(9) * t100 + t348;
t57 = -t276 * t362 - t280 * t89;
t55 = -t276 * t89 + t280 * t362;
t54 = -t270 * t87 + t272 * t88;
t52 = -t270 * t85 + t272 * t86;
t51 = -t270 * t83 + t272 * t84;
t50 = t270 * t86 + t272 * t85;
t49 = t270 * t84 + t272 * t83;
t48 = -pkin(9) * t87 + t350;
t44 = -pkin(5) * t362 + pkin(9) * t101 + t350;
t43 = -pkin(5) * t89 + pkin(9) * t88 - t348;
t42 = -t277 * t64 + t281 * t65;
t41 = t277 * t65 + t281 * t64;
t38 = -pkin(4) * t126 + qJ(5) * t40;
t37 = -t270 * t56 + t272 * t58;
t36 = -t270 * t55 + t272 * t57;
t34 = t270 * t57 + t272 * t55;
t33 = -qJ(5) * t110 - t39;
t32 = -t277 * t53 + t281 * t54;
t31 = t277 * t54 + t281 * t53;
t30 = t278 * t362 + t282 * t42;
t29 = -pkin(4) * t154 + qJ(5) * t112 + t40;
t25 = t278 * t89 + t282 * t32;
t24 = -qJ(5) * t64 - t270 * t44 + t272 * t59;
t23 = t281 * t40 - t349;
t22 = t277 * t40 + t347;
t21 = t126 * t278 + t23 * t282;
t20 = -qJ(5) * t53 - t270 * t43 + t272 * t48;
t19 = -pkin(4) * t362 + qJ(5) * t65 + t270 * t59 + t272 * t44;
t18 = -t277 * t35 + t281 * t37;
t17 = t277 * t37 + t281 * t35;
t16 = -pkin(4) * t89 + qJ(5) * t54 + t270 * t48 + t272 * t43;
t15 = t117 * t278 + t18 * t282;
t12 = -pkin(5) * t77 + pkin(9) * t14;
t11 = -pkin(9) * t56 - t13;
t10 = -pkin(5) * t117 + pkin(9) * t58 + t14;
t9 = t14 * t272 - t352;
t7 = -qJ(5) * t35 - t10 * t270 + t11 * t272;
t6 = -pkin(4) * t117 + qJ(5) * t37 + t10 * t272 + t11 * t270;
t5 = -t277 * t8 + t281 * t9;
t4 = t277 * t9 + t281 * t8;
t3 = -pkin(9) * t351 - qJ(5) * t8 - t12 * t270;
t2 = t278 * t77 + t282 * t5;
t1 = -pkin(4) * t77 - pkin(9) * t352 + qJ(5) * t9 + t12 * t272;
t26 = [0, 0, 0, 0, 0, qJDD(1), t310, t300, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t273 - t271 * t284) + t306, pkin(1) * (-qJDD(1) * t271 - t273 * t284) - t323, 0, pkin(1) * (t271 * t323 + t273 * t306), t298 * t278, t242 * t282 + t245 * t278, t327 + t282 * (-t262 + t356), -t297 * t282, t278 * (t264 - t356) + t326, 0, -t282 * t206 + pkin(2) * t245 + pkin(7) * t221 + pkin(1) * (t221 * t271 + t245 * t273), t278 * t206 - pkin(2) * t242 + pkin(7) * t222 + pkin(1) * (t222 * t271 - t242 * t273), pkin(2) * t247 + pkin(7) * t246 + pkin(1) * (t246 * t271 + t247 * t273) + t153, -pkin(2) * t206 + pkin(7) * t153 + pkin(1) * (t153 * t271 - t206 * t273), t278 * (t212 * t281 + t238 * t329) - t313, t278 * (t187 * t281 - t189 * t277) - t282 * t220, t278 * (-t225 * t277 + t363) - t282 * t190, t278 * (-t237 * t328 - t277 * t287) + t313, t278 * (t224 * t281 + t338) - t282 * t285, t227 + t278 * (t237 * t281 - t238 * t277) * t257, t278 * (-pkin(8) * t165 + t340) + t282 * (-pkin(3) * t165 + t129) - pkin(2) * t165 + pkin(7) * t133 + pkin(1) * (t133 * t271 - t165 * t273), t278 * (-pkin(8) * t175 + t339) + t282 * (-pkin(3) * t175 + t130) - pkin(2) * t175 + pkin(7) * t136 + pkin(1) * (t136 * t271 - t175 * t273), t278 * t296 + t312 * (t152 * t282 - t201 * t278) + t292 * (-t190 * t281 + t277 * t285), t312 * (t183 * t278 + t282 * t96) - t292 * t296, t278 * (-t143 * t277 + t144 * t281) - t314, t278 * (-t109 * t277 + t111 * t281) - t282 * t185, t278 * (-t137 * t277 + t139 * t281) - t282 * t150, t278 * (-t141 * t277 + t142 * t281) + t314, t278 * (-t138 * t277 + t140 * t281) - t282 * t293, t278 * (-t156 * t277 + t157 * t281) + t227, t278 * (-pkin(8) * t80 - t277 * t73 + t281 * t82) + t282 * (-pkin(3) * t80 - t369) - pkin(2) * t80 + pkin(7) * t72 + pkin(1) * (t271 * t72 - t273 * t80), t278 * (-pkin(8) * t98 - t277 * t74 + t281 * t97) + t282 * (-pkin(3) * t98 + t210 - t302) - pkin(2) * t98 + pkin(7) * t75 + pkin(1) * (t271 * t75 - t273 * t98), t278 * (-pkin(8) * t70 - t277 * t29 + t281 * t33) + t282 * (-pkin(3) * t70 - t354) - pkin(2) * t70 + pkin(7) * t60 + pkin(1) * (t271 * t60 - t273 * t70), t278 * (-pkin(8) * t22 - qJ(5) * t347 - t277 * t38) + t282 * (-pkin(3) * t22 - t355) - pkin(2) * t22 + pkin(7) * t21 + pkin(1) * (t21 * t271 - t22 * t273), t278 * (-t277 * t50 + t281 * t52) - t315, t278 * (-t277 * t34 + t281 * t36) - t282 * t131, t278 * (-t277 * t66 + t281 * t68) - t282 * t93, t278 * (-t277 * t49 + t281 * t51) + t315, t278 * (-t277 * t67 + t281 * t69) - t282 * t288, t278 * (-t277 * t78 + t281 * t79) + t282 * t232, t278 * (-pkin(8) * t31 - t16 * t277 + t20 * t281) + t282 * (-pkin(3) * t31 - t291) - pkin(2) * t31 + pkin(7) * t25 + pkin(1) * (t25 * t271 - t273 * t31), t278 * (-pkin(8) * t41 - t19 * t277 + t24 * t281) + t282 * (-pkin(3) * t41 - t286) - pkin(2) * t41 + pkin(7) * t30 + pkin(1) * (t271 * t30 - t273 * t41), t278 * (-pkin(8) * t17 - t277 * t6 + t281 * t7) + t282 * (-pkin(3) * t17 - t303) - pkin(2) * t17 + pkin(7) * t15 + pkin(1) * (t15 * t271 - t17 * t273), t278 * (-pkin(8) * t4 - t1 * t277 + t281 * t3) + t282 * (-pkin(3) * t4 - t304) - pkin(2) * t4 + pkin(7) * t2 + pkin(1) * (t2 * t271 - t273 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t325, 0, 0, 0, 0, 0, 0, t249 * t282 + t254 * t278, -t248 * t278 + t253 * t282, 0, -t194 * t282 + t195 * t278, 0, 0, 0, 0, 0, 0, t166 * t278 + t187 * t282, t176 * t278 + t191 * t282, t152 * t278 + t201 * t282, -t183 * t282 + t278 * t96, 0, 0, 0, 0, 0, 0, -t145 * t282 + t278 * t81, t278 * t99 - t282 * t360, -t154 * t282 + t278 * t71, -t126 * t282 + t23 * t278, 0, 0, 0, 0, 0, 0, t278 * t32 - t282 * t89, t278 * t42 - t282 * t362, -t117 * t282 + t18 * t278, t278 * t5 - t282 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t256, t262 - t264, t317, t256, t316, qJDD(3), -t194, -t195, 0, 0, t212 * t277 - t238 * t328, t187 * t277 + t189 * t281, t225 * t281 + t364, -t237 * t329 + t281 * t287, t224 * t277 - t337, (t237 * t277 + t238 * t281) * t257, pkin(3) * t187 + pkin(8) * t166 - t339, pkin(3) * t191 + pkin(8) * t176 + t340, pkin(3) * t201 + pkin(8) * t152 + t96, -pkin(3) * t183 + pkin(8) * t96, t143 * t281 + t144 * t277, t109 * t281 + t111 * t277, t137 * t281 + t139 * t277, t141 * t281 + t142 * t277, t138 * t281 + t140 * t277, t156 * t281 + t157 * t277, -pkin(3) * t145 + pkin(8) * t81 + t277 * t82 + t281 * t73, -pkin(3) * t360 + pkin(8) * t99 + t277 * t97 + t281 * t74, -pkin(3) * t154 + pkin(8) * t71 + t277 * t33 + t281 * t29, -pkin(3) * t126 + pkin(8) * t23 - qJ(5) * t349 + t281 * t38, t277 * t52 + t281 * t50, t277 * t36 + t281 * t34, t277 * t68 + t281 * t66, t277 * t51 + t281 * t49, t277 * t69 + t281 * t67, t277 * t79 + t281 * t78, -pkin(3) * t89 + pkin(8) * t32 + t16 * t281 + t20 * t277, -pkin(3) * t362 + pkin(8) * t42 + t19 * t281 + t24 * t277, -pkin(3) * t117 + pkin(8) * t18 + t277 * t7 + t281 * t6, -pkin(3) * t77 + pkin(8) * t5 + t1 * t281 + t277 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t334, t220, t190, -t334, t285, -t236, -t129, -t130, 0, 0, t336, t185, t150, -t336, t293, -t236, t369, t302 + 0.2e1 * t321, t354, t355, t132, t131, t93, -t132, t288, -t232, t291, t286, t303, t304; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145, t360, t154, t126, 0, 0, 0, 0, 0, 0, t89, t362, t117, t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t131, t93, -t132, t288, -t232, -t27, -t28, 0, 0;];
tauJ_reg  = t26;