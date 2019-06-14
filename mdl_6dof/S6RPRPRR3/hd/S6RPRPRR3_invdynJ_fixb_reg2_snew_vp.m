% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 18:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPRR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:34:58
% EndTime: 2019-05-05 18:35:14
% DurationCPUTime: 7.26s
% Computational Cost: add. (45634->493), mult. (96172->718), div. (0->0), fcn. (67340->12), ass. (0->294)
t274 = sin(qJ(6));
t269 = sin(pkin(11));
t271 = cos(pkin(11));
t276 = sin(qJ(3));
t310 = qJD(1) * t276;
t233 = -t271 * qJD(3) + t269 * t310;
t234 = t269 * qJD(3) + t271 * t310;
t275 = sin(qJ(5));
t279 = cos(qJ(5));
t208 = t279 * t233 + t234 * t275;
t210 = -t233 * t275 + t234 * t279;
t278 = cos(qJ(6));
t178 = t278 * t208 + t210 * t274;
t180 = -t208 * t274 + t210 * t278;
t130 = t180 * t178;
t308 = qJD(1) * qJD(3);
t258 = t276 * t308;
t280 = cos(qJ(3));
t306 = t280 * qJDD(1);
t241 = -t258 + t306;
t235 = -qJDD(5) + t241;
t231 = -qJDD(6) + t235;
t343 = -t130 - t231;
t350 = t274 * t343;
t349 = t278 * t343;
t282 = qJD(1) ^ 2;
t319 = t233 * t234;
t288 = -t241 - t319;
t348 = t269 * t288;
t347 = t271 * t288;
t320 = t210 * t208;
t286 = -t235 - t320;
t346 = t275 * t286;
t345 = t279 * t286;
t301 = t280 * t308;
t307 = t276 * qJDD(1);
t240 = t301 + t307;
t216 = t269 * qJDD(3) + t271 * t240;
t295 = -t271 * qJDD(3) + t240 * t269;
t297 = t275 * t216 + t279 * t295;
t160 = -qJD(5) * t210 - t297;
t161 = -t208 * qJD(5) + t279 * t216 - t275 * t295;
t114 = -qJD(6) * t178 + t160 * t274 + t161 * t278;
t309 = t280 * qJD(1);
t254 = -qJD(5) + t309;
t248 = -qJD(6) + t254;
t164 = t178 * t248;
t344 = t114 + t164;
t195 = t208 * t254;
t148 = -t195 + t161;
t342 = t195 + t161;
t221 = t233 * t309;
t196 = -t216 + t221;
t222 = t234 * t309;
t198 = -t295 - t222;
t175 = t178 ^ 2;
t176 = t180 ^ 2;
t341 = t208 ^ 2;
t207 = t210 ^ 2;
t340 = t233 ^ 2;
t232 = t234 ^ 2;
t247 = t248 ^ 2;
t252 = t254 ^ 2;
t339 = qJD(3) ^ 2;
t277 = sin(qJ(1));
t281 = cos(qJ(1));
t300 = t277 * g(1) - t281 * g(2);
t236 = qJDD(1) * pkin(1) + t300;
t293 = t281 * g(1) + t277 * g(2);
t238 = -t282 * pkin(1) - t293;
t270 = sin(pkin(10));
t272 = cos(pkin(10));
t296 = t272 * t236 - t270 * t238;
t203 = -qJDD(1) * pkin(2) - t282 * pkin(7) - t296;
t291 = t240 + t301;
t173 = -t291 * qJ(4) + (-t241 + t258) * pkin(3) + t203;
t311 = t270 * t236 + t272 * t238;
t204 = -pkin(2) * t282 + qJDD(1) * pkin(7) + t311;
t292 = -t280 * pkin(3) - t276 * qJ(4);
t294 = t282 * t292 + t204;
t312 = -g(3) + qJDD(2);
t299 = t276 * t312;
t182 = -pkin(3) * t339 + qJDD(3) * qJ(4) + t280 * t294 + t299;
t119 = 0.2e1 * qJD(4) * t234 - t271 * t173 + t269 * t182;
t112 = t288 * pkin(4) + pkin(8) * t196 - t119;
t120 = -0.2e1 * qJD(4) * t233 + t269 * t173 + t271 * t182;
t217 = -pkin(4) * t309 - t234 * pkin(8);
t116 = -pkin(4) * t340 - pkin(8) * t295 + t217 * t309 + t120;
t70 = -t279 * t112 + t275 * t116;
t46 = t286 * pkin(5) - pkin(9) * t148 - t70;
t190 = -pkin(5) * t254 - pkin(9) * t210;
t71 = t275 * t112 + t279 * t116;
t56 = -pkin(5) * t341 + pkin(9) * t160 + t190 * t254 + t71;
t28 = t274 * t56 - t278 * t46;
t29 = t274 * t46 + t278 * t56;
t14 = t274 * t29 - t278 * t28;
t338 = pkin(5) * t14;
t298 = -t278 * t160 + t274 * t161;
t285 = (-qJD(6) - t248) * t180 - t298;
t90 = t114 - t164;
t53 = t274 * t285 - t278 * t90;
t337 = pkin(5) * t53;
t336 = t14 * t275;
t335 = t14 * t279;
t42 = t275 * t71 - t279 * t70;
t334 = t269 * t42;
t333 = t271 * t42;
t257 = t280 * t312;
t181 = -qJDD(3) * pkin(3) - t339 * qJ(4) + t276 * t294 + qJDD(4) - t257;
t131 = t295 * pkin(4) - t340 * pkin(8) + t234 * t217 + t181;
t81 = -t160 * pkin(5) - pkin(9) * t341 + t210 * t190 + t131;
t332 = t274 * t81;
t331 = t278 * t81;
t123 = -t130 + t231;
t330 = t123 * t274;
t329 = t123 * t278;
t328 = t131 * t275;
t327 = t131 * t279;
t168 = t235 - t320;
t326 = t168 * t275;
t325 = t168 * t279;
t324 = t181 * t269;
t323 = t181 * t271;
t200 = t241 - t319;
t322 = t200 * t269;
t321 = t200 * t271;
t318 = t248 * t274;
t317 = t248 * t278;
t316 = t254 * t275;
t315 = t254 * t279;
t253 = t280 * t282 * t276;
t246 = qJDD(3) + t253;
t314 = t276 * t246;
t245 = -t253 + qJDD(3);
t313 = t280 * t245;
t305 = t280 * t130;
t304 = t280 * t320;
t303 = t280 * t319;
t302 = pkin(1) * t270 + pkin(7);
t15 = t274 * t28 + t278 * t29;
t43 = t275 * t70 + t279 * t71;
t78 = t119 * t269 + t271 * t120;
t188 = t204 * t276 - t257;
t189 = t280 * t204 + t299;
t151 = t276 * t188 + t280 * t189;
t290 = t119 * t271 - t120 * t269;
t128 = -t247 - t175;
t94 = t128 * t274 + t349;
t289 = pkin(5) * t94 - t28;
t153 = -t176 - t247;
t100 = t153 * t278 + t330;
t287 = pkin(5) * t100 - t29;
t284 = (-qJD(5) - t254) * t210 - t297;
t283 = -pkin(1) * t272 - pkin(2) + t292;
t266 = t280 ^ 2;
t265 = t276 ^ 2;
t262 = t266 * t282;
t261 = t265 * t282;
t251 = -t262 - t339;
t250 = -t261 - t339;
t244 = t261 + t262;
t243 = (t265 + t266) * qJDD(1);
t242 = -0.2e1 * t258 + t306;
t239 = 0.2e1 * t301 + t307;
t229 = t280 * t241;
t220 = -t232 - t262;
t219 = -t232 + t262;
t218 = -t262 + t340;
t214 = -t250 * t276 - t313;
t213 = t251 * t280 - t314;
t211 = -t262 - t340;
t199 = t216 + t221;
t197 = -t222 + t295;
t193 = -t232 - t340;
t192 = -t207 + t252;
t191 = -t252 + t341;
t187 = -t207 - t252;
t185 = -t220 * t269 + t321;
t184 = t220 * t271 + t322;
t183 = t207 - t341;
t174 = -t252 - t341;
t172 = t211 * t271 - t348;
t171 = t211 * t269 + t347;
t163 = -t176 + t247;
t162 = t175 - t247;
t159 = -t196 * t269 + t198 * t271;
t156 = (t208 * t279 - t210 * t275) * t254;
t155 = (t208 * t275 + t210 * t279) * t254;
t152 = -t207 - t341;
t150 = t185 * t280 + t199 * t276;
t149 = t172 * t280 + t197 * t276;
t143 = (qJD(5) - t254) * t210 + t297;
t142 = t191 * t279 + t326;
t141 = -t192 * t275 + t345;
t140 = t191 * t275 - t325;
t139 = t192 * t279 + t346;
t138 = t161 * t279 + t210 * t316;
t137 = t161 * t275 - t210 * t315;
t136 = -t160 * t275 - t208 * t315;
t135 = t160 * t279 - t208 * t316;
t133 = -t187 * t275 + t325;
t132 = t187 * t279 + t326;
t129 = t176 - t175;
t127 = t174 * t279 - t346;
t126 = t174 * t275 + t345;
t122 = (t178 * t278 - t180 * t274) * t248;
t121 = (t178 * t274 + t180 * t278) * t248;
t117 = -t175 - t176;
t113 = -qJD(6) * t180 - t298;
t111 = t162 * t278 + t330;
t110 = -t163 * t274 + t349;
t109 = t162 * t274 - t329;
t108 = t163 * t278 + t350;
t105 = t148 * t275 + t279 * t284;
t104 = -t143 * t279 - t275 * t342;
t103 = -t148 * t279 + t275 * t284;
t102 = -t143 * t275 + t279 * t342;
t101 = -t153 * t274 + t329;
t99 = -t132 * t269 + t133 * t271;
t98 = t132 * t271 + t133 * t269;
t97 = -pkin(8) * t132 + t327;
t96 = -pkin(8) * t126 + t328;
t95 = t128 * t278 - t350;
t93 = -t126 * t269 + t127 * t271;
t92 = t126 * t271 + t127 * t269;
t86 = (qJD(6) - t248) * t180 + t298;
t85 = t114 * t278 + t180 * t318;
t84 = t114 * t274 - t180 * t317;
t83 = -t113 * t274 - t178 * t317;
t82 = t113 * t278 - t178 * t318;
t80 = -t121 * t275 + t122 * t279;
t79 = t121 * t279 + t122 * t275;
t76 = -pkin(4) * t342 + pkin(8) * t133 + t328;
t75 = t276 * t342 + t280 * t99;
t73 = -pkin(4) * t143 + pkin(8) * t127 - t327;
t72 = t143 * t276 + t280 * t93;
t68 = -t109 * t275 + t111 * t279;
t67 = -t108 * t275 + t110 * t279;
t66 = t109 * t279 + t111 * t275;
t65 = t108 * t279 + t110 * t275;
t64 = -t103 * t269 + t105 * t271;
t63 = t103 * t271 + t105 * t269;
t62 = -t100 * t275 + t101 * t279;
t61 = t100 * t279 + t101 * t275;
t60 = -pkin(9) * t100 + t331;
t59 = t152 * t276 + t280 * t64;
t58 = -t275 * t94 + t279 * t95;
t57 = t275 * t95 + t279 * t94;
t55 = t274 * t90 + t278 * t285;
t54 = -t274 * t344 - t278 * t86;
t52 = -t274 * t86 + t278 * t344;
t51 = -pkin(9) * t94 + t332;
t50 = -t275 * t84 + t279 * t85;
t49 = -t275 * t82 + t279 * t83;
t48 = t275 * t85 + t279 * t84;
t47 = t275 * t83 + t279 * t82;
t44 = -pkin(5) * t344 + pkin(9) * t101 + t332;
t41 = -pkin(5) * t86 + pkin(9) * t95 - t331;
t40 = -t269 * t61 + t271 * t62;
t39 = t269 * t62 + t271 * t61;
t38 = -pkin(4) * t131 + pkin(8) * t43;
t37 = -pkin(8) * t103 - t42;
t36 = -t269 * t57 + t271 * t58;
t35 = t269 * t58 + t271 * t57;
t34 = -pkin(4) * t152 + pkin(8) * t105 + t43;
t33 = -t275 * t53 + t279 * t55;
t32 = -t275 * t52 + t279 * t54;
t31 = t275 * t55 + t279 * t53;
t30 = t275 * t54 + t279 * t52;
t26 = t276 * t344 + t280 * t40;
t25 = t276 * t86 + t280 * t36;
t24 = t271 * t43 - t334;
t23 = t269 * t43 + t333;
t22 = -pkin(8) * t61 - t275 * t44 + t279 * t60;
t21 = t131 * t276 + t24 * t280;
t20 = -pkin(8) * t57 - t275 * t41 + t279 * t51;
t19 = -pkin(4) * t344 + pkin(8) * t62 + t275 * t60 + t279 * t44;
t18 = -pkin(4) * t86 + pkin(8) * t58 + t275 * t51 + t279 * t41;
t17 = -t269 * t31 + t271 * t33;
t16 = t269 * t33 + t271 * t31;
t13 = t117 * t276 + t17 * t280;
t12 = -pkin(5) * t81 + pkin(9) * t15;
t11 = -pkin(9) * t53 - t14;
t10 = -pkin(5) * t117 + pkin(9) * t55 + t15;
t9 = t15 * t279 - t336;
t8 = t15 * t275 + t335;
t7 = -pkin(8) * t31 - t10 * t275 + t11 * t279;
t6 = -pkin(4) * t117 + pkin(8) * t33 + t10 * t279 + t11 * t275;
t5 = -t269 * t8 + t271 * t9;
t4 = t269 * t9 + t271 * t8;
t3 = -pkin(8) * t8 - pkin(9) * t335 - t12 * t275;
t2 = t276 * t81 + t280 * t5;
t1 = -pkin(4) * t81 + pkin(8) * t9 - pkin(9) * t336 + t12 * t279;
t27 = [0, 0, 0, 0, 0, qJDD(1), t300, t293, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t272 - t270 * t282) + t296, pkin(1) * (-qJDD(1) * t270 - t272 * t282) - t311, 0, pkin(1) * (t270 * t311 + t272 * t296), t291 * t276, t239 * t280 + t242 * t276, t314 + t280 * (-t261 + t339), -t276 * t301 + t229, t276 * (t262 - t339) + t313, 0, -t280 * t203 + pkin(2) * t242 + pkin(7) * t213 + pkin(1) * (t213 * t270 + t242 * t272), t276 * t203 - pkin(2) * t239 + pkin(7) * t214 + pkin(1) * (t214 * t270 - t239 * t272), pkin(2) * t244 + pkin(7) * t243 + pkin(1) * (t243 * t270 + t244 * t272) + t151, -pkin(2) * t203 + pkin(7) * t151 + pkin(1) * (t151 * t270 - t203 * t272), t276 * (t216 * t271 + t222 * t269) - t303, t276 * (-t197 * t271 - t199 * t269) + t280 * (-t232 + t340), t276 * (-t219 * t269 + t347) + t280 * t196, t276 * (-t221 * t271 + t269 * t295) + t303, t276 * (t218 * t271 + t322) - t280 * t198, t229 + t276 * (t233 * t271 - t234 * t269) * t309, t276 * (-qJ(4) * t171 + t324) + t280 * (-pkin(3) * t171 + t119) - pkin(2) * t171 + pkin(7) * t149 + pkin(1) * (t149 * t270 - t171 * t272), t276 * (-qJ(4) * t184 + t323) + t280 * (-pkin(3) * t184 + t120) - pkin(2) * t184 + pkin(7) * t150 + pkin(1) * (t150 * t270 - t184 * t272), t276 * t290 + t302 * (t159 * t280 + t193 * t276) + t283 * (t196 * t271 + t198 * t269), t302 * (t181 * t276 + t280 * t78) - t283 * t290, t276 * (-t137 * t269 + t138 * t271) - t304, t276 * (-t102 * t269 + t104 * t271) - t280 * t183, t276 * (-t139 * t269 + t141 * t271) - t280 * t148, t276 * (-t135 * t269 + t136 * t271) + t304, t276 * (-t140 * t269 + t142 * t271) - t280 * t284, t276 * (-t155 * t269 + t156 * t271) + t280 * t235, t276 * (-qJ(4) * t92 - t269 * t73 + t271 * t96) + t280 * (-pkin(3) * t92 - pkin(4) * t126 + t70) - pkin(2) * t92 + pkin(7) * t72 + pkin(1) * (t270 * t72 - t272 * t92), t276 * (-qJ(4) * t98 - t269 * t76 + t271 * t97) + t280 * (-pkin(3) * t98 - pkin(4) * t132 + t71) - pkin(2) * t98 + pkin(7) * t75 + pkin(1) * (t270 * t75 - t272 * t98), t276 * (-qJ(4) * t63 - t269 * t34 + t271 * t37) + t280 * (-pkin(3) * t63 - pkin(4) * t103) - pkin(2) * t63 + pkin(7) * t59 + pkin(1) * (t270 * t59 - t272 * t63), t276 * (-pkin(8) * t333 - qJ(4) * t23 - t269 * t38) + t280 * (-pkin(3) * t23 - pkin(4) * t42) - pkin(2) * t23 + pkin(7) * t21 + pkin(1) * (t21 * t270 - t23 * t272), t276 * (-t269 * t48 + t271 * t50) - t305, t276 * (-t269 * t30 + t271 * t32) - t280 * t129, t276 * (-t269 * t65 + t271 * t67) - t280 * t90, t276 * (-t269 * t47 + t271 * t49) + t305, t276 * (-t269 * t66 + t271 * t68) - t280 * t285, t276 * (-t269 * t79 + t271 * t80) + t280 * t231, t276 * (-qJ(4) * t35 - t18 * t269 + t20 * t271) + t280 * (-pkin(3) * t35 - pkin(4) * t57 - t289) - pkin(2) * t35 + pkin(7) * t25 + pkin(1) * (t25 * t270 - t272 * t35), t276 * (-qJ(4) * t39 - t19 * t269 + t22 * t271) + t280 * (-pkin(3) * t39 - pkin(4) * t61 - t287) - pkin(2) * t39 + pkin(7) * t26 + pkin(1) * (t26 * t270 - t272 * t39), t276 * (-qJ(4) * t16 - t269 * t6 + t271 * t7) + t280 * (-pkin(3) * t16 - pkin(4) * t31 - t337) - pkin(2) * t16 + pkin(7) * t13 + pkin(1) * (t13 * t270 - t16 * t272), t276 * (-qJ(4) * t4 - t1 * t269 + t271 * t3) + t280 * (-pkin(3) * t4 - pkin(4) * t8 - t338) - pkin(2) * t4 + pkin(7) * t2 + pkin(1) * (t2 * t270 - t272 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t312, 0, 0, 0, 0, 0, 0, t246 * t280 + t251 * t276, -t245 * t276 + t250 * t280, 0, -t188 * t280 + t189 * t276, 0, 0, 0, 0, 0, 0, t172 * t276 - t197 * t280, t185 * t276 - t199 * t280, t159 * t276 - t193 * t280, -t181 * t280 + t276 * t78, 0, 0, 0, 0, 0, 0, -t143 * t280 + t276 * t93, t276 * t99 - t280 * t342, -t152 * t280 + t276 * t64, -t131 * t280 + t24 * t276, 0, 0, 0, 0, 0, 0, t276 * t36 - t280 * t86, t276 * t40 - t280 * t344, -t117 * t280 + t17 * t276, t276 * t5 - t280 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t253, t261 - t262, t307, t253, t306, qJDD(3), -t188, -t189, 0, 0, t216 * t269 - t222 * t271, -t197 * t269 + t199 * t271, t219 * t271 + t348, -t221 * t269 - t271 * t295, t218 * t269 - t321, (t233 * t269 + t234 * t271) * t309, -pkin(3) * t197 + qJ(4) * t172 - t323, -pkin(3) * t199 + qJ(4) * t185 + t324, -pkin(3) * t193 + qJ(4) * t159 + t78, -pkin(3) * t181 + qJ(4) * t78, t137 * t271 + t138 * t269, t102 * t271 + t104 * t269, t139 * t271 + t141 * t269, t135 * t271 + t136 * t269, t140 * t271 + t142 * t269, t155 * t271 + t156 * t269, -pkin(3) * t143 + qJ(4) * t93 + t269 * t96 + t271 * t73, -pkin(3) * t342 + qJ(4) * t99 + t269 * t97 + t271 * t76, -pkin(3) * t152 + qJ(4) * t64 + t269 * t37 + t271 * t34, -pkin(3) * t131 - pkin(8) * t334 + qJ(4) * t24 + t271 * t38, t269 * t50 + t271 * t48, t269 * t32 + t271 * t30, t269 * t67 + t271 * t65, t269 * t49 + t271 * t47, t269 * t68 + t271 * t66, t269 * t80 + t271 * t79, -pkin(3) * t86 + qJ(4) * t36 + t18 * t271 + t20 * t269, -pkin(3) * t344 + qJ(4) * t40 + t19 * t271 + t22 * t269, -pkin(3) * t117 + qJ(4) * t17 + t269 * t7 + t271 * t6, -pkin(3) * t81 + qJ(4) * t5 + t1 * t271 + t269 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, t199, t193, t181, 0, 0, 0, 0, 0, 0, t143, t342, t152, t131, 0, 0, 0, 0, 0, 0, t86, t344, t117, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t320, t183, t148, -t320, t284, -t235, -t70, -t71, 0, 0, t130, t129, t90, -t130, t285, -t231, t289, t287, t337, t338; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t129, t90, -t130, t285, -t231, -t28, -t29, 0, 0;];
tauJ_reg  = t27;
