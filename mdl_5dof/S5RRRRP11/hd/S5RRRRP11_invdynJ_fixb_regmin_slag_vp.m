% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRP11
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP11_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP11_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP11_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP11_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP11_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP11_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:18:35
% EndTime: 2019-12-31 22:18:55
% DurationCPUTime: 8.30s
% Computational Cost: add. (7379->567), mult. (18381->766), div. (0->0), fcn. (14404->10), ass. (0->246)
t202 = cos(qJ(2));
t196 = sin(pkin(5));
t311 = qJD(1) * t196;
t186 = t202 * t311;
t247 = t186 - qJD(3);
t198 = sin(qJ(3));
t199 = sin(qJ(2));
t201 = cos(qJ(3));
t329 = cos(pkin(5));
t277 = t329 * qJD(1);
t241 = t277 + qJD(2);
t220 = qJD(3) * t241;
t270 = t329 * qJDD(1);
t236 = t270 + qJDD(2);
t309 = qJD(2) * t202;
t286 = t198 * t309;
t300 = qJDD(1) * t199;
t305 = qJD(3) * t201;
t55 = (qJD(1) * (t199 * t305 + t286) + t198 * t300) * t196 + t198 * t220 - t201 * t236;
t50 = qJDD(4) + t55;
t264 = pkin(1) * t277;
t145 = pkin(7) * t186 + t199 * t264;
t369 = -t145 - t247 * (pkin(3) * t198 - pkin(9) * t201);
t323 = t196 * t202;
t352 = cos(qJ(1));
t257 = t329 * t352;
t351 = sin(qJ(1));
t151 = t351 * t199 - t202 * t257;
t256 = t329 * t351;
t153 = t352 * t199 + t202 * t256;
t361 = g(1) * t153 + g(2) * t151;
t216 = -g(3) * t323 + t361;
t237 = qJD(2) * t264;
t259 = pkin(1) * t270;
t301 = qJD(1) * qJD(2);
t281 = t202 * t301;
t262 = t196 * t281;
t280 = t196 * t300;
t268 = t199 * t237 - t202 * t259 + (t262 + t280) * pkin(7);
t77 = -pkin(2) * t236 + t268;
t368 = t77 - t216;
t307 = qJD(3) * t198;
t367 = t198 * t186 - t307;
t288 = t199 * t311;
t122 = t198 * t288 - t201 * t241;
t116 = qJD(4) + t122;
t197 = sin(qJ(4));
t200 = cos(qJ(4));
t152 = t199 * t257 + t351 * t202;
t290 = t196 * t352;
t98 = t152 * t201 - t198 * t290;
t64 = -t151 * t200 + t197 * t98;
t65 = t151 * t197 + t200 * t98;
t193 = t196 ^ 2;
t366 = 0.2e1 * t193;
t124 = t198 * t241 + t201 * t288;
t365 = t247 * t124;
t291 = pkin(1) * t329;
t324 = t196 * t199;
t232 = -pkin(7) * t324 + t202 * t291;
t219 = t198 * t236;
t285 = t201 * t309;
t299 = qJDD(1) * t201;
t205 = t201 * t220 + (t199 * t299 + (-t199 * t307 + t285) * qJD(1)) * t196 + t219;
t282 = t199 * t301;
t263 = t196 * t282;
t298 = qJDD(1) * t202;
t185 = t196 * t298;
t296 = qJDD(3) - t185;
t218 = t263 + t296;
t88 = t200 * t124 - t197 * t247;
t328 = qJD(4) * t88;
t27 = t205 * t197 - t200 * t218 + t328;
t364 = (qJDD(2) + 0.2e1 * t270) * t196;
t265 = t201 * t186;
t363 = t265 - t305;
t238 = pkin(3) * t201 + pkin(9) * t198 + pkin(2);
t303 = qJD(4) * t200;
t142 = -pkin(7) * t288 + t202 * t264;
t254 = pkin(2) * t199 - pkin(8) * t202;
t143 = t254 * t311;
t316 = t201 * t142 + t198 * t143;
t63 = pkin(9) * t288 + t316;
t362 = -t369 * t197 + t200 * t63 + t238 * t303;
t154 = -t199 * t256 + t352 * t202;
t289 = t196 * t351;
t101 = t154 * t198 - t201 * t289;
t149 = t198 * t324 - t329 * t201;
t274 = -t152 * t198 - t201 * t290;
t223 = g(1) * t101 - g(2) * t274 + g(3) * t149;
t107 = pkin(8) * t241 + t145;
t350 = pkin(8) * t199;
t239 = -pkin(2) * t202 - pkin(1) - t350;
t135 = t239 * t196;
t115 = qJD(1) * t135;
t53 = -t198 * t107 + t201 * t115;
t45 = t247 * pkin(3) - t53;
t273 = t200 * t247;
t86 = t124 * t197 + t273;
t19 = t86 * pkin(4) - t88 * qJ(5) + t45;
t354 = pkin(9) * t50;
t360 = t116 * t19 - t354;
t133 = -t329 * pkin(2) - t232;
t150 = t329 * t198 + t201 * t324;
t58 = t149 * pkin(3) - t150 * pkin(9) + t133;
t314 = pkin(7) * t323 + t199 * t291;
t134 = t329 * pkin(8) + t314;
t317 = t201 * t134 + t198 * t135;
t60 = -pkin(9) * t323 + t317;
t243 = t197 * t58 + t200 * t60;
t229 = t254 * qJD(2);
t144 = t196 * t229;
t146 = t232 * qJD(2);
t221 = -t134 * t307 + t135 * t305 + t198 * t144 + t201 * t146;
t310 = qJD(2) * t199;
t287 = t196 * t310;
t34 = pkin(9) * t287 + t221;
t147 = t314 * qJD(2);
t93 = qJD(3) * t150 + t196 * t286;
t94 = -qJD(3) * t149 + t196 * t285;
t38 = t93 * pkin(3) - t94 * pkin(9) + t147;
t358 = -qJD(4) * t243 - t197 * t34 + t200 * t38;
t357 = t88 ^ 2;
t356 = t116 ^ 2;
t203 = qJD(1) ^ 2;
t355 = pkin(4) * t50;
t345 = t88 * t86;
t78 = pkin(3) * t124 + pkin(9) * t122;
t344 = t197 * t78 + t200 * t53;
t302 = qJD(4) * t201;
t306 = qJD(3) * t200;
t343 = -qJD(5) * t201 + (-t197 * t302 - t198 * t306) * pkin(8) - t362 - t367 * qJ(5);
t304 = qJD(4) * t197;
t341 = -t238 * t304 + (-t307 * pkin(8) - t63) * t197 + (t302 * pkin(8) - t369) * t200 + t367 * pkin(4);
t113 = t197 * t265 - t200 * t288;
t319 = t200 * t202;
t126 = (t197 * t199 + t201 * t319) * t196;
t114 = qJD(1) * t126;
t248 = pkin(4) * t197 - qJ(5) * t200;
t234 = pkin(8) + t248;
t249 = pkin(4) * t200 + qJ(5) * t197;
t127 = t198 * t142;
t62 = -pkin(3) * t288 - t143 * t201 + t127;
t340 = -pkin(4) * t113 + qJ(5) * t114 - t62 + (qJD(4) * t249 - qJD(5) * t200) * t198 + t234 * t305;
t339 = pkin(9) * qJD(4);
t338 = qJ(5) * t50;
t106 = -pkin(2) * t241 - t142;
t44 = t122 * pkin(3) - t124 * pkin(9) + t106;
t54 = t201 * t107 + t198 * t115;
t46 = -t247 * pkin(9) + t54;
t21 = t197 * t44 + t200 * t46;
t15 = qJ(5) * t116 + t21;
t337 = t116 * t15;
t336 = t116 * t21;
t335 = t116 * t86;
t334 = t197 * t50;
t333 = t200 * t50;
t26 = qJD(4) * t273 + t124 * t304 - t197 * t218 - t200 * t205;
t332 = t26 * t197;
t331 = t88 * t116;
t330 = -t197 * qJD(5) + t116 * t248 - t54;
t325 = t193 * t203;
t322 = t197 * t201;
t321 = t200 * t238;
t320 = t200 * t201;
t20 = -t197 * t46 + t200 * t44;
t318 = qJD(5) - t20;
t313 = pkin(8) * t320 - t197 * t238;
t194 = t199 ^ 2;
t312 = -t202 ^ 2 + t194;
t308 = qJD(3) * t197;
t294 = t202 * t325;
t293 = t197 * t323;
t292 = pkin(7) * t185 + t199 * t259 + t202 * t237;
t284 = t116 * t304;
t283 = pkin(1) * t366;
t215 = -pkin(7) * t263 + t292;
t76 = pkin(8) * t236 + t215;
t79 = (qJD(1) * t229 + qJDD(1) * t239) * t196;
t227 = t107 * t307 - t115 * t305 - t198 * t79 - t201 * t76;
t12 = pkin(9) * t218 - t227;
t18 = t55 * pkin(3) - pkin(9) * t205 + t77;
t278 = t197 * t12 - t200 * t18 + t46 * t303 + t44 * t304;
t276 = t107 * t305 + t115 * t307 + t198 * t76 - t201 * t79;
t275 = -t198 * t134 + t135 * t201;
t272 = t202 * t247;
t271 = t116 * t200;
t269 = qJD(3) * t247;
t102 = t154 * t201 + t198 * t289;
t68 = t102 * t197 - t153 * t200;
t261 = -g(1) * t64 + g(2) * t68;
t69 = t102 * t200 + t153 * t197;
t260 = g(1) * t65 - g(2) * t69;
t258 = g(1) * t274 + g(2) * t101;
t255 = t196 * t203 * t329;
t95 = t150 * t197 + t196 * t319;
t96 = t150 * t200 - t293;
t252 = -pkin(4) * t95 + qJ(5) * t96;
t59 = pkin(3) * t323 - t275;
t251 = t197 * t305 - t113;
t250 = t200 * t305 - t114;
t14 = -pkin(4) * t116 + t318;
t246 = t14 * t200 - t15 * t197;
t244 = -t197 * t60 + t200 * t58;
t240 = 0.2e1 * t277 + qJD(2);
t235 = pkin(3) + t249;
t233 = -t134 * t305 - t135 * t307 + t144 * t201 - t198 * t146;
t231 = -t116 * t303 - t334;
t3 = t200 * t12 + t197 * t18 + t44 * t303 - t46 * t304;
t230 = t197 * t38 + t200 * t34 + t58 * t303 - t60 * t304;
t228 = t116 * t45 - t354;
t125 = -t200 * t324 + t201 * t293;
t81 = -t151 * t322 - t152 * t200;
t83 = -t153 * t322 - t154 * t200;
t226 = g(1) * t83 + g(2) * t81 + g(3) * t125;
t82 = -t151 * t320 + t152 * t197;
t84 = -t153 * t320 + t154 * t197;
t225 = -g(1) * t84 - g(2) * t82 - g(3) * t126;
t224 = t281 + t300;
t222 = -g(1) * t102 - g(2) * t98 - g(3) * t150;
t213 = g(1) * t68 + g(2) * t64 + g(3) * t95 - t278;
t35 = -pkin(3) * t287 - t233;
t212 = t116 * t339 - t223;
t13 = -pkin(3) * t218 + t276;
t5 = t27 * pkin(4) + t26 * qJ(5) - t88 * qJD(5) + t13;
t211 = -t212 - t5;
t210 = -g(1) * t69 - g(2) * t65 - g(3) * t96 + t3;
t209 = t19 * t88 + qJDD(5) - t213;
t136 = t234 * t198;
t109 = t321 + (pkin(8) * t197 + pkin(4)) * t201;
t108 = -qJ(5) * t201 + t313;
t42 = -qJD(4) * t95 + t197 * t287 + t200 * t94;
t41 = -qJD(4) * t293 + t150 * t303 + t197 * t94 - t200 * t287;
t37 = pkin(4) * t88 + qJ(5) * t86;
t28 = -t252 + t59;
t25 = -pkin(4) * t149 - t244;
t24 = qJ(5) * t149 + t243;
t23 = -pkin(4) * t124 + t197 * t53 - t200 * t78;
t22 = qJ(5) * t124 + t344;
t9 = -t26 + t335;
t8 = pkin(4) * t41 - qJ(5) * t42 - qJD(5) * t96 + t35;
t7 = -pkin(4) * t93 - t358;
t6 = qJ(5) * t93 + qJD(5) * t149 + t230;
t2 = qJDD(5) + t278 - t355;
t1 = qJD(5) * t116 + t3 + t338;
t4 = [qJDD(1), g(1) * t351 - g(2) * t352, g(1) * t352 + g(2) * t351, (qJDD(1) * t194 + 0.2e1 * t199 * t281) * t193, (t199 * t298 - t312 * t301) * t366, t196 * t240 * t309 + t364 * t199, t364 * t202 - t240 * t287, t236 * t329, -t147 * t241 + t232 * t236 - t268 * t329 + g(1) * t152 - g(2) * t154 + (-t282 + t298) * t283, -g(1) * t151 + g(2) * t153 - t146 * t241 - t215 * t329 - t224 * t283 - t314 * t236, t124 * t94 + t205 * t150, -t94 * t122 - t124 * t93 - t205 * t149 - t150 * t55, -t94 * t247 + t150 * t296 + ((-t219 + (-t220 - t262) * t201) * t202 + (-(-qJD(1) * t307 + t299) * t323 + (qJD(1) * t150 + t124) * qJD(2)) * t199) * t196, t93 * t247 - t149 * t296 + (t55 * t202 + (-qJD(1) * t149 - t122) * t310) * t196, (-t296 * t202 + (-t186 - t247) * t310) * t196, g(1) * t98 - g(2) * t102 + t106 * t93 + t147 * t122 + t133 * t55 + t77 * t149 + t275 * t218 - t233 * t247 + t276 * t323 + t53 * t287, t106 * t94 + t147 * t124 + t133 * t205 + t77 * t150 - t317 * t218 + t221 * t247 - t227 * t323 - t54 * t287 + t258, -t26 * t96 + t42 * t88, t26 * t95 - t27 * t96 - t41 * t88 - t42 * t86, t116 * t42 - t149 * t26 + t50 * t96 + t88 * t93, -t116 * t41 - t149 * t27 - t50 * t95 - t86 * t93, t116 * t93 + t149 * t50, t358 * t116 + t13 * t95 - t149 * t278 + t20 * t93 + t244 * t50 + t59 * t27 + t35 * t86 + t45 * t41 + t260, -t230 * t116 + t13 * t96 - t3 * t149 - t21 * t93 - t243 * t50 - t59 * t26 + t35 * t88 + t45 * t42 + t261, -t116 * t7 - t14 * t93 - t149 * t2 + t19 * t41 - t25 * t50 + t27 * t28 + t5 * t95 + t8 * t86 + t260, -t1 * t95 + t14 * t42 - t15 * t41 + t2 * t96 - t24 * t27 - t25 * t26 - t6 * t86 + t7 * t88 - t258, t1 * t149 + t116 * t6 + t15 * t93 - t19 * t42 + t24 * t50 + t26 * t28 - t5 * t96 - t8 * t88 - t261, t1 * t24 + t15 * t6 + t5 * t28 + t19 * t8 + t2 * t25 + t14 * t7 - g(1) * (-t351 * pkin(1) - t152 * pkin(2) - pkin(3) * t98 - pkin(4) * t65 + pkin(7) * t290 - t151 * pkin(8) + pkin(9) * t274 - qJ(5) * t64) - g(2) * (t352 * pkin(1) + t154 * pkin(2) + t102 * pkin(3) + t69 * pkin(4) + pkin(7) * t289 + t153 * pkin(8) + t101 * pkin(9) + t68 * qJ(5)); 0, 0, 0, -t199 * t294, t312 * t325, -t202 * t255 + t280, t199 * t255 + t185, t236, pkin(1) * t199 * t325 + t145 * t241 + t216 - t268, pkin(1) * t294 + t142 * t241 + g(1) * t154 + g(2) * t152 + (pkin(7) * t301 + g(3)) * t324 - t292, (-qJD(3) * t288 + t236) * t198 ^ 2 + ((t196 * t224 + t220) * t198 - t365) * t201, t363 * t122 + t367 * t124 - t198 * t55 + t201 * t205, -t201 * t269 + t198 * t296 + (t201 * t272 + (qJD(2) * t198 - t124) * t199) * t311, t198 * t269 + t201 * t296 + (-t198 * t272 + (qJD(2) * t201 + t122) * t199) * t311, t247 * t288, -pkin(2) * t55 - t127 * t247 - t53 * t288 - t145 * t122 + (-pkin(8) * t218 - t106 * t247) * t198 + (pkin(8) * t269 + t143 * t247 - t368) * t201, -pkin(2) * t205 - t145 * t124 - t316 * t247 + t54 * t288 - t363 * t106 + (-t201 * t218 - t247 * t307) * pkin(8) + t368 * t198, -t198 * t200 * t26 + (-t198 * t304 + t250) * t88, t113 * t88 + t114 * t86 + (-t197 * t88 - t200 * t86) * t305 + (t332 - t200 * t27 + (t197 * t86 - t200 * t88) * qJD(4)) * t198, t201 * t26 + t250 * t116 + (-t247 * t88 - t284 + t333) * t198, t201 * t27 - t251 * t116 + (t247 * t86 + t231) * t198, -t116 * t198 * t247 - t50 * t201, -t50 * t321 - t45 * t113 - t62 * t86 + (t369 * t200 + (qJD(4) * t238 + t63) * t197) * t116 + (t45 * t308 + t278 + (qJD(3) * t86 + t231) * pkin(8)) * t201 + (t45 * t303 + t13 * t197 - t247 * t20 + (t116 * t308 + t27) * pkin(8)) * t198 + t225, -t313 * t50 - t62 * t88 - t45 * t114 + t362 * t116 + (t45 * t306 + t3 + (qJD(3) * t88 + t284) * pkin(8)) * t201 + (-t45 * t304 + t13 * t200 + t247 * t21 + (t116 * t306 - t26) * pkin(8)) * t198 + t226, -t109 * t50 + t136 * t27 + t2 * t201 + t340 * t86 + t251 * t19 - t341 * t116 + (t14 * t247 + t19 * t303 + t197 * t5) * t198 + t225, -t108 * t27 - t109 * t26 + t113 * t15 - t114 * t14 + t341 * t88 - t343 * t86 + t246 * t305 + (-t1 * t197 + t2 * t200 + (-t14 * t197 - t15 * t200) * qJD(4) + t216) * t198, -t1 * t201 + t108 * t50 + t136 * t26 - t340 * t88 - t250 * t19 + t343 * t116 + (-t15 * t247 + t19 * t304 - t200 * t5) * t198 - t226, t1 * t108 + t5 * t136 + t2 * t109 - g(1) * (pkin(4) * t84 + pkin(8) * t154 + qJ(5) * t83) - g(2) * (pkin(4) * t82 + pkin(8) * t152 + qJ(5) * t81) + t340 * t19 + t343 * t15 + t341 * t14 + t361 * t238 + (-pkin(4) * t126 - qJ(5) * t125 - (t202 * t238 + t350) * t196) * g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124 * t122, -t122 ^ 2 + t124 ^ 2, -t122 * t247 + t205, -t55 - t365, t218, -t106 * t124 - t54 * t247 + t223 - t276, t106 * t122 - t53 * t247 - t222 + t227, t271 * t88 - t332, (-t26 - t335) * t200 + (-t27 - t331) * t197, t116 * t271 - t88 * t124 + t334, t86 * t124 - t356 * t197 + t333, -t116 * t124, -pkin(3) * t27 - t20 * t124 - t54 * t86 + (t53 * t116 + t228) * t197 + (-t13 + (-t78 - t339) * t116 + t223) * t200, pkin(3) * t26 + t344 * t116 + t21 * t124 - t54 * t88 + t228 * t200 + (t13 + t212) * t197, t116 * t23 + t124 * t14 + t360 * t197 + t211 * t200 - t235 * t27 + t330 * t86, t22 * t86 - t23 * t88 + (t1 + t116 * t14 + (-t27 + t328) * pkin(9)) * t200 + (t2 - t337 + (qJD(4) * t86 - t26) * pkin(9)) * t197 + t222, -t116 * t22 - t124 * t15 + t211 * t197 - t360 * t200 - t235 * t26 - t330 * t88, -t14 * t23 - t15 * t22 + t330 * t19 + (qJD(4) * t246 + t1 * t200 + t2 * t197 + t222) * pkin(9) + (-t5 + t223) * t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t345, -t86 ^ 2 + t357, t9, t331 - t27, t50, -t45 * t88 + t213 + t336, t116 * t20 + t45 * t86 - t210, -t37 * t86 - t209 + t336 + 0.2e1 * t355, pkin(4) * t26 - qJ(5) * t27 + (t15 - t21) * t88 + (t14 - t318) * t86, 0.2e1 * t338 - t19 * t86 + t37 * t88 + (0.2e1 * qJD(5) - t20) * t116 + t210, t1 * qJ(5) - t2 * pkin(4) - t19 * t37 - t14 * t21 - g(1) * (-pkin(4) * t68 + qJ(5) * t69) - g(2) * (-pkin(4) * t64 + qJ(5) * t65) - g(3) * t252 + t318 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t345 - t50, t9, -t356 - t357, t209 - t337 - t355;];
tau_reg = t4;
