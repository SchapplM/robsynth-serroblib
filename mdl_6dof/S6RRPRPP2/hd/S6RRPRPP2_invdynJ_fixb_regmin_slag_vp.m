% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRPP2
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
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRPP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:52:32
% EndTime: 2019-03-09 09:52:42
% DurationCPUTime: 4.13s
% Computational Cost: add. (6797->530), mult. (15485->625), div. (0->0), fcn. (11135->10), ass. (0->260)
t201 = cos(qJ(2));
t325 = cos(pkin(9));
t263 = t325 * t201;
t173 = qJD(1) * t263;
t195 = sin(pkin(9));
t198 = sin(qJ(2));
t290 = qJD(1) * t198;
t136 = -t195 * t290 + t173;
t128 = qJD(4) - t136;
t151 = t195 * t201 + t198 * t325;
t137 = t151 * qJD(2);
t282 = t198 * qJDD(1);
t247 = qJDD(1) * t263 - t195 * t282;
t86 = qJD(1) * t137 + qJDD(4) - t247;
t371 = t86 * qJ(5) + t128 * qJD(5);
t138 = t151 * qJD(1);
t197 = sin(qJ(4));
t200 = cos(qJ(4));
t107 = qJD(2) * t197 + t138 * t200;
t315 = t128 * t197;
t370 = t107 * t315;
t284 = qJD(1) * qJD(2);
t268 = t198 * t284;
t210 = qJD(2) * t173 + qJDD(1) * t151 - t195 * t268;
t283 = qJD(2) * qJD(4);
t288 = qJD(4) * t197;
t228 = t197 * qJDD(2) - t138 * t288 + (t210 + t283) * t200;
t105 = -t200 * qJD(2) + t138 * t197;
t321 = t105 * t128;
t21 = t228 + t321;
t191 = qJ(2) + pkin(9);
t184 = sin(t191);
t202 = cos(qJ(1));
t309 = t184 * t202;
t199 = sin(qJ(1));
t345 = g(2) * t199;
t369 = g(1) * t309 + t184 * t345;
t189 = t201 * pkin(2);
t183 = t189 + pkin(1);
t207 = -t200 * qJDD(2) + t197 * t210;
t50 = qJD(4) * t107 + t207;
t368 = t50 * qJ(6) + t105 * qJD(6);
t367 = 0.2e1 * t371;
t196 = -qJ(3) - pkin(7);
t160 = t196 * t198;
t153 = qJD(1) * t160;
t161 = t196 * t201;
t154 = qJD(1) * t161;
t264 = t325 * t154;
t94 = t153 * t195 - t264;
t366 = qJD(5) * t197 + t94;
t252 = g(1) * t202 + t345;
t365 = t184 * t252;
t352 = pkin(4) + pkin(5);
t364 = t200 * t352;
t359 = g(1) * t199 - g(2) * t202;
t363 = t359 * t184;
t159 = -qJD(1) * t183 + qJD(3);
t57 = -pkin(3) * t136 - pkin(8) * t138 + t159;
t335 = qJD(2) * pkin(2);
t145 = t153 + t335;
t91 = t195 * t145 - t264;
t78 = qJD(2) * pkin(8) + t91;
t30 = -t197 * t78 + t200 * t57;
t295 = qJD(5) - t30;
t104 = t107 ^ 2;
t125 = t128 ^ 2;
t362 = -t125 - t104;
t229 = -t195 * t198 + t263;
t361 = t137 * qJ(5) - qJD(5) * t229;
t185 = cos(t191);
t360 = -t185 * pkin(3) - t184 * pkin(8);
t314 = t136 * t197;
t242 = t200 * t86 + (-t288 + t314) * t128;
t319 = t105 * t138;
t216 = t242 - t319;
t81 = t86 * pkin(4);
t358 = t81 - qJDD(5);
t18 = qJ(6) * t107 + t30;
t296 = qJD(5) - t18;
t12 = -t128 * t352 + t296;
t287 = qJD(4) * t200;
t280 = pkin(2) * t268 + qJDD(3);
t281 = t201 * qJDD(1);
t322 = qJDD(1) * pkin(1);
t92 = -qJD(2) * t138 + t247;
t33 = -pkin(2) * t281 - t92 * pkin(3) - pkin(8) * t210 + t280 - t322;
t265 = qJD(2) * t196;
t134 = -qJD(3) * t198 + t201 * t265;
t85 = qJDD(2) * pkin(2) + qJD(1) * t134 + qJDD(1) * t160;
t133 = qJD(3) * t201 + t198 * t265;
t93 = qJD(1) * t133 - qJDD(1) * t161;
t43 = t195 * t85 + t325 * t93;
t39 = qJDD(2) * pkin(8) + t43;
t233 = t197 * t33 + t200 * t39 + t57 * t287 - t288 * t78;
t3 = t233 + t371;
t2 = t3 + t368;
t357 = t12 * t128 + t2;
t141 = t195 * t154;
t90 = t325 * t145 + t141;
t255 = qJD(2) * pkin(3) + t90;
t223 = qJ(5) * t107 + t255;
t32 = pkin(4) * t105 - t223;
t350 = pkin(2) * t195;
t178 = pkin(8) + t350;
t331 = t178 * t86;
t356 = t128 * t32 - t331;
t17 = -t105 * t352 + qJD(6) + t223;
t300 = t200 * t202;
t303 = t197 * t199;
t129 = t185 * t303 + t300;
t298 = t202 * t197;
t302 = t199 * t200;
t131 = t185 * t298 - t302;
t262 = t197 * t39 - t200 * t33 + t78 * t287 + t57 * t288;
t311 = t184 * t197;
t219 = g(1) * t131 + g(2) * t129 + g(3) * t311 - t262;
t215 = t219 + t358;
t336 = qJ(6) * t228;
t355 = (qJD(6) + t17) * t107 + t215 + t336;
t354 = -g(3) * t184 - t252 * t185;
t353 = t105 ^ 2;
t351 = pkin(5) * t86;
t349 = pkin(2) * t198;
t346 = g(2) * t196;
t343 = g(3) * t185;
t342 = g(3) * t201;
t341 = t200 * pkin(4);
t340 = -t105 * t287 - t197 * t50;
t31 = t197 * t57 + t200 * t78;
t42 = -t195 * t93 + t325 * t85;
t67 = pkin(2) * t290 + pkin(3) * t138 - pkin(8) * t136;
t95 = t153 * t325 + t141;
t339 = t197 * t67 + t200 * t95;
t101 = t195 * t160 - t161 * t325;
t89 = -pkin(3) * t229 - pkin(8) * t151 - t183;
t338 = t200 * t101 + t197 * t89;
t337 = qJ(5) * t50;
t117 = t128 * qJ(5);
t19 = qJ(6) * t105 + t31;
t14 = t117 + t19;
t334 = t128 * t14;
t23 = t117 + t31;
t333 = t128 * t23;
t332 = t128 * t31;
t330 = t197 * t228;
t73 = t197 * t86;
t274 = t352 * t197;
t323 = qJ(5) * t200;
t237 = -t274 + t323;
t329 = t128 * t237 + t366;
t25 = t138 * qJ(5) + t339;
t297 = qJ(6) - t178;
t328 = -qJ(6) * t314 - qJD(6) * t200 + t288 * t297 - t25;
t149 = t297 * t200;
t87 = t197 * t95;
t327 = -qJD(4) * t149 - qJD(6) * t197 - t87 - (-qJ(6) * t136 - t67) * t200 + t352 * t138;
t248 = pkin(4) * t197 - t323;
t326 = t128 * t248 - t366;
t324 = qJ(5) * t105;
t320 = t105 * t136;
t318 = t107 * t105;
t317 = t107 * t128;
t316 = t107 * t138;
t259 = t128 * t200;
t313 = t151 * t197;
t312 = t151 * t200;
t310 = t184 * t200;
t308 = t185 * t199;
t307 = t185 * t200;
t306 = t185 * t202;
t305 = t197 * qJ(5);
t140 = t229 * qJD(2);
t304 = t197 * t140;
t301 = t200 * t140;
t299 = t202 * t196;
t293 = t369 * t197;
t292 = t369 * t200;
t192 = t198 ^ 2;
t291 = -t201 ^ 2 + t192;
t289 = qJD(4) * t178;
t285 = qJD(5) * t200;
t64 = t133 * t325 + t195 * t134;
t279 = t101 * t287 + t197 * t64 + t89 * t288;
t275 = t198 * t335;
t68 = pkin(3) * t137 - pkin(8) * t140 + t275;
t278 = t197 * t68 + t200 * t64 + t89 * t287;
t38 = -qJDD(2) * pkin(3) - t42;
t34 = -qJ(5) * t229 + t338;
t8 = t50 * pkin(4) - qJ(5) * t228 - t107 * qJD(5) + t38;
t4 = -pkin(5) * t50 + qJDD(6) - t8;
t272 = t4 - t343;
t271 = t325 * pkin(2);
t270 = t151 * t288;
t269 = t151 * t287;
t97 = t197 * t101;
t266 = t200 * t89 - t97;
t130 = t185 * t302 - t298;
t261 = -t129 * pkin(4) + qJ(5) * t130;
t132 = t185 * t300 + t303;
t260 = -t131 * pkin(4) + qJ(5) * t132;
t63 = t133 * t195 - t325 * t134;
t100 = -t325 * t160 - t161 * t195;
t257 = pkin(8) * t308 - t199 * t349;
t256 = pkin(8) * t306 - t202 * t349;
t179 = -t271 - pkin(3);
t254 = g(1) * t129 - g(2) * t131;
t253 = g(1) * t130 - g(2) * t132;
t5 = t262 - t358;
t1 = -qJD(6) * t107 - t336 - t351 + t5;
t250 = -t1 + t334;
t249 = t305 + t341;
t22 = -pkin(4) * t128 + t295;
t246 = -t197 * t23 + t200 * t22;
t245 = pkin(4) * t307 + t185 * t305 + t189 - t360;
t244 = t200 * t68 - t279;
t243 = -qJ(6) * t140 - qJD(6) * t151;
t241 = t128 * t287 - t136 * t259 + t73;
t239 = -t183 + t360;
t238 = t128 * t289 + t343;
t236 = -t305 - t364;
t235 = -0.2e1 * pkin(1) * t284 - pkin(7) * qJDD(2);
t234 = -t130 * pkin(4) - t129 * qJ(5) - t299;
t232 = -t270 + t301;
t231 = -t101 * t288 + t278;
t230 = -t128 * t255 - t331;
t171 = t202 * t183;
t227 = pkin(3) * t306 + t132 * pkin(4) + pkin(8) * t309 + t131 * qJ(5) + t171;
t226 = -t238 - t8;
t224 = t179 - t305;
t221 = -qJDD(1) * t183 + t280;
t220 = t241 + t316;
t204 = qJD(2) ^ 2;
t218 = -pkin(7) * t204 + 0.2e1 * t322 + t359;
t205 = qJD(1) ^ 2;
t217 = pkin(1) * t205 - pkin(7) * qJDD(1) + t252;
t213 = t107 * t32 - t215;
t212 = g(1) * t132 + g(2) * t130 + g(3) * t310 - t233;
t209 = t128 * t30 + t212;
t208 = -qJDD(4) + t92 + t318;
t206 = -t138 * t287 - t197 * t283 - t207;
t162 = qJ(5) * t310;
t148 = t297 * t197;
t147 = t224 - t341;
t118 = -t224 + t364;
t52 = pkin(4) * t107 + t324;
t51 = t151 * t248 + t100;
t41 = t151 * t237 - t100;
t40 = -t107 * t352 - t324;
t35 = pkin(4) * t229 - t266;
t26 = -pkin(4) * t138 - t200 * t67 + t87;
t24 = qJ(6) * t313 + t34;
t20 = t97 + (-qJ(6) * t151 - t89) * t200 + t352 * t229;
t13 = t248 * t140 + (qJD(4) * t249 - t285) * t151 + t63;
t11 = t237 * t140 + (qJD(4) * t236 + t285) * t151 - t63;
t10 = -pkin(4) * t137 - t244;
t9 = t231 + t361;
t7 = qJ(6) * t269 + (-qJD(4) * t101 - t243) * t197 + t278 + t361;
t6 = qJ(6) * t270 - t352 * t137 + (t243 - t68) * t200 + t279;
t15 = [qJDD(1), t359, t252, qJDD(1) * t192 + 0.2e1 * t201 * t268, 0.2e1 * t198 * t281 - 0.2e1 * t284 * t291, qJDD(2) * t198 + t201 * t204, qJDD(2) * t201 - t198 * t204, 0, t198 * t235 + t201 * t218, -t198 * t218 + t201 * t235, t100 * t210 + t101 * t92 + t64 * t136 - t91 * t137 + t63 * t138 - t90 * t140 - t42 * t151 + t229 * t43 - t252, t43 * t101 + t91 * t64 - t42 * t100 - t90 * t63 - t221 * t183 + t159 * t275 - g(1) * (-t199 * t183 - t299) - g(2) * (-t199 * t196 + t171) t107 * t232 + t228 * t312 (-t105 * t200 - t107 * t197) * t140 + (-t330 - t200 * t50 + (t105 * t197 - t107 * t200) * qJD(4)) * t151, t107 * t137 + t128 * t232 - t228 * t229 + t312 * t86, -t86 * t313 - t105 * t137 + t229 * t50 + (-t269 - t304) * t128, t128 * t137 - t229 * t86, t244 * t128 + t266 * t86 + t262 * t229 + t30 * t137 + t63 * t105 + t100 * t50 - t255 * t304 + (t38 * t197 - t255 * t287) * t151 + t253, -t231 * t128 - t338 * t86 + t233 * t229 - t31 * t137 + t63 * t107 + t100 * t228 - t255 * t301 + (t38 * t200 + t255 * t288) * t151 - t254, t32 * t304 - t10 * t128 + t105 * t13 - t137 * t22 + t229 * t5 - t35 * t86 + t50 * t51 + (t197 * t8 + t287 * t32) * t151 + t253, t10 * t107 - t105 * t9 - t34 * t50 + t35 * t228 + t363 + t246 * t140 + (-t197 * t3 + t200 * t5 + (-t197 * t22 - t200 * t23) * qJD(4)) * t151, -t32 * t301 - t107 * t13 + t128 * t9 + t137 * t23 - t229 * t3 + t34 * t86 - t228 * t51 + (-t200 * t8 + t288 * t32) * t151 + t254, t3 * t34 + t23 * t9 + t8 * t51 + t32 * t13 + t5 * t35 + t22 * t10 - g(1) * t234 - g(2) * t227 + (-g(1) * t239 + t346) * t199, -t17 * t304 + t1 * t229 - t105 * t11 - t12 * t137 - t128 * t6 - t20 * t86 - t41 * t50 + (-t17 * t287 - t197 * t4) * t151 + t253, t17 * t301 + t107 * t11 + t128 * t7 + t137 * t14 - t229 * t2 + t24 * t86 + t41 * t228 + (-t17 * t288 + t200 * t4) * t151 + t254, t105 * t7 - t107 * t6 - t20 * t228 + t24 * t50 - t363 + (-t12 * t200 + t14 * t197) * t140 + (-t1 * t200 + t197 * t2 + (t12 * t197 + t14 * t200) * qJD(4)) * t151, t2 * t24 + t14 * t7 + t1 * t20 + t12 * t6 + t4 * t41 + t17 * t11 - g(1) * (-pkin(5) * t130 + t234) - g(2) * (pkin(5) * t132 - qJ(6) * t309 + t227) + (-g(1) * (qJ(6) * t184 + t239) + t346) * t199; 0, 0, 0, -t198 * t205 * t201, t291 * t205, t282, t281, qJDD(2), t198 * t217 - t342, g(3) * t198 + t201 * t217, -t210 * t271 + t92 * t350 - (-t91 + t94) * t138 + (t90 - t95) * t136, t90 * t94 - t91 * t95 + (t325 * t42 - t342 + t195 * t43 + (-qJD(1) * t159 + t252) * t198) * pkin(2), t107 * t259 + t330 (t228 + t320) * t200 - t370 + t340, t241 - t316, t242 + t319, -t128 * t138, -t94 * t105 + t87 * t128 - t30 * t138 + t179 * t50 + (-t343 - t38 + (-t67 - t289) * t128) * t200 + t230 * t197 + t292, t179 * t228 + t339 * t128 + t31 * t138 - t94 * t107 + t230 * t200 + (t238 + t38) * t197 - t293, t326 * t105 + t128 * t26 + t138 * t22 + t147 * t50 + t356 * t197 + t226 * t200 + t292, t105 * t25 - t107 * t26 + (-t136 * t22 - t178 * t50 + t3 + (t107 * t178 + t22) * qJD(4)) * t200 + (t136 * t23 + t178 * t228 + t5 + (t105 * t178 - t23) * qJD(4)) * t197 + t354, -t326 * t107 - t128 * t25 - t138 * t23 - t147 * t228 + t226 * t197 - t200 * t356 + t293, t8 * t147 - t23 * t25 - t22 * t26 - g(1) * t256 - g(2) * t257 - g(3) * t245 + t326 * t32 + (qJD(4) * t246 + t5 * t197 + t3 * t200) * t178 + (pkin(3) + t249) * t365, -t105 * t329 - t118 * t50 + t12 * t138 - t128 * t327 + t148 * t86 - t17 * t315 + t200 * t272 + t292, t107 * t329 + t118 * t228 + t128 * t328 - t138 * t14 - t149 * t86 + t17 * t259 + t197 * t272 + t293, t328 * t105 - t327 * t107 + t148 * t228 - t149 * t50 + t250 * t197 - t357 * t200 - t354, -t2 * t149 - t1 * t148 + t4 * t118 - g(1) * (-qJ(6) * t306 + t256) - g(2) * (-qJ(6) * t308 + t257) - g(3) * (pkin(5) * t307 + t245) + t329 * t17 + t328 * t14 + t327 * t12 + (g(3) * qJ(6) + t252 * (pkin(3) - t236)) * t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136 ^ 2 - t138 ^ 2, -t136 * t91 + t138 * t90 + t221 - t359, 0, 0, 0, 0, 0, t216, -t125 * t200 - t316 - t73, t216 (-t228 + t320) * t200 + t370 + t340, t220, -t138 * t32 + (-t5 + t333) * t200 + (t128 * t22 + t3) * t197 - t359, t216, t220, t21 * t200 + (t50 - t317) * t197, t138 * t17 + t357 * t197 + t250 * t200 - t359; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t318, t104 - t353, t21, t206 + t317, t86, t107 * t255 + t219 + t332, -t105 * t255 + t209, -t105 * t52 - t213 + t332 + t81, -pkin(4) * t228 - t337 + (t23 - t31) * t107 + (t22 - t295) * t105, -t105 * t32 + t107 * t52 - t209 + t367, t3 * qJ(5) - t5 * pkin(4) - t32 * t52 - t22 * t31 - g(1) * t260 - g(2) * t261 - g(3) * (-pkin(4) * t311 + t162) + t295 * t23 (pkin(5) + t352) * t86 + t105 * t40 + t128 * t19 + t355, t105 * t17 - t107 * t40 - t128 * t18 - t212 + t367 + t368, t337 + t352 * t228 + (-t14 + t19) * t107 + (-t12 + t296) * t105, t2 * qJ(5) - t1 * t352 - t12 * t19 - t17 * t40 - g(1) * (-pkin(5) * t131 + t260) - g(2) * (-pkin(5) * t129 + t261) - g(3) * (-t184 * t274 + t162) + t296 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, t21, t362, t213 - t333, t208, t362, -t21, -t334 - t351 - t355; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t206 - t317, t228 - t321, -t104 - t353, -t105 * t14 + t107 * t12 + t272 + t365;];
tau_reg  = t15;
