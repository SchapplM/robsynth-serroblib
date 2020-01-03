% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRPP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:02:06
% EndTime: 2019-12-31 21:02:18
% DurationCPUTime: 6.45s
% Computational Cost: add. (5590->536), mult. (12858->676), div. (0->0), fcn. (8741->10), ass. (0->257)
t202 = cos(qJ(2));
t286 = qJD(1) * t202;
t165 = -qJD(3) + t286;
t199 = sin(qJ(2));
t276 = t199 * qJDD(1);
t178 = pkin(6) * t276;
t277 = qJD(1) * qJD(2);
t262 = t202 * t277;
t109 = -qJDD(2) * pkin(2) + pkin(6) * t262 + t178;
t192 = g(3) * t202;
t200 = sin(qJ(1));
t203 = cos(qJ(1));
t247 = g(1) * t203 + g(2) * t200;
t216 = t247 * t199 - t192;
t210 = -t109 + t216;
t359 = pkin(7) * qJD(3) * t165 + t210;
t186 = t202 * qJDD(1);
t353 = -t199 * t277 + t186;
t126 = qJDD(3) - t353;
t193 = qJ(3) + pkin(8);
t182 = sin(t193);
t197 = -qJ(4) - pkin(7);
t201 = cos(qJ(3));
t147 = t197 * t201;
t196 = sin(pkin(8));
t198 = sin(qJ(3));
t308 = t196 * t198;
t316 = cos(pkin(8));
t88 = -t316 * t147 + t197 * t308;
t358 = -t126 * t88 - t216 * t182;
t287 = qJD(1) * t199;
t271 = t198 * t287;
t278 = t201 * qJD(2);
t132 = t271 - t278;
t269 = t201 * t287;
t285 = qJD(2) * t198;
t134 = t269 + t285;
t73 = t316 * t132 + t134 * t196;
t319 = t73 * t165;
t281 = qJD(3) * t199;
t355 = qJD(1) * t281 - qJDD(2);
t67 = -qJD(3) * t278 + (-t262 - t276) * t201 + t355 * t198;
t68 = ((qJD(3) + t286) * qJD(2) + t276) * t198 + t355 * t201;
t36 = -t196 * t68 - t316 * t67;
t20 = t36 - t319;
t225 = -t196 * t132 + t316 * t134;
t357 = t73 * t225;
t259 = t316 * t198;
t128 = t196 * t201 + t259;
t112 = t128 * qJD(3);
t318 = -t128 * t286 + t112;
t258 = t316 * t201;
t282 = qJD(3) * t198;
t113 = qJD(3) * t258 - t196 * t282;
t270 = t198 * t286;
t317 = t196 * t270 - t258 * t286 + t113;
t283 = qJD(2) * t202;
t268 = t198 * t283;
t280 = qJD(3) * t201;
t356 = t199 * t280 + t268;
t336 = t225 ^ 2;
t260 = qJD(3) * t197;
t279 = qJD(4) * t201;
t107 = t198 * t260 + t279;
t220 = -qJD(4) * t198 + t201 * t260;
t297 = t201 * t202;
t234 = pkin(3) * t199 - qJ(4) * t297;
t249 = pkin(2) * t199 - pkin(7) * t202;
t135 = t249 * qJD(1);
t90 = pkin(6) * t271 + t201 * t135;
t66 = t234 * qJD(1) + t90;
t118 = t198 * t135;
t302 = t199 * t201;
t305 = t198 * t202;
t77 = t118 + (-pkin(6) * t302 - qJ(4) * t305) * qJD(1);
t331 = (t220 - t66) * t316 + (-t107 + t77) * t196;
t108 = t353 * pkin(6) + qJDD(2) * pkin(7);
t250 = pkin(2) * t202 + pkin(7) * t199;
t146 = -pkin(1) - t250;
t123 = t146 * qJD(1);
t180 = pkin(6) * t286;
t149 = qJD(2) * pkin(7) + t180;
t136 = t249 * qJD(2);
t84 = qJD(1) * t136 + t146 * qJDD(1);
t30 = t201 * t108 + t123 * t280 - t149 * t282 + t198 * t84;
t82 = t201 * t123 - t149 * t198;
t354 = t165 * t82 + t30;
t254 = -t180 + (-t270 + t282) * pkin(3);
t170 = pkin(6) * t297;
t96 = t198 * t146 + t170;
t35 = -t196 * t67 + t316 * t68;
t352 = pkin(4) * t35 - qJ(5) * t36 - qJD(5) * t225;
t127 = -t258 + t308;
t351 = -t126 * t127 + t165 * t318 + t73 * t287;
t338 = g(3) * t199;
t350 = t247 * t202 + t338;
t148 = -qJD(2) * pkin(2) + pkin(6) * t287;
t89 = pkin(3) * t132 + qJD(4) + t148;
t29 = pkin(4) * t73 - qJ(5) * t225 + t89;
t348 = -t29 * t225 - qJDD(5);
t347 = t127 * t36 + t128 * t35 + t225 * t318 + t317 * t73;
t346 = -0.2e1 * pkin(1);
t344 = pkin(4) * t126;
t343 = pkin(6) * t198;
t342 = g(1) * t200;
t339 = g(2) * t203;
t83 = t123 * t198 + t149 * t201;
t59 = -qJ(4) * t132 + t83;
t54 = t316 * t59;
t58 = -qJ(4) * t134 + t82;
t27 = t196 * t58 + t54;
t337 = t27 * t225;
t334 = t73 ^ 2;
t79 = t201 * t84;
t31 = -qJD(3) * t83 - t108 * t198 + t79;
t15 = pkin(3) * t126 + qJ(4) * t67 - qJD(4) * t134 + t31;
t18 = -qJ(4) * t68 - qJD(4) * t132 + t30;
t3 = t316 * t15 - t196 * t18;
t4 = t196 * t15 + t316 * t18;
t284 = qJD(2) * t199;
t291 = t201 * t136 + t284 * t343;
t40 = -t199 * t279 + t234 * qJD(2) + (-t170 + (qJ(4) * t199 - t146) * t198) * qJD(3) + t291;
t292 = t198 * t136 + t146 * t280;
t45 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t302 + (-qJD(4) * t199 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t202) * t198 + t292;
t14 = t196 * t40 + t316 * t45;
t332 = t318 * pkin(4) - t317 * qJ(5) - qJD(5) * t128 + t254;
t52 = -pkin(3) * t165 + t58;
t26 = t196 * t52 + t54;
t330 = pkin(4) * t287 - t331;
t39 = t196 * t66 + t316 * t77;
t33 = qJ(5) * t287 + t39;
t63 = t316 * t107 + t196 * t220;
t329 = t63 - t33;
t328 = t63 - t39;
t130 = t201 * t146;
t80 = -qJ(4) * t302 + t130 + (-pkin(3) - t343) * t202;
t307 = t198 * t199;
t86 = -qJ(4) * t307 + t96;
t47 = t196 * t80 + t316 * t86;
t325 = t165 * t225;
t323 = t165 * t83;
t322 = t196 * t59;
t321 = t198 * t67;
t320 = t201 * t68;
t315 = pkin(6) * qJDD(1);
t313 = t132 * t165;
t312 = t132 * t198;
t311 = t134 * t132;
t310 = t134 * t165;
t309 = t134 * t201;
t306 = t198 * t200;
t304 = t198 * t203;
t303 = t199 * t200;
t301 = t199 * t203;
t183 = cos(t193);
t300 = t200 * t183;
t299 = t200 * t201;
t298 = t200 * t202;
t296 = t201 * t203;
t177 = pkin(3) * t201 + pkin(2);
t151 = t202 * t177;
t295 = t202 * t203;
t294 = t203 * t182;
t28 = t316 * t58 - t322;
t293 = qJD(5) - t28;
t167 = pkin(3) * t307;
t137 = t199 * pkin(6) + t167;
t290 = t203 * pkin(1) + t200 * pkin(6);
t194 = t199 ^ 2;
t195 = t202 ^ 2;
t289 = t194 - t195;
t288 = t194 + t195;
t273 = t198 * t295;
t205 = qJD(1) ^ 2;
t272 = t199 * t205 * t202;
t92 = t356 * pkin(3) + pkin(6) * t283;
t267 = t202 * t278;
t265 = t165 * t287;
t256 = t199 * t262;
t172 = g(1) * t303;
t255 = -g(2) * t301 + t172;
t101 = t202 * t294 - t300;
t99 = t182 * t298 + t183 * t203;
t253 = g(1) * t99 - g(2) * t101;
t252 = t316 * t283;
t100 = t183 * t298 - t294;
t102 = t182 * t200 + t183 * t295;
t248 = g(1) * t100 - g(2) * t102;
t246 = -t334 - t336;
t245 = -t334 + t336;
t104 = t128 * t199;
t60 = -t113 * t199 - t196 * t267 - t198 * t252;
t244 = t104 * t35 - t60 * t73;
t243 = pkin(4) * t183 + qJ(5) * t182;
t242 = pkin(6) * t132 + t148 * t198;
t241 = pkin(6) * t134 + t148 * t201;
t239 = -t198 * t83 - t201 * t82;
t238 = -pkin(7) * t126 + qJD(3) * t148;
t236 = t35 - t325;
t235 = t35 + t325;
t230 = -pkin(6) * qJDD(2) + t277 * t346;
t114 = t198 * t298 + t296;
t228 = t127 * t35 + t318 * t73;
t13 = -t196 * t45 + t316 * t40;
t25 = t316 * t52 - t322;
t46 = -t196 * t86 + t316 * t80;
t227 = t126 * t198 - t165 * t280;
t226 = t126 * t201 + t165 * t282;
t224 = pkin(1) * t205 + t247;
t204 = qJD(2) ^ 2;
t223 = pkin(6) * t204 + qJDD(1) * t346 + t339;
t222 = g(1) * t101 + g(2) * t99 + t182 * t338 + t3;
t221 = pkin(3) * t306 + t177 * t295 - t197 * t301 + t290;
t87 = -t147 * t196 - t197 * t259;
t219 = g(2) * t199 * t300 - t126 * t87 + (g(1) * t301 - t192) * t183;
t51 = pkin(3) * t68 + qJDD(4) + t109;
t190 = t203 * pkin(6);
t218 = t197 * t303 + pkin(3) * t304 + t190 + (-pkin(1) - t151) * t200;
t217 = -t36 - t319;
t105 = -t196 * t307 + t199 * t258;
t61 = t199 * t112 + t196 * t268 - t201 * t252;
t213 = t104 * t36 + t105 * t35 - t225 * t60 - t61 * t73;
t212 = g(1) * t102 + g(2) * t100 + t183 * t338 - t4;
t211 = -t27 * t165 + t222;
t209 = t104 * t126 + t165 * t60 - t202 * t35 + t73 * t284;
t208 = -t88 * t35 + t36 * t87 - t63 * t73 - t350;
t207 = t51 - t216;
t175 = -t316 * pkin(3) - pkin(4);
t173 = pkin(3) * t196 + qJ(5);
t169 = pkin(3) * t299;
t122 = t126 * qJ(5);
t117 = t201 * t295 + t306;
t116 = -t273 + t299;
t115 = -t200 * t297 + t304;
t95 = -pkin(6) * t305 + t130;
t91 = -pkin(6) * t269 + t118;
t85 = -t126 * t202 - t165 * t284;
t70 = pkin(4) * t127 - qJ(5) * t128 - t177;
t57 = -t96 * qJD(3) + t291;
t56 = (-t199 * t278 - t202 * t282) * pkin(6) + t292;
t53 = pkin(4) * t104 - qJ(5) * t105 + t137;
t44 = t202 * pkin(4) - t46;
t43 = -qJ(5) * t202 + t47;
t32 = pkin(3) * t134 + pkin(4) * t225 + qJ(5) * t73;
t23 = -qJ(5) * t165 + t26;
t22 = t165 * pkin(4) + qJD(5) - t25;
t21 = t126 * t128 - t317 * t165 - t225 * t287;
t19 = -pkin(4) * t60 + qJ(5) * t61 - qJD(5) * t105 + t92;
t10 = -pkin(4) * t284 - t13;
t9 = qJ(5) * t284 - qJD(5) * t202 + t14;
t8 = t105 * t36 - t225 * t61;
t7 = t128 * t36 + t225 * t317;
t6 = t105 * t126 + t165 * t61 - t202 * t36 + t225 * t284;
t5 = t51 + t352;
t2 = qJDD(5) - t3 - t344;
t1 = -qJD(5) * t165 + t122 + t4;
t11 = [0, 0, 0, 0, 0, qJDD(1), -t339 + t342, t247, 0, 0, qJDD(1) * t194 + 0.2e1 * t256, 0.2e1 * t199 * t186 - 0.2e1 * t289 * t277, qJDD(2) * t199 + t202 * t204, qJDD(1) * t195 - 0.2e1 * t256, qJDD(2) * t202 - t199 * t204, 0, t230 * t199 + (-t223 + t342) * t202, t223 * t199 + t230 * t202 - t172, 0.2e1 * t288 * t315 - t247, -g(1) * (-pkin(1) * t200 + t190) - g(2) * t290 + (pkin(6) ^ 2 * t288 + pkin(1) ^ 2) * qJDD(1), -t67 * t302 + (-t198 * t281 + t267) * t134, (-t132 * t201 - t134 * t198) * t283 + (t321 - t320 + (-t309 + t312) * qJD(3)) * t199, (-t165 * t278 + t67) * t202 + (qJD(2) * t134 + t226) * t199, t356 * t132 + t68 * t307, (t165 * t285 + t68) * t202 + (-qJD(2) * t132 - t227) * t199, t85, -g(1) * t115 - g(2) * t117 + t126 * t95 - t165 * t57 + (qJD(2) * t242 - t31) * t202 + (pkin(6) * t68 + qJD(2) * t82 + t109 * t198 + t148 * t280) * t199, -g(1) * t114 - g(2) * t116 - t126 * t96 + t165 * t56 + (qJD(2) * t241 + t30) * t202 + (-pkin(6) * t67 - qJD(2) * t83 + t109 * t201 - t148 * t282) * t199, -t132 * t56 - t134 * t57 + t67 * t95 - t68 * t96 + t172 + t239 * t283 + (-t339 - t198 * t30 - t201 * t31 + (t198 * t82 - t201 * t83) * qJD(3)) * t199, t30 * t96 + t83 * t56 + t31 * t95 + t82 * t57 - g(1) * t190 - g(2) * (t203 * t250 + t290) - t146 * t342 + (t109 * t199 + t148 * t283) * pkin(6), t8, -t213, t6, t244, -t209, t85, t104 * t51 + t126 * t46 - t13 * t165 + t137 * t35 - t202 * t3 + t25 * t284 - t60 * t89 + t73 * t92 + t248, t105 * t51 - t126 * t47 + t137 * t36 + t14 * t165 + t202 * t4 + t225 * t92 - t26 * t284 - t61 * t89 - t253, -t104 * t4 - t105 * t3 - t13 * t225 - t14 * t73 + t25 * t61 + t26 * t60 - t35 * t47 - t36 * t46 + t255, -g(1) * t218 - g(2) * t221 + t25 * t13 + t51 * t137 + t26 * t14 + t3 * t46 + t4 * t47 + t89 * t92, t8, t6, t213, t85, t209, t244, t10 * t165 + t104 * t5 - t126 * t44 + t19 * t73 + t2 * t202 - t22 * t284 - t29 * t60 + t35 * t53 + t248, -t1 * t104 + t10 * t225 + t105 * t2 - t22 * t61 + t23 * t60 - t35 * t43 + t36 * t44 - t73 * t9 + t255, -t1 * t202 - t105 * t5 + t126 * t43 - t165 * t9 - t19 * t225 + t23 * t284 + t29 * t61 - t36 * t53 + t253, t1 * t43 + t23 * t9 + t5 * t53 + t29 * t19 + t2 * t44 + t22 * t10 - g(1) * (-pkin(4) * t100 - qJ(5) * t99 + t218) - g(2) * (pkin(4) * t102 + qJ(5) * t101 + t221); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t272, t289 * t205, t276, t272, t186, qJDD(2), t199 * t224 - t178 - t192, t338 + (t224 - t315) * t202, 0, 0, -t165 * t309 - t321, (-t67 + t313) * t201 + (-t68 + t310) * t198, (-t134 * t199 + t165 * t297) * qJD(1) + t227, -t165 * t312 - t320, (t132 * t199 - t165 * t305) * qJD(1) + t226, t265, -pkin(2) * t68 + t165 * t90 + t238 * t198 + (-t199 * t82 - t202 * t242) * qJD(1) + t359 * t201, pkin(2) * t67 - t165 * t91 + t238 * t201 + (t199 * t83 - t202 * t241) * qJD(1) - t359 * t198, t132 * t91 + t134 * t90 + ((qJD(3) * t134 - t68) * pkin(7) + t354) * t201 + (-t31 + t323 + (qJD(3) * t132 - t67) * pkin(7)) * t198 - t350, -t148 * t180 - t82 * t90 - t83 * t91 + t210 * pkin(2) + (qJD(3) * t239 - t31 * t198 + t30 * t201 - t350) * pkin(7), t7, -t347, t21, t228, t351, t265, t127 * t51 - t331 * t165 - t177 * t35 - t25 * t287 + t254 * t73 + t318 * t89 + t219, t128 * t51 + t328 * t165 - t177 * t36 + t225 * t254 + t26 * t287 + t317 * t89 + t358, -t127 * t4 - t128 * t3 - t225 * t331 - t317 * t25 - t318 * t26 + t39 * t73 + t208, t4 * t88 - t3 * t87 - t51 * t177 - g(3) * (-t197 * t199 + t151) + t254 * t89 + t328 * t26 + t331 * t25 + t247 * (t177 * t199 + t197 * t202), t7, t21, t347, t265, -t351, t228, t127 * t5 + t330 * t165 + t22 * t287 + t318 * t29 + t332 * t73 + t35 * t70 + t219, -t1 * t127 + t128 * t2 + t317 * t22 + t225 * t330 - t318 * t23 + t33 * t73 + t208, -t128 * t5 - t329 * t165 - t225 * t332 - t23 * t287 - t317 * t29 - t36 * t70 - t358, -g(3) * t151 + t1 * t88 + t2 * t87 + t5 * t70 + t332 * t29 + t329 * t23 + t330 * t22 + (-g(3) * t243 + t197 * t247) * t202 + (g(3) * t197 + t247 * (t177 + t243)) * t199; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t311, -t132 ^ 2 + t134 ^ 2, -t67 - t313, -t311, -t310 - t68, t126, -t149 * t280 - g(1) * t116 + g(2) * t114 - t134 * t148 - t323 + t79 + (-qJD(3) * t123 - t108 + t338) * t198, g(1) * t117 - g(2) * t115 + g(3) * t302 + t132 * t148 - t354, 0, 0, t357, t245, t20, -t357, -t235, t126, -t89 * t225 + (t126 * t316 - t134 * t73) * pkin(3) + t211, -t165 * t28 + t73 * t89 + (-t126 * t196 - t134 * t225) * pkin(3) + t212, t26 * t225 - t337 + (-t196 * t35 - t316 * t36) * pkin(3) + (-t25 + t28) * t73, -g(1) * t169 + t25 * t27 - t26 * t28 + (g(2) * t296 - t89 * t134 + t4 * t196 + t350 * t198 + t3 * t316) * pkin(3), t357, t20, -t245, t126, t235, -t357, -t32 * t73 + (pkin(4) - t175) * t126 + t211 + t348, -t173 * t35 + t175 * t36 + t225 * t23 - t337 + (t22 - t293) * t73, t126 * t173 - t29 * t73 + t32 * t225 + t122 + (-0.2e1 * qJD(5) + t28) * t165 - t212, t1 * t173 + t2 * t175 - t29 * t32 - t22 * t27 - g(1) * (-pkin(3) * t273 - pkin(4) * t101 + qJ(5) * t102 + t169) - g(2) * (-pkin(3) * t114 - pkin(4) * t99 + qJ(5) * t100) - g(3) * (-t167 + (-pkin(4) * t182 + qJ(5) * t183) * t199) + t293 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236, -t217, t246, t225 * t25 + t26 * t73 + t207, 0, 0, 0, 0, 0, 0, t236, t246, t217, -t22 * t225 + t23 * t73 + t207 + t352; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126 + t357, t20, -t165 ^ 2 - t336, t165 * t23 - t222 - t344 - t348;];
tau_reg = t11;
