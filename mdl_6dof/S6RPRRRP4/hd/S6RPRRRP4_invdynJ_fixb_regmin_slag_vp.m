% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRP4
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
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:08:49
% EndTime: 2019-03-09 06:08:59
% DurationCPUTime: 4.52s
% Computational Cost: add. (8740->430), mult. (21295->532), div. (0->0), fcn. (16944->14), ass. (0->242)
t186 = pkin(10) + qJ(3);
t179 = qJ(4) + t186;
t170 = sin(t179);
t193 = sin(qJ(1));
t197 = cos(qJ(1));
t224 = g(1) * t197 + g(2) * t193;
t359 = t224 * t170;
t190 = sin(qJ(5));
t194 = cos(qJ(5));
t262 = qJD(3) + qJD(4);
t191 = sin(qJ(4));
t195 = cos(qJ(4));
t188 = cos(pkin(10));
t196 = cos(qJ(3));
t276 = t196 * t188;
t187 = sin(pkin(10));
t192 = sin(qJ(3));
t281 = t187 * t192;
t139 = -t276 + t281;
t129 = t139 * qJD(1);
t307 = pkin(7) + qJ(2);
t154 = t307 * t187;
t141 = qJD(1) * t154;
t155 = t307 * t188;
t142 = qJD(1) * t155;
t216 = t141 * t192 - t142 * t196;
t88 = -pkin(8) * t129 - t216;
t84 = t195 * t88;
t140 = t187 * t196 + t188 * t192;
t130 = t140 * qJD(1);
t334 = -t196 * t141 - t142 * t192;
t87 = -pkin(8) * t130 + t334;
t85 = qJD(3) * pkin(3) + t87;
t55 = t191 * t85 + t84;
t49 = pkin(9) * t262 + t55;
t249 = pkin(2) * t188 + pkin(1);
t148 = -qJD(1) * t249 + qJD(2);
t112 = pkin(3) * t129 + t148;
t217 = -t129 * t191 + t195 * t130;
t99 = t195 * t129 + t130 * t191;
t56 = pkin(4) * t99 - pkin(9) * t217 + t112;
t22 = -t190 * t49 + t194 * t56;
t91 = t190 * t262 + t194 * t217;
t15 = -qJ(6) * t91 + t22;
t348 = qJD(5) + t99;
t11 = pkin(5) * t348 + t15;
t164 = g(3) * t170;
t171 = cos(t179);
t251 = t224 * t171 + t164;
t248 = qJD(1) * t281;
t263 = t188 * qJDD(1);
t264 = t187 * qJDD(1);
t271 = qJD(3) * t196;
t250 = t188 * qJD(1) * t271 + t192 * t263 + t196 * t264;
t106 = -qJD(3) * t248 + t250;
t132 = t140 * qJD(3);
t222 = t192 * t264 - t196 * t263;
t107 = qJD(1) * t132 + t222;
t269 = qJD(4) * t195;
t270 = qJD(4) * t191;
t208 = -t195 * t106 + t191 * t107 + t129 * t269 + t130 * t270;
t236 = t191 * t106 + t195 * t107;
t330 = qJD(4) * t217;
t47 = t236 + t330;
t147 = -qJDD(1) * t249 + qJDD(2);
t93 = t107 * pkin(3) + t147;
t18 = t47 * pkin(4) + pkin(9) * t208 + t93;
t267 = qJD(5) * t194;
t268 = qJD(5) * t190;
t261 = qJDD(3) + qJDD(4);
t265 = qJD(1) * qJD(2);
t322 = t307 * qJDD(1) + t265;
t116 = t322 * t187;
t117 = t322 * t188;
t235 = -t196 * t116 - t192 * t117;
t43 = qJDD(3) * pkin(3) - pkin(8) * t106 + qJD(3) * t216 + t235;
t218 = -t192 * t116 + t196 * t117;
t46 = -pkin(8) * t107 + t334 * qJD(3) + t218;
t323 = -(qJD(4) * t85 + t46) * t195 - t191 * t43 + t88 * t270;
t9 = pkin(9) * t261 - t323;
t211 = t190 * t18 + t194 * t9 + t56 * t267 - t49 * t268;
t32 = qJD(5) * t91 - t190 * t208 - t194 * t261;
t229 = t194 * t262;
t89 = t190 * t217 - t229;
t3 = -qJ(6) * t32 - qJD(6) * t89 + t211;
t14 = t194 * t18;
t23 = t190 * t56 + t194 * t49;
t31 = -qJD(5) * t229 - t190 * t261 + t194 * t208 + t217 * t268;
t45 = qJDD(5) + t47;
t1 = pkin(5) * t45 + qJ(6) * t31 - qJD(5) * t23 - qJD(6) * t91 - t190 * t9 + t14;
t16 = -qJ(6) * t89 + t23;
t327 = t16 * t348 + t1;
t357 = t194 * t348;
t358 = -t11 * t357 - t327 * t190 + t3 * t194 - t251;
t336 = t262 * t99;
t355 = -t208 + t336;
t29 = t31 * t190;
t351 = t194 * t99;
t354 = -t29 + (t267 + t351) * t91;
t38 = t190 * t45;
t304 = t267 * t348 + t38;
t308 = t91 * t217;
t353 = t348 * t351 + t304 - t308;
t83 = t191 * t88;
t54 = t195 * t85 - t83;
t48 = -pkin(4) * t262 - t54;
t352 = t48 * t99;
t350 = t217 * t99;
t256 = pkin(3) * t269;
t58 = t195 * t87 - t83;
t349 = t58 - t256;
t178 = cos(t186);
t174 = pkin(5) * t194 + pkin(4);
t189 = -qJ(6) - pkin(9);
t233 = -t170 * t189 + t171 * t174;
t346 = pkin(3) * t178 + t233;
t287 = qJDD(1) * pkin(1);
t333 = g(1) * t193 - g(2) * t197;
t213 = -qJDD(2) + t287 + t333;
t344 = t217 ^ 2 - t99 ^ 2;
t72 = pkin(4) * t217 + pkin(9) * t99;
t343 = t112 * t99 + t251 + t323;
t181 = t194 * qJ(6);
t342 = -pkin(5) * t217 - t99 * t181;
t309 = t89 * t217;
t338 = t348 * t217;
t39 = t194 * t45;
t337 = -t268 * t348 + t39;
t335 = t333 * t170;
t273 = -t192 * t154 + t196 * t155;
t275 = t217 * qJD(3);
t332 = t275 - t236;
t277 = t194 * t197;
t280 = t190 * t193;
t121 = t171 * t280 + t277;
t278 = t193 * t194;
t279 = t190 * t197;
t123 = -t171 * t279 + t278;
t331 = -g(1) * t123 + g(2) * t121;
t329 = qJ(2) * qJDD(1);
t312 = g(3) * t171;
t328 = -t312 + t359;
t40 = t48 * t268;
t326 = t194 * t359 - t22 * t217 + t40;
t240 = t191 * t46 - t195 * t43 + t88 * t269 + t85 * t270;
t10 = -pkin(4) * t261 + t240;
t311 = g(3) * t190;
t325 = t10 * t190 + t171 * t311 + t23 * t217 + t48 * t267;
t324 = -t112 * t217 - t240 + t328;
t321 = t91 ^ 2;
t320 = pkin(3) * t132;
t319 = pkin(3) * t195;
t318 = pkin(5) * t190;
t108 = t195 * t139 + t140 * t191;
t131 = t139 * qJD(3);
t73 = -qJD(4) * t108 - t131 * t195 - t132 * t191;
t310 = t48 * t73;
t306 = -t15 + t11;
t305 = -t190 * t32 - t89 * t267;
t303 = t190 * t72 + t194 * t54;
t61 = pkin(3) * t130 + t72;
t302 = t190 * t61 + t194 * t58;
t234 = -t196 * t154 - t155 * t192;
t95 = -pkin(8) * t140 + t234;
t96 = -pkin(8) * t139 + t273;
t67 = t191 * t95 + t195 * t96;
t63 = t194 * t67;
t109 = -t139 * t191 + t140 * t195;
t114 = pkin(3) * t139 - t249;
t68 = pkin(4) * t108 - pkin(9) * t109 + t114;
t301 = t190 * t68 + t63;
t299 = t190 * t91;
t298 = t190 * t348;
t30 = t194 * t31;
t297 = t194 * t89;
t296 = t194 * t91;
t180 = t194 * qJD(6);
t173 = t191 * pkin(3) + pkin(9);
t274 = -qJ(6) - t173;
t232 = qJD(5) * t274;
t286 = t99 * t190;
t293 = -qJ(6) * t286 + t190 * t232 + t194 * t256 + t180 - t302;
t60 = t194 * t61;
t292 = t194 * t232 - t60 + (-qJD(6) + t349) * t190 + t342;
t241 = qJD(5) * t189;
t291 = t180 - t303 + (-qJ(6) * t99 + t241) * t190;
t70 = t194 * t72;
t290 = t194 * t241 - t70 + (-qJD(6) + t54) * t190 + t342;
t288 = qJD(5) * t348;
t284 = t109 * t190;
t272 = t187 ^ 2 + t188 ^ 2;
t66 = t191 * t96 - t195 * t95;
t204 = -t154 * t271 + qJD(2) * t276 + (-qJD(2) * t187 - qJD(3) * t155) * t192;
t75 = -pkin(8) * t132 + t204;
t200 = -t140 * qJD(2) - t273 * qJD(3);
t76 = pkin(8) * t131 + t200;
t26 = -qJD(4) * t66 + t191 * t76 + t195 * t75;
t74 = qJD(4) * t109 - t131 * t191 + t195 * t132;
t35 = pkin(4) * t74 - pkin(9) * t73 + t320;
t259 = t190 * t35 + t194 * t26 + t68 * t267;
t258 = -t30 + t305;
t257 = pkin(9) * t288;
t254 = t173 * t288;
t253 = t91 * t268;
t247 = t109 * t267;
t246 = -t10 - t312;
t245 = pkin(8) + t307 + t318;
t243 = -qJD(5) * t56 - t9;
t237 = t272 * qJD(1) ^ 2;
t228 = 0.2e1 * t272;
t57 = t191 * t87 + t84;
t227 = pkin(3) * t270 - t57;
t226 = -pkin(9) * t45 + t352;
t225 = -t49 * t267 + t14;
t221 = -t173 * t45 + t352;
t220 = t297 + t299;
t219 = -qJ(6) * t73 - qJD(6) * t109;
t215 = t170 * t174 + t171 * t189;
t214 = -t286 * t348 + t337;
t212 = -t253 - t30;
t209 = -t249 - t346;
t205 = t213 + t287;
t6 = t32 * pkin(5) + qJDD(6) + t10;
t202 = t228 * t265 - t224;
t27 = qJD(4) * t67 + t191 * t75 - t195 * t76;
t177 = sin(t186);
t175 = -pkin(4) - t319;
t157 = pkin(9) * t194 + t181;
t156 = t189 * t190;
t135 = t173 * t194 + t181;
t134 = t274 * t190;
t124 = t171 * t277 + t280;
t122 = -t171 * t278 + t279;
t86 = t89 ^ 2;
t65 = t194 * t68;
t36 = t89 * pkin(5) + qJD(6) + t48;
t34 = t194 * t35;
t25 = -qJ(6) * t284 + t301;
t19 = pkin(5) * t108 - t109 * t181 - t190 * t67 + t65;
t5 = -qJ(6) * t247 + (-qJD(5) * t67 + t219) * t190 + t259;
t4 = pkin(5) * t74 - t190 * t26 + t34 + t219 * t194 + (-t63 + (qJ(6) * t109 - t68) * t190) * qJD(5);
t2 = [qJDD(1), t333, t224, t205 * t188, -t205 * t187, t228 * t329 + t202, t213 * pkin(1) + (t272 * t329 + t202) * qJ(2), t106 * t140 - t130 * t131, -t106 * t139 - t107 * t140 + t129 * t131 - t130 * t132, -qJD(3) * t131 + qJDD(3) * t140, -qJD(3) * t132 - qJDD(3) * t139, 0, t200 * qJD(3) + t234 * qJDD(3) - t107 * t249 + t148 * t132 + t147 * t139 + t178 * t333, -t204 * qJD(3) - t273 * qJDD(3) - t106 * t249 - t148 * t131 + t147 * t140 - t177 * t333, -t109 * t208 + t217 * t73, t108 * t208 - t109 * t47 - t217 * t74 - t73 * t99, t109 * t261 + t262 * t73, -t108 * t261 - t262 * t74, 0, t93 * t108 + t112 * t74 + t114 * t47 + t171 * t333 - t261 * t66 - t262 * t27 + t99 * t320, t93 * t109 + t112 * t73 - t114 * t208 + t217 * t320 - t26 * t262 - t261 * t67 - t335, t109 * t212 + t73 * t296, -t220 * t73 + (t29 - t194 * t32 + (t190 * t89 - t296) * qJD(5)) * t109, -t108 * t31 + t337 * t109 + t357 * t73 + t74 * t91, -t108 * t32 - t109 * t304 - t298 * t73 - t74 * t89, t108 * t45 + t348 * t74 (-t267 * t67 + t34) * t348 + t65 * t45 + t225 * t108 + t22 * t74 + t27 * t89 + t66 * t32 + t48 * t247 - g(1) * t122 - g(2) * t124 + ((-qJD(5) * t68 - t26) * t348 - t67 * t45 + t243 * t108 + t10 * t109 + t310) * t190 -(-t268 * t67 + t259) * t348 - t301 * t45 - t211 * t108 - t23 * t74 + t27 * t91 - t66 * t31 + t194 * t310 - g(1) * t121 - g(2) * t123 + (t10 * t194 - t40) * t109, t19 * t31 - t25 * t32 - t4 * t91 - t5 * t89 + (-t11 * t194 - t16 * t190) * t73 + t335 + (-t1 * t194 - t190 * t3 + (t11 * t190 - t16 * t194) * qJD(5)) * t109, t3 * t25 + t16 * t5 + t1 * t19 + t11 * t4 + t6 * (pkin(5) * t284 + t66) + t36 * ((t190 * t73 + t247) * pkin(5) + t27) + (-g(1) * t245 + g(2) * t209) * t197 + (-g(1) * t209 - g(2) * t245) * t193; 0, 0, 0, -t263, t264, -t237, -qJ(2) * t237 - t213, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t130 + t222 (-t129 - t248) * qJD(3) + t250, 0, 0, 0, 0, 0, t236 + t275 + 0.2e1 * t330, -t208 - t336, 0, 0, 0, 0, 0, t214 - t309, -t348 * t357 - t308 - t38 -(t297 - t299) * t99 - t212 + t305, -t36 * t217 + t327 * t194 + (-t11 * t348 + t3) * t190 - t333; 0, 0, 0, 0, 0, 0, 0, t129 * t130, -t129 ^ 2 + t130 ^ 2 (t129 - t248) * qJD(3) + t250, -t222, qJDD(3), -g(3) * t178 - t148 * t130 + t177 * t224 + t235, g(3) * t177 + t148 * t129 + t224 * t178 - t218, t350, t344, t355, t332, t261, t57 * t262 + (-t130 * t99 + t195 * t261 - t262 * t270) * pkin(3) + t324, t58 * t262 + (-t130 * t217 - t191 * t261 - t262 * t269) * pkin(3) + t343, t354, -t220 * t99 - t253 + t258, t353, t214 + t309, -t338, t175 * t32 - t60 * t348 + t227 * t89 + (t246 - t254) * t194 + (t348 * t349 + t221) * t190 + t326, -t175 * t31 + t302 * t348 + t227 * t91 + (-t256 * t348 + t221) * t194 + (-t359 + t254) * t190 + t325, t134 * t31 - t135 * t32 - t292 * t91 - t293 * t89 + t358, t3 * t135 + t1 * t134 + t6 * (-t174 - t319) - g(3) * t346 + (-t84 + (pkin(3) * qJD(4) - t87) * t191 + t348 * t318) * t36 + t293 * t16 + t292 * t11 + t224 * (pkin(3) * t177 + t215); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t350, t344, t355, t332, t261, t262 * t55 + t324, t262 * t54 + t343, t354, -t298 * t91 - t351 * t89 + t258, t353, -t298 * t348 + t309 + t39, -t338, -pkin(4) * t32 - t55 * t89 - t70 * t348 + (t348 * t54 + t226) * t190 + (t246 - t257) * t194 + t326, pkin(4) * t31 + t303 * t348 - t55 * t91 + t226 * t194 + (-t359 + t257) * t190 + t325, t156 * t31 - t157 * t32 - t290 * t91 - t291 * t89 + t358, t3 * t157 + t1 * t156 - t6 * t174 - g(3) * t233 + (pkin(5) * t298 - t55) * t36 + t291 * t16 + t290 * t11 + t224 * t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91 * t89, -t86 + t321, t348 * t89 - t31, t348 * t91 - t32, t45, t23 * t348 - t48 * t91 + (t243 + t164) * t190 + t225 + t331, g(1) * t124 - g(2) * t122 + t164 * t194 + t22 * t348 + t48 * t89 - t211, pkin(5) * t31 - t306 * t89, t306 * t16 + (t170 * t311 - t36 * t91 + t1 + t331) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86 - t321, t11 * t91 + t16 * t89 - t328 + t6;];
tau_reg  = t2;
