% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:16:29
% EndTime: 2019-03-08 21:16:41
% DurationCPUTime: 4.46s
% Computational Cost: add. (3353->520), mult. (8063->705), div. (0->0), fcn. (6473->12), ass. (0->250)
t197 = sin(qJ(3));
t350 = qJ(4) * t197 + pkin(2);
t200 = cos(qJ(3));
t155 = -pkin(3) * t200 - t350;
t245 = pkin(3) * t197 - qJ(4) * t200;
t124 = qJD(3) * t245 - qJD(4) * t197;
t195 = cos(pkin(11));
t193 = sin(pkin(11));
t194 = sin(pkin(6));
t201 = cos(qJ(2));
t308 = t200 * t201;
t283 = t194 * t308;
t257 = t193 * t283;
t198 = sin(qJ(2));
t302 = qJD(1) * t194;
t276 = t198 * t302;
t349 = qJD(1) * t257 + (t124 - t276) * t195;
t150 = qJD(2) * pkin(8) + t276;
t322 = cos(pkin(6));
t261 = qJD(1) * t322;
t217 = t197 * t150 - t200 * t261;
t298 = qJD(3) * t193;
t300 = qJD(2) * t197;
t144 = t195 * t300 + t298;
t140 = t144 ^ 2;
t289 = t195 * qJD(3);
t142 = t193 * t300 - t289;
t348 = -t142 ^ 2 - t140;
t189 = t200 * qJDD(2);
t287 = qJD(2) * qJD(3);
t268 = t197 * t287;
t347 = -t268 + t189;
t299 = qJD(2) * t200;
t181 = qJD(6) + t299;
t346 = t181 - qJD(6);
t345 = qJDD(3) * pkin(3) - qJDD(4);
t297 = qJD(3) * t197;
t101 = (t193 * t198 + t195 * t308) * t194;
t93 = qJD(1) * t101;
t344 = qJ(5) * t297 - qJD(5) * t200 - t93;
t312 = t194 * t198;
t133 = t197 * t312 - t200 * t322;
t321 = cos(pkin(10));
t247 = t322 * t321;
t320 = sin(pkin(10));
t127 = t198 * t247 + t201 * t320;
t263 = t194 * t321;
t74 = t127 * t197 + t200 * t263;
t246 = t322 * t320;
t129 = -t198 * t246 + t201 * t321;
t262 = t194 * t320;
t76 = t129 * t197 - t200 * t262;
t224 = g(1) * t76 + g(2) * t74 + g(3) * t133;
t196 = sin(qJ(6));
t199 = cos(qJ(6));
t63 = t142 * t196 + t144 * t199;
t187 = t195 * qJDD(3);
t267 = t200 * t287;
t286 = t197 * qJDD(2);
t221 = t267 + t286;
t98 = t193 * t221 - t187;
t99 = qJDD(3) * t193 + t195 * t221;
t13 = qJD(6) * t63 + t196 * t99 - t199 * t98;
t202 = qJD(3) ^ 2;
t288 = qJD(1) * qJD(2);
t269 = t198 * t288;
t311 = t194 * t201;
t242 = -qJDD(1) * t311 + t194 * t269;
t126 = t198 * t320 - t201 * t247;
t128 = t198 * t321 + t201 * t246;
t251 = g(1) * t128 + g(2) * t126;
t343 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t202 + t194 * (-g(3) * t201 + t269) - t242 + t251;
t107 = qJDD(2) * pkin(8) + (qJDD(1) * t198 + t201 * t288) * t194;
t173 = t197 * t261;
t259 = qJDD(1) * t322;
t296 = qJD(3) * t200;
t258 = -qJD(3) * t173 - t197 * t107 - t150 * t296 + t200 * t259;
t216 = t258 + t345;
t294 = qJD(5) * t144;
t326 = qJ(5) * t99;
t206 = t216 + t294 + t326;
t340 = pkin(4) + pkin(5);
t6 = -t340 * t98 + t206;
t342 = t6 + t224;
t301 = qJD(2) * t194;
t273 = t198 * t301;
t272 = t201 * t301;
t78 = -qJD(3) * t133 + t200 * t272;
t44 = t193 * t273 + t195 * t78;
t134 = t197 * t322 + t200 * t312;
t73 = t134 * t195 - t193 * t311;
t79 = qJD(3) * t134 + t197 * t272;
t341 = qJD(2) * (-t200 * t44 + t297 * t73) - t133 * t99 - t144 * t79 - t189 * t73;
t339 = pkin(4) * t98;
t335 = pkin(9) * t197;
t333 = -pkin(9) + qJ(4);
t248 = t197 * t259;
t21 = t248 + qJDD(3) * qJ(4) + t200 * t107 + (qJD(4) - t217) * qJD(3);
t36 = qJD(2) * t124 + qJDD(2) * t155 + t242;
t11 = t193 * t36 + t195 * t21;
t279 = -pkin(8) * t193 - pkin(4);
t309 = t195 * t200;
t285 = pkin(9) * t309;
t332 = (-t285 + (-pkin(5) + t279) * t197) * qJD(3) - t349;
t110 = t193 * t124;
t310 = t195 * t197;
t313 = t193 * t200;
t331 = -t110 - (-pkin(8) * t310 + pkin(9) * t313) * qJD(3) - t344;
t284 = pkin(8) * t297;
t71 = -t195 * t284 + t110;
t330 = t71 + t344;
t329 = t279 * t297 - t349;
t96 = t200 * t150 + t173;
t87 = qJD(3) * qJ(4) + t96;
t275 = t201 * t302;
t97 = qJD(2) * t155 - t275;
t30 = t193 * t97 + t195 * t87;
t328 = t193 * t284 + t349;
t327 = t71 - t93;
t325 = qJD(2) * pkin(2);
t61 = -t199 * t142 + t144 * t196;
t324 = t181 * t61;
t323 = t181 * t63;
t149 = t245 * qJD(2);
t46 = t193 * t149 - t195 * t217;
t318 = qJ(5) * t195;
t314 = t193 * t199;
t307 = qJDD(1) - g(3);
t147 = t193 * t196 + t195 * t199;
t222 = t147 * t200;
t306 = -qJD(2) * t222 - t147 * qJD(6);
t274 = t193 * t299;
t282 = t196 * t309;
t290 = qJD(6) * t199;
t291 = qJD(6) * t196;
t305 = -qJD(2) * t282 + t193 * t290 - t195 * t291 + t199 * t274;
t270 = qJ(4) * t189;
t304 = qJD(4) * t274 + t193 * t270;
t103 = pkin(8) * t309 + t193 * t155;
t191 = t197 ^ 2;
t192 = t200 ^ 2;
t303 = t191 - t192;
t295 = qJD(4) * t195;
t293 = qJD(5) * t193;
t39 = qJ(5) * t300 + t46;
t281 = t155 * t126;
t280 = t155 * t128;
t10 = -t193 * t21 + t195 * t36;
t234 = pkin(4) * t189 + qJDD(5) - t10;
t8 = -pkin(4) * t268 + t234;
t2 = t347 * pkin(5) - pkin(9) * t99 + t8;
t7 = qJ(5) * t268 + (-qJ(5) * qJDD(2) - qJD(2) * qJD(5)) * t200 + t11;
t5 = pkin(9) * t98 + t7;
t278 = -t196 * t5 + t199 * t2;
t277 = qJ(4) * t297;
t271 = qJD(5) * t310;
t266 = t201 * t287;
t265 = qJ(5) * t193 + pkin(3);
t29 = -t193 * t87 + t195 * t97;
t45 = t149 * t195 + t193 * t217;
t232 = -t193 * t340 + t318;
t260 = -t232 * t299 + t293 + t96;
t177 = pkin(8) * t313;
t102 = t155 * t195 - t177;
t256 = pkin(3) * t283 + pkin(8) * t312 + t311 * t350;
t255 = t142 * t275;
t254 = t144 * t275;
t253 = t197 * t275;
t252 = t195 * t270;
t250 = g(1) * t129 + g(2) * t127;
t249 = t196 * t2 + t199 * t5;
t86 = -qJ(5) * t200 + t103;
t244 = pkin(4) * t193 - t318;
t22 = pkin(4) * t299 + qJD(5) - t29;
t14 = pkin(5) * t299 - pkin(9) * t144 + t22;
t23 = -qJ(5) * t299 + t30;
t15 = pkin(9) * t142 + t23;
t3 = t14 * t199 - t15 * t196;
t4 = t14 * t196 + t15 * t199;
t190 = t200 * pkin(4);
t59 = pkin(5) * t200 + t177 + t190 + (-t155 - t335) * t195;
t64 = t193 * t335 + t86;
t241 = t196 * t64 - t199 * t59;
t240 = t196 * t59 + t199 * t64;
t72 = t134 * t193 + t195 * t311;
t19 = -t196 * t73 + t199 * t72;
t20 = t196 * t72 + t199 * t73;
t239 = qJ(4) * t99 + qJD(4) * t144;
t148 = -t195 * t196 + t314;
t238 = t181 ^ 2;
t203 = qJD(2) ^ 2;
t237 = qJDD(2) * t201 - t198 * t203;
t236 = pkin(8) + t244;
t235 = qJD(3) * pkin(3) - qJD(4) - t217;
t158 = t333 * t193;
t230 = pkin(9) * t274 - qJD(6) * t158 - t295 + t39;
t159 = t333 * t195;
t229 = qJD(4) * t193 - qJD(6) * t159 - (-t197 * t340 - t285) * qJD(2) + t45;
t228 = -t142 * t290 + t144 * t291 - t196 * t98 - t199 * t99;
t43 = t193 * t78 - t195 * t273;
t227 = -t44 * t142 + t144 * t43 + t72 * t99 - t73 * t98;
t100 = -t195 * t312 + t257;
t47 = -t126 * t313 - t127 * t195;
t49 = -t128 * t313 - t129 * t195;
t226 = g(1) * t49 + g(2) * t47 + g(3) * t100;
t48 = -t126 * t309 + t127 * t193;
t50 = -t128 * t309 + t129 * t193;
t225 = -g(1) * t50 - g(2) * t48 - g(3) * t101;
t75 = t127 * t200 - t197 * t263;
t77 = t129 * t200 + t197 * t262;
t223 = g(1) * t77 + g(2) * t75 + g(3) * t134;
t117 = t147 * t197;
t220 = -pkin(8) + t232;
t9 = -t206 + t339;
t219 = t224 - t9;
t218 = t224 + t216;
t215 = -g(3) * t311 + t251;
t214 = qJ(5) * t144 + t235;
t213 = -t195 * qJ(4) * t98 - t142 * t295 - t223;
t210 = t224 + t258;
t208 = t193 * t286 - t187 + (-t144 + t298) * t299;
t151 = -t275 - t325;
t207 = -pkin(8) * qJDD(3) + (t151 + t275 - t325) * qJD(3);
t205 = t72 * t189 + t79 * t142 + t133 * t98 + (t200 * t43 - t297 * t72) * qJD(2);
t204 = -t210 - t345;
t152 = -pkin(4) * t195 - t265;
t146 = -qJDD(6) - t347;
t136 = t195 * t340 + t265;
t123 = t133 * pkin(3);
t119 = t142 * t299;
t116 = t196 * t310 - t197 * t314;
t115 = t236 * t197;
t91 = -t102 + t190;
t82 = t220 * t197;
t67 = t76 * pkin(3);
t66 = t74 * pkin(3);
t65 = t236 * t296 - t271;
t57 = t220 * t296 + t271;
t56 = t244 * t299 + t96;
t54 = t119 + t99;
t53 = qJD(3) * t282 + qJD(6) * t117 - t296 * t314;
t52 = qJD(6) * t148 * t197 + qJD(3) * t222;
t42 = -pkin(4) * t300 - t45;
t35 = t128 * t193 + t195 * t77;
t34 = -t128 * t195 + t193 * t77;
t33 = t126 * t193 + t195 * t75;
t32 = -t126 * t195 + t193 * t75;
t27 = pkin(4) * t142 - t214;
t16 = -t142 * t340 + t214;
t1 = [t307, 0, t237 * t194 (-qJDD(2) * t198 - t201 * t203) * t194, 0, 0, 0, 0, 0, -qJD(3) * t79 - qJDD(3) * t133 + (-t197 * t266 + t200 * t237) * t194, -qJD(3) * t78 - qJDD(3) * t134 + (-t197 * t237 - t200 * t266) * t194, t205, -t341, t227, -t10 * t72 + t11 * t73 - t133 * t216 - t235 * t79 - t29 * t43 + t30 * t44 - g(3), t205, t227, t341, t133 * t9 + t22 * t43 + t23 * t44 + t27 * t79 + t7 * t73 + t72 * t8 - g(3), 0, 0, 0, 0, 0 (-qJD(6) * t20 - t196 * t44 + t199 * t43) * t181 - t19 * t146 - t79 * t61 - t133 * t13 -(qJD(6) * t19 + t196 * t43 + t199 * t44) * t181 + t20 * t146 - t79 * t63 + t133 * t228; 0, qJDD(2), t307 * t311 + t251, -t307 * t312 + t250, qJDD(2) * t191 + 0.2e1 * t197 * t267, 0.2e1 * t189 * t197 - 0.2e1 * t287 * t303, qJDD(3) * t197 + t200 * t202, qJDD(3) * t200 - t197 * t202, 0, t207 * t197 + t343 * t200, -t343 * t197 + t207 * t200 (-t255 + pkin(8) * t98 - t193 * t216 + (qJD(2) * t102 + t29) * qJD(3)) * t197 + (-qJDD(2) * t102 - t10 + (pkin(8) * t142 - t193 * t235) * qJD(3) - t328 * qJD(2)) * t200 + t225 (-t254 + pkin(8) * t99 - t195 * t216 + (-qJD(2) * t103 - t30) * qJD(3)) * t197 + (qJDD(2) * t103 + t11 + (pkin(8) * t144 - t195 * t235) * qJD(3) + t327 * qJD(2)) * t200 + t226, -t102 * t99 - t103 * t98 - t328 * t144 - t327 * t142 + (-t193 * t30 - t195 * t29) * t296 + (-t10 * t195 - t11 * t193 + t215) * t197, t11 * t103 + t10 * t102 + t235 * t253 - g(1) * t280 - g(2) * t281 - g(3) * t256 + t327 * t30 + t328 * t29 + (-t197 * t216 - t235 * t296 - t250) * pkin(8), t115 * t98 + t142 * t65 + (-t255 + t193 * t9 + (-qJD(2) * t91 - t22) * qJD(3)) * t197 + (qJD(2) * t329 + qJDD(2) * t91 + t27 * t298 + t8) * t200 + t225, -t86 * t98 + t91 * t99 + t329 * t144 - t330 * t142 + (-t193 * t23 + t195 * t22) * t296 + (-t193 * t7 + t195 * t8 + t215) * t197, -t115 * t99 - t144 * t65 + (t254 - t195 * t9 + (qJD(2) * t86 + t23) * qJD(3)) * t197 + (-qJD(2) * t330 - qJDD(2) * t86 - t27 * t289 - t7) * t200 - t226, t7 * t86 + t9 * t115 + t8 * t91 - g(1) * (pkin(4) * t50 + pkin(8) * t129 + qJ(5) * t49 + t280) - g(2) * (pkin(4) * t48 + pkin(8) * t127 + qJ(5) * t47 + t281) - g(3) * (pkin(4) * t101 + qJ(5) * t100 + t256) + (t65 - t253) * t27 + t330 * t23 + t329 * t22, -t117 * t228 + t52 * t63, t116 * t228 - t117 * t13 - t52 * t61 - t53 * t63, -t117 * t146 + t181 * t52 - t200 * t228 - t297 * t63, t116 * t146 - t13 * t200 - t181 * t53 + t297 * t61, -t146 * t200 - t181 * t297, t241 * t146 + t278 * t200 + t57 * t61 + t82 * t13 + t6 * t116 + t16 * t53 - g(1) * (t196 * t49 + t199 * t50) - g(2) * (t196 * t47 + t199 * t48) - g(3) * (t100 * t196 + t101 * t199) + (-qJD(3) * t3 + t275 * t61) * t197 + (t196 * t331 + t199 * t332) * t181 + (-t181 * t240 - t200 * t4) * qJD(6), t240 * t146 - t249 * t200 + t57 * t63 - t82 * t228 + t6 * t117 + t16 * t52 - g(1) * (-t196 * t50 + t199 * t49) - g(2) * (-t196 * t48 + t199 * t47) - g(3) * (t100 * t199 - t101 * t196) + (qJD(3) * t4 + t275 * t63) * t197 + (-t196 * t332 + t199 * t331) * t181 + (t181 * t241 - t200 * t3) * qJD(6); 0, 0, 0, 0, -t197 * t203 * t200, t303 * t203, t286, t189, qJDD(3), qJD(3) * t96 - t151 * t300 + t210, -t248 + (-qJD(2) * t151 - t107) * t200 + t223, -pkin(3) * t98 - t142 * t96 + t218 * t195 + (-t197 * t29 + t200 * t45 + (t200 * t235 - t277) * t193) * qJD(2) + t304, t252 - pkin(3) * t99 - t144 * t96 - t218 * t193 + (t197 * t30 - t200 * t46 + (-t277 + (qJD(4) + t235) * t200) * t195) * qJD(2), t142 * t46 + t144 * t45 + (t29 * t299 + t11) * t195 + (t299 * t30 - t10 + t239) * t193 + t213, t216 * pkin(3) + g(1) * t67 + g(2) * t66 + g(3) * t123 - t29 * t45 - t30 * t46 + t235 * t96 + (-t193 * t29 + t195 * t30) * qJD(4) + (-t10 * t193 + t11 * t195 - t223) * qJ(4), t152 * t98 + (-t56 - t293) * t142 + t219 * t195 + (t197 * t22 - t200 * t42 + (-t200 * t27 - t277) * t193) * qJD(2) + t304, t142 * t39 - t144 * t42 + (-t22 * t299 + t7) * t195 + (t23 * t299 + t239 + t8) * t193 + t213, -t252 + t144 * t56 - t152 * t99 + (t219 + t294) * t193 + (-t197 * t23 + t200 * t39 + (t277 + (-qJD(4) + t27) * t200) * t195) * qJD(2), t9 * t152 - t23 * t39 - t27 * t56 - t22 * t42 - g(1) * (qJ(4) * t77 - t67) - g(2) * (qJ(4) * t75 - t66) - g(3) * (qJ(4) * t134 - t123) + (pkin(4) * t224 + t7 * qJ(4) + t23 * qJD(4)) * t195 + (t8 * qJ(4) + qJ(5) * t224 + t22 * qJD(4) - t27 * qJD(5)) * t193, -t148 * t228 + t306 * t63, -t13 * t148 + t147 * t228 - t305 * t63 - t306 * t61, -t146 * t148 + t181 * t306 + t300 * t63, t146 * t147 - t181 * t305 - t300 * t61, t181 * t300 -(t158 * t199 - t159 * t196) * t146 + t136 * t13 + t3 * t300 + t260 * t61 + (t196 * t230 + t199 * t229) * t181 + t305 * t16 + t342 * t147 (t158 * t196 + t159 * t199) * t146 - t136 * t228 - t4 * t300 + t260 * t63 + (-t196 * t229 + t199 * t230) * t181 + t306 * t16 + t342 * t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, t54, t348, t142 * t30 + t144 * t29 + t204, t208, t348, -t54, t339 - t326 + t142 * t23 + (-qJD(5) - t22) * t144 + t204, 0, 0, 0, 0, 0, -t13 - t323, t228 + t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142 * t144 + t347, -t119 + t99, -t192 * t203 - t140, -g(1) * t34 - g(2) * t32 - g(3) * t72 + t144 * t27 + (-pkin(4) * t297 + t200 * t23) * qJD(2) + t234, 0, 0, 0, 0, 0, -t144 * t61 - t199 * t146 - t196 * t238, -t144 * t63 + t196 * t146 - t199 * t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63 * t61, -t61 ^ 2 + t63 ^ 2, -t228 + t324, -t13 + t323, -t146, -t16 * t63 - g(1) * (-t196 * t35 + t199 * t34) - g(2) * (-t196 * t33 + t199 * t32) - g(3) * t19 + t278 + t346 * t4, t16 * t61 - g(1) * (-t196 * t34 - t199 * t35) - g(2) * (-t196 * t32 - t199 * t33) + g(3) * t20 - t249 + t346 * t3;];
tau_reg  = t1;
