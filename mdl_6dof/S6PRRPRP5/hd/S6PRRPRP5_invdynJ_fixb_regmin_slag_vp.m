% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:49:28
% EndTime: 2019-03-08 21:49:34
% DurationCPUTime: 4.27s
% Computational Cost: add. (3450->470), mult. (7553->607), div. (0->0), fcn. (5526->10), ass. (0->244)
t161 = sin(pkin(6));
t165 = sin(qJ(2));
t168 = cos(qJ(2));
t259 = qJD(1) * qJD(2);
t240 = t168 * t259;
t162 = cos(pkin(6));
t268 = qJD(3) * t162;
t300 = qJDD(2) * pkin(8);
t345 = qJD(1) * t268 + t300 + (qJDD(1) * t165 + t240) * t161;
t164 = sin(qJ(3));
t258 = qJD(2) * qJD(3);
t238 = t164 * t258;
t167 = cos(qJ(3));
t254 = t167 * qJDD(2);
t343 = -t238 + t254;
t323 = pkin(4) + pkin(8);
t234 = -qJ(4) * t164 - pkin(2);
t277 = qJD(1) * t161;
t243 = t165 * t277;
t120 = qJD(2) * pkin(8) + t243;
t276 = qJD(1) * t162;
t73 = t164 * t120 - t167 * t276;
t344 = qJD(4) + t73;
t292 = t161 * t167;
t104 = t162 * t164 + t165 * t292;
t171 = qJD(2) ^ 2;
t207 = qJDD(2) * t168 - t165 * t171;
t236 = t168 * t258;
t291 = t161 * t168;
t245 = qJD(2) * t291;
t293 = t161 * t165;
t250 = t164 * t293;
t57 = -qJD(3) * t250 + (t245 + t268) * t167;
t342 = qJD(3) * t57 + qJDD(3) * t104 + t161 * (t164 * t207 + t167 * t236);
t103 = -t162 * t167 + t250;
t58 = qJD(3) * t104 + t164 * t245;
t341 = -qJD(3) * t58 - qJDD(3) * t103 + t161 * (-t164 * t236 + t167 * t207);
t157 = qJD(3) * qJ(4);
t74 = t167 * t120 + t164 * t276;
t66 = -t157 - t74;
t169 = -pkin(3) - pkin(9);
t108 = t169 * t167 + t234;
t264 = qJD(3) * t167;
t118 = t323 * t264;
t131 = t323 * t164;
t163 = sin(qJ(5));
t166 = cos(qJ(5));
t261 = qJD(5) * t166;
t262 = qJD(5) * t163;
t266 = qJD(3) * t164;
t150 = pkin(3) * t266;
t304 = qJ(4) * t167;
t217 = pkin(9) * t164 - t304;
t263 = qJD(4) * t164;
t183 = qJD(3) * t217 - t263;
t76 = t150 + t183;
t286 = t164 * t168;
t80 = (t163 * t286 + t165 * t166) * t161;
t340 = qJD(1) * t80 + t108 * t262 - t163 * t118 - t131 * t261 - t166 * t76;
t65 = -qJD(3) * pkin(3) + t344;
t267 = qJD(3) * t163;
t271 = qJD(2) * t167;
t111 = t166 * t271 + t267;
t273 = qJD(2) * t164;
t146 = qJD(5) + t273;
t298 = t111 * t146;
t40 = qJD(5) * t111 - t166 * qJDD(3) + t163 * t343;
t339 = t40 - t298;
t242 = t163 * t271;
t265 = qJD(3) * t166;
t113 = -t242 + t265;
t296 = t113 * t146;
t41 = -qJD(5) * t242 + qJDD(3) * t163 + (qJD(3) * qJD(5) + t343) * t166;
t338 = -t41 + t296;
t210 = t146 * t163;
t270 = qJD(3) * t111;
t237 = t167 * t258;
t255 = t164 * qJDD(2);
t191 = t237 + t255;
t110 = qJDD(5) + t191;
t96 = t166 * t110;
t337 = -t146 * t210 - t270 + t96;
t336 = -t146 * t262 + t96;
t285 = t166 * t168;
t334 = t163 * t165 - t164 * t285;
t282 = pkin(4) * t273 + t344;
t275 = qJD(1) * t168;
t247 = t161 * t275;
t123 = -pkin(3) * t167 + t234;
t274 = qJD(2) * t123;
t75 = -t247 + t274;
t333 = t75 * t273 + qJDD(4);
t332 = t120 * t264 + t345 * t164;
t307 = cos(pkin(10));
t233 = t161 * t307;
t232 = t307 * t165;
t160 = sin(pkin(10));
t294 = t160 * t168;
t98 = t162 * t232 + t294;
t54 = -t164 * t233 + t167 * t98;
t231 = t307 * t168;
t295 = t160 * t165;
t100 = -t162 * t295 + t231;
t56 = t160 * t161 * t164 + t100 * t167;
t194 = g(1) * t56 + g(2) * t54 + g(3) * t104;
t241 = t165 * t259;
t216 = -qJDD(1) * t291 + t161 * t241;
t97 = -t162 * t231 + t295;
t99 = t162 * t294 + t232;
t223 = g(1) * t99 + g(2) * t97;
t170 = qJD(3) ^ 2;
t317 = pkin(8) * t170;
t331 = 0.2e1 * qJDD(2) * pkin(2) + t161 * (-g(3) * t168 + t241) - t216 + t223 - t317;
t269 = qJD(3) * t113;
t290 = t163 * t110;
t324 = t146 ^ 2;
t329 = t324 * t166 + t269 + t290;
t189 = -g(3) * t291 + t223;
t256 = qJDD(2) * t123;
t198 = -qJ(4) * t264 - t263;
t203 = pkin(3) * t238 + t216;
t26 = qJD(2) * t198 + t203 + t256;
t95 = t150 + t198;
t327 = qJD(2) * (-t95 + t243) + t189 - t256 - t26 - t317;
t213 = t166 * t108 + t163 * t131;
t326 = -t213 * qJD(5) + t118 * t166 - t163 * t76 + t334 * t277;
t325 = t113 ^ 2;
t320 = qJ(6) * t264 + qJD(6) * t164 - t340;
t319 = -pkin(5) * t264 - t326;
t318 = pkin(5) * t110;
t42 = t169 * qJD(3) + t282;
t63 = qJD(2) * t108 - t247;
t15 = t163 * t42 + t166 * t63;
t8 = qJ(6) * t146 + t15;
t315 = t146 * t8;
t62 = pkin(4) * t271 + t74;
t151 = pkin(3) * t273;
t90 = qJD(2) * t217 + t151;
t314 = t163 * t62 + t166 * t90;
t219 = pkin(5) * t166 + qJ(6) * t163;
t205 = -pkin(4) - t219;
t313 = qJD(5) * t219 - qJD(6) * t166 - t205 * t273 + t344;
t312 = qJD(2) * pkin(2);
t311 = t146 * t15;
t310 = t166 * t40;
t309 = t167 * t97;
t308 = t167 * t99;
t306 = pkin(8) * qJDD(3);
t303 = qJ(6) * t110;
t302 = qJD(3) * t74;
t299 = qJDD(3) * pkin(3);
t297 = t113 * t111;
t289 = t163 * t164;
t287 = t164 * t166;
t284 = t167 * t168;
t283 = t169 * t110;
t14 = -t163 * t63 + t166 * t42;
t281 = qJD(6) - t14;
t280 = qJDD(1) - g(3);
t132 = t323 * t167;
t158 = t164 ^ 2;
t159 = t167 ^ 2;
t279 = t158 - t159;
t278 = t158 + t159;
t272 = qJD(2) * t165;
t260 = qJD(5) * t169;
t257 = qJDD(1) * t162;
t190 = -t167 * t257 + t332;
t188 = qJDD(4) + t190;
t10 = t191 * pkin(4) + t169 * qJDD(3) + t188;
t18 = qJD(2) * t183 + qJDD(2) * t108 + t203;
t253 = -t163 * t10 - t166 * t18 - t42 * t261;
t252 = -pkin(3) * t309 + t234 * t97;
t251 = -pkin(3) * t308 + t234 * t99;
t248 = t164 * t171 * t167;
t246 = t161 * t272;
t235 = t166 * t10 - t163 * t18 - t63 * t261 - t42 * t262;
t229 = -t120 * t266 + t164 * t257 + t345 * t167;
t48 = t157 + t62;
t228 = pkin(2) * t291 + pkin(8) * t293 + (pkin(3) * t284 + qJ(4) * t286) * t161;
t227 = t111 * t247;
t226 = t113 * t247;
t222 = g(1) * t100 + g(2) * t98;
t7 = -pkin(5) * t146 + t281;
t220 = t163 * t7 + t166 * t8;
t218 = -pkin(5) * t163 + qJ(6) * t166;
t212 = -t108 * t163 + t131 * t166;
t155 = qJDD(3) * qJ(4);
t156 = qJD(3) * qJD(4);
t19 = -t155 - t156 - t229;
t122 = qJ(4) - t218;
t201 = -t262 * t63 - t253;
t60 = -t103 * t163 + t161 * t285;
t59 = t103 * t166 + t163 * t291;
t199 = -t146 * t261 - t290;
t33 = t163 * t98 + t287 * t97;
t35 = t100 * t163 + t287 * t99;
t79 = t334 * t161;
t197 = -g(1) * t35 - g(2) * t33 - g(3) * t79;
t34 = t166 * t98 - t289 * t97;
t36 = t100 * t166 - t289 * t99;
t196 = -g(1) * t36 - g(2) * t34 - g(3) * t80;
t53 = t98 * t164 + t167 * t233;
t55 = t100 * t164 - t160 * t292;
t195 = g(1) * t55 + g(2) * t53 + g(3) * t103;
t12 = -qJD(5) * t60 + t163 * t246 - t166 * t58;
t192 = t104 * t41 + t110 * t59 + t57 * t111 - t12 * t146;
t21 = pkin(5) * t111 - qJ(6) * t113 + t48;
t187 = t146 * t21 + t283;
t13 = qJD(5) * t59 + t163 * t58 + t166 * t246;
t185 = t104 * t40 - t110 * t60 - t113 * t57 + t13 * t146;
t27 = t163 * t97 - t53 * t166;
t29 = t163 * t99 - t55 * t166;
t184 = g(1) * t29 + g(2) * t27 - g(3) * t59 + t235;
t181 = -t146 * t260 - t194;
t11 = pkin(4) * t343 - t19;
t3 = pkin(5) * t41 + qJ(6) * t40 - qJD(6) * t113 + t11;
t180 = -t181 - t3;
t121 = -t247 - t312;
t179 = -t306 + (t121 + t247 - t312) * qJD(3);
t178 = t306 + (-t247 - t75 - t274) * qJD(3);
t28 = t163 * t53 + t166 * t97;
t30 = t163 * t55 + t166 * t99;
t177 = -g(1) * t30 - g(2) * t28 + g(3) * t60 + t201;
t176 = qJD(3) * t73 - t194 + t229;
t175 = -t190 + t195;
t174 = t113 * t21 + qJDD(6) - t184;
t20 = t188 - t299;
t172 = t164 * t20 - t167 * t19 + (t164 * t66 + t167 * t65) * qJD(3) - t222;
t117 = t323 * t266;
t116 = -qJ(4) * t271 + t151;
t94 = t103 * pkin(3);
t72 = t167 * t219 + t132;
t52 = pkin(5) * t113 + qJ(6) * t111;
t51 = t55 * pkin(3);
t50 = t53 * pkin(3);
t44 = -pkin(5) * t164 - t212;
t43 = qJ(6) * t164 + t213;
t31 = (qJD(5) * t218 + qJD(6) * t163) * t167 + (-pkin(8) + t205) * t266;
t23 = -pkin(5) * t271 + t163 * t90 - t166 * t62;
t22 = qJ(6) * t271 + t314;
t2 = qJDD(6) - t235 - t318;
t1 = qJD(6) * t146 + t201 + t303;
t4 = [t280, 0, t207 * t161 (-qJDD(2) * t165 - t168 * t171) * t161, 0, 0, 0, 0, 0, t341, -t342 (t103 * t164 + t104 * t167) * qJDD(2) + (t164 * t58 + t167 * t57 + (t103 * t167 - t104 * t164) * qJD(3)) * qJD(2), -t341, t342, t103 * t20 - t104 * t19 - t57 * t66 + t58 * t65 - g(3) + (-t168 * t26 + t272 * t75) * t161, 0, 0, 0, 0, 0, t192, -t185, t192, -t111 * t13 + t113 * t12 + t40 * t59 + t41 * t60, t185, -t1 * t60 + t104 * t3 + t12 * t7 + t13 * t8 - t2 * t59 + t21 * t57 - g(3); 0, qJDD(2), t280 * t291 + t223, -t280 * t293 + t222, qJDD(2) * t158 + 0.2e1 * t164 * t237, 0.2e1 * t164 * t254 - 0.2e1 * t258 * t279, qJDD(3) * t164 + t167 * t170, qJDD(3) * t167 - t164 * t170, 0, t179 * t164 + t167 * t331, -t164 * t331 + t179 * t167, t278 * t300 + (-g(3) * t165 - t240 * t278) * t161 + t172, t178 * t164 - t167 * t327, t164 * t327 + t178 * t167, t26 * t123 + t75 * t95 - g(1) * t251 - g(2) * t252 - g(3) * t228 + (-t165 * t75 + (-t164 * t65 + t167 * t66) * t168) * t277 + t172 * pkin(8), t163 * t167 * t40 + (t163 * t266 - t167 * t261) * t113 (-t111 * t163 + t113 * t166) * t266 + (t163 * t41 + t310 + (t111 * t166 + t113 * t163) * qJD(5)) * t167 (t146 * t267 - t40) * t164 + (t199 + t269) * t167 (t146 * t265 - t41) * t164 + (-t270 - t336) * t167, t110 * t164 + t146 * t264, t212 * t110 - t117 * t111 + t132 * t41 + (-t265 * t48 + t235) * t164 + t326 * t146 + (t14 * qJD(3) + t11 * t166 - t262 * t48 - t227) * t167 + t196, -t213 * t110 - t117 * t113 - t132 * t40 + ((qJD(3) * t48 + qJD(5) * t63) * t163 + t253) * t164 + t340 * t146 + (-t15 * qJD(3) - t11 * t163 - t261 * t48 - t226) * t167 - t197, -t110 * t44 + t111 * t31 + t41 * t72 + (-t21 * t265 - t2) * t164 - t319 * t146 + (-qJD(3) * t7 + t166 * t3 - t21 * t262 - t227) * t167 + t196, -t40 * t44 - t41 * t43 + t319 * t113 - t320 * t111 + t220 * t266 + (-t1 * t166 - t163 * t2 + (t163 * t8 - t166 * t7) * qJD(5) + t189) * t167, t110 * t43 - t113 * t31 + t40 * t72 + (-t21 * t267 + t1) * t164 + t320 * t146 + (qJD(3) * t8 + t163 * t3 + t21 * t261 + t226) * t167 + t197, t1 * t43 + t3 * t72 + t21 * t31 + t2 * t44 - g(1) * (pkin(5) * t36 - pkin(9) * t308 + qJ(6) * t35 + t323 * t100 + t251) - g(2) * (pkin(5) * t34 - pkin(9) * t309 + qJ(6) * t33 + t323 * t98 + t252) - g(3) * (pkin(5) * t80 + qJ(6) * t79 + t228) + t320 * t8 + t319 * t7 + (-t21 * t167 * t275 - g(3) * (pkin(4) * t165 + pkin(9) * t284)) * t161; 0, 0, 0, 0, -t248, t279 * t171, t255, t254, qJDD(3), -t121 * t273 + t175 + t302, -t121 * t271 - t176 (-pkin(3) * t164 + t304) * qJDD(2), -0.2e1 * t299 - t302 + (-qJD(2) * t116 - t257) * t167 - t195 + t332 + t333, 0.2e1 * t155 + 0.2e1 * t156 + (t116 * t164 + t167 * t75) * qJD(2) + t176, -t19 * qJ(4) - t20 * pkin(3) - t75 * t116 - t65 * t74 - g(1) * (qJ(4) * t56 - t51) - g(2) * (qJ(4) * t54 - t50) - g(3) * (qJ(4) * t104 - t94) - t344 * t66, -t113 * t210 - t310 (-t41 - t296) * t166 + (t40 + t298) * t163 (-t113 * t167 - t146 * t289) * qJD(2) + t336 (t111 * t167 - t146 * t287) * qJD(2) + t199, -t146 * t271, -t14 * t271 + qJ(4) * t41 + t282 * t111 + (t283 + (t48 - t62) * t146) * t166 + (t11 + (t90 - t260) * t146 - t194) * t163, -qJ(4) * t40 + t314 * t146 + t15 * t271 + t282 * t113 + (-t146 * t48 - t283) * t163 + (t11 + t181) * t166, t111 * t313 + t122 * t41 + t146 * t23 - t163 * t180 + t166 * t187 + t271 * t7, t111 * t22 - t113 * t23 + (-t8 * t273 + t169 * t40 + t2 + (-t111 * t169 - t8) * qJD(5)) * t166 + (-t7 * t273 - t169 * t41 - t1 + (t113 * t169 - t7) * qJD(5)) * t163 + t195, -t113 * t313 + t122 * t40 - t146 * t22 + t163 * t187 + t166 * t180 - t271 * t8, -t8 * t22 - t7 * t23 - g(1) * (-pkin(9) * t55 - t51) - g(2) * (-pkin(9) * t53 - t50) - g(3) * (-pkin(9) * t103 - t94) + t313 * t21 + (qJD(5) * t220 + t1 * t163 - t2 * t166) * t169 + (t3 - t194) * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t255, qJDD(3) + t248, -t158 * t171 - t170, qJD(3) * t66 - t175 - t299 + t333, 0, 0, 0, 0, 0, t337, -t329, t337, t163 * t338 + t166 * t339, t329, -qJD(3) * t21 + (-t2 + t315) * t166 + (t146 * t7 + t1) * t163 - t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t297, -t111 ^ 2 + t325, -t339, t338, t110, -t113 * t48 + t184 + t311, t111 * t48 + t14 * t146 - t177, -t111 * t52 - t174 + t311 + 0.2e1 * t318, pkin(5) * t40 - qJ(6) * t41 + (-t15 + t8) * t113 + (t7 - t281) * t111, 0.2e1 * t303 - t111 * t21 + t113 * t52 + (0.2e1 * qJD(6) - t14) * t146 + t177, t1 * qJ(6) - t2 * pkin(5) - t21 * t52 - t7 * t15 - g(1) * (-pkin(5) * t29 + qJ(6) * t30) - g(2) * (-pkin(5) * t27 + qJ(6) * t28) - g(3) * (pkin(5) * t59 - qJ(6) * t60) + t281 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110 + t297, -t339, -t324 - t325, t174 - t315 - t318;];
tau_reg  = t4;
