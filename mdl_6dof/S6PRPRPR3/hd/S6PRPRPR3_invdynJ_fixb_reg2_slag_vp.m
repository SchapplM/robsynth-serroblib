% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:37:45
% EndTime: 2019-03-08 19:37:52
% DurationCPUTime: 4.26s
% Computational Cost: add. (4110->462), mult. (9439->595), div. (0->0), fcn. (7608->12), ass. (0->249)
t169 = sin(qJ(4));
t268 = qJD(4) * t169;
t164 = sin(pkin(6));
t162 = sin(pkin(11));
t165 = cos(pkin(11));
t170 = sin(qJ(2));
t173 = cos(qJ(2));
t217 = t162 * t173 + t165 * t170;
t92 = t217 * t164;
t84 = qJD(1) * t92;
t343 = pkin(4) * t268 - t84;
t261 = qJD(2) * qJD(4);
t243 = t169 * t261;
t172 = cos(qJ(4));
t254 = t172 * qJDD(2);
t339 = -t243 + t254;
t306 = qJ(5) * t172;
t225 = pkin(9) * t169 - t306;
t262 = t169 * qJD(5);
t189 = qJD(4) * t225 - t262;
t342 = t189 + t343;
t167 = cos(pkin(6));
t145 = qJD(1) * t167 + qJD(3);
t127 = t172 * t145;
t275 = qJD(1) * t164;
t249 = t173 * t275;
t119 = qJD(2) * pkin(2) + t249;
t250 = t170 * t275;
t74 = t162 * t119 + t165 * t250;
t72 = qJD(2) * pkin(8) + t74;
t50 = t169 * t72 - t127;
t341 = -qJD(5) - t50;
t274 = qJD(2) * t164;
t244 = qJD(1) * t274;
t260 = qJDD(1) * t164;
t340 = t170 * t260 + t173 * t244;
t168 = sin(qJ(6));
t171 = cos(qJ(6));
t267 = qJD(4) * t171;
t247 = t169 * t267;
t263 = qJD(6) * t172;
t338 = t168 * t263 + t247;
t174 = -pkin(4) - pkin(9);
t284 = t169 * qJ(5);
t240 = -pkin(3) - t284;
t199 = t172 * t174 + t240;
t140 = t173 * t260;
t303 = qJDD(2) * pkin(2);
t88 = -t170 * t244 + t140 + t303;
t253 = t162 * t340 - t165 * t88;
t232 = pkin(4) * t243 + t253;
t14 = qJD(2) * t189 + qJDD(2) * t199 + t232;
t280 = qJD(5) - t127 + (pkin(5) * qJD(2) + t72) * t169;
t22 = qJD(4) * t174 + t280;
t126 = t162 * t250;
t73 = t165 * t119 - t126;
t44 = qJD(2) * t199 - t73;
t222 = t168 * t44 - t171 * t22;
t242 = t172 * t261;
t255 = t169 * qJDD(2);
t200 = t242 + t255;
t141 = t167 * qJDD(1) + qJDD(3);
t266 = qJD(4) * t172;
t42 = t162 * t88 + t165 * t340;
t37 = qJDD(2) * pkin(8) + t42;
t233 = -t172 * t141 + t145 * t268 + t169 * t37 + t72 * t266;
t218 = qJDD(5) + t233;
t6 = t200 * pkin(5) + qJDD(4) * t174 + t218;
t1 = -t222 * qJD(6) + t171 * t14 + t168 * t6;
t273 = qJD(2) * t169;
t148 = qJD(6) + t273;
t337 = t148 * t222 + t1;
t13 = t168 * t22 + t171 * t44;
t2 = -qJD(6) * t13 - t168 * t14 + t171 * t6;
t336 = t13 * t148 + t2;
t163 = sin(pkin(10));
t166 = cos(pkin(10));
t282 = t173 * t165;
t109 = t162 * t170 - t282;
t201 = t167 * t109;
t58 = -t163 * t217 - t166 * t201;
t311 = t172 * t58;
t335 = pkin(4) * t311 + t58 * t284;
t61 = t163 * t201 - t166 * t217;
t310 = t172 * t61;
t334 = pkin(4) * t310 + t61 * t284;
t293 = t164 * t170;
t91 = t162 * t293 - t164 * t282;
t309 = t172 * t91;
t333 = -pkin(4) * t309 - t91 * t284;
t40 = -qJD(4) * pkin(4) - t341;
t269 = qJD(4) * t168;
t272 = qJD(2) * t172;
t111 = t171 * t272 + t269;
t270 = qJD(4) * t111;
t108 = qJDD(6) + t200;
t99 = t171 * t108;
t332 = t270 - t99;
t51 = t145 * t169 + t172 * t72;
t43 = -qJD(4) * qJ(5) - t51;
t214 = t111 * t148;
t64 = qJD(6) * t111 - t171 * qJDD(4) + t168 * t339;
t331 = t64 - t214;
t245 = t168 * t272;
t113 = -t245 + t267;
t298 = t113 * t148;
t65 = -qJD(6) * t245 + t168 * qJDD(4) + (qJD(4) * qJD(6) + t339) * t171;
t330 = -t65 + t298;
t213 = -pkin(4) * t172 + t240;
t54 = qJD(2) * t213 - t73;
t329 = t54 * t273 + qJDD(5);
t153 = pkin(5) * t272;
t35 = t153 - t43;
t328 = t108 * t174 + t148 * t35;
t248 = t168 * t268;
t197 = -t171 * t263 + t248;
t285 = t168 * t172;
t327 = t108 * t285 - t148 * t197;
t288 = t167 * t173;
t289 = t167 * t170;
t277 = -t162 * t288 - t165 * t289;
t62 = -t166 * t109 + t163 * t277;
t57 = t163 * t109 + t166 * t277;
t321 = t165 * pkin(2);
t101 = t213 - t321;
t207 = -qJ(5) * t266 - t262;
t17 = qJD(2) * t207 + qJDD(2) * t213 + t232;
t151 = pkin(2) * t162 + pkin(8);
t175 = qJD(4) ^ 2;
t203 = g(1) * t61 + g(2) * t58 - g(3) * t91;
t193 = t151 * t175 + t203;
t316 = t207 + t343;
t326 = qJD(2) * t316 + qJDD(2) * t101 + t17 + t193;
t70 = t167 * t169 + t172 * t92;
t86 = t109 * t274;
t20 = qJD(4) * t70 - t86 * t169;
t69 = -t167 * t172 + t169 * t92;
t85 = qJD(2) * t92;
t325 = qJD(2) * (-t172 * t85 + t268 * t91) - t20 * qJD(4) - t69 * qJDD(4) - t254 * t91;
t21 = -t92 * t268 + (qJD(4) * t167 - t86) * t172;
t324 = qJD(2) * (t169 * t85 + t266 * t91) - t21 * qJD(4) - t70 * qJDD(4) + t255 * t91;
t223 = -t13 * t171 - t168 * t222;
t323 = -qJD(6) * t223 + t1 * t168 + t2 * t171;
t322 = pkin(5) + pkin(8);
t320 = pkin(5) + t151;
t286 = t168 * t169;
t104 = t320 * t169;
t93 = t199 - t321;
t55 = t104 * t171 - t168 * t93;
t87 = t165 * t249 - t126;
t105 = t320 * t172;
t97 = qJD(4) * t105;
t319 = qJD(6) * t55 + t168 * t97 + t171 * t342 - t286 * t87;
t283 = t169 * t171;
t56 = t104 * t168 + t171 * t93;
t318 = -qJD(6) * t56 - t168 * t342 + t171 * t97 - t283 * t87;
t317 = t113 * t266 - t64 * t169;
t313 = t169 * t65;
t312 = t171 * t64;
t308 = t65 * t168;
t291 = t164 * t173;
t307 = pkin(2) * t291 - t91 * pkin(3);
t305 = qJD(2) * t84;
t304 = qJD(2) * t87;
t302 = qJDD(4) * pkin(4);
t300 = t111 * t171;
t299 = t113 * t111;
t297 = t113 * t168;
t295 = t163 * t170;
t294 = t164 * t169;
t292 = t164 * t172;
t287 = t168 * t108;
t281 = t35 * qJD(4);
t279 = qJDD(1) - g(3);
t278 = t338 * t148;
t160 = t169 ^ 2;
t161 = t172 ^ 2;
t276 = t160 - t161;
t265 = qJD(6) * t168;
t264 = qJD(6) * t171;
t259 = qJDD(2) * t160;
t258 = qJDD(2) * t161;
t257 = qJDD(4) * qJ(5);
t256 = qJDD(4) * t151;
t252 = t166 * t288;
t31 = t166 * t292 - t169 * t57;
t32 = -t166 * t294 - t172 * t57;
t238 = -t31 * pkin(4) + qJ(5) * t32;
t33 = -t163 * t292 + t169 * t62;
t34 = t163 * t294 + t172 * t62;
t237 = -t33 * pkin(4) + qJ(5) * t34;
t236 = -t69 * pkin(4) + qJ(5) * t70;
t234 = -t169 * t141 - t145 * t266 - t172 * t37 + t72 * t268;
t231 = pkin(8) * t92 + t307;
t230 = t113 * t247;
t228 = t169 * t242;
t224 = t13 * t168 - t171 * t222;
t23 = -t168 * t91 + t171 * t69;
t24 = t168 * t69 + t171 * t91;
t221 = t169 * t40 - t172 * t43;
t220 = t169 * t50 + t172 * t51;
t131 = pkin(2) * t252;
t216 = -pkin(2) * t295 + t58 * pkin(3) + t131;
t215 = t148 * t168;
t208 = -t163 * t288 - t166 * t170;
t206 = g(1) * t33 + g(2) * t31 + g(3) * t69;
t205 = -g(1) * t34 - g(2) * t32 - g(3) * t70;
t204 = -g(1) * t62 + g(2) * t57 - g(3) * t92;
t202 = -g(3) * t167 + (-g(1) * t163 + g(2) * t166) * t164;
t198 = -pkin(8) * t57 + t216;
t196 = pkin(2) * t208 + t61 * pkin(3);
t192 = qJD(4) * qJD(5) - t234 + t257;
t152 = -pkin(3) - t321;
t71 = -qJD(2) * pkin(3) - t73;
t191 = -t256 + (qJD(2) * t152 + t71 + t87) * qJD(4);
t190 = t256 + (-qJD(2) * t101 - t54 - t87) * qJD(4);
t188 = t206 - t233;
t187 = t205 - t234;
t7 = pkin(5) * t339 + t192;
t186 = -qJD(6) * t148 * t174 + t205 + t7;
t185 = pkin(8) * t62 + t196;
t184 = -g(1) * t208 - g(3) * t291;
t9 = t218 - t302;
t183 = t9 * t169 + t192 * t172 + (t169 * t43 + t172 * t40) * qJD(4);
t182 = -t234 * t172 + t233 * t169 + (-t169 * t51 + t172 * t50) * qJD(4);
t180 = qJD(4) * t51 + t188;
t179 = (-t160 - t161) * t304 + t204 + (t259 + t258) * t151;
t36 = -qJDD(2) * pkin(3) + t253;
t178 = qJDD(2) * t152 + t193 - t305 + t36;
t177 = (t169 * t69 + t172 * t70) * qJDD(2) + (t169 * t20 + t172 * t21 + (-t169 * t70 + t172 * t69) * qJD(4)) * qJD(2);
t176 = qJD(2) ^ 2;
t155 = pkin(4) * t273;
t144 = t169 * t176 * t172;
t130 = t276 * t176;
t129 = qJDD(4) * t172 - t175 * t169;
t128 = qJDD(4) * t169 + t172 * t175;
t114 = -qJ(5) * t272 + t155;
t103 = -0.2e1 * t228 + t258;
t102 = 0.2e1 * t228 + t259;
t96 = t320 * t268;
t94 = qJD(2) * t225 + t155;
t90 = t111 * t248;
t80 = 0.2e1 * t169 * t254 - 0.2e1 * t261 * t276;
t39 = t153 + t51;
t19 = t168 * t39 + t171 * t94;
t18 = -t168 * t94 + t171 * t39;
t4 = qJD(6) * t23 + t20 * t168 + t85 * t171;
t3 = -qJD(6) * t24 - t85 * t168 + t20 * t171;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t279, 0, 0, 0, 0, 0, 0 (qJDD(2) * t173 - t170 * t176) * t164 (-qJDD(2) * t170 - t173 * t176) * t164, 0, -g(3) + (t167 ^ 2 + (t170 ^ 2 + t173 ^ 2) * t164 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -qJD(2) * t85 - qJDD(2) * t91, qJD(2) * t86 - qJDD(2) * t92, 0, t141 * t167 + t253 * t91 + t42 * t92 - t73 * t85 - t74 * t86 - g(3), 0, 0, 0, 0, 0, 0, t325, t324, t177, t20 * t50 + t21 * t51 + t233 * t69 - t234 * t70 + t36 * t91 + t71 * t85 - g(3), 0, 0, 0, 0, 0, 0, t177, -t325, -t324, t17 * t91 + t192 * t70 + t20 * t40 - t21 * t43 + t54 * t85 + t69 * t9 - g(3), 0, 0, 0, 0, 0, 0, t108 * t23 + t111 * t21 + t148 * t3 + t65 * t70, -t108 * t24 + t113 * t21 - t148 * t4 - t64 * t70, -t111 * t4 - t113 * t3 + t23 * t64 - t24 * t65, t1 * t24 + t13 * t4 + t2 * t23 + t21 * t35 - t222 * t3 + t7 * t70 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t140 - g(2) * (t252 - t295) + t184, -g(1) * (t163 * t289 - t166 * t173) - g(2) * (-t163 * t173 - t166 * t289) - t279 * t293, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t165 * t303 - t203 - t253 + t305, -t162 * t303 - t204 + t304 - t42, 0, -g(2) * t131 + t73 * t84 - t74 * t87 + (g(2) * t295 + t42 * t162 - t165 * t253 + t184) * pkin(2), t102, t80, t128, t103, t129, 0, t169 * t191 - t172 * t178, t169 * t178 + t172 * t191, t179 + t182, -g(1) * t185 - g(2) * t198 - g(3) * t231 + t151 * t182 + t36 * t152 - t220 * t87 - t71 * t84, 0, -t128, -t129, t102, t80, t103, t179 + t183, t190 * t169 + t172 * t326, -t169 * t326 + t190 * t172, t17 * t101 - g(1) * (t185 + t334) - g(2) * (t198 + t335) - g(3) * (t231 + t333) - t221 * t87 + t316 * t54 + t183 * t151, t113 * t197 + t285 * t64, t230 - t90 + (t308 + t312 + (t297 + t300) * qJD(6)) * t172, t317 - t327, t65 * t171 * t172 - t111 * t338, -t313 + (-t270 - t99) * t172 + t278, t108 * t169 + t148 * t266, t105 * t65 + t55 * t108 - t96 * t111 + t204 * t171 + (-t168 * t203 - t267 * t35 + t2) * t169 + t318 * t148 + (-qJD(4) * t222 - t111 * t87 + t7 * t171 - t265 * t35) * t172, -t105 * t64 - t56 * t108 - t96 * t113 - t204 * t168 + (-t171 * t203 + t269 * t35 - t1) * t169 - t319 * t148 + (-qJD(4) * t13 - t113 * t87 - t7 * t168 - t264 * t35) * t172, t55 * t64 - t56 * t65 - t318 * t113 - t319 * t111 - t223 * t268 + (qJD(6) * t224 - t1 * t171 + t168 * t2 - t203) * t172, t1 * t56 + t2 * t55 + t7 * t105 - g(1) * (pkin(9) * t310 + t322 * t62 + t196 + t334) - g(2) * (pkin(9) * t311 - t322 * t57 + t216 + t335) - g(3) * (-pkin(9) * t309 + t322 * t92 + t307 + t333) + (-t172 * t87 - t96) * t35 + t319 * t13 - t318 * t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202 + t141, 0, 0, 0, 0, 0, 0, t129, -t128, 0, qJD(4) * t220 - t169 * t234 - t172 * t233 + t202, 0, 0, 0, 0, 0, 0, 0, -t129, t128, qJD(4) * t221 + t169 * t192 - t9 * t172 + t202, 0, 0, 0, 0, 0, 0, t172 * t332 + t278 + t313, t317 + t327, -t230 - t90 + (t308 - t312 + (-t297 + t300) * qJD(6)) * t172 (qJD(4) * t224 + t7) * t169 + (t281 - t323) * t172 + t202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, t130, t255, t144, t254, qJDD(4), -t273 * t71 + t180, -t50 * qJD(4) - t272 * t71 - t187, 0, 0, qJDD(4), -t255, -t254, -t144, t130, t144 (-pkin(4) * t169 + t306) * qJDD(2), -t114 * t272 - t180 - 0.2e1 * t302 + t329, 0.2e1 * t257 + (0.2e1 * qJD(5) + t50) * qJD(4) + (t114 * t169 + t172 * t54) * qJD(2) + t187, -t9 * pkin(4) - g(1) * t237 - g(2) * t238 - g(3) * t236 + t192 * qJ(5) - t54 * t114 + t341 * t43 - t40 * t51, -t113 * t215 - t312 (-t65 - t298) * t171 + (t64 + t214) * t168, -t148 * t265 + t99 + (-t113 * t172 - t148 * t286) * qJD(2), t171 * t214 + t308, -t148 * t264 - t287 + (t111 * t172 - t148 * t283) * qJD(2), -t148 * t272, qJ(5) * t65 + t280 * t111 - t18 * t148 + t186 * t168 + t171 * t328 + t222 * t272, -qJ(5) * t64 + t280 * t113 + t13 * t272 + t19 * t148 - t168 * t328 + t186 * t171, t19 * t111 + t18 * t113 + (-t13 * t273 + t174 * t64 - t2 + (-t111 * t174 - t13) * qJD(6)) * t171 + (-t222 * t273 - t174 * t65 - t1 + (t113 * t174 - t222) * qJD(6)) * t168 + t206, t7 * qJ(5) - t13 * t19 + t222 * t18 - g(1) * (-pkin(9) * t33 + t237) - g(2) * (-pkin(9) * t31 + t238) - g(3) * (-pkin(9) * t69 + t236) + t280 * t35 + t323 * t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t255, qJDD(4) + t144, -t160 * t176 - t175, qJD(4) * t43 - t188 - t302 + t329, 0, 0, 0, 0, 0, 0, -t148 * t215 - t332, -t148 ^ 2 * t171 - qJD(4) * t113 - t287, t168 * t330 + t171 * t331, t168 * t337 + t336 * t171 - t206 - t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t299, -t111 ^ 2 + t113 ^ 2, -t331, -t299, t330, t108, -t35 * t113 - g(1) * (t168 * t61 + t171 * t33) - g(2) * (t168 * t58 + t171 * t31) - g(3) * t23 + t336, t35 * t111 - g(1) * (-t168 * t33 + t171 * t61) - g(2) * (-t168 * t31 + t171 * t58) + g(3) * t24 - t337, 0, 0;];
tau_reg  = t5;
