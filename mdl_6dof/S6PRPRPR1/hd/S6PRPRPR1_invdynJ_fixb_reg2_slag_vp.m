% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6PRPRPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:28:04
% EndTime: 2019-03-08 19:28:13
% DurationCPUTime: 5.96s
% Computational Cost: add. (6519->529), mult. (15210->715), div. (0->0), fcn. (12565->16), ass. (0->259)
t193 = sin(pkin(6));
t191 = sin(pkin(11));
t194 = cos(pkin(11));
t200 = sin(qJ(2));
t203 = cos(qJ(2));
t231 = t191 * t203 + t194 * t200;
t127 = t231 * t193;
t286 = t203 * t194;
t140 = t191 * t200 - t286;
t192 = sin(pkin(10));
t195 = cos(pkin(10));
t196 = cos(pkin(6));
t289 = t196 * t203;
t291 = t196 * t200;
t283 = -t191 * t289 - t194 * t291;
t80 = t192 * t140 + t195 * t283;
t85 = -t195 * t140 + t192 * t283;
t241 = g(1) * t85 - g(2) * t80;
t220 = -g(3) * t127 - t241;
t121 = qJD(1) * t127;
t199 = sin(qJ(4));
t275 = qJD(4) * t199;
t343 = pkin(4) * t275 - t121;
t202 = cos(qJ(4));
t175 = pkin(2) * t191 + pkin(8);
t285 = qJ(5) + t175;
t248 = qJD(4) * t285;
t120 = qJD(5) * t202 - t199 * t248;
t280 = qJD(1) * t193;
t256 = t200 * t280;
t152 = t191 * t256;
t255 = t203 * t280;
t124 = t194 * t255 - t152;
t190 = sin(pkin(12));
t216 = -qJD(5) * t199 - t202 * t248;
t317 = cos(pkin(12));
t251 = t317 * t202;
t222 = -t190 * t199 + t251;
t326 = t120 * t317 - t222 * t124 + t190 * t216;
t141 = t190 * t202 + t199 * t317;
t133 = t141 * qJD(4);
t136 = t222 * qJD(4);
t342 = pkin(5) * t133 - pkin(9) * t136 + t343;
t277 = qJD(2) * t193;
t254 = qJD(1) * t277;
t268 = qJDD(1) * t193;
t341 = t200 * t268 + t203 * t254;
t179 = pkin(4) * t202 + pkin(3);
t270 = qJD(2) * qJD(4);
t253 = t199 * t270;
t163 = t203 * t268;
t316 = qJDD(2) * pkin(2);
t125 = -t200 * t254 + t163 + t316;
t257 = -t194 * t125 + t191 * t341;
t46 = pkin(4) * t253 - qJDD(2) * t179 + qJDD(5) + t257;
t267 = t199 * qJDD(2);
t236 = -qJDD(2) * t251 + t190 * t267;
t94 = qJD(2) * t133 + t236;
t165 = qJD(2) * t251;
t212 = qJDD(2) * t141 - t190 * t253;
t95 = qJD(4) * t165 + t212;
t18 = pkin(5) * t94 - pkin(9) * t95 + t46;
t198 = sin(qJ(6));
t201 = cos(qJ(6));
t147 = qJD(2) * pkin(2) + t255;
t108 = t191 * t147 + t194 * t256;
t106 = qJD(2) * pkin(8) + t108;
t245 = qJ(5) * qJD(2) + t106;
t168 = qJD(1) * t196 + qJD(3);
t303 = t168 * t199;
t66 = t202 * t245 + t303;
t58 = t317 * t66;
t153 = t202 * t168;
t65 = -t199 * t245 + t153;
t60 = qJD(4) * pkin(4) + t65;
t29 = t190 * t60 + t58;
t25 = qJD(4) * pkin(9) + t29;
t276 = qJD(2) * t199;
t131 = t190 * t276 - t165;
t134 = t141 * qJD(2);
t107 = t147 * t194 - t152;
t92 = -qJD(2) * t179 + qJD(5) - t107;
t41 = pkin(5) * t131 - pkin(9) * t134 + t92;
t235 = t198 * t25 - t201 * t41;
t166 = t196 * qJDD(1) + qJDD(3);
t151 = t202 * t166;
t269 = qJD(2) * qJD(5);
t274 = qJD(4) * t202;
t69 = t191 * t125 + t194 * t341;
t64 = qJDD(2) * pkin(8) + t69;
t20 = qJDD(4) * pkin(4) + t151 - t245 * t274 + (-qJ(5) * qJDD(2) - qJD(4) * t168 - t269 - t64) * t199;
t263 = t199 * t166 + t168 * t274 + t202 * t64;
t26 = -t106 * t275 + t263;
t266 = t202 * qJDD(2);
t21 = t202 * t269 + (-t253 + t266) * qJ(5) + t26;
t6 = t190 * t20 + t317 * t21;
t4 = qJDD(4) * pkin(9) + t6;
t1 = -t235 * qJD(6) + t198 * t18 + t201 * t4;
t129 = qJD(6) + t131;
t340 = t129 * t235 + t1;
t10 = t198 * t41 + t201 * t25;
t2 = -qJD(6) * t10 + t201 * t18 - t198 * t4;
t339 = t10 * t129 + t2;
t113 = qJD(4) * t198 + t134 * t201;
t247 = t129 * t198;
t338 = t113 * t247;
t187 = qJ(4) + pkin(12);
t181 = sin(t187);
t295 = t193 * t200;
t126 = t191 * t295 - t193 * t286;
t218 = t140 * t196;
t81 = -t192 * t231 - t195 * t218;
t84 = t192 * t218 - t195 * t231;
t221 = g(1) * t84 + g(2) * t81 - g(3) * t126;
t337 = t221 * t181;
t119 = t126 * t201;
t290 = t196 * t202;
t307 = t127 * t199;
t102 = t290 - t307;
t103 = t127 * t202 + t196 * t199;
t45 = t190 * t102 + t103 * t317;
t32 = -t198 * t45 + t119;
t273 = qJD(6) * t198;
t287 = t201 * t136;
t226 = t141 * t273 - t287;
t304 = t141 * t201;
t93 = qJDD(6) + t94;
t336 = -t129 * t226 + t93 * t304;
t335 = t134 ^ 2;
t334 = pkin(2) * t194;
t333 = pkin(4) * t199;
t154 = -t179 - t334;
t86 = -pkin(5) * t222 - pkin(9) * t141 + t154;
t138 = t285 * t202;
t249 = t285 * t199;
t90 = t138 * t317 - t190 * t249;
t37 = -t198 * t90 + t201 * t86;
t330 = qJD(6) * t37 + t342 * t198 + t201 * t326;
t38 = t198 * t86 + t201 * t90;
t329 = -qJD(6) * t38 - t198 * t326 + t342 * t201;
t271 = t201 * qJD(4);
t111 = t134 * t198 - t271;
t250 = -t201 * qJDD(4) + t198 * t95;
t43 = qJD(6) * t113 + t250;
t328 = -t111 * t287 - t43 * t304;
t327 = t120 * t190 - t141 * t124 - t216 * t317;
t42 = -qJD(6) * t271 - t198 * qJDD(4) + t134 * t273 - t201 * t95;
t325 = t113 * t133 + t222 * t42;
t323 = t190 * t66;
t322 = t198 * t42;
t320 = t198 * t93;
t319 = -t136 * t131 - t141 * t94;
t272 = qJD(6) * t201;
t318 = -t111 * t272 - t198 * t43;
t315 = t106 * t199;
t314 = t111 * t131;
t313 = t111 * t134;
t312 = t111 * t198;
t311 = t113 * t111;
t310 = t113 * t134;
t309 = t113 * t201;
t308 = t126 * t198;
t306 = t134 * t131;
t305 = t141 * t198;
t182 = cos(t187);
t302 = t182 * t198;
t301 = t182 * t201;
t299 = t192 * t193;
t298 = t192 * t200;
t297 = t193 * t195;
t296 = t193 * t199;
t294 = t193 * t202;
t293 = t193 * t203;
t288 = t198 * t136;
t284 = qJDD(1) - g(3);
t188 = t199 ^ 2;
t189 = t202 ^ 2;
t282 = t188 - t189;
t105 = -qJD(2) * pkin(3) - t107;
t279 = qJD(2) * t105;
t278 = qJD(2) * t124;
t262 = t113 * t288;
t261 = t195 * t294;
t260 = t195 * t289;
t205 = qJD(2) ^ 2;
t259 = t199 * t205 * t202;
t169 = pkin(2) * t293;
t197 = -qJ(5) - pkin(8);
t258 = -t126 * t179 - t127 * t197 + t169;
t5 = -t190 * t21 + t317 * t20;
t246 = t129 * t201;
t243 = t202 * t253;
t158 = pkin(2) * t260;
t242 = -pkin(2) * t298 + t158;
t240 = pkin(5) * t182 + pkin(9) * t181;
t239 = t10 * t201 + t198 * t235;
t238 = -t10 * t198 + t201 * t235;
t75 = t153 - t315;
t76 = t106 * t202 + t303;
t234 = t75 * t199 - t76 * t202;
t233 = -t111 * t133 + t222 * t43;
t33 = t201 * t45 + t308;
t232 = t133 * t134 - t222 * t95;
t230 = -t129 * t273 - t131 * t247 + t201 * t93;
t229 = t81 * t179 + t80 * t197 + t242;
t228 = -t192 * t289 - t195 * t200;
t28 = t317 * t60 - t323;
t227 = t141 * t272 + t288;
t48 = t181 * t80 - t182 * t297;
t50 = -t181 * t85 + t182 * t299;
t98 = -t127 * t181 + t182 * t196;
t225 = -g(1) * t50 - g(2) * t48 - g(3) * t98;
t49 = -t181 * t297 - t182 * t80;
t51 = t181 * t299 + t182 * t85;
t99 = t127 * t182 + t181 * t196;
t224 = -g(1) * t51 - g(2) * t49 - g(3) * t99;
t174 = pkin(4) * t190 + pkin(9);
t24 = -qJD(4) * pkin(5) - t28;
t223 = t129 * t24 - t174 * t93;
t219 = -g(1) * t299 + g(2) * t297 - g(3) * t196;
t217 = t228 * pkin(2);
t215 = -qJD(2) * t121 + t221;
t214 = t84 * t179 - t197 * t85 + t217;
t177 = -pkin(3) - t334;
t213 = -qJDD(4) * t175 + (qJD(2) * t177 + t105 + t124) * qJD(4);
t3 = -qJDD(4) * pkin(5) - t5;
t211 = qJD(6) * t129 * t174 - t225 + t3;
t210 = -g(1) * t228 - g(3) * t293;
t209 = qJD(6) * t238 + t1 * t201 - t2 * t198;
t208 = -t129 * t227 - t305 * t93;
t27 = -t76 * qJD(4) - t199 * t64 + t151;
t207 = -t27 * t199 + t26 * t202 + (-t199 * t76 - t202 * t75) * qJD(4);
t204 = qJD(4) ^ 2;
t63 = -qJDD(2) * pkin(3) + t257;
t206 = qJDD(2) * t177 + t175 * t204 + t215 + t63;
t176 = -pkin(4) * t317 - pkin(5);
t170 = pkin(4) * t290;
t157 = t192 * pkin(4) * t294;
t156 = qJDD(4) * t202 - t199 * t204;
t155 = qJDD(4) * t199 + t202 * t204;
t130 = t131 ^ 2;
t123 = t140 * t277;
t122 = qJD(2) * t127;
t97 = qJD(4) * t136 + qJDD(4) * t141;
t96 = -qJD(4) * t133 + qJDD(4) * t222;
t89 = t138 * t190 + t249 * t317;
t77 = pkin(4) * t276 + pkin(5) * t134 + pkin(9) * t131;
t53 = qJD(4) * t102 - t123 * t202;
t52 = -qJD(4) * t103 + t123 * t199;
t44 = -t102 * t317 + t103 * t190;
t31 = t317 * t65 - t323;
t30 = t190 * t65 + t58;
t23 = t190 * t52 + t317 * t53;
t22 = t190 * t53 - t317 * t52;
t14 = t198 * t77 + t201 * t31;
t13 = -t198 * t31 + t201 * t77;
t8 = -qJD(6) * t33 + t122 * t201 - t198 * t23;
t7 = t32 * qJD(6) + t122 * t198 + t201 * t23;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t284, 0, 0, 0, 0, 0, 0 (qJDD(2) * t203 - t200 * t205) * t193 (-qJDD(2) * t200 - t203 * t205) * t193, 0, -g(3) + (t196 ^ 2 + (t200 ^ 2 + t203 ^ 2) * t193 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -qJD(2) * t122 - qJDD(2) * t126, qJD(2) * t123 - qJDD(2) * t127, 0, -t107 * t122 - t108 * t123 + t126 * t257 + t127 * t69 + t166 * t196 - g(3), 0, 0, 0, 0, 0, 0, -t126 * t266 + qJD(4) * t52 + qJDD(4) * t102 + (-t122 * t202 + t126 * t275) * qJD(2), t126 * t267 - qJD(4) * t53 - qJDD(4) * t103 + (t122 * t199 + t126 * t274) * qJD(2) (-t102 * t199 + t103 * t202) * qJDD(2) + (-t199 * t52 + t202 * t53 + (-t102 * t202 - t103 * t199) * qJD(4)) * qJD(2), t102 * t27 + t103 * t26 + t105 * t122 + t126 * t63 + t52 * t75 + t53 * t76 - g(3), 0, 0, 0, 0, 0, 0, -qJD(4) * t22 - qJDD(4) * t44 + t122 * t131 + t126 * t94, -qJD(4) * t23 - qJDD(4) * t45 + t122 * t134 + t126 * t95, -t131 * t23 + t134 * t22 + t44 * t95 - t45 * t94, t122 * t92 + t126 * t46 - t22 * t28 + t23 * t29 - t44 * t5 + t45 * t6 - g(3), 0, 0, 0, 0, 0, 0, t111 * t22 + t129 * t8 + t32 * t93 + t43 * t44, t113 * t22 - t129 * t7 - t33 * t93 - t42 * t44, -t111 * t7 - t113 * t8 + t32 * t42 - t33 * t43, t1 * t33 + t10 * t7 + t2 * t32 + t22 * t24 - t235 * t8 + t3 * t44 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t163 - g(2) * (t260 - t298) + t210, -g(1) * (t192 * t291 - t195 * t203) - g(2) * (-t192 * t203 - t195 * t291) - t284 * t295, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t194 * t316 - t215 - t257, -t191 * t316 - t220 + t278 - t69, 0, -g(2) * t158 + t107 * t121 - t108 * t124 + (g(2) * t298 + t69 * t191 - t194 * t257 + t210) * pkin(2), qJDD(2) * t188 + 0.2e1 * t243, 0.2e1 * t199 * t266 - 0.2e1 * t270 * t282, t155, qJDD(2) * t189 - 0.2e1 * t243, t156, 0, t199 * t213 - t202 * t206, t199 * t206 + t202 * t213, t207 + t220 + (qJDD(2) * t175 - t278) * (t188 + t189) t63 * t177 - t105 * t121 - g(1) * (pkin(3) * t84 + pkin(8) * t85 + t217) - g(2) * (pkin(3) * t81 - pkin(8) * t80 + t242) - g(3) * (-pkin(3) * t126 + pkin(8) * t127 + t169) + t234 * t124 + t207 * t175, t134 * t136 + t141 * t95, -t232 + t319, t97, t131 * t133 - t222 * t94, t96, 0, -qJDD(4) * t89 - t121 * t131 + t133 * t92 - t222 * t46 + t154 * t94 - t221 * t182 + (t131 * t333 - t327) * qJD(4), -qJDD(4) * t90 - t121 * t134 + t136 * t92 + t141 * t46 + t154 * t95 + t337 + (t134 * t333 - t326) * qJD(4), -t131 * t326 - t133 * t29 + t134 * t327 - t136 * t28 - t141 * t5 + t222 * t6 + t89 * t95 - t90 * t94 + t220, -g(1) * t214 - g(2) * t229 - g(3) * t258 + t46 * t154 - t327 * t28 + t326 * t29 + t343 * t92 - t5 * t89 + t6 * t90, -t113 * t226 - t304 * t42, -t262 + (t322 + (-t309 + t312) * qJD(6)) * t141 + t328, t325 + t336, t111 * t227 + t305 * t43, t208 + t233, t129 * t133 - t222 * t93, t37 * t93 - t2 * t222 - t235 * t133 + t89 * t43 + t24 * t288 - g(1) * (t198 * t85 + t301 * t84) - g(2) * (-t198 * t80 + t301 * t81) - g(3) * (-t126 * t301 + t127 * t198) + (t3 * t198 + t24 * t272) * t141 + t329 * t129 + t327 * t111, -t38 * t93 + t1 * t222 - t10 * t133 - t89 * t42 + t24 * t287 - g(1) * (t201 * t85 - t302 * t84) - g(2) * (-t201 * t80 - t302 * t81) - g(3) * (t126 * t302 + t127 * t201) + (t3 * t201 - t24 * t273) * t141 - t330 * t129 + t327 * t113, t37 * t42 - t38 * t43 + t238 * t136 - t329 * t113 - t330 * t111 - t337 + (-qJD(6) * t239 - t1 * t198 - t2 * t201) * t141, t1 * t38 + t2 * t37 + t3 * t89 - g(1) * (t240 * t84 + t214) - g(2) * (t240 * t81 + t229) - g(3) * (-t126 * t240 + t258) - t329 * t235 + t327 * t24 + t330 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t219 + t166, 0, 0, 0, 0, 0, 0, t156, -t155, 0, -qJD(4) * t234 + t199 * t26 + t202 * t27 + t219, 0, 0, 0, 0, 0, 0, t96, -t97, t232 + t319, -t133 * t28 + t136 * t29 + t141 * t6 + t222 * t5 + t219, 0, 0, 0, 0, 0, 0, t208 - t233, t325 - t336, t262 + (-t322 + (t309 + t312) * qJD(6)) * t141 + t328, t133 * t24 + t136 * t239 + t141 * t209 - t222 * t3 + t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t259, t282 * t205, t267, t259, t266, qJDD(4), -g(3) * t102 + t151 + (-g(1) * t192 + g(2) * t195) * t294 + (t241 - t64 - t279) * t199, -t202 * t279 - g(1) * (-t192 * t296 - t202 * t85) - g(2) * (t195 * t296 + t202 * t80) + g(3) * t103 + (t75 + t315) * qJD(4) - t263, 0, 0, t306, -t130 + t335 (t165 + t131) * qJD(4) + t212, -t306, -t236, qJDD(4), t30 * qJD(4) - t92 * t134 + (qJDD(4) * t317 - t131 * t276) * pkin(4) + t225 + t5, qJD(4) * t31 + t131 * t92 + (-qJDD(4) * t190 - t134 * t276) * pkin(4) - t224 - t6 (t29 - t30) * t134 + (-t28 + t31) * t131 + (-t190 * t94 - t317 * t95) * pkin(4), -g(1) * t157 - g(3) * t170 + t28 * t30 - t29 * t31 + (g(2) * t261 + t5 * t317 + t6 * t190 + (-qJD(2) * t92 - t220) * t199) * pkin(4), t113 * t246 - t322 (-t42 - t314) * t201 - t338 + t318, t129 * t246 - t310 + t320, t111 * t247 - t201 * t43, t230 + t313, -t129 * t134, -t111 * t30 - t129 * t13 + t134 * t235 + t176 * t43 + t198 * t223 - t201 * t211, t10 * t134 - t113 * t30 + t129 * t14 - t176 * t42 + t198 * t211 + t201 * t223, t111 * t14 + t113 * t13 + (t131 * t235 - t174 * t43 + t1 + (t113 * t174 + t235) * qJD(6)) * t201 + (-t10 * t131 - t174 * t42 - t2 + (t111 * t174 - t10) * qJD(6)) * t198 + t224, t3 * t176 - t10 * t14 + t235 * t13 - t24 * t30 - g(1) * (pkin(5) * t50 + pkin(9) * t51 - t333 * t85 + t157) - g(2) * (pkin(5) * t48 + pkin(9) * t49 + (t199 * t80 - t261) * pkin(4)) - g(3) * (-pkin(4) * t307 + pkin(5) * t98 + pkin(9) * t99 + t170) + t209 * t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t134 * qJD(4) + t236 (t165 - t131) * qJD(4) + t212, -t130 - t335, t131 * t29 + t134 * t28 + t221 + t46, 0, 0, 0, 0, 0, 0, t230 - t313, -t129 ^ 2 * t201 - t310 - t320 (t42 - t314) * t201 + t338 + t318, -t24 * t134 + t340 * t198 + t339 * t201 + t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t311, -t111 ^ 2 + t113 ^ 2, t111 * t129 - t42, -t311, -t250 + (-qJD(6) + t129) * t113, t93, -t24 * t113 - g(1) * (-t198 * t51 - t201 * t84) - g(2) * (-t198 * t49 - t201 * t81) - g(3) * (-t198 * t99 + t119) + t339, t24 * t111 - g(1) * (t198 * t84 - t201 * t51) - g(2) * (t198 * t81 - t201 * t49) - g(3) * (-t201 * t99 - t308) - t340, 0, 0;];
tau_reg  = t9;
