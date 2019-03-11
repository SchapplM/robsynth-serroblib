% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:08:21
% EndTime: 2019-03-09 08:08:32
% DurationCPUTime: 4.15s
% Computational Cost: add. (5131->477), mult. (12126->607), div. (0->0), fcn. (9048->12), ass. (0->238)
t197 = sin(pkin(9));
t201 = sin(qJ(2));
t204 = cos(qJ(2));
t300 = cos(pkin(9));
t157 = t197 * t204 + t300 * t201;
t196 = sin(pkin(10));
t198 = cos(pkin(10));
t200 = sin(qJ(6));
t203 = cos(qJ(6));
t328 = -t196 * t203 + t198 * t200;
t84 = t328 * t157;
t193 = qJ(2) + pkin(9);
t190 = cos(t193);
t181 = g(3) * t190;
t309 = qJ(3) + pkin(7);
t263 = qJD(2) * t309;
t135 = t204 * qJD(3) - t201 * t263;
t166 = t309 * t204;
t103 = t135 * qJD(1) + qJDD(1) * t166;
t136 = -t201 * qJD(3) - t204 * t263;
t165 = t309 * t201;
t93 = qJDD(2) * pkin(2) + t136 * qJD(1) - qJDD(1) * t165;
t50 = -t197 * t103 + t300 * t93;
t230 = qJDD(2) * pkin(3) - qJDD(4) + t50;
t267 = t230 - t181;
t189 = sin(t193);
t202 = sin(qJ(1));
t205 = cos(qJ(1));
t251 = g(1) * t205 + g(2) * t202;
t331 = t251 * t189;
t208 = -t331 - t267;
t261 = t300 * t204;
t175 = qJD(1) * t261;
t282 = qJD(1) * t201;
t138 = t197 * t282 - t175;
t277 = qJD(6) - t138;
t141 = t157 * qJD(1);
t117 = t198 * qJD(2) - t141 * t196;
t233 = qJD(2) * t196 + t198 * t141;
t57 = t117 * t203 + t200 * t233;
t338 = t277 * t57;
t59 = -t117 * t200 + t203 * t233;
t337 = t277 * t59;
t278 = qJD(6) * t203;
t279 = qJD(6) * t200;
t301 = t328 * t138 + t196 * t278 - t198 * t279;
t336 = t117 * t138;
t335 = t233 ^ 2;
t191 = t204 * pkin(2);
t186 = t191 + pkin(1);
t164 = -qJD(1) * t186 + qJD(3);
t68 = pkin(3) * t138 - qJ(4) * t141 + t164;
t160 = qJD(1) * t165;
t307 = qJD(2) * pkin(2);
t153 = -t160 + t307;
t161 = qJD(1) * t166;
t262 = t300 * t161;
t99 = t197 * t153 + t262;
t89 = qJD(2) * qJ(4) + t99;
t38 = t196 * t68 + t198 * t89;
t27 = t138 * qJ(5) + t38;
t140 = t157 * qJD(2);
t275 = t201 * qJDD(1);
t245 = -qJDD(1) * t261 + t197 * t275;
t100 = qJD(1) * t140 + t245;
t276 = qJD(1) * qJD(2);
t266 = t201 * t276;
t101 = qJD(2) * t175 + t157 * qJDD(1) - t197 * t266;
t218 = pkin(2) * t266 - qJDD(1) * t186 + qJDD(3);
t28 = pkin(3) * t100 - qJ(4) * t101 - qJD(4) * t141 + t218;
t51 = t300 * t103 + t197 * t93;
t45 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t51;
t10 = -t196 * t45 + t198 * t28;
t248 = qJDD(5) - t10;
t8 = -pkin(4) * t100 + t248;
t334 = -t138 * t27 + t8;
t330 = -t190 * pkin(3) - t189 * qJ(4);
t316 = g(1) * t202;
t265 = -g(2) * t205 + t316;
t280 = qJD(5) * t233;
t75 = qJDD(2) * t196 + t101 * t198;
t329 = qJ(5) * t75 + t280;
t327 = -qJD(6) + t277;
t74 = -t198 * qJDD(2) + t101 * t196;
t16 = t59 * qJD(6) + t200 * t75 - t203 * t74;
t156 = t196 * t200 + t198 * t203;
t302 = t277 * t156;
t94 = -qJDD(6) + t100;
t325 = t277 * t302 - t328 * t94;
t137 = t138 ^ 2;
t324 = t100 * t198 + t117 * t141 - t137 * t196;
t298 = t100 * t196;
t323 = t137 * t198 + t141 * t233 + t298;
t299 = qJ(5) * t198;
t320 = pkin(4) + pkin(5);
t322 = t320 * t196 - t299;
t211 = t230 + t329;
t7 = -t320 * t74 + t211;
t321 = t7 + t331;
t319 = pkin(4) * t74;
t318 = pkin(2) * t201;
t317 = pkin(8) * t196;
t313 = g(3) * t189;
t312 = g(3) * t204;
t311 = t198 * pkin(4);
t179 = pkin(2) * t197 + qJ(4);
t310 = -pkin(8) + t179;
t11 = t196 * t28 + t198 * t45;
t221 = -t197 * t201 + t261;
t143 = t221 * qJD(2);
t273 = t201 * t307;
t60 = pkin(3) * t140 - qJ(4) * t143 - qJD(4) * t157 + t273;
t82 = t300 * t135 + t197 * t136;
t30 = t196 * t60 + t198 * t82;
t146 = t197 * t161;
t105 = -t300 * t160 - t146;
t80 = pkin(2) * t282 + pkin(3) * t141 + qJ(4) * t138;
t47 = t198 * t105 + t196 * t80;
t305 = t141 * t57;
t303 = t59 * t141;
t112 = -t197 * t165 + t300 * t166;
t97 = -pkin(3) * t221 - qJ(4) * t157 - t186;
t54 = t198 * t112 + t196 * t97;
t297 = t179 * t198;
t296 = t189 * t205;
t295 = t190 * t205;
t294 = t196 * qJ(5);
t293 = t196 * t202;
t290 = t198 * t205;
t289 = t202 * t198;
t288 = t205 * t196;
t287 = t205 * t309;
t98 = t300 * t153 - t146;
t231 = qJD(2) * pkin(3) - qJD(4) + t98;
t286 = -qJD(4) - t231;
t216 = qJ(5) * t233 + t231;
t36 = -pkin(4) * t117 - t216;
t285 = qJD(4) - t36;
t284 = (g(1) * t290 + g(2) * t289) * t189;
t194 = t201 ^ 2;
t283 = -t204 ^ 2 + t194;
t281 = qJD(4) * t198;
t274 = t204 * qJDD(1);
t33 = t141 * qJ(5) + t47;
t43 = -qJ(5) * t221 + t54;
t272 = t100 * t297;
t173 = t205 * t186;
t271 = pkin(3) * t295 + qJ(4) * t296 + t173;
t9 = -t211 + t319;
t269 = -t9 - t181;
t2 = -pkin(8) * t75 - t320 * t100 + t248;
t6 = t100 * qJ(5) + t138 * qJD(5) + t11;
t5 = pkin(8) * t74 + t6;
t268 = t203 * t2 - t200 * t5;
t70 = t196 * t82;
t29 = t198 * t60 - t70;
t37 = -t196 * t89 + t198 * t68;
t95 = t196 * t105;
t46 = t198 * t80 - t95;
t107 = t196 * t112;
t53 = t198 * t97 - t107;
t104 = -t197 * t160 + t262;
t260 = qJD(5) * t196 - t322 * t138 + t104;
t81 = t135 * t197 - t300 * t136;
t111 = t300 * t165 + t197 * t166;
t259 = t277 ^ 2;
t17 = t140 * qJ(5) - qJD(5) * t221 + t30;
t258 = g(3) * (t191 - t330);
t257 = -g(2) * t296 + t189 * t316;
t256 = t156 * t94 - t301 * t277;
t185 = -t300 * pkin(2) - pkin(3);
t255 = qJD(5) - t37;
t254 = -pkin(3) * t189 - t318;
t126 = t190 * t293 + t290;
t128 = t190 * t288 - t289;
t253 = -g(1) * t126 + g(2) * t128;
t127 = t190 * t289 - t288;
t129 = t190 * t290 + t293;
t252 = g(1) * t127 - g(2) * t129;
t250 = -t196 * t6 + t198 * t8;
t249 = t200 * t2 + t203 * t5;
t247 = t143 * t36 + t157 * t9;
t246 = pkin(4) * t196 - t299;
t244 = -t10 * t198 - t11 * t196;
t14 = -pkin(8) * t233 - t320 * t138 + t255;
t19 = -pkin(8) * t117 + t27;
t3 = t14 * t203 - t19 * t200;
t4 = t14 * t200 + t19 * t203;
t243 = -t143 * t231 - t157 * t230;
t242 = -t196 * t37 + t198 * t38;
t24 = t107 + (-pkin(8) * t157 - t97) * t198 + t320 * t221;
t31 = t157 * t317 + t43;
t241 = -t200 * t31 + t203 * t24;
t240 = t200 * t24 + t203 * t31;
t239 = qJD(4) * t233 + t179 * t75;
t236 = -t117 ^ 2 - t335;
t235 = t126 * t203 - t127 * t200;
t234 = t126 * t200 + t127 * t203;
t227 = -0.2e1 * pkin(1) * t276 - pkin(7) * qJDD(2);
t226 = t138 * t233 + t74;
t150 = t310 * t196;
t225 = qJD(6) * t150 + t138 * t317 + t281 - t33;
t151 = t310 * t198;
t224 = qJD(4) * t196 - qJD(6) * t151 - t95 - (pkin(8) * t138 - t80) * t198 + t320 * t141;
t223 = g(3) * t328;
t222 = t117 * t278 - t200 * t74 - t203 * t75 + t233 * t279;
t85 = t156 * t157;
t220 = t185 - t294;
t219 = -qJD(5) * t157 * t198 + t81;
t217 = t75 + t336;
t206 = qJD(2) ^ 2;
t215 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t206 + t265;
t207 = qJD(1) ^ 2;
t214 = pkin(1) * t207 - pkin(7) * qJDD(1) + t251;
t213 = -t198 * t75 - t196 * t74 + (t117 * t198 + t196 * t233) * t138;
t210 = t117 * t281 - t251 * t190 - t74 * t297 - t313;
t209 = (-g(1) * (-t186 + t330) - g(2) * t309) * t202;
t169 = qJ(4) * t295;
t167 = t202 * t190 * qJ(4);
t149 = t220 - t311;
t121 = t320 * t198 - t220;
t77 = t128 * t200 + t129 * t203;
t76 = t128 * t203 - t129 * t200;
t55 = t246 * t157 + t111;
t52 = -t246 * t138 + t104;
t49 = -t157 * t322 - t111;
t44 = pkin(4) * t221 - t53;
t42 = qJD(6) * t85 + t328 * t143;
t41 = -qJD(6) * t84 + t156 * t143;
t34 = -t141 * pkin(4) - t46;
t32 = t246 * t143 + t219;
t26 = -pkin(4) * t138 + t255;
t22 = t117 * t320 + t216;
t21 = -t143 * t322 - t219;
t20 = -t140 * pkin(4) - t29;
t13 = t143 * t317 + t17;
t12 = t70 + (-pkin(8) * t143 - t60) * t198 - t320 * t140;
t1 = [qJDD(1), t265, t251, qJDD(1) * t194 + 0.2e1 * t204 * t266, 0.2e1 * t201 * t274 - 0.2e1 * t283 * t276, qJDD(2) * t201 + t204 * t206, qJDD(2) * t204 - t201 * t206, 0, t201 * t227 + t204 * t215, -t201 * t215 + t204 * t227, -t100 * t112 + t101 * t111 - t138 * t82 - t140 * t99 + t141 * t81 - t143 * t98 - t157 * t50 + t221 * t51 - t251, t51 * t112 + t99 * t82 - t50 * t111 - t98 * t81 - t218 * t186 + t164 * t273 - g(1) * (-t186 * t202 + t287) - g(2) * (t202 * t309 + t173) -t10 * t221 + t100 * t53 + t111 * t74 - t117 * t81 + t138 * t29 + t140 * t37 + t196 * t243 + t252, -t100 * t54 + t11 * t221 + t111 * t75 - t138 * t30 - t140 * t38 + t198 * t243 + t233 * t81 + t253, t117 * t30 - t233 * t29 - t53 * t75 - t54 * t74 + t244 * t157 + (-t196 * t38 - t198 * t37) * t143 + t257, -g(1) * t287 - g(2) * t271 + t10 * t53 + t11 * t54 - t111 * t230 - t231 * t81 + t37 * t29 + t38 * t30 + t209, -t100 * t44 - t117 * t32 - t138 * t20 - t140 * t26 + t196 * t247 + t221 * t8 + t55 * t74 + t252, t117 * t17 + t233 * t20 - t43 * t74 + t44 * t75 + t250 * t157 + (-t196 * t27 + t198 * t26) * t143 + t257, t100 * t43 + t138 * t17 + t140 * t27 - t198 * t247 - t221 * t6 - t233 * t32 - t55 * t75 - t253, t6 * t43 + t27 * t17 + t9 * t55 + t36 * t32 + t8 * t44 + t26 * t20 - g(1) * (-pkin(4) * t127 - qJ(5) * t126 + t287) - g(2) * (pkin(4) * t129 + qJ(5) * t128 + t271) + t209, -t222 * t85 + t41 * t59, -t16 * t85 + t222 * t84 - t41 * t57 - t42 * t59, -t140 * t59 - t221 * t222 + t277 * t41 - t85 * t94, t140 * t57 - t16 * t221 - t277 * t42 + t84 * t94, -t140 * t277 - t221 * t94 (t203 * t12 - t200 * t13) * t277 - t241 * t94 + t268 * t221 - t3 * t140 + t21 * t57 + t49 * t16 + t7 * t84 + t22 * t42 + g(1) * t234 - g(2) * t77 + (-t221 * t4 - t240 * t277) * qJD(6) -(t200 * t12 + t203 * t13) * t277 + t240 * t94 - t249 * t221 + t4 * t140 + t21 * t59 - t49 * t222 + t7 * t85 + t22 * t41 + g(1) * t235 - g(2) * t76 + (-t221 * t3 - t241 * t277) * qJD(6); 0, 0, 0, -t201 * t207 * t204, t283 * t207, t275, t274, qJDD(2), t201 * t214 - t312, g(3) * t201 + t204 * t214 (-t104 + t99) * t141 + (t105 - t98) * t138 + (-t100 * t197 - t300 * t101) * pkin(2), t98 * t104 - t99 * t105 + (t300 * t50 - t312 + t197 * t51 + (-qJD(1) * t164 + t251) * t201) * pkin(2), -t179 * t298 + t104 * t117 - t141 * t37 + t185 * t74 + t267 * t198 + (t286 * t196 - t46) * t138 + t284, -t272 - t104 * t233 + t141 * t38 + t185 * t75 + (t286 * t198 + t47) * t138 + t208 * t196, -t117 * t47 + t233 * t46 + (-t138 * t37 + t11) * t198 + (-t138 * t38 - t10 + t239) * t196 + t210, -t230 * t185 - t38 * t47 - t37 * t46 + t231 * t104 - g(1) * (t205 * t254 + t169) - g(2) * (t202 * t254 + t167) - t258 + (-t10 * t196 + t11 * t198) * t179 + t242 * qJD(4), t117 * t52 + t138 * t34 + t141 * t26 + t149 * t74 + t269 * t198 + (qJD(5) * t117 - t100 * t179 - t138 * t285) * t196 + t284, -t117 * t33 - t233 * t34 + (t138 * t26 + t6) * t198 + (t239 + t334) * t196 + t210, t272 + t233 * t52 - t141 * t27 - t149 * t75 + (t198 * t285 - t33) * t138 + (t269 + t280 + t331) * t196, t9 * t149 - t27 * t33 - t36 * t52 - t26 * t34 - g(1) * (-t205 * t318 + t169) - g(2) * (-t202 * t318 + t167) - t258 + (-pkin(4) * t181 + qJD(4) * t27 + t179 * t6) * t198 + (-qJ(5) * t181 + qJD(4) * t26 - qJD(5) * t36 + t179 * t8) * t196 + (pkin(3) + t294 + t311) * t331, t222 * t328 - t302 * t59, t156 * t222 + t16 * t328 - t301 * t59 + t302 * t57, t303 - t325, t256 - t305, t277 * t141 -(t150 * t203 - t151 * t200) * t94 + t121 * t16 + t3 * t141 + t260 * t57 + t301 * t22 - (t200 * t225 - t203 * t224) * t277 + (-t181 + t321) * t156 (t150 * t200 + t151 * t203) * t94 - t121 * t222 - t4 * t141 + t260 * t59 - t302 * t22 + t190 * t223 - (t200 * t224 + t203 * t225) * t277 - t321 * t328; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141 ^ 2 - t137, t138 * t99 + t141 * t98 + t218 - t265, t324, -t323, t213, t138 * t242 + t141 * t231 - t244 - t265, t324, t213, t323, -t141 * t36 + (t196 * t26 + t198 * t27) * t138 - t250 - t265, 0, 0, 0, 0, 0, t256 + t305, t303 + t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226, t217, t236, -t117 * t38 + t233 * t37 + t208, t226, t236, -t217, -t117 * t27 - t233 * t26 + t208 + t319 - t329, 0, 0, 0, 0, 0, -t16 - t337, t222 + t338; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t141 - t117 * t233 - t245, t75 - t336, -t137 - t335, -g(1) * t128 - g(2) * t126 - t196 * t313 + t233 * t36 + t334, 0, 0, 0, 0, 0, -t200 * t259 - t203 * t94 - t233 * t57, t200 * t94 - t203 * t259 - t233 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59 * t57, -t57 ^ 2 + t59 ^ 2, -t222 + t338, -t16 + t337, -t94, -g(1) * t76 - g(2) * t235 + t189 * t223 - t22 * t59 + t327 * t4 + t268, g(1) * t77 + g(2) * t234 + t156 * t313 + t22 * t57 + t327 * t3 - t249;];
tau_reg  = t1;
