% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR13
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR13_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR13_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:46:37
% EndTime: 2019-12-31 21:46:51
% DurationCPUTime: 5.11s
% Computational Cost: add. (4360->479), mult. (10992->652), div. (0->0), fcn. (8574->10), ass. (0->233)
t173 = sin(qJ(3));
t177 = cos(qJ(3));
t304 = cos(pkin(5));
t250 = t304 * qJD(1);
t218 = t250 + qJD(2);
t174 = sin(qJ(2));
t171 = sin(pkin(5));
t285 = qJD(1) * t171;
t266 = t174 * t285;
t100 = t173 * t218 + t177 * t266;
t92 = qJD(5) + t100;
t241 = pkin(1) * t250;
t178 = cos(qJ(2));
t265 = t178 * t285;
t117 = pkin(7) * t265 + t174 * t241;
t83 = pkin(8) * t218 + t117;
t216 = -pkin(2) * t178 - pkin(8) * t174 - pkin(1);
t111 = t216 * t171;
t91 = qJD(1) * t111;
t40 = t173 * t83 - t177 * t91;
t290 = -qJD(4) - t40;
t276 = qJDD(1) * t174;
t259 = t171 * t276;
t277 = qJD(1) * qJD(2);
t260 = t171 * t277;
t339 = t178 * t260 + t259;
t145 = -qJD(3) + t265;
t175 = sin(qJ(1));
t326 = cos(qJ(1));
t234 = t304 * t326;
t123 = t174 * t175 - t178 * t234;
t172 = sin(qJ(5));
t176 = cos(qJ(5));
t124 = t174 * t234 + t175 * t178;
t267 = t171 * t326;
t72 = t124 * t173 + t177 * t267;
t338 = t123 * t176 + t172 * t72;
t337 = -t123 * t172 + t176 * t72;
t203 = qJD(3) * t218;
t247 = t304 * qJDD(1);
t213 = t247 + qJDD(2);
t280 = qJD(3) * t177;
t283 = qJD(2) * t178;
t43 = (qJD(1) * (t173 * t283 + t174 * t280) + t173 * t276) * t171 + t173 * t203 - t177 * t213;
t168 = t171 ^ 2;
t336 = 0.2e1 * t168;
t255 = t178 * t304;
t297 = t171 * t174;
t210 = pkin(1) * t255 - pkin(7) * t297;
t335 = (qJDD(2) + 0.2e1 * t247) * t171;
t125 = t174 * t326 + t175 * t255;
t231 = g(1) * t125 + g(2) * t123;
t294 = t171 * t178;
t197 = -g(3) * t294 + t231;
t243 = t173 * t265;
t282 = qJD(3) * t173;
t334 = -t173 * qJD(4) - t117 + (-t243 + t282) * pkin(3);
t289 = t100 * pkin(4) - t290;
t98 = t173 * t266 - t177 * t218;
t327 = pkin(4) * t98;
t41 = t173 * t91 + t177 * t83;
t34 = qJ(4) * t145 - t41;
t22 = -t34 - t327;
t329 = pkin(3) + pkin(9);
t271 = t173 * t297;
t242 = qJD(3) * t271;
t42 = qJD(1) * t242 - t173 * t213 + (-t203 - t339) * t177;
t39 = -qJDD(5) + t42;
t333 = t329 * t39 + (t22 - t41 + t327) * t92;
t257 = -t173 * qJ(4) - pkin(2);
t127 = -t177 * t329 + t257;
t328 = pkin(4) + pkin(8);
t146 = t328 * t173;
t147 = t328 * t177;
t256 = t174 * t304;
t126 = -t175 * t256 + t178 * t326;
t229 = g(1) * t126 + g(2) * t124;
t114 = -pkin(7) * t266 + t178 * t241;
t102 = t173 * t114;
t232 = pkin(2) * t174 - pkin(8) * t178;
t115 = t232 * t285;
t249 = -t177 * t115 + t102;
t291 = t177 * t178;
t332 = -t146 * t39 - (qJD(5) * t127 - qJD(3) * t147 + (pkin(4) * t291 - t174 * t329) * t285 + t249) * t92 - t229;
t287 = pkin(1) * t256 + pkin(7) * t294;
t110 = pkin(8) * t304 + t287;
t208 = t232 * qJD(2);
t116 = t171 * t208;
t118 = t210 * qJD(2);
t204 = t110 * t282 - t111 * t280 - t116 * t173 - t118 * t177;
t284 = qJD(2) * t174;
t17 = -t171 * (qJ(4) * t284 - qJD(4) * t178) + t204;
t275 = qJDD(1) * t178;
t161 = t171 * t275;
t239 = t174 * t260;
t112 = qJDD(3) - t161 + t239;
t220 = t145 * t176 - t172 * t98;
t16 = -qJD(5) * t220 + t172 * t112 - t176 * t43;
t18 = t145 * t329 + t289;
t82 = -pkin(2) * t218 - t114;
t182 = -t100 * qJ(4) + t82;
t21 = t329 * t98 + t182;
t223 = t172 * t21 - t176 * t18;
t240 = pkin(1) * qJD(2) * t304;
t214 = qJD(1) * t240;
t235 = pkin(1) * t247;
t268 = pkin(7) * t161 + t174 * t235 + t178 * t214;
t195 = -pkin(7) * t239 + t268;
t56 = pkin(8) * t213 + t195;
t59 = (qJD(1) * t208 + qJDD(1) * t216) * t171;
t253 = t173 * t56 - t177 * t59 + t280 * t83 + t282 * t91;
t226 = qJDD(4) + t253;
t4 = -t42 * pkin(4) - t112 * t329 + t226;
t245 = pkin(7) * t339 + t174 * t214 - t178 * t235;
t57 = -pkin(2) * t213 + t245;
t181 = qJ(4) * t42 - qJD(4) * t100 + t57;
t6 = t329 * t43 + t181;
t1 = -qJD(5) * t223 + t172 * t4 + t176 * t6;
t330 = t100 ^ 2;
t180 = qJD(1) ^ 2;
t325 = pkin(3) * t112;
t324 = pkin(8) * t112;
t121 = -t177 * t304 + t271;
t321 = t121 * pkin(3);
t61 = -t145 * t172 - t176 * t98;
t320 = t61 * t92;
t319 = t220 * t92;
t318 = t110 * t177 + t111 * t173;
t317 = pkin(8) * qJD(3);
t316 = t100 * t98;
t315 = t145 * t34;
t314 = t145 * t41;
t313 = t145 * t61;
t312 = t145 * t220;
t311 = t145 * t98;
t278 = qJD(5) * t176;
t279 = qJD(5) * t172;
t15 = t112 * t176 + t145 * t279 + t172 * t43 + t278 * t98;
t310 = t15 * t176;
t309 = t172 * t39;
t308 = t172 * t92;
t36 = t176 * t39;
t307 = t98 * qJ(4);
t303 = qJ(4) * t177;
t306 = -qJ(4) * t280 + t265 * t303 + t334;
t288 = t114 * t177 + t115 * t173;
t293 = t173 * t178;
t305 = -t328 * t282 - (-pkin(4) * t293 + qJ(4) * t174) * t285 - t288;
t302 = t100 * t145;
t108 = t112 * qJ(4);
t295 = t171 * t177;
t122 = t173 * t304 + t174 * t295;
t301 = t122 * qJ(4);
t298 = t168 * t180;
t296 = t171 * t175;
t292 = t176 * t178;
t263 = t171 * t283;
t119 = pkin(7) * t263 + t174 * t240;
t169 = t174 ^ 2;
t286 = -t178 ^ 2 + t169;
t281 = qJD(3) * t176;
t272 = t178 * t298;
t270 = t171 * t292;
t264 = t171 * t284;
t262 = pkin(1) * t336;
t261 = t174 * t277;
t254 = -t173 * t59 - t177 * t56 - t280 * t91 + t282 * t83;
t252 = -t110 * t173 + t111 * t177;
t73 = t124 * t177 - t173 * t267;
t76 = t126 * t173 - t175 * t295;
t237 = g(1) * t72 - g(2) * t76;
t77 = t126 * t177 + t173 * t296;
t236 = -g(1) * t73 + g(2) * t77;
t46 = pkin(3) * t294 - t252;
t233 = t171 * t180 * t304;
t230 = -g(1) * t123 + g(2) * t125;
t8 = t172 * t18 + t176 * t21;
t25 = pkin(4) * t122 + pkin(9) * t294 + t46;
t109 = -pkin(2) * t304 - t210;
t187 = t109 - t301;
t29 = t121 * t329 + t187;
t222 = -t172 * t29 + t176 * t25;
t221 = t172 * t25 + t176 * t29;
t217 = 0.2e1 * t250 + qJD(2);
t45 = qJ(4) * t294 - t318;
t215 = -t110 * t280 - t111 * t282 + t177 * t116 - t118 * t173;
t132 = qJD(4) * t145;
t9 = -t108 + t132 + t254;
t212 = pkin(3) * t177 - t257;
t209 = -t308 * t92 - t36;
t70 = t121 * t176 + t172 * t294;
t207 = t145 * t177;
t205 = -g(1) * t77 - g(2) * t73 - g(3) * t122;
t69 = -t242 + (qJD(3) * t304 + t263) * t177;
t201 = -t69 * qJ(4) - t122 * qJD(4) + t119;
t199 = -t176 * t92 ^ 2 + t309;
t198 = t171 * (t172 * t293 + t174 * t176);
t196 = -g(3) * t297 - t229;
t194 = t127 * t39 + (-qJD(5) * t146 - t334 + t145 * (pkin(9) * t173 - t303)) * t92;
t193 = -t145 * t82 - t324;
t28 = t98 * pkin(3) + t182;
t192 = t145 * t28 + t324;
t190 = g(1) * t76 + g(2) * t72 + g(3) * t121 - t253;
t189 = t205 - t254;
t2 = -qJD(5) * t8 - t172 * t6 + t176 * t4;
t188 = t145 * t317 + t197;
t10 = pkin(3) * t43 + t181;
t186 = -t10 + t188;
t185 = -t42 - t311;
t5 = -pkin(4) * t43 - t9;
t184 = t5 + (t329 * t92 + t307) * t92 + t205;
t183 = t100 * t28 + qJDD(4) - t190;
t90 = qJD(1) * t198;
t89 = t172 * t266 - t176 * t243;
t71 = -t121 * t172 + t270;
t68 = qJD(3) * t122 + t173 * t263;
t58 = pkin(3) * t100 + t307;
t51 = t125 * t176 + t172 * t76;
t50 = -t125 * t172 + t176 * t76;
t49 = -pkin(3) * t266 + t249;
t47 = -qJ(4) * t266 - t288;
t44 = t187 + t321;
t32 = pkin(3) * t145 - t290;
t30 = -pkin(4) * t121 - t45;
t27 = qJD(5) * t70 + t68 * t172 + t176 * t264;
t26 = -qJD(5) * t270 - t68 * t176 + (qJD(5) * t121 + t264) * t172;
t20 = pkin(3) * t68 + t201;
t19 = -pkin(3) * t264 - t215;
t14 = t329 * t68 + t201;
t13 = -t68 * pkin(4) - t17;
t12 = t69 * pkin(4) - t264 * t329 - t215;
t11 = t226 - t325;
t3 = [qJDD(1), g(1) * t175 - g(2) * t326, g(1) * t326 + g(2) * t175, (qJDD(1) * t169 + 0.2e1 * t178 * t261) * t168, (t174 * t275 - t277 * t286) * t336, t174 * t335 + t217 * t263, t178 * t335 - t217 * t264, t213 * t304, -t119 * t218 + t210 * t213 - t245 * t304 + g(1) * t124 - g(2) * t126 + (-t261 + t275) * t262, -t118 * t218 - t287 * t213 - t195 * t304 + (-t178 * t277 - t276) * t262 + t230, t100 * t69 - t122 * t42, -t100 * t68 + t121 * t42 - t122 * t43 - t69 * t98, t122 * t112 - t69 * t145 + (t100 * t284 + t178 * t42) * t171, -t121 * t112 + t68 * t145 + (t178 * t43 - t284 * t98) * t171, (-t112 * t178 - t145 * t284) * t171, -t215 * t145 + t252 * t112 + t119 * t98 + t109 * t43 + t57 * t121 + t82 * t68 + (t178 * t253 - t284 * t40) * t171 - t236, -t204 * t145 - t318 * t112 + t119 * t100 - t109 * t42 + t57 * t122 + t82 * t69 + (-t178 * t254 - t284 * t41) * t171 - t237, t100 * t19 + t11 * t122 + t121 * t9 + t17 * t98 + t32 * t69 + t34 * t68 - t42 * t46 + t43 * t45 - t230, -t10 * t121 + t112 * t46 - t145 * t19 - t20 * t98 - t28 * t68 - t43 * t44 + (-t11 * t178 + t284 * t32) * t171 + t236, -t10 * t122 - t20 * t100 - t45 * t112 + t17 * t145 - t28 * t69 + t44 * t42 + (t178 * t9 - t284 * t34) * t171 + t237, t10 * t44 + t28 * t20 + t9 * t45 + t34 * t17 + t11 * t46 + t32 * t19 - g(1) * (-pkin(1) * t175 - pkin(2) * t124 - pkin(3) * t73 + pkin(7) * t267 - pkin(8) * t123 - qJ(4) * t72) - g(2) * (pkin(1) * t326 + pkin(2) * t126 + pkin(3) * t77 + pkin(7) * t296 + pkin(8) * t125 + qJ(4) * t76), -t15 * t71 - t220 * t27, t15 * t70 + t16 * t71 + t220 * t26 - t27 * t61, t122 * t15 - t220 * t69 + t27 * t92 + t39 * t71, -t122 * t16 - t26 * t92 - t39 * t70 - t61 * t69, -t122 * t39 + t69 * t92, (-qJD(5) * t221 + t176 * t12 - t172 * t14) * t92 - t222 * t39 + t2 * t122 - t223 * t69 + t13 * t61 + t30 * t16 - t5 * t70 + t22 * t26 + g(1) * t338 - g(2) * t51, -(qJD(5) * t222 + t172 * t12 + t176 * t14) * t92 + t221 * t39 - t1 * t122 - t8 * t69 - t13 * t220 + t30 * t15 - t5 * t71 + t22 * t27 + g(1) * t337 - g(2) * t50; 0, 0, 0, -t174 * t272, t286 * t298, -t178 * t233 + t259, t174 * t233 + t161, t213, pkin(1) * t174 * t298 + t117 * t218 + t197 - t245, pkin(1) * t272 + t114 * t218 + (pkin(7) * t277 + g(3)) * t297 + t229 - t268, -t100 * t207 - t42 * t173, (-t42 + t311) * t177 + (-t43 + t302) * t173, -t145 * t280 + t173 * t112 + (-t100 * t174 + t145 * t291) * t285, t145 * t282 + t177 * t112 + (-t145 * t293 + t174 * t98) * t285, t145 * t266, t40 * t266 - pkin(2) * t43 - t102 * t145 - t117 * t98 + t193 * t173 + (-t57 + (t115 + t317) * t145 + t197) * t177, pkin(2) * t42 - t288 * t145 + t41 * t266 - t117 * t100 + t193 * t177 + (-t188 + t57) * t173, -t100 * t49 - t47 * t98 + (-t9 - t145 * t32 + (qJD(3) * t100 - t43) * pkin(8)) * t177 + (t11 - t315 + (qJD(3) * t98 - t42) * pkin(8)) * t173 + t196, t145 * t49 + t173 * t192 - t177 * t186 + t212 * t43 - t266 * t32 - t306 * t98, -t100 * t306 - t145 * t47 + t173 * t186 + t177 * t192 - t212 * t42 + t266 * t34, -t32 * t49 - t34 * t47 + t306 * t28 + (t11 * t173 - t9 * t177 + (t173 * t34 + t177 * t32) * qJD(3) + t196) * pkin(8) + (-t10 + t197) * t212, -t15 * t172 * t177 - (t172 * t282 - t177 * t278 - t90) * t220, t90 * t61 - t220 * t89 + (-t172 * t61 - t176 * t220) * t282 + (-t310 + t16 * t172 + (-t172 * t220 + t176 * t61) * qJD(5)) * t177, -t90 * t92 + (qJD(3) * t308 + t15) * t173 + (-t278 * t92 + t309 + t312) * t177, t89 * t92 + (t281 * t92 - t16) * t173 + (t279 * t92 + t313 + t36) * t177, -t39 * t173 - t207 * t92, t147 * t16 - t22 * t89 + t305 * t61 + t194 * t172 + t332 * t176 + (t172 * t231 - t22 * t281 + t2) * t173 - g(3) * t198 + (t145 * t223 + t5 * t176 - t22 * t279) * t177, -t1 * t173 + t147 * t15 - t22 * t90 - t305 * t220 + (t173 * t231 + t194) * t176 + (t22 * t282 - t332) * t172 - g(3) * (-t172 * t174 + t173 * t292) * t171 + (t145 * t8 - t5 * t172 - t22 * t278) * t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t316, -t98 ^ 2 + t330, t185, -t302 - t43, t112, -t100 * t82 + t190 - t314, t145 * t40 + t82 * t98 - t189, pkin(3) * t42 - qJ(4) * t43 + (-t34 - t41) * t100 + (t32 + t290) * t98, t58 * t98 + t183 + t314 - 0.2e1 * t325, t100 * t58 + t145 * t290 - t28 * t98 + 0.2e1 * t108 - t132 + t189, -t9 * qJ(4) - t11 * pkin(3) - t28 * t58 - t32 * t41 - g(1) * (-pkin(3) * t76 + qJ(4) * t77) - g(2) * (-pkin(3) * t72 + qJ(4) * t73) - g(3) * (t301 - t321) + t290 * t34, t220 * t308 + t310, (-t16 + t319) * t176 + (-t15 + t320) * t172, -t220 * t98 + t209, -t61 * t98 + t199, t92 * t98, qJ(4) * t16 + t172 * t184 + t176 * t333 - t223 * t98 + t289 * t61, qJ(4) * t15 - t172 * t333 + t176 * t184 - t220 * t289 - t8 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, t112 - t316, -t145 ^ 2 - t330, t183 - t315 - t325, 0, 0, 0, 0, 0, t209 + t313, t199 - t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220 * t61, t220 ^ 2 - t61 ^ 2, t15 + t320, -t16 - t319, -t39, -g(1) * t50 - g(2) * t337 - g(3) * t70 + t22 * t220 + t8 * t92 + t2, g(1) * t51 + g(2) * t338 - g(3) * t71 + t22 * t61 - t223 * t92 - t1;];
tau_reg = t3;
