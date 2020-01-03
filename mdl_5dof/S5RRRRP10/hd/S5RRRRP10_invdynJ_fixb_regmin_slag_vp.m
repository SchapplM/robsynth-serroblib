% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRP10
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
% tau_reg [5x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:12:24
% EndTime: 2019-12-31 22:12:41
% DurationCPUTime: 6.64s
% Computational Cost: add. (5436->492), mult. (13644->683), div. (0->0), fcn. (10669->10), ass. (0->225)
t191 = cos(qJ(2));
t182 = sin(pkin(5));
t280 = qJD(1) * t182;
t170 = t191 * t280;
t220 = t170 - qJD(3);
t187 = sin(qJ(2));
t183 = cos(pkin(5));
t279 = qJD(1) * t183;
t262 = pkin(1) * t279;
t126 = pkin(7) * t170 + t187 * t262;
t186 = sin(qJ(3));
t190 = cos(qJ(3));
t341 = -t126 - t220 * (pkin(3) * t186 - pkin(9) * t190);
t322 = cos(qJ(1));
t253 = t322 * t191;
t188 = sin(qJ(1));
t291 = t187 * t188;
t133 = -t183 * t253 + t291;
t254 = t322 * t187;
t290 = t188 * t191;
t135 = t183 * t290 + t254;
t295 = t182 * t191;
t203 = g(1) * t135 + g(2) * t133 - g(3) * t295;
t269 = qJDD(1) * t183;
t171 = qJDD(2) + t269;
t270 = qJD(1) * qJD(2);
t246 = t191 * t270;
t226 = t182 * t246;
t231 = qJD(2) * t262;
t268 = qJDD(1) * t187;
t245 = t182 * t268;
t260 = pkin(1) * t269;
t230 = t187 * t231 - t191 * t260 + (t226 + t245) * pkin(7);
t68 = -t171 * pkin(2) + t230;
t340 = t68 - t203;
t134 = t183 * t254 + t290;
t136 = -t183 * t291 + t253;
t222 = g(1) * t136 + g(2) * t134;
t298 = t182 * t187;
t202 = -g(3) * t298 - t222;
t252 = t187 * t280;
t123 = -pkin(7) * t252 + t191 * t262;
t114 = t186 * t123;
t224 = pkin(2) * t187 - pkin(8) * t191;
t124 = t224 * t280;
t57 = -pkin(3) * t252 - t124 * t190 + t114;
t234 = qJD(2) + t279;
t111 = t186 * t252 - t190 * t234;
t184 = -qJ(5) - pkin(9);
t339 = -qJ(5) * t111 + qJD(4) * t184;
t185 = sin(qJ(4));
t189 = cos(qJ(4));
t255 = t182 * t322;
t88 = t134 * t190 - t186 * t255;
t338 = -t133 * t189 + t185 * t88;
t337 = t133 * t185 + t189 * t88;
t217 = qJD(3) * t234;
t277 = qJD(2) * t191;
t250 = t186 * t277;
t275 = qJD(3) * t190;
t51 = -t190 * t171 + (qJD(1) * (t187 * t275 + t250) + t186 * t268) * t182 + t186 * t217;
t179 = t182 ^ 2;
t264 = 0.2e1 * t179;
t247 = t187 * t270;
t227 = t182 * t247;
t266 = qJDD(1) * t191;
t169 = t182 * t266;
t265 = qJDD(3) - t169;
t205 = t227 + t265;
t335 = pkin(8) * t205;
t276 = qJD(3) * t186;
t261 = pkin(8) * t276;
t333 = t185 * t261 + t341 * t189;
t113 = t186 * t234 + t190 * t252;
t332 = t220 * t113;
t249 = t190 * t277;
t267 = qJDD(1) * t190;
t293 = t186 * t171;
t194 = (t187 * t267 + (-t187 * t276 + t249) * qJD(1)) * t182 + t190 * t217 + t293;
t228 = t190 * t170;
t331 = t228 - t275;
t150 = -pkin(3) * t190 - pkin(9) * t186 - pkin(2);
t272 = qJD(4) * t189;
t286 = t190 * t123 + t186 * t124;
t58 = pkin(9) * t252 + t286;
t330 = -t150 * t272 - t341 * t185 + t189 * t58;
t297 = t182 * t188;
t92 = t136 * t190 + t186 * t297;
t59 = t135 * t189 - t185 * t92;
t296 = t182 * t190;
t132 = t183 * t186 + t187 * t296;
t288 = t189 * t191;
t85 = t132 * t185 + t182 * t288;
t328 = -g(1) * t59 + g(2) * t338 + g(3) * t85;
t77 = t189 * t113 - t185 * t220;
t26 = qJD(4) * t77 + t194 * t185 - t189 * t205;
t326 = t77 ^ 2;
t192 = qJD(1) ^ 2;
t105 = qJD(4) + t111;
t97 = -t234 * pkin(2) - t123;
t39 = t111 * pkin(3) - t113 * pkin(9) + t97;
t219 = -pkin(2) * t191 - pkin(8) * t187 - pkin(1);
t122 = t219 * t182;
t104 = qJD(1) * t122;
t98 = t234 * pkin(8) + t126;
t50 = t186 * t104 + t190 * t98;
t42 = -t220 * pkin(9) + t50;
t19 = -t185 * t42 + t189 * t39;
t9 = -qJ(5) * t77 + t19;
t8 = pkin(4) * t105 + t9;
t325 = -t9 + t8;
t321 = pkin(1) * t187;
t320 = pkin(4) * t185;
t103 = (t185 * t187 + t190 * t288) * t280;
t289 = t189 * t190;
t176 = pkin(8) * t289;
t229 = t186 * t170;
t271 = qJD(5) * t189;
t302 = qJ(5) * t186;
t315 = -pkin(4) * t229 + qJ(5) * t103 + t185 * t58 - t186 * t271 + (pkin(4) * t186 - qJ(5) * t289) * qJD(3) + (-t176 + (-t150 + t302) * t185) * qJD(4) + t333;
t102 = t185 * t228 - t189 * t252;
t292 = t186 * t189;
t314 = qJ(5) * t102 + (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t292 + (-qJD(5) * t186 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t190) * t185 - t330;
t49 = t190 * t104 - t186 * t98;
t69 = pkin(3) * t113 + pkin(9) * t111;
t313 = t185 * t69 + t189 * t49;
t172 = pkin(7) * t298;
t120 = t172 + (-pkin(1) * t191 - pkin(2)) * t183;
t131 = -t183 * t190 + t186 * t298;
t54 = pkin(3) * t131 - pkin(9) * t132 + t120;
t283 = pkin(7) * t295 + t183 * t321;
t121 = pkin(8) * t183 + t283;
t287 = t190 * t121 + t186 * t122;
t56 = -pkin(9) * t295 + t287;
t312 = t185 * t54 + t189 * t56;
t238 = t189 * t220;
t75 = t113 * t185 + t238;
t310 = t105 * t75;
t273 = qJD(4) * t185;
t25 = qJD(4) * t238 + t113 * t273 - t185 * t205 - t189 * t194;
t309 = t185 * t25;
t46 = qJDD(4) + t51;
t308 = t185 * t46;
t307 = t189 * t46;
t306 = t77 * t105;
t305 = t339 * t185 + t271 - t313;
t64 = t189 * t69;
t304 = -pkin(4) * t113 - t64 + t339 * t189 + (-qJD(5) + t49) * t185;
t236 = t105 * t189;
t299 = t179 * t192;
t294 = t183 * t191;
t282 = t185 * t150 + t176;
t180 = t187 ^ 2;
t281 = -t191 ^ 2 + t180;
t278 = qJD(2) * t187;
t274 = qJD(4) * t105;
t259 = t191 * t299;
t258 = t185 * t295;
t257 = -pkin(7) * t169 - t187 * t260 - t191 * t231;
t256 = pkin(8) + t320;
t251 = t182 * t278;
t248 = t182 * t183 * t192;
t243 = -t185 * t56 + t189 * t54;
t200 = -pkin(7) * t227 - t257;
t67 = pkin(8) * t171 + t200;
t211 = t224 * qJD(2);
t70 = (qJD(1) * t211 + qJDD(1) * t219) * t182;
t241 = t104 * t276 + t186 * t67 - t190 * t70 + t98 * t275;
t240 = -t186 * t121 + t122 * t190;
t237 = t191 * t220;
t235 = qJD(3) * t220;
t233 = qJD(2) + 0.2e1 * t279;
t232 = t171 + t269;
t87 = t134 * t186 + t190 * t255;
t91 = t136 * t186 - t188 * t296;
t225 = -g(1) * t87 + g(2) * t91;
t55 = pkin(3) * t295 - t240;
t221 = t189 * t275 - t103;
t20 = t185 * t39 + t189 * t42;
t125 = t182 * t211;
t127 = (pkin(1) * t294 - t172) * qJD(2);
t218 = -t121 * t275 - t122 * t276 + t125 * t190 - t186 * t127;
t215 = -t104 * t275 - t186 * t70 - t190 * t67 + t98 * t276;
t12 = pkin(9) * t205 - t215;
t16 = t51 * pkin(3) - pkin(9) * t194 + t68;
t214 = -t189 * t12 - t185 * t16 - t39 * t272 + t42 * t273;
t206 = -t121 * t276 + t122 * t275 + t186 * t125 + t190 * t127;
t29 = pkin(9) * t251 + t206;
t128 = t283 * qJD(2);
t83 = qJD(3) * t132 + t182 * t250;
t84 = -qJD(3) * t131 + t182 * t249;
t33 = pkin(3) * t83 - pkin(9) * t84 + t128;
t213 = t185 * t33 + t189 * t29 + t54 * t272 - t56 * t273;
t41 = t220 * pkin(3) - t49;
t210 = -pkin(9) * t46 + t105 * t41;
t209 = g(1) * t91 + g(2) * t87 + g(3) * t131;
t208 = -g(1) * t92 - g(2) * t88 - g(3) * t132;
t207 = t246 + t268;
t30 = -pkin(3) * t251 - t218;
t4 = -qJD(4) * t20 - t185 * t12 + t189 * t16;
t198 = -t312 * qJD(4) - t185 * t29 + t189 * t33;
t13 = -pkin(3) * t205 + t241;
t197 = pkin(9) * t274 + t13 - t209;
t196 = pkin(8) * t274 - t203;
t7 = t26 * pkin(4) + qJDD(5) + t13;
t178 = pkin(4) * t189 + pkin(3);
t155 = t184 * t189;
t154 = t184 * t185;
t142 = t189 * t150;
t95 = -t185 * t302 + t282;
t86 = t132 * t189 - t258;
t79 = -qJ(5) * t292 + t142 + (-pkin(8) * t185 - pkin(4)) * t190;
t74 = t75 ^ 2;
t60 = t135 * t185 + t189 * t92;
t38 = -qJD(4) * t85 + t185 * t251 + t189 * t84;
t37 = -qJD(4) * t258 + t132 * t272 + t185 * t84 - t189 * t251;
t27 = t75 * pkin(4) + qJD(5) + t41;
t21 = -qJ(5) * t85 + t312;
t17 = pkin(4) * t131 - qJ(5) * t86 + t243;
t10 = -qJ(5) * t75 + t20;
t6 = -qJ(5) * t37 - qJD(5) * t85 + t213;
t5 = pkin(4) * t83 - qJ(5) * t38 - qJD(5) * t86 + t198;
t2 = -qJ(5) * t26 - qJD(5) * t75 - t214;
t1 = pkin(4) * t46 + qJ(5) * t25 - qJD(5) * t77 + t4;
t3 = [qJDD(1), g(1) * t188 - g(2) * t322, g(1) * t322 + g(2) * t188, (qJDD(1) * t180 + 0.2e1 * t187 * t246) * t179, (t187 * t266 - t281 * t270) * t264, (t187 * t232 + t233 * t277) * t182, (t191 * t232 - t233 * t278) * t182, t171 * t183, -t128 * t234 - t172 * t171 - t230 * t183 + g(1) * t134 - g(2) * t136 + (t171 * t294 + (-t247 + t266) * t264) * pkin(1), -t207 * pkin(1) * t264 - g(1) * t133 + g(2) * t135 - t127 * t234 - t283 * t171 - t200 * t183, t113 * t84 + t194 * t132, -t84 * t111 - t113 * t83 - t194 * t131 - t132 * t51, -t84 * t220 + t132 * t265 + ((-t293 + (-t217 - t226) * t190) * t191 + (-(-qJD(1) * t276 + t267) * t295 + (qJD(1) * t132 + t113) * qJD(2)) * t187) * t182, t83 * t220 - t131 * t265 + (t51 * t191 + (-qJD(1) * t131 - t111) * t278) * t182, (-t265 * t191 + (-t170 - t220) * t278) * t182, g(1) * t88 - g(2) * t92 + t128 * t111 + t120 * t51 + t68 * t131 + t240 * t205 - t218 * t220 + t241 * t295 + t49 * t251 + t97 * t83, t128 * t113 + t120 * t194 + t68 * t132 - t287 * t205 + t206 * t220 - t215 * t295 - t50 * t251 + t97 * t84 + t225, -t25 * t86 + t38 * t77, t25 * t85 - t26 * t86 - t37 * t77 - t38 * t75, t105 * t38 - t131 * t25 + t46 * t86 + t77 * t83, -t105 * t37 - t131 * t26 - t46 * t85 - t75 * t83, t105 * t83 + t131 * t46, g(1) * t337 - g(2) * t60 + t198 * t105 + t13 * t85 + t4 * t131 + t19 * t83 + t243 * t46 + t55 * t26 + t30 * t75 + t41 * t37, -g(1) * t338 - g(2) * t59 - t213 * t105 + t13 * t86 + t214 * t131 - t20 * t83 - t55 * t25 + t30 * t77 - t312 * t46 + t41 * t38, -t1 * t86 - t10 * t37 + t17 * t25 - t2 * t85 - t21 * t26 - t38 * t8 - t5 * t77 - t6 * t75 - t225, t2 * t21 + t10 * t6 + t1 * t17 + t8 * t5 + t7 * (pkin(4) * t85 + t55) + t27 * (pkin(4) * t37 + t30) - g(1) * (-t188 * pkin(1) - t134 * pkin(2) + pkin(7) * t255 - t133 * t256 - t178 * t88 + t184 * t87) - g(2) * (t322 * pkin(1) + t136 * pkin(2) + pkin(7) * t297 + t256 * t135 + t92 * t178 - t91 * t184); 0, 0, 0, -t187 * t259, t281 * t299, -t191 * t248 + t245, t187 * t248 + t169, t171, t126 * t234 + t299 * t321 + t203 - t230, pkin(1) * t259 + t123 * t234 + (pkin(7) * t270 + g(3)) * t298 + t222 + t257, (-qJD(3) * t252 + t171) * t186 ^ 2 + ((t182 * t207 + t217) * t186 - t332) * t190, -t186 * t51 + t190 * t194 + (t229 - t276) * t113 + t331 * t111, -t190 * t235 + t186 * t265 + (t190 * t237 + (qJD(2) * t186 - t113) * t187) * t280, t186 * t235 + t190 * t265 + (-t186 * t237 + (qJD(2) * t190 + t111) * t187) * t280, t220 * t252, -pkin(2) * t51 - t114 * t220 - t49 * t252 - t126 * t111 + (-t220 * t97 - t335) * t186 + (pkin(8) * t235 + t124 * t220 - t340) * t190, -t190 * t335 - pkin(2) * t194 - t126 * t113 + t50 * t252 - t331 * t97 + (-t261 - t286) * t220 + t340 * t186, -t25 * t292 + (-t186 * t273 + t221) * t77, t102 * t77 + t103 * t75 + (-t185 * t77 - t189 * t75) * t275 + (t309 - t189 * t26 + (t185 * t75 - t189 * t77) * qJD(4)) * t186, t190 * t25 + t221 * t105 + (-t105 * t273 - t220 * t77 + t307) * t186, t190 * t26 + (-t185 * t275 + t102) * t105 + (-t105 * t272 + t220 * t75 - t308) * t186, -t105 * t186 * t220 - t190 * t46, -t41 * t102 + t142 * t46 - t57 * t75 + t333 * t105 + ((-qJD(4) * t150 + t58) * t105 + t202) * t185 + (t41 * t185 * qJD(3) - t4 + (qJD(3) * t75 - t308) * pkin(8) - t196 * t189) * t190 + (pkin(8) * t26 + t13 * t185 - t19 * t220 + t41 * t272) * t186, -t282 * t46 - t57 * t77 - t41 * t103 + t330 * t105 + t202 * t189 + (-t214 + (pkin(8) * t77 + t41 * t189) * qJD(3) + t196 * t185) * t190 + (-t41 * t273 + t13 * t189 + t220 * t20 + (qJD(3) * t236 - t25) * pkin(8)) * t186, t10 * t102 + t103 * t8 + t25 * t79 - t26 * t95 - t315 * t77 - t314 * t75 + (-t10 * t185 - t189 * t8) * t275 + (-t1 * t189 - t185 * t2 + (-t10 * t189 + t185 * t8) * qJD(4) + t203) * t186, t1 * t79 + t314 * t10 + t2 * t95 + t315 * t8 + ((t272 * t186 - t102) * pkin(4) - t57) * t27 + (t7 * t186 + t27 * t275 + t202) * t256 + t203 * (t178 * t190 - t184 * t186 + pkin(2)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113 * t111, -t111 ^ 2 + t113 ^ 2, -t111 * t220 + t194, -t332 - t51, t205, -t97 * t113 - t50 * t220 + t209 - t241, t97 * t111 - t49 * t220 - t208 + t215, t236 * t77 - t309, (-t25 - t310) * t189 + (-t26 - t306) * t185, t105 * t236 - t113 * t77 + t308, -t105 ^ 2 * t185 + t113 * t75 + t307, -t105 * t113, -pkin(3) * t26 - t64 * t105 - t19 * t113 - t50 * t75 + (t49 * t105 + t210) * t185 - t197 * t189, pkin(3) * t25 + t313 * t105 + t20 * t113 + t197 * t185 + t210 * t189 - t50 * t77, t154 * t25 + t155 * t26 - t304 * t77 - t305 * t75 + (-t105 * t8 + t2) * t189 + (-t10 * t105 - t1) * t185 + t208, -t2 * t155 + t1 * t154 - t7 * t178 - g(1) * (-t178 * t91 - t184 * t92) - g(2) * (-t178 * t87 - t184 * t88) - g(3) * (-t131 * t178 - t132 * t184) + t304 * t8 + (t105 * t320 - t50) * t27 + t305 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77 * t75, -t74 + t326, -t25 + t310, -t26 + t306, t46, t20 * t105 - t41 * t77 + t328 + t4, g(1) * t60 + g(2) * t337 + g(3) * t86 + t19 * t105 + t41 * t75 + t214, pkin(4) * t25 - t325 * t75, t325 * t10 + (-t27 * t77 + t1 + t328) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74 - t326, t10 * t75 + t8 * t77 - t209 + t7;];
tau_reg = t3;
