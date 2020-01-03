% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRR16
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRR16_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR16_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:47:28
% EndTime: 2019-12-31 20:47:39
% DurationCPUTime: 4.02s
% Computational Cost: add. (16278->401), mult. (37287->549), div. (0->0), fcn. (26910->10), ass. (0->251)
t206 = sin(qJ(5));
t205 = cos(pkin(5));
t199 = t205 * qJD(1) + qJD(2);
t207 = sin(qJ(4));
t211 = cos(qJ(4));
t204 = sin(pkin(5));
t212 = cos(qJ(2));
t265 = qJD(1) * t212;
t251 = t204 * t265;
t159 = t211 * t199 - t207 * t251;
t208 = sin(qJ(2));
t266 = qJD(1) * t208;
t250 = t204 * t266;
t187 = qJD(4) + t250;
t210 = cos(qJ(5));
t138 = t206 * t159 - t210 * t187;
t140 = t210 * t159 + t206 * t187;
t113 = t140 * t138;
t190 = qJD(2) * t250;
t261 = t212 * qJDD(1);
t171 = t204 * t261 - t190;
t198 = t205 * qJDD(1) + qJDD(2);
t244 = t211 * t171 + t207 * t198;
t125 = -t159 * qJD(4) - t244;
t124 = qJDD(5) - t125;
t310 = -t113 + t124;
t319 = t206 * t310;
t157 = t207 * t199 + t211 * t251;
t134 = t159 * t157;
t170 = (qJD(2) * t265 + qJDD(1) * t208) * t204;
t162 = qJDD(4) + t170;
t309 = -t134 + t162;
t318 = t207 * t309;
t317 = t210 * t310;
t316 = t211 * t309;
t299 = pkin(3) + pkin(7);
t195 = t199 ^ 2;
t202 = t204 ^ 2;
t214 = qJD(1) ^ 2;
t279 = t202 * t214;
t197 = t212 ^ 2 * t279;
t174 = -t197 - t195;
t257 = t212 * t279;
t184 = t208 * t257;
t306 = -t198 - t184;
t315 = pkin(1) * (t174 * t208 - t212 * t306);
t196 = t208 ^ 2 * t279;
t152 = -t196 - t195;
t169 = -t184 + t198;
t314 = pkin(1) * (t152 * t212 - t169 * t208);
t271 = t208 * t306;
t313 = pkin(7) * (-t212 * t174 - t271);
t269 = t212 * t169;
t312 = pkin(7) * (t208 * t152 + t269);
t132 = t157 * pkin(4) - t159 * pkin(9);
t301 = t187 ^ 2;
t209 = sin(qJ(1));
t213 = cos(qJ(1));
t247 = t209 * g(1) - t213 * g(2);
t293 = t204 * pkin(7);
t164 = qJDD(1) * pkin(1) + t214 * t293 + t247;
t167 = pkin(3) * t250 - t199 * pkin(8);
t249 = qJD(3) * t266;
t188 = -0.2e1 * t204 * t249;
t177 = t199 * pkin(2) * t250;
t292 = t205 * g(3);
t223 = -t170 * qJ(3) + t177 - t292;
t280 = t199 * t212;
t258 = qJ(3) * t280;
t300 = -pkin(2) - pkin(8);
t78 = -pkin(3) * t197 + t188 + t300 * t171 + (-t164 + (-t167 * t208 - t258) * qJD(1)) * t204 + t223;
t286 = t211 * t78;
t294 = t198 * pkin(8);
t273 = t208 * qJ(3);
t238 = -t212 * pkin(2) - t273;
t267 = qJD(1) * t204;
t166 = t238 * t267;
t283 = t164 * t205;
t150 = t212 * t283;
t304 = -t198 * pkin(2) - t195 * qJ(3) + qJDD(3);
t242 = -t150 + t304;
t240 = t213 * g(1) + t209 * g(2);
t165 = -t214 * pkin(1) + qJDD(1) * t293 - t240;
t272 = t208 * t165;
t220 = t242 + t272;
t296 = g(3) * t212;
t305 = (t296 + (-pkin(3) * t280 + t166 * t208) * qJD(1)) * t204 + t170 * pkin(3) + t220;
t47 = t286 + (-pkin(8) * t184 - t294 + t305) * t207;
t34 = -t301 * pkin(4) + t162 * pkin(9) - t157 * t132 + t47;
t228 = -t207 * t171 + t211 * t198;
t126 = -t157 * qJD(4) + t228;
t147 = t187 * t157;
t108 = t126 - t147;
t278 = t204 * t208;
t243 = -g(3) * t278 + t212 * t165;
t219 = -t195 * pkin(2) + t198 * qJ(3) + 0.2e1 * qJD(3) * t199 + t166 * t251 + t243;
t256 = t208 * t283;
t95 = t219 + t256;
t77 = t171 * pkin(3) - pkin(8) * t197 + t199 * t167 + t95;
t45 = -t108 * pkin(9) + (t187 * t159 - t125) * pkin(4) + t77;
t17 = t206 * t34 - t210 * t45;
t18 = t206 * t45 + t210 * t34;
t10 = t206 * t17 + t210 * t18;
t311 = t202 * (qJD(1) * t199 - t205 * t214);
t179 = t199 * t251;
t308 = t170 + t179;
t307 = t179 - t170;
t153 = qJD(5) + t157;
t245 = t206 * t126 - t210 * t162;
t68 = (qJD(5) - t153) * t140 + t245;
t252 = t199 * t266;
t145 = -t190 + (t252 + t261) * t204;
t303 = ((t197 - t195) * t208 + t269) * t204 + t205 * t145;
t302 = ((-t196 + t195) * t212 - t271) * t204 - t205 * t307;
t136 = t138 ^ 2;
t137 = t140 ^ 2;
t151 = t153 ^ 2;
t155 = t157 ^ 2;
t156 = t159 ^ 2;
t46 = t207 * t78 - t211 * (t306 * pkin(8) + t305);
t33 = -t162 * pkin(4) - t301 * pkin(9) + t159 * t132 + t46;
t298 = -pkin(4) * t33 + pkin(9) * t10;
t297 = pkin(4) * t211;
t295 = t171 * pkin(2);
t30 = t206 * t33;
t83 = t113 + t124;
t291 = t206 * t83;
t290 = t207 * t77;
t109 = t126 + t147;
t218 = (-qJD(4) + t187) * t159 - t244;
t74 = -t211 * t109 + t207 * t218;
t289 = t208 * t74;
t31 = t210 * t33;
t288 = t210 * t83;
t287 = t211 * t77;
t285 = t153 * t206;
t284 = t153 * t210;
t282 = t187 * t207;
t281 = t187 * t211;
t277 = t204 * t212;
t121 = t134 + t162;
t274 = t207 * t121;
t270 = t211 * t121;
t176 = -t196 - t197;
t268 = pkin(1) * (-t204 * t176 + (t145 * t208 + t212 * t307) * t205) + (t212 * t145 - t208 * t307) * t293;
t264 = qJD(4) + t187;
t262 = qJD(5) + t153;
t111 = -t137 - t151;
t55 = -t206 * t111 - t288;
t235 = -t210 * t126 - t206 * t162;
t73 = t262 * t138 + t235;
t260 = pkin(4) * t73 + pkin(9) * t55 + t30;
t99 = -t151 - t136;
t52 = t210 * t99 - t319;
t69 = -t262 * t140 - t245;
t259 = pkin(4) * t69 + pkin(9) * t52 - t31;
t255 = t207 * t113;
t254 = t208 * t134;
t253 = t211 * t113;
t248 = pkin(4) * t207 + qJ(3);
t123 = t153 * t138;
t90 = -t138 * qJD(5) - t235;
t72 = t123 + t90;
t43 = t206 * t72 - t210 * t68;
t92 = t136 + t137;
t246 = pkin(4) * t92 + pkin(9) * t43 + t10;
t241 = g(3) * t277 - t150;
t239 = t166 * t267 + t165;
t148 = t204 * t164 + t292;
t9 = -t210 * t17 + t206 * t18;
t21 = t207 * t47 - t211 * t46;
t237 = t207 * t46 + t211 * t47;
t226 = pkin(1) - t238;
t225 = -qJD(1) * t258 - t164;
t130 = t241 + t272;
t131 = t243 + t256;
t222 = (t208 * t130 + t212 * t131) * t204;
t221 = t212 * t300 - pkin(1) - t273;
t217 = -t223 + t295;
t101 = (t166 * t266 + t296) * t204 + t220;
t182 = t205 * t198;
t175 = t196 - t197;
t146 = -t190 + (-t252 + t261) * t204;
t142 = -t156 + t301;
t141 = t155 - t301;
t135 = -t156 - t301;
t133 = t156 - t155;
t129 = -t301 - t155;
t128 = (t170 * t204 + t212 * t311) * t208;
t127 = (t171 * t204 - t208 * t311) * t212;
t119 = -t137 + t151;
t118 = t136 - t151;
t117 = -t155 - t156;
t116 = pkin(2) * t307 + qJ(3) * t145;
t115 = (-t157 * t211 + t159 * t207) * t187;
t112 = t137 - t136;
t110 = -t264 * t157 + t228;
t105 = t264 * t159 + t244;
t104 = t211 * t126 - t159 * t282;
t103 = -t207 * t125 + t157 * t281;
t102 = t205 * t175 + (t208 * t146 + t308 * t212) * t204;
t98 = t211 * t141 - t274;
t97 = -t207 * t142 + t316;
t94 = -t207 * t135 - t270;
t93 = t211 * t135 - t274;
t89 = -t140 * qJD(5) - t245;
t87 = t207 * t129 + t316;
t86 = (-t138 * t210 + t140 * t206) * t153;
t85 = (-t138 * t206 - t140 * t210) * t153;
t81 = pkin(2) * t306 - qJ(3) * t174 + t101;
t80 = -pkin(2) * t152 + qJ(3) * t169 + t95;
t75 = -t211 * t105 - t207 * t108;
t71 = -t123 + t90;
t65 = -t140 * t285 + t210 * t90;
t64 = t140 * t284 + t206 * t90;
t63 = t138 * t284 - t206 * t89;
t62 = t138 * t285 + t210 * t89;
t61 = t207 * t124 + t211 * t86;
t60 = t210 * t118 - t291;
t59 = -t206 * t119 + t317;
t58 = t206 * t118 + t288;
t57 = t210 * t119 + t319;
t56 = -pkin(2) * t101 + qJ(3) * t95;
t54 = t210 * t111 - t291;
t51 = t206 * t99 + t317;
t49 = t211 * t65 + t255;
t48 = t211 * t63 - t255;
t42 = -t206 * t71 + t210 * t69;
t41 = -t206 * t68 - t210 * t72;
t40 = t206 * t69 + t210 * t71;
t38 = -t207 * t68 + t211 * t60;
t37 = t207 * t72 + t211 * t59;
t35 = t207 * t55 + t211 * t73;
t28 = t207 * t52 + t211 * t69;
t27 = t207 * t112 + t211 * t42;
t26 = qJ(3) * t110 + t300 * t93 + t287;
t24 = t207 * t43 + t211 * t92;
t23 = qJ(3) * t105 + t300 * t87 + t290;
t20 = -pkin(9) * t54 + t31;
t19 = -pkin(9) * t51 + t30;
t14 = -pkin(4) * t54 + t18;
t13 = -pkin(4) * t51 + t17;
t12 = qJ(3) * t117 + t300 * t74 - t21;
t11 = qJ(3) * t77 + t300 * t21;
t7 = -pkin(9) * t41 - t9;
t5 = t207 * t10 - t211 * t33;
t4 = qJ(3) * t54 - t207 * t14 + t211 * t20 + t300 * t35;
t3 = qJ(3) * t51 - t207 * t13 + t211 * t19 + t300 * t28;
t2 = t211 * t7 + t300 * t24 + t248 * t41;
t1 = t300 * t5 + (-pkin(9) * t211 + t248) * t9;
t6 = [0, 0, 0, 0, 0, qJDD(1), t247, t240, 0, 0, t128, t102, t302, t127, t303, t182, (-t130 + t315) * t205 + (pkin(1) * t146 + t212 * t148 - t313) * t204, (-t131 + t314) * t205 + (-pkin(1) * t308 - t208 * t148 - t312) * t204, t222 + t268, pkin(1) * (t204 * t148 + (-t130 * t212 + t131 * t208) * t205) + pkin(7) * t222, t182, -t302, -t303, t128, t102, t127, t205 * t116 + (t212 * (-pkin(2) * t176 + t219) + (-qJ(3) * t176 + t239 * t208 + (g(3) * t204 + t283) * t212 + t242) * t208) * t204 + t268, (t81 - t315) * t205 + (t212 * (t188 - t217) + t313 + t225 * t277 - t226 * t146) * t204, (t80 - t314) * t205 + (t208 * t217 + t312 + (-t225 + 0.2e1 * t249) * t278 + t226 * t308) * t204, (t56 + pkin(1) * (-t101 * t212 + t208 * t95)) * t205 + (pkin(7) * (t208 * t101 + t212 * t95) - t226 * (-qJ(3) * t308 - t148 + t177 + t188 - t295)) * t204, t205 * t104 + (t254 + t212 * (-t207 * t126 - t159 * t281)) * t204, t205 * t75 + (t208 * t133 + t212 * (t207 * t105 - t211 * t108)) * t204, t205 * t97 + (t208 * t109 + t212 * (-t211 * t142 - t318)) * t204, t205 * t103 + (-t254 + t212 * (-t211 * t125 - t157 * t282)) * t204, t205 * t98 + (t208 * t218 + t212 * (-t207 * t141 - t270)) * t204, t162 * t278 + t205 * t115 + (t157 * t207 + t159 * t211) * t187 * t277, (t23 + pkin(1) * (t105 * t208 - t212 * t87)) * t205 + (t208 * (pkin(3) * t87 - t46) + t212 * (pkin(3) * t105 + t287) + pkin(7) * (t212 * t105 + t208 * t87) + t221 * (t211 * t129 - t318)) * t204, (pkin(3) * t93 - qJ(3) * t94 - t286 + (-t241 + t294 - (-pkin(8) * t257 + t239) * t208 + t307 * pkin(3) - t304) * t207) * t278 + (pkin(3) * t110 + t300 * t94 - t290) * t277 + t205 * t26 + pkin(1) * (-t204 * t94 + (t110 * t208 - t212 * t93) * t205) + (t212 * t110 + t208 * t93) * t293, (t12 + pkin(1) * (t117 * t208 - t212 * t74)) * t205 + (pkin(3) * t289 + t212 * (pkin(3) * t117 - t237) + pkin(7) * (t212 * t117 + t289) + t221 * (t207 * t109 + t211 * t218)) * t204, (t11 + pkin(1) * (t208 * t77 - t21 * t212)) * t205 + (t221 * t237 + t299 * (t208 * t21 + t212 * t77)) * t204, t205 * t49 + (t208 * t64 + t212 * (-t207 * t65 + t253)) * t204, t205 * t27 + (t208 * t40 + t212 * (t211 * t112 - t207 * t42)) * t204, t205 * t37 + (t208 * t57 + t212 * (-t207 * t59 + t211 * t72)) * t204, t205 * t48 + (t208 * t62 + t212 * (-t207 * t63 - t253)) * t204, t205 * t38 + (t208 * t58 + t212 * (-t207 * t60 - t211 * t68)) * t204, t205 * t61 + (t208 * t85 + t212 * (t211 * t124 - t207 * t86)) * t204, (t3 + pkin(1) * (t208 * t51 - t212 * t28)) * t205 + (t208 * (pkin(3) * t28 + t259) + t212 * (pkin(3) * t51 - t211 * t13 - t207 * t19) + pkin(7) * (t208 * t28 + t212 * t51) + t221 * (-t207 * t69 + t211 * t52)) * t204, (t4 + pkin(1) * (t208 * t54 - t212 * t35)) * t205 + (t208 * (pkin(3) * t35 + t260) + t212 * (pkin(3) * t54 - t211 * t14 - t207 * t20) + pkin(7) * (t208 * t35 + t212 * t54) + t221 * (-t207 * t73 + t211 * t55)) * t204, (t2 + pkin(1) * (t208 * t41 - t212 * t24)) * t205 + (t208 * (pkin(3) * t24 + t246) + t212 * (-t207 * t7 + (pkin(3) + t297) * t41) + pkin(7) * (t208 * t24 + t212 * t41) + t221 * (-t207 * t92 + t211 * t43)) * t204, (t1 + pkin(1) * (t208 * t9 - t212 * t5)) * t205 + ((t299 * t5 + t298) * t208 + (pkin(9) * t207 + t297 + t299) * t212 * t9 + t221 * (t211 * t10 + t207 * t33)) * t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t184, t175, -t307, t184, t145, t198, -t130, -t131, 0, 0, t198, t307, -t145, -t184, t175, t184, t116, t81, t80, t56, t104, t75, t97, t103, t98, t115, t23, t26, t12, t11, t49, t27, t37, t48, t38, t61, t3, t4, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t307, -t306, t152, t101, 0, 0, 0, 0, 0, 0, t87, t93, t74, t21, 0, 0, 0, 0, 0, 0, t28, t35, t24, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, t133, t109, -t134, t218, t162, -t46, -t47, 0, 0, t64, t40, t57, t62, t58, t85, t259, t260, t246, t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, t112, t72, -t113, -t68, t124, -t17, -t18, 0, 0;];
tauJ_reg = t6;
