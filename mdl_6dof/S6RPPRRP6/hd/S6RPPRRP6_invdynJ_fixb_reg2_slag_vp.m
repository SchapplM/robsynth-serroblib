% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRRP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRP6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:11:24
% EndTime: 2019-03-09 02:11:31
% DurationCPUTime: 3.56s
% Computational Cost: add. (3454->443), mult. (6035->505), div. (0->0), fcn. (3374->6), ass. (0->226)
t125 = sin(qJ(5));
t128 = cos(qJ(5));
t126 = sin(qJ(4));
t236 = qJD(1) * t126;
t86 = qJD(5) + t236;
t165 = t128 * t86;
t226 = t128 * qJD(4);
t129 = cos(qJ(4));
t235 = qJD(1) * t129;
t68 = t125 * t235 - t226;
t159 = t68 * t165;
t227 = qJD(5) * t129;
t206 = t125 * t227;
t219 = t129 * qJDD(1);
t27 = (t126 * t226 + t206) * qJD(1) - qJD(5) * t226 - t125 * qJDD(4) - t128 * t219;
t25 = t128 * t27;
t233 = qJD(4) * t125;
t70 = t128 * t235 + t233;
t268 = t70 * t86;
t225 = qJD(1) * qJD(4);
t205 = t126 * t225;
t28 = t70 * qJD(5) - t128 * qJDD(4) + (-t205 + t219) * t125;
t301 = t125 * (t28 + t268) + t159 + t25;
t294 = t28 - t268;
t127 = sin(qJ(1));
t130 = cos(qJ(1));
t289 = g(1) * t130 + g(2) * t127;
t298 = t129 * t289;
t124 = pkin(1) + qJ(3);
t288 = qJD(1) * t124;
t83 = -qJD(2) + t288;
t182 = -qJD(1) * t83 - t289;
t232 = qJD(4) * t126;
t261 = t128 * t70;
t262 = t128 * t68;
t264 = t128 * t28;
t265 = t125 * t27;
t297 = t129 * ((t125 * t68 - t261) * qJD(5) - t264 + t265) + (t125 * t70 + t262) * t232;
t120 = t126 ^ 2;
t121 = t129 ^ 2;
t238 = t120 + t121;
t99 = qJ(2) * qJD(1) + qJD(3);
t82 = -pkin(7) * qJD(1) + t99;
t296 = t238 * t82;
t230 = qJD(5) * t125;
t204 = t129 * t225;
t220 = t126 * qJDD(1);
t65 = qJDD(5) + t204 + t220;
t263 = t128 * t65;
t160 = -t86 * t230 + t263;
t166 = t125 * t86;
t213 = qJD(1) * t166;
t193 = t126 * t213 - t160;
t50 = t68 * t235;
t295 = t193 - t50;
t14 = t68 * t86 - t27;
t144 = -pkin(8) * qJD(5) * t86 - t298;
t117 = qJD(1) * qJD(2);
t118 = qJ(2) * qJDD(1);
t201 = qJDD(3) + t117 + t118;
t72 = -pkin(7) * qJDD(1) + t201;
t33 = -qJDD(4) * pkin(4) - t129 * t72 + t82 * t232;
t293 = -t33 + t144;
t292 = t129 * t28 - t68 * t232;
t110 = t129 * pkin(8);
t277 = pkin(4) * t126;
t290 = t110 - t277;
t56 = -qJD(4) * pkin(4) - t129 * t82;
t18 = pkin(5) * t68 - qJ(6) * t70 + t56;
t279 = pkin(8) * t65;
t287 = t18 * t86 - t279;
t228 = qJD(5) * t128;
t161 = t125 * t65 + t86 * t228;
t250 = t126 * t128;
t215 = t86 * t250;
t256 = t129 * t70;
t286 = qJD(1) * (t215 + t256) + t161;
t123 = -pkin(7) + qJ(2);
t285 = (qJD(2) + t83 + t288) * qJD(4) + qJDD(4) * t123;
t229 = qJD(5) * t126;
t207 = t123 * t229;
t231 = qJD(4) * t129;
t234 = qJD(2) * t126;
t275 = pkin(8) * t126;
t186 = pkin(4) * t129 + t275;
t66 = t186 * qJD(4) + qJD(3);
t195 = t124 + t277;
t77 = -t110 + t195;
t12 = -(qJD(5) * t77 + t123 * t231 + t234) * t125 - (-t66 + t207) * t128;
t283 = (t86 * t233 - t28) * t126 - t129 * (qJD(4) * t68 + t161);
t212 = t86 * t226;
t254 = qJD(4) * t70;
t282 = t126 * (-t160 + t254) - t129 * (-t27 + t212) + t213;
t281 = t70 ^ 2;
t108 = 0.2e1 * t117;
t280 = pkin(5) * t65;
t278 = pkin(8) * t70;
t276 = pkin(7) * t130;
t114 = g(1) * t127;
t113 = g(2) * t130;
t274 = g(3) * t126;
t273 = g(3) * t129;
t49 = -qJD(2) + (-t290 + t124) * qJD(1);
t75 = t126 * t82;
t55 = qJD(4) * pkin(8) + t75;
t21 = t125 * t49 + t128 * t55;
t17 = qJ(6) * t86 + t21;
t272 = t17 * t86;
t271 = t21 * t86;
t269 = t70 * t68;
t177 = pkin(5) * t125 - qJ(6) * t128;
t267 = -qJD(6) * t125 + t177 * t86 - t75;
t245 = t128 * t129;
t74 = t186 * qJD(1);
t36 = t125 * t74 + t82 * t245;
t45 = t123 * t250 + t125 * t77;
t266 = qJ(6) * t65;
t260 = t128 * t74;
t259 = t128 * t77;
t257 = t129 * t56;
t132 = qJD(1) ^ 2;
t253 = t120 * t132;
t252 = t125 * t126;
t251 = t125 * t129;
t249 = t126 * t130;
t248 = t127 * t125;
t247 = t127 * t128;
t246 = t127 * t129;
t244 = t129 * t130;
t243 = t130 * t128;
t20 = -t125 * t55 + t128 * t49;
t242 = qJD(6) - t20;
t241 = t130 * pkin(1) + t127 * qJ(2);
t240 = t114 - t113;
t131 = qJD(4) ^ 2;
t237 = -t131 - t132;
t224 = qJD(3) * qJD(1);
t223 = qJDD(1) * t124;
t221 = qJDD(4) * t126;
t122 = qJDD(1) * pkin(1);
t218 = t122 - qJDD(2);
t214 = t68 ^ 2 - t281;
t211 = t129 * t132 * t126;
t210 = t86 * t235;
t209 = t130 * qJ(3) + t241;
t208 = t129 * t226;
t203 = qJDD(2) - t240;
t116 = qJDD(1) * qJ(3);
t202 = -t116 - t218;
t200 = t238 * t72;
t26 = t66 * qJD(1) - qJDD(1) * t290 - t202;
t34 = qJDD(4) * pkin(8) + t126 * t72 + t82 * t231;
t5 = -t125 * t34 + t128 * t26 - t55 * t228 - t49 * t230;
t199 = t123 * t208 + t125 * t66 + t128 * t234 + t77 * t228;
t198 = qJD(5) * t68 - t27;
t197 = qJD(1) + t229;
t196 = 0.2e1 * t118 + t108 - t289;
t107 = t130 * qJ(2);
t194 = pkin(8) * t246 + t107 - t276;
t192 = t126 * t204;
t191 = -t122 + t203;
t57 = t126 * t248 - t243;
t59 = t125 * t249 + t247;
t190 = -g(1) * t57 + g(2) * t59;
t58 = t125 * t130 + t126 * t247;
t60 = t126 * t243 - t248;
t189 = g(1) * t58 - g(2) * t60;
t181 = t198 * pkin(8);
t179 = -t124 * t127 + t107;
t178 = pkin(5) * t128 + qJ(6) * t125;
t16 = -pkin(5) * t86 + t242;
t176 = t125 * t17 - t128 * t16;
t175 = t125 * t16 + t128 * t17;
t174 = t125 * t21 + t128 * t20;
t173 = t125 * t20 - t128 * t21;
t73 = -t202 + t224;
t169 = t83 * qJD(3) + t73 * t124;
t167 = -t116 + t191;
t78 = -pkin(4) - t178;
t162 = -t123 + t177;
t4 = t125 * t26 + t128 * t34 + t49 * t228 - t55 * t230;
t157 = -g(1) * (pkin(4) * t244 + pkin(8) * t249) - g(2) * (pkin(4) * t246 + t127 * t275);
t156 = pkin(4) * t249 - pkin(8) * t244 + t209;
t154 = t86 * t56 - t279;
t153 = -t72 - t182;
t152 = t294 * t125;
t149 = g(1) * t59 + g(2) * t57 + g(3) * t251 + t5;
t147 = t68 * t166 - t264;
t146 = (g(2) * pkin(7) + g(1) * t195) * t127;
t145 = -pkin(8) * t264 - t126 * t289 - t273;
t3 = pkin(5) * t28 + qJ(6) * t27 - qJD(6) * t70 + t33;
t143 = -t3 + t144;
t142 = -t123 * t131 + t114 + t223 + t224 + t73;
t1 = qJD(6) * t86 + t266 + t4;
t2 = qJDD(6) - t5 - t280;
t141 = -t176 * qJD(5) + t1 * t128 + t125 * t2;
t140 = -t175 * qJD(5) - t1 * t125 + t128 * t2;
t139 = -t174 * qJD(5) - t5 * t125 + t4 * t128;
t138 = t173 * qJD(5) - t125 * t4 - t128 * t5;
t137 = t18 * t70 + qJDD(6) - t149;
t136 = t292 * t125 + t227 * t262;
t135 = -g(1) * t60 - g(2) * t58 - g(3) * t245 + t4;
t134 = -t65 * t252 + (-t125 * t231 - t197 * t128) * t86 - t292;
t133 = -t28 * t250 - t68 * t208 + t197 * t261 + (qJD(1) * t68 + t198 * t126 + t70 * t231) * t125;
t104 = t121 * t132;
t102 = qJDD(4) * t129;
t95 = g(2) * t244;
t93 = g(3) * t250;
t46 = t162 * t129;
t44 = -t123 * t252 + t259;
t41 = -t259 + (t123 * t125 - pkin(5)) * t126;
t40 = qJ(6) * t126 + t45;
t39 = t126 * t65 + t86 * t231;
t37 = pkin(5) * t70 + qJ(6) * t68;
t35 = -t82 * t251 + t260;
t32 = -t260 + (-pkin(5) * qJD(1) + t125 * t82) * t129;
t29 = qJ(6) * t235 + t36;
t15 = -t162 * t232 + (t178 * qJD(5) - qJD(6) * t128 - qJD(2)) * t129;
t13 = (t215 - t256) * qJD(1) + t161;
t11 = -t125 * t207 + t199;
t10 = -pkin(5) * t231 - t12;
t9 = qJ(6) * t231 + (-t123 * t230 + qJD(6)) * t126 + t199;
t8 = t70 * t165 - t265;
t7 = -t70 * t206 + (-t129 * t27 - t70 * t232) * t128;
t6 = (-t27 - t212) * t126 + (t160 + t254) * t129;
t19 = [0, 0, 0, 0, 0, qJDD(1), t240, t289, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -0.2e1 * t122 + t203, t196, t218 * pkin(1) - g(1) * (-pkin(1) * t127 + t107) - g(2) * t241 + (t108 + t118) * qJ(2), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) + t196, -t167 + t223 + 0.2e1 * t224, -g(1) * t179 - g(2) * t209 + t201 * qJ(2) + t99 * qJD(2) + t169, qJDD(1) * t121 - 0.2e1 * t192, -0.2e1 * t126 * t219 + 0.2e1 * (t120 - t121) * t225, -t126 * t131 + t102, qJDD(1) * t120 + 0.2e1 * t192, -t129 * t131 - t221, 0, t285 * t129 + (t142 - t113) * t126, -t126 * t285 + t142 * t129 - t95, t289 + t238 * (-qJDD(1) * t123 - t117 - t72) -g(1) * (t179 - t276) - g(2) * (-pkin(7) * t127 + t209) + t123 * t200 + qJD(2) * t296 + t169, t7, t297, t6, t136, t283, t39, t12 * t86 + t44 * t65 + (t5 + (t123 * t68 - t125 * t56) * qJD(4)) * t126 + (-qJD(2) * t68 + qJD(4) * t20 - t123 * t28 + t125 * t33 + t228 * t56) * t129 + t189, -t11 * t86 - t45 * t65 + (-t4 + (t123 * t70 - t128 * t56) * qJD(4)) * t126 + (-qJD(2) * t70 - qJD(4) * t21 + t123 * t27 + t128 * t33 - t230 * t56) * t129 + t190, -t11 * t68 - t12 * t70 + t27 * t44 - t28 * t45 + t95 + t174 * t232 + (t138 - t114) * t129, t4 * t45 + t21 * t11 + t5 * t44 + t20 * t12 - qJD(2) * t257 - g(1) * t194 - g(2) * t156 + t146 + (-t129 * t33 + t232 * t56) * t123, t7, t6, -t297, t39, -t283, t136, -t10 * t86 + t15 * t68 + t28 * t46 - t41 * t65 + (-t18 * t233 - t2) * t126 + (-qJD(4) * t16 + t125 * t3 + t18 * t228) * t129 + t189, t10 * t70 - t27 * t41 - t28 * t40 - t68 * t9 + t95 + t176 * t232 + (t140 - t114) * t129, -t15 * t70 + t27 * t46 + t40 * t65 + t86 * t9 + (t18 * t226 + t1) * t126 + (qJD(4) * t17 - t128 * t3 + t18 * t230) * t129 - t190, t1 * t40 + t17 * t9 + t3 * t46 + t18 * t15 + t2 * t41 + t16 * t10 - g(1) * (-pkin(5) * t58 - qJ(6) * t57 + t194) - g(2) * (pkin(5) * t60 + qJ(6) * t59 + t156) + t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t132, -qJ(2) * t132 + t191, 0, 0, 0, 0, 0, 0, 0, -t132, -qJDD(1) (-qJD(3) - t99) * qJD(1) + t167, 0, 0, 0, 0, 0, 0, -0.2e1 * t204 - t220, 0.2e1 * t205 - t219, t104 + t253 (-qJD(3) - t296) * qJD(1) + t167, 0, 0, 0, 0, 0, 0, t193 + t50, t286, t128 * t14 + t152 (t126 * t173 + t257) * qJD(1) + t138 - t240, 0, 0, 0, 0, 0, 0, t166 * t86 - t263 + t50, -t25 + t159 + t152, -t286 (-t126 * t175 + t129 * t18) * qJD(1) + t140 - t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t132, t182 + t201, 0, 0, 0, 0, 0, 0, t237 * t126 + t102, t237 * t129 - t221, -t238 * qJDD(1), t200 + t182, 0, 0, 0, 0, 0, 0, t134, t282, t133, -t174 * qJD(1) + (-qJD(4) * t173 - t33) * t129 + (qJD(4) * t56 + t139) * t126 - t289, 0, 0, 0, 0, 0, 0, t134, t133, -t282, -t176 * qJD(1) + (qJD(4) * t175 - t3) * t129 + (qJD(4) * t18 + t141) * t126 - t289; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, t104 - t253, t219, -t211, -t220, qJDD(4), -t153 * t129 + t274, t153 * t126 + t273, 0, 0, t8, -t301, t13, t147, -t295, -t210, -pkin(4) * t28 + t154 * t125 + t128 * t293 - t20 * t235 - t35 * t86 - t68 * t75 + t93, t21 * t235 - t70 * t75 + pkin(4) * t27 + t36 * t86 + t154 * t128 + (-t274 - t293) * t125, t35 * t70 + t36 * t68 + (-t20 * t236 + t4 + (-t20 + t278) * qJD(5)) * t128 + (t181 - t5 - t271) * t125 + t145, -t33 * pkin(4) + t139 * pkin(8) - g(3) * t290 - t20 * t35 - t21 * t36 - t56 * t75 + t157, t8, t13, t301, -t210, t295, t147, t125 * t287 + t143 * t128 + t16 * t235 + t267 * t68 + t28 * t78 + t32 * t86 + t93, t29 * t68 - t32 * t70 + (t16 * t236 + t1 + (t16 + t278) * qJD(5)) * t128 + (t181 + t2 - t272) * t125 + t145, -t17 * t235 + t27 * t78 - t29 * t86 - t267 * t70 - t287 * t128 + (t143 + t274) * t125, t141 * pkin(8) - g(3) * t110 - t16 * t32 - t17 * t29 + t267 * t18 + t157 + (t3 - t274) * t78 - t178 * t298; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t269, -t214, t14, -t269, -t294, t65, -t56 * t70 + t149 + t271, t20 * t86 + t56 * t68 - t135, 0, 0, t269, t14, t214, t65, t294, -t269, -t37 * t68 - t137 + t271 + 0.2e1 * t280, pkin(5) * t27 - qJ(6) * t28 + (t17 - t21) * t70 + (t16 - t242) * t68, 0.2e1 * t266 - t18 * t68 + t37 * t70 + (0.2e1 * qJD(6) - t20) * t86 + t135, t1 * qJ(6) - t2 * pkin(5) - t18 * t37 - t16 * t21 - g(1) * (-pkin(5) * t59 + qJ(6) * t60) - g(2) * (-pkin(5) * t57 + qJ(6) * t58) + t242 * t17 + t177 * t273; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65 + t269, t14, -t86 ^ 2 - t281, t137 - t272 - t280;];
tau_reg  = t19;
