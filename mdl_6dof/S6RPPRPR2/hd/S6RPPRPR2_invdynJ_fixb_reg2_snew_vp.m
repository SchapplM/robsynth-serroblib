% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRPR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:02:50
% EndTime: 2019-05-05 14:03:03
% DurationCPUTime: 5.13s
% Computational Cost: add. (10296->367), mult. (23687->500), div. (0->0), fcn. (16470->10), ass. (0->202)
t209 = qJD(4) ^ 2;
t200 = cos(pkin(10));
t204 = sin(qJ(4));
t198 = sin(pkin(10));
t207 = cos(qJ(4));
t246 = t198 * t207;
t223 = t200 * t204 + t246;
t168 = t223 * qJD(1);
t268 = t168 ^ 2;
t148 = t268 + t209;
t247 = t198 * t204;
t166 = (-t200 * t207 + t247) * qJD(1);
t249 = t168 * t166;
t286 = qJDD(4) + t249;
t297 = t286 * t204;
t91 = t148 * t207 + t297;
t309 = pkin(7) * t91;
t296 = t286 * t207;
t93 = -t148 * t204 + t296;
t308 = pkin(7) * t93;
t307 = t198 * t93 + t200 * t91;
t60 = t198 * t91 - t200 * t93;
t149 = t268 - t209;
t287 = qJDD(4) - t249;
t299 = t207 * t287;
t300 = t204 * t287;
t306 = t198 * (t149 * t204 + t299) - t200 * (t149 * t207 - t300);
t269 = t166 ^ 2;
t144 = t269 - t209;
t303 = t198 * (-t144 * t207 + t297) - t200 * (t144 * t204 + t296);
t155 = t166 * qJD(4);
t165 = t223 * qJDD(1);
t133 = t165 - t155;
t277 = t155 - t133;
t298 = t277 * qJ(5);
t123 = -t209 - t269;
t75 = t123 * t204 + t299;
t295 = pkin(7) * t75;
t78 = -t123 * t207 + t300;
t294 = pkin(7) * t78;
t191 = t198 ^ 2;
t192 = t200 ^ 2;
t210 = qJD(1) ^ 2;
t177 = (t191 + t192) * t210;
t205 = sin(qJ(1));
t208 = cos(qJ(1));
t228 = g(1) * t208 + g(2) * t205;
t175 = -pkin(1) * t210 - t228;
t199 = sin(pkin(9));
t201 = cos(pkin(9));
t227 = g(1) * t205 - t208 * g(2);
t222 = qJDD(1) * pkin(1) + t227;
t245 = t201 * t175 + t199 * t222;
t278 = -pkin(2) * t210 + qJDD(1) * qJ(3) + 0.2e1 * qJD(1) * qJD(3) + t245;
t285 = t198 * t78 - t200 * t75;
t284 = t198 * t75 + t200 * t78;
t271 = -t269 - t268;
t283 = pkin(2) * t271;
t282 = pkin(3) * t271;
t281 = t201 * t271;
t203 = sin(qJ(6));
t206 = cos(qJ(6));
t140 = qJD(4) * t203 - t206 * t166;
t142 = qJD(4) * t206 + t166 * t203;
t106 = t142 * t140;
t120 = qJDD(6) + t133;
t274 = -t106 + t120;
t280 = t203 * t274;
t279 = t206 * t274;
t231 = -t199 * t175 + t201 * t222;
t111 = -qJDD(1) * pkin(2) - t210 * qJ(3) + qJDD(3) - t231;
t234 = pkin(1) * t199 + qJ(3);
t235 = pkin(1) * t201 + pkin(2);
t276 = -t235 * qJDD(1) + t234 * t177 + t111;
t159 = qJD(6) + t168;
t109 = t159 * t140;
t240 = t200 * qJDD(1);
t241 = t198 * qJDD(1);
t225 = t204 * t241 - t207 * t240;
t244 = t168 * qJD(4);
t131 = t225 + t244;
t87 = -t140 * qJD(6) + t206 * qJDD(4) + t203 * t131;
t275 = -t109 + t87;
t101 = t155 + t133;
t272 = t268 - t269;
t138 = t140 ^ 2;
t139 = t142 ^ 2;
t157 = t159 ^ 2;
t267 = 2 * qJD(5);
t266 = -pkin(4) - pkin(8);
t265 = pkin(4) * t207;
t195 = -g(3) + qJDD(2);
t183 = t200 * t195;
t85 = t183 + (pkin(3) * t200 * t210 - pkin(7) * qJDD(1) - t278) * t198;
t248 = t192 * t210;
t98 = t198 * t195 + t278 * t200;
t89 = -pkin(3) * t248 + pkin(7) * t240 + t98;
t55 = t204 * t89 - t207 * t85;
t56 = t204 * t85 + t207 * t89;
t27 = t204 * t56 - t207 * t55;
t264 = t198 * t27;
t143 = pkin(5) * t168 - qJD(4) * pkin(8);
t119 = pkin(4) * t166 - qJ(5) * t168;
t220 = -t209 * pkin(4) - t119 * t166 + t56;
t242 = qJDD(4) * qJ(5);
t30 = t242 - t131 * pkin(5) - t269 * pkin(8) + (t267 + t143) * qJD(4) + t220;
t263 = t203 * t30;
t74 = t106 + t120;
t262 = t203 * t74;
t96 = -pkin(3) * t240 + t111 + (-t191 * t210 - t248) * pkin(7);
t261 = t204 * t96;
t260 = t206 * t30;
t259 = t206 * t74;
t258 = t207 * t96;
t251 = t159 * t203;
t250 = t159 * t206;
t239 = -t139 - t157;
t237 = t204 * t106;
t236 = t207 * t106;
t233 = qJ(5) * t204 + pkin(3);
t97 = t278 * t198 - t183;
t61 = t198 * t97 + t200 * t98;
t44 = -qJDD(4) * pkin(4) - t209 * qJ(5) + t119 * t168 + qJDD(5) + t55;
t29 = pkin(5) * t101 - pkin(8) * t287 + t44;
t218 = t131 * pkin(4) + t298 + t96;
t232 = pkin(4) * qJD(4) - (2 * qJD(5));
t33 = (-t143 + t232) * t168 - pkin(5) * t269 + pkin(8) * t131 + t218;
t17 = t203 * t33 - t206 * t29;
t28 = t204 * t55 + t207 * t56;
t230 = qJDD(4) * t203 - t206 * t131;
t18 = t203 * t29 + t206 * t33;
t7 = -t17 * t206 + t18 * t203;
t8 = t203 * t17 + t18 * t206;
t221 = (-qJD(6) + t159) * t142 - t230;
t219 = qJD(4) * t267 + t220;
t217 = t198 * (t207 * t133 - t204 * t244) + t200 * (t204 * t133 + t207 * t244);
t43 = t219 + t242;
t215 = t198 * (t131 * t204 + t207 * t155) + t200 * (-t207 * t131 + t204 * t155);
t214 = -t168 * t267 + t218;
t213 = (t198 * (-t166 * t207 + t168 * t204) + t200 * (-t166 * t204 - t168 * t207)) * qJD(4);
t188 = t192 * qJDD(1);
t187 = t191 * qJDD(1);
t176 = t188 + t187;
t132 = t165 - 0.2e1 * t155;
t130 = t225 + 0.2e1 * t244;
t108 = -t139 + t157;
t107 = t138 - t157;
t103 = t139 - t138;
t100 = t131 + t244;
t99 = t131 - t244;
t88 = -t157 - t138;
t86 = -qJD(6) * t142 - t230;
t83 = -t138 - t139;
t82 = t101 * t204 - t207 * t225;
t81 = t165 * t204 - t207 * t99;
t80 = -t101 * t207 - t204 * t225;
t79 = -t165 * t207 - t204 * t99;
t71 = (t140 * t203 + t142 * t206) * t159;
t69 = t109 + t87;
t65 = (qJD(6) + t159) * t142 + t230;
t63 = -t142 * t250 - t203 * t87;
t62 = -t140 * t251 - t206 * t86;
t58 = -t107 * t203 - t259;
t57 = -t108 * t206 - t280;
t53 = -t203 * t239 - t259;
t52 = t206 * t239 - t262;
t51 = -t198 * t80 + t200 * t82;
t50 = -t198 * t79 + t200 * t81;
t49 = t206 * t88 - t280;
t48 = t203 * t88 + t279;
t45 = t232 * t168 + t218;
t42 = t203 * t69 + t206 * t221;
t41 = t203 * t221 - t206 * t69;
t40 = t203 * t65 - t206 * t275;
t39 = (t100 + t244) * pkin(4) + t214;
t38 = -pkin(4) * t244 - t214 - t298;
t37 = -qJ(5) * t271 + t44;
t36 = -pkin(4) * t271 + t43;
t35 = t204 * t52 + t207 * t275;
t34 = t204 * t275 - t207 * t52;
t32 = t204 * t48 + t207 * t65;
t31 = t204 * t65 - t207 * t48;
t25 = t204 * t41 + t207 * t83;
t24 = t204 * t83 - t207 * t41;
t23 = t204 * t44 + t207 * t43;
t22 = t204 * t43 - t207 * t44;
t21 = pkin(5) * t41 - qJ(5) * t42;
t20 = -t198 * t34 + t200 * t35;
t19 = -t198 * t31 + t200 * t32;
t15 = t200 * t28 - t264;
t14 = -t198 * t24 + t200 * t25;
t13 = pkin(5) * t275 + t266 * t53 - t263;
t12 = pkin(5) * t65 + t266 * t49 + t260;
t10 = pkin(5) * t52 - qJ(5) * t53 - t18;
t9 = pkin(5) * t48 - qJ(5) * t49 - t17;
t6 = t204 * t7 + t207 * t30;
t5 = t204 * t30 - t207 * t7;
t4 = pkin(5) * t83 + t266 * t42 - t8;
t3 = pkin(5) * t7 - qJ(5) * t8;
t2 = pkin(5) * t30 + t266 * t8;
t1 = -t198 * t5 + t200 * t6;
t11 = [0, 0, 0, 0, 0, qJDD(1), t227, t228, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t201 - t199 * t210) + t231, pkin(1) * (-qJDD(1) * t199 - t201 * t210) - t245, 0, pkin(1) * (t199 * t245 + t201 * t231), t187, 0.2e1 * t198 * t240, 0, t188, 0, 0, -t276 * t200, t276 * t198, pkin(2) * t177 + qJ(3) * t176 + pkin(1) * (t176 * t199 + t177 * t201) + t61, -pkin(2) * t111 + qJ(3) * t61 + pkin(1) * (-t111 * t201 + t199 * t61), t217, t198 * (-t130 * t207 - t132 * t204) + t200 * (-t130 * t204 + t132 * t207), t306, t215, -t303, t213, t198 * (t261 - t295) + t200 * (-pkin(3) * t130 - t258 - t294) - pkin(2) * t130 - qJ(3) * t284 + pkin(1) * (-t130 * t201 - t199 * t284), t198 * (t258 + t309) + t200 * (-pkin(3) * t132 + t261 - t308) - pkin(2) * t132 + qJ(3) * t60 + pkin(1) * (-t132 * t201 + t199 * t60), t198 * (-pkin(7) * t79 - t27) + t200 * (pkin(7) * t81 + t28 - t282) - t283 + qJ(3) * t50 + pkin(1) * (t199 * t50 - t281), -pkin(7) * t264 + t200 * (-pkin(3) * t96 + pkin(7) * t28) - pkin(2) * t96 + qJ(3) * t15 + pkin(1) * (t15 * t199 - t201 * t96), t213, -t306, t303, t217, t198 * (-t100 * t207 + t204 * t277) + t200 * (-t100 * t204 - t207 * t277), t215, t198 * (-pkin(7) * t80 - t204 * t36 + t207 * t37) + t200 * (pkin(7) * t82 + t204 * t37 + t207 * t36 - t282) - t283 + qJ(3) * t51 + pkin(1) * (t199 * t51 - t281), t198 * (-t204 * t39 + t295) + t200 * (t207 * t39 + t294) + t234 * t284 + (qJ(5) * t246 + t200 * t233 + t235) * t100, t198 * (t207 * t38 - t309) + t200 * (t204 * t38 + t308) - t234 * t60 - (-pkin(4) * t247 + t200 * (pkin(3) + t265) + t235) * t277, (t198 * (pkin(4) * t204 - qJ(5) * t207) + t200 * (-t233 - t265) - t235) * t45 + (t234 + pkin(7)) * (-t198 * t22 + t200 * t23), t198 * (-t204 * t63 + t236) + t200 * (t207 * t63 + t237), t198 * (t103 * t207 - t204 * t40) + t200 * (t103 * t204 + t207 * t40), t198 * (-t204 * t57 + t207 * t69) + t200 * (t204 * t69 + t207 * t57), t198 * (-t204 * t62 - t236) + t200 * (t207 * t62 - t237), t198 * (-t204 * t58 + t207 * t221) + t200 * (t204 * t221 + t207 * t58), t198 * (t120 * t207 - t204 * t71) + t200 * (t120 * t204 + t207 * t71), t198 * (-pkin(7) * t31 - t12 * t204 + t207 * t9) + t200 * (-pkin(3) * t49 + pkin(7) * t32 + t12 * t207 + t204 * t9) - pkin(2) * t49 + qJ(3) * t19 + pkin(1) * (t19 * t199 - t201 * t49), t198 * (-pkin(7) * t34 + t10 * t207 - t13 * t204) + t200 * (-pkin(3) * t53 + pkin(7) * t35 + t10 * t204 + t13 * t207) - pkin(2) * t53 + qJ(3) * t20 + pkin(1) * (t199 * t20 - t201 * t53), t198 * (-pkin(7) * t24 - t204 * t4 + t207 * t21) + t200 * (-pkin(3) * t42 + pkin(7) * t25 + t204 * t21 + t207 * t4) - pkin(2) * t42 + qJ(3) * t14 + pkin(1) * (t14 * t199 - t201 * t42), t198 * (-pkin(7) * t5 - t2 * t204 + t207 * t3) + t200 * (-pkin(3) * t8 + pkin(7) * t6 + t2 * t207 + t204 * t3) - pkin(2) * t8 + qJ(3) * t1 + pkin(1) * (t1 * t199 - t201 * t8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198 * t98 - t200 * t97, 0, 0, 0, 0, 0, 0, -t285, -t307, t198 * t81 + t200 * t79, t198 * t28 + t200 * t27, 0, 0, 0, 0, 0, 0, t198 * t82 + t200 * t80, t285, t307, t198 * t23 + t200 * t22, 0, 0, 0, 0, 0, 0, t198 * t32 + t200 * t31, t198 * t35 + t200 * t34, t198 * t25 + t200 * t24, t198 * t6 + t200 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t240, t241, -t177, t111, 0, 0, 0, 0, 0, 0, t130, t132, t271, t96, 0, 0, 0, 0, 0, 0, t271, -t100, t277, t45, 0, 0, 0, 0, 0, 0, t49, t53, t42, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t249, t272, t165, -t249, -t225, qJDD(4), -t55, -t56, 0, 0, qJDD(4), -t101, t99, t249, t272, -t249, -pkin(4) * t101 - qJ(5) * t225, -pkin(4) * t287 - qJ(5) * t123 + t44, pkin(4) * t148 + (qJDD(4) + t286) * qJ(5) + t219, -pkin(4) * t44 + qJ(5) * t43, -t142 * t251 + t206 * t87, -t203 * t275 - t206 * t65, -t108 * t203 + t279, t140 * t250 - t203 * t86, t107 * t206 - t262, (-t140 * t206 + t142 * t203) * t159, qJ(5) * t65 + t266 * t48 + t263, qJ(5) * t275 + t266 * t52 + t260, qJ(5) * t83 + t266 * t41 - t7, qJ(5) * t30 + t266 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t287, -t148, t44, 0, 0, 0, 0, 0, 0, t48, t52, t41, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t103, t69, -t106, t221, t120, -t17, -t18, 0, 0;];
tauJ_reg  = t11;
