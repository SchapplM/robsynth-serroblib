% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:21:05
% EndTime: 2019-03-09 05:21:11
% DurationCPUTime: 2.97s
% Computational Cost: add. (5940->373), mult. (12154->483), div. (0->0), fcn. (8728->14), ass. (0->220)
t199 = qJD(1) ^ 2;
t196 = cos(qJ(1));
t176 = g(2) * t196;
t192 = sin(qJ(1));
t177 = g(1) * t192;
t270 = t177 - t176;
t208 = -qJ(2) * t199 - t270;
t186 = qJ(3) + qJ(4);
t170 = sin(t186);
t171 = cos(t186);
t318 = g(3) * t170 - t171 * t270;
t169 = pkin(10) + t186;
t156 = cos(t169);
t178 = qJDD(3) + qJDD(4);
t190 = sin(qJ(4));
t194 = cos(qJ(4));
t195 = cos(qJ(3));
t266 = qJD(1) * t195;
t191 = sin(qJ(3));
t267 = qJD(1) * t191;
t116 = t190 * t267 - t194 * t266;
t197 = -pkin(1) - pkin(7);
t144 = t197 * qJD(1) + qJD(2);
t105 = -pkin(8) * t267 + t144 * t191;
t91 = t194 * t105;
t106 = -pkin(8) * t266 + t195 * t144;
t93 = qJD(3) * pkin(3) + t106;
t223 = -t190 * t93 - t91;
t142 = t197 * qJDD(1) + qJDD(2);
t129 = t195 * t142;
t265 = qJD(3) * t191;
t256 = t195 * qJDD(1);
t259 = qJD(1) * qJD(3);
t314 = t191 * t259 - t256;
t71 = qJDD(3) * pkin(3) + pkin(8) * t314 - t144 * t265 + t129;
t264 = qJD(3) * t195;
t249 = t195 * t259;
t257 = t191 * qJDD(1);
t315 = t249 + t257;
t77 = -pkin(8) * t315 + t142 * t191 + t144 * t264;
t205 = qJD(4) * t223 - t190 * t77 + t194 * t71;
t127 = t190 * t195 + t191 * t194;
t180 = qJD(3) + qJD(4);
t85 = t180 * t127;
t62 = -qJD(1) * t85 - t190 * t257 + t194 * t256;
t17 = pkin(4) * t178 - qJ(5) * t62 + qJD(5) * t116 + t205;
t187 = sin(pkin(10));
t188 = cos(pkin(10));
t117 = t127 * qJD(1);
t263 = qJD(4) * t190;
t304 = (qJD(4) * t93 + t77) * t194 - t105 * t263 + t190 * t71;
t251 = t191 * t263;
t215 = -qJD(1) * t251 - t190 * t314;
t234 = t180 * t195;
t63 = (qJD(1) * t234 + t257) * t194 + t215;
t19 = -qJ(5) * t63 - qJD(5) * t117 + t304;
t5 = t17 * t188 - t187 * t19;
t3 = -pkin(5) * t178 - t5;
t317 = t156 * t177 + t3;
t241 = t116 * t187 - t188 * t117;
t306 = qJD(6) - t241;
t189 = sin(qJ(6));
t193 = cos(qJ(6));
t221 = -t188 * t116 - t117 * t187;
t64 = -t193 * t180 + t189 * t221;
t316 = t306 * t64;
t155 = sin(t169);
t313 = g(3) * t155 + t156 * t176;
t239 = t193 * t306;
t38 = -t187 * t62 - t188 * t63;
t37 = qJDD(6) - t38;
t287 = t189 * t37;
t312 = -t239 * t306 - t287;
t260 = qJD(6) * t193;
t261 = qJD(6) * t189;
t39 = -t187 * t63 + t188 * t62;
t21 = t189 * t178 + t180 * t260 + t193 * t39 - t221 * t261;
t86 = -t190 * t265 + t194 * t234 - t251;
t51 = t187 * t86 + t188 * t85;
t66 = t180 * t189 + t193 * t221;
t128 = -t190 * t191 + t194 * t195;
t81 = t127 * t187 - t188 * t128;
t311 = -t21 * t81 - t51 * t66;
t220 = t188 * t127 + t128 * t187;
t225 = -t187 * t85 + t188 * t86;
t281 = qJ(5) * t117;
t59 = -t223 - t281;
t289 = t187 * t59;
t109 = t116 * qJ(5);
t90 = t190 * t105;
t243 = t194 * t93 - t90;
t58 = t109 + t243;
t54 = pkin(4) * t180 + t58;
t27 = t188 * t54 - t289;
t55 = t188 * t59;
t28 = t187 * t54 + t55;
t6 = t187 * t17 + t188 * t19;
t310 = t220 * t6 + t225 * t28 - t27 * t51 - t5 * t81 - t270;
t157 = pkin(4) * t187 + pkin(9);
t300 = pkin(4) * t116;
t42 = pkin(5) * t221 - pkin(9) * t241 - t300;
t309 = (qJD(6) * t157 + t42) * t306;
t162 = pkin(3) * t194 + pkin(4);
t278 = t188 * t190;
t108 = pkin(3) * t278 + t187 * t162;
t102 = pkin(9) + t108;
t166 = pkin(3) * t266;
t308 = (qJD(6) * t102 + t166 + t42) * t306;
t292 = pkin(8) - t197;
t133 = t292 * t191;
t134 = t292 * t195;
t272 = -t194 * t133 - t190 * t134;
t181 = qJDD(1) * qJ(2);
t229 = g(1) * t196 + g(2) * t192;
t182 = qJD(1) * qJD(2);
t253 = 0.2e1 * t182;
t305 = 0.2e1 * t181 + t253 - t229;
t123 = t292 * t265;
t124 = qJD(3) * t134;
t204 = -qJD(4) * t272 + t194 * t123 + t124 * t190;
t202 = qJ(5) * t85 - qJD(5) * t128 + t204;
t262 = qJD(4) * t194;
t210 = t190 * t123 - t194 * t124 + t133 * t263 - t134 * t262;
t36 = -qJ(5) * t86 - qJD(5) * t127 + t210;
t15 = t187 * t202 + t188 * t36;
t136 = pkin(3) * t267 + qJD(1) * qJ(2);
t87 = pkin(4) * t117 + qJD(5) + t136;
t40 = -pkin(5) * t241 - pkin(9) * t221 + t87;
t248 = pkin(9) * t178 + qJD(6) * t40 + t6;
t25 = -pkin(5) * t180 - t27;
t236 = t133 * t190 - t194 * t134;
t213 = -qJ(5) * t128 + t236;
t70 = -qJ(5) * t127 + t272;
t44 = t187 * t213 + t188 * t70;
t174 = t191 * pkin(3);
t159 = qJ(2) + t174;
t231 = pkin(4) * t127 + t159;
t45 = pkin(5) * t220 + pkin(9) * t81 + t231;
t302 = -t220 * t248 - (qJD(6) * t45 + t15) * t306 - t25 * t51 - t3 * t81 - t44 * t37;
t298 = t25 * t241;
t297 = t25 * t81;
t296 = t45 * t37;
t295 = t64 * t221;
t294 = t66 * t221;
t293 = t306 * t221;
t291 = t194 * t106 - t90;
t290 = pkin(3) * qJD(4);
t288 = t189 * t21;
t34 = t193 * t37;
t242 = -t106 * t190 - t91;
t214 = t242 + t281;
t60 = t109 + t291;
t286 = -t187 * t60 + t188 * t214 + (t187 * t194 + t278) * t290;
t279 = t187 * t190;
t285 = -t187 * t214 - t188 * t60 + (t188 * t194 - t279) * t290;
t284 = t128 * t178 - t85 * t180;
t283 = pkin(1) * qJDD(1);
t280 = t116 * t117;
t277 = t189 * t192;
t276 = t189 * t196;
t275 = t192 * t193;
t274 = t193 * t196;
t271 = t196 * pkin(1) + t192 * qJ(2);
t185 = t195 ^ 2;
t269 = t191 ^ 2 - t185;
t198 = qJD(3) ^ 2;
t268 = -t198 - t199;
t145 = pkin(3) * t264 + qJD(2);
t258 = qJDD(3) * t191;
t252 = t81 * t261;
t94 = pkin(3) * t315 + t181 + t182;
t207 = pkin(4) * t63 + qJDD(5) + t94;
t11 = -pkin(5) * t38 - pkin(9) * t39 + t207;
t26 = pkin(9) * t180 + t28;
t245 = qJD(6) * t26 - t11;
t240 = -t193 * t178 + t189 * t39;
t235 = qJD(6) * t220 + qJD(1);
t233 = qJDD(2) - t283;
t230 = pkin(4) * t86 + t145;
t227 = t221 * t28 + t241 * t27;
t226 = -t306 * t51 - t37 * t81;
t13 = t189 * t40 + t193 * t26;
t224 = t13 * t221 + t189 * t317 + t25 * t260;
t222 = -t127 * t178 - t180 * t86;
t12 = -t189 * t26 + t193 * t40;
t219 = -t12 * t221 + t193 * t313 + t25 * t261;
t218 = t34 + (t189 * t241 - t261) * t306;
t216 = g(3) * t156 - t248;
t107 = -pkin(3) * t279 + t162 * t188;
t212 = 0.2e1 * qJ(2) * t259 + qJDD(3) * t197;
t30 = t188 * t58 - t289;
t209 = -t157 * t37 + t30 * t306 - t298;
t206 = -t102 * t37 - t285 * t306 - t298;
t203 = -t197 * t198 + t305;
t201 = g(3) * t171 + t136 * t117 + t170 * t270 - t304;
t200 = t136 * t116 + t205 + t318;
t179 = -qJ(5) - pkin(8) - pkin(7);
t173 = t196 * qJ(2);
t168 = qJDD(3) * t195;
t158 = -pkin(4) * t188 - pkin(5);
t131 = pkin(4) * t170 + t174;
t101 = -pkin(5) - t107;
t100 = t155 * t274 - t277;
t99 = t155 * t276 + t275;
t98 = t155 * t275 + t276;
t97 = -t155 * t277 + t274;
t67 = t116 ^ 2 - t117 ^ 2;
t49 = -t116 * t180 + (-t180 * t266 - t257) * t194 - t215;
t48 = t117 * t180 + t62;
t43 = t187 * t70 - t188 * t213;
t29 = t187 * t58 + t55;
t22 = qJD(6) * t66 + t240;
t20 = pkin(5) * t225 + pkin(9) * t51 + t230;
t14 = t187 * t36 - t188 * t202;
t10 = t193 * t11;
t9 = t239 * t66 + t288;
t8 = -t294 - t312;
t7 = t218 + t295;
t1 = (t21 - t316) * t193 + (-t306 * t66 - t22) * t189;
t2 = [qJDD(1), t270, t229, qJDD(2) - t270 - 0.2e1 * t283, t305, -t233 * pkin(1) - g(1) * (-pkin(1) * t192 + t173) - g(2) * t271 + (t253 + t181) * qJ(2), qJDD(1) * t185 - 0.2e1 * t191 * t249, -0.2e1 * t191 * t256 + 0.2e1 * t259 * t269, -t191 * t198 + t168, -t195 * t198 - t258, 0, t191 * t203 + t195 * t212, -t191 * t212 + t195 * t203, t116 * t85 + t128 * t62, t116 * t86 + t117 * t85 - t127 * t62 - t128 * t63, t284, t222, 0, t145 * t117 + t94 * t127 + t136 * t86 + t159 * t63 - t170 * t229 + t178 * t236 + t180 * t204, -t145 * t116 + t94 * t128 - t136 * t85 + t159 * t62 - t171 * t229 - t178 * t272 - t180 * t210, t14 * t221 + t15 * t241 + t38 * t44 + t39 * t43 - t310, t6 * t44 + t28 * t15 - t5 * t43 - t27 * t14 + t207 * t231 + t87 * t230 - g(1) * (t131 * t196 + t173 + (-pkin(1) + t179) * t192) - g(2) * (t131 * t192 - t179 * t196 + t271) t193 * t311 + t66 * t252 -(-t189 * t66 - t193 * t64) * t51 - (-t288 - t193 * t22 + (t189 * t64 - t193 * t66) * qJD(6)) * t81, t193 * t226 + t21 * t220 + t225 * t66 + t252 * t306, t260 * t306 * t81 - t189 * t226 - t22 * t220 - t225 * t64, t220 * t37 + t225 * t306, -g(1) * t100 - g(2) * t98 + t10 * t220 + t12 * t225 + t14 * t64 + t43 * t22 + (t20 * t306 + t296 + (-t220 * t26 - t306 * t44 - t297) * qJD(6)) * t193 + t302 * t189, g(1) * t99 - g(2) * t97 - t13 * t225 + t14 * t66 + t43 * t21 + (-(-qJD(6) * t44 + t20) * t306 - t296 + t245 * t220 + qJD(6) * t297) * t189 + t302 * t193; 0, 0, 0, qJDD(1), -t199, t233 + t208, 0, 0, 0, 0, 0, t191 * t268 + t168, t195 * t268 - t258, 0, 0, 0, 0, 0, -qJD(1) * t117 + t284, qJD(1) * t116 + t222, t220 * t38 + t221 * t51 + t225 * t241 + t39 * t81, -qJD(1) * t87 + t310, 0, 0, 0, 0, 0, -t220 * t287 + t22 * t81 + t51 * t64 + (-t189 * t225 - t193 * t235) * t306, -t220 * t34 + (t189 * t235 - t193 * t225) * t306 - t311; 0, 0, 0, 0, 0, 0, t195 * t199 * t191, -t269 * t199, t256, -t257, qJDD(3), g(3) * t191 + t195 * t208 + t129, g(3) * t195 + (-t142 - t208) * t191, -t280, t67, t48, t49, t178, -t242 * t180 + (-t117 * t266 + t178 * t194 - t180 * t263) * pkin(3) + t200, t291 * t180 + (t116 * t266 - t178 * t190 - t180 * t262) * pkin(3) + t201, -t107 * t39 + t108 * t38 + t221 * t286 + t241 * t285 + t227, t6 * t108 + t5 * t107 - t87 * (t166 - t300) + g(3) * t131 + t285 * t28 - t286 * t27 - t270 * (pkin(3) * t195 + pkin(4) * t171) t9, t1, t8, t7, -t293, t101 * t22 + t286 * t64 + (-t317 - t308) * t193 + t206 * t189 + t219, t101 * t21 + t286 * t66 + t206 * t193 + (-t313 + t308) * t189 + t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t280, t67, t48, t49, t178, -t180 * t223 + t200, t180 * t243 + t201, -t29 * t221 - t30 * t241 + (t187 * t38 - t188 * t39) * pkin(4) + t227, t27 * t29 - t28 * t30 + (t116 * t87 + t187 * t6 + t188 * t5 + t318) * pkin(4), t9, t1, t8, t7, -t293, t158 * t22 - t29 * t64 + t209 * t189 + (-t317 - t309) * t193 + t219, t158 * t21 - t29 * t66 + t209 * t193 + (-t313 + t309) * t189 + t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t221 ^ 2 - t241 ^ 2, t221 * t27 - t241 * t28 + t207 - t229, 0, 0, 0, 0, 0, t218 - t295, -t294 + t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t64, -t64 ^ 2 + t66 ^ 2, t21 + t316, -t240 + (-qJD(6) + t306) * t66, t37, -g(1) * t97 - g(2) * t99 + t13 * t306 + t189 * t216 - t25 * t66 - t26 * t260 + t10, g(1) * t98 - g(2) * t100 + t12 * t306 + t189 * t245 + t193 * t216 + t25 * t64;];
tau_reg  = t2;
