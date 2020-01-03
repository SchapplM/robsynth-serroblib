% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:26:42
% EndTime: 2019-12-31 20:26:53
% DurationCPUTime: 4.03s
% Computational Cost: add. (5194->438), mult. (15047->634), div. (0->0), fcn. (12484->12), ass. (0->239)
t175 = sin(pkin(5));
t184 = cos(qJ(2));
t170 = pkin(2) * t184 + pkin(1);
t220 = t170 * qJDD(1);
t180 = sin(qJ(2));
t259 = qJD(1) * qJD(2);
t243 = t180 * t259;
t226 = t175 * t243;
t255 = pkin(2) * t226 + qJDD(3);
t318 = t175 * t220 - t255;
t174 = sin(pkin(10));
t267 = qJD(1) * t180;
t249 = t175 * t267;
t176 = cos(pkin(10));
t272 = t184 * t176;
t251 = t175 * t272;
t128 = qJD(1) * t251 - t174 * t249;
t124 = qJD(4) - t128;
t211 = t174 * t184 + t176 * t180;
t198 = qJD(1) * t211;
t131 = t175 * t198;
t183 = cos(qJ(4));
t177 = cos(pkin(5));
t268 = qJD(1) * t177;
t231 = qJD(2) + t268;
t152 = t183 * t231;
t179 = sin(qJ(4));
t103 = t131 * t179 - t152;
t102 = qJD(5) + t103;
t181 = sin(qJ(1));
t273 = t181 * t184;
t185 = cos(qJ(1));
t274 = t180 * t185;
t140 = -t177 * t273 - t274;
t271 = t184 * t185;
t275 = t180 * t181;
t206 = t177 * t271 - t275;
t317 = -g(1) * t140 - g(2) * t206;
t178 = sin(qJ(5));
t182 = cos(qJ(5));
t135 = t211 * t177;
t143 = t174 * t180 - t272;
t213 = t135 * t185 - t143 * t181;
t278 = t175 * t185;
t82 = -t179 * t278 + t183 * t213;
t199 = t143 * t177;
t95 = -t181 * t211 - t185 * t199;
t316 = t178 * t82 + t182 * t95;
t315 = -t178 * t95 + t182 * t82;
t190 = -t183 * t131 - t179 * t231;
t55 = -t182 * t124 - t178 * t190;
t314 = t124 * t55;
t261 = qJD(5) * t178;
t262 = qJD(4) * t183;
t285 = t128 * t183;
t80 = t131 * t178 + t182 * t285;
t191 = -t179 * t261 + t182 * t262 - t80;
t312 = t191 * t102;
t258 = qJDD(1) * t177;
t162 = qJDD(2) + t258;
t263 = qJD(4) * t179;
t244 = t184 * t259;
t93 = -t174 * t226 + (qJDD(1) * t211 + t176 * t244) * t175;
t37 = qJD(4) * t152 - t131 * t263 + t179 * t162 + t183 * t93;
t57 = t124 * t178 - t182 * t190;
t256 = t184 * qJDD(1);
t240 = t175 * t256;
t148 = t176 * t240;
t197 = t211 * qJD(2);
t257 = qJDD(1) * t180;
t242 = t174 * t257;
t91 = qJDD(4) - t148 + (qJD(1) * t197 + t242) * t175;
t12 = qJD(5) * t57 + t178 * t37 - t182 * t91;
t134 = t211 * t175;
t113 = t134 * t179 - t177 * t183;
t209 = t179 * t213 + t183 * t278;
t212 = -t135 * t181 - t143 * t185;
t280 = t175 * t181;
t84 = -t179 * t212 + t183 * t280;
t201 = g(1) * t84 - g(2) * t209 - g(3) * t113;
t307 = pkin(1) * t184;
t165 = t177 * t307;
t160 = qJD(1) * t165;
t304 = pkin(7) + qJ(3);
t245 = t304 * t180;
t227 = t175 * t245;
t196 = t177 * pkin(2) - t227;
t106 = qJD(2) * pkin(2) + qJD(1) * t196 + t160;
t276 = t177 * t180;
t164 = pkin(1) * t276;
t279 = t175 * t184;
t118 = (t304 * t279 + t164) * qJD(1);
t277 = t176 * t118;
t51 = t174 * t106 + t277;
t44 = pkin(8) * t231 + t51;
t225 = t170 * t175;
t137 = -qJD(1) * t225 + qJD(3);
t61 = -t128 * pkin(3) - t131 * pkin(8) + t137;
t23 = t179 * t61 + t183 * t44;
t253 = pkin(1) * t256;
t159 = t177 * t253;
t254 = pkin(1) * qJD(2) * t177;
t228 = qJD(1) * t254;
t235 = qJD(2) * t304;
t265 = qJD(3) * t180;
t48 = -t180 * t228 + pkin(2) * t162 + t159 + (-qJDD(1) * t245 + (-t184 * t235 - t265) * qJD(1)) * t175;
t193 = qJD(3) * t184 - t180 * t235;
t250 = pkin(7) * t240 + qJDD(1) * t164 + t184 * t228;
t58 = (qJ(3) * t256 + qJD(1) * t193) * t175 + t250;
t21 = t174 * t48 + t176 * t58;
t18 = pkin(8) * t162 + t21;
t92 = t148 + (-qJD(2) * t198 - t242) * t175;
t33 = -t92 * pkin(3) - t93 * pkin(8) - t318;
t236 = t179 * t18 - t183 * t33;
t4 = -t91 * pkin(4) + qJD(4) * t23 + t236;
t311 = t102 * (-pkin(4) * t190 + t102 * pkin(9)) + t201 + t4;
t169 = -pkin(2) * t176 - pkin(3);
t142 = -pkin(4) * t183 - pkin(9) * t179 + t169;
t38 = -qJD(4) * t190 - t183 * t162 + t179 * t93;
t36 = qJDD(5) + t38;
t117 = -qJD(1) * t227 + t160;
t64 = t174 * t117 + t277;
t310 = t102 * (-t64 + t124 * (pkin(4) * t179 - pkin(9) * t183)) + t142 * t36;
t15 = pkin(9) * t124 + t23;
t110 = t174 * t118;
t50 = t176 * t106 - t110;
t43 = -pkin(3) * t231 - t50;
t19 = t103 * pkin(4) + pkin(9) * t190 + t43;
t219 = t15 * t178 - t182 * t19;
t205 = t179 * t33 + t183 * t18 + t61 * t262 - t263 * t44;
t3 = pkin(9) * t91 + t205;
t20 = -t174 * t58 + t176 * t48;
t17 = -pkin(3) * t162 - t20;
t6 = pkin(4) * t38 - pkin(9) * t37 + t17;
t1 = -t219 * qJD(5) + t178 * t6 + t182 * t3;
t171 = t175 ^ 2;
t309 = 0.2e1 * t171;
t186 = qJD(1) ^ 2;
t308 = pkin(1) * t171;
t305 = g(3) * t184;
t116 = t165 + t196;
t270 = pkin(7) * t279 + t164;
t125 = qJ(3) * t279 + t270;
t72 = t174 * t116 + t176 * t125;
t63 = pkin(8) * t177 + t72;
t281 = t175 * t180;
t133 = t174 * t281 - t251;
t78 = t133 * pkin(3) - t134 * pkin(8) - t225;
t215 = t179 * t78 + t183 * t63;
t65 = t117 * t176 - t110;
t74 = pkin(2) * t249 + pkin(3) * t131 - pkin(8) * t128;
t303 = t179 * t74 + t183 * t65;
t302 = t102 * t55;
t260 = qJD(5) * t182;
t11 = t124 * t260 + t178 * t91 + t182 * t37 + t190 * t261;
t301 = t11 * t178;
t300 = t128 * t57;
t298 = t178 * t36;
t296 = t179 * t91;
t295 = t182 * t36;
t293 = t183 * t12;
t292 = t57 * t102;
t290 = t103 * t124;
t289 = t103 * t131;
t288 = t190 * t124;
t287 = t190 * t131;
t286 = t124 * t179;
t168 = pkin(2) * t174 + pkin(8);
t284 = t168 * t178;
t283 = t168 * t182;
t282 = t171 * t186;
t172 = t180 ^ 2;
t269 = -t184 ^ 2 + t172;
t266 = qJD(2) * t180;
t264 = qJD(4) * t168;
t252 = t180 * t282;
t248 = t175 * t266;
t247 = t175 * t177 * t186;
t246 = t304 * t175;
t241 = t175 * t257;
t238 = -t11 * t183 + t57 * t263;
t161 = t184 * t254;
t107 = t175 * t193 + t161;
t108 = -t175 * t265 + (-t184 * t246 - t164) * qJD(2);
t45 = t107 * t174 - t176 * t108;
t71 = t116 * t176 - t174 * t125;
t233 = t124 * t183;
t232 = t102 * t182;
t230 = qJD(2) + 0.2e1 * t268;
t229 = t162 + t258;
t223 = g(1) * t185 + g(2) * t181;
t222 = g(1) * t181 - g(2) * t185;
t8 = t15 * t182 + t178 * t19;
t29 = pkin(9) * t133 + t215;
t114 = t134 * t183 + t177 * t179;
t62 = -pkin(3) * t177 - t71;
t32 = pkin(4) * t113 - pkin(9) * t114 + t62;
t218 = t178 * t32 + t182 * t29;
t217 = -t178 * t29 + t182 * t32;
t22 = -t179 * t44 + t183 * t61;
t46 = t107 * t176 + t108 * t174;
t129 = t175 * t197;
t130 = t143 * t175 * qJD(2);
t75 = pkin(2) * t248 + pkin(3) * t129 + pkin(8) * t130;
t216 = -t179 * t46 + t183 * t75;
t214 = -t179 * t63 + t183 * t78;
t77 = t114 * t182 + t133 * t178;
t76 = t114 * t178 - t133 * t182;
t210 = -t124 * t263 + t128 * t286 + t183 * t91;
t79 = -t182 * t131 + t178 * t285;
t208 = (-t178 * t262 + t79) * t102;
t207 = -t102 * t260 - t298;
t204 = t179 * t75 + t183 * t46 + t78 * t262 - t263 * t63;
t203 = t124 * t43 - t168 * t91;
t98 = t181 * t199 - t185 * t211;
t200 = g(1) * t98 + g(2) * t95 - g(3) * t133;
t194 = -t17 - t200;
t192 = t270 * t177;
t14 = -pkin(4) * t124 - t22;
t189 = -pkin(9) * t36 + (t14 + t22) * t102;
t2 = -qJD(5) * t8 - t178 * t3 + t182 * t6;
t188 = qJD(5) * t102 * t168 + t200;
t187 = -g(1) * t212 - g(2) * t213 - g(3) * t134 + (pkin(9) * t131 - qJD(5) * t142 + t303) * t102;
t141 = -t177 * t275 + t271;
t139 = -t177 * t274 - t273;
t136 = pkin(2) * t276 - t246;
t85 = t179 * t280 + t183 * t212;
t70 = -qJD(4) * t113 - t130 * t183;
t69 = qJD(4) * t114 - t130 * t179;
t40 = -t178 * t98 + t182 * t85;
t39 = -t178 * t85 - t182 * t98;
t28 = -pkin(4) * t133 - t214;
t27 = -qJD(5) * t76 + t129 * t178 + t70 * t182;
t26 = qJD(5) * t77 - t129 * t182 + t70 * t178;
t24 = -pkin(4) * t131 + t179 * t65 - t183 * t74;
t13 = pkin(4) * t69 - pkin(9) * t70 + t45;
t10 = -t129 * pkin(4) + qJD(4) * t215 - t216;
t9 = pkin(9) * t129 + t204;
t5 = [qJDD(1), t222, t223, (qJDD(1) * t172 + 0.2e1 * t184 * t243) * t171, (t180 * t256 - t259 * t269) * t309, (qJD(2) * t184 * t230 + t180 * t229) * t175, (t184 * t229 - t230 * t266) * t175, t162 * t177, t253 * t309 + (-pkin(7) * t281 + t165) * t162 + (-pkin(7) * t241 + t159) * t177 - g(1) * t139 - g(2) * t141 - t270 * qJD(2) ^ 2 + 0.2e1 * (-t180 * t308 - t192) * t259, -(-pkin(7) * t248 + t161) * t231 - t270 * t162 - (-pkin(7) * t226 + t250) * t177 + g(1) * t206 - g(2) * t140 + 0.2e1 * (-t244 - t257) * t308, t128 * t46 - t129 * t51 + t130 * t50 + t131 * t45 - t133 * t21 - t134 * t20 - t175 * t223 - t71 * t93 + t72 * t92, t21 * t72 + t51 * t46 + t20 * t71 - t50 * t45 - g(1) * (-t136 * t185 - t170 * t181) - g(2) * (-t136 * t181 + t170 * t185) + (pkin(2) * t137 * t266 + t318 * t170) * t175, t114 * t37 - t190 * t70, -t103 * t70 - t113 * t37 - t114 * t38 + t190 * t69, t114 * t91 + t124 * t70 - t129 * t190 + t133 * t37, -t103 * t129 - t113 * t91 - t124 * t69 - t133 * t38, t124 * t129 + t133 * t91, t216 * t124 + t214 * t91 - t236 * t133 + t22 * t129 + t45 * t103 + t62 * t38 + t17 * t113 + t43 * t69 + g(1) * t82 - g(2) * t85 + (-t124 * t215 - t133 * t23) * qJD(4), -g(1) * t209 - g(2) * t84 + t17 * t114 - t204 * t124 - t23 * t129 - t205 * t133 - t190 * t45 - t215 * t91 + t62 * t37 + t43 * t70, t11 * t77 + t27 * t57, -t11 * t76 - t12 * t77 - t26 * t57 - t27 * t55, t102 * t27 + t11 * t113 + t36 * t77 + t57 * t69, -t102 * t26 - t113 * t12 - t36 * t76 - t55 * t69, t102 * t69 + t113 * t36, (-qJD(5) * t218 + t182 * t13 - t178 * t9) * t102 + t217 * t36 + t2 * t113 - t219 * t69 + t10 * t55 + t28 * t12 + t4 * t76 + t14 * t26 + g(1) * t315 - g(2) * t40, -(qJD(5) * t217 + t178 * t13 + t182 * t9) * t102 - t218 * t36 - t1 * t113 - t8 * t69 + t10 * t57 + t28 * t11 + t4 * t77 + t14 * t27 - g(1) * t316 - g(2) * t39; 0, 0, 0, -t184 * t252, t269 * t282, -t184 * t247 + t241, t180 * t247 + t240, t162, pkin(1) * t252 + t159 + (-pkin(7) * t257 - t305) * t175 + t186 * t192 + t317, t282 * t307 + (-pkin(7) * t249 + t160) * t268 + g(1) * t141 - g(2) * t139 + g(3) * t281 + t160 * qJD(2) - t250, (t51 - t64) * t131 + (t50 - t65) * t128 + (t174 * t92 - t176 * t93) * pkin(2), t50 * t64 - t51 * t65 + (t21 * t174 + t20 * t176 + (-t137 * t267 - t305) * t175 + t317) * pkin(2), t37 * t179 - t190 * t233, (t37 - t290) * t183 + (-t38 + t288) * t179, t124 * t233 + t287 + t296, t210 + t289, -t124 * t131, -t64 * t103 - t22 * t131 + t169 * t38 + (t65 * t124 + t203) * t179 + ((-t74 - t264) * t124 + t194) * t183, t169 * t37 + t303 * t124 + t23 * t131 + t64 * t190 + t203 * t183 + (t124 * t264 - t194) * t179, t11 * t182 * t179 + t191 * t57, t55 * t80 + t57 * t79 + (-t178 * t57 - t182 * t55) * t262 + (-t301 - t12 * t182 + (t178 * t55 - t182 * t57) * qJD(5)) * t179, (t295 - t300) * t179 + t312 + t238, t293 + t208 + (t207 - t314) * t179, t102 * t286 - t36 * t183, -t14 * t79 - t24 * t55 + t310 * t182 + t187 * t178 + (-t36 * t284 - t2 + (t14 * t178 + t168 * t55) * qJD(4) - t188 * t182) * t183 + (t14 * t260 + t168 * t12 + t219 * t128 + t4 * t178 + (t102 * t284 - t219) * qJD(4)) * t179, -t14 * t80 - t24 * t57 - t310 * t178 + t187 * t182 + (-t36 * t283 + t1 + (t14 * t182 + t168 * t57) * qJD(4) + t188 * t178) * t183 + (-t14 * t261 + t168 * t11 + t8 * t128 + t4 * t182 + (t102 * t283 - t8) * qJD(4)) * t179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128 ^ 2 - t131 ^ 2, -g(3) * t177 - t128 * t51 + t131 * t50 + (-t220 - t222) * t175 + t255, 0, 0, 0, 0, 0, t210 - t289, -t124 ^ 2 * t183 + t287 - t296, 0, 0, 0, 0, 0, -t293 + t208 + (t207 + t314) * t179, (-t295 - t300) * t179 - t312 + t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t190 * t103, -t103 ^ 2 + t190 ^ 2, t37 + t290, -t38 - t288, t91, t190 * t43 - t201 - t236 + (-qJD(4) + t124) * t23, g(1) * t85 + g(2) * t82 + g(3) * t114 + t103 * t43 + t124 * t22 - t205, t232 * t57 + t301, (t11 - t302) * t182 + (-t12 - t292) * t178, t102 * t232 + t190 * t57 + t298, -t102 ^ 2 * t178 - t190 * t55 + t295, t102 * t190, -pkin(4) * t12 + t189 * t178 - t311 * t182 - t190 * t219 - t23 * t55, -pkin(4) * t11 + t311 * t178 + t189 * t182 - t190 * t8 - t23 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 * t55, -t55 ^ 2 + t57 ^ 2, t11 + t302, -t12 + t292, t36, -g(1) * t39 + g(2) * t316 + g(3) * t76 + t8 * t102 - t14 * t57 + t2, g(1) * t40 + g(2) * t315 + g(3) * t77 - t219 * t102 + t14 * t55 - t1;];
tau_reg = t5;
