% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% tau_reg [6x31]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:36:23
% EndTime: 2019-03-09 02:36:31
% DurationCPUTime: 3.32s
% Computational Cost: add. (4338->406), mult. (9065->532), div. (0->0), fcn. (6976->14), ass. (0->213)
t168 = sin(pkin(10));
t169 = cos(pkin(10));
t173 = sin(qJ(4));
t177 = cos(qJ(4));
t119 = t168 * t177 + t169 * t173;
t111 = t119 * qJD(1);
t296 = qJD(5) + t111;
t238 = qJD(1) * t168;
t217 = t173 * t238;
t243 = t177 * t169;
t222 = qJD(1) * t243;
t113 = -t217 + t222;
t172 = sin(qJ(5));
t176 = cos(qJ(5));
t230 = t176 * qJD(4);
t86 = t113 * t172 - t230;
t301 = t296 * t86;
t163 = pkin(10) + qJ(4);
t148 = sin(t163);
t149 = cos(t163);
t174 = sin(qJ(1));
t178 = cos(qJ(1));
t285 = g(1) * t174 - g(2) * t178;
t184 = -g(3) * t148 + t149 * t285;
t170 = -pkin(1) - qJ(3);
t134 = t170 * qJD(1) + qJD(2);
t214 = -pkin(7) * qJD(1) + t134;
t103 = t214 * t168;
t104 = t214 * t169;
t236 = qJD(4) * t177;
t237 = qJD(4) * t173;
t278 = -qJD(1) * qJD(3) + qJDD(1) * t170;
t122 = qJDD(2) + t278;
t210 = -pkin(7) * qJDD(1) + t122;
t91 = t210 * t168;
t92 = t210 * t169;
t196 = -t103 * t236 - t104 * t237 - t173 * t91 + t177 * t92;
t22 = -qJDD(4) * pkin(4) - t196;
t300 = pkin(8) * qJD(5) * t296 + t184 + t22;
t175 = cos(qJ(6));
t171 = sin(qJ(6));
t249 = t171 * t176;
t121 = t172 * t175 + t249;
t280 = qJD(5) + qJD(6);
t81 = t280 * t121;
t269 = t121 * t111 + t81;
t88 = qJD(4) * t172 + t113 * t176;
t198 = t171 * t86 - t175 * t88;
t38 = t171 * t88 + t175 * t86;
t299 = t198 * t38;
t235 = qJD(5) * t172;
t248 = t172 * t111;
t298 = t235 + t248;
t209 = t176 * t296;
t185 = -qJD(4) * t217 + t119 * qJDD(1);
t221 = t169 * t236;
t76 = qJD(1) * t221 + t185;
t72 = qJDD(5) + t76;
t261 = t172 * t72;
t297 = -t209 * t296 - t261;
t295 = t198 ^ 2 - t38 ^ 2;
t102 = qJD(6) + t296;
t232 = qJD(6) * t175;
t233 = qJD(6) * t171;
t226 = t169 * qJDD(1);
t227 = t168 * qJDD(1);
t199 = -t173 * t227 + t177 * t226;
t281 = qJD(4) * t111;
t75 = t199 - t281;
t34 = qJD(5) * t230 + t172 * qJDD(4) - t113 * t235 + t176 * t75;
t212 = -t176 * qJDD(4) + t172 * t75;
t35 = qJD(5) * t88 + t212;
t8 = -t171 * t35 + t175 * t34 - t86 * t232 - t233 * t88;
t294 = t102 * t38 + t8;
t57 = t177 * t103 + t173 * t104;
t52 = qJD(4) * pkin(8) + t57;
t147 = qJD(1) * qJ(2) + qJD(3);
t127 = pkin(3) * t238 + t147;
t53 = pkin(4) * t111 - pkin(8) * t113 + t127;
t24 = t172 * t53 + t176 * t52;
t14 = -pkin(9) * t86 + t24;
t12 = t14 * t233;
t167 = qJ(5) + qJ(6);
t154 = cos(t167);
t274 = g(3) * t149;
t289 = -t173 * t103 + t104 * t177;
t51 = -qJD(4) * pkin(4) - t289;
t30 = pkin(5) * t86 + t51;
t251 = t154 * t174;
t153 = sin(t167);
t252 = t153 * t178;
t97 = t148 * t251 + t252;
t250 = t154 * t178;
t253 = t153 * t174;
t99 = t148 * t250 - t253;
t293 = g(1) * t97 - g(2) * t99 + t154 * t274 + t30 * t38 + t12;
t197 = t173 * t92 + t177 * t91;
t21 = qJDD(4) * pkin(8) + qJD(4) * t289 + t197;
t164 = qJDD(1) * qJ(2);
t165 = qJD(1) * qJD(2);
t284 = t164 + t165;
t130 = qJDD(3) + t284;
t117 = pkin(3) * t227 + t130;
t29 = pkin(4) * t76 - pkin(8) * t75 + t117;
t27 = t176 * t29;
t2 = pkin(5) * t72 - pkin(9) * t34 - qJD(5) * t24 - t172 * t21 + t27;
t234 = qJD(5) * t176;
t191 = t172 * t29 + t176 * t21 + t53 * t234 - t235 * t52;
t3 = -pkin(9) * t35 + t191;
t224 = -t171 * t3 + t175 * t2;
t23 = -t172 * t52 + t176 * t53;
t13 = -pkin(9) * t88 + t23;
t11 = pkin(5) * t296 + t13;
t263 = t14 * t175;
t5 = t171 * t11 + t263;
t96 = -t148 * t253 + t250;
t98 = t148 * t252 + t251;
t292 = -g(1) * t96 - g(2) * t98 - qJD(6) * t5 + t153 * t274 + t30 * t198 + t224;
t182 = qJD(6) * t198 - t171 * t34 - t175 * t35;
t291 = -t102 * t198 + t182;
t239 = t168 ^ 2 + t169 ^ 2;
t287 = t134 * t239;
t203 = g(1) * t178 + g(2) * t174;
t286 = t130 - t203;
t283 = t171 * t235 + t172 * t233;
t282 = -qJD(6) * t176 - t234;
t69 = qJDD(6) + t72;
t264 = t121 * t69;
t120 = t171 * t172 - t175 * t176;
t270 = (-t111 - t280) * t120;
t279 = -t270 * t102 - t264;
t277 = 0.2e1 * t165;
t276 = pkin(8) + pkin(9);
t273 = -pkin(7) + t170;
t74 = pkin(4) * t113 + pkin(8) * t111;
t272 = t172 * t74 + t176 * t289;
t118 = t168 * t173 - t243;
t140 = t168 * pkin(3) + qJ(2);
t73 = pkin(4) * t119 + pkin(8) * t118 + t140;
t125 = t273 * t168;
t126 = t273 * t169;
t79 = t125 * t177 + t126 * t173;
t77 = t176 * t79;
t271 = t172 * t73 + t77;
t268 = t113 * t38;
t267 = t113 * t198;
t266 = t113 * t86;
t265 = t113 * t88;
t44 = t120 * t69;
t262 = t172 * t34;
t62 = t176 * t72;
t260 = pkin(1) * qJDD(1);
t115 = -t168 * t237 + t221;
t259 = t102 * t115;
t114 = -t168 * t236 - t169 * t237;
t257 = t114 * t172;
t256 = t114 * t176;
t255 = t118 * t172;
t254 = t118 * t176;
t247 = t172 * t174;
t246 = t172 * t178;
t245 = t174 * t176;
t244 = t176 * t178;
t242 = t114 * qJD(4) - t118 * qJDD(4);
t241 = t178 * pkin(1) + t174 * qJ(2);
t223 = qJD(5) * t276;
t220 = t118 * t234;
t216 = qJD(6) * t11 + t3;
t213 = -qJD(5) * t53 - t21;
t211 = t239 * t122;
t208 = qJDD(2) - t260;
t207 = qJD(5) * t119 + qJD(1);
t206 = -t269 * t102 - t44;
t205 = t298 * pkin(5) - t57;
t204 = -t234 * t52 + t27;
t131 = t276 * t172;
t201 = pkin(9) * t248 + qJD(6) * t131 + t172 * t223 + t272;
t132 = t276 * t176;
t65 = t176 * t74;
t200 = pkin(5) * t113 + qJD(6) * t132 - t172 * t289 + t65 + (pkin(9) * t111 + t223) * t176;
t78 = t125 * t173 - t126 * t177;
t194 = -t296 * t298 + t62;
t193 = -qJD(4) * t115 - qJDD(4) * t119;
t47 = -qJD(3) * t119 - qJD(4) * t78;
t70 = pkin(4) * t115 - pkin(8) * t114 + qJD(2);
t190 = t172 * t70 + t176 * t47 + t73 * t234 - t235 * t79;
t189 = -pkin(8) * t72 + t296 * t51;
t188 = t220 - t257;
t187 = t102 * t120;
t181 = t284 + t286;
t48 = -qJD(3) * t118 + qJD(4) * t79;
t179 = qJD(1) ^ 2;
t157 = t178 * qJ(2);
t146 = -pkin(5) * t176 - pkin(4);
t110 = t148 * t244 - t247;
t109 = t148 * t246 + t245;
t108 = t148 * t245 + t246;
t107 = -t148 * t247 + t244;
t67 = t120 * t118;
t66 = t121 * t118;
t63 = t176 * t73;
t60 = t176 * t70;
t49 = -pkin(5) * t255 + t78;
t28 = -pkin(5) * t188 + t48;
t25 = pkin(9) * t255 + t271;
t19 = pkin(5) * t119 + pkin(9) * t254 - t172 * t79 + t63;
t16 = t114 * t249 + (-t280 * t254 + t257) * t175 + t283 * t118;
t15 = -t114 * t120 + t118 * t81;
t10 = pkin(5) * t35 + t22;
t7 = pkin(9) * t188 + t190;
t6 = -pkin(9) * t256 + pkin(5) * t115 - t172 * t47 + t60 + (-t77 + (-pkin(9) * t118 - t73) * t172) * qJD(5);
t4 = t175 * t11 - t14 * t171;
t1 = [qJDD(1), t285, t203, qJDD(2) - t285 - 0.2e1 * t260, 0.2e1 * t164 + t277 - t203, -t208 * pkin(1) - g(1) * (-pkin(1) * t174 + t157) - g(2) * t241 + (t164 + t277) * qJ(2), t181 * t168, t181 * t169, t285 + t239 * (-t122 - t278) t130 * qJ(2) + t147 * qJD(2) - g(1) * (t170 * t174 + t157) - g(2) * (qJ(3) * t178 + t241) + t170 * t211 - qJD(3) * t287, t113 * t114 - t118 * t75, -t111 * t114 - t113 * t115 + t118 * t76 - t119 * t75, t242, t193, 0, qJD(2) * t111 - qJD(4) * t48 - qJDD(4) * t78 + t115 * t127 + t117 * t119 + t140 * t76 - t148 * t203, qJD(2) * t113 - qJD(4) * t47 - qJDD(4) * t79 + t114 * t127 - t117 * t118 + t140 * t75 - t149 * t203, t88 * t256 + (-t176 * t34 + t235 * t88) * t118 (-t172 * t88 - t176 * t86) * t114 + (t262 + t176 * t35 + (-t172 * t86 + t176 * t88) * qJD(5)) * t118, -t72 * t254 + t115 * t88 + t119 * t34 + (t118 * t235 + t256) * t296, -t115 * t86 - t119 * t35 + t188 * t296 + t72 * t255, t115 * t296 + t119 * t72 (-t234 * t79 + t60) * t296 + t63 * t72 + t204 * t119 + t23 * t115 + t48 * t86 + t78 * t35 - t51 * t220 - g(1) * t110 - g(2) * t108 + ((-qJD(5) * t73 - t47) * t296 - t79 * t72 + t213 * t119 - t22 * t118 + t51 * t114) * t172, -t190 * t296 - t271 * t72 - t191 * t119 - t24 * t115 + t48 * t88 + t78 * t34 + t51 * t256 + g(1) * t109 - g(2) * t107 + (-t22 * t176 + t235 * t51) * t118, -t15 * t198 + t67 * t8, -t15 * t38 + t16 * t198 + t182 * t67 + t66 * t8, t102 * t15 - t115 * t198 + t119 * t8 + t67 * t69, -t102 * t16 - t115 * t38 + t119 * t182 + t66 * t69, t119 * t69 + t259 (-t171 * t7 + t175 * t6) * t102 + (-t171 * t25 + t175 * t19) * t69 + t224 * t119 + t4 * t115 + t28 * t38 - t49 * t182 - t10 * t66 + t30 * t16 - g(1) * t99 - g(2) * t97 + ((-t171 * t19 - t175 * t25) * t102 - t5 * t119) * qJD(6), g(1) * t98 - g(2) * t96 + t10 * t67 - t5 * t115 + t12 * t119 + t30 * t15 - t28 * t198 + t49 * t8 + (-(-qJD(6) * t25 + t6) * t102 - t19 * t69 - t2 * t119) * t171 + (-(qJD(6) * t19 + t7) * t102 - t25 * t69 - t216 * t119) * t175; 0, 0, 0, qJDD(1), -t179, -qJ(2) * t179 + t208 - t285, -t179 * t168, -t179 * t169, -t239 * qJDD(1), -qJD(1) * t147 + t211 - t285, 0, 0, 0, 0, 0, -qJD(1) * t111 + t242, -qJD(1) * t113 + t193, 0, 0, 0, 0, 0, -t119 * t261 - t114 * t86 + t118 * t35 + (-t115 * t172 - t176 * t207) * t296, -t119 * t62 - t114 * t88 + t118 * t34 + (-t115 * t176 + t172 * t207) * t296, 0, 0, 0, 0, 0, -t114 * t38 - t118 * t182 - t121 * t259 + qJD(1) * t187 + ((t282 * t175 + t283) * t102 - t264) * t119, t114 * t198 + t118 * t8 + t115 * t187 + t121 * t102 * qJD(1) + (-(t282 * t171 - t172 * t232 - t175 * t235) * t102 + t44) * t119; 0, 0, 0, 0, 0, 0, t227, t226, -t239 * t179, qJD(1) * t287 + t286, 0, 0, 0, 0, 0 (t113 + t222) * qJD(4) + t185, t199 - 0.2e1 * t281, 0, 0, 0, 0, 0, t194 - t266, -t265 + t297, 0, 0, 0, 0, 0, t206 - t268, t267 + t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113 * t111, -t111 ^ 2 + t113 ^ 2, t199 (t113 - t222) * qJD(4) - t185, qJDD(4), qJD(4) * t57 - t113 * t127 - t184 + t196, t111 * t127 + t148 * t285 - t197 + t274, t209 * t88 + t262 (t34 - t301) * t176 + (-t296 * t88 - t35) * t172, -t265 - t297, t194 + t266, -t296 * t113, -pkin(4) * t35 - t65 * t296 - t23 * t113 - t57 * t86 + (t289 * t296 + t189) * t172 - t300 * t176, -pkin(4) * t34 + t24 * t113 + t172 * t300 + t189 * t176 + t272 * t296 - t57 * t88, t121 * t8 - t198 * t270, -t120 * t8 + t121 * t182 + t198 * t269 - t270 * t38, t267 - t279, t206 + t268, -t102 * t113 (-t131 * t175 - t132 * t171) * t69 - t146 * t182 + t10 * t120 - t4 * t113 + t205 * t38 + t269 * t30 + (t171 * t201 - t175 * t200) * t102 - t184 * t154 -(-t131 * t171 + t132 * t175) * t69 + t146 * t8 + t10 * t121 + t5 * t113 - t205 * t198 + t270 * t30 + (t171 * t200 + t175 * t201) * t102 + t184 * t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88 * t86, -t86 ^ 2 + t88 ^ 2, t34 + t301, -t212 + (-qJD(5) + t296) * t88, t72, -g(1) * t107 - g(2) * t109 + t296 * t24 - t51 * t88 + (t213 + t274) * t172 + t204, g(1) * t108 - g(2) * t110 + t176 * t274 + t23 * t296 + t51 * t86 - t191, -t299, t295, t294, t291, t69 -(-t13 * t171 - t263) * t102 + (-t102 * t233 + t175 * t69 - t88 * t38) * pkin(5) + t292 (-t14 * t102 - t2) * t171 + (t13 * t102 - t216) * t175 + (-t102 * t232 - t171 * t69 + t198 * t88) * pkin(5) + t293; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299, t295, t294, t291, t69, t102 * t5 + t292, t102 * t4 - t171 * t2 - t175 * t216 + t293;];
tau_reg  = t1;
