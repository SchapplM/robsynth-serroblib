% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRRR7
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:13:05
% EndTime: 2019-12-05 17:13:13
% DurationCPUTime: 3.42s
% Computational Cost: add. (4652->385), mult. (10588->522), div. (0->0), fcn. (7730->14), ass. (0->192)
t172 = sin(qJ(3));
t280 = pkin(7) + pkin(6);
t227 = qJD(3) * t280;
t129 = t172 * t227;
t176 = cos(qJ(3));
t130 = t176 * t227;
t135 = t280 * t172;
t136 = t280 * t176;
t171 = sin(qJ(4));
t175 = cos(qJ(4));
t251 = t175 * t176;
t253 = t171 * t172;
t123 = -t251 + t253;
t177 = cos(qJ(2));
t197 = t123 * t177;
t240 = qJD(4) * t175;
t241 = qJD(4) * t171;
t266 = qJD(1) * t197 - t175 * t129 - t171 * t130 - t135 * t240 - t136 * t241;
t124 = t171 * t176 + t172 * t175;
t198 = t124 * t177;
t82 = -t171 * t135 + t175 * t136;
t265 = qJD(1) * t198 - t82 * qJD(4) + t171 * t129 - t175 * t130;
t230 = t176 * qJDD(2);
t231 = t172 * qJDD(2);
t163 = qJD(3) + qJD(4);
t286 = t163 * t124;
t48 = qJD(2) * t286 + t171 * t231 - t175 * t230;
t289 = -pkin(8) * t286 + t266;
t207 = t163 * t253;
t242 = qJD(3) * t176;
t77 = -t175 * t242 - t176 * t240 + t207;
t288 = t77 * pkin(8) + t265;
t225 = qJD(2) * t251;
t245 = qJD(2) * t172;
t226 = t171 * t245;
t113 = -t225 + t226;
t115 = t124 * qJD(2);
t170 = sin(qJ(5));
t174 = cos(qJ(5));
t205 = t113 * t170 - t174 * t115;
t64 = t174 * t113 + t115 * t170;
t287 = t64 * t205;
t168 = sin(pkin(9));
t169 = cos(pkin(9));
t210 = g(1) * t169 + g(2) * t168;
t201 = t210 * t177;
t173 = sin(qJ(2));
t271 = g(3) * t173;
t193 = t201 + t271;
t202 = t210 * t173;
t270 = g(3) * t177;
t192 = t202 - t270;
t234 = qJD(2) * qJD(3);
t222 = t176 * t234;
t285 = -t222 - t231;
t246 = qJD(1) * t173;
t141 = qJD(2) * pkin(6) + t246;
t236 = t177 * qJD(1);
t262 = qJD(2) * pkin(2);
t142 = -t236 - t262;
t165 = t172 ^ 2;
t166 = t176 ^ 2;
t248 = t165 + t166;
t216 = t248 * t177;
t284 = t141 * t216 + t142 * t173;
t12 = t205 ^ 2 - t64 ^ 2;
t238 = qJD(5) * t174;
t239 = qJD(5) * t170;
t212 = -qJD(4) * t225 - t171 * t230 + t285 * t175;
t47 = qJD(2) * t207 + t212;
t13 = t113 * t238 + t115 * t239 + t170 * t48 + t174 * t47;
t158 = qJD(5) + t163;
t8 = t158 * t64 - t13;
t109 = t115 * pkin(8);
t218 = pkin(7) * qJD(2) + t141;
t105 = t218 * t176;
t90 = t171 * t105;
t104 = t218 * t172;
t93 = qJD(3) * pkin(3) - t104;
t54 = t175 * t93 - t90;
t34 = -t109 + t54;
t32 = pkin(4) * t163 + t34;
t275 = pkin(8) * t113;
t92 = t175 * t105;
t55 = t171 * t93 + t92;
t35 = t55 - t275;
t235 = qJD(1) * qJD(2);
t120 = qJDD(2) * pkin(6) + t173 * qJDD(1) + t177 * t235;
t53 = qJDD(3) * pkin(3) + t285 * pkin(7) - t172 * t120 - t141 * t242;
t223 = t172 * t234;
t243 = qJD(3) * t172;
t56 = -t141 * t243 + t176 * t120 + (-t223 + t230) * pkin(7);
t16 = -qJD(4) * t55 - t171 * t56 + t175 * t53;
t162 = qJDD(3) + qJDD(4);
t6 = pkin(4) * t162 + pkin(8) * t47 + t16;
t219 = t105 * t241 - t171 * t53 - t175 * t56 - t93 * t240;
t7 = -pkin(8) * t48 - t219;
t1 = t174 * (qJD(5) * t32 + t7) + t170 * t6 - t35 * t239;
t167 = qJ(3) + qJ(4);
t161 = qJ(5) + t167;
t153 = sin(t161);
t154 = cos(t161);
t255 = t169 * t177;
t256 = t168 * t177;
t276 = pkin(3) * t176;
t156 = pkin(2) + t276;
t121 = -qJD(2) * t156 - t236;
t76 = t113 * pkin(4) + t121;
t183 = t64 * t76 - g(1) * (-t153 * t168 - t154 * t255) - g(2) * (t153 * t169 - t154 * t256) + t154 * t271 - t1;
t260 = t174 * t35;
t18 = t170 * t32 + t260;
t2 = -qJD(5) * t18 - t170 * t7 + t174 * t6;
t182 = t205 * t76 - g(1) * (-t153 * t255 + t154 * t168) - g(2) * (-t153 * t256 - t154 * t169) + t153 * t271 + t2;
t190 = qJD(5) * t205 + t170 * t47 - t174 * t48;
t9 = -t158 * t205 + t190;
t107 = t123 * t173;
t199 = pkin(3) * t243 - t246;
t159 = sin(t167);
t160 = cos(t167);
t282 = t159 * t271 - g(1) * (-t159 * t255 + t160 * t168) - g(2) * (-t159 * t256 - t160 * t169);
t213 = -t177 * qJDD(1) + t173 * t235;
t259 = qJDD(2) * pkin(2);
t119 = t213 - t259;
t179 = qJD(3) ^ 2;
t281 = -pkin(6) * t179 + t173 * (t210 + t235) - t119 + t259 - t270;
t81 = -t175 * t135 - t136 * t171;
t61 = -pkin(8) * t124 + t81;
t62 = -pkin(8) * t123 + t82;
t27 = -t170 * t62 + t174 * t61;
t278 = qJD(5) * t27 + t288 * t170 + t289 * t174;
t28 = t170 * t61 + t174 * t62;
t277 = -qJD(5) * t28 - t289 * t170 + t288 * t174;
t155 = pkin(3) * t175 + pkin(4);
t254 = t170 * t171;
t57 = t104 * t171 - t92;
t38 = t57 + t275;
t58 = -t175 * t104 - t90;
t39 = -t109 + t58;
t264 = -t170 * t38 - t174 * t39 + t155 * t238 + (-t171 * t239 + (t174 * t175 - t254) * qJD(4)) * pkin(3);
t252 = t171 * t174;
t263 = t170 * t39 - t174 * t38 - t155 * t239 + (-t171 * t238 + (-t170 * t175 - t252) * qJD(4)) * pkin(3);
t261 = t170 * t35;
t258 = t115 * t113;
t250 = qJDD(1) - g(3);
t249 = t165 - t166;
t180 = qJD(2) ^ 2;
t247 = t179 + t180;
t244 = qJD(2) * t173;
t237 = t121 * qJD(2);
t233 = qJDD(2) * t177;
t232 = qJDD(3) * t172;
t228 = t172 * t180 * t176;
t132 = pkin(4) * t160 + t276;
t217 = t248 * t120;
t215 = pkin(4) * t286 + t199;
t214 = t248 * qJDD(2);
t211 = t172 * t222;
t209 = g(1) * t168 - g(2) * t169;
t17 = t174 * t32 - t261;
t208 = -t17 * t64 - t18 * t205;
t106 = t124 * t173;
t59 = -t106 * t174 + t107 * t170;
t60 = -t106 * t170 - t107 * t174;
t75 = -t123 * t170 + t124 * t174;
t200 = t209 * t176;
t79 = pkin(3) * t223 - qJDD(2) * t156 + t213;
t189 = -pkin(6) * qJDD(3) + (t142 + t236 - t262) * qJD(3);
t188 = -g(1) * (-t159 * t168 - t160 * t255) - g(2) * (t159 * t169 - t160 * t256) + t113 * t121 + t160 * t271 + t219;
t187 = -t142 * qJD(2) - t120 + t193;
t181 = -t115 * t121 + t16 + t282;
t164 = -pkin(8) - t280;
t157 = qJDD(5) + t162;
t131 = -pkin(3) * t172 - pkin(4) * t159;
t128 = pkin(2) + t132;
t111 = pkin(3) * t252 + t155 * t170;
t110 = -pkin(3) * t254 + t155 * t174;
t97 = pkin(4) * t123 - t156;
t84 = pkin(3) * t245 + pkin(4) * t115;
t74 = t174 * t123 + t124 * t170;
t51 = -t113 ^ 2 + t115 ^ 2;
t37 = -qJD(2) * t198 + t163 * t107;
t36 = -qJD(2) * t197 - t173 * t286;
t31 = t115 * t163 - t48;
t30 = -t212 + (t113 - t226) * t163;
t29 = pkin(4) * t48 + t79;
t24 = qJD(5) * t75 - t170 * t77 + t174 * t286;
t23 = t123 * t238 + t124 * t239 + t170 * t286 + t174 * t77;
t20 = t174 * t34 - t261;
t19 = -t170 * t34 - t260;
t11 = -qJD(5) * t60 - t170 * t36 + t174 * t37;
t10 = qJD(5) * t59 + t170 * t37 + t174 * t36;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t250, 0, 0, 0, 0, 0, 0, -t173 * t180 + t233, -qJDD(2) * t173 - t177 * t180, 0, -g(3) + (t173 ^ 2 + t177 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, (-0.2e1 * t223 + t230) * t177 + (-t176 * t247 - t232) * t173, (-qJDD(3) * t173 - 0.2e1 * t177 * t234) * t176 + (t173 * t247 - t233) * t172, t173 * t214 + t180 * t216, t284 * qJD(2) - t119 * t177 + t173 * t217 - g(3), 0, 0, 0, 0, 0, 0, -t106 * t162 + t113 * t244 + t163 * t37 - t177 * t48, t107 * t162 + t115 * t244 - t163 * t36 + t177 * t47, -t106 * t47 + t107 * t48 - t113 * t36 - t115 * t37, -t106 * t16 + t107 * t219 + t173 * t237 - t177 * t79 + t36 * t55 + t37 * t54 - g(3), 0, 0, 0, 0, 0, 0, t11 * t158 + t157 * t59 + t177 * t190 + t244 * t64, -t10 * t158 + t13 * t177 - t157 * t60 - t205 * t244, -t10 * t64 + t11 * t205 + t13 * t59 + t190 * t60, t1 * t60 + t10 * t18 + t11 * t17 - t177 * t29 + t2 * t59 + t244 * t76 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t250 * t177 + t202, -t173 * t250 + t201, 0, 0, qJDD(2) * t165 + 0.2e1 * t211, 0.2e1 * t172 * t230 - 0.2e1 * t234 * t249, t176 * t179 + t232, qJDD(2) * t166 - 0.2e1 * t211, qJDD(3) * t176 - t172 * t179, 0, t189 * t172 + t281 * t176, -t281 * t172 + t189 * t176, -t271 + t217 + pkin(6) * t214 + (-t235 * t248 - t210) * t177, (-t119 + t192) * pkin(2) + (t217 - t193) * pkin(6) - t284 * qJD(1), -t115 * t77 - t124 * t47, t113 * t77 - t115 * t286 + t123 * t47 - t124 * t48, t124 * t162 - t163 * t77, t113 * t286 + t123 * t48, -t123 * t162 - t163 * t286, 0, t113 * t199 + t121 * t286 + t123 * t79 - t156 * t48 + t160 * t192 + t162 * t81 + t163 * t265, t115 * t199 - t121 * t77 + t124 * t79 + t156 * t47 - t159 * t192 - t162 * t82 - t163 * t266, -t113 * t266 - t115 * t265 + t123 * t219 - t124 * t16 - t286 * t55 + t47 * t81 - t48 * t82 + t54 * t77 - t193, -t219 * t82 + t16 * t81 - t79 * t156 - g(3) * (t156 * t177 + t173 * t280) + t266 * t55 + t265 * t54 + t199 * t121 + t210 * (t156 * t173 - t177 * t280), -t13 * t75 + t205 * t23, t13 * t74 + t190 * t75 + t205 * t24 + t23 * t64, t157 * t75 - t158 * t23, -t190 * t74 + t24 * t64, -t157 * t74 - t158 * t24, 0, t192 * t154 + t157 * t27 + t277 * t158 - t190 * t97 + t215 * t64 + t24 * t76 + t29 * t74, -t13 * t97 - t153 * t192 - t157 * t28 - t278 * t158 - t205 * t215 - t23 * t76 + t29 * t75, -t1 * t74 + t13 * t27 + t17 * t23 - t18 * t24 + t190 * t28 - t2 * t75 + t205 * t277 - t278 * t64 - t193, t1 * t28 + t2 * t27 + t29 * t97 - g(3) * (t128 * t177 - t164 * t173) + t215 * t76 + t278 * t18 + t277 * t17 + t210 * (t128 * t173 + t164 * t177); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t228, t249 * t180, t231, t228, t230, qJDD(3), t172 * t187 - t200, t172 * t209 + t176 * t187, 0, 0, t258, t51, t30, -t258, t31, t162, -t163 * t57 + (-t113 * t245 + t162 * t175 - t163 * t241) * pkin(3) + t181, t163 * t58 + (-t115 * t245 - t162 * t171 - t163 * t240) * pkin(3) + t188, (t55 + t57) * t115 + (-t54 + t58) * t113 + (-t171 * t48 + t175 * t47 + (-t113 * t175 + t115 * t171) * qJD(4)) * pkin(3), -t54 * t57 - t55 * t58 + (-t219 * t171 + t16 * t175 - t200 + (t193 - t237) * t172 + (-t171 * t54 + t175 * t55) * qJD(4)) * pkin(3), -t287, t12, t8, t287, t9, t157, t110 * t157 + t158 * t263 - t64 * t84 + t182, -t111 * t157 - t158 * t264 + t205 * t84 + t183, t110 * t13 + t111 * t190 + t205 * t263 - t264 * t64 + t208, t1 * t111 + t2 * t110 - t76 * t84 - g(1) * (t131 * t255 + t132 * t168) - g(2) * (t131 * t256 - t132 * t169) - t131 * t271 + t264 * t18 + t263 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t258, t51, t30, -t258, t31, t162, t163 * t55 + t181, t163 * t54 + t188, 0, 0, -t287, t12, t8, t287, t9, t157, -t158 * t19 + (-t115 * t64 + t157 * t174 - t158 * t239) * pkin(4) + t182, t158 * t20 + (t115 * t205 - t157 * t170 - t158 * t238) * pkin(4) + t183, -t19 * t205 + t20 * t64 + (t13 * t174 + t190 * t170 + (-t170 * t205 - t174 * t64) * qJD(5)) * pkin(4) + t208, -t17 * t19 - t18 * t20 + (t1 * t170 + t2 * t174 - t76 * t115 + (-t17 * t170 + t174 * t18) * qJD(5) + t282) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t287, t12, t8, t287, t9, t157, t158 * t18 + t182, t158 * t17 + t183, 0, 0;];
tau_reg = t3;
