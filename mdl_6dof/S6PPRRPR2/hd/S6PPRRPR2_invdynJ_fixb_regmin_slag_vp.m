% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PPRRPR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% tau_reg [6x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PPRRPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:50:57
% EndTime: 2019-03-08 18:51:04
% DurationCPUTime: 3.19s
% Computational Cost: add. (2492->358), mult. (6437->511), div. (0->0), fcn. (6063->14), ass. (0->201)
t123 = sin(qJ(3));
t126 = cos(qJ(3));
t117 = sin(pkin(7));
t120 = cos(pkin(7));
t116 = sin(pkin(12));
t119 = cos(pkin(12));
t245 = cos(pkin(11));
t246 = cos(pkin(6));
t182 = t246 * t245;
t244 = sin(pkin(11));
t141 = t244 * t116 - t119 * t182;
t118 = sin(pkin(6));
t194 = t118 * t245;
t273 = t117 * t194 + t120 * t141;
t72 = t116 * t182 + t244 * t119;
t32 = t72 * t123 + t126 * t273;
t181 = t246 * t244;
t142 = t245 * t116 + t119 * t181;
t193 = t118 * t244;
t274 = -t117 * t193 + t142 * t120;
t73 = -t116 * t181 + t245 * t119;
t34 = t73 * t123 + t126 * t274;
t195 = t117 * t246;
t185 = t126 * t195;
t231 = t119 * t120;
t208 = t126 * t231;
t234 = t116 * t123;
t209 = t118 * t234;
t54 = -t118 * t208 - t185 + t209;
t160 = g(1) * t34 + g(2) * t32 + g(3) * t54;
t102 = t246 * qJDD(1) + qJDD(2);
t105 = t246 * qJD(1) + qJD(2);
t227 = qJD(1) * t118;
t204 = t119 * t227;
t189 = t120 * t204;
t213 = qJDD(1) * t118;
t197 = t119 * t213;
t226 = qJD(3) * t123;
t203 = t117 * t226;
t205 = t116 * t227;
t224 = qJD(3) * t126;
t135 = -(t102 * t117 + t120 * t197) * t126 + qJDD(1) * t209 + t105 * t203 + t189 * t226 + t205 * t224;
t122 = sin(qJ(4));
t214 = qJD(3) * qJD(4);
t199 = t122 * t214;
t134 = pkin(4) * t199 + t135;
t216 = t122 * qJD(5);
t125 = cos(qJ(4));
t221 = qJD(4) * t125;
t162 = -qJ(5) * t221 - t216;
t196 = -t122 * qJ(5) - pkin(3);
t172 = t125 * pkin(4) - t196;
t279 = t172 * qJDD(3);
t9 = t162 * qJD(3) + t134 - t279;
t285 = -t9 + t160;
t210 = t125 * qJDD(3);
t282 = t199 - t210;
t164 = t116 * t126 + t123 * t231;
t154 = t164 * t118;
t233 = t117 * t123;
t46 = qJD(1) * t154 + t105 * t233;
t44 = qJD(3) * pkin(9) + t46;
t67 = t120 * t105 - t117 * t204;
t255 = -t122 * t44 + t125 * t67;
t280 = -qJD(5) + t255;
t217 = t122 * qJD(3);
t107 = qJD(6) + t217;
t284 = t107 - qJD(6);
t283 = t46 * qJD(3) - t135 + t160;
t26 = -qJD(4) * pkin(4) - t280;
t225 = qJD(3) * t125;
t281 = t172 * qJD(3);
t129 = qJD(3) ^ 2;
t173 = qJDD(3) * t126 - t123 * t129;
t202 = t117 * t224;
t76 = t122 * t120 + t125 * t233;
t57 = t76 * qJD(4) + t122 * t202;
t207 = t122 * t233;
t75 = -t125 * t120 + t207;
t278 = -t57 * qJD(4) - t75 * qJDD(4) + t117 * (t173 * t125 - t126 * t199);
t198 = t125 * t214;
t56 = qJD(4) * t207 - t120 * t221 - t125 * t202;
t277 = -t56 * qJD(4) + t76 * qJDD(4) + t117 * (t173 * t122 + t126 * t198);
t222 = qJD(4) * t122;
t276 = pkin(4) * t222 - t46;
t29 = t122 * t67 + t125 * t44;
t27 = -qJD(4) * qJ(5) - t29;
t121 = sin(qJ(6));
t220 = qJD(6) * t121;
t124 = cos(qJ(6));
t211 = t122 * qJDD(3);
t158 = t198 + t211;
t80 = qJDD(6) + t158;
t74 = t124 * t80;
t275 = -t107 * t220 + t74;
t230 = pkin(5) * t217 - t280;
t165 = t208 - t234;
t238 = qJDD(3) * pkin(9);
t20 = t238 + (t102 * t123 + t105 * t224) * t117 + (t165 * qJD(3) * qJD(1) + t164 * qJDD(1)) * t118;
t272 = t122 * t20 + t44 * t221 + t67 * t222;
t131 = t141 * t117 - t120 * t194;
t33 = -t123 * t273 + t72 * t126;
t12 = t33 * t122 - t131 * t125;
t132 = t142 * t117 + t120 * t193;
t35 = -t123 * t274 + t73 * t126;
t14 = t35 * t122 - t132 * t125;
t55 = t123 * t195 + t154;
t71 = -t118 * t119 * t117 + t246 * t120;
t36 = t55 * t122 - t71 * t125;
t271 = -g(1) * t14 - g(2) * t12 - g(3) * t36;
t108 = pkin(5) * t225;
t127 = -pkin(4) - pkin(10);
t23 = t108 - t27;
t270 = t127 * t80 + (t23 - t108 - t29) * t107;
t45 = -t123 * t205 + (t105 * t117 + t189) * t126;
t200 = t121 * t225;
t48 = -qJD(6) * t200 + t121 * qJDD(4) + (qJD(4) * qJD(6) - t282) * t124;
t254 = t162 + t276;
t128 = qJD(4) ^ 2;
t257 = pkin(9) * t128;
t267 = t254 * qJD(3) + t257 - t279 - t285;
t50 = (t165 * t118 + t185) * qJD(3);
t11 = -t55 * t222 + (qJD(4) * t71 + t50) * t125;
t37 = t71 * t122 + t55 * t125;
t51 = t55 * qJD(3);
t266 = (t122 * t51 + t54 * t221) * qJD(3) - t11 * qJD(4) - t37 * qJDD(4) + t54 * t211;
t10 = t37 * qJD(4) + t50 * t122;
t265 = (-t125 * t51 + t54 * t222) * qJD(3) - t10 * qJD(4) - t36 * qJDD(4) - t54 * t210;
t264 = pkin(5) + pkin(9);
t242 = qJ(5) * t125;
t180 = pkin(10) * t122 - t242;
t145 = t180 * qJD(4) - t216;
t256 = -t145 - t276;
t253 = qJD(3) * pkin(3);
t252 = t121 * t80;
t215 = t124 * qJD(4);
t83 = -t200 + t215;
t251 = t124 * t83;
t65 = t120 * t102 - t117 * t197;
t250 = t125 * t65;
t218 = t121 * qJD(4);
t81 = t124 * t225 + t218;
t47 = -t81 * qJD(6) + t124 * qJDD(4) + t282 * t121;
t249 = t47 * t124;
t248 = t81 * t107;
t247 = t83 * t107;
t243 = pkin(9) * qJDD(4);
t241 = qJD(4) * t81;
t240 = qJD(4) * t83;
t237 = qJDD(4) * pkin(4);
t236 = t107 * t121;
t235 = t107 * t124;
t232 = t117 * t126;
t114 = t122 ^ 2;
t115 = t125 ^ 2;
t229 = t114 - t115;
t219 = qJD(6) * t125;
t212 = qJDD(4) * qJ(5);
t97 = t264 * t125;
t206 = t122 * t129 * t125;
t192 = -t122 * t65 - t125 * t20 - t67 * t221 + t44 * t222;
t22 = t127 * qJD(4) + t230;
t79 = t127 * t125 + t196;
t30 = t79 * qJD(3) - t45;
t178 = t121 * t30 - t124 * t22;
t7 = t121 * t22 + t124 * t30;
t16 = -t54 * t121 + t36 * t124;
t17 = t36 * t121 + t54 * t124;
t168 = -t121 * t75 + t124 * t232;
t167 = t121 * t232 + t124 * t75;
t166 = -qJD(6) * t235 - t252;
t163 = qJDD(5) - t250 + t272;
t13 = t131 * t122 + t33 * t125;
t15 = t132 * t122 + t35 * t125;
t161 = -g(1) * t15 - g(2) * t13 - g(3) * t37;
t159 = g(1) * t35 + g(2) * t33 + g(3) * t55;
t8 = t145 * qJD(3) + t79 * qJDD(3) + t134;
t157 = t160 - t8;
t96 = t264 * t122;
t153 = -t96 * t80 + t159;
t5 = t163 - t237;
t43 = -t45 - t253;
t152 = -t243 + (t43 + t45 - t253) * qJD(4);
t38 = -t45 - t281;
t148 = t243 + (-t38 - t45 + t281) * qJD(4);
t147 = qJD(4) * qJD(5) - t192 + t212;
t146 = -g(1) * t193 + g(2) * t194 - g(3) * t246;
t144 = t161 - t192;
t140 = -t29 * qJD(4) + t271 + t272;
t110 = pkin(4) * t217;
t3 = -t282 * pkin(5) + t147;
t136 = t3 + (t180 * qJD(3) - qJD(6) * t127 + t110) * t107 + t161;
t133 = 0.2e1 * qJDD(3) * pkin(3) - t257 + t283;
t130 = t5 * t122 + t147 * t125 + (t122 * t27 + t125 * t26) * qJD(4) - t159;
t90 = qJD(4) * t97;
t89 = t264 * t222;
t84 = -qJ(5) * t225 + t110;
t31 = t38 * t217;
t2 = t158 * pkin(5) + t127 * qJDD(4) + t163;
t1 = t124 * t2;
t4 = [qJDD(1) - g(3), t102 * t246 - g(3) + (t116 ^ 2 + t119 ^ 2) * t118 ^ 2 * qJDD(1), 0, -t51 * qJD(3) - t54 * qJDD(3), -t50 * qJD(3) - t55 * qJDD(3), 0, 0, 0, 0, 0, t265, t266 (t122 * t36 + t125 * t37) * qJDD(3) + (t10 * t122 + t11 * t125 + (-t122 * t37 + t125 * t36) * qJD(4)) * qJD(3), -t265, -t266, t26 * t10 - t27 * t11 + t147 * t37 + t5 * t36 + t38 * t51 + t9 * t54 - g(3), 0, 0, 0, 0, 0 (-qJD(6) * t17 + t10 * t124 - t51 * t121) * t107 + t16 * t80 + t11 * t81 + t37 * t48 -(qJD(6) * t16 + t10 * t121 + t51 * t124) * t107 - t17 * t80 + t11 * t83 + t37 * t47; 0, t146 + t102, 0, t173 * t117 (-qJDD(3) * t123 - t126 * t129) * t117, 0, 0, 0, 0, 0, t278, -t277 (t122 * t75 + t125 * t76) * qJDD(3) + (t122 * t57 - t125 * t56 + (-t122 * t76 + t125 * t75) * qJD(4)) * qJD(3), -t278, t277, t26 * t57 + t27 * t56 + t147 * t76 + t5 * t75 + (-t126 * t9 + t226 * t38) * t117 + t146, 0, 0, 0, 0, 0 (qJD(6) * t168 - t121 * t203 + t124 * t57) * t107 + t167 * t80 - t56 * t81 + t76 * t48 -(qJD(6) * t167 + t121 * t57 + t124 * t203) * t107 + t168 * t80 - t56 * t83 + t76 * t47; 0, 0, qJDD(3), t283, -t102 * t233 - t164 * t213 + (-t105 * t232 - t165 * t227 + t45) * qJD(3) + t159, t114 * qJDD(3) + 0.2e1 * t122 * t198, 0.2e1 * t122 * t210 - 0.2e1 * t229 * t214, qJDD(4) * t122 + t125 * t128, qJDD(4) * t125 - t122 * t128, 0, t152 * t122 + t133 * t125, -t133 * t122 + t152 * t125, t130 + (-t45 * qJD(3) + t238) * (t114 + t115) t148 * t122 + t125 * t267, -t122 * t267 + t148 * t125 (-t122 * t26 + t125 * t27) * t45 + t254 * t38 + t130 * pkin(9) + t285 * t172, -t219 * t251 + (-t125 * t47 + t222 * t83) * t121 (-t121 * t81 + t251) * t222 + (t121 * t48 - t249 + (t121 * t83 + t124 * t81) * qJD(6)) * t125 (t107 * t218 + t47) * t122 + (t166 + t240) * t125 (t107 * t215 - t48) * t122 + (-t241 - t275) * t125, t107 * t221 + t122 * t80, -t79 * t252 + t97 * t48 - t89 * t81 - t153 * t124 + (t121 * t157 - t215 * t23 + t1) * t122 + ((-t122 * t45 + t90) * t124 + t256 * t121) * t107 + ((-t121 * t96 - t124 * t79) * t107 - t7 * t122) * qJD(6) + (-qJD(4) * t178 + t3 * t124 - t220 * t23 - t45 * t81) * t125, t97 * t47 - t89 * t83 + (-qJD(4) * t7 - t45 * t83) * t125 + (-t23 * t219 - t79 * t80 + (-qJD(6) * t96 + t256) * t107 + (-qJD(6) * t22 + t157) * t122) * t124 + (-(-qJD(6) * t79 + t90) * t107 - t3 * t125 + (qJD(4) * t23 + qJD(6) * t30 + t107 * t45 - t2) * t122 + t153) * t121; 0, 0, 0, 0, 0, -t206, t229 * t129, t211, t210, qJDD(4), -t43 * t217 - t140 + t250, qJD(4) * t255 - t43 * t225 - t144 (-pkin(4) * t122 + t242) * qJDD(3), -0.2e1 * t237 + qJDD(5) + t31 + (-qJD(3) * t84 - t65) * t125 + t140, 0.2e1 * t212 + (0.2e1 * qJD(5) - t255) * qJD(4) + (t122 * t84 + t125 * t38) * qJD(3) + t144, t147 * qJ(5) - t5 * pkin(4) - t38 * t84 - t26 * t29 - g(1) * (-pkin(4) * t14 + qJ(5) * t15) - g(2) * (-pkin(4) * t12 + qJ(5) * t13) - g(3) * (-pkin(4) * t36 + qJ(5) * t37) + t280 * t27, -t236 * t83 + t249 (-t48 - t247) * t124 + (-t47 + t248) * t121 (-t122 * t236 - t125 * t83) * qJD(3) + t275 (-t122 * t235 + t125 * t81) * qJD(3) + t166, -t107 * t225, qJ(5) * t48 + t136 * t121 + t124 * t270 + t178 * t225 + t230 * t81, qJ(5) * t47 - t121 * t270 + t136 * t124 + t7 * t225 + t230 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, qJDD(4) + t206, -t114 * t129 - t128, qJD(4) * t27 + t271 + t31 + t5, 0, 0, 0, 0, 0, -t107 * t236 - t241 + t74, -t107 ^ 2 * t124 - t240 - t252; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83 * t81, -t81 ^ 2 + t83 ^ 2, t47 + t248, t247 - t48, t80, -t121 * t8 + t1 - t23 * t83 - g(1) * (-t121 * t34 + t124 * t14) - g(2) * (t12 * t124 - t121 * t32) - g(3) * t16 + t284 * t7, -t124 * t8 - t121 * t2 + t23 * t81 - g(1) * (-t121 * t14 - t124 * t34) - g(2) * (-t12 * t121 - t124 * t32) + g(3) * t17 - t284 * t178;];
tau_reg  = t4;
