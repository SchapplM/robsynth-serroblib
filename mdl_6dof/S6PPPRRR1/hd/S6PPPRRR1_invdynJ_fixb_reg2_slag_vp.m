% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PPPRRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPPRRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:40:42
% EndTime: 2019-03-08 18:40:50
% DurationCPUTime: 4.72s
% Computational Cost: add. (9338->450), mult. (24571->673), div. (0->0), fcn. (25242->18), ass. (0->214)
t142 = sin(pkin(8));
t145 = cos(pkin(8));
t141 = sin(pkin(13));
t143 = sin(pkin(6));
t144 = cos(pkin(13));
t260 = cos(pkin(7));
t216 = t144 * t260;
t255 = sin(pkin(14));
t258 = cos(pkin(14));
t173 = (-t255 * t141 + t258 * t216) * t143;
t257 = sin(pkin(7));
t205 = t258 * t257;
t261 = cos(pkin(6));
t164 = t261 * t205 + t173;
t212 = t143 * t144 * t257;
t177 = -t261 * t260 + t212;
t276 = -t177 * t142 + t164 * t145;
t147 = sin(qJ(5));
t150 = cos(qJ(5));
t210 = pkin(5) * t147 - pkin(11) * t150;
t121 = t210 * qJD(5);
t275 = -pkin(5) * t150 - pkin(11) * t147;
t123 = -pkin(4) + t275;
t126 = t261 * qJDD(1) + qJDD(2);
t106 = -qJDD(1) * t212 + t260 * t126 + qJDD(3);
t129 = t261 * qJD(1) + qJD(2);
t107 = -qJD(1) * t212 + t260 * t129 + qJD(3);
t148 = sin(qJ(4));
t151 = cos(qJ(4));
t240 = qJD(4) * t148;
t224 = t142 * t240;
t238 = qJD(4) * t151;
t246 = t142 * t151;
t84 = qJDD(1) * t173 + t126 * t205;
t266 = t145 * t84;
t86 = qJD(1) * t173 + t129 * t205;
t265 = t145 * t86;
t172 = (t258 * t141 + t255 * t216) * t143;
t204 = t257 * t255;
t85 = qJDD(1) * t172 + t126 * t204;
t284 = -qJD(4) * t265 - t85;
t87 = qJD(1) * t172 + t129 * t204;
t27 = t106 * t246 - t107 * t224 + t284 * t148 + t151 * t266 - t87 * t238;
t12 = qJD(4) * t121 + qJDD(4) * t123 - t27;
t146 = sin(qJ(6));
t149 = cos(qJ(6));
t202 = t107 * t142 + t265;
t48 = t148 * t202 + t151 * t87;
t46 = qJD(4) * pkin(10) + t48;
t66 = t107 * t145 - t142 * t86;
t35 = t147 * t66 + t150 * t46;
t33 = qJD(5) * pkin(11) + t35;
t83 = t148 * t87;
t47 = t151 * t202 - t83;
t44 = qJD(4) * t123 - t47;
t203 = t146 * t33 - t149 * t44;
t235 = qJD(5) * t150;
t253 = qJDD(4) * pkin(10);
t223 = t142 * t238;
t247 = t142 * t148;
t211 = -t106 * t247 - t107 * t223 - t148 * t266 + t284 * t151;
t26 = -t87 * t240 - t211;
t24 = t26 + t253;
t65 = t106 * t145 - t142 * t84;
t227 = -t147 * t65 - t150 * t24 - t66 * t235;
t236 = qJD(5) * t147;
t7 = -t46 * t236 - t227;
t5 = qJDD(5) * pkin(11) + t7;
t1 = -t203 * qJD(6) + t146 * t12 + t149 * t5;
t239 = qJD(4) * t150;
t130 = -qJD(6) + t239;
t285 = -t203 * t130 + t1;
t176 = t260 * t142 + t145 * t205;
t98 = t148 * t204 - t151 * t176;
t233 = qJD(6) * t147;
t283 = qJD(4) * t233 - qJDD(5);
t228 = t147 * qJDD(4);
t92 = ((qJD(6) + t239) * qJD(5) + t228) * t146 + t283 * t149;
t282 = -qJD(6) * t150 * pkin(10) + t121 - t48;
t10 = t146 * t44 + t149 * t33;
t2 = -qJD(6) * t10 + t149 * t12 - t146 * t5;
t281 = t10 * t130 - t2;
t259 = cos(pkin(12));
t207 = t261 * t259;
t256 = sin(pkin(12));
t179 = -t256 * t141 + t144 * t207;
t169 = t179 * t260;
t178 = t141 * t207 + t256 * t144;
t199 = t143 * t205;
t157 = -t258 * t169 + t178 * t255 + t259 * t199;
t218 = t143 * t259;
t165 = t179 * t257 + t260 * t218;
t279 = t165 * t142 + t157 * t145;
t206 = t261 * t256;
t181 = -t259 * t141 - t144 * t206;
t170 = t181 * t260;
t180 = -t141 * t206 + t259 * t144;
t158 = -t258 * t170 + t180 * t255 - t256 * t199;
t217 = t143 * t256;
t166 = t181 * t257 - t260 * t217;
t278 = t166 * t142 + t158 * t145;
t277 = t276 * qJD(4);
t62 = t150 * t65;
t8 = -qJD(5) * t35 - t147 * t24 + t62;
t230 = t149 * qJD(5);
t232 = qJD(6) * t149;
t244 = t149 * t150;
t270 = -t147 * t230 * pkin(10) + t123 * t232 + t282 * t146 - t47 * t244;
t234 = qJD(6) * t146;
t245 = t146 * t150;
t269 = t146 * t236 * pkin(10) - t123 * t234 + t282 * t149 + t47 * t245;
t268 = qJD(4) * pkin(4);
t229 = qJD(4) * qJD(5);
t222 = t150 * t229;
t91 = -qJD(6) * t230 + (-t222 - t228) * t149 + t283 * t146;
t264 = t146 * t91;
t263 = t147 * t46;
t262 = t149 * t92;
t254 = qJDD(4) * pkin(4);
t241 = qJD(4) * t147;
t116 = t146 * t241 - t230;
t237 = qJD(5) * t146;
t118 = t149 * t241 + t237;
t252 = t116 * t118;
t251 = t116 * t130;
t250 = t116 * t146;
t249 = t118 * t130;
t248 = t118 * t149;
t139 = t147 ^ 2;
t140 = t150 ^ 2;
t243 = t139 - t140;
t136 = t150 * qJDD(4);
t153 = qJD(4) ^ 2;
t225 = t147 * t153 * t150;
t221 = t151 * t229;
t213 = t147 * t222;
t209 = -t10 * t146 + t149 * t203;
t159 = -t142 * t164 - t145 * t177;
t95 = t261 * t204 + t172;
t53 = t276 * t148 + t95 * t151;
t39 = t147 * t159 + t53 * t150;
t52 = t95 * t148 - t276 * t151;
t16 = t146 * t52 + t149 * t39;
t15 = -t146 * t39 + t149 * t52;
t112 = -t142 * t205 + t145 * t260;
t99 = t176 * t148 + t151 * t204;
t73 = t112 * t147 + t150 * t99;
t59 = t146 * t98 + t149 * t73;
t58 = -t146 * t73 + t149 * t98;
t34 = t150 * t66 - t263;
t201 = t150 * t112 - t147 * t99;
t200 = qJDD(4) * t151 - t148 * t153;
t198 = t143 * t204;
t114 = t145 * t147 + t150 * t247;
t104 = -t114 * t146 - t149 * t246;
t194 = -t114 * t149 + t146 * t246;
t113 = -t150 * t145 + t147 * t247;
t115 = t147 * t229 + qJDD(6) - t136;
t193 = t115 * t146 - t130 * t232;
t192 = t115 * t149 + t130 * t234;
t38 = t53 * t147 - t150 * t159;
t70 = t255 * t169 + t178 * t258 - t259 * t198;
t41 = -t279 * t148 + t70 * t151;
t71 = t255 * t170 + t180 * t258 + t256 * t198;
t43 = -t278 * t148 + t71 * t151;
t56 = t142 * t157 - t145 * t165;
t57 = t142 * t158 - t145 * t166;
t191 = -g(1) * (-t147 * t43 + t150 * t57) - g(2) * (-t147 * t41 + t150 * t56) + g(3) * t38;
t18 = t147 * t56 + t150 * t41;
t20 = t147 * t57 + t150 * t43;
t190 = -g(1) * t20 - g(2) * t18 - g(3) * t39;
t40 = t148 * t70 + t279 * t151;
t42 = t148 * t71 + t278 * t151;
t189 = g(1) * t42 + g(2) * t40 + g(3) * t52;
t188 = g(1) * t43 + g(2) * t41 + g(3) * t53;
t6 = -qJDD(5) * pkin(5) - t8;
t187 = t191 - t6;
t32 = -qJD(5) * pkin(5) - t34;
t185 = -pkin(11) * t115 - t130 * t32;
t45 = -t47 - t268;
t184 = -pkin(10) * qJDD(5) + (t45 + t47 - t268) * qJD(5);
t183 = qJD(4) * t48 + t189;
t182 = -g(1) * t217 + g(2) * t218 - g(3) * t261;
t175 = pkin(11) * qJD(6) * t130 + t187;
t152 = qJD(5) ^ 2;
t25 = -t27 - t254;
t168 = -pkin(10) * t152 + t183 - t25 + t254;
t167 = -t8 * t147 + t7 * t150 + (-t147 * t35 - t150 * t34) * qJD(5) - t188;
t154 = g(1) * t166 + g(2) * t165 + g(3) * t177;
t120 = t210 * qJD(4);
t110 = pkin(10) * t244 + t123 * t146;
t109 = -pkin(10) * t245 + t123 * t149;
t103 = qJD(5) * t114 + t147 * t223;
t102 = -qJD(5) * t113 + t150 * t223;
t94 = t99 * qJD(4);
t93 = t98 * qJD(4);
t64 = qJD(6) * t194 - t102 * t146 + t149 * t224;
t63 = qJD(6) * t104 + t102 * t149 + t146 * t224;
t55 = qJD(5) * t73 - t147 * t93;
t54 = qJD(5) * t201 - t150 * t93;
t51 = t52 * pkin(4);
t50 = t277 * t148 + t95 * t238;
t49 = t277 * t151 - t95 * t240;
t37 = t42 * pkin(4);
t36 = t40 * pkin(4);
t31 = t120 * t146 + t149 * t34;
t30 = t120 * t149 - t146 * t34;
t29 = -t59 * qJD(6) - t146 * t54 + t149 * t94;
t28 = t58 * qJD(6) + t146 * t94 + t149 * t54;
t14 = -qJD(5) * t38 + t49 * t150;
t13 = qJD(5) * t39 + t49 * t147;
t4 = t15 * qJD(6) + t14 * t149 + t146 * t50;
t3 = -t16 * qJD(6) - t14 * t146 + t149 * t50;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t126 * t261 - g(3) + (t141 ^ 2 + t144 ^ 2) * t143 ^ 2 * qJDD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106 * t177 + t164 * t84 + t85 * t95 - g(3), 0, 0, 0, 0, 0, 0, -qJD(4) * t50 - qJDD(4) * t52, -qJD(4) * t49 - qJDD(4) * t53, 0, t159 * t65 + t26 * t53 - t27 * t52 - t47 * t50 + t48 * t49 - g(3), 0, 0, 0, 0, 0, 0, -t52 * t136 - qJD(5) * t13 - qJDD(5) * t38 + (-t150 * t50 + t52 * t236) * qJD(4), t52 * t228 - qJD(5) * t14 - qJDD(5) * t39 + (t147 * t50 + t52 * t235) * qJD(4) (t147 * t38 + t150 * t39) * qJDD(4) + (t13 * t147 + t14 * t150 + (-t147 * t39 + t150 * t38) * qJD(5)) * qJD(4), -t13 * t34 + t14 * t35 + t25 * t52 - t38 * t8 + t39 * t7 + t45 * t50 - g(3), 0, 0, 0, 0, 0, 0, t115 * t15 + t116 * t13 - t130 * t3 + t38 * t92, -t115 * t16 + t118 * t13 + t130 * t4 - t38 * t91, -t116 * t4 - t118 * t3 + t15 * t91 - t16 * t92, t1 * t16 + t10 * t4 + t13 * t32 + t15 * t2 - t203 * t3 + t38 * t6 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182 + t126, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106 * t260 + t85 * t204 + t84 * t205 + t182, 0, 0, 0, 0, 0, 0, -qJD(4) * t94 - qJDD(4) * t98, qJD(4) * t93 - qJDD(4) * t99, 0, t65 * t112 + t26 * t99 - t27 * t98 - t47 * t94 - t48 * t93 + t182, 0, 0, 0, 0, 0, 0, -t98 * t136 - qJD(5) * t55 + qJDD(5) * t201 + (-t150 * t94 + t98 * t236) * qJD(4), t98 * t228 - qJD(5) * t54 - qJDD(5) * t73 + (t147 * t94 + t98 * t235) * qJD(4) (-t147 * t201 + t150 * t73) * qJDD(4) + (t147 * t55 + t150 * t54 + (-t147 * t73 - t150 * t201) * qJD(5)) * qJD(4), t201 * t8 + t25 * t98 - t34 * t55 + t35 * t54 + t45 * t94 + t7 * t73 + t182, 0, 0, 0, 0, 0, 0, t115 * t58 + t116 * t55 - t130 * t29 - t201 * t92, -t115 * t59 + t118 * t55 + t130 * t28 + t201 * t91, -t116 * t28 - t118 * t29 + t58 * t91 - t59 * t92, t1 * t59 + t10 * t28 + t2 * t58 - t201 * t6 - t203 * t29 + t32 * t55 + t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154 + t106, 0, 0, 0, 0, 0, 0, t200 * t142 (-qJDD(4) * t148 - t151 * t153) * t142, 0, t65 * t145 + (t26 * t148 + t27 * t151 + (-t148 * t47 + t151 * t48) * qJD(4)) * t142 + t154, 0, 0, 0, 0, 0, 0, -qJD(5) * t103 - qJDD(5) * t113 + (-t147 * t221 + t150 * t200) * t142, -qJD(5) * t102 - qJDD(5) * t114 + (-t147 * t200 - t150 * t221) * t142 (t113 * t147 + t114 * t150) * qJDD(4) + (t102 * t150 + t103 * t147 + (t113 * t150 - t114 * t147) * qJD(5)) * qJD(4), t7 * t114 + t35 * t102 - t8 * t113 - t34 * t103 + (-t151 * t25 + t240 * t45) * t142 + t154, 0, 0, 0, 0, 0, 0, t103 * t116 + t104 * t115 + t113 * t92 - t130 * t64, t103 * t118 - t113 * t91 + t115 * t194 + t130 * t63, t104 * t91 - t116 * t63 - t118 * t64 + t194 * t92, -t1 * t194 + t10 * t63 + t32 * t103 + t2 * t104 + t6 * t113 - t203 * t64 + t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(4), t183 + t27 (t47 + t83) * qJD(4) + t188 + t211, 0, 0, qJDD(4) * t139 + 0.2e1 * t213, 0.2e1 * t147 * t136 - 0.2e1 * t243 * t229, qJDD(5) * t147 + t150 * t152, qJDD(4) * t140 - 0.2e1 * t213, qJDD(5) * t150 - t147 * t152, 0, t147 * t184 + t150 * t168, -t147 * t168 + t150 * t184, t167 + (-t47 * qJD(4) + t253) * (t139 + t140) -t25 * pkin(4) + g(1) * t37 + g(2) * t36 + g(3) * t51 - t45 * t48 + (t34 * t147 - t35 * t150) * t47 + t167 * pkin(10), -t147 * t149 * t91 + (-t146 * t233 + t150 * t230) * t118 (-t116 * t149 - t118 * t146) * t235 + (t264 - t262 + (-t248 + t250) * qJD(6)) * t147 (-t130 * t230 + t91) * t150 + (qJD(5) * t118 + t192) * t147, t146 * t147 * t92 + (t146 * t235 + t147 * t232) * t116 (t130 * t237 + t92) * t150 + (-qJD(5) * t116 - t193) * t147, -t115 * t150 - t130 * t236, t109 * t115 - t269 * t130 - t188 * t146 + (-t2 + (pkin(10) * t116 + t146 * t32) * qJD(5) + t189 * t149) * t150 + (pkin(10) * t92 - qJD(5) * t203 - t47 * t116 + t6 * t146 + t232 * t32) * t147, -t110 * t115 + t270 * t130 - t188 * t149 + (t1 + (pkin(10) * t118 + t149 * t32) * qJD(5) - t189 * t146) * t150 + (-pkin(10) * t91 - t10 * qJD(5) - t47 * t118 + t6 * t149 - t234 * t32) * t147, t109 * t91 - t110 * t92 - t269 * t118 - t270 * t116 + t209 * t235 + (-t1 * t146 - t149 * t2 + (-t10 * t149 - t146 * t203) * qJD(6) + t189) * t147, t1 * t110 + t2 * t109 - t32 * t147 * t47 - g(1) * (t275 * t42 - t37) - g(2) * (t275 * t40 - t36) - g(3) * (t275 * t52 - t51) - t269 * t203 + t270 * t10 + (t6 * t147 + t235 * t32 - t188) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t225, t243 * t153, t228, t225, t136, qJDD(5), t62 + (-qJD(4) * t45 - t24) * t147 + t191, -t45 * t239 + (t34 + t263) * qJD(5) - t190 + t227, 0, 0, -t130 * t248 - t264 (-t91 + t251) * t149 + (-t92 + t249) * t146 (-t118 * t147 + t130 * t244) * qJD(4) + t193, -t130 * t250 - t262 (t116 * t147 - t130 * t245) * qJD(4) + t192, t130 * t241, -pkin(5) * t92 - t116 * t35 + t130 * t30 + t146 * t185 + t149 * t175 + t203 * t241, pkin(5) * t91 + t10 * t241 - t118 * t35 - t130 * t31 - t146 * t175 + t149 * t185, t116 * t31 + t118 * t30 + ((qJD(6) * t118 - t92) * pkin(11) + t285) * t149 + ((qJD(6) * t116 - t91) * pkin(11) + t281) * t146 + t190, -t10 * t31 + t203 * t30 - t32 * t35 + t187 * pkin(5) + (qJD(6) * t209 + t1 * t149 - t2 * t146 + t190) * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t252, -t116 ^ 2 + t118 ^ 2, -t91 - t251, -t252, -t249 - t92, t115, -t32 * t118 - g(1) * (-t146 * t20 + t149 * t42) - g(2) * (-t146 * t18 + t149 * t40) - g(3) * t15 - t281, t32 * t116 - g(1) * (-t146 * t42 - t149 * t20) - g(2) * (-t146 * t40 - t149 * t18) + g(3) * t16 - t285, 0, 0;];
tau_reg  = t9;
