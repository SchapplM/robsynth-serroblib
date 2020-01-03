% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRP10
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP10_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP10_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:11
% EndTime: 2019-12-31 20:11:18
% DurationCPUTime: 3.65s
% Computational Cost: add. (2992->432), mult. (6333->514), div. (0->0), fcn. (3682->6), ass. (0->240)
t151 = sin(qJ(2));
t222 = t151 * qJD(1);
t119 = qJD(4) + t222;
t153 = cos(qJ(4));
t150 = sin(qJ(4));
t223 = t150 * qJD(2);
t154 = cos(qJ(2));
t229 = qJD(1) * t154;
t87 = t153 * t229 + t223;
t182 = t87 * t119;
t218 = qJD(1) * qJD(2);
t204 = t151 * t218;
t216 = t154 * qJDD(1);
t295 = -t204 + t216;
t41 = qJD(4) * t87 - t153 * qJDD(2) + t295 * t150;
t297 = t41 - t182;
t296 = t41 + t182;
t283 = pkin(3) + pkin(6);
t259 = t41 * qJ(5);
t203 = t154 * t218;
t217 = t151 * qJDD(1);
t172 = t203 + t217;
t86 = qJDD(4) + t172;
t279 = t86 * pkin(4);
t118 = pkin(2) * t204;
t253 = qJ(3) * t154;
t186 = pkin(7) * t151 - t253;
t221 = t151 * qJD(3);
t167 = qJD(2) * t186 - t221;
t136 = t151 * qJ(3);
t201 = -pkin(1) - t136;
t284 = pkin(2) + pkin(7);
t171 = -t284 * t154 + t201;
t30 = qJD(1) * t167 + qJDD(1) * t171 + t118;
t58 = t171 * qJD(1);
t130 = pkin(6) * t222;
t290 = qJD(3) + t130;
t234 = pkin(3) * t222 + t290;
t62 = -t284 * qJD(2) + t234;
t32 = t150 * t62 + t153 * t58;
t117 = pkin(6) * t203;
t127 = pkin(6) * t217;
t202 = qJDD(3) + t117 + t127;
t44 = t172 * pkin(3) - t284 * qJDD(2) + t202;
t7 = -qJD(4) * t32 - t150 * t30 + t153 * t44;
t205 = t150 * t229;
t220 = t153 * qJD(2);
t89 = -t205 + t220;
t1 = -t89 * qJD(5) + t259 + t279 + t7;
t31 = -t150 * t58 + t153 * t62;
t20 = -t89 * qJ(5) + t31;
t19 = t119 * pkin(4) + t20;
t42 = -qJD(4) * t205 + t150 * qJDD(2) + (qJD(2) * qJD(4) + t295) * t153;
t258 = t42 * qJ(5);
t225 = qJD(4) * t153;
t226 = qJD(4) * t150;
t6 = t150 * t44 + t153 * t30 + t62 * t225 - t58 * t226;
t2 = -t87 * qJD(5) - t258 + t6;
t155 = cos(qJ(1));
t245 = t151 * t155;
t152 = sin(qJ(1));
t247 = t151 * t152;
t210 = -g(1) * t245 - g(2) * t247 + g(3) * t154;
t21 = -t87 * qJ(5) + t32;
t262 = t21 * t119;
t294 = -(t119 * t19 - t2) * t150 + (t1 + t262) * t153 + t210;
t256 = t89 * t119;
t27 = -t42 + t256;
t292 = t42 + t256;
t69 = t153 * t86;
t291 = -t119 * t226 + t69;
t187 = g(1) * t155 + g(2) * t152;
t289 = t42 * pkin(4) + qJDD(5);
t240 = t153 * t154;
t236 = t155 * t153;
t243 = t152 * t150;
t73 = t151 * t236 - t243;
t237 = t155 * t150;
t242 = t152 * t153;
t75 = t151 * t242 + t237;
t288 = -g(1) * t73 - g(2) * t75 + g(3) * t240;
t146 = qJD(2) * qJ(3);
t131 = pkin(6) * t229;
t96 = pkin(3) * t229 + t131;
t71 = t146 + t96;
t286 = t119 * t71 - t284 * t86;
t285 = t89 ^ 2;
t278 = g(1) * t152;
t275 = g(2) * t155;
t274 = g(3) * t151;
t273 = t150 * pkin(4);
t140 = t154 * pkin(2);
t272 = t154 * pkin(7);
t271 = t89 * t87;
t269 = -t20 + t19;
t134 = pkin(2) * t222;
t66 = qJD(1) * t186 + t134;
t40 = t150 * t96 + t153 * t66;
t219 = t153 * qJD(5);
t235 = qJ(5) + t284;
t249 = t150 * t151;
t39 = -t150 * t66 + t153 * t96;
t268 = t235 * t226 - t219 - (pkin(4) * t154 - qJ(5) * t249) * qJD(1) - t39;
t99 = t235 * t153;
t267 = -t153 * qJ(5) * t222 - qJD(4) * t99 - t150 * qJD(5) - t40;
t106 = t283 * t151;
t233 = t140 + t136;
t190 = t233 + t272;
t81 = -pkin(1) - t190;
t48 = t150 * t106 + t153 * t81;
t266 = t150 * t86;
t265 = t153 * t41;
t264 = t153 * t89;
t261 = t31 * t119;
t260 = t32 * t119;
t257 = t42 * t150;
t126 = t153 * pkin(4) + pkin(3);
t255 = pkin(4) * t225 + t126 * t222 + t290;
t254 = pkin(6) * qJDD(2);
t252 = qJD(2) * t87;
t251 = qJD(2) * t89;
t250 = qJDD(2) * pkin(2);
t248 = t150 * t154;
t246 = t151 * t153;
t158 = qJD(1) ^ 2;
t244 = t151 * t158;
t241 = t152 * t154;
t149 = -qJ(5) - pkin(7);
t239 = t154 * t149;
t238 = t154 * t155;
t107 = t283 * t154;
t232 = t155 * pkin(1) + t152 * pkin(6);
t147 = t151 ^ 2;
t148 = t154 ^ 2;
t231 = t147 - t148;
t230 = t147 + t148;
t228 = qJD(2) * t151;
t227 = qJD(2) * t154;
t224 = qJD(4) * t154;
t215 = pkin(4) * t248;
t214 = t151 * t237;
t200 = -t87 * pkin(4) - qJD(5);
t46 = -t200 + t71;
t213 = t46 * t226;
t212 = t46 * t225;
t128 = pkin(6) * t216;
t144 = qJDD(2) * qJ(3);
t145 = qJD(2) * qJD(3);
t211 = t128 + t144 + t145;
t209 = t283 * qJD(2);
t207 = t150 * t224;
t206 = t119 * t229;
t199 = qJ(5) * t154 - t81;
t198 = -qJD(2) * pkin(2) + qJD(3);
t197 = pkin(2) * t238 + qJ(3) * t245 + t232;
t196 = -t127 - t210;
t195 = pkin(3) * t216 + t211;
t194 = qJD(1) * t209;
t193 = t151 * t203;
t192 = g(1) * t75 - g(2) * t73;
t74 = t214 + t242;
t76 = -t151 * t243 + t236;
t191 = -g(1) * t76 - g(2) * t74;
t189 = t230 * qJDD(1) * pkin(6);
t157 = qJD(2) ^ 2;
t188 = pkin(6) * t157 + t275;
t185 = -t150 * t31 + t153 * t32;
t100 = t130 + t198;
t105 = -t131 - t146;
t183 = t100 * t154 + t105 * t151;
t181 = t119 * t150;
t180 = t201 - t140;
t177 = t187 * t154;
t72 = t180 * qJD(1);
t176 = t72 * t222 + qJDD(3) - t196;
t175 = -0.2e1 * pkin(1) * t218 - t254;
t174 = -t119 * t225 - t266;
t133 = pkin(2) * t228;
t55 = t133 + t167;
t97 = t283 * t227;
t16 = t106 * t225 + t150 * t97 + t153 * t55 - t81 * t226;
t173 = -qJ(3) * t227 - t221;
t170 = 0.2e1 * qJDD(1) * pkin(1) - t188;
t101 = -pkin(1) - t233;
t169 = t254 + (-qJD(1) * t101 - t72) * qJD(2);
t168 = -t177 - t274;
t45 = -t151 * t194 + t195;
t18 = t45 + t289;
t166 = t18 + t168;
t43 = qJD(1) * t173 + qJDD(1) * t180 + t118;
t68 = t133 + t173;
t165 = qJD(1) * t68 + qJDD(1) * t101 + t188 + t43;
t164 = g(1) * t74 - g(2) * t76 - g(3) * t248 - t6;
t59 = pkin(6) * t204 - t211;
t65 = t202 - t250;
t163 = qJD(2) * t183 + t65 * t151 - t59 * t154;
t162 = qJD(4) * t119 * t284 + t168 + t45;
t160 = t7 + t288;
t141 = t155 * pkin(6);
t125 = qJ(3) + t273;
t123 = g(1) * t241;
t116 = qJ(3) * t238;
t114 = qJ(3) * t241;
t113 = t154 * t244;
t104 = t231 * t158;
t103 = qJDD(2) * t154 - t157 * t151;
t102 = qJDD(2) * t151 + t157 * t154;
t98 = t235 * t150;
t95 = t151 * t209;
t93 = -qJ(3) * t229 + t134;
t92 = t153 * t106;
t85 = t87 ^ 2;
t83 = t148 * qJDD(1) - 0.2e1 * t193;
t82 = t147 * qJDD(1) + 0.2e1 * t193;
t80 = t153 * t97;
t70 = pkin(4) * t240 + t107;
t54 = 0.2e1 * t151 * t216 - 0.2e1 * t231 * t218;
t50 = -pkin(4) * t207 + (-pkin(6) - t126) * t228;
t49 = t119 * t227 + t86 * t151;
t47 = -t150 * t81 + t92;
t36 = -qJ(5) * t240 + t48;
t35 = -t85 + t285;
t34 = t151 * pkin(4) + t199 * t150 + t92;
t25 = -t119 ^ 2 * t153 - t251 - t266;
t24 = -t119 * t181 - t252 + t69;
t23 = (-t119 * t246 + t154 * t87) * qJD(1) + t174;
t22 = (-t119 * t249 - t154 * t89) * qJD(1) + t291;
t17 = -t48 * qJD(4) - t150 * t55 + t80;
t15 = t153 * t182 + t257;
t14 = -t181 * t89 - t265;
t13 = -t87 * t207 + (t154 * t42 - t87 * t228) * t153;
t12 = -t224 * t264 + (t154 * t41 + t89 * t228) * t150;
t11 = -t154 * t219 + (t151 * t220 + t207) * qJ(5) + t16;
t10 = pkin(4) * t227 + t80 + t199 * t225 + (-qJ(5) * t228 - qJD(4) * t106 + qJD(5) * t154 - t55) * t150;
t9 = (t119 * t223 - t41) * t151 + (t174 + t251) * t154;
t8 = (t119 * t220 - t42) * t151 + (-t252 - t291) * t154;
t5 = t27 * t150 + t153 * t297;
t4 = t150 * t296 - t292 * t153;
t3 = (-t150 * t87 + t264) * t228 + (t257 + t265 + (t150 * t89 + t153 * t87) * qJD(4)) * t154;
t26 = [0, 0, 0, 0, 0, qJDD(1), -t275 + t278, t187, 0, 0, t82, t54, t102, t83, t103, 0, t151 * t175 + t154 * t170 + t123, t175 * t154 + (-t170 - t278) * t151, -t187 + 0.2e1 * t189, -g(1) * (-t152 * pkin(1) + t141) - g(2) * t232 + (t230 * pkin(6) ^ 2 + pkin(1) ^ 2) * qJDD(1), 0, -t102, -t103, t82, t54, t83, t189 + t163 - t187, t151 * t169 + t154 * t165 - t123, t169 * t154 + (-t165 + t278) * t151, t163 * pkin(6) - g(1) * t141 - g(2) * t197 + t43 * t101 - t180 * t278 + t72 * t68, t12, t3, t9, t13, t8, t49, t107 * t42 + t17 * t119 + t47 * t86 - t95 * t87 + (-t71 * t220 + t7) * t151 + (qJD(2) * t31 + t45 * t153 - t71 * t226) * t154 + t191, -t107 * t41 - t16 * t119 - t48 * t86 - t95 * t89 + (t71 * t223 - t6) * t151 + (-qJD(2) * t32 - t45 * t150 - t71 * t225) * t154 + t192, -t16 * t87 - t17 * t89 + t47 * t41 - t48 * t42 + t123 + t185 * t228 + (-t275 + t150 * t7 - t153 * t6 + (t150 * t32 + t153 * t31) * qJD(4)) * t154, t6 * t48 + t32 * t16 + t7 * t47 + t31 * t17 + t45 * t107 - t71 * t95 - g(1) * (t155 * pkin(3) + t141) - g(2) * (pkin(7) * t238 + t197) + (-g(1) * (t180 - t272) - g(2) * pkin(3)) * t152, t12, t3, t9, t13, t8, t49, t10 * t119 + t34 * t86 + t70 * t42 + t50 * t87 + (-t46 * t220 + t1) * t151 + (qJD(2) * t19 + t18 * t153 - t213) * t154 + t191, -t11 * t119 - t36 * t86 - t70 * t41 + t50 * t89 + (t46 * t223 - t2) * t151 + (-qJD(2) * t21 - t18 * t150 - t212) * t154 + t192, -t10 * t89 - t11 * t87 + t34 * t41 - t36 * t42 + t123 + (-t150 * t19 + t153 * t21) * t228 + (-t275 + t1 * t150 - t153 * t2 + (t150 * t21 + t153 * t19) * qJD(4)) * t154, t2 * t36 + t21 * t11 + t1 * t34 + t19 * t10 + t18 * t70 + t46 * t50 - g(1) * (t155 * t126 + t141) - g(2) * (pkin(4) * t214 - t149 * t238 + t197) + (-g(1) * (-pkin(4) * t249 + t180 + t239) - g(2) * t126) * t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, t104, t217, t113, t216, qJDD(2), pkin(1) * t244 + t196, t274 - t128 + (pkin(1) * t158 + t187) * t154, 0, 0, qJDD(2), -t217, -t216, -t113, t104, t113, (-pkin(2) * t151 + t253) * qJDD(1) + ((-t105 - t146) * t151 + (-t100 + t198) * t154) * qJD(1), -t93 * t229 + t176 - 0.2e1 * t250, t128 + 0.2e1 * t144 + 0.2e1 * t145 + (qJD(1) * t93 - g(3)) * t151 + (qJD(1) * t72 - t187) * t154, -t59 * qJ(3) - t105 * qJD(3) - t65 * pkin(2) - t72 * t93 - g(1) * (-pkin(2) * t245 + t116) - g(2) * (-pkin(2) * t247 + t114) - g(3) * t233 - t183 * qJD(1) * pkin(6), t14, t4, t22, t15, t23, -t206, qJ(3) * t42 - t39 * t119 + t162 * t150 + t286 * t153 - t31 * t229 + t234 * t87, -qJ(3) * t41 + t40 * t119 - t286 * t150 + t162 * t153 + t32 * t229 + t234 * t89, t39 * t89 + t40 * t87 + (-t32 * t222 - t284 * t41 - t7 + (t284 * t87 - t32) * qJD(4)) * t153 + (t31 * t222 + t284 * t42 - t6 + (-t284 * t89 + t31) * qJD(4)) * t150 - t210, -g(1) * t116 - g(2) * t114 - g(3) * t190 + t45 * qJ(3) + t234 * t71 - t31 * t39 - t32 * t40 + (-qJD(4) * t185 - t6 * t150 + t187 * t151 - t7 * t153) * t284, t14, t4, t22, t15, t23, -t206, t212 + t125 * t42 - t99 * t86 + t255 * t87 + t268 * t119 + (-t154 * t19 + t46 * t246) * qJD(1) + t166 * t150, -t213 - t125 * t41 + t98 * t86 + t255 * t89 - t267 * t119 + (t154 * t21 - t46 * t249) * qJD(1) + t166 * t153, -t267 * t87 - t268 * t89 - t99 * t41 + t98 * t42 - t294, -t2 * t98 - t1 * t99 + t18 * t125 - g(1) * (t155 * t215 + t116) - g(2) * (t152 * t215 + t114) - g(3) * (t233 - t239) + t255 * t46 + t267 * t21 + t268 * t19 + (-g(3) * t273 + t187 * (pkin(2) - t149)) * t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t217, qJDD(2) + t113, -t147 * t158 - t157, t105 * qJD(2) + t117 + t176 - t250, 0, 0, 0, 0, 0, 0, t24, t25, t5, -t71 * qJD(2) + (t7 + t260) * t153 + (t6 - t261) * t150 + t210, 0, 0, 0, 0, 0, 0, t24, t25, t5, -t46 * qJD(2) + t294; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t271, t35, -t297, -t271, t27, t86, -t71 * t89 + t160 + t260, t71 * t87 + t164 + t261, 0, 0, t271, t35, -t297, -t271, t27, t86, 0.2e1 * t279 + t259 + t262 + (t200 - t46) * t89 + t160, -t285 * pkin(4) + t258 + t20 * t119 + (qJD(5) + t46) * t87 + t164, t41 * pkin(4) - t269 * t87, t269 * t21 + (-t46 * t89 + t1 + t288) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t292, -t296, -t85 - t285, t19 * t89 + t21 * t87 - t177 + (-g(3) - t194) * t151 + t195 + t289;];
tau_reg = t26;
