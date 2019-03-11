% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:03:27
% EndTime: 2019-03-08 23:03:34
% DurationCPUTime: 2.90s
% Computational Cost: add. (5139->316), mult. (13424->460), div. (0->0), fcn. (10639->12), ass. (0->205)
t159 = sin(qJ(4));
t160 = sin(qJ(3));
t163 = cos(qJ(4));
t164 = cos(qJ(3));
t121 = t159 * t160 - t163 * t164;
t274 = pkin(8) + pkin(9);
t225 = qJD(3) * t274;
t123 = t160 * t225;
t124 = t164 * t225;
t129 = t274 * t164;
t165 = cos(qJ(2));
t156 = sin(pkin(6));
t240 = qJD(1) * t156;
t220 = t165 * t240;
t234 = qJD(4) * t159;
t128 = t274 * t160;
t249 = t128 * t163;
t290 = qJD(4) * t249 - t121 * t220 + t163 * t123 + t159 * t124 + t129 * t234;
t122 = t159 * t164 + t160 * t163;
t186 = t128 * t159 - t129 * t163;
t289 = t186 * qJD(4) + t122 * t220 + t159 * t123 - t163 * t124;
t237 = qJD(2) * t160;
t221 = t159 * t237;
t235 = qJD(2) * t164;
t114 = -t163 * t235 + t221;
t218 = t163 * t237;
t116 = -t159 * t235 - t218;
t155 = sin(pkin(12));
t253 = cos(pkin(12));
t202 = -t253 * t114 + t116 * t155;
t276 = qJD(6) - t202;
t288 = qJD(6) - t276;
t152 = qJD(3) + qJD(4);
t93 = t152 * t122;
t287 = -qJ(5) * t93 - qJD(5) * t121 - t290;
t158 = sin(qJ(6));
t162 = cos(qJ(6));
t178 = -t155 * t114 - t253 * t116;
t64 = -t162 * t152 + t158 * t178;
t286 = t276 * t64;
t92 = t152 * t121;
t285 = t92 * qJ(5) - t122 * qJD(5) + t289;
t231 = qJD(2) * qJD(3);
t213 = t164 * t231;
t284 = qJD(4) * t235 + t213;
t161 = sin(qJ(2));
t224 = t161 * t240;
t261 = qJD(3) * pkin(3);
t176 = t160 * t261 - t224;
t204 = t162 * t276;
t214 = t160 * t231;
t200 = qJD(4) * t218 + t284 * t159 + t163 * t214;
t80 = -t152 * t221 + t284 * t163;
t55 = t155 * t80 + t253 * t200;
t259 = t158 * t55;
t283 = -t276 * t204 - t259;
t148 = pkin(3) * t163 + pkin(4);
t207 = t253 * t159;
t107 = pkin(3) * t207 + t155 * t148;
t103 = pkin(10) + t107;
t150 = pkin(3) * t237;
t273 = pkin(4) * t116;
t38 = pkin(5) * t178 - pkin(10) * t202 - t273;
t282 = (qJD(6) * t103 + t150 + t38) * t276;
t233 = qJD(6) * t158;
t280 = -t162 * t55 + t233 * t276;
t238 = qJD(2) * t156;
t215 = qJD(1) * t238;
t197 = t165 * t215;
t206 = t274 * qJD(2) + t224;
t157 = cos(pkin(6));
t239 = qJD(1) * t157;
t96 = -t206 * t160 + t164 * t239;
t68 = t96 * qJD(3) + t164 * t197;
t90 = t96 + t261;
t279 = (qJD(4) * t90 + t68) * t163;
t108 = t116 * qJ(5);
t97 = t160 * t239 + t206 * t164;
t87 = t159 * t97;
t209 = t163 * t90 - t87;
t49 = t108 + t209;
t277 = pkin(4) * t93 + t176;
t69 = -t97 * qJD(3) - t160 * t197;
t210 = t159 * t69 - t97 * t234;
t11 = -t200 * qJ(5) - t114 * qJD(5) + t210 + t279;
t89 = t163 * t97;
t189 = -t159 * t90 - t89;
t190 = -t159 * t68 + t163 * t69;
t171 = t189 * qJD(4) + t190;
t168 = -t80 * qJ(5) + t116 * qJD(5) + t171;
t3 = t155 * t11 - t253 * t168;
t175 = -qJ(5) * t122 - t129 * t159 - t249;
t73 = -qJ(5) * t121 - t186;
t45 = t155 * t175 + t253 * t73;
t86 = -t155 * t121 + t253 * t122;
t196 = t3 * t86 - t45 * t55;
t252 = qJ(5) * t114;
t50 = -t189 - t252;
t260 = t155 * t50;
t40 = pkin(4) * t152 + t49;
t22 = t253 * t40 - t260;
t20 = -t152 * pkin(5) - t22;
t265 = t285 * t155 + t287 * t253;
t149 = -pkin(3) * t164 - pkin(2);
t109 = t149 * qJD(2) - t220;
t81 = t114 * pkin(4) + qJD(5) + t109;
t31 = -pkin(5) * t202 - pkin(10) * t178 + t81;
t4 = t253 * t11 + t155 * t168;
t185 = pkin(4) * t121 + t149;
t85 = t253 * t121 + t122 * t155;
t41 = pkin(5) * t85 - pkin(10) * t86 + t185;
t58 = -t155 * t93 - t253 * t92;
t275 = -(qJD(6) * t31 + t4) * t85 + t20 * t58 + (-qJD(6) * t41 - t265) * t276 + t196;
t272 = t20 * t202;
t271 = t20 * t86;
t270 = t41 * t55;
t269 = t64 * t178;
t66 = t152 * t158 + t162 * t178;
t268 = t66 * t178;
t267 = t276 * t178;
t266 = t287 * t155 - t285 * t253;
t42 = t253 * t50;
t23 = t155 * t40 + t42;
t264 = t163 * t96 - t87;
t263 = pkin(3) * qJD(4);
t262 = qJD(2) * pkin(2);
t56 = -t155 * t200 + t253 * t80;
t258 = t158 * t56;
t257 = t158 * t202;
t232 = qJD(6) * t162;
t29 = t152 * t232 + t162 * t56 - t178 * t233;
t256 = t29 * t158;
t208 = -t159 * t96 - t89;
t183 = t208 + t252;
t54 = t108 + t264;
t255 = -t155 * t54 + t253 * t183 + (t155 * t163 + t207) * t263;
t248 = t155 * t159;
t254 = -t155 * t183 - t253 * t54 + (t253 * t163 - t248) * t263;
t251 = t109 * t116;
t250 = t116 * t114;
t247 = t156 * t161;
t246 = t156 * t165;
t167 = qJD(2) ^ 2;
t245 = t156 * t167;
t166 = qJD(3) ^ 2;
t244 = t166 * t160;
t243 = t166 * t164;
t110 = pkin(3) * t214 + t161 * t215;
t241 = t160 ^ 2 - t164 ^ 2;
t236 = qJD(2) * t161;
t229 = t86 * t233;
t228 = t276 * t232;
t227 = t161 * t245;
t21 = pkin(10) * t152 + t23;
t191 = t158 * t21 - t162 * t31;
t226 = t178 * t191 + t20 * t233;
t223 = t156 * t236;
t222 = t165 * t238;
t217 = -pkin(3) * t152 - t90;
t6 = t158 * t31 + t162 * t21;
t201 = t3 * t158 + t178 * t6 + t20 * t232;
t199 = t160 * t222;
t198 = t164 * t222;
t57 = -t155 * t92 + t253 * t93;
t195 = pkin(5) * t57 - pkin(10) * t58 + t277;
t194 = t178 * t23 + t202 * t22;
t193 = t276 * t58 + t55 * t86;
t112 = t157 * t164 - t160 * t247;
t113 = t157 * t160 + t164 * t247;
t187 = t112 * t163 - t113 * t159;
t75 = t112 * t159 + t113 * t163;
t184 = t257 * t276 - t280;
t182 = t109 * t114 - t210;
t47 = t155 * t187 + t253 * t75;
t181 = -t158 * t47 - t162 * t246;
t180 = t158 * t246 - t162 * t47;
t179 = t262 * qJD(2);
t106 = -pkin(3) * t248 + t253 * t148;
t63 = t200 * pkin(4) + t110;
t174 = -t103 * t55 - t254 * t276 - t272;
t173 = -0.2e1 * qJD(3) * t262;
t146 = -t253 * pkin(4) - pkin(5);
t145 = pkin(4) * t155 + pkin(10);
t102 = -pkin(5) - t106;
t95 = -qJD(3) * t113 - t199;
t94 = qJD(3) * t112 + t198;
t67 = -t114 ^ 2 + t116 ^ 2;
t62 = -t116 * t152 - t200;
t61 = t114 * t152 + t80;
t46 = t155 * t75 - t253 * t187;
t44 = t155 * t73 - t253 * t175;
t35 = -t75 * qJD(4) - t159 * t94 + t163 * t95;
t34 = t187 * qJD(4) + t159 * t95 + t163 * t94;
t30 = t66 * qJD(6) + t258;
t25 = t253 * t49 - t260;
t24 = t155 * t49 + t42;
t17 = t55 * pkin(5) - t56 * pkin(10) + t63;
t14 = t162 * t17;
t13 = t155 * t35 + t253 * t34;
t12 = t155 * t34 - t253 * t35;
t9 = t66 * t204 + t256;
t8 = -t268 - t283;
t7 = t184 + t269;
t1 = (t29 - t286) * t162 + (-t276 * t66 - t30) * t158;
t2 = [0, 0, -t227, -t165 * t245, 0, 0, 0, 0, 0, -t164 * t227 + (t95 - t199) * qJD(3), t160 * t227 + (-t94 - t198) * qJD(3), 0, 0, 0, 0, 0, t35 * t152 + (t114 * t236 - t165 * t200) * t156, -t34 * t152 + (-t116 * t236 - t165 * t80) * t156, t12 * t178 + t13 * t202 + t46 * t56 - t47 * t55, -t22 * t12 + t23 * t13 + t3 * t46 + t4 * t47 + (-t165 * t63 + t236 * t81) * t156, 0, 0, 0, 0, 0 (qJD(6) * t180 - t158 * t13 + t162 * t223) * t276 + t181 * t55 + t12 * t64 + t46 * t30 -(qJD(6) * t181 + t162 * t13 + t158 * t223) * t276 + t180 * t55 + t12 * t66 + t46 * t29; 0, 0, 0, 0, 0.2e1 * t160 * t213, -0.2e1 * t241 * t231, t243, -t244, 0, -pkin(8) * t243 + t160 * t173, pkin(8) * t244 + t164 * t173, t116 * t92 + t122 * t80, t92 * t114 + t116 * t93 - t80 * t121 - t122 * t200, -t92 * t152, -t93 * t152, 0, t109 * t93 + t110 * t121 + t176 * t114 + t149 * t200 + t152 * t289, -t109 * t92 + t110 * t122 - t176 * t116 + t149 * t80 + t152 * t290, t178 * t266 + t202 * t265 - t22 * t58 - t23 * t57 - t4 * t85 + t44 * t56 + t196, t63 * t185 - t266 * t22 + t265 * t23 + t277 * t81 + t3 * t44 + t4 * t45, -t66 * t229 + (t29 * t86 + t58 * t66) * t162 (-t158 * t66 - t162 * t64) * t58 + (-t256 - t162 * t30 + (t158 * t64 - t162 * t66) * qJD(6)) * t86, t162 * t193 - t229 * t276 + t29 * t85 + t66 * t57, -t158 * t193 - t228 * t86 - t30 * t85 - t64 * t57, t276 * t57 + t55 * t85, t14 * t85 + t44 * t30 - t191 * t57 + t266 * t64 + (t270 + t195 * t276 + (-t21 * t85 - t276 * t45 + t271) * qJD(6)) * t162 + t275 * t158, t44 * t29 - t6 * t57 + t266 * t66 + (-t270 - (-qJD(6) * t21 + t17) * t85 - qJD(6) * t271 + (qJD(6) * t45 - t195) * t276) * t158 + t275 * t162; 0, 0, 0, 0, -t160 * t167 * t164, t241 * t167, 0, 0, 0, t160 * t179, t164 * t179, -t250, t67, t61, t62, 0, -t208 * t152 - t114 * t150 + t251 + (t217 * t159 - t89) * qJD(4) + t190, t264 * t152 + t116 * t150 + (t217 * qJD(4) - t68) * t163 + t182, -t106 * t56 - t107 * t55 + t178 * t255 + t202 * t254 + t194, t4 * t107 - t3 * t106 - t81 * (t150 - t273) + t254 * t23 - t255 * t22, t9, t1, t8, t7, -t267, t102 * t30 + t255 * t64 + (-t3 - t282) * t162 + t174 * t158 + t226, t102 * t29 + t158 * t282 + t174 * t162 + t255 * t66 + t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t250, t67, t61, t62, 0, -t189 * t152 + t171 + t251, t209 * t152 + t182 - t279, -t24 * t178 - t25 * t202 + (-t155 * t55 - t253 * t56) * pkin(4) + t194, t22 * t24 - t23 * t25 + (t116 * t81 + t155 * t4 - t253 * t3) * pkin(4), t9, t1, t8, t7, -t267, t146 * t30 - t3 * t162 - (-t158 * t25 + t162 * t38) * t276 - t24 * t64 - t20 * t257 + (-t228 - t259) * t145 + t226, t146 * t29 + (t158 * t38 + t162 * t25) * t276 - t24 * t66 - t162 * t272 + t280 * t145 + t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178 ^ 2 - t202 ^ 2, t178 * t22 - t202 * t23 + t63, 0, 0, 0, 0, 0, t184 - t269, -t268 + t283; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t64, -t64 ^ 2 + t66 ^ 2, t29 + t286, -t288 * t66 - t258, t55, -t158 * t4 - t20 * t66 - t288 * t6 + t14, -t158 * t17 - t162 * t4 + t191 * t288 + t20 * t64;];
tauc_reg  = t2;
