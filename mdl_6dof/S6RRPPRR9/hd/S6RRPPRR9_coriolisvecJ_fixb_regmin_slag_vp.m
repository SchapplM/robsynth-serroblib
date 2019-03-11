% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRR9_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:32:07
% EndTime: 2019-03-09 09:32:20
% DurationCPUTime: 4.27s
% Computational Cost: add. (3379->435), mult. (8947->594), div. (0->0), fcn. (6438->8), ass. (0->222)
t158 = sin(pkin(6));
t167 = cos(qJ(2));
t239 = qJD(1) * t167;
t217 = t158 * t239;
t122 = qJD(5) + t217;
t159 = cos(pkin(6));
t229 = t159 * qJD(1);
t143 = qJD(2) + t229;
t163 = sin(qJ(5));
t166 = cos(qJ(5));
t164 = sin(qJ(2));
t240 = qJD(1) * t158;
t221 = t164 * t240;
t94 = t166 * t143 + t163 * t221;
t263 = t94 * t122;
t227 = qJD(1) * qJD(2);
t215 = t158 * t227;
t125 = t167 * t215;
t62 = qJD(5) * t94 - t166 * t125;
t286 = -t62 - t263;
t115 = t166 * t221;
t92 = t163 * t143 - t115;
t91 = qJD(6) + t92;
t226 = pkin(1) * t229;
t244 = -pkin(8) * t221 + t167 * t226;
t228 = qJD(3) - t244;
t162 = sin(qJ(6));
t165 = cos(qJ(6));
t55 = -t165 * t122 + t162 * t94;
t285 = t122 * t55;
t268 = t122 * t92;
t234 = qJD(5) * t163;
t61 = qJD(5) * t115 + t163 * t125 - t143 * t234;
t284 = -t61 + t268;
t101 = pkin(8) * t217 + t164 * t226;
t84 = pkin(3) * t217 + t101;
t248 = qJD(4) + t84;
t133 = t143 * pkin(2);
t283 = -t143 * qJ(4) - t133;
t142 = t143 ^ 2;
t156 = t164 ^ 2;
t155 = t158 ^ 2;
t259 = t155 * qJD(1) ^ 2;
t282 = -t156 * t259 - t142;
t274 = pkin(1) * t164;
t148 = t159 * t274;
t257 = t158 * t167;
t281 = pkin(8) * t257 + t148;
t128 = t143 * qJ(3);
t249 = pkin(4) * t217 + t248;
t32 = -t143 * pkin(9) + t128 + t249;
t160 = qJ(3) - pkin(9);
t161 = pkin(2) + qJ(4);
t200 = -t161 * t167 - pkin(1);
t66 = (-t160 * t164 + t200) * t158;
t45 = qJD(1) * t66;
t14 = t163 * t32 + t166 * t45;
t235 = qJD(3) * t164;
t169 = (-t235 + (-t160 * qJD(2) - qJD(4)) * t167) * t158;
t202 = t164 * t215;
t116 = pkin(2) * t202;
t247 = qJ(4) * t202 + t116;
t27 = qJD(1) * t169 + t247;
t258 = t158 * t164;
t206 = (-pkin(3) - pkin(4)) * t258;
t189 = qJD(1) * t206;
t124 = t143 * qJD(3);
t275 = pkin(1) * t159;
t225 = qJD(2) * t275;
t205 = qJD(1) * t225;
t246 = pkin(8) * t202 - t167 * t205;
t73 = -t124 + t246;
t36 = qJD(2) * t189 - t73;
t171 = -qJD(5) * t14 - t163 * t27 + t166 * t36;
t4 = pkin(5) * t202 - t171;
t280 = (t94 * pkin(5) + t91 * pkin(10)) * t91 + t4;
t265 = t165 * t62;
t279 = -t162 * t91 ^ 2 + t265;
t57 = t162 * t122 + t165 * t94;
t18 = qJD(6) * t57 + t162 * t61 + t165 * t202;
t238 = qJD(2) * t164;
t64 = t189 + t244;
t251 = qJD(3) - t64;
t31 = -t251 - t283;
t278 = (t160 * t238 - t167 * t31) * t240 - t31 * qJD(5);
t96 = -t159 * qJ(3) - t281;
t80 = pkin(3) * t257 - t96;
t54 = pkin(4) * t257 - t159 * pkin(9) + t80;
t191 = t163 * t54 + t166 * t66;
t220 = t158 * t238;
t136 = pkin(2) * t220;
t245 = qJ(4) * t220 + t136;
t35 = t169 + t245;
t276 = pkin(3) + pkin(8);
t210 = t158 * (-pkin(4) - t276);
t140 = t167 * t225;
t151 = t159 * qJD(3);
t243 = t140 + t151;
t47 = t210 * t238 + t243;
t277 = -qJD(5) * t191 - t163 * t35 + t166 * t47;
t95 = pkin(8) * t125 + t164 * t205;
t50 = pkin(3) * t125 - t143 * qJD(4) + t95;
t37 = -pkin(4) * t125 - t50;
t10 = t62 * pkin(5) - t61 * pkin(10) + t37;
t12 = t122 * pkin(10) + t14;
t16 = t92 * pkin(5) - t94 * pkin(10) + t31;
t196 = t162 * t12 - t165 * t16;
t232 = qJD(5) * t166;
t181 = t163 * t36 + t166 * t27 + t32 * t232 - t234 * t45;
t3 = -pkin(10) * t202 + t181;
t1 = -qJD(6) * t196 + t162 * t10 + t165 * t3;
t273 = t55 * t91;
t272 = t55 * t94;
t271 = t57 * t91;
t270 = t57 * t94;
t129 = qJ(4) * t221;
t141 = pkin(2) * t221;
t63 = -t160 * t217 + t129 + t141;
t269 = t163 * t64 + t166 * t63;
t267 = t162 * t62;
t266 = t162 * t91;
t212 = t165 * t91;
t230 = qJD(6) * t165;
t17 = t122 * t230 + t165 * t61 + (-qJD(6) * t94 - t202) * t162;
t264 = t17 * t162;
t199 = pkin(5) * t166 + pkin(10) * t163;
t262 = qJD(5) * t199 - (-pkin(4) - t199) * t217 + t248;
t261 = qJ(3) * t164;
t260 = t122 * t166;
t256 = t160 * t162;
t255 = t160 * t165;
t254 = t163 * t167;
t253 = t165 * t167;
t83 = -pkin(3) * t221 + t244;
t250 = qJD(3) - t83;
t157 = t167 ^ 2;
t242 = t156 - t157;
t241 = qJ(3) * qJD(2);
t237 = qJD(2) * t166;
t236 = qJD(2) * t167;
t233 = qJD(5) * t165;
t231 = qJD(6) * t162;
t224 = t122 * t254;
t222 = t167 * t259;
t98 = -t159 * pkin(2) + pkin(8) * t258 - t167 * t275;
t219 = t158 * t236;
t218 = t122 * t234;
t216 = t155 * t227;
t112 = t163 * pkin(5) - t166 * pkin(10) + t161;
t211 = -pkin(10) * t221 - qJD(6) * t112 + t269;
t209 = t143 + t229;
t208 = 0.2e1 * t216;
t207 = -qJD(4) - t241;
t204 = t164 * t222;
t203 = t159 * qJ(4) - t98;
t201 = -0.2e1 * pkin(1) * t216;
t198 = t101 * t143 - t95;
t6 = t165 * t12 + t162 * t16;
t20 = pkin(10) * t257 + t191;
t104 = t159 * t166 + t163 * t258;
t182 = -t159 * t163 + t166 * t258;
t53 = t206 + t203;
t25 = -pkin(5) * t182 - t104 * pkin(10) + t53;
t195 = t162 * t25 + t165 * t20;
t194 = -t162 * t20 + t165 * t25;
t13 = -t163 * t45 + t166 * t32;
t192 = -t163 * t66 + t166 * t54;
t102 = t281 * qJD(2);
t190 = t102 * t143 + t95 * t159;
t188 = -pkin(8) * t220 + t140;
t99 = -qJ(3) * t217 + t141;
t86 = -t143 * t217 + t125;
t187 = (-qJD(2) + t143) * t221;
t186 = -t230 * t91 - t267;
t185 = -t231 * t91 + t265;
t183 = -t104 * t162 + t158 * t253;
t77 = t104 * t165 + t162 * t257;
t180 = t163 * t47 + t166 * t35 + t54 * t232 - t234 * t66;
t179 = t91 * t122;
t177 = t122 * t57;
t97 = (-pkin(2) * t167 - pkin(1) - t261) * t158;
t176 = t212 * t91 + t267;
t175 = (-qJ(3) * t236 - t235) * t158;
t11 = -t122 * pkin(5) - t13;
t173 = -pkin(10) * t62 + (t11 + t13) * t91;
t49 = -pkin(3) * t202 - t73;
t79 = (t200 - t261) * t158;
t2 = -qJD(6) * t6 + t165 * t10 - t162 * t3;
t172 = (t167 * t207 - t235) * t158;
t150 = t159 * qJD(4);
t48 = t150 + (t167 * t210 - t148) * qJD(2);
t126 = qJD(3) * t217;
t108 = t163 * t202;
t90 = -t151 - t188;
t89 = (-t162 * t164 + t163 * t253) * t240;
t88 = (t162 * t254 + t164 * t165) * t240;
t87 = qJD(1) * t97;
t85 = t136 + t175;
t82 = -t128 - t101;
t81 = t129 + t99;
t78 = -t133 + t228;
t75 = qJD(5) * t104 - t166 * t219;
t74 = qJD(5) * t182 + t163 * t219;
t72 = qJD(1) * t175 + t116;
t71 = -t150 + (t276 * t257 + t148) * qJD(2);
t70 = -t276 * t220 + t243;
t69 = pkin(3) * t258 - t203;
t68 = t87 * t221;
t67 = qJD(1) * t79;
t60 = t128 + t248;
t46 = t67 * t217;
t43 = t172 + t245;
t40 = t250 + t283;
t38 = qJD(1) * t172 + t247;
t24 = qJD(6) * t183 - t162 * t220 + t74 * t165;
t23 = qJD(6) * t77 + t74 * t162 + t165 * t220;
t21 = pkin(5) * t221 + t163 * t63 - t166 * t64;
t19 = -pkin(5) * t257 - t192;
t15 = t75 * pkin(5) - t74 * pkin(10) + t48;
t8 = pkin(5) * t220 - t277;
t7 = -pkin(10) * t220 + t180;
t5 = [0, 0, 0, t164 * t167 * t208, -t242 * t208, t209 * t219, -t209 * t220, 0, t164 * t201 - t190, -t143 * t188 + t159 * t246 + t167 * t201 (t164 * t95 - t167 * t73 + (t164 * t82 + t167 * t78) * qJD(2) + (t102 * t164 - t167 * t90 + (t164 * t96 + t167 * t98) * qJD(2)) * qJD(1)) * t158 (-t87 * t238 + t167 * t72 + (t167 * t85 - t238 * t97) * qJD(1)) * t158 + t190, -t90 * t143 - t73 * t159 + (-t87 * t236 - t164 * t72 + (-t164 * t85 - t236 * t97) * qJD(1)) * t158, t102 * t78 + t72 * t97 + t73 * t96 + t82 * t90 + t85 * t87 + t95 * t98 (t164 * t50 + t167 * t49 + (-t164 * t60 + t167 * t40) * qJD(2) + (t164 * t71 + t167 * t70 + (-t164 * t80 + t167 * t69) * qJD(2)) * qJD(1)) * t158, t70 * t143 + t49 * t159 + (-t67 * t236 - t164 * t38 + (-t164 * t43 - t236 * t79) * qJD(1)) * t158, -t71 * t143 - t50 * t159 + (t67 * t238 - t167 * t38 + (-t167 * t43 + t238 * t79) * qJD(1)) * t158, t38 * t79 + t40 * t71 + t43 * t67 + t49 * t80 + t50 * t69 + t60 * t70, t104 * t61 + t74 * t94, -t104 * t62 + t182 * t61 - t74 * t92 - t75 * t94, t74 * t122 + (t167 * t61 + (-qJD(1) * t104 - t94) * t238) * t158, -t75 * t122 + (-t167 * t62 + (-qJD(1) * t182 + t92) * t238) * t158 (-t122 * t158 - t155 * t239) * t238, t277 * t122 + t48 * t92 + t53 * t62 - t37 * t182 + t31 * t75 + (t171 * t167 + (-qJD(1) * t192 - t13) * t238) * t158, -t180 * t122 + t48 * t94 + t53 * t61 + t37 * t104 + t31 * t74 + (-t181 * t167 + (qJD(1) * t191 + t14) * t238) * t158, t17 * t77 + t24 * t57, t17 * t183 - t18 * t77 - t23 * t57 - t24 * t55, -t17 * t182 + t24 * t91 + t57 * t75 + t62 * t77, t18 * t182 + t183 * t62 - t23 * t91 - t55 * t75, -t182 * t62 + t75 * t91 (-qJD(6) * t195 + t165 * t15 - t162 * t7) * t91 + t194 * t62 - t2 * t182 - t196 * t75 + t8 * t55 + t19 * t18 - t4 * t183 + t11 * t23 -(qJD(6) * t194 + t162 * t15 + t165 * t7) * t91 - t195 * t62 + t1 * t182 - t6 * t75 + t8 * t57 + t19 * t17 + t4 * t77 + t11 * t24; 0, 0, 0, -t204, t242 * t259, t86, t187, 0, t259 * t274 + t198, pkin(1) * t222 + t143 * t244 + t246, t126 + ((-pkin(2) * qJD(2) - t244 - t78) * t167 + (-t101 - t82 - t241) * t164) * t240, -t217 * t99 - t198 + t68, t228 * t143 + (t164 * t99 + t167 * t87) * t240 - t73, -t95 * pkin(2) - t73 * qJ(3) - t78 * t101 - t228 * t82 - t87 * t99, t126 + ((-qJD(2) * t161 - t40 - t83) * t167 + (t207 + t60 - t84) * t164) * t240, -t83 * t143 + 0.2e1 * t124 + t46 + (-pkin(3) * qJD(2) + t81) * t221 - t246, t248 * t143 + (-t164 * t67 + t167 * t81) * t240 - t50, t49 * qJ(3) - t50 * t161 - t248 * t40 + t250 * t60 - t67 * t81, -t163 * t263 + t61 * t166, t284 * t163 + t286 * t166, -t218 + (-t224 + (t94 - t237) * t164) * t240, -t122 * t232 + t108 + (-t164 * t92 - t167 * t260) * t240, t122 * t221, t13 * t221 + t161 * t62 + t249 * t92 + (t37 + (-qJD(5) * t160 + t63) * t122) * t163 + (t122 * t251 - t278) * t166, -t14 * t221 + t161 * t61 + t37 * t166 + t249 * t94 + (-t160 * t232 + t269) * t122 + (-qJD(3) * t122 + t278) * t163, t17 * t165 * t166 + (-t163 * t233 - t166 * t231 - t89) * t57, t89 * t55 + t57 * t88 + (t162 * t57 + t165 * t55) * t234 + (-t264 - t165 * t18 + (t162 * t55 - t165 * t57) * qJD(6)) * t166, -t89 * t91 + (-t233 * t91 + t17) * t163 + (t177 + t185) * t166, t88 * t91 + (qJD(5) * t266 - t18) * t163 + (t186 - t285) * t166, t62 * t163 + t260 * t91, t112 * t265 - t11 * t88 - t21 * t55 + (t162 * t211 + t165 * t262) * t91 + ((-qJD(3) * t162 - t160 * t230) * t91 - t62 * t256 + t2 + (-t11 * t162 + t160 * t55) * qJD(5)) * t163 + (-t196 * t217 + t11 * t230 - qJD(3) * t55 - t160 * t18 + t4 * t162 + (-t256 * t91 - t196) * qJD(5)) * t166, -t112 * t267 - t11 * t89 - t21 * t57 + (-t162 * t262 + t165 * t211) * t91 + (-(qJD(3) * t165 - t160 * t231) * t91 - t62 * t255 - t1 + (-t11 * t165 + t160 * t57) * qJD(5)) * t163 + (-t6 * t217 - t11 * t231 - qJD(3) * t57 - t160 * t17 + t4 * t165 + (-t255 * t91 - t6) * qJD(5)) * t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t204, t282, t143 * t82 + t68 + t95, t86, t282, -t204, -t143 * t60 + t221 * t67 + t50, 0, 0, 0, 0, 0, t286, t284, 0, 0, 0, 0, 0, t272 - t279, t176 + t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187, -t204, -t157 * t259 - t142, t143 * t40 + t46 + t49, 0, 0, 0, 0, 0, -t218 - t143 * t92 + (-t164 * t237 - t224) * t240, -t122 ^ 2 * t166 - t143 * t94 + t108, 0, 0, 0, 0, 0, -t143 * t212 + (-t162 * t179 - t18) * t166 + (t186 + t285) * t163, t143 * t266 + (-t165 * t179 - t17) * t166 + (t177 - t185) * t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94 * t92, -t92 ^ 2 + t94 ^ 2, t61 + t268, -t62 + t263, -t202, t14 * t122 - t31 * t94 + t171, t122 * t13 + t31 * t92 - t181, t212 * t57 + t264 (t17 - t273) * t165 + (-t18 - t271) * t162, t176 - t270, t272 + t279, -t91 * t94, -pkin(5) * t18 - t14 * t55 + t173 * t162 - t280 * t165 + t196 * t94, -pkin(5) * t17 - t14 * t57 + t280 * t162 + t173 * t165 + t6 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 * t55, -t55 ^ 2 + t57 ^ 2, t17 + t273, -t18 + t271, t62, -t11 * t57 + t6 * t91 + t2, t11 * t55 - t196 * t91 - t1;];
tauc_reg  = t5;
