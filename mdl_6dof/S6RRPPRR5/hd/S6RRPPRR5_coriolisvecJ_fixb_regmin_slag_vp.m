% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPRR5
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
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:11:04
% EndTime: 2019-03-09 09:11:15
% DurationCPUTime: 4.31s
% Computational Cost: add. (3353->438), mult. (8861->592), div. (0->0), fcn. (6386->8), ass. (0->215)
t162 = sin(qJ(6));
t160 = cos(pkin(6));
t233 = t160 * qJD(1);
t140 = qJD(2) + t233;
t166 = cos(qJ(5));
t163 = sin(qJ(5));
t164 = sin(qJ(2));
t159 = sin(pkin(6));
t243 = qJD(1) * t159;
t224 = t164 * t243;
t204 = t163 * t224;
t282 = t166 * t140 + t204;
t82 = qJD(6) + t282;
t212 = t82 ^ 2;
t165 = cos(qJ(6));
t236 = qJD(5) * t166;
t167 = cos(qJ(2));
t239 = qJD(2) * t167;
t179 = t163 * t239 + t164 * t236;
t174 = t179 * t159;
t237 = qJD(5) * t163;
t220 = t140 * t237;
t55 = qJD(1) * t174 - t220;
t267 = t165 * t55;
t288 = -t162 * t212 + t267;
t242 = qJD(1) * t167;
t219 = t159 * t242;
t284 = -t219 - qJD(5);
t86 = -t163 * t140 + t166 * t224;
t50 = t162 * t86 + t165 * t284;
t287 = t50 * t86;
t52 = -t162 * t284 + t165 * t86;
t286 = t52 * t86;
t264 = t86 * t284;
t265 = t282 * t284;
t232 = qJD(1) * qJD(2);
t217 = t159 * t232;
t202 = t167 * t217;
t54 = -qJD(5) * t204 - t140 * t236 + t166 * t202;
t285 = -t54 - t265;
t257 = t159 * t164;
t135 = qJD(4) * t257;
t250 = qJ(4) * t202 + qJD(1) * t135;
t256 = t159 * t167;
t276 = pkin(1) * t164;
t181 = pkin(8) * t256 + t160 * t276;
t280 = t181 * qJD(2);
t89 = qJD(1) * t280;
t47 = t89 - t250;
t222 = t159 * t239;
t59 = -qJ(4) * t222 - t135 + t280;
t283 = t284 * t52;
t168 = -pkin(2) - pkin(3);
t281 = qJD(2) * t168;
t185 = (pkin(4) * t167 - pkin(9) * t164) * t159;
t79 = -pkin(1) * t243 - pkin(2) * t219 - qJ(3) * t224;
t60 = pkin(3) * t219 + qJD(4) - t79;
t34 = qJD(1) * t185 + t60;
t118 = t140 * qJ(3);
t231 = pkin(1) * t233;
t95 = pkin(8) * t219 + t164 * t231;
t76 = -qJ(4) * t219 + t95;
t56 = t118 + t76;
t42 = -t140 * pkin(9) + t56;
t14 = t163 * t34 + t166 * t42;
t156 = pkin(4) - t168;
t178 = t159 * (-pkin(9) * t167 - t156 * t164);
t173 = qJD(2) * t178;
t136 = qJD(3) * t257;
t249 = qJ(3) * t202 + qJD(1) * t136;
t28 = qJD(1) * t173 + t249;
t203 = t164 * t217;
t109 = qJ(4) * t203;
t133 = t167 * t231;
t113 = qJD(2) * t133;
t115 = t140 * qJD(3);
t238 = qJD(4) * t167;
t241 = qJD(2) * t164;
t35 = t109 + t113 + t115 + (-pkin(8) * t241 - t238) * t243;
t172 = -qJD(5) * t14 - t163 * t35 + t166 * t28;
t4 = pkin(5) * t203 - t172;
t279 = (t86 * pkin(5) + t82 * pkin(10)) * t82 + t4;
t18 = qJD(6) * t52 + t162 * t54 + t165 * t203;
t91 = -t159 * pkin(1) - pkin(2) * t256 - qJ(3) * t257;
t72 = pkin(3) * t256 - t91;
t45 = t185 + t72;
t90 = t160 * qJ(3) + t181;
t71 = -qJ(4) * t256 + t90;
t61 = -t160 * pkin(9) + t71;
t191 = t163 * t45 + t166 * t61;
t246 = qJ(3) * t222 + t136;
t33 = t173 + t246;
t277 = pkin(1) * t160;
t230 = qJD(2) * t277;
t134 = t167 * t230;
t148 = t160 * qJD(3);
t48 = t134 + t148 + (-t238 + (-pkin(8) + qJ(4)) * t241) * t159;
t278 = -qJD(5) * t191 - t163 * t48 + t166 * t33;
t10 = t55 * pkin(5) - t54 * pkin(10) - t47;
t12 = -pkin(10) * t284 + t14;
t119 = qJ(4) * t224;
t94 = -pkin(8) * t224 + t133;
t251 = qJD(3) - t94;
t70 = -t140 * pkin(2) + t251;
t40 = -t140 * pkin(3) - t119 + t70;
t31 = t140 * pkin(4) - t40;
t15 = pkin(5) * t282 - t86 * pkin(10) + t31;
t197 = t162 * t12 - t165 * t15;
t184 = t163 * t28 + t166 * t35 + t34 * t236 - t42 * t237;
t3 = -pkin(10) * t203 + t184;
t1 = -t197 * qJD(6) + t162 * t10 + t165 * t3;
t275 = t50 * t82;
t274 = t52 * t82;
t253 = t166 * t167;
t80 = (t162 * t253 + t164 * t165) * t243;
t273 = t80 * t82;
t81 = (-t162 * t164 + t165 * t253) * t243;
t272 = t81 * t82;
t122 = qJ(3) * t219;
t43 = qJD(1) * t178 + t122;
t74 = t119 + t94;
t271 = t163 * t43 + t166 * t74;
t270 = pkin(8) * qJD(2);
t269 = t162 * t55;
t268 = t164 * t86;
t234 = qJD(6) * t165;
t17 = -t284 * t234 + t165 * t54 + (-qJD(6) * t86 - t203) * t162;
t266 = t17 * t162;
t263 = t94 * t140;
t262 = t95 * t140;
t261 = t76 + t284 * (pkin(5) * t163 - pkin(10) * t166);
t260 = qJD(5) * t82;
t259 = t284 * t163;
t155 = t159 ^ 2;
t258 = t155 * qJD(1) ^ 2;
t161 = qJ(3) - pkin(9);
t255 = t161 * t162;
t254 = t161 * t165;
t252 = qJD(3) - t74;
t248 = 0.2e1 * t115 + t113;
t157 = t164 ^ 2;
t158 = t167 ^ 2;
t245 = t157 - t158;
t244 = qJ(3) * qJD(2);
t240 = qJD(2) * t166;
t235 = qJD(6) * t162;
t229 = t162 * t260;
t228 = t165 * t260;
t227 = t167 * t259;
t226 = t284 * t253;
t225 = t167 * t258;
t92 = -t160 * pkin(2) + pkin(8) * t257 - t167 * t277;
t223 = t159 * t241;
t221 = t284 * t237;
t218 = t155 * t232;
t206 = t168 * t257;
t75 = qJD(1) * t206 + t122;
t216 = -t75 - t270;
t93 = pkin(2) * t224 - t122;
t215 = t93 - t270;
t211 = t165 * t82;
t105 = t166 * pkin(5) + t163 * pkin(10) + t156;
t210 = -pkin(10) * t224 - qJD(6) * t105 + t271;
t209 = t140 + t233;
t208 = t164 * t230;
t207 = 0.2e1 * t218;
t205 = t164 * t225;
t201 = -0.2e1 * pkin(1) * t218;
t199 = t163 * t203 + t236 * t284;
t6 = t165 * t12 + t162 * t15;
t196 = -t140 * t280 - t89 * t160;
t20 = pkin(10) * t256 + t191;
t62 = -t160 * pkin(3) - qJ(4) * t257 + t92;
t53 = t160 * pkin(4) - t62;
t98 = t160 * t166 + t163 * t257;
t99 = -t160 * t163 + t166 * t257;
t25 = t98 * pkin(5) - t99 * pkin(10) + t53;
t195 = t162 * t25 + t165 * t20;
t194 = -t162 * t20 + t165 * t25;
t13 = -t163 * t42 + t166 * t34;
t192 = -t163 * t61 + t166 * t45;
t190 = qJD(2) * t206;
t189 = -pkin(8) * t223 + t134;
t188 = -t82 * t234 - t269;
t187 = t82 * t235 - t267;
t186 = -t99 * t162 + t165 * t256;
t69 = t162 * t256 + t99 * t165;
t183 = t163 * t33 + t166 * t48 + t45 * t236 - t61 * t237;
t182 = t284 * t50;
t180 = -pkin(8) * t203 + t113;
t11 = pkin(5) * t284 - t13;
t175 = -pkin(10) * t55 + (t11 + t13) * t82;
t2 = -qJD(6) * t6 + t165 * t10 - t162 * t3;
t171 = (t161 * t241 - t167 * t31) * t243 - t31 * qJD(5);
t170 = -t179 * t243 + t220;
t100 = t140 * t219;
t97 = -t140 ^ 2 - t157 * t258;
t83 = t148 + t189;
t78 = -t100 + t202;
t77 = pkin(2) * t223 - t246;
t73 = t118 + t95;
t67 = -qJD(5) * t98 + t166 * t222;
t66 = -t160 * t237 + t174;
t65 = t115 + t180;
t63 = pkin(2) * t203 - t249;
t58 = t190 + t246;
t46 = qJD(1) * t190 + t249;
t24 = qJD(6) * t186 - t162 * t223 + t67 * t165;
t23 = qJD(6) * t69 + t67 * t162 + t165 * t223;
t21 = pkin(5) * t224 + t163 * t74 - t166 * t43;
t19 = -pkin(5) * t256 - t192;
t16 = t66 * pkin(5) - t67 * pkin(10) - t59;
t8 = pkin(5) * t223 - t278;
t7 = -pkin(10) * t223 + t183;
t5 = [0, 0, 0, t164 * t167 * t207, -t245 * t207, t209 * t222, -t209 * t223, 0, t164 * t201 + t196, -t140 * t189 - t160 * t180 + t167 * t201 (t79 * t241 - t167 * t63 + (-t167 * t77 + t91 * t241) * qJD(1)) * t159 + t196 (t164 * t89 + t167 * t65 + (-t164 * t73 + t167 * t70) * qJD(2) + (t164 * t280 + t167 * t83 + (-t164 * t90 + t167 * t92) * qJD(2)) * qJD(1)) * t159, t83 * t140 + t65 * t160 + (-t79 * t239 - t164 * t63 + (-t164 * t77 - t91 * t239) * qJD(1)) * t159, t280 * t70 + t63 * t91 + t65 * t90 + t73 * t83 + t79 * t77 + t89 * t92, -t59 * t140 - t47 * t160 + (-t60 * t241 + t167 * t46 + (t167 * t58 - t72 * t241) * qJD(1)) * t159, t48 * t140 + t35 * t160 + (t60 * t239 + t164 * t46 + (t164 * t58 + t72 * t239) * qJD(1)) * t159 (-t164 * t47 - t167 * t35 + (t164 * t56 - t167 * t40) * qJD(2) + (-t164 * t59 - t167 * t48 + (t164 * t71 - t167 * t62) * qJD(2)) * qJD(1)) * t159, t35 * t71 + t40 * t59 + t46 * t72 + t47 * t62 + t56 * t48 + t60 * t58, t54 * t99 + t86 * t67, -t282 * t67 - t54 * t98 - t99 * t55 - t86 * t66, -t67 * t284 + (t167 * t54 + (-qJD(1) * t99 - t86) * t241) * t159, t66 * t284 + (-t167 * t55 + (qJD(1) * t98 + t282) * t241) * t159 (-t155 * t242 + t159 * t284) * t241, -t278 * t284 - t59 * t282 + t53 * t55 - t47 * t98 + t31 * t66 + (t172 * t167 + (-qJD(1) * t192 - t13) * t241) * t159, t183 * t284 - t59 * t86 + t53 * t54 - t47 * t99 + t31 * t67 + (-t184 * t167 + (qJD(1) * t191 + t14) * t241) * t159, t17 * t69 + t52 * t24, t17 * t186 - t69 * t18 - t52 * t23 - t24 * t50, t17 * t98 + t24 * t82 + t52 * t66 + t69 * t55, -t18 * t98 + t186 * t55 - t23 * t82 - t50 * t66, t55 * t98 + t82 * t66 (-qJD(6) * t195 + t165 * t16 - t162 * t7) * t82 + t194 * t55 + t2 * t98 - t197 * t66 + t8 * t50 + t19 * t18 - t4 * t186 + t11 * t23 -(qJD(6) * t194 + t162 * t16 + t165 * t7) * t82 - t195 * t55 - t1 * t98 - t6 * t66 + t8 * t52 + t19 * t17 + t4 * t69 + t11 * t24; 0, 0, 0, -t205, t245 * t258, t78 (-qJD(2) + t140) * t224, 0, -pkin(8) * t202 + t262 + (-t160 * t232 + t258) * t276, pkin(1) * t225 - t180 + t263, t262 + (-t208 + (-t164 * t79 + t167 * t215) * t159) * qJD(1) ((t73 - t95 - t244) * t164 + (-pkin(2) * qJD(2) + t251 - t70) * t167) * t243, -t263 + (t164 * t215 + t167 * t79) * t243 + t248, -t89 * pkin(2) + t65 * qJ(3) + t251 * t73 - t70 * t95 - t79 * t93, t76 * t140 + (-t208 + (t164 * t60 + t167 * t216) * t159) * qJD(1) + t250, -t74 * t140 + t109 + ((-qJD(4) - t60) * t167 + t216 * t164) * t243 + t248 ((-t56 + t76 + t244) * t164 + (-t252 + t40 - t281) * t167) * t243, t35 * qJ(3) + t47 * t168 + t252 * t56 - t40 * t76 - t60 * t75, -t54 * t163 + t166 * t264, t285 * t166 + (t55 - t264) * t163 (t226 + t268) * t243 + t199, -t221 + (-t227 + (-t282 + t240) * t164) * t243, -t284 * t224, t13 * t224 + t156 * t55 + t76 * t282 + (-t47 - (-qJD(5) * t161 - t43) * t284) * t166 + (t252 * t284 + t171) * t163, -t14 * t224 + t156 * t54 + t47 * t163 + t76 * t86 - (t161 * t237 + t271) * t284 + (qJD(3) * t284 + t171) * t166, -t17 * t165 * t163 + (t163 * t235 - t165 * t236 - t81) * t52, t81 * t50 + t52 * t80 + (t162 * t52 + t165 * t50) * t236 + (t266 + t165 * t18 + (-t162 * t50 + t165 * t52) * qJD(6)) * t163, -t272 + (t17 - t228) * t166 + (t187 + t283) * t163, t273 + (-t18 + t229) * t166 + (-t182 - t188) * t163, t55 * t166 + t259 * t82, t105 * t267 - t11 * t80 - t21 * t50 + (t162 * t210 + t165 * t261) * t82 + ((-qJD(3) * t162 - t161 * t234) * t82 - t55 * t255 + t2 + (-t11 * t162 + t161 * t50) * qJD(5)) * t166 + (t197 * t219 - t11 * t234 + qJD(3) * t50 + t161 * t18 - t4 * t162 + (t255 * t82 + t197) * qJD(5)) * t163, -t105 * t269 - t11 * t81 - t21 * t52 + (-t162 * t261 + t165 * t210) * t82 + (-(qJD(3) * t165 - t161 * t235) * t82 - t55 * t254 - t1 + (-t11 * t165 + t161 * t52) * qJD(5)) * t166 + (t6 * t219 + t11 * t235 + qJD(3) * t52 + t161 * t17 - t4 * t165 + (t254 * t82 + t6) * qJD(5)) * t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t205, t78, t97, -t73 * t140 + (t79 * t257 + t280) * qJD(1), -t205, t97, -t78, -t56 * t140 + (-t257 * t60 + t280) * qJD(1) - t250, 0, 0, 0, 0, 0, t170 + t264, t285, 0, 0, 0, 0, 0, t287 - t288, t165 * t212 + t269 + t286; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-qJD(2) - t140) * t224, t100 + t202 (-t157 - t158) * t258 (t167 * t56 + (t40 + t281) * t164) * t243 + t249, 0, 0, 0, 0, 0, t221 + (t227 + (-t282 - t240) * t164) * t243 (t226 - t268) * t243 + t199, 0, 0, 0, 0, 0, -t273 + (-t18 - t229) * t166 + (-t182 + t188) * t163, -t272 + (-t17 - t228) * t166 + (t187 - t283) * t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86 * t282, -t282 ^ 2 + t86 ^ 2, t54 - t265, t170 - t264, -t203, -t14 * t284 - t31 * t86 + t172, -t13 * t284 + t282 * t31 - t184, t211 * t52 + t266 (t17 - t275) * t165 + (-t18 - t274) * t162, t211 * t82 + t269 - t286, t287 + t288, -t82 * t86, -pkin(5) * t18 - t14 * t50 + t175 * t162 - t165 * t279 + t197 * t86, -pkin(5) * t17 - t14 * t52 + t162 * t279 + t175 * t165 + t6 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t50, -t50 ^ 2 + t52 ^ 2, t17 + t275, -t18 + t274, t55, -t11 * t52 + t6 * t82 + t2, t11 * t50 - t197 * t82 - t1;];
tauc_reg  = t5;
