% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPPR5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:23:44
% EndTime: 2019-03-09 08:23:55
% DurationCPUTime: 3.59s
% Computational Cost: add. (2224->421), mult. (5618->552), div. (0->0), fcn. (3602->6), ass. (0->225)
t159 = cos(qJ(2));
t230 = qJD(1) * t159;
t133 = -qJD(6) + t230;
t281 = qJD(6) + t133;
t280 = pkin(4) + qJ(3);
t221 = qJD(1) * qJD(2);
t279 = -0.2e1 * t221;
t265 = pkin(3) + qJ(5);
t154 = sin(pkin(9));
t157 = sin(qJ(2));
t231 = qJD(1) * t157;
t213 = t154 * t231;
t155 = cos(pkin(9));
t222 = t155 * qJD(2);
t95 = t213 - t222;
t278 = t265 * t95;
t204 = t159 * t221;
t121 = t154 * t204;
t277 = qJ(5) * t121 + t95 * qJD(5);
t211 = t155 * t231;
t229 = qJD(2) * t154;
t97 = t211 + t229;
t268 = t97 ^ 2;
t91 = t95 ^ 2;
t276 = t268 + t91;
t197 = t95 + t222;
t59 = t197 * t230;
t275 = t265 * t155;
t226 = qJD(4) * t154;
t140 = pkin(7) * t230;
t212 = t154 * t230;
t235 = pkin(3) * t212 + t140;
t274 = -qJD(5) * t155 - t226 - t235;
t242 = t155 * t157;
t273 = pkin(4) * t242 + t159 * qJ(5);
t225 = qJD(4) * t159;
t228 = qJD(2) * t157;
t272 = -qJ(4) * t228 + t225;
t153 = t159 ^ 2;
t161 = qJD(1) ^ 2;
t246 = t153 * t161;
t271 = -t268 - t246;
t156 = sin(qJ(6));
t158 = cos(qJ(6));
t194 = t155 * t204;
t47 = t156 * t97 + t158 * t95;
t17 = qJD(6) * t47 + t121 * t156 - t158 * t194;
t227 = qJD(2) * t159;
t142 = pkin(7) * t227;
t209 = t154 * t227;
t234 = pkin(3) * t209 + t142;
t270 = -t157 * (qJD(4) * t155 - qJD(5) * t154) + t234;
t224 = qJD(6) * t156;
t247 = t133 * t159;
t269 = qJD(1) * (-t156 * t247 + t158 * t228) + t133 * t224;
t267 = pkin(4) + pkin(8);
t266 = t95 * t97;
t264 = -pkin(5) - qJ(4);
t99 = t154 * t158 + t155 * t156;
t170 = t159 * t99;
t81 = t99 * qJD(6);
t263 = -qJD(1) * t170 + t81;
t210 = t155 * t230;
t241 = t155 * t158;
t100 = t154 * t156 - t241;
t82 = t100 * qJD(6);
t262 = -t156 * t212 + t158 * t210 + t82;
t139 = pkin(7) * t231;
t107 = (qJD(3) - t139) * qJD(2);
t184 = pkin(2) * t157 - qJ(3) * t159;
t79 = qJD(2) * t184 - t157 * qJD(3);
t71 = t79 * qJD(1);
t38 = t155 * t107 + t154 * t71;
t261 = qJD(2) * pkin(2);
t45 = t156 * t95 - t158 * t97;
t260 = t133 * t45;
t259 = t133 * t47;
t258 = t155 * t79;
t119 = qJD(2) * qJ(3) + t140;
t252 = qJ(3) * t157;
t111 = -pkin(2) * t159 - pkin(1) - t252;
t86 = t111 * qJD(1);
t49 = -t154 * t119 + t155 * t86;
t39 = pkin(3) * t230 + qJD(4) - t49;
t257 = t157 * t39;
t50 = t155 * t119 + t154 * t86;
t41 = qJ(4) * t230 - t50;
t256 = t157 * t41;
t251 = qJ(5) * t154;
t174 = t264 * t155 + t251;
t165 = t159 * t174;
t255 = -qJD(1) * t165 + t274;
t180 = qJ(4) * t155 - t251;
t171 = t159 * t180;
t254 = -qJD(1) * t171 - t274;
t240 = t155 * t159;
t131 = pkin(7) * t240;
t68 = t154 * t111 + t131;
t136 = qJ(4) * t231;
t103 = t184 * qJD(1);
t84 = t154 * t103;
t253 = t136 + t84;
t250 = qJD(3) * t95;
t249 = qJD(3) * t97;
t248 = qJD(4) * t97;
t245 = t154 * t157;
t244 = t154 * t159;
t243 = t155 * t103;
t239 = t159 * t161;
t160 = qJD(2) ^ 2;
t238 = t160 * t157;
t237 = t160 * t159;
t132 = pkin(7) * t204;
t236 = pkin(3) * t121 + t132;
t233 = pkin(3) * t245 + t157 * pkin(7);
t115 = t280 * t154;
t116 = t280 * t155;
t232 = t157 ^ 2 - t153;
t223 = qJD(6) * t158;
t146 = t159 * qJD(5);
t220 = pkin(8) * t240;
t129 = pkin(7) * t244;
t134 = t157 * t221;
t219 = qJ(4) * t134 + t38;
t218 = pkin(7) * t228;
t217 = t154 * t267;
t216 = t97 * t230;
t215 = -pkin(7) * t154 - pkin(3);
t37 = -t154 * t107 + t155 * t71;
t179 = pkin(4) * t194 + qJD(1) * t146 - t37;
t206 = t265 * t157;
t4 = (-t206 + t220) * t221 + t179;
t195 = t159 * t217;
t5 = (-t225 + (pkin(5) * t157 - t195) * qJD(2)) * qJD(1) + t219;
t214 = -t156 * t4 + t158 * t5;
t208 = t159 * t222;
t205 = t264 * t159;
t203 = t154 * qJ(4) + pkin(2);
t202 = -t134 + t266;
t198 = -qJD(3) + t261;
t186 = -t139 + t198;
t201 = t186 - t261;
t168 = qJ(4) * t97 + t186;
t19 = t168 - t278;
t83 = t203 + t275;
t200 = -qJD(2) * t83 + t19;
t199 = pkin(1) * t279;
t196 = pkin(4) * t210 - t243;
t67 = t111 * t155 - t129;
t193 = qJD(2) * t206;
t192 = qJD(1) * t205;
t61 = qJ(4) * t159 - t68;
t191 = pkin(4) * t208 + t146 - t258;
t190 = t157 * t215;
t189 = -pkin(4) * t229 - qJD(4);
t57 = t121 - t216;
t188 = t156 * t5 + t158 * t4;
t172 = qJ(5) * t230 + t39;
t7 = t267 * t97 + t172;
t8 = -t267 * t95 + qJD(5) + t192 + t50;
t2 = t156 * t8 + t158 * t7;
t187 = t156 * t7 - t158 * t8;
t72 = t154 * t79;
t185 = t72 - t272;
t151 = t159 * pkin(3);
t62 = t151 - t67;
t15 = pkin(4) * t97 + t172;
t20 = -pkin(4) * t95 + qJD(5) - t41;
t183 = -t15 * t155 + t154 * t20;
t29 = t129 + t151 + (pkin(8) * t157 - t111) * t155 + t273;
t30 = -t157 * t217 + t205 + t68;
t182 = t156 * t30 + t158 * t29;
t181 = -t156 * t29 + t158 * t30;
t64 = -pkin(7) * t211 + t84;
t56 = -t155 * t218 + t72;
t177 = (-qJ(5) + t215) * t157;
t164 = t177 + t220;
t90 = pkin(8) * t155 + t116;
t176 = -qJD(1) * t164 + qJD(3) * t154 + qJD(6) * t90 - t196;
t163 = (-pkin(7) * t155 + pkin(5)) * t157 - t195;
t89 = pkin(8) * t154 + t115;
t175 = -qJD(1) * t163 + qJD(3) * t155 - qJD(6) * t89 - t253;
t16 = t158 * t121 + t156 * t194 + t97 * t223 - t224 * t95;
t173 = -pkin(4) * t244 - pkin(7) * t242;
t166 = -qJ(4) * t194 + t236;
t31 = t166 - t248;
t162 = (-t156 * t228 - t158 * t247) * qJD(1) + t133 * t223;
t120 = qJD(3) * t212;
t108 = -pkin(3) * t155 - t203;
t77 = t99 * t157;
t76 = t156 * t245 - t157 * t241;
t74 = -qJ(4) * t242 + t233;
t73 = t264 * t154 - pkin(2) - t275;
t66 = -qJ(4) * t210 + t235;
t63 = pkin(7) * t213 + t243;
t60 = t157 * t180 - t233;
t58 = (-t95 + t222) * t230;
t55 = t154 * t218 + t258;
t54 = qJD(1) * t190 - t243;
t53 = -t136 - t64;
t52 = (-qJ(4) * t227 - qJD(4) * t157) * t155 + t234;
t48 = -pkin(4) * t245 - t61;
t44 = t157 * t174 + t233;
t43 = qJD(2) * t190 - t258;
t42 = t62 + t273;
t40 = qJD(1) * t173 + t253;
t36 = -t56 + t272;
t35 = pkin(3) * t95 - t168;
t33 = qJD(2) * t170 - t157 * t82;
t32 = t156 * t209 + t157 * t81 - t158 * t208;
t28 = qJD(1) * t177 + t196;
t27 = -pkin(3) * t134 - t37;
t26 = qJD(2) * t171 - t270;
t25 = qJD(2) * t173 + t185;
t24 = qJD(1) * t225 - t219;
t21 = qJD(2) * t177 + t191;
t18 = qJD(2) * t165 + t270;
t14 = qJD(2) * t163 + t185;
t13 = qJD(2) * t164 + t191;
t12 = t189 * t230 + t219;
t11 = -t31 - t277;
t10 = -qJD(1) * t193 + t179;
t9 = t264 * t97 - t186 + t278;
t6 = t192 * t222 + t236 - t248 + t277;
t1 = [0, 0, 0, 0.2e1 * t159 * t134, t232 * t279, t237, -t238, 0, -pkin(7) * t237 + t157 * t199, pkin(7) * t238 + t159 * t199 (-qJD(1) * t55 - t37) * t159 + ((pkin(7) * t95 - t154 * t186) * t159 + (t49 + (t67 + 0.2e1 * t129) * qJD(1)) * t157) * qJD(2) (qJD(1) * t56 + t38) * t159 + ((pkin(7) * t97 - t155 * t186) * t159 + (-t50 + (-t68 + 0.2e1 * t131) * qJD(1)) * t157) * qJD(2), -t55 * t97 - t56 * t95 + (-t154 * t38 - t155 * t37) * t157 + (-t154 * t50 - t155 * t49 + (-t154 * t68 - t155 * t67) * qJD(1)) * t227, t37 * t67 + t38 * t68 + t49 * t55 + t50 * t56 + (-t186 + t139) * t142, t36 * t95 + t43 * t97 + (t154 * t24 + t155 * t27) * t157 + (t154 * t41 + t155 * t39 + (t154 * t61 + t155 * t62) * qJD(1)) * t227, -t31 * t245 - t52 * t95 + (-qJD(1) * t43 - t27) * t159 + (-t35 * t244 + t257 + (t157 * t62 - t244 * t74) * qJD(1)) * qJD(2), -t31 * t242 - t52 * t97 + (qJD(1) * t36 + t24) * t159 + (-t35 * t240 - t256 + (-t157 * t61 - t240 * t74) * qJD(1)) * qJD(2), t24 * t61 + t27 * t62 + t31 * t74 + t35 * t52 + t36 * t41 + t39 * t43, t11 * t242 + t26 * t97 + (-qJD(1) * t25 - t12) * t159 + (t19 * t240 + t157 * t20 + (t157 * t48 + t240 * t60) * qJD(1)) * qJD(2), -t21 * t97 + t25 * t95 + (-t10 * t155 + t12 * t154) * t157 + ((t154 * t48 - t155 * t42) * qJD(1) + t183) * t227, -t11 * t245 - t26 * t95 + (qJD(1) * t21 + t10) * t159 + (-t19 * t244 - t15 * t157 + (-t157 * t42 - t244 * t60) * qJD(1)) * qJD(2), t10 * t42 + t11 * t60 + t12 * t48 + t15 * t21 + t19 * t26 + t20 * t25, t16 * t77 + t33 * t47, -t16 * t76 - t17 * t77 - t32 * t47 - t33 * t45, -t133 * t33 - t159 * t16 + (qJD(1) * t77 + t47) * t228, t133 * t32 + t159 * t17 + (-qJD(1) * t76 - t45) * t228 (-t133 - t230) * t228 -(-t156 * t13 + t158 * t14) * t133 - t214 * t159 + t18 * t45 + t44 * t17 + t6 * t76 + t9 * t32 + (t133 * t182 + t159 * t2) * qJD(6) + (qJD(1) * t181 - t187) * t228 (t158 * t13 + t156 * t14) * t133 + t188 * t159 + t18 * t47 + t44 * t16 + t6 * t77 + t9 * t33 + (t133 * t181 - t159 * t187) * qJD(6) + (-qJD(1) * t182 - t2) * t228; 0, 0, 0, -t157 * t239, t232 * t161, 0, 0, 0, t161 * pkin(1) * t157, pkin(1) * t239, t120 + ((-qJ(3) * t229 - t49) * t157 + (-pkin(7) * t197 + t154 * t201 + t63) * t159) * qJD(1) ((-qJ(3) * t222 + t50) * t157 + (-t64 + (-t97 + t229) * pkin(7) + (t186 - t198) * t155) * t159) * qJD(1), t63 * t97 + t64 * t95 + (t230 * t49 - t250 + t38) * t155 + (t230 * t50 + t249 - t37) * t154, -t49 * t63 - t50 * t64 + (-t154 * t49 + t155 * t50) * qJD(3) + (-t154 * t37 + t155 * t38) * qJ(3) + t201 * t140, -t53 * t95 - t54 * t97 + (-t230 * t39 - t24 - t250) * t155 + (-t230 * t41 + t249 + t27) * t154, t155 * t31 - t120 + (t66 + t226) * t95 + (-t257 + t159 * t54 + (t159 * t35 + (-t108 * t159 + t252) * qJD(2)) * t154) * qJD(1), t66 * t97 + (-t31 + t248) * t154 + (t256 - t159 * t53 + (qJ(3) * t228 + (-qJD(2) * t108 - qJD(3) + t35) * t159) * t155) * qJD(1), t108 * t31 - t35 * t66 - t39 * t54 - t41 * t53 + (-qJ(3) * t24 - qJD(3) * t41) * t155 + (qJ(3) * t27 + qJD(3) * t39 - qJD(4) * t35) * t154, t11 * t154 + t254 * t97 + ((qJD(2) * t116 - t20) * t157 + (t40 + (-qJD(3) - t200) * t155) * t159) * qJD(1), -t10 * t154 - t12 * t155 + t28 * t97 - t40 * t95 + (-t154 * t97 + t155 * t95) * qJD(3) + ((-t115 * t155 + t116 * t154) * qJD(2) - t183) * t230, t11 * t155 + t120 - t254 * t95 + ((-qJD(2) * t115 + t15) * t157 + (t154 * t200 - t28) * t159) * qJD(1), t10 * t115 + t11 * t83 + t116 * t12 - t15 * t28 - t20 * t40 + t254 * t19 + (t15 * t154 + t155 * t20) * qJD(3), t16 * t100 + t263 * t47, -t100 * t17 + t16 * t99 - t262 * t47 - t263 * t45, -t263 * t133 + (qJD(2) * t100 - t47) * t231, t262 * t133 + (qJD(2) * t99 + t45) * t231, t133 * t231, t73 * t17 - t6 * t99 + t262 * t9 + t255 * t45 + (t156 * t176 - t158 * t175) * t133 + ((-t156 * t89 + t158 * t90) * qJD(2) + t187) * t231, t6 * t100 + t73 * t16 + t263 * t9 + t255 * t47 + (t156 * t175 + t158 * t176) * t133 + (-(t156 * t90 + t158 * t89) * qJD(2) + t2) * t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t59, -t276, t49 * t97 + t50 * t95 + t132, -t276, -t57, -t59, -t41 * t95 + (-qJD(4) - t39) * t97 + t166, -t59, t276, t57, t20 * t95 + (-qJD(4) - t15) * t97 + t166 + t277, 0, 0, 0, 0, 0, t17 - t259, t16 + t260; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t202, t271, t35 * t97 + (-pkin(3) * t228 - t159 * t41) * qJD(1) - t37, t271, -t58, t202, -t19 * t97 + (t159 * t20 - t193) * qJD(1) + t179, 0, 0, 0, 0, 0, t45 * t97 + t162, t47 * t97 - t269; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134 + t266, t121 + t216, -t91 - t246, t19 * t95 + (-t15 + t189) * t230 + t219, 0, 0, 0, 0, 0, -t45 * t95 + t269, -t47 * t95 + t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t45, -t45 ^ 2 + t47 ^ 2, t16 - t260, -t17 - t259, t134, -t281 * t2 - t47 * t9 + t214, t281 * t187 + t45 * t9 - t188;];
tauc_reg  = t1;
