% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:12:25
% EndTime: 2019-03-09 06:12:34
% DurationCPUTime: 3.89s
% Computational Cost: add. (9510->393), mult. (25040->492), div. (0->0), fcn. (19771->8), ass. (0->203)
t147 = cos(pkin(10));
t152 = cos(qJ(3));
t220 = t147 * t152;
t146 = sin(pkin(10));
t150 = sin(qJ(3));
t221 = t146 * t150;
t174 = -t220 + t221;
t120 = t174 * qJD(1);
t127 = t146 * t152 + t147 * t150;
t121 = t127 * qJD(1);
t149 = sin(qJ(4));
t258 = cos(qJ(4));
t166 = -t258 * t120 - t149 * t121;
t208 = qJD(5) - t166;
t187 = t208 ^ 2;
t204 = qJD(3) + qJD(4);
t252 = pkin(7) + qJ(2);
t132 = t252 * t146;
t128 = qJD(1) * t132;
t133 = t252 * t147;
t129 = qJD(1) * t133;
t176 = t128 * t150 - t129 * t152;
t92 = -pkin(8) * t120 - t176;
t88 = t149 * t92;
t269 = -t152 * t128 - t129 * t150;
t91 = -pkin(8) * t121 + t269;
t90 = qJD(3) * pkin(3) + t91;
t53 = t258 * t90 - t88;
t47 = -pkin(4) * t204 - t53;
t148 = sin(qJ(5));
t151 = cos(qJ(5));
t165 = -t149 * t120 + t258 * t121;
t93 = t148 * t165 - t151 * t204;
t95 = t148 * t204 + t151 * t165;
t30 = t93 * pkin(5) - t95 * qJ(6) + t47;
t289 = t208 * t30;
t89 = t258 * t92;
t54 = t149 * t90 + t89;
t48 = pkin(9) * t204 + t54;
t140 = -pkin(2) * t147 - pkin(1);
t130 = qJD(1) * t140 + qJD(2);
t109 = pkin(3) * t120 + t130;
t55 = -pkin(4) * t166 - pkin(9) * t165 + t109;
t22 = t148 * t55 + t151 * t48;
t11 = qJ(6) * t208 + t22;
t195 = qJD(4) * t258;
t214 = qJD(4) * t149;
t205 = qJD(1) * qJD(3);
t192 = t152 * t205;
t193 = t150 * t205;
t217 = t146 * t192 + t147 * t193;
t78 = -t217 * pkin(8) - qJD(2) * t120 + t269 * qJD(3);
t136 = t147 * t192;
t116 = -t146 * t193 + t136;
t163 = t127 * qJD(2);
t160 = qJD(1) * t163;
t79 = -pkin(8) * t116 + qJD(3) * t176 - t160;
t158 = -t149 * t79 - t90 * t195 + t92 * t214 - t258 * t78;
t212 = qJD(5) * t151;
t213 = qJD(5) * t148;
t161 = t258 * t116 - t149 * t217;
t285 = qJD(4) * t166;
t155 = t161 + t285;
t194 = t217 * pkin(3);
t188 = t149 * t116 + t258 * t217;
t266 = qJD(4) * t165;
t70 = t188 + t266;
t34 = t70 * pkin(4) - pkin(9) * t155 + t194;
t190 = -t148 * t158 - t151 * t34 + t48 * t212 + t55 * t213;
t259 = pkin(5) * t70;
t4 = t190 - t259;
t288 = t11 * t208 - t4;
t67 = t148 * t70;
t246 = t208 * t212 + t67;
t283 = t166 * t151;
t173 = -t208 * t283 + t246;
t241 = t165 * t95;
t287 = t173 - t241;
t210 = t166 * qJD(3);
t286 = -t210 + t161;
t186 = qJD(5) * t204;
t42 = t165 * t213 + (-t155 - t186) * t151;
t40 = t42 * t148;
t284 = -t40 + (t212 - t283) * t95;
t282 = t165 * t166;
t281 = t165 ^ 2 - t166 ^ 2;
t73 = pkin(4) * t165 - pkin(9) * t166;
t280 = -t109 * t166 + t158;
t181 = pkin(5) * t148 - qJ(6) * t151;
t279 = pkin(5) * t213 - qJ(6) * t212 - qJD(6) * t148 - t166 * t181;
t233 = t151 * t93;
t236 = t148 * t95;
t177 = t233 + t236;
t183 = -t42 * t151 - t95 * t213;
t154 = t148 * t155;
t228 = qJD(5) * t95;
t43 = t154 + t228;
t250 = -t148 * t43 - t93 * t212;
t278 = t166 * t177 + t183 + t250;
t169 = t148 * t34 - t151 * t158 + t55 * t212 - t48 * t213;
t244 = qJ(6) * t70;
t2 = qJD(6) * t208 + t169 + t244;
t275 = t4 * t148 + t2 * t151;
t185 = pkin(3) * t195;
t273 = pkin(5) * t165;
t242 = t165 * t93;
t69 = t151 * t70;
t272 = t208 * t213 - t69;
t16 = t149 * t78 + t92 * t195 + t90 * t214 - t258 * t79;
t7 = pkin(5) * t43 + qJ(6) * t42 - qJD(6) * t95 + t16;
t196 = -t7 * t151 + t30 * t213;
t191 = -t16 * t151 + t47 * t213;
t56 = t149 * t91 + t89;
t184 = pkin(3) * t214 - t56;
t175 = t132 * t150 - t133 * t152;
t100 = -pkin(8) * t174 - t175;
t222 = t132 * t152;
t99 = -pkin(8) * t127 - t133 * t150 - t222;
t64 = t149 * t100 - t258 * t99;
t271 = qJ(6) * t165;
t270 = t208 * t165;
t211 = t165 * qJD(3);
t268 = t211 - t188;
t21 = -t148 * t48 + t151 * t55;
t218 = qJD(6) - t21;
t10 = -pkin(5) * t208 + t218;
t267 = t10 * t148 + t11 * t151;
t265 = t10 * t165 + t196;
t264 = -t109 * t165 - t16;
t255 = t7 * t148;
t263 = -t11 * t165 - t255;
t262 = -t21 * t165 + t191;
t261 = t16 * t148 + t22 * t165 + t47 * t212;
t260 = t95 ^ 2;
t123 = t127 * qJD(3);
t257 = pkin(3) * t123;
t256 = t30 * t95;
t253 = t95 * t93;
t249 = t148 * t73 + t151 * t53;
t57 = t258 * t91 - t88;
t60 = pkin(3) * t121 + t73;
t248 = t148 * t60 + t151 * t57;
t65 = t258 * t100 + t149 * t99;
t108 = t258 * t127 - t149 * t174;
t112 = pkin(3) * t174 + t140;
t164 = -t149 * t127 - t174 * t258;
t66 = -pkin(4) * t164 - pkin(9) * t108 + t112;
t247 = t148 * t66 + t151 * t65;
t141 = pkin(3) * t149 + pkin(9);
t239 = t141 * t70;
t122 = t174 * qJD(3);
t76 = t164 * qJD(4) - t258 * t122 - t149 * t123;
t238 = t148 * t76;
t237 = t148 * t93;
t235 = t151 * t43;
t234 = t151 * t76;
t232 = t151 * t95;
t230 = t184 + t279;
t229 = t279 - t54;
t227 = t208 * t148;
t226 = t166 * t148;
t216 = t146 ^ 2 + t147 ^ 2;
t215 = qJD(3) * t121;
t207 = t10 * t212 + t275;
t206 = qJD(1) * qJD(2);
t203 = t258 * pkin(3);
t199 = qJD(1) * t221;
t197 = t141 * t213;
t189 = t216 * qJD(1) ^ 2;
t182 = t151 * pkin(5) + t148 * qJ(6);
t180 = t10 * t151 - t11 * t148;
t179 = -t166 * t47 - t239;
t178 = -t148 * t53 + t151 * t73;
t172 = 0.2e1 * t216 * t206;
t171 = t166 * t227 - t272;
t131 = -pkin(4) - t182;
t170 = t208 * t22 - t190;
t159 = -qJD(3) * t222 + qJD(2) * t220 + (-qJD(2) * t146 - qJD(3) * t133) * t150;
t80 = -pkin(8) * t123 + t159;
t156 = qJD(3) * t175 - t163;
t81 = pkin(8) * t122 + t156;
t25 = -t64 * qJD(4) + t149 * t81 + t258 * t80;
t77 = t108 * qJD(4) - t149 * t122 + t258 * t123;
t37 = pkin(4) * t77 - pkin(9) * t76 + t257;
t168 = t148 * t37 + t151 * t25 + t66 * t212 - t65 * t213;
t162 = t246 * pkin(9);
t157 = qJD(5) * t180 + t275;
t26 = t65 * qJD(4) + t149 * t80 - t258 * t81;
t142 = -t203 - pkin(4);
t125 = -t203 + t131;
t58 = pkin(5) * t95 + qJ(6) * t93;
t38 = t108 * t181 + t64;
t28 = pkin(5) * t164 + t148 * t65 - t151 * t66;
t27 = -qJ(6) * t164 + t247;
t23 = t208 * t93 - t42;
t20 = -t178 - t273;
t19 = t249 + t271;
t18 = t148 * t57 - t151 * t60 - t273;
t17 = t248 + t271;
t8 = t181 * t76 + (qJD(5) * t182 - qJD(6) * t151) * t108 + t26;
t6 = -pkin(5) * t77 + t247 * qJD(5) + t148 * t25 - t151 * t37;
t5 = qJ(6) * t77 - qJD(6) * t164 + t168;
t1 = [0, 0, 0, 0, 0, t172, qJ(2) * t172, t116 * t127 - t121 * t122, -t116 * t174 + t122 * t120 - t121 * t123 - t127 * t217, -t122 * qJD(3), -t123 * qJD(3), 0, t156 * qJD(3) + t130 * t123 + t140 * t217, -qJD(3) * t159 + t140 * t116 - t130 * t122, t108 * t155 + t165 * t76, -t108 * t70 + t155 * t164 - t165 * t77 + t166 * t76, t76 * t204, -t77 * t204, 0, t112 * t70 + t109 * t77 - t26 * t204 + (-t123 * t166 - t164 * t217) * pkin(3), t108 * t194 + t109 * t76 + t112 * t155 + t165 * t257 - t204 * t25, t183 * t108 + t76 * t232, -t177 * t76 + (t40 - t235 + (-t232 + t237) * qJD(5)) * t108, t108 * t69 + t164 * t42 + t77 * t95 + (-t108 * t213 + t234) * t208, -t108 * t67 + t164 * t43 - t77 * t93 + (-t108 * t212 - t238) * t208, -t164 * t70 + t208 * t77, t190 * t164 + t21 * t77 + t26 * t93 + t64 * t43 + ((-qJD(5) * t65 + t37) * t208 + t66 * t70 + t47 * qJD(5) * t108) * t151 + ((-qJD(5) * t66 - t25) * t208 - t65 * t70 + t16 * t108 + t47 * t76) * t148, -t191 * t108 + t164 * t169 - t168 * t208 - t22 * t77 + t47 * t234 - t247 * t70 + t26 * t95 - t64 * t42, t30 * t238 - t10 * t77 - t208 * t6 + t164 * t4 - t28 * t70 + t38 * t43 + t8 * t93 + (t212 * t30 + t255) * t108, -t27 * t43 - t28 * t42 - t5 * t93 + t6 * t95 + t180 * t76 + (-t267 * qJD(5) - t148 * t2 + t151 * t4) * t108, t196 * t108 + t11 * t77 - t164 * t2 + t208 * t5 - t30 * t234 + t27 * t70 + t38 * t42 - t8 * t95, t10 * t6 + t11 * t5 + t2 * t27 + t28 * t4 + t30 * t8 + t38 * t7; 0, 0, 0, 0, 0, -t189, -qJ(2) * t189, 0, 0, 0, 0, 0, t215 + t217, t136 + (-t120 - t199) * qJD(3), 0, 0, 0, 0, 0, t188 + t211 + 0.2e1 * t266, t161 + t210 + 0.2e1 * t285, 0, 0, 0, 0, 0, t171 - t242, -t151 * t187 - t241 - t67, -t148 * t187 - t242 + t69 (t233 - t236) * t166 - t183 + t250, t173 + t241, -t165 * t30 + t288 * t151 + (t10 * t208 + t2) * t148; 0, 0, 0, 0, 0, 0, 0, t121 * t120, -t120 ^ 2 + t121 ^ 2, t136 + (t120 - t199) * qJD(3), t215 - t217, 0, -t130 * t121 - t160, t130 * t120 + t174 * t206, -t282, t281, t286, t268, 0, t56 * t204 + (t121 * t166 - t204 * t214) * pkin(3) + t264, t57 * t204 + (-t121 * t165 - t195 * t204) * pkin(3) + t280, t284, t278, t287, t171 + t242, -t270, t142 * t43 + t184 * t93 + t179 * t148 + ((-qJD(5) * t141 - t60) * t151 + (-t185 + t57) * t148) * t208 + t262, -t142 * t42 + t184 * t95 + t179 * t151 + (-t151 * t185 + t197 + t248) * t208 + t261, t125 * t43 + t230 * t93 + (-t166 * t30 - t239) * t148 + (-t141 * t212 - t148 * t185 + t18) * t208 + t265, t17 * t93 - t18 * t95 + (-t93 * t185 - t10 * t166 + (-t43 + t228) * t141) * t151 + (t95 * t185 + t166 * t11 - t141 * t42 + (t141 * t93 - t11) * qJD(5)) * t148 + t207, t125 * t42 - t230 * t95 + (-t17 - t197) * t208 + (t185 * t208 + t239 - t289) * t151 + t263, -t10 * t18 - t11 * t17 + t7 * t125 + t157 * t141 + t267 * t185 + t230 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t282, t281, t286, t268, 0, t204 * t54 + t264, t204 * t53 + t280, t284, t278, t287, -t208 * t227 + t242 + t69, -t270, -pkin(4) * t43 - t178 * t208 - t226 * t47 - t54 * t93 - t162 + t262, pkin(4) * t42 + t272 * pkin(9) + t208 * t249 - t283 * t47 - t54 * t95 + t261, t131 * t43 + t20 * t208 - t226 * t30 + t229 * t93 - t162 + t265, -t11 * t213 + t19 * t93 - t20 * t95 - t180 * t166 + (-t40 - t235 + (t232 + t237) * qJD(5)) * pkin(9) + t207, t131 * t42 - t229 * t95 + (-pkin(9) * t213 - t19) * t208 + (pkin(9) * t70 - t289) * t151 + t263, pkin(9) * t157 - t10 * t20 - t11 * t19 + t131 * t7 + t229 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, -t93 ^ 2 + t260, t23, -t148 * t186 - t165 * t212 + t208 * t95 - t154, t70, -t47 * t95 + t170, t208 * t21 + t47 * t93 - t169, -t58 * t93 + t170 - t256 + 0.2e1 * t259, pkin(5) * t42 - qJ(6) * t43 + (t11 - t22) * t95 + (t10 - t218) * t93, 0.2e1 * t244 - t30 * t93 + t58 * t95 + (0.2e1 * qJD(6) - t21) * t208 + t169, -pkin(5) * t4 + qJ(6) * t2 - t10 * t22 + t11 * t218 - t30 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253 - t70, t23, -t260 - t187, t256 - t288;];
tauc_reg  = t1;
