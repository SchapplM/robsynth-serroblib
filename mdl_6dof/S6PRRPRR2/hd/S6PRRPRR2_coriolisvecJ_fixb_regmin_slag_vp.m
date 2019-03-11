% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:00:16
% EndTime: 2019-03-08 22:00:25
% DurationCPUTime: 3.27s
% Computational Cost: add. (3942->352), mult. (10271->511), div. (0->0), fcn. (8216->12), ass. (0->203)
t156 = sin(pkin(12));
t162 = sin(qJ(3));
t226 = qJD(2) * t162;
t158 = cos(pkin(12));
t166 = cos(qJ(3));
t236 = t158 * t166;
t125 = qJD(2) * t236 - t156 * t226;
t268 = qJD(5) + qJD(6);
t281 = t125 - t268;
t163 = sin(qJ(2));
t157 = sin(pkin(6));
t228 = qJD(1) * t157;
t209 = t163 * t228;
t224 = qJD(3) * t162;
t280 = pkin(3) * t224 - t209;
t134 = t156 * t166 + t158 * t162;
t126 = t134 * qJD(3);
t133 = t156 * t162 - t236;
t129 = t133 * qJD(3);
t279 = -pkin(4) * t126 - pkin(9) * t129 - t280;
t167 = cos(qJ(2));
t208 = t167 * t228;
t103 = t133 * t208;
t262 = -qJ(4) - pkin(8);
t201 = qJD(3) * t262;
t122 = t166 * qJD(4) + t162 * t201;
t123 = -t162 * qJD(4) + t166 * t201;
t70 = t122 * t158 + t123 * t156;
t248 = t70 + t103;
t127 = t134 * qJD(2);
t161 = sin(qJ(5));
t165 = cos(qJ(5));
t219 = t165 * qJD(3);
t106 = t127 * t161 - t219;
t138 = qJD(2) * pkin(8) + t209;
t198 = qJ(4) * qJD(2) + t138;
t159 = cos(pkin(6));
t227 = qJD(1) * t159;
t207 = t162 * t227;
t99 = t166 * t198 + t207;
t86 = t156 * t99;
t146 = t166 * t227;
t98 = -t162 * t198 + t146;
t91 = qJD(3) * pkin(3) + t98;
t34 = t158 * t91 - t86;
t31 = -qJD(3) * pkin(4) - t34;
t27 = pkin(5) * t106 + t31;
t108 = qJD(3) * t161 + t127 * t165;
t160 = sin(qJ(6));
t164 = cos(qJ(6));
t45 = t106 * t164 + t108 * t160;
t278 = t27 * t45;
t137 = t160 * t165 + t161 * t164;
t258 = t281 * t137;
t185 = t106 * t160 - t108 * t164;
t277 = t185 * t45;
t223 = qJD(5) * t161;
t244 = t125 * t161;
t276 = t223 - t244;
t275 = t185 ^ 2 - t45 ^ 2;
t119 = qJD(5) - t125;
t115 = qJD(6) + t119;
t220 = qJD(6) * t164;
t221 = qJD(6) * t160;
t218 = qJD(2) * qJD(3);
t203 = t166 * t218;
t204 = t162 * t218;
t117 = -t156 * t204 + t158 * t203;
t58 = qJD(5) * t219 + t117 * t165 - t127 * t223;
t59 = qJD(5) * t108 + t161 * t117;
t9 = -t106 * t220 - t108 * t221 - t160 * t59 + t164 * t58;
t274 = t115 * t45 + t9;
t116 = qJD(2) * t126;
t255 = t158 * t99;
t35 = t156 * t91 + t255;
t32 = qJD(3) * pkin(9) + t35;
t215 = -pkin(3) * t166 - pkin(2);
t118 = qJD(2) * t215 + qJD(4) - t208;
t51 = -t125 * pkin(4) - t127 * pkin(9) + t118;
t20 = t161 * t51 + t165 * t32;
t190 = qJD(4) + t208;
t269 = (-t138 * t166 - t207) * qJD(3) + (-qJ(4) * qJD(3) * t166 - t162 * t190) * qJD(2);
t57 = (-t138 * t162 + t146) * qJD(3) + (-qJ(4) * t224 + t166 * t190) * qJD(2);
t25 = t156 * t269 + t158 * t57;
t124 = pkin(3) * t204 + qJD(2) * t209;
t43 = pkin(4) * t116 - pkin(9) * t117 + t124;
t41 = t165 * t43;
t172 = -qJD(5) * t20 - t161 * t25 + t41;
t2 = t116 * pkin(5) - t58 * pkin(10) + t172;
t222 = qJD(5) * t165;
t181 = t161 * t43 + t165 * t25 + t222 * t51 - t223 * t32;
t3 = -pkin(10) * t59 + t181;
t214 = -t160 * t3 + t164 * t2;
t14 = -pkin(10) * t106 + t20;
t254 = t164 * t14;
t19 = -t161 * t32 + t165 * t51;
t13 = -pkin(10) * t108 + t19;
t8 = pkin(5) * t119 + t13;
t5 = t160 * t8 + t254;
t273 = -qJD(6) * t5 + t185 * t27 + t214;
t171 = qJD(6) * t185 - t160 * t58 - t164 * t59;
t272 = -t115 * t185 + t171;
t271 = -t161 * t103 - t165 * t279;
t141 = t262 * t162;
t142 = t262 * t166;
t105 = t141 * t156 - t142 * t158;
t84 = pkin(4) * t133 - pkin(9) * t134 + t215;
t270 = t105 * t223 + t161 * t279 - t165 * t248 - t84 * t222;
t249 = t156 * t122 - t123 * t158 - t134 * t208;
t78 = t137 * t134;
t232 = t165 * t129;
t179 = -t134 * t223 - t232;
t136 = t160 * t161 - t164 * t165;
t259 = t281 * t136;
t267 = -t115 * t259 - t137 * t116;
t11 = t14 * t221;
t206 = qJD(6) * t8 + t3;
t266 = t160 * t2 + t164 * t206 - t11;
t92 = t165 * t105;
t265 = pkin(10) * t232 + t126 * pkin(5) - t161 * t70 + (-t92 + (pkin(10) * t134 - t84) * t161) * qJD(5) + t271;
t210 = t134 * t222;
t233 = t161 * t129;
t180 = t210 - t233;
t264 = pkin(10) * t180 + t270;
t149 = pkin(3) * t156 + pkin(9);
t263 = pkin(10) + t149;
t39 = t158 * t98 - t86;
t71 = pkin(3) * t226 + pkin(4) * t127 - pkin(9) * t125;
t261 = t161 * t71 + t165 * t39;
t260 = t161 * t84 + t92;
t257 = qJD(2) * pkin(2);
t256 = t127 * t45;
t24 = t156 * t57 - t158 * t269;
t253 = t24 * t165;
t252 = t185 * t127;
t251 = t58 * t161;
t250 = pkin(5) * t180 + t249;
t247 = t106 * t119;
t246 = t108 * t119;
t245 = t108 * t127;
t243 = t127 * t106;
t242 = t134 * t161;
t241 = t134 * t165;
t239 = t157 * t163;
t238 = t157 * t167;
t169 = qJD(2) ^ 2;
t237 = t157 * t169;
t235 = t161 * t116;
t109 = t165 * t116;
t168 = qJD(3) ^ 2;
t231 = t168 * t162;
t230 = t168 * t166;
t229 = t162 ^ 2 - t166 ^ 2;
t225 = qJD(2) * t163;
t216 = t163 * t237;
t150 = -pkin(3) * t158 - pkin(4);
t213 = t157 * t225;
t212 = qJD(2) * t238;
t37 = t156 * t98 + t255;
t200 = qJD(5) * t263;
t104 = -t141 * t158 - t142 * t156;
t199 = t119 * t165;
t197 = t162 * t212;
t196 = t166 * t212;
t195 = t115 * t258 - t136 * t116;
t194 = pkin(5) * t276 - t37;
t131 = t263 * t161;
t193 = -pkin(10) * t244 + qJD(6) * t131 + t161 * t200 + t261;
t132 = t263 * t165;
t64 = t165 * t71;
t192 = t127 * pkin(5) + qJD(6) * t132 - t161 * t39 + t64 + (-pkin(10) * t125 + t200) * t165;
t81 = t165 * t84;
t26 = pkin(5) * t133 - pkin(10) * t241 - t105 * t161 + t81;
t28 = -pkin(10) * t242 + t260;
t189 = t160 * t26 + t164 * t28;
t130 = t159 * t162 + t166 * t239;
t182 = t159 * t166 - t162 * t239;
t77 = t130 * t158 + t156 * t182;
t183 = t161 * t238 - t165 * t77;
t55 = -t161 * t77 - t165 * t238;
t188 = t160 * t183 + t164 * t55;
t187 = t160 * t55 - t164 * t183;
t186 = -t105 * t116 + t24 * t134;
t184 = -t119 * t276 + t109;
t177 = t257 * qJD(2);
t176 = -t116 * t149 + t119 * t31;
t174 = -0.2e1 * qJD(3) * t257;
t140 = -pkin(5) * t165 + t150;
t97 = -qJD(3) * t130 - t197;
t96 = qJD(3) * t182 + t196;
t85 = t116 * t133;
t79 = t136 * t134;
t76 = t130 * t156 - t158 * t182;
t68 = pkin(5) * t242 + t104;
t38 = t156 * t97 + t158 * t96;
t36 = t156 * t96 - t158 * t97;
t22 = -t221 * t242 + (t241 * t268 - t233) * t164 + t179 * t160;
t21 = t136 * t129 - t268 * t78;
t17 = qJD(5) * t183 - t161 * t38 + t165 * t213;
t16 = qJD(5) * t55 + t161 * t213 + t165 * t38;
t12 = pkin(5) * t59 + t24;
t4 = -t14 * t160 + t164 * t8;
t1 = [0, 0, -t216, -t167 * t237, 0, 0, 0, 0, 0, -t166 * t216 + (t97 - t197) * qJD(3), t162 * t216 + (-t96 - t196) * qJD(3), -t116 * t77 + t117 * t76 + t125 * t38 + t127 * t36, t24 * t76 + t25 * t77 - t34 * t36 + t35 * t38 + (t118 * t225 - t124 * t167) * t157, 0, 0, 0, 0, 0, t106 * t36 + t116 * t55 + t119 * t17 + t59 * t76, t108 * t36 + t116 * t183 - t119 * t16 + t58 * t76, 0, 0, 0, 0, 0 (-qJD(6) * t187 - t160 * t16 + t164 * t17) * t115 + t188 * t116 + t36 * t45 - t76 * t171 -(qJD(6) * t188 + t164 * t16 + t160 * t17) * t115 - t187 * t116 - t36 * t185 + t76 * t9; 0, 0, 0, 0, 0.2e1 * t162 * t203, -0.2e1 * t229 * t218, t230, -t231, 0, -pkin(8) * t230 + t162 * t174, pkin(8) * t231 + t166 * t174, t104 * t117 + t125 * t248 - t35 * t126 + t127 * t249 + t34 * t129 - t25 * t133 + t186, t24 * t104 + t25 * t105 + t118 * t280 + t124 * t215 + t248 * t35 - t249 * t34, t108 * t179 + t241 * t58 -(-t106 * t165 - t108 * t161) * t129 + (-t251 - t165 * t59 + (t106 * t161 - t108 * t165) * qJD(5)) * t134, t108 * t126 + t109 * t134 + t119 * t179 + t58 * t133, -t106 * t126 - t119 * t180 - t59 * t133 - t134 * t235, t119 * t126 + t85, t81 * t116 + (-t222 * t32 + t41) * t133 + t19 * t126 + t104 * t59 + t31 * t210 + (-t105 * t222 + t271) * t119 + t249 * t106 + ((-qJD(5) * t84 - t70) * t119 + (-qJD(5) * t51 - t25) * t133 - t31 * t129 + t186) * t161, -t260 * t116 - t181 * t133 - t20 * t126 + t104 * t58 - t31 * t232 + (-t223 * t31 + t253) * t134 + t270 * t119 + t249 * t108, -t185 * t21 - t79 * t9, -t171 * t79 + t185 * t22 - t21 * t45 - t78 * t9, t115 * t21 - t116 * t79 - t126 * t185 + t133 * t9, -t115 * t22 - t116 * t78 - t126 * t45 + t133 * t171, t115 * t126 + t85 (-t160 * t28 + t164 * t26) * t116 + t214 * t133 + t4 * t126 - t68 * t171 + t12 * t78 + t27 * t22 + t250 * t45 + (t160 * t264 + t164 * t265) * t115 + (-t115 * t189 - t133 * t5) * qJD(6), -t189 * t116 - t266 * t133 - t5 * t126 + t68 * t9 - t12 * t79 + t27 * t21 - t250 * t185 + ((-qJD(6) * t26 + t264) * t164 + (qJD(6) * t28 - t265) * t160) * t115; 0, 0, 0, 0, -t162 * t169 * t166, t229 * t169, 0, 0, 0, t162 * t177, t166 * t177 (t35 - t37) * t127 + (t34 - t39) * t125 + (-t116 * t156 - t117 * t158) * pkin(3), t34 * t37 - t35 * t39 + (-t118 * t226 + t156 * t25 - t158 * t24) * pkin(3), t108 * t199 + t251 (t58 - t247) * t165 + (-t59 - t246) * t161, t119 * t199 + t235 - t245, t184 + t243, -t119 * t127, -t37 * t106 - t19 * t127 + t150 * t59 - t253 + (-t149 * t222 - t64) * t119 + (t119 * t39 + t176) * t161, -t37 * t108 + t20 * t127 + t150 * t58 + t24 * t161 + (t149 * t223 + t261) * t119 + t176 * t165, t9 * t137 - t185 * t259, -t9 * t136 + t137 * t171 - t185 * t258 - t259 * t45, t252 - t267, t195 + t256, -t115 * t127 (-t131 * t164 - t132 * t160) * t116 - t140 * t171 + t12 * t136 - t4 * t127 + t194 * t45 - t258 * t27 + (t160 * t193 - t164 * t192) * t115 -(-t131 * t160 + t132 * t164) * t116 + t140 * t9 + t12 * t137 + t5 * t127 - t194 * t185 + t259 * t27 + (t160 * t192 + t164 * t193) * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125 ^ 2 - t127 ^ 2, -t125 * t35 + t127 * t34 + t124, 0, 0, 0, 0, 0, t184 - t243, -t119 ^ 2 * t165 - t235 - t245, 0, 0, 0, 0, 0, t195 - t256, t252 + t267; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108 * t106, -t106 ^ 2 + t108 ^ 2, t58 + t247, t246 - t59, t116, -t31 * t108 + t20 * t119 + t172, t106 * t31 + t119 * t19 - t181, -t277, t275, t274, t272, t116 -(-t13 * t160 - t254) * t115 + (-t108 * t45 - t115 * t221 + t116 * t164) * pkin(5) + t273, t278 + t11 + (-t115 * t14 - t2) * t160 + (t115 * t13 - t206) * t164 + (t108 * t185 - t115 * t220 - t116 * t160) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t277, t275, t274, t272, t116, t5 * t115 + t273, t4 * t115 - t266 + t278;];
tauc_reg  = t1;
