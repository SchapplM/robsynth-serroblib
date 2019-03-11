% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:38:21
% EndTime: 2019-03-08 21:38:30
% DurationCPUTime: 3.27s
% Computational Cost: add. (4191->373), mult. (10687->520), div. (0->0), fcn. (8026->10), ass. (0->189)
t164 = sin(qJ(3));
t167 = cos(qJ(3));
t191 = pkin(3) * t164 - qJ(4) * t167;
t111 = t191 * qJD(3) - t164 * qJD(4);
t159 = sin(pkin(11));
t226 = qJD(3) * t164;
t218 = pkin(8) * t226;
t144 = t159 * t218;
t161 = cos(pkin(11));
t165 = sin(qJ(2));
t160 = sin(pkin(6));
t230 = qJD(1) * t160;
t168 = cos(qJ(2));
t238 = t167 * t168;
t253 = t161 * t111 + t144 - (-t159 * t238 + t161 * t165) * t230;
t277 = t159 * t111 - (t159 * t165 + t161 * t238) * t230;
t239 = t161 * t167;
t186 = pkin(4) * t164 - pkin(9) * t239;
t180 = t186 * qJD(3);
t276 = t180 + t253;
t240 = t161 * t164;
t244 = t159 * t167;
t275 = (-pkin(8) * t240 - pkin(9) * t244) * qJD(3) + t277;
t220 = t167 * qJD(2);
t151 = -qJD(5) + t220;
t221 = t161 * qJD(3);
t227 = qJD(2) * t164;
t125 = -t159 * t227 + t221;
t222 = t159 * qJD(3);
t126 = t161 * t227 + t222;
t163 = sin(qJ(5));
t166 = cos(qJ(5));
t71 = -t166 * t125 + t163 * t126;
t248 = t71 * t151;
t219 = qJD(2) * qJD(3);
t207 = t167 * t219;
t193 = t159 * t207;
t196 = t166 * t167 * t221;
t223 = qJD(5) * t166;
t39 = (qJD(5) * t126 + t193) * t163 - qJD(2) * t196 - t125 * t223;
t274 = t39 - t248;
t273 = t71 ^ 2;
t128 = t163 * t159 - t166 * t161;
t224 = qJD(5) * t163;
t266 = -t159 * t224 + t161 * t223;
t233 = t128 * t220 + t266;
t129 = t166 * t159 + t163 * t161;
t114 = t129 * qJD(5);
t181 = t129 * t167;
t232 = -qJD(2) * t181 + t114;
t188 = t163 * t125 + t166 * t126;
t264 = t188 ^ 2;
t136 = -t167 * pkin(3) - t164 * qJ(4) - pkin(2);
t124 = t161 * t136;
t78 = -pkin(9) * t240 + t124 + (-pkin(8) * t159 - pkin(4)) * t167;
t103 = pkin(8) * t239 + t159 * t136;
t245 = t159 * t164;
t89 = -pkin(9) * t245 + t103;
t272 = t276 * t163 + t275 * t166 + t78 * t223 - t89 * t224;
t130 = t191 * qJD(2);
t210 = t165 * t230;
t134 = qJD(2) * pkin(8) + t210;
t162 = cos(pkin(6));
t229 = qJD(1) * t162;
t267 = -t164 * t134 + t167 * t229;
t58 = t161 * t130 - t159 * t267;
t41 = t186 * qJD(2) + t58;
t216 = t159 * t220;
t59 = t159 * t130 + t161 * t267;
t48 = -pkin(9) * t216 + t59;
t259 = pkin(9) + qJ(4);
t139 = t259 * t159;
t140 = t259 * t161;
t91 = -t163 * t139 + t166 * t140;
t271 = t129 * qJD(4) + t91 * qJD(5) - t163 * t48 + t166 * t41;
t254 = t163 * t78 + t166 * t89;
t270 = t254 * qJD(5) + t275 * t163 - t276 * t166;
t187 = -t166 * t139 - t163 * t140;
t269 = -t128 * qJD(4) + t187 * qJD(5) - t163 * t41 - t166 * t48;
t268 = t151 * t188;
t265 = t188 * qJD(5);
t263 = qJ(6) * t226 - t167 * qJD(6) + t272;
t262 = -pkin(5) * t226 + t270;
t204 = -qJD(3) * pkin(3) + qJD(4);
t92 = t204 - t267;
t60 = -t125 * pkin(4) + t92;
t15 = t71 * pkin(5) - qJ(6) * t188 + t60;
t261 = t15 * t188;
t260 = t188 * t71;
t258 = qJ(6) * t227 - t269;
t257 = pkin(5) * t227 + t271;
t147 = t164 * t229;
t100 = t167 * t134 + t147;
t80 = pkin(4) * t216 + t100;
t256 = -pkin(5) * t232 + qJ(6) * t233 + t129 * qJD(6) + t80;
t209 = t168 * t230;
t195 = qJD(2) * t209;
t64 = t167 * t195 + (qJD(4) + t267) * qJD(3);
t81 = (t111 + t210) * qJD(2);
t26 = t159 * t81 + t161 * t64;
t101 = t136 * qJD(2) - t209;
t96 = qJD(3) * qJ(4) + t100;
t43 = t159 * t101 + t161 * t96;
t203 = t161 * t218;
t252 = -t203 + t277;
t251 = qJD(2) * pkin(2);
t225 = qJD(3) * t167;
t67 = qJD(3) * t147 + t134 * t225 + t164 * t195;
t250 = t67 * t159;
t249 = t67 * t161;
t247 = qJD(3) * t187;
t246 = qJD(3) * t91;
t243 = t160 * t165;
t242 = t160 * t168;
t170 = qJD(2) ^ 2;
t241 = t160 * t170;
t169 = qJD(3) ^ 2;
t237 = t169 * t164;
t236 = t169 * t167;
t42 = t161 * t101 - t159 * t96;
t27 = -pkin(4) * t220 - t126 * pkin(9) + t42;
t34 = t125 * pkin(9) + t43;
t10 = -t163 * t34 + t166 * t27;
t235 = qJD(6) - t10;
t213 = t167 * t222;
t120 = pkin(4) * t213 + pkin(8) * t225;
t131 = pkin(4) * t245 + t164 * pkin(8);
t231 = t164 ^ 2 - t167 ^ 2;
t228 = qJD(2) * t160;
t217 = t165 * t241;
t154 = -t161 * pkin(4) - pkin(3);
t215 = t165 * t228;
t214 = t168 * t228;
t206 = t164 * t219;
t25 = -t159 * t64 + t161 * t81;
t19 = qJD(2) * t180 + t25;
t22 = -pkin(9) * t193 + t26;
t205 = t163 * t22 - t166 * t19 + t34 * t223 + t27 * t224;
t202 = pkin(5) * t206;
t201 = t71 * t209;
t200 = t188 * t209;
t50 = pkin(4) * t193 + t67;
t199 = t164 * t214;
t198 = t167 * t214;
t197 = t164 * t209;
t194 = qJ(6) * t206;
t135 = -t209 - t251;
t192 = -t135 - t209;
t11 = t163 * t27 + t166 * t34;
t189 = -t163 * t89 + t166 * t78;
t117 = t162 * t164 + t167 * t243;
t85 = -t117 * t159 - t161 * t242;
t86 = t117 * t161 - t159 * t242;
t35 = t163 * t86 - t166 * t85;
t36 = t163 * t85 + t166 * t86;
t185 = -t11 * t151 - t205;
t116 = -t162 * t167 + t164 * t243;
t184 = -t163 * t19 - t166 * t22 - t27 * t223 + t34 * t224;
t178 = qJD(3) * t181;
t2 = -t202 + t205;
t40 = qJD(2) * t178 + t265;
t88 = t117 * qJD(3) + t199;
t87 = -t116 * qJD(3) + t198;
t56 = -t87 * t159 + t161 * t215;
t57 = t159 * t215 + t87 * t161;
t9 = t36 * qJD(5) + t163 * t57 - t166 * t56;
t177 = t116 * t40 + t9 * t151 - t35 * t206 + t88 * t71;
t176 = qJD(3) * (-t192 - t251);
t175 = -qJ(4) * t226 + (t204 - t92) * t167;
t3 = t40 * pkin(5) + t39 * qJ(6) - qJD(6) * t188 + t50;
t8 = -t35 * qJD(5) + t163 * t56 + t166 * t57;
t173 = t116 * t39 - t8 * t151 - t188 * t88 + t36 * t206;
t171 = t40 - t268;
t109 = t128 * t164;
t108 = t129 * t164;
t102 = -pkin(8) * t244 + t124;
t69 = t128 * pkin(5) - t129 * qJ(6) + t154;
t62 = t164 * t266 + t178;
t61 = t164 * t114 + t163 * t213 - t196;
t47 = t108 * pkin(5) + t109 * qJ(6) + t131;
t33 = pkin(5) * t188 + t71 * qJ(6);
t32 = t167 * pkin(5) - t189;
t31 = -t167 * qJ(6) + t254;
t16 = -t39 - t248;
t12 = t62 * pkin(5) + t61 * qJ(6) + t109 * qJD(6) + t120;
t7 = -t151 * qJ(6) + t11;
t5 = t151 * pkin(5) + t235;
t1 = -t151 * qJD(6) - t184 + t194;
t4 = [0, 0, -t217, -t168 * t241, 0, 0, 0, 0, 0, -t167 * t217 + (-t88 - t199) * qJD(3), t164 * t217 + (-t87 - t198) * qJD(3), -t88 * t125 + (-t167 * t56 + (t116 * t244 + t164 * t85) * qJD(3)) * qJD(2), t88 * t126 + (t167 * t57 + (t116 * t239 - t164 * t86) * qJD(3)) * qJD(2), t57 * t125 - t56 * t126 + (-t159 * t86 - t161 * t85) * t207, t67 * t116 + t25 * t85 + t26 * t86 + t42 * t56 + t43 * t57 + t92 * t88, 0, 0, 0, 0, 0, t177, -t173, t177, t188 * t9 - t35 * t39 - t36 * t40 - t8 * t71, t173, t1 * t36 + t3 * t116 + t15 * t88 + t2 * t35 + t5 * t9 + t7 * t8; 0, 0, 0, 0, 0.2e1 * t167 * t206, -0.2e1 * t231 * t219, t236, -t237, 0, -pkin(8) * t236 + t164 * t176, pkin(8) * t237 + t167 * t176 (t125 * t209 + t250 + (qJD(2) * t102 + t42) * qJD(3)) * t164 + (-t25 + (-pkin(8) * t125 + t159 * t92) * qJD(3) + (t144 - t253) * qJD(2)) * t167 (-t126 * t209 + t249 + (-qJD(2) * t103 - t43) * qJD(3)) * t164 + (t26 + (pkin(8) * t126 + t161 * t92) * qJD(3) + (t203 + t252) * qJD(2)) * t167 (-t159 * t26 - t161 * t25) * t164 - t253 * t126 + t252 * t125 + (-t159 * t43 - t161 * t42 + (-t102 * t161 - t103 * t159) * qJD(2)) * t225, -t92 * t197 + t25 * t102 + t26 * t103 + t252 * t43 + t253 * t42 + (t164 * t67 + t92 * t225) * pkin(8), t39 * t109 - t188 * t61, t39 * t108 + t109 * t40 - t188 * t62 + t61 * t71, t61 * t151 + t39 * t167 + (-qJD(2) * t109 + t188) * t226, t62 * t151 + t40 * t167 + (-qJD(2) * t108 - t71) * t226 (-t151 - t220) * t226, t205 * t167 + t120 * t71 + t131 * t40 + t50 * t108 + t60 * t62 + t270 * t151 + (-t201 + (qJD(2) * t189 + t10) * qJD(3)) * t164, -t184 * t167 + t120 * t188 - t131 * t39 - t50 * t109 - t60 * t61 + t272 * t151 + (-t200 + (-t254 * qJD(2) - t11) * qJD(3)) * t164, t3 * t108 + t12 * t71 + t15 * t62 + t2 * t167 + t47 * t40 + t262 * t151 + (-t201 + (-qJD(2) * t32 - t5) * qJD(3)) * t164, -t1 * t108 - t2 * t109 + t188 * t262 - t263 * t71 - t31 * t40 - t32 * t39 - t5 * t61 - t7 * t62, -t1 * t167 + t3 * t109 - t12 * t188 + t15 * t61 + t47 * t39 - t263 * t151 + (t200 + (qJD(2) * t31 + t7) * qJD(3)) * t164, t1 * t31 + t2 * t32 + t3 * t47 + t263 * t7 + t262 * t5 + (t12 - t197) * t15; 0, 0, 0, 0, -t164 * t170 * t167, t231 * t170, 0, 0, 0, t100 * qJD(3) - t135 * t227 - t67, t192 * t220, t100 * t125 - t249 + (t175 * t159 - t164 * t42 + t167 * t58) * qJD(2), -t100 * t126 + t250 + (t175 * t161 + t164 * t43 - t167 * t59) * qJD(2), -t59 * t125 + t58 * t126 + (qJD(4) * t125 + t42 * t220 + t26) * t161 + (qJD(4) * t126 + t43 * t220 - t25) * t159, -t67 * pkin(3) - t92 * t100 - t42 * t58 - t43 * t59 + (-t159 * t42 + t161 * t43) * qJD(4) + (-t25 * t159 + t26 * t161) * qJ(4), -t39 * t129 + t188 * t233, t39 * t128 - t129 * t40 - t188 * t232 - t233 * t71, -t233 * t151 + (qJD(3) * t129 - t188) * t227, t232 * t151 + (-qJD(3) * t128 + t71) * t227, t151 * t227, t50 * t128 + t154 * t40 - t80 * t71 + t232 * t60 + t271 * t151 + (-t10 + t247) * t227, t50 * t129 - t154 * t39 - t80 * t188 + t233 * t60 + t269 * t151 + (t11 - t246) * t227, t3 * t128 + t69 * t40 - t256 * t71 + t257 * t151 + t232 * t15 + (t5 + t247) * t227, -t1 * t128 + t2 * t129 + t187 * t39 + t188 * t257 - t232 * t7 + t233 * t5 + t258 * t71 - t91 * t40, -t3 * t129 + t69 * t39 + t256 * t188 + t258 * t151 - t233 * t15 + (-t7 + t246) * t227, t1 * t91 - t256 * t15 - t187 * t2 + t257 * t5 - t258 * t7 + t3 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t126 + t222) * t220 (-t125 + t221) * t220, -t125 ^ 2 - t126 ^ 2, -t43 * t125 + t42 * t126 + t67, 0, 0, 0, 0, 0, t171, -t274, t171, -t264 - t273, t274, -t188 * t5 + t7 * t71 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t260, t264 - t273, t16, -t129 * t207 - t265 - t268, t206, -t188 * t60 + t185, -t10 * t151 + t60 * t71 + t184, -t33 * t71 + t185 + 0.2e1 * t202 - t261, pkin(5) * t39 - t40 * qJ(6) + (-t11 + t7) * t188 + (t5 - t235) * t71, 0.2e1 * t194 - t15 * t71 + t33 * t188 + (-0.2e1 * qJD(6) + t10) * t151 - t184, -t2 * pkin(5) + t1 * qJ(6) - t5 * t11 - t15 * t33 + t235 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t206 + t260, t16, -t151 ^ 2 - t264, t7 * t151 + t2 + t261;];
tauc_reg  = t4;
