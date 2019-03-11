% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:52:29
% EndTime: 2019-03-09 09:52:37
% DurationCPUTime: 2.91s
% Computational Cost: add. (5269->401), mult. (13216->488), div. (0->0), fcn. (9289->6), ass. (0->200)
t157 = cos(qJ(2));
t228 = cos(pkin(9));
t194 = t228 * t157;
t144 = qJD(1) * t194;
t153 = sin(pkin(9));
t155 = sin(qJ(2));
t212 = qJD(1) * t155;
t121 = -t153 * t212 + t144;
t117 = qJD(4) - t121;
t136 = t153 * t157 + t155 * t228;
t123 = t136 * qJD(1);
t154 = sin(qJ(4));
t156 = cos(qJ(4));
t91 = qJD(2) * t154 + t123 * t156;
t231 = t91 * t117;
t206 = qJD(1) * qJD(2);
t198 = t155 * t206;
t169 = qJD(2) * t144 - t153 * t198;
t51 = qJD(4) * t91 + t154 * t169;
t271 = t51 - t231;
t122 = t136 * qJD(2);
t111 = qJD(1) * t122;
t270 = qJ(5) * t111 + qJD(5) * t117;
t207 = t156 * qJD(2);
t211 = qJD(4) * t154;
t180 = qJD(4) * t207 - t123 * t211 + t156 * t169;
t89 = t123 * t154 - t207;
t240 = t117 * t89;
t27 = t180 + t240;
t269 = -0.2e1 * t206;
t268 = qJ(6) * t51 + qJD(6) * t89;
t253 = -qJ(3) - pkin(7);
t142 = t253 * t157;
t139 = qJD(1) * t142;
t126 = t153 * t139;
t141 = t253 * t155;
t138 = qJD(1) * t141;
t245 = qJD(2) * pkin(2);
t130 = t138 + t245;
t76 = t130 * t228 + t126;
t190 = qJD(2) * pkin(3) + t76;
t172 = t91 * qJ(5) + t190;
t258 = pkin(4) + pkin(5);
t17 = -t258 * t89 + qJD(6) + t172;
t267 = (qJD(6) + t17) * t91;
t114 = t117 ^ 2;
t88 = t91 ^ 2;
t266 = -t114 - t88;
t195 = t228 * t139;
t78 = t153 * t138 - t195;
t265 = -t154 * qJD(5) - t78;
t264 = t258 * t156;
t202 = -pkin(2) * t157 - pkin(1);
t189 = t202 * qJD(1);
t140 = qJD(3) + t189;
t55 = -pkin(3) * t121 - pkin(8) * t123 + t140;
t77 = t130 * t153 - t195;
t71 = qJD(2) * pkin(8) + t77;
t28 = -t154 * t71 + t156 * t55;
t215 = qJD(5) - t28;
t263 = 0.2e1 * t270;
t175 = -t153 * t155 + t194;
t262 = qJ(5) * t122 - qJD(5) * t175;
t100 = t156 * t111;
t225 = t121 * t154;
t181 = t100 + (-t211 + t225) * t117;
t239 = t123 * t89;
t163 = t181 - t239;
t18 = t91 * qJ(6) + t28;
t216 = qJD(5) - t18;
t12 = -t117 * t258 + t216;
t210 = qJD(4) * t156;
t196 = qJD(2) * t253;
t118 = t157 * qJD(3) + t155 * t196;
t103 = t118 * qJD(1);
t119 = -t155 * qJD(3) + t157 * t196;
t164 = qJD(1) * t119;
t49 = t103 * t228 + t153 * t164;
t145 = pkin(2) * t198;
t54 = pkin(3) * t111 - pkin(8) * t169 + t145;
t177 = -t154 * t54 - t156 * t49 - t210 * t55 + t211 * t71;
t5 = -t177 + t270;
t2 = t5 + t268;
t261 = t117 * t12 + t2;
t147 = pkin(2) * t153 + pkin(8);
t222 = t147 * t111;
t30 = pkin(4) * t89 - t172;
t260 = t117 * t30 - t222;
t259 = t89 ^ 2;
t106 = t111 * pkin(4);
t257 = t156 * pkin(4);
t48 = t103 * t153 - t164 * t228;
t10 = pkin(4) * t51 - qJ(5) * t180 - qJD(5) * t91 + t48;
t6 = -pkin(5) * t51 - t10;
t256 = t6 * t154;
t255 = t6 * t156;
t254 = t91 * t89;
t217 = qJ(6) - t147;
t134 = t217 * t156;
t65 = pkin(2) * t212 + pkin(3) * t123 - pkin(8) * t121;
t79 = t138 * t228 + t126;
t73 = t154 * t79;
t252 = t73 + (-qJ(6) * t121 - t65) * t156 - t258 * t123 + qJD(4) * t134 + t154 * qJD(6);
t248 = t154 * t65 + t156 * t79;
t24 = qJ(5) * t123 + t248;
t251 = -qJ(6) * t225 - t156 * qJD(6) + t211 * t217 - t24;
t227 = qJ(5) * t156;
t178 = -t154 * t258 + t227;
t250 = -t117 * t178 + t265;
t249 = -t154 * t51 - t210 * t89;
t29 = t154 * t55 + t156 * t71;
t75 = -pkin(3) * t175 - pkin(8) * t136 + t202;
t85 = t141 * t153 - t142 * t228;
t247 = t154 * t75 + t156 * t85;
t246 = qJ(5) * t51;
t244 = t10 * t154;
t243 = t10 * t156;
t109 = t117 * qJ(5);
t22 = t109 + t29;
t242 = t117 * t22;
t241 = t117 * t29;
t19 = qJ(6) * t89 + t29;
t14 = t109 + t19;
t238 = t14 * t117;
t237 = t154 * t91;
t236 = t156 * t89;
t235 = t48 * t154;
t234 = t48 * t156;
t233 = t180 * t154;
t232 = t89 * qJ(5);
t230 = t91 * t123;
t187 = pkin(4) * t154 - t227;
t229 = t117 * t187 + t265;
t226 = qJ(6) * t136;
t191 = t117 * t156;
t125 = t175 * qJD(2);
t224 = t125 * t154;
t223 = t125 * t156;
t221 = t154 * qJ(5);
t99 = t154 * t111;
t160 = qJD(1) ^ 2;
t220 = t157 * t160;
t159 = qJD(2) ^ 2;
t219 = t159 * t155;
t218 = t159 * t157;
t213 = t155 ^ 2 - t157 ^ 2;
t209 = qJD(5) * t156;
t62 = t118 * t228 + t119 * t153;
t205 = t154 * t62 + t210 * t85 + t211 * t75;
t203 = t155 * t245;
t66 = pkin(3) * t122 - pkin(8) * t125 + t203;
t204 = t154 * t66 + t156 * t62 + t210 * t75;
t31 = -qJ(5) * t175 + t247;
t201 = t136 * t211;
t200 = t136 * t210;
t199 = t147 * t211;
t81 = t154 * t85;
t197 = t156 * t75 - t81;
t193 = t154 * t49 - t156 * t54 + t210 * t71 + t211 * t55;
t192 = pkin(1) * t269;
t61 = t153 * t118 - t119 * t228;
t84 = -t141 * t228 - t153 * t142;
t148 = -pkin(2) * t228 - pkin(3);
t7 = -t106 + t193;
t171 = -qJ(6) * t180 + t7;
t165 = -t111 * pkin(5) + t171;
t1 = -t91 * qJD(6) + t165;
t188 = -t1 + t238;
t21 = -pkin(4) * t117 + t215;
t186 = -t154 * t22 + t156 * t21;
t185 = t236 + t237;
t184 = t156 * t66 - t205;
t183 = -qJ(6) * t125 - qJD(6) * t136;
t182 = t117 * t210 - t121 * t191 + t99;
t179 = t156 * t180 - t211 * t91;
t176 = -t211 * t85 + t204;
t174 = -t117 * t190 - t222;
t173 = t30 * t91 + t7;
t170 = t148 - t221;
t168 = t182 + t230;
t167 = t117 * t28 + t177;
t162 = -qJD(2) * t123 + t254;
t133 = t217 * t154;
t132 = t170 - t257;
t110 = -t170 + t264;
t37 = pkin(4) * t91 + t232;
t36 = t136 * t187 + t84;
t34 = t136 * t178 - t84;
t33 = -t258 * t91 - t232;
t32 = pkin(4) * t175 - t197;
t25 = -t123 * pkin(4) - t156 * t65 + t73;
t23 = t154 * t226 + t31;
t20 = t81 + (-t75 - t226) * t156 + t258 * t175;
t13 = t187 * t125 + (-t209 + (t221 + t257) * qJD(4)) * t136 + t61;
t11 = t178 * t125 + (t209 + (-t221 - t264) * qJD(4)) * t136 - t61;
t9 = -t122 * pkin(4) - t184;
t8 = t176 + t262;
t4 = qJ(6) * t200 + (-qJD(4) * t85 - t183) * t154 + t204 + t262;
t3 = qJ(6) * t201 - t258 * t122 + (t183 - t66) * t156 + t205;
t15 = [0, 0, 0, 0.2e1 * t157 * t198, t213 * t269, t218, -t219, 0, -pkin(7) * t218 + t155 * t192, pkin(7) * t219 + t157 * t192, -t111 * t85 + t121 * t62 - t122 * t77 + t123 * t61 - t125 * t76 + t136 * t48 + t169 * t84 + t175 * t49, t48 * t84 + t49 * t85 - t76 * t61 + t77 * t62 + (t140 + t189) * t203, t136 * t179 + t223 * t91, -t185 * t125 + (-t233 - t156 * t51 + (t154 * t89 - t156 * t91) * qJD(4)) * t136, t136 * t100 + t91 * t122 - t180 * t175 + (-t201 + t223) * t117, -t136 * t99 - t89 * t122 + t51 * t175 + (-t200 - t224) * t117, -t111 * t175 + t117 * t122, t184 * t117 + t197 * t111 + t193 * t175 + t28 * t122 + t61 * t89 + t84 * t51 - t190 * t224 + (-t190 * t210 + t235) * t136, -t176 * t117 - t247 * t111 - t177 * t175 - t29 * t122 + t61 * t91 + t84 * t180 - t190 * t223 + (t190 * t211 + t234) * t136, t30 * t224 - t32 * t111 - t9 * t117 - t21 * t122 + t13 * t89 + t7 * t175 + t36 * t51 + (t210 * t30 + t244) * t136, -t31 * t51 + t32 * t180 - t8 * t89 + t9 * t91 + t186 * t125 + (-t5 * t154 + t7 * t156 + (-t154 * t21 - t156 * t22) * qJD(4)) * t136, -t30 * t223 + t31 * t111 + t8 * t117 + t22 * t122 - t13 * t91 - t5 * t175 - t36 * t180 + (t211 * t30 - t243) * t136, t10 * t36 + t13 * t30 + t21 * t9 + t22 * t8 + t31 * t5 + t32 * t7, -t17 * t224 + t1 * t175 - t11 * t89 - t20 * t111 - t3 * t117 - t12 * t122 - t34 * t51 + (-t17 * t210 - t256) * t136, t17 * t223 + t11 * t91 + t23 * t111 + t4 * t117 + t14 * t122 - t2 * t175 + t34 * t180 + (-t17 * t211 + t255) * t136, -t20 * t180 + t23 * t51 - t3 * t91 + t4 * t89 + (-t12 * t156 + t14 * t154) * t125 + (-t1 * t156 + t2 * t154 + (t12 * t154 + t14 * t156) * qJD(4)) * t136, t1 * t20 + t11 * t17 + t12 * t3 + t14 * t4 + t2 * t23 + t34 * t6; 0, 0, 0, -t155 * t220, t213 * t160, 0, 0, 0, t160 * pkin(1) * t155, pkin(1) * t220 (t77 - t78) * t123 + (-t79 + t76) * t121 + (-t111 * t153 - t169 * t228) * pkin(2), t76 * t78 - t77 * t79 + (-t140 * t212 + t153 * t49 - t228 * t48) * pkin(2), t191 * t91 + t233, t121 * t185 + t179 + t249, t182 - t230, t181 + t239, -t117 * t123, -t28 * t123 + t148 * t51 - t234 - t78 * t89 + (t73 + (-qJD(4) * t147 - t65) * t156) * t117 + t174 * t154, t29 * t123 + t148 * t180 + t235 - t78 * t91 + (t199 + t248) * t117 + t174 * t156, -t243 + t21 * t123 + t132 * t51 + t229 * t89 + (-t147 * t210 + t25) * t117 + t260 * t154, t24 * t89 - t25 * t91 + (-t121 * t21 - t147 * t51 + t5 + (t147 * t91 + t21) * qJD(4)) * t156 + (t121 * t22 + t147 * t180 + t7 + (t147 * t89 - t22) * qJD(4)) * t154, -t244 - t22 * t123 - t132 * t180 - t229 * t91 + (-t24 - t199) * t117 - t260 * t156, t10 * t132 - t21 * t25 - t22 * t24 + t229 * t30 + (qJD(4) * t186 + t7 * t154 + t5 * t156) * t147, -t110 * t51 + t133 * t111 + t12 * t123 + t250 * t89 + t255 + (-t154 * t17 + t252) * t117, t110 * t180 - t134 * t111 + t117 * t251 - t14 * t123 + t17 * t191 - t250 * t91 + t256, t133 * t180 - t134 * t51 + t154 * t188 - t156 * t261 + t251 * t89 + t252 * t91, -t1 * t133 + t6 * t110 - t12 * t252 - t2 * t134 + t14 * t251 - t17 * t250; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121 ^ 2 - t123 ^ 2, -t121 * t77 + t123 * t76 + t145, 0, 0, 0, 0, 0, t163, -t114 * t156 - t230 - t99, t163 (t236 - t237) * t121 - t179 + t249, t168, -t30 * t123 + (-t7 + t242) * t156 + (t117 * t21 + t5) * t154, t163, t168, t154 * t271 + t156 * t27, t17 * t123 + t154 * t261 + t156 * t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t254, t88 - t259, t27, -t271, t111, t190 * t91 - t193 + t241, -t190 * t89 + t167, -t37 * t89 + t106 - t173 + t241, -pkin(4) * t180 - t246 + (t22 - t29) * t91 + (t21 - t215) * t89, -t30 * t89 + t37 * t91 - t167 + t263, -t7 * pkin(4) + t5 * qJ(5) - t21 * t29 + t215 * t22 - t30 * t37, t19 * t117 + t33 * t89 + t267 + (pkin(5) + t258) * t111 - t171, -t117 * t18 + t17 * t89 - t33 * t91 - t177 + t263 + t268, t246 + t258 * t180 + (-t14 + t19) * t91 + (-t12 + t216) * t89, t2 * qJ(5) - t1 * t258 - t12 * t19 + t14 * t216 - t17 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t162, t27, t266, t173 - t242, t162, t266, -t27, t165 - t238 - t267; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51 - t231, t180 - t240, -t88 - t259, t12 * t91 - t14 * t89 + t6;];
tauc_reg  = t15;
