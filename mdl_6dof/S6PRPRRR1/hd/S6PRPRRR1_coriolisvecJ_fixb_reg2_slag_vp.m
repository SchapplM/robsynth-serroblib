% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:25:06
% EndTime: 2019-03-08 20:25:16
% DurationCPUTime: 3.55s
% Computational Cost: add. (6937->384), mult. (16822->532), div. (0->0), fcn. (13183->12), ass. (0->213)
t139 = sin(pkin(6));
t138 = sin(pkin(12));
t140 = cos(pkin(12));
t145 = sin(qJ(2));
t148 = cos(qJ(2));
t168 = t138 * t148 + t140 * t145;
t107 = t168 * t139;
t101 = qJD(1) * t107;
t144 = sin(qJ(4));
t209 = qJD(4) * t144;
t178 = pkin(4) * t209 - t101;
t212 = qJD(1) * t139;
t192 = t145 * t212;
t121 = t138 * t192;
t191 = t148 * t212;
t104 = t140 * t191 - t121;
t131 = t138 * pkin(2) + pkin(8);
t252 = pkin(9) + t131;
t186 = qJD(4) * t252;
t109 = t144 * t186;
t143 = sin(qJ(5));
t147 = cos(qJ(4));
t257 = cos(qJ(5));
t193 = t257 * t147;
t218 = t143 * t144;
t160 = t193 - t218;
t115 = t252 * t144;
t116 = t252 * t147;
t161 = -t257 * t115 - t143 * t116;
t179 = t147 * t186;
t249 = t161 * qJD(5) - t160 * t104 - t257 * t109 - t143 * t179;
t190 = t257 * qJD(5);
t203 = qJD(4) + qJD(5);
t87 = -qJD(4) * t193 - t147 * t190 + t203 * t218;
t118 = t143 * t147 + t257 * t144;
t88 = t203 * t118;
t266 = t88 * pkin(5) + t87 * pkin(10) + t178;
t142 = sin(qJ(6));
t146 = cos(qJ(6));
t205 = qJD(6) * t146;
t141 = cos(pkin(6));
t128 = t141 * qJD(1) + qJD(3);
t122 = t147 * t128;
t120 = qJD(2) * pkin(2) + t191;
t95 = t138 * t120 + t140 * t192;
t92 = qJD(2) * pkin(8) + t95;
t188 = pkin(9) * qJD(2) + t92;
t175 = t188 * t144;
t68 = t122 - t175;
t67 = qJD(4) * pkin(4) + t68;
t197 = t257 * t67;
t216 = t144 * t128;
t69 = t188 * t147 + t216;
t238 = t143 * t69;
t31 = t197 - t238;
t28 = -t203 * pkin(5) - t31;
t156 = t69 * qJD(4);
t106 = (t138 * t145 - t140 * t148) * t139;
t103 = qJD(2) * t106;
t100 = qJD(1) * t103;
t217 = t144 * t100;
t40 = t257 * (-t156 + t217);
t208 = qJD(4) * t147;
t225 = -t147 * t100 + t128 * t208;
t44 = -qJD(4) * t175 + t225;
t187 = t143 * t44 - t40;
t196 = t257 * t69;
t32 = t143 * t67 + t196;
t9 = t32 * qJD(5) + t187;
t265 = t9 * t142 + t28 * t205;
t79 = -t143 * t115 + t257 * t116;
t248 = t79 * qJD(5) - t118 * t104 - t143 * t109 + t257 * t179;
t264 = t87 * t203;
t206 = qJD(6) * t142;
t263 = t9 * t146 - t28 * t206;
t29 = t203 * pkin(10) + t32;
t211 = qJD(2) * t144;
t111 = -qJD(2) * t193 + t143 * t211;
t113 = t118 * qJD(2);
t194 = -t147 * pkin(4) - pkin(3);
t94 = t140 * t120 - t121;
t82 = t194 * qJD(2) - t94;
t51 = t111 * pkin(5) - t113 * pkin(10) + t82;
t17 = t142 * t51 + t146 * t29;
t170 = t142 * t29 - t146 * t51;
t262 = t142 * t170 + t146 * t17;
t260 = t203 * qJD(2);
t80 = t160 * t260;
t98 = t146 * t113 + t142 * t203;
t53 = t98 * qJD(6) + t142 * t80;
t110 = qJD(6) + t111;
t220 = t118 * t146;
t234 = t146 * t87;
t81 = t88 * qJD(2);
t261 = -t110 * (t118 * t206 + t234) + t81 * t220;
t204 = qJD(2) * qJD(4);
t189 = t144 * t204;
t102 = qJD(2) * t107;
t99 = qJD(1) * t102;
t86 = pkin(4) * t189 + t99;
t33 = t81 * pkin(5) - t80 * pkin(10) + t86;
t207 = qJD(5) * t143;
t159 = t143 * t217 - t69 * t207 + t257 * t44;
t8 = -t143 * t156 + t67 * t190 + t159;
t2 = -t170 * qJD(6) + t142 * t33 + t146 * t8;
t89 = -t107 * t144 + t141 * t147;
t90 = t107 * t147 + t141 * t144;
t165 = -t143 * t90 + t257 * t89;
t259 = t9 * t165;
t258 = t9 * t161;
t256 = t140 * pkin(2);
t1 = t2 * t146;
t255 = t9 * t160;
t184 = t146 * t203;
t96 = t142 * t113 - t184;
t253 = t98 * t96;
t123 = t194 - t256;
t77 = -pkin(5) * t160 - t118 * pkin(10) + t123;
t42 = -t142 * t79 + t146 * t77;
t251 = t42 * qJD(6) + t266 * t142 + t146 * t249;
t43 = t142 * t77 + t146 * t79;
t250 = -t43 * qJD(6) - t142 * t249 + t266 * t146;
t247 = -t53 * t220 + t96 * t234;
t52 = -qJD(6) * t184 + t113 * t206 - t146 * t80;
t246 = t160 * t52 + t98 * t88;
t245 = t87 * t111 - t118 * t81;
t244 = pkin(4) * qJD(5);
t241 = t142 * t81;
t240 = t142 * t87;
t239 = t142 * t96;
t237 = t144 * t92;
t235 = t146 * t81;
t233 = t146 * t98;
t232 = t52 * t142;
t231 = t53 * t146;
t230 = t81 * t160;
t229 = t82 * t113;
t228 = t96 * t110;
t227 = t98 * t110;
t226 = t99 * t106;
t224 = t110 * t113;
t223 = t111 * t142;
t222 = t111 * t146;
t221 = t113 * t111;
t150 = qJD(2) ^ 2;
t219 = t139 * t150;
t149 = qJD(4) ^ 2;
t215 = t149 * t144;
t214 = t149 * t147;
t136 = t144 ^ 2;
t137 = t147 ^ 2;
t213 = t136 - t137;
t210 = qJD(2) * t147;
t202 = -t17 * t223 + t170 * t222 + t1;
t201 = t98 * t240;
t199 = pkin(4) * t211;
t195 = t144 * t150 * t147;
t185 = t110 * t146;
t183 = pkin(4) * t190;
t182 = t17 * t113 + t265;
t181 = t147 * t189;
t34 = t143 * t68 + t196;
t180 = pkin(4) * t207 - t34;
t83 = t113 * pkin(5) + t111 * pkin(10);
t177 = t160 * t53 - t88 * t96;
t134 = t143 * pkin(4) + pkin(10);
t174 = t111 * t28 - t134 * t81;
t173 = -t113 * t88 + t160 * t80;
t172 = t142 * t17 - t146 * t170;
t73 = t122 - t237;
t74 = t147 * t92 + t216;
t169 = t144 * t73 - t147 * t74;
t55 = t143 * t89 + t257 * t90;
t36 = t106 * t146 - t142 * t55;
t37 = t106 * t142 + t146 * t55;
t166 = t113 * t170 - t263;
t163 = t101 * qJD(2) - t131 * t149 - t99;
t132 = -pkin(3) - t256;
t91 = -qJD(2) * pkin(3) - t94;
t162 = qJD(4) * (qJD(2) * t132 + t104 + t91);
t3 = -qJD(6) * t17 - t142 * t8 + t146 * t33;
t158 = -t172 * qJD(6) - t3 * t142;
t157 = -t232 + (t233 + t239) * qJD(6);
t155 = t143 * (-pkin(9) * t210 - t74);
t154 = t82 * t111 - t159;
t153 = (-t118 * t205 + t240) * t110 - t118 * t241;
t152 = t158 + t1;
t48 = -t92 * t209 + t225;
t49 = -t74 * qJD(4) + t217;
t151 = -t49 * t144 + t48 * t147 + (-t144 * t74 - t147 * t73) * qJD(4);
t135 = -t257 * pkin(4) - pkin(5);
t84 = t88 * t203;
t75 = t83 + t199;
t71 = -t111 ^ 2 + t113 ^ 2;
t64 = t113 * t203 - t118 * t260;
t63 = (qJD(2) * t160 + t111) * t203;
t59 = t89 * qJD(4) - t103 * t147;
t58 = -t90 * qJD(4) + t103 * t144;
t35 = t257 * t68 - t238;
t25 = t110 * t185 - t98 * t113 + t241;
t24 = -t110 ^ 2 * t142 + t96 * t113 + t235;
t23 = t142 * t228 - t231;
t22 = t98 * t185 - t232;
t21 = t142 * t83 + t146 * t31;
t20 = -t142 * t31 + t146 * t83;
t19 = t142 * t75 + t146 * t35;
t18 = -t142 * t35 + t146 * t75;
t15 = t55 * qJD(5) + t143 * t59 - t257 * t58;
t14 = t165 * qJD(5) + t143 * t58 + t257 * t59;
t6 = (-t52 - t228) * t146 + (-t53 - t227) * t142;
t5 = -t37 * qJD(6) + t102 * t146 - t142 * t14;
t4 = t36 * qJD(6) + t102 * t142 + t146 * t14;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145 * t219, -t148 * t219, 0, 0, 0, 0, 0, 0, 0, 0, -t102 * qJD(2), t103 * qJD(2), 0, -t100 * t107 - t94 * t102 - t95 * t103 + t226, 0, 0, 0, 0, 0, 0, t58 * qJD(4) + (-t102 * t147 + t106 * t209) * qJD(2), -t59 * qJD(4) + (t102 * t144 + t106 * t208) * qJD(2) (-t144 * t58 + t147 * t59 + (-t144 * t90 - t147 * t89) * qJD(4)) * qJD(2), t91 * t102 + t48 * t90 + t49 * t89 + t73 * t58 + t74 * t59 + t226, 0, 0, 0, 0, 0, 0, t102 * t111 + t106 * t81 - t15 * t203, t102 * t113 + t106 * t80 - t14 * t203, -t14 * t111 + t15 * t113 - t165 * t80 - t55 * t81, t82 * t102 + t86 * t106 + t32 * t14 - t31 * t15 + t8 * t55 - t259, 0, 0, 0, 0, 0, 0, t5 * t110 + t15 * t96 - t165 * t53 + t36 * t81, -t4 * t110 + t15 * t98 + t165 * t52 - t37 * t81, t36 * t52 - t37 * t53 - t4 * t96 - t5 * t98, t15 * t28 + t17 * t4 - t170 * t5 + t2 * t37 + t3 * t36 - t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t168 * t212 + t101) * qJD(2) (qJD(1) * t106 + t104) * qJD(2), 0, t94 * t101 - t95 * t104 + (-t100 * t138 - t140 * t99) * pkin(2), 0.2e1 * t181, -0.2e1 * t213 * t204, t214, -0.2e1 * t181, -t215, 0, t144 * t162 + t163 * t147, -t163 * t144 + t147 * t162 (-t136 - t137) * t104 * qJD(2) + t151, -t91 * t101 + t169 * t104 + t151 * t131 + t99 * t132, -t113 * t87 + t80 * t118, t173 + t245, -t264, t111 * t88 - t230, -t84, 0, t178 * t111 + t123 * t81 - t160 * t86 - t203 * t248 + t82 * t88, t178 * t113 + t86 * t118 + t123 * t80 - t203 * t249 - t82 * t87, -t249 * t111 + t248 * t113 + t9 * t118 + t160 * t8 - t161 * t80 + t31 * t87 - t32 * t88 - t79 * t81, t86 * t123 + t178 * t82 - t248 * t31 + t249 * t32 + t8 * t79 - t258, -t87 * t233 + (-t52 * t146 - t98 * t206) * t118, t201 + (t232 + (-t233 + t239) * qJD(6)) * t118 + t247, t246 + t261, -t87 * t239 + (t142 * t53 + t96 * t205) * t118, t153 + t177, t110 * t88 - t230, t250 * t110 + t118 * t265 - t160 * t3 - t161 * t53 - t170 * t88 - t28 * t240 + t248 * t96 + t42 * t81, -t251 * t110 + t118 * t263 + t160 * t2 + t161 * t52 - t17 * t88 - t28 * t234 + t248 * t98 - t43 * t81, t42 * t52 - t43 * t53 - t250 * t98 - t251 * t96 + t172 * t87 + (-qJD(6) * t262 - t142 * t2 - t146 * t3) * t118, t17 * t251 - t170 * t250 + t2 * t43 + t248 * t28 + t3 * t42 - t258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t215, -t214, 0, -t169 * qJD(4) + t48 * t144 + t49 * t147, 0, 0, 0, 0, 0, 0, -t84, t264, -t173 + t245, t8 * t118 - t31 * t88 - t32 * t87 - t255, 0, 0, 0, 0, 0, 0, t153 - t177, t246 - t261, t118 * t157 - t201 + t247, t118 * t152 - t262 * t87 + t28 * t88 - t255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t195, t213 * t150, 0, t195, 0, 0 (-qJD(2) * t91 + t100) * t144, -t91 * t210 + (t73 + t237) * qJD(4) - t225, 0, 0, t221, t71, t63, -t221, t64, 0, -t69 * t190 + t40 + t34 * t203 - t111 * t199 - t229 + (-qJD(5) * t67 - t203 * t244 - t44) * t143 (-t197 + t35) * qJD(5) + (-t155 + t35) * qJD(4) + (-t113 * t211 - t203 * t190) * pkin(4) + t154 (t32 - t34) * t113 + (-t31 + t35) * t111 + (-t257 * t80 - t143 * t81 + (-t257 * t111 + t113 * t143) * qJD(5)) * pkin(4), t31 * t34 - t32 * t35 + (-t82 * t211 - t257 * t9 + t143 * t8 + (-t143 * t31 + t257 * t32) * qJD(5)) * pkin(4), t22, t6, t25, t23, t24, -t224, t135 * t53 + t180 * t96 + t174 * t142 + (-t134 * t205 - t142 * t183 - t18) * t110 + t166, -t135 * t52 + t180 * t98 + t174 * t146 + (t134 * t206 - t146 * t183 + t19) * t110 + t182, t18 * t98 + t19 * t96 + (-t96 * t183 - t134 * t53 + (t134 * t98 + t170) * qJD(6)) * t146 + (t98 * t183 - t134 * t52 - t3 + (t134 * t96 - t17) * qJD(6)) * t142 + t202, t9 * t135 + t170 * t18 - t17 * t19 - t28 * t34 + (t143 * t28 + t257 * t262) * t244 + t152 * t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, t71, t63, -t221, t64, 0, t32 * qJD(4) - t187 - t229 (-t197 + t31) * qJD(5) + (-t155 + t31) * qJD(4) + t154, 0, 0, t22, t6, t25, t23, t24, -t224, t28 * t223 - pkin(5) * t53 - t20 * t110 - t32 * t96 + (-t110 * t205 - t241) * pkin(10) + t166, t28 * t222 + pkin(5) * t52 + t21 * t110 - t32 * t98 + (t110 * t206 - t235) * pkin(10) + t182, t20 * t98 + t21 * t96 + (t157 - t231) * pkin(10) + t158 + t202, -t9 * pkin(5) + pkin(10) * t152 - t17 * t21 + t170 * t20 - t28 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, -t96 ^ 2 + t98 ^ 2, -t52 + t228, -t253, t227 - t53, t81, t17 * t110 - t28 * t98 + t3, -t110 * t170 + t28 * t96 - t2, 0, 0;];
tauc_reg  = t7;
