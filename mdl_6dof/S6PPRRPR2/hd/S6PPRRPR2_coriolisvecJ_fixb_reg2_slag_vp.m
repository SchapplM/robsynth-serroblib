% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PPRRPR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:50:55
% EndTime: 2019-03-08 18:51:01
% DurationCPUTime: 2.94s
% Computational Cost: add. (3654->332), mult. (9701->474), div. (0->0), fcn. (8128->12), ass. (0->190)
t110 = sin(qJ(4));
t179 = t110 * qJD(5);
t186 = qJD(4) * t110;
t105 = sin(pkin(6));
t103 = sin(pkin(12));
t111 = sin(qJ(3));
t114 = cos(qJ(3));
t106 = cos(pkin(12));
t107 = cos(pkin(7));
t199 = t106 * t107;
t129 = t103 * t114 + t111 * t199;
t124 = t129 * t105;
t104 = sin(pkin(7));
t202 = t104 * t111;
t108 = cos(pkin(6));
t94 = t108 * qJD(1) + qJD(2);
t175 = t94 * t202;
t47 = qJD(1) * t124 + t175;
t241 = -pkin(4) * t186 + t179 + t47;
t113 = cos(qJ(4));
t145 = pkin(10) * t110 - qJ(5) * t113;
t128 = t145 * qJD(4);
t240 = -t128 + t241;
t45 = qJD(3) * pkin(9) + t47;
t190 = qJD(1) * t105;
t168 = t106 * t190;
t64 = -t104 * t168 + t107 * t94;
t60 = t113 * t64;
t31 = t110 * t45 - t60;
t239 = -qJD(5) - t31;
t109 = sin(qJ(6));
t112 = cos(qJ(6));
t178 = t112 * qJD(4);
t188 = qJD(3) * t113;
t75 = -t109 * t188 + t178;
t201 = t104 * t114;
t203 = t103 * t111;
t232 = (t114 * t199 - t203) * t105;
t238 = t108 * t201 + t232;
t237 = qJD(1) * t232 + t94 * t201;
t115 = -pkin(4) - pkin(10);
t157 = pkin(5) * qJD(3) + t45;
t192 = t157 * t110 + qJD(5) - t60;
t21 = t115 * qJD(4) + t192;
t46 = (t104 * t94 + t107 * t168) * t114 - t190 * t203;
t158 = -t110 * qJ(5) - pkin(3);
t72 = t115 * t113 + t158;
t35 = t72 * qJD(3) - t46;
t144 = t109 * t35 - t112 * t21;
t118 = t47 - t179;
t177 = qJD(3) * qJD(4);
t160 = t110 * t177;
t95 = pkin(4) * t160;
t33 = t95 + (t128 + t118) * qJD(3);
t211 = t110 * t64;
t42 = t237 * qJD(3);
t212 = t110 * t42;
t9 = t212 + (t157 * t113 + t211) * qJD(4);
t1 = -t144 * qJD(6) + t109 * t9 + t112 * t33;
t180 = t110 * qJD(3);
t96 = qJD(6) + t180;
t236 = t144 * t96 + t1;
t6 = t109 * t21 + t112 * t35;
t2 = -qJD(6) * t6 - t109 * t33 + t112 * t9;
t235 = t6 * t96 + t2;
t25 = -qJD(4) * pkin(4) - t239;
t32 = t113 * t45 + t211;
t28 = -qJD(4) * qJ(5) - t32;
t181 = t109 * qJD(4);
t73 = t112 * t188 + t181;
t142 = t73 * t96;
t61 = t73 * qJD(6) - t109 * t160;
t234 = t61 - t142;
t217 = t75 * t96;
t62 = t75 * qJD(6) - t112 * t160;
t233 = -t62 + t217;
t52 = t108 * t202 + t124;
t68 = -t105 * t106 * t104 + t108 * t107;
t38 = t68 * t110 + t52 * t113;
t49 = t238 * qJD(3);
t15 = t38 * qJD(4) + t49 * t110;
t50 = t52 * qJD(3);
t231 = (-t113 * t50 - t186 * t238) * qJD(3) - t15 * qJD(4);
t16 = -t52 * t186 + (qJD(4) * t68 + t49) * t113;
t185 = qJD(4) * t113;
t230 = qJD(3) * (t110 * t50 - t185 * t238) - t16 * qJD(4);
t116 = qJD(4) ^ 2;
t138 = pkin(9) * t116;
t165 = qJ(5) * t185;
t214 = t165 + t241;
t34 = t95 + (t118 - t165) * qJD(3);
t229 = qJD(3) * t214 - t138 - t34;
t117 = qJD(3) ^ 2;
t166 = qJD(3) * t201;
t151 = t113 * t166;
t170 = t110 * t202;
t55 = qJD(4) * t170 - t107 * t185 - t151;
t228 = qJD(4) * (-t55 + t151) - t117 * t170;
t150 = t110 * t166;
t169 = t113 * t202;
t70 = t110 * t107 + t169;
t56 = t70 * qJD(4) + t150;
t227 = (t56 + t150) * qJD(4) + t117 * t169;
t226 = pkin(5) + pkin(9);
t12 = t32 * qJD(4) + t212;
t37 = t52 * t110 - t68 * t113;
t223 = t12 * t37;
t69 = -t113 * t107 + t170;
t222 = t12 * t69;
t43 = t47 * qJD(3);
t221 = t43 * t238;
t176 = -t113 * t42 - t64 * t185 + t45 * t186;
t7 = (-pkin(5) * t180 + qJD(5)) * qJD(4) - t176;
t220 = t7 * t109;
t219 = t7 * t112;
t218 = t75 * t73;
t197 = t110 * t112;
t85 = t226 * t110;
t54 = t109 * t85 + t112 * t72;
t86 = t226 * t113;
t78 = qJD(4) * t86;
t216 = t54 * qJD(6) - t240 * t109 - t112 * t78 + t46 * t197;
t198 = t109 * t110;
t53 = -t109 * t72 + t112 * t85;
t215 = -t53 * qJD(6) - t109 * t78 + t240 * t112 + t46 * t198;
t213 = qJD(3) * pkin(3);
t210 = t112 * t61;
t209 = t113 * t75;
t208 = t114 * t43;
t207 = t115 * t96;
t206 = t12 * t110;
t205 = t62 * t109;
t80 = -t113 * pkin(4) + t158;
t204 = qJD(3) * t80;
t200 = t104 * t117;
t196 = t116 * t110;
t195 = t116 * t113;
t101 = t110 ^ 2;
t102 = t113 ^ 2;
t191 = t101 - t102;
t189 = qJD(3) * t111;
t184 = qJD(6) * t109;
t183 = qJD(6) * t112;
t182 = qJD(6) * t113;
t173 = t96 * t197;
t172 = t96 * t183;
t167 = t104 * t189;
t164 = t109 * t182;
t163 = t112 * t182;
t159 = t113 * t177;
t39 = -t46 + t204;
t156 = t39 * t180 + t212;
t149 = t110 * t159;
t148 = (-t101 - t102) * t46 * qJD(3);
t147 = t109 * t144 + t112 * t6;
t17 = t109 * t238 + t37 * t112;
t18 = t37 * t109 - t112 * t238;
t143 = -qJD(3) * t102 + t110 * t96;
t141 = t109 * t96;
t44 = -t46 - t213;
t136 = qJD(4) * (t44 + t46 - t213);
t135 = qJD(4) * (-t39 - t46 - t204);
t132 = -t109 * t69 + t112 * t201;
t57 = t109 * t201 + t112 * t69;
t98 = pkin(5) * t188;
t22 = -t28 + t98;
t131 = t110 * t22 + t115 * t185;
t10 = -qJD(4) * qJD(5) + t176;
t122 = -t10 * t113 + t206 + (t110 * t28 + t113 * t25) * qJD(4);
t121 = -t176 * t113 + t206 + (-t110 * t32 + t113 * t31) * qJD(4);
t120 = (t110 * t15 + t113 * t16 + (-t110 * t38 + t113 * t37) * qJD(4)) * qJD(3);
t119 = (t110 * t56 - t113 * t55 + (-t110 * t70 + t113 * t69) * qJD(4)) * qJD(3);
t100 = pkin(4) * t180;
t93 = t110 * t117 * t113;
t89 = t112 * t159;
t84 = -0.2e1 * t149;
t83 = 0.2e1 * t149;
t82 = t191 * t117;
t77 = t226 * t186;
t76 = -qJ(5) * t188 + t100;
t71 = -0.2e1 * t191 * t177;
t66 = t145 * qJD(3) + t100;
t30 = t132 * qJD(6) - t109 * t167 + t112 * t56;
t29 = t57 * qJD(6) + t109 * t56 + t112 * t167;
t24 = t32 + t98;
t14 = t109 * t24 + t112 * t66;
t13 = -t109 * t66 + t112 * t24;
t4 = t17 * qJD(6) + t15 * t109 + t50 * t112;
t3 = -t18 * qJD(6) - t50 * t109 + t15 * t112;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50 * qJD(3), -t49 * qJD(3), 0, t42 * t52 - t46 * t50 + t47 * t49 - t221, 0, 0, 0, 0, 0, 0, t231, t230, t120, t31 * t15 + t32 * t16 - t176 * t38 + t44 * t50 - t221 + t223, 0, 0, 0, 0, 0, 0, t120, -t231, -t230, -t10 * t38 + t25 * t15 - t28 * t16 - t238 * t34 + t39 * t50 + t223, 0, 0, 0, 0, 0, 0, t159 * t17 + t16 * t73 + t3 * t96 + t38 * t62, -t159 * t18 + t16 * t75 - t38 * t61 - t4 * t96, t17 * t61 - t18 * t62 - t3 * t75 - t4 * t73, t1 * t18 - t144 * t3 + t22 * t16 + t2 * t17 + t7 * t38 + t6 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111 * t200, -t114 * t200, 0 (t111 * t42 - t208 + (-t111 * t46 + t114 * t47) * qJD(3)) * t104, 0, 0, 0, 0, 0, 0, -t227, -t228, t119, -t176 * t70 + t222 + t31 * t56 - t32 * t55 + (t189 * t44 - t208) * t104, 0, 0, 0, 0, 0, 0, t119, t227, t228, -t10 * t70 + t222 + t25 * t56 + t28 * t55 + (-t114 * t34 + t189 * t39) * t104, 0, 0, 0, 0, 0, 0, t159 * t57 + t30 * t96 - t55 * t73 + t70 * t62, t132 * t159 - t29 * t96 - t55 * t75 - t70 * t61, t132 * t62 - t29 * t73 - t30 * t75 + t57 * t61, -t1 * t132 - t144 * t30 + t2 * t57 - t22 * t55 + t6 * t29 + t7 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t129 * t190 - t175 + t47) * qJD(3) (-t237 + t46) * qJD(3), 0, 0, t83, t71, t195, t84, -t196, 0, t110 * t136 - t138 * t113, t110 * t138 + t113 * t136, t148 + t121, -t43 * pkin(3) - t44 * t47 + (-t110 * t31 - t113 * t32) * t46 + t121 * pkin(9), 0, -t195, t196, t83, t71, t84, t148 + t122, t110 * t135 - t113 * t229, t110 * t229 + t113 * t135, t34 * t80 + (-t110 * t25 + t113 * t28) * t46 - t214 * t39 + t122 * pkin(9), -t75 * t163 + (t113 * t61 + t186 * t75) * t109 (-t109 * t73 + t112 * t75) * t186 + (t205 + t210 + (t109 * t75 + t112 * t73) * qJD(6)) * t113, -t96 * t163 - t61 * t110 + (t109 * t143 + t209) * qJD(4), -t73 * t164 + (t113 * t62 - t186 * t73) * t112, t96 * t164 - t62 * t110 + (t112 * t143 - t113 * t73) * qJD(4) (t96 + t180) * t185, t86 * t62 - t77 * t73 - t216 * t96 + (-t178 * t22 + t2) * t110 + (-t22 * t184 + t219 - t46 * t73 + (qJD(3) * t53 - t144) * qJD(4)) * t113, -t86 * t61 - t77 * t75 + t215 * t96 + (t181 * t22 - t1) * t110 + (-t22 * t183 - t220 - t46 * t75 + (-qJD(3) * t54 - t6) * qJD(4)) * t113, t53 * t61 - t54 * t62 + t216 * t75 + t215 * t73 + t147 * t186 + (-t1 * t112 + t109 * t2 + (t109 * t6 - t112 * t144) * qJD(6)) * t113, t1 * t54 + t2 * t53 + t7 * t86 - t215 * t6 + t216 * t144 + (-t113 * t46 - t77) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, t82, 0, t93, 0, 0 (-qJD(3) * t44 - t42) * t110, -t31 * qJD(4) - t188 * t44 + t176, 0, 0, 0, 0, 0, -t93, t82, t93, 0, -t188 * t76 + t156 (0.2e1 * qJD(5) + t31) * qJD(4) + (t110 * t76 + t113 * t39) * qJD(3) - t176, -t12 * pkin(4) - t10 * qJ(5) + t239 * t28 - t25 * t32 - t39 * t76, -t141 * t75 - t210 (-t62 - t217) * t112 + (t61 + t142) * t109, -t96 * t184 + t89 + (-t198 * t96 - t209) * qJD(3), t112 * t142 + t205, -t172 + (-t173 + (t73 - t181) * t113) * qJD(3), -t96 * t188, qJ(5) * t62 + t220 - t13 * t96 + t192 * t73 + (-t109 * t207 + t112 * t22) * qJD(6) + (t112 * t131 + t113 * t144) * qJD(3), -qJ(5) * t61 + t219 + t14 * t96 + t192 * t75 + (-t109 * t22 - t112 * t207) * qJD(6) + (-t109 * t131 + t113 * t6) * qJD(3), t13 * t75 + t14 * t73 + (-t6 * t180 + t115 * t61 - t2 + (-t115 * t73 - t6) * qJD(6)) * t112 + (-t144 * t180 - t115 * t62 - t1 + (t115 * t75 - t144) * qJD(6)) * t109, t7 * qJ(5) + t144 * t13 - t6 * t14 + t192 * t22 + (qJD(6) * t147 + t1 * t109 + t2 * t112) * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, -t101 * t117 - t116 (t28 + t32) * qJD(4) + t156, 0, 0, 0, 0, 0, 0, -qJD(4) * t73 - t141 * t96 + t89, -t172 - qJD(4) * t75 + (-t113 * t181 - t173) * qJD(3), t233 * t109 + t234 * t112, -t22 * qJD(4) + t236 * t109 + t235 * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t218, -t73 ^ 2 + t75 ^ 2, -t234, -t218, t233, t159, -t22 * t75 + t235, t22 * t73 - t236, 0, 0;];
tauc_reg  = t5;
