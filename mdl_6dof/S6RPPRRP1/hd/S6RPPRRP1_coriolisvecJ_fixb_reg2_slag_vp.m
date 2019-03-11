% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:58:43
% EndTime: 2019-03-09 01:58:50
% DurationCPUTime: 2.35s
% Computational Cost: add. (5277->325), mult. (12700->411), div. (0->0), fcn. (9230->8), ass. (0->181)
t138 = cos(pkin(10));
t231 = cos(qJ(4));
t177 = t231 * t138;
t125 = qJD(1) * t177;
t136 = sin(pkin(10));
t141 = sin(qJ(4));
t187 = t141 * t136;
t109 = qJD(1) * t187 - t125;
t104 = qJD(5) + t109;
t140 = sin(qJ(5));
t142 = cos(qJ(5));
t184 = qJD(5) * t142;
t185 = qJD(5) * t140;
t154 = t177 - t187;
t148 = t154 * qJD(3);
t237 = qJD(1) * t148;
t127 = sin(pkin(9)) * pkin(1) + qJ(3);
t121 = t127 * qJD(1);
t132 = t138 * qJD(2);
t220 = pkin(7) * qJD(1);
t91 = t132 + (-t121 - t220) * t136;
t100 = t136 * qJD(2) + t138 * t121;
t92 = t138 * t220 + t100;
t239 = -t141 * t92 + t231 * t91;
t44 = qJD(4) * t239 + t237;
t57 = t141 * t91 + t231 * t92;
t54 = qJD(4) * pkin(8) + t57;
t120 = -cos(pkin(9)) * pkin(1) - pkin(3) * t138 - pkin(2);
t105 = t120 * qJD(1) + qJD(3);
t119 = t231 * t136 + t141 * t138;
t111 = t119 * qJD(1);
t64 = t109 * pkin(4) - t111 * pkin(8) + t105;
t114 = t119 * qJD(4);
t101 = qJD(1) * t114;
t176 = qJD(4) * t187;
t124 = qJD(1) * t176;
t150 = -qJD(4) * t125 + t124;
t73 = t101 * pkin(4) + t150 * pkin(8);
t171 = -t140 * t73 - t142 * t44 - t64 * t184 + t54 * t185;
t31 = -t140 * t54 + t142 * t64;
t167 = -t104 * t31 - t171;
t182 = t142 * qJD(4);
t158 = qJD(5) * t182 - t111 * t185 - t142 * t150;
t32 = t140 * t64 + t142 * t54;
t9 = -qJD(5) * t32 - t140 * t44 + t142 * t73;
t144 = -qJ(6) * t158 + t9;
t230 = t101 * pkin(5);
t95 = qJD(4) * t140 + t111 * t142;
t1 = -t95 * qJD(6) + t144 + t230;
t93 = t111 * t140 - t182;
t21 = -qJ(6) * t93 + t32;
t212 = t21 * t104;
t241 = t1 + t212;
t240 = t32 * t104 + t9;
t227 = pkin(7) + t127;
t115 = t227 * t136;
t116 = t227 * t138;
t238 = -t231 * t115 - t141 * t116;
t200 = qJD(5) * t95;
t66 = -t140 * t150 + t200;
t113 = -qJD(4) * t177 + t176;
t192 = t113 * t142;
t98 = t142 * t101;
t236 = -(t119 * t185 + t192) * t104 + t119 * t98;
t235 = -t111 * t114 - t150 * t154;
t234 = t95 ^ 2;
t233 = t111 ^ 2;
t232 = t93 * pkin(5);
t149 = t119 * qJD(3);
t45 = qJD(1) * t149 + t57 * qJD(4);
t229 = t45 * t238;
t228 = t95 * t93;
t226 = -qJ(6) - pkin(8);
t20 = -qJ(6) * t95 + t31;
t17 = pkin(5) * t104 + t20;
t225 = t17 - t20;
t84 = pkin(4) * t111 + pkin(8) * t109;
t36 = t140 * t84 + t142 * t239;
t224 = -t140 * t66 - t93 * t184;
t190 = t119 * t142;
t223 = -t66 * t190 + t93 * t192;
t222 = t95 * t114 - t154 * t158;
t82 = -t141 * t115 + t231 * t116;
t75 = t142 * t82;
t76 = -pkin(4) * t154 - pkin(8) * t119 + t120;
t42 = t140 * t76 + t75;
t221 = -t119 * t101 + t113 * t109;
t218 = t104 * t93;
t217 = t111 * t93;
t216 = t140 * t93;
t215 = t140 * t95;
t214 = t142 * t93;
t213 = t142 * t95;
t23 = t66 * pkin(5) + t45;
t211 = t23 * t140;
t210 = t23 * t142;
t208 = t45 * t154;
t207 = t45 * t140;
t206 = t45 * t142;
t205 = t158 * t140;
t204 = t95 * t104;
t203 = t95 * t111;
t172 = qJD(5) * t226;
t197 = t109 * t140;
t202 = -qJ(6) * t197 + t142 * qJD(6) + t140 * t172 - t36;
t196 = t109 * t142;
t35 = -t140 * t239 + t142 * t84;
t201 = -pkin(5) * t111 - qJ(6) * t196 - t140 * qJD(6) + t142 * t172 - t35;
t199 = t101 * t154;
t198 = t104 * t111;
t195 = t111 * t109;
t193 = t113 * t140;
t191 = t119 * t140;
t189 = t140 * t101;
t186 = t136 ^ 2 + t138 ^ 2;
t183 = t113 * qJD(4);
t59 = t238 * qJD(4) + t148;
t85 = pkin(4) * t114 + pkin(8) * t113;
t181 = t140 * t85 + t142 * t59 + t76 * t184;
t180 = t95 * t193;
t175 = t119 * t184;
t174 = qJD(6) + t232;
t173 = -t140 * t59 + t142 * t85;
t41 = -t140 * t82 + t142 * t76;
t170 = qJD(1) * t186;
t169 = t104 * t142;
t157 = t66 * qJ(6) + t171;
t2 = -qJD(6) * t93 - t157;
t168 = -t104 * t17 + t2;
t166 = -t114 * t93 + t154 * t66;
t165 = t140 * t21 + t142 * t17;
t164 = t140 * t17 - t142 * t21;
t163 = t140 * t32 + t142 * t31;
t162 = t140 * t31 - t142 * t32;
t161 = -t100 * t138 + t136 * (-t121 * t136 + t132);
t160 = qJ(6) * t113 - qJD(6) * t119;
t159 = t98 + (-t185 - t197) * t104;
t155 = t142 * t158 - t95 * t185;
t153 = t175 - t193;
t53 = -qJD(4) * pkin(4) - t239;
t151 = -pkin(8) * t101 + t104 * t53;
t146 = -t163 * qJD(5) - t9 * t140 - t142 * t171;
t145 = -t153 * t104 - t119 * t189;
t60 = t82 * qJD(4) + t149;
t130 = -pkin(5) * t142 - pkin(4);
t123 = t226 * t142;
t122 = t226 * t140;
t108 = t109 ^ 2;
t103 = t114 * qJD(4);
t90 = t93 ^ 2;
t61 = pkin(5) * t191 - t238;
t48 = t104 * t114 - t199;
t47 = -t90 + t234;
t46 = -pkin(5) * t197 + t57;
t40 = t174 + t53;
t39 = t204 - t66;
t38 = t158 + t218;
t37 = t153 * pkin(5) + t60;
t34 = -qJ(6) * t191 + t42;
t33 = -pkin(5) * t154 - qJ(6) * t190 + t41;
t29 = -t104 ^ 2 * t142 - t189 - t203;
t28 = t104 * t169 + t189 - t203;
t27 = t159 + t217;
t26 = t159 - t217;
t25 = t104 * t216 - t66 * t142;
t24 = t95 * t169 + t205;
t19 = -t224 * t119 - t93 * t193;
t18 = t155 * t119 - t95 * t192;
t16 = -t42 * qJD(5) + t173;
t15 = -t82 * t185 + t181;
t14 = t145 - t166;
t13 = t145 + t166;
t12 = t222 - t236;
t11 = t222 + t236;
t10 = -qJ(6) * t175 + (-qJD(5) * t82 + t160) * t140 + t181;
t7 = (-t214 - t215) * t109 + t155 + t224;
t6 = (-t214 + t215) * t109 - t155 + t224;
t5 = t114 * pkin(5) + t160 * t142 + (-t75 + (qJ(6) * t119 - t76) * t140) * qJD(5) + t173;
t4 = t180 + (-t205 + (-t213 + t216) * qJD(5)) * t119 + t223;
t3 = -t180 + (t205 + (t213 + t216) * qJD(5)) * t119 + t223;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t170 (t127 * t170 - t161) * qJD(3), -t111 * t113 - t119 * t150, t221 + t235, -t183, t109 * t114 - t199, -t103, 0, -qJD(4) * t60 + t101 * t120 + t105 * t114, -t59 * qJD(4) - t105 * t113 - t120 * t150, -t82 * t101 - t59 * t109 + t60 * t111 + t113 * t239 - t57 * t114 + t45 * t119 + t150 * t238 + t154 * t44, -t239 * t60 + t44 * t82 + t57 * t59 - t229, t18, t4, t11, t19, t13, t48, -t53 * t193 + t41 * t101 + t16 * t104 + t31 * t114 - t9 * t154 + t60 * t93 - t238 * t66 + (t184 * t53 + t207) * t119, -t53 * t192 - t42 * t101 - t15 * t104 - t32 * t114 - t171 * t154 + t60 * t95 - t238 * t158 + (-t185 * t53 + t206) * t119, -t15 * t93 - t16 * t95 - t41 * t158 - t42 * t66 + t163 * t113 + (qJD(5) * t162 + t140 * t171 - t9 * t142) * t119, t15 * t32 + t16 * t31 - t171 * t42 + t41 * t9 + t53 * t60 - t229, t18, t4, t11, t19, t13, t48, -t40 * t193 - t1 * t154 + t33 * t101 + t5 * t104 + t17 * t114 + t37 * t93 + t61 * t66 + (t184 * t40 + t211) * t119, -t40 * t192 - t10 * t104 - t34 * t101 - t21 * t114 + t2 * t154 + t37 * t95 + t61 * t158 + (-t185 * t40 + t210) * t119, -t10 * t93 - t33 * t158 - t34 * t66 - t5 * t95 + t165 * t113 + (qJD(5) * t164 - t1 * t142 - t2 * t140) * t119, t1 * t33 + t10 * t21 + t17 * t5 + t2 * t34 + t23 * t61 + t37 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, t183, t221 - t235, -t113 * t57 - t114 * t239 + t119 * t44 - t208, 0, 0, 0, 0, 0, 0, t14, t12, t3, t113 * t162 + t53 * t114 + t119 * t146 - t208, 0, 0, 0, 0, 0, 0, t14, t12, t3, t40 * t114 - t23 * t154 + t164 * t113 + (-qJD(5) * t165 - t1 * t140 + t2 * t142) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t186 * qJD(1) ^ 2, t161 * qJD(1), 0, 0, 0, 0, 0, 0, 0.2e1 * t111 * qJD(4), -t124 + (t125 - t109) * qJD(4), -t108 - t233, t109 * t57 + t111 * t239, 0, 0, 0, 0, 0, 0, t26, t29, t6, -t53 * t111 + t167 * t140 + t240 * t142, 0, 0, 0, 0, 0, 0, t26, t29, t6, -t40 * t111 + t168 * t140 + t241 * t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, -t108 + t233, -t124 + (t125 + t109) * qJD(4), -t195, 0, 0 -(qJD(3) + t105) * t111, t105 * t109 - t237, 0, 0, t24, t7, t28, t25, t27, -t198, -pkin(4) * t66 - t31 * t111 - t206 - t57 * t93 + (-pkin(8) * t184 - t35) * t104 + t151 * t140, -pkin(4) * t158 + t32 * t111 + t207 - t57 * t95 + (pkin(8) * t185 + t36) * t104 + t151 * t142, t35 * t95 + t36 * t93 + ((-t66 + t200) * pkin(8) + t167) * t142 + ((qJD(5) * t93 + t158) * pkin(8) - t240) * t140, -t45 * pkin(4) + pkin(8) * t146 - t31 * t35 - t32 * t36 - t53 * t57, t24, t7, t28, t25, t27, -t198, t122 * t101 - t17 * t111 + t130 * t66 - t210 - t46 * t93 + t201 * t104 + (t109 * t40 + (t40 + t232) * qJD(5)) * t140, t40 * t196 + t123 * t101 + t21 * t111 + t130 * t158 + t211 - t46 * t95 - t202 * t104 + (pkin(5) * t215 + t142 * t40) * qJD(5), -t122 * t158 + t123 * t66 - t241 * t140 + t168 * t142 - t201 * t95 - t202 * t93, t1 * t122 - t2 * t123 + t23 * t130 + (pkin(5) * t185 - t46) * t40 + t202 * t21 + t201 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, t47, t38, -t228, t39, t101, -t53 * t95 + t240, t53 * t93 - t167, 0, 0, t228, t47, t38, -t228, t39, t101, 0.2e1 * t230 + t212 + (-t174 - t40) * t95 + t144, -t234 * pkin(5) + t20 * t104 + (qJD(6) + t40) * t93 + t157, -pkin(5) * t158 - t225 * t93, t225 * t21 + (-t40 * t95 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 + t204, t158 - t218, -t90 - t234, t17 * t95 + t21 * t93 + t23;];
tauc_reg  = t8;
