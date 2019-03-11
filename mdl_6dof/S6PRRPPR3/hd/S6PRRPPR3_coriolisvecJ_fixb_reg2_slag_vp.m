% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPPR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:10:57
% EndTime: 2019-03-08 21:11:05
% DurationCPUTime: 2.89s
% Computational Cost: add. (2547->365), mult. (5994->472), div. (0->0), fcn. (3789->8), ass. (0->190)
t115 = sin(qJ(3));
t116 = sin(qJ(2));
t110 = sin(pkin(6));
t191 = qJD(1) * t110;
t165 = t116 * t191;
t67 = qJD(2) * pkin(8) + t165;
t207 = t115 * t67;
t118 = cos(qJ(3));
t111 = cos(pkin(6));
t190 = qJD(1) * t111;
t80 = t118 * t190;
t152 = -t80 + t207;
t232 = qJD(4) + t152;
t39 = -qJD(3) * pkin(3) + t232;
t186 = qJD(2) * t118;
t114 = sin(qJ(6));
t117 = cos(qJ(6));
t177 = t117 * qJD(3);
t60 = t114 * t186 - t177;
t188 = qJD(2) * t115;
t85 = qJD(6) + t188;
t133 = t60 * t85;
t180 = qJD(6) * t118;
t162 = t114 * t180;
t41 = -qJD(6) * t177 + (t115 * t177 + t162) * qJD(2);
t228 = t41 - t133;
t179 = t114 * qJD(3);
t61 = t117 * t186 + t179;
t217 = t61 * t85;
t176 = qJD(2) * qJD(3);
t156 = t115 * t176;
t42 = qJD(6) * t61 - t114 * t156;
t129 = t42 - t217;
t216 = pkin(8) - qJ(5);
t143 = pkin(5) * t115 + pkin(9) * t118;
t112 = qJD(2) * pkin(2);
t119 = cos(qJ(2));
t79 = t119 * t191;
t68 = -t79 - t112;
t48 = -pkin(3) * t186 - qJ(4) * t188 + t68;
t31 = pkin(4) * t186 + qJD(5) - t48;
t17 = qJD(2) * t143 + t31;
t120 = -pkin(3) - pkin(4);
t106 = -pkin(9) + t120;
t89 = qJ(5) * t188;
t132 = t232 - t89;
t21 = qJD(3) * t106 + t132;
t139 = t114 * t21 - t117 * t17;
t130 = pkin(5) * t118 + t106 * t115;
t127 = t130 * qJD(3);
t155 = t118 * t176;
t97 = t115 * qJD(4);
t214 = qJ(4) * t155 + qJD(2) * t97;
t15 = (t127 - t165) * qJD(2) + t214;
t184 = qJD(3) * t118;
t166 = qJ(5) * t184;
t178 = t115 * qJD(5);
t189 = qJD(2) * t110;
t158 = qJD(1) * t189;
t146 = t119 * t158;
t185 = qJD(3) * t115;
t157 = qJD(1) * t185;
t23 = t111 * t157 + t115 * t146 + t184 * t67;
t16 = (-t166 - t178) * qJD(2) + t23;
t1 = -t139 * qJD(6) + t114 * t15 + t117 * t16;
t231 = t139 * t85 + t1;
t6 = t114 * t17 + t117 * t21;
t2 = -qJD(6) * t6 - t114 * t16 + t117 * t15;
t230 = -t6 * t85 - t2;
t107 = qJD(3) * qJ(4);
t160 = t115 * t190;
t36 = -(qJ(5) * qJD(2) - t67) * t118 + t160;
t229 = t107 + t36;
t46 = t118 * t67 + t160;
t40 = t107 + t46;
t167 = qJD(3) * t120;
t227 = (-t152 + t207) * qJD(3);
t163 = t119 * t189;
t147 = t115 * t163;
t200 = t110 * t116;
t56 = t111 * t115 + t118 * t200;
t34 = qJD(3) * t56 + t147;
t169 = t115 * t200;
t55 = -t111 * t118 + t169;
t226 = -(t115 * t56 - t118 * t55) * qJD(3) + t115 * t34;
t69 = -pkin(3) * t118 - qJ(4) * t115 - pkin(2);
t58 = pkin(4) * t118 - t69;
t47 = t143 + t58;
t72 = t216 * t115;
t19 = -t114 * t72 + t117 * t47;
t197 = t117 * t119;
t213 = qJ(4) * t184 + t97;
t27 = t127 + t213;
t52 = t184 * t216 - t178;
t223 = t19 * qJD(6) + t114 * t27 + t117 * t52 - (-t114 * t116 + t115 * t197) * t191;
t198 = t114 * t119;
t20 = t114 * t47 + t117 * t72;
t222 = -t20 * qJD(6) - t114 * t52 + t117 * t27 - (-t115 * t198 - t116 * t117) * t191;
t105 = qJD(3) * qJD(4);
t215 = qJD(3) * t80 + t118 * t146;
t22 = -t185 * t67 + t215;
t18 = t105 + t22;
t183 = qJD(5) * t118;
t82 = qJ(5) * t156;
t11 = qJD(2) * t183 - t18 - t82;
t221 = t11 * t56;
t73 = t216 * t118;
t220 = t11 * t73;
t219 = t23 * t55;
t218 = t60 * t61;
t211 = t11 * t114;
t210 = t11 * t117;
t209 = t115 * t23;
t206 = t115 * t85;
t205 = t117 * t85;
t204 = t118 * t229;
t203 = t41 * t114;
t202 = t42 * t117;
t109 = t118 ^ 2;
t122 = qJD(2) ^ 2;
t201 = t109 * t122;
t199 = t110 * t122;
t121 = qJD(3) ^ 2;
t196 = t121 * t115;
t99 = t121 * t118;
t35 = t152 - t89;
t195 = qJD(4) + t35;
t193 = -qJD(5) - t31;
t108 = t115 ^ 2;
t192 = t108 + t109;
t187 = qJD(2) * t116;
t182 = qJD(6) * t114;
t181 = qJD(6) * t117;
t175 = 0.2e1 * t105 + t215;
t174 = pkin(3) * t185;
t173 = t114 * t206;
t172 = t115 * t205;
t171 = t85 * t182;
t170 = t85 * t181;
t168 = t116 * t199;
t164 = t110 * t187;
t161 = t117 * t180;
t30 = (t165 + t174) * qJD(2) - t214;
t159 = -pkin(8) * t121 - t30;
t154 = t68 - t112;
t151 = qJD(2) * t69 + t48;
t150 = t193 * t115;
t148 = t115 * t167;
t145 = t115 * t155;
t142 = -t114 * t6 + t117 * t139;
t141 = t114 * t139 + t117 * t6;
t137 = qJD(2) * t109 - t206;
t136 = qJD(3) * t46 - t23;
t37 = t110 * t197 - t114 * t55;
t38 = t110 * t198 + t117 * t55;
t26 = qJD(3) * pkin(5) + t229;
t131 = -t106 * t184 - t115 * t26;
t128 = t192 * t146;
t126 = qJD(6) * t141 + t1 * t114 + t117 * t2;
t125 = t209 + t118 * t18 + (-t115 * t40 + t118 * t39) * qJD(3);
t124 = t209 + t118 * t22 + (-t115 * t46 + t118 * t152) * qJD(3);
t33 = -qJD(3) * t169 + (qJD(3) * t111 + t163) * t118;
t14 = -t115 * t168 + (t118 * t163 + t33) * qJD(3);
t13 = -t118 * t168 + (-t34 - t147) * qJD(3);
t123 = qJD(2) * t226 + t186 * t33;
t113 = qJ(4) + pkin(5);
t96 = t108 * t122;
t91 = qJ(4) * t186;
t84 = t115 * t122 * t118;
t81 = -t96 - t121;
t75 = -0.2e1 * t145;
t74 = 0.2e1 * t145;
t71 = -t96 + t201;
t65 = t116 * t115 * t158;
t63 = t110 * t119 * t157;
t62 = pkin(3) * t188 - t91;
t59 = (-t108 + t109) * t176;
t57 = 0.2e1 * t59;
t53 = t174 - t213;
t51 = t185 * t216 + t183;
t50 = t120 * t188 + t91;
t49 = t148 + t213;
t32 = qJD(2) * t130 + t91;
t25 = t167 + t132;
t24 = (t148 - t165) * qJD(2) + t214;
t10 = t114 * t32 + t117 * t36;
t9 = -t114 * t36 + t117 * t32;
t8 = qJD(6) * t37 - t114 * t164 + t117 * t34;
t7 = -qJD(6) * t38 - t114 * t34 - t117 * t164;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t168, -t119 * t199, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t14, t123, t22 * t56 + t219 + t33 * t46 + t34 * t152 + (t68 - t79) * t164, 0, 0, 0, 0, 0, 0, t13, t123, t14, t18 * t56 + t219 + t33 * t40 + t34 * t39 + (-t119 * t30 + t187 * t48) * t110, 0, 0, 0, 0, 0, 0, t14, -t13 (-t118 * t33 - t226) * qJD(2), -t221 + t16 * t55 + t25 * t34 + t229 * t33 + (t119 * t24 - t187 * t31) * t110, 0, 0, 0, 0, 0, 0, t155 * t37 - t33 * t60 - t42 * t56 + t7 * t85, -t155 * t38 - t33 * t61 + t41 * t56 - t8 * t85, -t37 * t41 + t38 * t42 + t60 * t8 + t61 * t7, t1 * t38 - t139 * t7 + t2 * t37 + t26 * t33 + t6 * t8 - t221; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t57, t99, t75, -t196, 0, -pkin(8) * t99 + t154 * t185 + t63, pkin(8) * t196 + (t154 + t79) * t184, -t128 + t124 ((-t115 * t152 - t118 * t46) * t119 + (-t68 - t112) * t116) * t191 + t124 * pkin(8), t74, t99, -0.2e1 * t59, 0, t196, t75, t63 + t151 * t185 + ((-t53 + t165) * qJD(2) + t159) * t118, -t128 + t125, t65 + (-qJD(2) * t53 + t159) * t115 + (-t151 - t79) * t184, t30 * t69 + t48 * t53 + (-t116 * t48 + (-t115 * t39 - t118 * t40) * t119) * t191 + t125 * pkin(8), t75, t57, -t196, t74, t99, 0, t65 + (qJD(2) * t49 + t24) * t115 + (-t51 + (qJD(2) * t58 + t31 - t79) * t118) * qJD(3), -t118 * t24 - t63 + (t115 * t31 + t52) * qJD(3) + (t58 * t185 + (-t49 - t165) * t118) * qJD(2), t11 * t118 - t115 * t16 + (t115 * t229 - t118 * t25) * qJD(3) + (-t115 * t52 + t118 * t51 + (t115 * t73 - t118 * t72) * qJD(3) + t192 * t79) * qJD(2), -t220 + t16 * t72 + t24 * t58 + t25 * t52 - t229 * t51 + t31 * t49 + (t116 * t31 + (-t115 * t25 - t204) * t119) * t191, -t61 * t162 + (-t118 * t41 - t185 * t61) * t117 (t114 * t61 + t117 * t60) * t185 + (t203 - t202 + (t114 * t60 - t117 * t61) * qJD(6)) * t118, t85 * t162 + t115 * t41 + (-t117 * t137 - t118 * t61) * qJD(3), t60 * t161 + (t118 * t42 - t185 * t60) * t114, t85 * t161 + t115 * t42 + (t114 * t137 + t118 * t60) * qJD(3) (t85 + t188) * t184, -t42 * t73 + t51 * t60 + t222 * t85 + (t179 * t26 + t2) * t115 + (t60 * t79 - t26 * t181 + t211 + (qJD(2) * t19 - t139) * qJD(3)) * t118, t41 * t73 + t51 * t61 - t223 * t85 + (t177 * t26 - t1) * t115 + (t61 * t79 + t26 * t182 + t210 + (-qJD(2) * t20 - t6) * qJD(3)) * t118, t118 * t126 + t142 * t185 - t19 * t41 + t20 * t42 + t222 * t61 + t223 * t60, t1 * t20 - t220 + t19 * t2 + t223 * t6 - t222 * t139 + (-t118 * t79 - t51) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, -t71, 0, t84, 0, 0, -t188 * t68 + t136, -t186 * t68 - t215 + t227, 0, 0, -t84, 0, t71, 0, 0, t84 (-t115 * t48 + t118 * t62) * qJD(2) + t136, 0, -t227 + (t115 * t62 + t118 * t48) * qJD(2) + t175, -pkin(3) * t23 + qJ(4) * t18 + t232 * t40 - t39 * t46 - t48 * t62, t84, -t71, 0, -t84, 0, 0, t82 + (t35 - t207) * qJD(3) + (-t115 * t50 + t118 * t193) * qJD(2) + t175, -qJD(3) * t36 + ((-qJ(5) * qJD(3) + t50) * t118 + t150) * qJD(2) + t23 (-t195 + t25 - t167) * t186, -qJ(4) * t11 + t120 * t16 + t195 * t229 - t25 * t36 - t31 * t50, t205 * t61 - t203 (-t41 - t133) * t117 + (-t42 - t217) * t114, -t170 + (-t172 + (t61 - t179) * t118) * qJD(2), t114 * t133 - t202, t171 + (t173 + (-t60 - t177) * t118) * qJD(2), -t85 * t186, -t210 - t113 * t42 - t85 * t9 - t195 * t60 + (-t106 * t205 - t114 * t26) * qJD(6) + (t114 * t131 + t118 * t139) * qJD(2), t10 * t85 + t211 + t113 * t41 - t195 * t61 + (t106 * t114 * t85 - t117 * t26) * qJD(6) + (t117 * t131 + t118 * t6) * qJD(2), -t10 * t60 - t61 * t9 + (-t139 * t188 + t106 * t42 - t1 + (-t106 * t61 - t139) * qJD(6)) * t117 + (t6 * t188 + t106 * t41 + t2 + (-t106 * t60 + t6) * qJD(6)) * t114, -t10 * t6 - t11 * t113 + t139 * t9 + t195 * t26 + (qJD(6) * t142 + t1 * t117 - t114 * t2) * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, 0, t81, -qJD(3) * t40 + t188 * t48 + t23, 0, 0, 0, 0, 0, 0, t81, t84, 0, -qJD(3) * t229 + (t150 - t166) * qJD(2) + t23, 0, 0, 0, 0, 0, 0, -t170 + qJD(3) * t60 + (-t118 * t179 - t172) * qJD(2), t171 + qJD(3) * t61 + (-t118 * t177 + t173) * qJD(2), t114 * t228 + t117 * t129, -qJD(3) * t26 + t114 * t230 + t117 * t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t155, 0.2e1 * t156, -t96 - t201 (-t165 + t204 + (t25 + t167) * t115) * qJD(2) + t214, 0, 0, 0, 0, 0, 0, -t171 + (-t173 + (-t60 + t177) * t118) * qJD(2), -t170 + (-t172 + (-t61 - t179) * t118) * qJD(2), t114 * t129 - t117 * t228 (t115 * t141 + t118 * t26) * qJD(2) + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t218, -t60 ^ 2 + t61 ^ 2, t228, -t218, t129, t155, t26 * t61 - t230, -t26 * t60 - t231, 0, 0;];
tauc_reg  = t3;
