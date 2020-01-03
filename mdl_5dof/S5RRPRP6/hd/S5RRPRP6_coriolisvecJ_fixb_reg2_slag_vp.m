% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:58:27
% EndTime: 2019-12-31 19:58:34
% DurationCPUTime: 2.15s
% Computational Cost: add. (3915->329), mult. (10154->441), div. (0->0), fcn. (7134->6), ass. (0->180)
t133 = sin(qJ(4));
t135 = cos(qJ(4));
t174 = qJD(4) * t135;
t175 = qJD(4) * t133;
t132 = sin(pkin(8));
t190 = cos(pkin(8));
t134 = sin(qJ(2));
t136 = cos(qJ(2));
t209 = -qJ(3) - pkin(6);
t162 = qJD(2) * t209;
t98 = t136 * qJD(3) + t134 * t162;
t92 = t98 * qJD(1);
t141 = -t134 * qJD(3) + t136 * t162;
t93 = t141 * qJD(1);
t49 = t132 * t93 + t190 * t92;
t172 = qJD(1) * qJD(2);
t165 = t134 * t172;
t123 = pkin(2) * t165;
t121 = t132 * t165;
t160 = t190 * t136;
t122 = qJD(1) * t160;
t142 = qJD(2) * t122 - t121;
t113 = t132 * t136 + t190 * t134;
t102 = t113 * qJD(2);
t95 = qJD(1) * t102;
t53 = t95 * pkin(3) - t142 * pkin(7) + t123;
t176 = qJD(1) * t134;
t100 = t132 * t176 - t122;
t128 = -t136 * pkin(2) - pkin(1);
t177 = qJD(1) * t128;
t117 = qJD(3) + t177;
t178 = qJD(1) * t113;
t54 = t100 * pkin(3) - pkin(7) * t178 + t117;
t119 = t209 * t134;
t115 = qJD(1) * t119;
t204 = qJD(2) * pkin(2);
t109 = t115 + t204;
t120 = t209 * t136;
t116 = qJD(1) * t120;
t161 = t190 * t116;
t75 = t132 * t109 - t161;
t70 = qJD(2) * pkin(7) + t75;
t159 = -t133 * t53 - t135 * t49 - t54 * t174 + t70 * t175;
t31 = -t133 * t70 + t135 * t54;
t97 = qJD(4) + t100;
t226 = -t31 * t97 - t159;
t225 = -0.2e1 * t172;
t32 = t133 * t54 + t135 * t70;
t11 = -qJD(4) * t32 - t133 * t49 + t135 * t53;
t173 = t135 * qJD(2);
t147 = qJD(4) * t173 + t135 * t142 - t175 * t178;
t139 = -qJ(5) * t147 + t11;
t218 = t95 * pkin(4);
t85 = t133 * qJD(2) + t135 * t178;
t2 = -t85 * qJD(5) + t139 + t218;
t83 = t133 * t178 - t173;
t25 = -t83 * qJ(5) + t32;
t216 = t25 * t97;
t224 = t2 + t216;
t223 = t32 * t97 + t11;
t222 = t135 * t95 - t97 * t175;
t51 = t85 * qJD(4) + t133 * t142;
t221 = t85 ^ 2;
t220 = t178 ^ 2;
t219 = t83 * pkin(4);
t217 = pkin(2) * t134;
t48 = t132 * t92 - t190 * t93;
t80 = -t190 * t119 - t132 * t120;
t213 = t48 * t80;
t212 = t83 * t97;
t211 = t85 * t83;
t210 = t85 * t97;
t24 = -t85 * qJ(5) + t31;
t16 = t97 * pkin(4) + t24;
t208 = t16 - t24;
t125 = t132 * pkin(2) + pkin(7);
t180 = qJ(5) + t125;
t157 = qJD(4) * t180;
t184 = t135 * qJ(5);
t169 = pkin(2) * t176;
t64 = pkin(3) * t178 + t100 * pkin(7) + t169;
t106 = t132 * t116;
t77 = t190 * t115 + t106;
t34 = -t133 * t77 + t135 * t64;
t207 = pkin(4) * t178 + t133 * qJD(5) + t100 * t184 + t135 * t157 + t34;
t185 = t133 * qJ(5);
t35 = t133 * t64 + t135 * t77;
t206 = -t135 * qJD(5) + t100 * t185 + t133 * t157 + t35;
t205 = -t133 * t51 - t83 * t174;
t143 = -t132 * t134 + t160;
t73 = -pkin(3) * t143 - t113 * pkin(7) + t128;
t81 = t132 * t119 - t190 * t120;
t78 = t135 * t81;
t39 = t133 * t73 + t78;
t203 = t178 * t83;
t202 = t133 * t83;
t201 = t133 * t85;
t200 = t133 * t95;
t164 = qJD(5) + t219;
t74 = t190 * t109 + t106;
t69 = -qJD(2) * pkin(3) - t74;
t41 = t164 + t69;
t199 = t135 * t41;
t198 = t135 * t83;
t28 = t51 * pkin(4) + t48;
t197 = t28 * t133;
t196 = t28 * t135;
t195 = t147 * t133;
t194 = t51 * t135;
t193 = t85 * t178;
t192 = t95 * t143;
t191 = t97 * t178;
t189 = t100 * t133;
t188 = t178 * t100;
t105 = t143 * qJD(2);
t187 = t105 * t133;
t186 = t105 * t135;
t138 = qJD(1) ^ 2;
t183 = t136 * t138;
t137 = qJD(2) ^ 2;
t182 = t137 * t134;
t181 = t137 * t136;
t179 = t134 ^ 2 - t136 ^ 2;
t63 = t132 * t141 + t190 * t98;
t170 = t134 * t204;
t65 = t102 * pkin(3) - t105 * pkin(7) + t170;
t171 = t133 * t65 + t135 * t63 + t73 * t174;
t167 = t134 * t183;
t166 = t113 * t174;
t62 = t132 * t98 - t190 * t141;
t163 = -t133 * t63 + t135 * t65;
t38 = -t133 * t81 + t135 * t73;
t158 = pkin(1) * t225;
t76 = t132 * t115 - t161;
t156 = t135 * t97;
t155 = 0.2e1 * t178;
t154 = t136 * t165;
t126 = -t190 * pkin(2) - pkin(3);
t153 = qJD(4) * t125 * t97 + t48;
t146 = t51 * qJ(5) + t159;
t3 = -t83 * qJD(5) - t146;
t152 = -t97 * t16 + t3;
t151 = -t133 * t32 - t135 * t31;
t150 = -t198 - t201;
t149 = -qJ(5) * t105 - qJD(5) * t113;
t148 = -t97 * t189 + t222;
t145 = t135 * t147 - t85 * t175;
t144 = -t125 * t95 + t97 * t69;
t118 = -t135 * pkin(4) + t126;
t111 = t180 * t135;
t110 = t180 * t133;
t99 = t100 ^ 2;
t82 = t83 ^ 2;
t60 = t133 * t113 * pkin(4) + t80;
t44 = -pkin(4) * t189 + t76;
t40 = t97 * t102 - t192;
t37 = -t82 + t221;
t36 = (t166 + t187) * pkin(4) + t62;
t33 = -t113 * t185 + t39;
t30 = t210 - t51;
t29 = t147 + t212;
t26 = -pkin(4) * t143 - t113 * t184 + t38;
t22 = -t97 ^ 2 * t135 - t193 - t200;
t21 = t97 * t156 - t193 + t200;
t20 = t148 + t203;
t19 = t148 - t203;
t18 = t97 * t202 - t194;
t17 = t85 * t156 + t195;
t15 = -t205 * t113 + t83 * t187;
t14 = t145 * t113 + t85 * t186;
t13 = -t39 * qJD(4) + t163;
t12 = -t81 * t175 + t171;
t9 = -t97 * t187 - t83 * t102 + t51 * t143 + (-t97 * t174 - t200) * t113;
t8 = t85 * t102 + t113 * t222 - t143 * t147 + t97 * t186;
t7 = -qJ(5) * t166 + (-qJD(4) * t81 + t149) * t133 + t171;
t6 = t102 * pkin(4) + t149 * t135 + (-t78 + (qJ(5) * t113 - t73) * t133) * qJD(4) + t163;
t5 = t150 * t100 + t145 + t205;
t4 = (-t198 + t201) * t100 - t145 + t205;
t1 = t150 * t105 + (-t195 - t194 + (-t135 * t85 + t202) * qJD(4)) * t113;
t10 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t154, t179 * t225, t181, -0.2e1 * t154, -t182, 0, -pkin(6) * t181 + t134 * t158, pkin(6) * t182 + t136 * t158, 0, 0, t105 * t178 + t113 * t142, -t105 * t100 - t102 * t178 - t113 * t95 + t142 * t143, t105 * qJD(2), t100 * t102 - t192, -t102 * qJD(2), 0, t117 * t102 + t128 * t95 + (-t62 + (-qJD(1) * t143 + t100) * t217) * qJD(2), t117 * t105 - t128 * t121 + (t122 * t128 + t155 * t217 - t63) * qJD(2), -t63 * t100 - t75 * t102 - t74 * t105 + t48 * t113 + t142 * t80 + t143 * t49 + t178 * t62 - t81 * t95, t213 + t49 * t81 - t74 * t62 + t75 * t63 + (t117 + t177) * t170, t14, t1, t8, t15, t9, t40, t69 * t187 + t31 * t102 - t11 * t143 + t13 * t97 + t38 * t95 + t80 * t51 + t62 * t83 + (t48 * t133 + t174 * t69) * t113, t69 * t186 - t159 * t143 - t32 * t102 - t12 * t97 - t39 * t95 + t80 * t147 + t62 * t85 + (t48 * t135 - t175 * t69) * t113, -t12 * t83 - t13 * t85 - t38 * t147 - t39 * t51 + t151 * t105 + (t159 * t133 - t11 * t135 + (t133 * t31 - t135 * t32) * qJD(4)) * t113, t11 * t38 + t32 * t12 + t31 * t13 - t159 * t39 + t69 * t62 + t213, t14, t1, t8, t15, t9, t40, t41 * t187 + t16 * t102 - t2 * t143 + t26 * t95 + t36 * t83 + t60 * t51 + t6 * t97 + (t174 * t41 + t197) * t113, t41 * t186 - t25 * t102 + t3 * t143 - t33 * t95 + t36 * t85 + t60 * t147 - t7 * t97 + (-t175 * t41 + t196) * t113, -t26 * t147 - t33 * t51 - t6 * t85 - t7 * t83 + (-t133 * t25 - t135 * t16) * t105 + (-t3 * t133 - t2 * t135 + (t133 * t16 - t135 * t25) * qJD(4)) * t113, t16 * t6 + t2 * t26 + t25 * t7 + t28 * t60 + t3 * t33 + t41 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t167, t179 * t138, 0, t167, 0, 0, t138 * pkin(1) * t134, pkin(1) * t183, 0, 0, t188, -t99 + t220, -t121 + (t122 + t100) * qJD(2), -t188, 0, 0, t76 * qJD(2) - t100 * t169 - t117 * t178 - t48, t77 * qJD(2) + t117 * t100 - t169 * t178 - t49, (t75 - t76) * t178 + (t77 - t74) * t100 + (-t132 * t95 - t190 * t142) * pkin(2), t74 * t76 - t75 * t77 + (-t117 * t176 + t132 * t49 - t190 * t48) * pkin(2), t17, t5, t21, t18, t20, -t191, t126 * t51 + t133 * t144 - t135 * t153 - t178 * t31 - t34 * t97 - t76 * t83, t126 * t147 + t133 * t153 + t135 * t144 + t178 * t32 + t35 * t97 - t76 * t85, t34 * t85 + t35 * t83 + (-t100 * t31 - t125 * t51 - t159 + (t125 * t85 - t31) * qJD(4)) * t135 + (-t100 * t32 + t125 * t147 - t11 + (t125 * t83 - t32) * qJD(4)) * t133, t48 * t126 - t31 * t34 - t32 * t35 - t69 * t76 + (qJD(4) * t151 - t11 * t133 - t135 * t159) * t125, t17, t5, t21, t18, t20, -t191, -t16 * t178 - t110 * t95 + t118 * t51 - t196 - t44 * t83 - t207 * t97 + (t100 * t41 + (t41 + t219) * qJD(4)) * t133, t100 * t199 + t25 * t178 - t111 * t95 + t118 * t147 + t197 - t44 * t85 + t206 * t97 + (pkin(4) * t201 + t199) * qJD(4), t110 * t147 - t111 * t51 - t133 * t224 + t152 * t135 + t206 * t83 + t207 * t85, -t2 * t110 + t3 * t111 + t28 * t118 + (pkin(4) * t175 - t44) * t41 - t206 * t25 - t207 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155 * qJD(2), -t121 + (t122 - t100) * qJD(2), -t99 - t220, t75 * t100 + t178 * t74 + t123, 0, 0, 0, 0, 0, 0, t19, t22, t4, t226 * t133 + t223 * t135 - t69 * t178, 0, 0, 0, 0, 0, 0, t19, t22, t4, t152 * t133 + t135 * t224 - t178 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, t37, t29, -t211, t30, t95, -t69 * t85 + t223, t69 * t83 - t226, 0, 0, t211, t37, t29, -t211, t30, t95, 0.2e1 * t218 + t216 + (-t164 - t41) * t85 + t139, -t221 * pkin(4) + t24 * t97 + (qJD(5) + t41) * t83 + t146, -pkin(4) * t147 - t208 * t83, t208 * t25 + (-t41 * t85 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51 + t210, t147 - t212, -t82 - t221, t16 * t85 + t25 * t83 + t28;];
tauc_reg = t10;
