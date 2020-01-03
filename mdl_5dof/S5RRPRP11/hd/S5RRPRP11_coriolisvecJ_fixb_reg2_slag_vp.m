% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP11_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:14:04
% EndTime: 2019-12-31 20:14:10
% DurationCPUTime: 2.02s
% Computational Cost: add. (2242->309), mult. (5140->402), div. (0->0), fcn. (2809->4), ass. (0->172)
t206 = pkin(3) + pkin(6);
t114 = sin(qJ(2));
t177 = qJD(1) * t114;
t100 = qJD(4) + t177;
t113 = sin(qJ(4));
t115 = cos(qJ(4));
t170 = qJD(4) * t115;
t171 = qJD(4) * t113;
t116 = cos(qJ(2));
t145 = pkin(7) * t114 - qJ(3) * t116;
t168 = t114 * qJD(3);
t124 = t145 * qJD(2) - t168;
t167 = qJD(1) * qJD(2);
t156 = t114 * t167;
t99 = pkin(2) * t156;
t34 = t124 * qJD(1) + t99;
t117 = -pkin(2) - pkin(7);
t155 = -t114 * qJ(3) - pkin(1);
t66 = t117 * t116 + t155;
t47 = t66 * qJD(1);
t102 = pkin(6) * t177;
t211 = qJD(3) + t102;
t180 = pkin(3) * t177 + t211;
t51 = t117 * qJD(2) + t180;
t101 = t116 * t167;
t98 = pkin(6) * t101;
t65 = pkin(3) * t101 + t98;
t128 = -t113 * t65 - t115 * t34 - t51 * t170 + t47 * t171;
t17 = -t113 * t47 + t115 * t51;
t223 = -t100 * t17 - t128;
t154 = t113 * t34 - t115 * t65 + t47 * t170 + t51 * t171;
t18 = t113 * t51 + t115 * t47;
t222 = t100 * t18 - t154;
t174 = qJD(2) * t114;
t175 = qJD(2) * t113;
t176 = qJD(1) * t116;
t69 = t115 * t176 + t175;
t38 = qJD(4) * t69 - t113 * t156;
t193 = t115 * t38;
t157 = t113 * t176;
t39 = qJD(2) * t170 - qJD(4) * t157 - t115 * t156;
t35 = t39 * t113;
t173 = qJD(2) * t115;
t71 = -t157 + t173;
t221 = ((t113 * t71 + t115 * t69) * qJD(4) + t35 + t193) * t116 - (t113 * t69 - t115 * t71) * t174;
t112 = t116 ^ 2;
t186 = t100 * t114;
t134 = qJD(1) * t112 - t186;
t169 = qJD(4) * t116;
t160 = t113 * t169;
t220 = (t115 * t134 + t116 * t69) * qJD(2) - t100 * t160 + t114 * t39;
t135 = t100 * t71;
t136 = t69 * t100;
t219 = (t38 + t136) * t113 - (t39 + t135) * t115;
t218 = -0.2e1 * t167;
t14 = qJ(5) * t100 + t18;
t150 = pkin(4) * t101;
t2 = -t150 + t154;
t217 = t100 * t14 - t2;
t162 = t115 * t177;
t215 = (t162 + t170) * t100;
t89 = t206 * t114;
t199 = t113 * t89 + t115 * t66;
t105 = pkin(2) * t174;
t43 = t105 + t124;
t172 = qJD(2) * t116;
t78 = t206 * t172;
t10 = -qJD(4) * t199 - t113 * t43 + t115 * t78;
t208 = t71 ^ 2;
t207 = t100 ^ 2;
t110 = qJD(2) * qJ(3);
t103 = pkin(6) * t176;
t77 = pkin(3) * t176 + t103;
t61 = t110 + t77;
t20 = pkin(4) * t69 - qJ(5) * t71 + t61;
t205 = t20 * t71;
t109 = qJD(2) * qJD(3);
t76 = t206 * t174;
t53 = -qJD(1) * t76 + t109;
t6 = t39 * pkin(4) + t38 * qJ(5) - t71 * qJD(5) + t53;
t204 = t6 * t113;
t203 = t6 * t115;
t202 = t71 * t69;
t185 = t113 * t117;
t58 = t69 * t170;
t201 = -t117 * t58 - t39 * t185;
t147 = pkin(4) * t115 + qJ(5) * t113;
t132 = -pkin(3) - t147;
t200 = -t147 * qJD(4) + t115 * qJD(5) + t132 * t177 - t211;
t106 = pkin(2) * t177;
t56 = t145 * qJD(1) + t106;
t24 = t113 * t77 + t115 * t56;
t198 = qJD(2) * pkin(2);
t194 = t113 * t53;
t191 = t116 * t71;
t190 = t117 * t38;
t189 = t117 * t71;
t188 = t53 * t115;
t84 = -pkin(2) * t116 + t155;
t62 = qJD(1) * t84;
t184 = t114 * t115;
t119 = qJD(1) ^ 2;
t183 = t116 * t119;
t118 = qJD(2) ^ 2;
t182 = t118 * t114;
t181 = t118 * t116;
t179 = qJD(5) - t17;
t90 = t206 * t116;
t111 = t114 ^ 2;
t178 = t111 - t112;
t166 = t69 * t162 + t35 + t58;
t165 = t69 ^ 2 - t208;
t164 = t100 * t185;
t163 = t100 * t115 * t117;
t161 = t117 * t172;
t159 = t115 * t169;
t158 = t100 * t176;
t91 = t113 * t101;
t153 = -qJD(2) * t71 - t91;
t152 = pkin(1) * t218;
t151 = qJD(3) - t198;
t149 = t114 * t101;
t93 = t115 * t101;
t148 = qJ(5) * t101;
t146 = -pkin(4) * t113 + qJ(5) * t115;
t13 = -pkin(4) * t100 + t179;
t144 = t113 * t13 + t115 * t14;
t143 = -t113 * t17 + t115 * t18;
t23 = -t113 * t56 + t115 * t77;
t29 = -t113 * t66 + t115 * t89;
t137 = -0.2e1 * qJD(2) * t62;
t133 = t100 * t113;
t127 = -qJ(3) * t172 - t168;
t44 = t127 * qJD(1) + t99;
t59 = t105 + t127;
t131 = pkin(6) * t118 + qJD(1) * t59 + t44;
t129 = t69 * t176 - t215 - t91;
t9 = t113 * t78 + t115 * t43 + t89 * t170 - t66 * t171;
t125 = t39 - t135;
t123 = -qJD(2) * t69 - t100 * t133 + t93;
t122 = t113 * t135 - t166 + t193;
t121 = -t69 * t160 + (t116 * t39 - t69 * t174) * t115;
t80 = pkin(6) * t156 - t109;
t82 = t102 + t151;
t86 = -t103 - t110;
t120 = -t80 * t116 + (t116 * t82 + (t86 + t103) * t114) * qJD(2);
t97 = t114 * t183;
t88 = -0.2e1 * t149;
t87 = 0.2e1 * t149;
t85 = t178 * t119;
t83 = qJ(3) - t146;
t81 = t117 * t93;
t74 = -qJ(3) * t176 + t106;
t60 = t178 * t218;
t50 = t62 * t177;
t46 = (t100 + t177) * t172;
t42 = t147 * t116 + t90;
t32 = pkin(4) * t71 + qJ(5) * t69;
t26 = -pkin(4) * t114 - t29;
t25 = qJ(5) * t114 + t199;
t22 = -pkin(4) * t176 - t23;
t21 = qJ(5) * t176 + t24;
t19 = t136 - t38;
t16 = -t100 * t171 + t93 + (-t113 * t186 - t191) * qJD(1);
t15 = (t146 * qJD(4) + qJD(5) * t113) * t116 + (-pkin(6) + t132) * t174;
t12 = -t71 * t133 - t193;
t11 = -t71 * t159 + (t116 * t38 + t71 * t174) * t113;
t8 = -pkin(4) * t172 - t10;
t7 = -t100 * t159 - t38 * t114 + (-t134 * t113 + t191) * qJD(2);
t5 = qJ(5) * t172 + qJD(5) * t114 + t9;
t1 = qJD(5) * t100 - t128 + t148;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t60, t181, t88, -t182, 0, -pkin(6) * t181 + t114 * t152, pkin(6) * t182 + t116 * t152, 0, 0, 0, -t181, t182, t87, t60, t88, t120, t114 * t137 + t131 * t116, -t131 * t114 + t116 * t137, t120 * pkin(6) + t44 * t84 + t62 * t59, t11, t221, t7, t121, -t220, t46, t10 * t100 + t90 * t39 - t76 * t69 + (-t61 * t173 - t154) * t114 + (-t61 * t171 + t188 + (qJD(1) * t29 + t17) * qJD(2)) * t116, -t9 * t100 - t90 * t38 - t76 * t71 + (t61 * t175 + t128) * t114 + (-t61 * t170 - t194 + (-qJD(1) * t199 - t18) * qJD(2)) * t116, -t10 * t71 + t29 * t38 - t199 * t39 - t69 * t9 + t143 * t174 + (-t113 * t154 + t115 * t128 + (t113 * t18 + t115 * t17) * qJD(4)) * t116, t10 * t17 - t128 * t199 - t154 * t29 + t18 * t9 + t53 * t90 - t61 * t76, t11, t7, -t221, t46, t220, t121, -t8 * t100 + t15 * t69 + t42 * t39 + (-t173 * t20 - t2) * t114 + (-t20 * t171 + t203 + (-qJD(1) * t26 - t13) * qJD(2)) * t116, -t25 * t39 - t26 * t38 - t5 * t69 + t71 * t8 + t144 * t174 + (-t1 * t115 - t113 * t2 + (t113 * t14 - t115 * t13) * qJD(4)) * t116, t5 * t100 - t15 * t71 + t42 * t38 + (-t175 * t20 + t1) * t114 + (t20 * t170 + t204 + (qJD(1) * t25 + t14) * qJD(2)) * t116, t1 * t25 + t13 * t8 + t14 * t5 + t15 * t20 + t2 * t26 + t42 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, t85, 0, t97, 0, 0, t119 * pkin(1) * t114, pkin(1) * t183, 0, 0, 0, 0, 0, -t97, t85, t97, ((-t86 - t110) * t114 + (t151 - t82) * t116) * qJD(1), -t74 * t176 + t50, 0.2e1 * t109 + (t114 * t74 + t116 * t62) * qJD(1), -t80 * qJ(3) - t86 * qJD(3) - t62 * t74 + (-t114 * t86 + (-t82 - t198) * t116) * qJD(1) * pkin(6), t12, t219, t16, t166, t129, -t158, qJ(3) * t39 - t100 * t23 + t194 + t81 + t180 * t69 + (t115 * t61 - t164) * qJD(4) + (-t116 * t17 + t61 * t184) * qJD(1), -qJ(3) * t38 + t100 * t24 + t188 + t180 * t71 + (-t113 * t61 - t163) * qJD(4) + (t116 * t18 + (-t114 * t61 - t161) * t113) * qJD(1), t23 * t71 + t24 * t69 + (t190 - t222) * t115 + (t17 * t177 + t128 + (t17 + t189) * qJD(4)) * t113 + t201, t53 * qJ(3) - t17 * t23 - t18 * t24 + t180 * t61 + (qJD(4) * t143 - t113 * t128 - t115 * t154) * t117, t12, t16, -t219, -t158, -t129, t115 * t136 + t35, t100 * t22 + t204 + t39 * t83 + t81 - t200 * t69 + (t115 * t20 - t164) * qJD(4) + (t116 * t13 + t184 * t20) * qJD(1), t21 * t69 - t22 * t71 + (t190 - t217) * t115 + (-t13 * t177 - t1 + (-t13 + t189) * qJD(4)) * t113 + t201, -t100 * t21 - t203 + t38 * t83 + t200 * t71 + (t113 * t20 + t163) * qJD(4) + (-t116 * t14 + (t114 * t20 + t161) * t113) * qJD(1), -t13 * t22 - t14 * t21 + t6 * t83 - t200 * t20 + (qJD(4) * t144 + t1 * t113 - t2 * t115) * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, -t111 * t119 - t118, qJD(2) * t86 + t50 + t98, 0, 0, 0, 0, 0, 0, t123, -t207 * t115 + t153, t122, -t61 * qJD(2) + t223 * t113 + t115 * t222, 0, 0, 0, 0, 0, 0, t123, t122, -t153 + t215, -t20 * qJD(2) + t217 * t115 + (t100 * t13 + t1) * t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, -t165, t19, -t202, -t125, t101, -t61 * t71 + t222, t61 * t69 - t223, 0, 0, t202, t19, t165, t101, t125, -t202, -t32 * t69 + 0.2e1 * t150 - t205 + t222, pkin(4) * t38 - t39 * qJ(5) + (t14 - t18) * t71 + (t13 - t179) * t69, 0.2e1 * t148 - t20 * t69 + t32 * t71 + (0.2e1 * qJD(5) - t17) * t100 - t128, -t2 * pkin(4) + t1 * qJ(5) - t13 * t18 + t14 * t179 - t20 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101 + t202, t19, -t207 - t208, t205 - t217;];
tauc_reg = t3;
