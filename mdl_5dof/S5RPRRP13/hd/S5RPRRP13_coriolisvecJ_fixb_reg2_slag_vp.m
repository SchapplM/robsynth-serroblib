% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP13_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP13_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:37
% EndTime: 2019-12-31 18:59:44
% DurationCPUTime: 1.95s
% Computational Cost: add. (2012->275), mult. (4272->373), div. (0->0), fcn. (2366->4), ass. (0->154)
t79 = cos(qJ(3));
t149 = qJD(3) * t79;
t76 = sin(qJ(4));
t131 = t76 * t149;
t78 = cos(qJ(4));
t147 = qJD(4) * t78;
t148 = qJD(4) * t76;
t77 = sin(qJ(3));
t112 = pkin(3) * t79 + pkin(7) * t77;
t56 = t112 * qJD(3) + qJD(2);
t43 = t56 * qJD(1);
t63 = pkin(3) * t77 - pkin(7) * t79 + qJ(2);
t49 = t63 * qJD(1);
t80 = -pkin(1) - pkin(6);
t68 = t80 * qJD(1) + qJD(2);
t62 = t77 * t68;
t51 = qJD(3) * pkin(7) + t62;
t119 = t68 * t131 + t51 * t147 + t49 * t148 - t78 * t43;
t20 = t49 * t76 + t51 * t78;
t153 = qJD(1) * t77;
t202 = -t153 - qJD(4);
t92 = -t20 * t202 - t119;
t150 = qJD(3) * t77;
t152 = qJD(1) * t79;
t133 = t78 * t152;
t151 = qJD(3) * t76;
t59 = t133 + t151;
t172 = t59 * t78;
t132 = t76 * t150;
t66 = qJD(1) * t132;
t32 = t59 * qJD(4) - t66;
t175 = t32 * t78;
t145 = t78 * qJD(3);
t146 = qJD(4) * t79;
t127 = t76 * t146;
t89 = t77 * t145 + t127;
t31 = t89 * qJD(1) - qJD(4) * t145;
t177 = t31 * t76;
t57 = t76 * t152 - t145;
t201 = (qJD(4) * (t57 * t76 - t172) - t175 + t177) * t79 + (t57 * t78 + t59 * t76) * t150;
t125 = t78 * t146;
t169 = t202 * t77;
t75 = t79 ^ 2;
t99 = qJD(1) * t75 + t169;
t200 = (t57 * t79 + t76 * t99) * qJD(3) - t202 * t125 + t77 * t32;
t118 = t57 + t145;
t199 = (t118 * t79 + t76 * t169) * qJD(1) + t202 * t148;
t116 = qJD(4) * t77 + qJD(1);
t161 = t79 * t31;
t162 = t78 * t79;
t170 = t202 * t76;
t198 = ((-t59 + t133) * t77 - t202 * t162) * qJD(3) + t116 * t170 - t161;
t97 = t59 * t202;
t98 = t57 * t202;
t197 = (t31 - t98) * t78 + (t32 - t97) * t76;
t135 = 0.2e1 * qJD(1);
t16 = -qJ(5) * t202 + t20;
t142 = qJD(1) * qJD(3);
t70 = t79 * t142;
t115 = pkin(4) * t70;
t2 = -t115 + t119;
t196 = t16 * t202 + t2;
t194 = t32 + t97;
t166 = t77 * t80;
t192 = t78 * t166 + t76 * t63;
t187 = qJD(4) * t192 - t78 * t56;
t184 = t59 ^ 2;
t183 = pkin(7) * t59;
t4 = pkin(4) * t32 + qJ(5) * t31 - qJD(5) * t59 + t68 * t150;
t182 = t4 * t76;
t181 = t4 * t78;
t155 = qJD(3) * pkin(3);
t171 = t68 * t79;
t52 = -t155 - t171;
t18 = pkin(4) * t57 - qJ(5) * t59 + t52;
t179 = t18 * t59;
t176 = t31 * t77;
t173 = t59 * t57;
t168 = t76 * t79;
t167 = t77 * t78;
t61 = t112 * qJD(1);
t164 = t78 * t61;
t163 = t78 * t63;
t81 = qJD(3) ^ 2;
t160 = t81 * t77;
t159 = t81 * t79;
t109 = pkin(4) * t76 - qJ(5) * t78;
t158 = t76 * qJD(5) + t202 * t109 + t62;
t25 = t68 * t162 + t76 * t61;
t157 = t77 ^ 2 - t75;
t82 = qJD(1) ^ 2;
t156 = -t81 - t82;
t154 = t82 * qJ(2);
t19 = t49 * t78 - t51 * t76;
t144 = qJD(5) - t19;
t143 = qJ(2) * qJD(3);
t141 = pkin(7) * t170;
t140 = pkin(7) * t202 * t78;
t137 = t79 * t82 * t77;
t129 = t80 * t149;
t136 = t78 * t129 + t63 * t147 + t76 * t56;
t134 = pkin(7) * t149;
t130 = t79 * t145;
t126 = t80 * t148;
t124 = t202 * t152;
t123 = qJD(2) * t135;
t122 = t57 ^ 2 - t184;
t121 = t76 * t80 - pkin(4);
t120 = -t52 + t171;
t117 = -t59 + t151;
t114 = t77 * t70;
t113 = qJ(5) * t70;
t110 = pkin(4) * t78 + qJ(5) * t76;
t15 = pkin(4) * t202 + t144;
t108 = t15 * t78 - t16 * t76;
t107 = t15 * t76 + t16 * t78;
t106 = t19 * t78 + t20 * t76;
t105 = t19 * t76 - t20 * t78;
t101 = (qJD(4) * t57 - t31) * pkin(7);
t100 = t120 * qJD(3);
t95 = t109 - t80;
t94 = -t18 * t77 + t134;
t93 = t52 * t77 - t134;
t91 = -t68 * t130 - t49 * t147 + t51 * t148 - t76 * t43;
t88 = -t76 * t98 - t175;
t1 = -qJD(5) * t202 + t113 - t91;
t87 = t108 * qJD(4) + t1 * t78 + t2 * t76;
t86 = -t106 * qJD(4) + t119 * t76 - t78 * t91;
t85 = t32 * t168 + (t125 - t132) * t57;
t84 = t57 * t150 + (-t32 - t66) * t79 - (-t116 * t78 - t131) * t202;
t83 = -t32 * t167 - t57 * t130 + t116 * t172 + (t116 * t57 + t59 * t149 - t176) * t76;
t72 = qJ(2) * t123;
t64 = -pkin(3) - t110;
t38 = t95 * t79;
t37 = (-t202 + t153) * t149;
t33 = -t76 * t166 + t163;
t30 = t121 * t77 - t163;
t29 = qJ(5) * t77 + t192;
t28 = pkin(7) * t175;
t26 = pkin(4) * t59 + qJ(5) * t57;
t24 = -t68 * t168 + t164;
t22 = -t164 + (-pkin(4) * qJD(1) + t68 * t76) * t79;
t21 = qJ(5) * t152 + t25;
t17 = -t31 - t98;
t14 = (t110 * qJD(4) - qJD(5) * t78) * t79 - t95 * t150;
t13 = -t76 * t129 - t187;
t12 = -t77 * t126 + t136;
t11 = -t202 * t147 + (t117 * t79 - t167 * t202) * qJD(1);
t10 = t121 * t149 + t187;
t9 = qJ(5) * t149 + (qJD(5) - t126) * t77 + t136;
t8 = -t78 * t97 - t177;
t7 = -t78 * t161 - t89 * t59;
t3 = t202 * t127 - t176 + (t59 * t79 + t99 * t78) * qJD(3);
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, t72, -0.2e1 * t114, 0.2e1 * t157 * t142, -t160, 0.2e1 * t114, -t159, 0, -t80 * t160 + (qJD(2) * t77 + t143 * t79) * t135, -t80 * t159 + (qJD(2) * t79 - t143 * t77) * t135, 0, t72, t7, t201, t3, t85, -t200, t37, -t13 * t202 - t119 * t77 + (t147 * t52 - t32 * t80) * t79 + ((qJD(1) * t33 + t19) * t79 + (t120 * t76 + t57 * t80) * t77) * qJD(3), t12 * t202 + t91 * t77 + (-t148 * t52 + t31 * t80) * t79 + ((-qJD(1) * t192 - t20) * t79 + (t120 * t78 + t59 * t80) * t77) * qJD(3), -t12 * t57 - t13 * t59 + t31 * t33 - t32 * t192 + t106 * t150 + (qJD(4) * t105 + t119 * t78 + t76 * t91) * t79, -t100 * t166 - t119 * t33 + t20 * t12 + t19 * t13 - t192 * t91, t7, t3, -t201, t37, t200, t85, t10 * t202 + t14 * t57 + t38 * t32 + (-t151 * t18 - t2) * t77 + (t18 * t147 + t182 + (-qJD(1) * t30 - t15) * qJD(3)) * t79, t10 * t59 - t29 * t32 - t30 * t31 - t57 * t9 - t108 * t150 + (-qJD(4) * t107 - t1 * t76 + t2 * t78) * t79, -t14 * t59 + t38 * t31 - t9 * t202 + (t145 * t18 + t1) * t77 + (t18 * t148 - t181 + (qJD(1) * t29 + t16) * qJD(3)) * t79, t1 * t29 + t10 * t15 + t14 * t18 + t16 * t9 + t2 * t30 + t38 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, -t154, 0, 0, 0, 0, 0, 0, t156 * t77, t156 * t79, 0, -t154, 0, 0, 0, 0, 0, 0, t84, -t198, t83, -t105 * t149 - t106 * qJD(1) + (-t100 + t86) * t77, 0, 0, 0, 0, 0, 0, t84, t83, t198, t108 * qJD(1) + (qJD(3) * t107 - t4) * t79 + (qJD(3) * t18 + t87) * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, -t157 * t82, 0, -t137, 0, 0, -t79 * t154, t77 * t154, 0, 0, t8, -t197, t11, t88, t199, t124, -pkin(3) * t32 + t24 * t202 - t118 * t62 + (t52 * t76 + t140) * qJD(4) + (-t19 * t79 + t76 * t93) * qJD(1), pkin(3) * t31 - t25 * t202 + t117 * t62 + (t52 * t78 - t141) * qJD(4) + (t20 * t79 + t78 * t93) * qJD(1), t24 * t59 + t25 * t57 - t28 + (-t19 * t153 - t91 + (-t19 + t183) * qJD(4)) * t78 + (t101 - t92) * t76, -t19 * t24 - t20 * t25 + (-t52 - t155) * t62 + t86 * pkin(7), t8, t11, t197, t124, -t199, t88, -t22 * t202 + t32 * t64 - t181 - t158 * t57 + (t18 * t76 + t140) * qJD(4) + (t15 * t79 - t76 * t94) * qJD(1), t21 * t57 - t22 * t59 - t28 + (t15 * t153 + t1 + (t15 + t183) * qJD(4)) * t78 + (t101 + t196) * t76, t21 * t202 + t31 * t64 - t182 + t158 * t59 + (-t18 * t78 + t141) * qJD(4) + (-t16 * t79 + t78 * t94) * qJD(1), pkin(7) * t87 - t15 * t22 - t158 * t18 - t16 * t21 + t4 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, -t122, t17, -t173, -t194, t70, -t52 * t59 + t92, -t19 * t202 + t52 * t57 + t91, 0, 0, t173, t17, t122, t70, t194, -t173, -t26 * t57 + 0.2e1 * t115 - t179 + t92, pkin(4) * t31 - t32 * qJ(5) + (t16 - t20) * t59 + (t15 - t144) * t57, 0.2e1 * t113 - t18 * t57 + t26 * t59 - (0.2e1 * qJD(5) - t19) * t202 - t91, -t2 * pkin(4) + t1 * qJ(5) + t144 * t16 - t15 * t20 - t18 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70 + t173, t17, -t202 ^ 2 - t184, t179 + t196;];
tauc_reg = t5;
