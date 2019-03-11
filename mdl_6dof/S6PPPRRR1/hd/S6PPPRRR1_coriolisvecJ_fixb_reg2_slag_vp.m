% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PPPRRR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:40:35
% EndTime: 2019-03-08 18:40:42
% DurationCPUTime: 2.65s
% Computational Cost: add. (6595->297), mult. (17688->476), div. (0->0), fcn. (17335->16), ass. (0->166)
t99 = sin(qJ(5));
t195 = pkin(10) * t99;
t101 = cos(qJ(6));
t100 = sin(qJ(4));
t103 = cos(qJ(4));
t94 = cos(pkin(13));
t96 = cos(pkin(7));
t168 = t94 * t96;
t88 = sin(pkin(14));
t89 = sin(pkin(13));
t92 = sin(pkin(6));
t93 = cos(pkin(14));
t109 = (t93 * t168 - t88 * t89) * t92;
t91 = sin(pkin(7));
t170 = t91 * t93;
t97 = cos(pkin(6));
t81 = t97 * qJD(1) + qJD(2);
t52 = qJD(1) * t109 + t81 * t170;
t140 = t92 * t94 * t91;
t66 = -qJD(1) * t140 + t96 * t81 + qJD(3);
t90 = sin(pkin(8));
t95 = cos(pkin(8));
t121 = t52 * t95 + t66 * t90;
t171 = t89 * t93;
t172 = t88 * t91;
t53 = t81 * t172 + (t88 * t168 + t171) * t92 * qJD(1);
t33 = t121 * t100 + t103 * t53;
t102 = cos(qJ(5));
t124 = pkin(5) * t99 - pkin(11) * t102;
t78 = t124 * qJD(5);
t26 = (t78 + t33) * qJD(4);
t31 = qJD(4) * pkin(10) + t33;
t44 = -t90 * t52 + t95 * t66;
t22 = t102 * t31 + t99 * t44;
t20 = qJD(5) * pkin(11) + t22;
t107 = t100 * t53 - t121 * t103;
t79 = -t102 * pkin(5) - t99 * pkin(11) - pkin(4);
t27 = t79 * qJD(4) + t107;
t98 = sin(qJ(6));
t5 = t101 * t27 - t98 * t20;
t21 = t102 * t44 - t99 * t31;
t28 = t107 * qJD(4);
t7 = t21 * qJD(5) - t102 * t28;
t1 = qJD(6) * t5 + t101 * t7 + t98 * t26;
t143 = t102 * qJD(4);
t82 = -qJD(6) + t143;
t194 = t5 * t82 + t1;
t193 = qJD(6) * t102 * pkin(10) + t33 - t78;
t6 = t101 * t20 + t98 * t27;
t2 = -qJD(6) * t6 + t101 * t26 - t98 * t7;
t192 = t6 * t82 - t2;
t116 = t95 * t170 + t90 * t96;
t190 = -t100 * t172 + t116 * t103;
t169 = t91 * t97;
t56 = t93 * t169 + t109;
t71 = t97 * t96 - t140;
t120 = t56 * t95 + t71 * t90;
t57 = t92 * t171 + (t92 * t168 + t169) * t88;
t189 = -t57 * t100 + t120 * t103;
t146 = qJD(6) * t101;
t147 = qJD(5) * t102;
t108 = t99 * t146 + t98 * t147;
t141 = qJD(5) * qJD(6);
t65 = qJD(4) * t108 + t98 * t141;
t37 = t120 * t100 + t57 * t103;
t46 = -t56 * t90 + t71 * t95;
t119 = t46 * t102 - t37 * t99;
t8 = t22 * qJD(5) - t99 * t28;
t188 = t119 * t8;
t59 = t116 * t100 + t103 * t172;
t70 = -t90 * t170 + t95 * t96;
t118 = t102 * t70 - t99 * t59;
t185 = t8 * t118;
t163 = t100 * t90;
t72 = -t102 * t95 + t99 * t163;
t184 = t8 * t72;
t183 = t8 * t98;
t182 = t8 * t99;
t19 = -qJD(5) * pkin(5) - t21;
t181 = t19 * t98;
t29 = t33 * qJD(4);
t180 = t29 * t189;
t179 = t29 * t190;
t142 = qJD(4) * qJD(5);
t129 = t102 * t142;
t153 = qJD(4) * t99;
t136 = t98 * t153;
t64 = qJD(6) * t136 + (-t129 - t141) * t101;
t178 = t64 * t98;
t144 = t101 * qJD(5);
t74 = t136 - t144;
t149 = t98 * qJD(5);
t76 = t101 * t153 + t149;
t177 = t74 * t76;
t176 = t74 * t82;
t175 = t76 * t82;
t174 = t8 * t101;
t173 = t82 * t98;
t151 = qJD(6) * t98;
t160 = t102 * t98;
t167 = t193 * t101 + t107 * t160 - t149 * t195 + t79 * t151;
t150 = t101 * t102;
t166 = -t107 * t150 + t144 * t195 - t79 * t146 + t193 * t98;
t86 = t99 ^ 2;
t87 = t102 ^ 2;
t165 = t86 - t87;
t164 = qJD(4) * pkin(4);
t162 = t101 * t19;
t161 = t101 * t82;
t159 = t103 * t29;
t157 = t103 * t90;
t105 = qJD(4) ^ 2;
t156 = t105 * t90;
t154 = t65 * t101;
t152 = qJD(5) * t99;
t148 = qJD(4) * t100;
t138 = t99 * t151;
t137 = t100 * t156;
t135 = t99 * t105 * t102;
t134 = t90 * t148;
t133 = qJD(4) * t157;
t132 = t82 * t146;
t130 = t99 * t142;
t30 = t107 - t164;
t128 = -qJD(4) * t30 + t28;
t127 = t99 * t133;
t126 = t102 * t133;
t125 = t99 * t129;
t123 = -t101 * t5 - t6 * t98;
t24 = t37 * t102 + t46 * t99;
t11 = -t101 * t189 - t24 * t98;
t12 = t24 * t101 - t189 * t98;
t48 = t102 * t59 + t99 * t70;
t40 = -t101 * t190 - t98 * t48;
t41 = t101 * t48 - t190 * t98;
t117 = qJD(4) * t86 - t102 * t82;
t115 = -pkin(11) * t152 - t102 * t19;
t73 = t102 * t163 + t99 * t95;
t114 = -t101 * t73 + t98 * t157;
t62 = -t101 * t157 - t98 * t73;
t104 = qJD(5) ^ 2;
t113 = pkin(10) * t104;
t112 = qJD(5) * (t30 - t107 - t164);
t106 = t7 * t102 + t182 + (-t102 * t21 - t22 * t99) * qJD(5);
t77 = t124 * qJD(4);
t68 = pkin(10) * t150 + t98 * t79;
t67 = -pkin(10) * t160 + t101 * t79;
t61 = qJD(5) * t73 + t127;
t60 = -t72 * qJD(5) + t126;
t55 = t59 * qJD(4);
t54 = t190 * qJD(4);
t43 = qJD(6) * t114 + t101 * t134 - t98 * t60;
t42 = qJD(6) * t62 + t101 * t60 + t98 * t134;
t39 = t48 * qJD(5) + t99 * t54;
t38 = t118 * qJD(5) + t102 * t54;
t35 = t37 * qJD(4);
t34 = t189 * qJD(4);
t18 = t101 * t21 + t98 * t77;
t17 = t101 * t77 - t98 * t21;
t16 = -t41 * qJD(6) + t101 * t55 - t98 * t38;
t15 = t40 * qJD(6) + t101 * t38 + t98 * t55;
t10 = t119 * qJD(5) + t34 * t102;
t9 = t24 * qJD(5) + t34 * t99;
t4 = t11 * qJD(6) + t10 * t101 + t35 * t98;
t3 = -t12 * qJD(6) - t10 * t98 + t35 * t101;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35 * qJD(4), -t34 * qJD(4), 0, t107 * t35 - t28 * t37 + t33 * t34 - t180, 0, 0, 0, 0, 0, 0, -t9 * qJD(5) + (-t102 * t35 - t152 * t189) * qJD(4), -t10 * qJD(5) + (-t147 * t189 + t35 * t99) * qJD(4) (t10 * t102 + t9 * t99 + (-t102 * t119 - t24 * t99) * qJD(5)) * qJD(4), t10 * t22 - t21 * t9 + t24 * t7 + t30 * t35 - t180 - t188, 0, 0, 0, 0, 0, 0, t11 * t130 - t119 * t65 - t3 * t82 + t9 * t74, t119 * t64 - t12 * t130 + t4 * t82 + t9 * t76, t11 * t64 - t12 * t65 - t3 * t76 - t4 * t74, t1 * t12 + t11 * t2 + t19 * t9 + t3 * t5 + t4 * t6 - t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55 * qJD(4), -t54 * qJD(4), 0, t107 * t55 - t28 * t59 + t33 * t54 - t179, 0, 0, 0, 0, 0, 0, -t39 * qJD(5) + (-t102 * t55 - t152 * t190) * qJD(4), -t38 * qJD(5) + (-t147 * t190 + t55 * t99) * qJD(4) (t102 * t38 + t39 * t99 + (-t102 * t118 - t48 * t99) * qJD(5)) * qJD(4), -t21 * t39 + t22 * t38 + t30 * t55 + t7 * t48 - t179 - t185, 0, 0, 0, 0, 0, 0, -t118 * t65 + t130 * t40 - t16 * t82 + t39 * t74, t118 * t64 - t130 * t41 + t15 * t82 + t39 * t76, -t15 * t74 - t16 * t76 + t40 * t64 - t41 * t65, t1 * t41 + t6 * t15 + t5 * t16 + t19 * t39 + t2 * t40 - t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137, -t103 * t156, 0 (-t100 * t28 - t159 + (t100 * t107 + t103 * t33) * qJD(4)) * t90, 0, 0, 0, 0, 0, 0, -t102 * t137 + (-t61 - t127) * qJD(5), t99 * t137 + (-t60 - t126) * qJD(5) (t102 * t60 + t61 * t99 + (t102 * t72 - t73 * t99) * qJD(5)) * qJD(4), -t21 * t61 + t22 * t60 + t7 * t73 + t184 + (t30 * t148 - t159) * t90, 0, 0, 0, 0, 0, 0, t130 * t62 - t43 * t82 + t61 * t74 + t72 * t65, t114 * t130 + t42 * t82 + t61 * t76 - t72 * t64, t114 * t65 - t42 * t74 - t43 * t76 + t62 * t64, -t1 * t114 + t19 * t61 + t2 * t62 + t6 * t42 + t5 * t43 + t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t125, -0.2e1 * t165 * t142, t104 * t102, -0.2e1 * t125, -t104 * t99, 0, -t102 * t113 + t112 * t99, t102 * t112 + t113 * t99 -(-t86 - t87) * t28 + t106, -t29 * pkin(4) - t30 * t33 - (-t102 * t22 + t21 * t99) * t107 + t106 * pkin(10), -t76 * t138 + (t76 * t147 - t64 * t99) * t101 (-t101 * t74 - t76 * t98) * t147 + (-t154 + t178 + (-t101 * t76 + t74 * t98) * qJD(6)) * t99, t82 * t138 + t102 * t64 + (t117 * t101 + t76 * t99) * qJD(5), t65 * t98 * t99 + t108 * t74, t99 * t132 + t65 * t102 + (-t117 * t98 - t74 * t99) * qJD(5) (-t82 - t143) * t152, t167 * t82 + (-t2 + (pkin(10) * t74 + t181) * qJD(5)) * t102 + (t19 * t146 + pkin(10) * t65 + t107 * t74 + t183 + (qJD(4) * t67 + t5) * qJD(5)) * t99, -t166 * t82 + (t1 + (pkin(10) * t76 + t162) * qJD(5)) * t102 + (-t19 * t151 - pkin(10) * t64 + t174 + t107 * t76 + (-qJD(4) * t68 - t6) * qJD(5)) * t99, t67 * t64 - t68 * t65 + t167 * t76 + t166 * t74 + t123 * t147 + (-t1 * t98 - t101 * t2 + (-t101 * t6 + t5 * t98) * qJD(6)) * t99, t19 * t99 * t107 + t1 * t68 + t2 * t67 - t166 * t6 - t167 * t5 + (t19 * t147 + t182) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, t165 * t105, 0, t135, 0, 0, t128 * t99, t128 * t102, 0, 0, -t76 * t161 - t178 (-t65 + t175) * t98 + (-t64 + t176) * t101, -t132 + (t82 * t150 + (-t76 + t149) * t99) * qJD(4), -t74 * t173 - t154, t82 * t151 + (-t82 * t160 + (t74 + t144) * t99) * qJD(4), t82 * t153, -pkin(5) * t65 - t174 + t17 * t82 - t22 * t74 + (pkin(11) * t161 + t181) * qJD(6) + (t115 * t98 - t5 * t99) * qJD(4), pkin(5) * t64 - t18 * t82 - t22 * t76 + t183 + (-pkin(11) * t173 + t162) * qJD(6) + (t101 * t115 + t6 * t99) * qJD(4), t17 * t76 + t18 * t74 + ((qJD(6) * t74 - t64) * pkin(11) + t192) * t98 + ((qJD(6) * t76 - t65) * pkin(11) + t194) * t101, -t8 * pkin(5) - t5 * t17 - t6 * t18 - t19 * t22 + (qJD(6) * t123 + t1 * t101 - t2 * t98) * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, -t74 ^ 2 + t76 ^ 2, -t64 - t176, -t177, -t175 - t65, t130, -t19 * t76 - t192, t19 * t74 - t194, 0, 0;];
tauc_reg  = t13;
