% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:05:03
% EndTime: 2019-03-09 09:05:09
% DurationCPUTime: 2.00s
% Computational Cost: add. (2548->236), mult. (7362->449), div. (0->0), fcn. (7255->10), ass. (0->143)
t101 = cos(qJ(5));
t102 = cos(qJ(2));
t94 = sin(pkin(6));
t156 = t102 * t94;
t99 = sin(qJ(2));
t168 = t94 * t99;
t93 = sin(pkin(11));
t95 = cos(pkin(11));
t68 = -t95 * t156 + t93 * t168;
t96 = cos(pkin(6));
t98 = sin(qJ(5));
t114 = t68 * t101 - t96 * t98;
t179 = qJD(5) * t114;
t173 = pkin(3) + pkin(9);
t165 = pkin(8) + qJ(3);
t128 = t165 * t99;
t57 = (pkin(1) * t102 + pkin(2)) * t96 - t94 * t128;
t172 = pkin(1) * t96;
t144 = t99 * t172;
t62 = t165 * t156 + t144;
t39 = t95 * t57 - t93 * t62;
t69 = (t102 * t93 + t95 * t99) * t94;
t24 = t69 * pkin(4) - t173 * t96 - t39;
t129 = -pkin(2) * t102 - pkin(1);
t106 = -t69 * qJ(4) + t129 * t94;
t31 = t173 * t68 + t106;
t178 = t101 * t31 + t98 * t24;
t154 = qJD(2) * t99;
t88 = t94 ^ 2;
t177 = t154 * t88;
t100 = cos(qJ(6));
t120 = pkin(5) * t101 + pkin(10) * t98;
t176 = t100 * t120;
t147 = qJD(5) * t101;
t97 = sin(qJ(6));
t150 = qJD(6) * t97;
t73 = -t100 * t147 + t98 * t150;
t145 = qJD(6) * t101;
t134 = t97 * t145;
t152 = qJD(5) * t98;
t56 = t96 * t101 + t68 * t98;
t65 = qJD(2) * t69;
t38 = qJD(5) * t56 - t65 * t101;
t158 = t101 * t38;
t175 = -t100 * (-t114 * t152 - t158) + t114 * t134;
t149 = qJD(2) * t102;
t131 = t94 * t149;
t137 = t94 * t154;
t66 = t131 * t95 - t137 * t93;
t81 = pkin(2) * t137;
t107 = -t66 * qJ(4) - t69 * qJD(4) + t81;
t20 = t173 * t65 + t107;
t82 = t149 * t172;
t50 = t82 + (-qJD(2) * t128 + qJD(3) * t102) * t94;
t51 = -t62 * qJD(2) - qJD(3) * t168;
t29 = t93 * t50 - t95 * t51;
t21 = t66 * pkin(4) + t29;
t6 = -qJD(5) * t178 + t101 * t21 - t98 * t20;
t174 = 0.2e1 * qJD(4);
t118 = t69 * t100 - t56 * t97;
t37 = t65 * t98 + t179;
t16 = qJD(6) * t118 + t37 * t100 + t66 * t97;
t171 = t16 * t97;
t170 = t29 * t69;
t85 = -t95 * pkin(2) - pkin(3);
t83 = -pkin(9) + t85;
t169 = t83 * t98;
t167 = t97 * t98;
t42 = t56 * t100 + t69 * t97;
t15 = qJD(6) * t42 - t66 * t100 + t37 * t97;
t166 = t98 * t15;
t164 = t42 * t147 + t16 * t98;
t30 = t95 * t50 + t93 * t51;
t40 = t93 * t57 + t95 * t62;
t91 = t100 ^ 2;
t162 = t97 ^ 2 - t91;
t90 = t98 ^ 2;
t92 = t101 ^ 2;
t161 = t90 - t92;
t160 = t90 + t92;
t159 = t100 * t83;
t157 = t101 * t83;
t155 = t66 * t101;
t153 = qJD(5) * t118;
t151 = qJD(6) * t92;
t148 = qJD(5) * t100;
t146 = qJD(6) * t100;
t143 = -0.2e1 * pkin(5) * qJD(6);
t142 = t83 * t167;
t86 = t96 * qJD(4);
t26 = -t86 - t30;
t34 = -t96 * qJ(4) - t40;
t141 = t98 * t159;
t140 = t97 * t157;
t139 = t83 * t151;
t136 = t97 * t152;
t135 = t88 * t149;
t133 = t98 * t148;
t132 = t98 * t147;
t130 = t97 * t146;
t127 = t100 * t145;
t84 = t93 * pkin(2) + qJ(4);
t125 = t162 * qJD(6);
t124 = t161 * qJD(5);
t122 = t97 * t133;
t121 = t98 * pkin(5) - t101 * pkin(10);
t11 = t69 * pkin(10) + t178;
t25 = -t68 * pkin(4) - t34;
t14 = -pkin(5) * t114 - t56 * pkin(10) + t25;
t8 = t100 * t11 + t97 * t14;
t119 = -t100 * t118 + t42 * t97;
t115 = t101 * t24 - t98 * t31;
t17 = -t65 * pkin(4) - t26;
t77 = t121 + t84;
t53 = t97 * t77 + t141;
t10 = -t69 * pkin(5) - t115;
t4 = -t66 * pkin(5) - t6;
t113 = t10 * t150 - t4 * t100;
t112 = t10 * t146 + t4 * t97;
t110 = t152 * t69 - t155;
t44 = -t147 * t69 - t66 * t98;
t109 = -t100 * t38 - t114 * t150;
t108 = -t114 * t146 + t97 * t38;
t5 = -t101 * t20 - t24 * t147 + t152 * t31 - t98 * t21;
t105 = -t66 * pkin(10) + t5;
t104 = -t108 - t153;
t103 = t38 * pkin(5) - t37 * pkin(10) + t17;
t76 = -t127 + t136;
t75 = t146 * t98 + t147 * t97;
t74 = t133 + t134;
t71 = (-pkin(8) * t156 - t144) * qJD(2);
t70 = pkin(8) * t137 - t82;
t52 = t100 * t77 - t142;
t49 = t114 * t136;
t43 = t68 * pkin(3) + t106;
t35 = -t96 * pkin(3) - t39;
t33 = t100 * qJD(4) - t53 * qJD(6) + (-t140 + t176) * qJD(5);
t32 = -t77 * t146 - t97 * (qJD(5) * t120 + qJD(4)) + t73 * t83;
t28 = t65 * pkin(3) + t107;
t7 = t100 * t14 - t97 * t11;
t2 = -qJD(6) * t8 + t100 * t103 + t105 * t97;
t1 = t100 * t105 - t103 * t97 + t11 * t150 - t14 * t146;
t3 = [0, 0, 0, 0.2e1 * t99 * t135, 0.2e1 * (t102 ^ 2 - t99 ^ 2) * t88 * qJD(2), 0.2e1 * t96 * t131, -0.2e1 * t96 * t137, 0, -0.2e1 * pkin(1) * t177 + 0.2e1 * t71 * t96, -0.2e1 * pkin(1) * t135 + 0.2e1 * t70 * t96, -0.2e1 * t30 * t68 - 0.2e1 * t39 * t66 - 0.2e1 * t40 * t65 + 0.2e1 * t170, 0.2e1 * pkin(2) * t129 * t177 - 0.2e1 * t39 * t29 + 0.2e1 * t40 * t30, 0.2e1 * t26 * t68 + 0.2e1 * t34 * t65 + 0.2e1 * t35 * t66 + 0.2e1 * t170, -0.2e1 * t28 * t68 + 0.2e1 * t29 * t96 - 0.2e1 * t43 * t65, -0.2e1 * t26 * t96 - 0.2e1 * t28 * t69 - 0.2e1 * t43 * t66, 0.2e1 * t34 * t26 + 0.2e1 * t43 * t28 + 0.2e1 * t35 * t29, 0.2e1 * t56 * t37, 0.2e1 * t114 * t37 - 0.2e1 * t56 * t38, 0.2e1 * t37 * t69 + 0.2e1 * t56 * t66, 0.2e1 * t114 * t66 - 0.2e1 * t38 * t69, 0.2e1 * t69 * t66, -0.2e1 * t114 * t17 + 0.2e1 * t115 * t66 + 0.2e1 * t25 * t38 + 0.2e1 * t6 * t69, 0.2e1 * t17 * t56 - 0.2e1 * t178 * t66 + 0.2e1 * t25 * t37 + 0.2e1 * t5 * t69, 0.2e1 * t42 * t16, 0.2e1 * t118 * t16 - 0.2e1 * t42 * t15, -0.2e1 * t114 * t16 + 0.2e1 * t42 * t38, 0.2e1 * t114 * t15 + 0.2e1 * t118 * t38, -0.2e1 * t114 * t38, 0.2e1 * t10 * t15 - 0.2e1 * t114 * t2 - 0.2e1 * t118 * t4 + 0.2e1 * t7 * t38, -0.2e1 * t1 * t114 + 0.2e1 * t10 * t16 - 0.2e1 * t8 * t38 + 0.2e1 * t4 * t42; 0, 0, 0, 0, 0, t131, -t137, 0, t71, t70 (-t65 * t93 - t66 * t95) * pkin(2) (-t29 * t95 + t30 * t93) * pkin(2), -qJD(4) * t68 - t84 * t65 + t85 * t66, t29, 0.2e1 * t86 + t30, -t34 * qJD(4) - t26 * t84 + t29 * t85, t37 * t101 - t152 * t56, -t158 - t37 * t98 + (-t101 * t56 - t114 * t98) * qJD(5), -t110, t44, 0, t83 * t155 - qJD(4) * t114 + t17 * t98 + t84 * t38 + (t101 * t25 - t69 * t169) * qJD(5), -t66 * t169 + qJD(4) * t56 + t17 * t101 + t84 * t37 + (-t69 * t157 - t25 * t98) * qJD(5), -t42 * t134 + (t101 * t16 - t152 * t42) * t100, t119 * t152 + (-t100 * t15 - t171 + (-t100 * t42 - t118 * t97) * qJD(6)) * t101, t164 + t175, -t166 - t49 + (-t108 + t153) * t101, -t114 * t147 + t38 * t98, -t33 * t114 + t52 * t38 + (t2 + (-t10 * t97 - t118 * t83) * qJD(5)) * t98 + (qJD(5) * t7 - t15 * t83 + t112) * t101, -t32 * t114 - t53 * t38 + (t1 + (-t10 * t100 + t42 * t83) * qJD(5)) * t98 + (-qJD(5) * t8 - t16 * t83 - t113) * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, t84 * t174, -0.2e1 * t132, 0.2e1 * t124, 0, 0, 0, 0.2e1 * qJD(4) * t98 + 0.2e1 * t147 * t84, 0.2e1 * qJD(4) * t101 - 0.2e1 * t152 * t84, -0.2e1 * t130 * t92 - 0.2e1 * t132 * t91, 0.4e1 * t101 * t122 + 0.2e1 * t162 * t151, -0.2e1 * t98 * t134 - 0.2e1 * t161 * t148, 0.2e1 * t124 * t97 - 0.2e1 * t127 * t98, 0.2e1 * t132, -0.2e1 * t100 * t139 + 0.2e1 * t33 * t98 + 0.2e1 * (t52 + 0.2e1 * t142) * t147, 0.2e1 * t97 * t139 + 0.2e1 * t32 * t98 + 0.2e1 * (-t53 + 0.2e1 * t141) * t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0, -t65, -t66, t28, 0, 0, 0, 0, 0, t44, t110, 0, 0, 0, 0, 0, t101 * t104 + t166 - t49, t164 - t175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, t29, 0, 0, 0, 0, 0, -t110, t44, 0, 0, 0, 0, 0 (t97 * t179 - t15) * t101 + t104 * t98 (t114 * t148 - t16) * t101 + (qJD(5) * t42 + t109) * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160 * t146, t160 * t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t38, t66, t6, t5, t146 * t42 + t171, -qJD(6) * t119 + t16 * t100 - t97 * t15, t108, -t109, 0, -pkin(5) * t15 - pkin(10) * t108 + t113, -pkin(5) * t16 + pkin(10) * t109 + t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, -t147, 0, -t83 * t152, -t83 * t147, -t101 * t125 - t122, -0.4e1 * t97 * t127 + t162 * t152, t75, -t73, 0 (-t140 - t176) * qJD(6) + (t121 * t97 - t141) * qJD(5) (pkin(10) * t167 + (pkin(5) * t97 - t159) * t101) * qJD(6) + (t100 * t121 + t142) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, t152, 0, 0, 0, 0, 0, t73, t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, -t147, 0, 0, 0, 0, 0, -t74, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t130, -0.2e1 * t125, 0, 0, 0, t97 * t143, t100 * t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t38, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t76, t147, t33, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, -t150, 0, -pkin(10) * t146, pkin(10) * t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
