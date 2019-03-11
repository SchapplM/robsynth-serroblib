% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRRRP3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:11:07
% EndTime: 2019-03-09 00:11:13
% DurationCPUTime: 2.11s
% Computational Cost: add. (2137->241), mult. (5663->438), div. (0->0), fcn. (5236->10), ass. (0->131)
t94 = cos(qJ(3));
t148 = qJD(3) * t94;
t91 = sin(qJ(4));
t130 = t91 * t148;
t93 = cos(qJ(4));
t146 = qJD(4) * t93;
t92 = sin(qJ(3));
t172 = t92 * t146 + t130;
t160 = t92 * t93;
t164 = pkin(8) * t91;
t113 = -pkin(3) * t94 - t92 * pkin(9);
t65 = -pkin(2) + t113;
t55 = t93 * t65;
t34 = -pkin(10) * t160 + t55 + (-pkin(4) - t164) * t94;
t163 = cos(qJ(5));
t159 = t93 * t94;
t79 = pkin(8) * t159;
t154 = t91 * t65 + t79;
t161 = t91 * t92;
t41 = -pkin(10) * t161 + t154;
t36 = t163 * t41;
t90 = sin(qJ(5));
t158 = t90 * t34 + t36;
t166 = pkin(9) + pkin(10);
t69 = t166 * t91;
t70 = t166 * t93;
t155 = t163 * t70 - t90 * t69;
t150 = sin(pkin(6));
t111 = t150 * cos(qJ(2));
t108 = qJD(2) * t111;
t110 = t150 * sin(qJ(2));
t151 = cos(pkin(6));
t46 = t92 * t110 - t151 * t94;
t167 = t46 * qJD(3) - t94 * t108;
t171 = -qJD(4) * t111 - t167;
t102 = t94 * t110 + t151 * t92;
t106 = qJD(2) * t110;
t170 = t102 * qJD(4) - t106;
t119 = t163 * qJD(5);
t169 = t163 * qJD(4) + t119;
t88 = t93 ^ 2;
t153 = t91 ^ 2 - t88;
t118 = t153 * qJD(4);
t168 = qJD(4) + qJD(5);
t145 = qJD(4) * t94;
t136 = t91 * t145;
t149 = qJD(3) * t93;
t104 = t92 * t149 + t136;
t112 = pkin(3) * t92 - pkin(9) * t94;
t63 = t112 * qJD(3);
t25 = t104 * pkin(8) - t65 * t146 - t91 * t63;
t165 = pkin(4) * t92;
t162 = t90 * t91;
t84 = qJD(3) * t92;
t131 = t91 * t84;
t156 = pkin(8) * t131 + t93 * t63;
t64 = pkin(4) * t161 + t92 * pkin(8);
t87 = t92 ^ 2;
t152 = -t94 ^ 2 + t87;
t147 = qJD(4) * t91;
t144 = qJD(5) * t90;
t143 = -0.2e1 * pkin(2) * qJD(3);
t142 = -0.2e1 * pkin(3) * qJD(4);
t141 = t90 * t161;
t83 = pkin(8) * t148;
t42 = pkin(4) * t172 + t83;
t140 = pkin(4) * t147;
t139 = pkin(4) * t144;
t138 = t163 * pkin(4);
t137 = t92 * t147;
t134 = t93 * t145;
t133 = t46 * t147;
t132 = t46 * t146;
t129 = t91 * t146;
t128 = t92 * t148;
t127 = t93 * t148;
t82 = -pkin(4) * t93 - pkin(3);
t126 = t163 * t93;
t125 = qJD(4) * t166;
t16 = (-pkin(10) * t159 + t165) * qJD(3) + (-t79 + (pkin(10) * t92 - t65) * t91) * qJD(4) + t156;
t20 = -pkin(10) * t172 - t25;
t124 = t163 * t16 - t90 * t20;
t123 = t163 * t34 - t90 * t41;
t122 = -t163 * t69 - t70 * t90;
t121 = qJD(3) * t163;
t117 = t152 * qJD(3);
t116 = t92 * t127;
t115 = pkin(4) * t119;
t114 = t94 * t121;
t109 = qJD(3) * t111;
t57 = t163 * t91 + t90 * t93;
t5 = -t34 * t119 + t41 * t144 - t90 * t16 - t163 * t20;
t61 = t91 * t125;
t62 = t93 * t125;
t23 = t69 * t119 + t70 * t144 + t163 * t61 + t90 * t62;
t6 = -t158 * qJD(5) + t124;
t24 = -t155 * qJD(5) - t163 * t62 + t90 * t61;
t38 = t168 * t57;
t100 = -t102 * t91 - t93 * t111;
t99 = t90 * t100;
t97 = t163 * t100;
t96 = t170 * t93 + t171 * t91;
t95 = -t170 * t91 + t171 * t93;
t81 = t138 + pkin(5);
t75 = -0.2e1 * t128;
t56 = -t126 + t162;
t45 = t92 * t126 - t141;
t44 = t57 * t92;
t43 = pkin(5) * t56 + t82;
t40 = t102 * t93 - t91 * t111;
t39 = t102 * qJD(3) + t92 * t108;
t37 = t168 * t162 - t169 * t93;
t35 = pkin(5) * t44 + t64;
t29 = pkin(5) * t38 + t140;
t28 = -qJ(6) * t56 + t155;
t27 = -qJ(6) * t57 + t122;
t26 = -t154 * qJD(4) + t156;
t22 = t91 * t114 - t90 * t137 - qJD(5) * t141 + (t90 * t148 + t169 * t92) * t93;
t21 = -t93 * t114 + t90 * t130 + t38 * t92;
t19 = t163 * t40 + t99;
t18 = -t90 * t40 + t97;
t11 = pkin(5) * t22 + t42;
t10 = -qJ(6) * t44 + t158;
t9 = -pkin(5) * t94 - t45 * qJ(6) + t123;
t8 = t37 * qJ(6) - t57 * qJD(6) + t24;
t7 = -qJ(6) * t38 - qJD(6) * t56 - t23;
t4 = -qJD(5) * t99 - t40 * t119 - t163 * t96 - t90 * t95;
t3 = -qJD(5) * t97 + t40 * t144 - t163 * t95 + t90 * t96;
t2 = -qJ(6) * t22 - qJD(6) * t44 - t5;
t1 = pkin(5) * t84 + t21 * qJ(6) - t45 * qJD(6) + t6;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t18 * t4 - 0.2e1 * t19 * t3 + 0.2e1 * t39 * t46; 0, 0, -t106, -t108, 0, 0, 0, 0, 0, -t94 * t106 - t92 * t109, t92 * t106 - t94 * t109, 0, 0, 0, 0, 0, t100 * t84 + t46 * t130 + t92 * t132 + t39 * t161 + t96 * t94, t46 * t127 - t92 * t133 + t39 * t160 - t40 * t84 + t95 * t94, 0, 0, 0, 0, 0, t18 * t84 + t46 * t22 + t39 * t44 - t4 * t94, -t19 * t84 - t46 * t21 - t3 * t94 + t39 * t45, t18 * t21 - t19 * t22 + t3 * t44 - t4 * t45, t1 * t18 - t10 * t3 + t11 * t46 + t19 * t2 + t35 * t39 + t4 * t9; 0, 0, 0, 0, 0.2e1 * t128, -0.2e1 * t117, 0, 0, 0, t92 * t143, t94 * t143, 0.2e1 * t88 * t128 - 0.2e1 * t87 * t129, -0.4e1 * t91 * t116 + 0.2e1 * t87 * t118, 0.2e1 * t92 * t136 + 0.2e1 * t152 * t149, -0.2e1 * t91 * t117 + 0.2e1 * t92 * t134, t75, 0.2e1 * t55 * t84 - 0.2e1 * t26 * t94 + 0.2e1 * (t91 * t128 + t87 * t146) * pkin(8), -0.2e1 * t25 * t94 - 0.2e1 * t154 * t84 + 0.2e1 * (-t87 * t147 + 0.2e1 * t116) * pkin(8), -0.2e1 * t45 * t21, 0.2e1 * t21 * t44 - 0.2e1 * t22 * t45, 0.2e1 * t21 * t94 + 0.2e1 * t45 * t84, 0.2e1 * t22 * t94 - 0.2e1 * t44 * t84, t75, 0.2e1 * t123 * t84 + 0.2e1 * t64 * t22 + 0.2e1 * t42 * t44 - 0.2e1 * t6 * t94, -0.2e1 * t158 * t84 - 0.2e1 * t64 * t21 + 0.2e1 * t42 * t45 - 0.2e1 * t5 * t94, -0.2e1 * t1 * t45 - 0.2e1 * t10 * t22 - 0.2e1 * t2 * t44 + 0.2e1 * t21 * t9, 0.2e1 * t1 * t9 + 0.2e1 * t10 * t2 + 0.2e1 * t11 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t167, 0, 0, 0, 0, 0, -t39 * t93 + t133, t39 * t91 + t132, 0, 0, 0, 0, 0, t38 * t46 + t39 * t56, -t37 * t46 + t39 * t57, t18 * t37 - t19 * t38 + t3 * t56 - t4 * t57, t18 * t8 + t19 * t7 + t27 * t4 - t28 * t3 + t29 * t46 + t39 * t43; 0, 0, 0, 0, 0, 0, t148, -t84, 0, -t83, pkin(8) * t84, -t92 * t118 + t91 * t127, -0.4e1 * t92 * t129 - t153 * t148, t131 - t134, t104, 0 (pkin(9) * t159 + (-pkin(3) * t93 + t164) * t92) * qJD(4) + (t113 * t91 - t79) * qJD(3) (pkin(8) * t160 + t112 * t91) * qJD(4) + (t113 * t93 + t94 * t164) * qJD(3), -t21 * t57 - t37 * t45, t21 * t56 - t22 * t57 + t37 * t44 - t38 * t45, t37 * t94 + t57 * t84, t38 * t94 - t56 * t84, 0, t122 * t84 + t44 * t140 + t82 * t22 - t24 * t94 + t64 * t38 + t42 * t56, t45 * t140 - t155 * t84 - t82 * t21 - t23 * t94 - t64 * t37 + t42 * t57, -t1 * t57 - t10 * t38 - t2 * t56 + t21 * t27 - t22 * t28 + t37 * t9 - t44 * t7 - t45 * t8, t1 * t27 + t10 * t7 + t11 * t43 + t2 * t28 + t29 * t35 + t8 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t129, -0.2e1 * t118, 0, 0, 0, t91 * t142, t93 * t142, -0.2e1 * t57 * t37, 0.2e1 * t37 * t56 - 0.2e1 * t38 * t57, 0, 0, 0, 0.2e1 * t56 * t140 + 0.2e1 * t38 * t82, 0.2e1 * t57 * t140 - 0.2e1 * t37 * t82, 0.2e1 * t27 * t37 - 0.2e1 * t28 * t38 - 0.2e1 * t56 * t7 - 0.2e1 * t57 * t8, 0.2e1 * t27 * t8 + 0.2e1 * t28 * t7 + 0.2e1 * t29 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t95, 0, 0, 0, 0, 0, t4, t3, 0, t4 * t81 + (-t3 * t90 + (t163 * t19 - t18 * t90) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127 - t137, -t172, t84, t26, t25, 0, 0, -t21, -t22, t84, t121 * t165 + (-t36 + (pkin(4) * t94 - t34) * t90) * qJD(5) + t124 (t119 * t94 - t90 * t84) * pkin(4) + t5, t81 * t21 + (-t22 * t90 + (-t163 * t44 + t45 * t90) * qJD(5)) * pkin(4), t1 * t81 + (t2 * t90 + (t163 * t10 - t9 * t90) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, -t147, 0, -pkin(9) * t146, pkin(9) * t147, 0, 0, -t37, -t38, 0, t24, t23, t81 * t37 + (-t38 * t90 + (-t163 * t56 + t57 * t90) * qJD(5)) * pkin(4), t8 * t81 + (t7 * t90 + (t163 * t28 - t27 * t90) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t139, -0.2e1 * t115, 0, 0.2e1 * (t138 - t81) * t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, t4 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t22, t84, t6, t5, pkin(5) * t21, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t38, 0, t24, t23, pkin(5) * t37, t8 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, -t115, 0, -pkin(5) * t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t12;
