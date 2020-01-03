% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:40
% EndTime: 2019-12-31 19:01:44
% DurationCPUTime: 1.28s
% Computational Cost: add. (1647->180), mult. (3902->266), div. (0->0), fcn. (2734->8), ass. (0->116)
t148 = cos(qJ(4));
t115 = qJD(1) * t148;
t86 = sin(qJ(3));
t129 = qJD(1) * t86;
t85 = sin(qJ(4));
t88 = cos(qJ(3));
t154 = t88 * t115 - t85 * t129;
t79 = qJD(3) + qJD(4);
t36 = t154 * t79;
t54 = qJD(5) - t154;
t153 = qJD(5) - t54;
t73 = sin(pkin(9)) * pkin(1) + pkin(6);
t149 = pkin(7) + t73;
t112 = t149 * qJD(1);
t50 = t88 * qJD(2) - t112 * t86;
t62 = t148 * t86 + t85 * t88;
t152 = qJD(1) * t62;
t42 = t79 * t62;
t37 = t42 * qJD(1);
t87 = cos(qJ(5));
t135 = t87 * t37;
t84 = sin(qJ(5));
t126 = qJD(5) * t84;
t94 = t148 * t88 - t85 * t86;
t41 = t79 * t94;
t96 = t62 * t126 - t87 * t41;
t151 = t62 * t135 - t96 * t54;
t51 = t86 * qJD(2) + t112 * t88;
t128 = qJD(1) * t88;
t57 = -t86 * t115 - t85 * t128;
t102 = t87 * t57 - t84 * t79;
t17 = -t102 * qJD(5) + t84 * t36;
t114 = qJD(3) * t149;
t52 = t86 * t114;
t53 = t88 * t114;
t59 = t149 * t86;
t60 = t149 * t88;
t95 = -t148 * t59 - t85 * t60;
t13 = t95 * qJD(4) - t148 * t52 - t85 * t53;
t136 = t85 * t51;
t130 = qJD(3) * pkin(3);
t46 = t50 + t130;
t22 = t148 * t46 - t136;
t20 = -t79 * pkin(4) - t22;
t74 = -cos(pkin(9)) * pkin(1) - pkin(2);
t63 = -t88 * pkin(3) + t74;
t58 = t63 * qJD(1);
t29 = -pkin(4) * t154 + t57 * pkin(8) + t58;
t33 = -pkin(4) * t94 - t62 * pkin(8) + t63;
t35 = t148 * t60 - t85 * t59;
t113 = t148 * qJD(4);
t127 = qJD(4) * t85;
t44 = t50 * qJD(3);
t45 = t51 * qJD(3);
t6 = t46 * t113 - t51 * t127 + t148 * t44 - t85 * t45;
t116 = t148 * t45 + t85 * t44;
t118 = t148 * t51;
t23 = t85 * t46 + t118;
t7 = t23 * qJD(4) + t116;
t150 = -(qJD(5) * t33 + t13) * t54 + (qJD(5) * t29 + t6) * t94 + t20 * t41 - t35 * t37 + t7 * t62;
t125 = qJD(5) * t87;
t16 = t79 * t125 + t57 * t126 + t87 * t36;
t147 = t16 * t84;
t146 = t20 * t154;
t145 = t20 * t62;
t144 = t33 * t37;
t143 = t41 * t79;
t47 = -t84 * t57 - t87 * t79;
t142 = t47 * t54;
t141 = t102 * t54;
t140 = t54 * t57;
t139 = t57 * t154;
t137 = t84 * t37;
t89 = qJD(3) ^ 2;
t134 = t89 * t86;
t133 = t89 * t88;
t132 = -t102 * t42 - t16 * t94;
t131 = t86 ^ 2 - t88 ^ 2;
t65 = qJD(1) * t74;
t123 = qJD(1) * qJD(3);
t121 = t86 * t130;
t120 = pkin(3) * t129;
t117 = t86 * t123;
t110 = t54 * t87;
t38 = -t57 * pkin(4) - pkin(8) * t154;
t76 = t85 * pkin(3) + pkin(8);
t108 = qJD(5) * t76 + t120 + t38;
t21 = t79 * pkin(8) + t23;
t9 = t87 * t21 + t84 * t29;
t107 = t20 * t125 - t9 * t57 + t7 * t84;
t27 = t85 * t50 + t118;
t106 = pkin(3) * t127 - t27;
t105 = t17 * t94 - t42 * t47;
t104 = -t76 * t37 - t146;
t103 = t84 * t21 - t87 * t29;
t100 = 0.2e1 * qJD(3) * t65;
t99 = -t103 * t57 + t20 * t126 - t7 * t87;
t28 = t148 * t50 - t136;
t98 = -pkin(3) * t113 + t28;
t97 = t58 * t57 - t116;
t92 = (-t62 * t125 - t84 * t41) * t54 - t62 * t137;
t91 = -t154 * t58 - t6;
t90 = qJD(1) ^ 2;
t77 = -t148 * pkin(3) - pkin(4);
t39 = t42 * t79;
t30 = -t154 ^ 2 + t57 ^ 2;
t26 = (-t152 - t57) * t79;
t15 = t42 * pkin(4) - t41 * pkin(8) + t121;
t14 = t35 * qJD(4) + t148 * t53 - t85 * t52;
t11 = pkin(3) * t117 + t37 * pkin(4) - t36 * pkin(8);
t10 = t87 * t11;
t4 = -t102 * t57 + t54 * t110 + t137;
t3 = -t54 ^ 2 * t84 - t47 * t57 + t135;
t2 = -t102 * t110 + t147;
t1 = (t16 - t142) * t87 + (-t17 + t141) * t84;
t5 = [0, 0, 0, 0, 0.2e1 * t88 * t117, -0.2e1 * t131 * t123, t133, -t134, 0, t86 * t100 - t73 * t133, t88 * t100 + t73 * t134, t36 * t62 - t57 * t41, t154 * t41 + t36 * t94 - t62 * t37 + t57 * t42, t143, -t39, 0, -t14 * t79 + t63 * t37 + t58 * t42 + (-qJD(1) * t94 - t154) * t121, -t13 * t79 + t63 * t36 + t58 * t41 + (-t57 + t152) * t121, t16 * t87 * t62 + t102 * t96, (t102 * t84 - t47 * t87) * t41 + (-t147 - t17 * t87 + (t102 * t87 + t47 * t84) * qJD(5)) * t62, t132 + t151, t105 + t92, -t37 * t94 + t54 * t42, -t10 * t94 + t14 * t47 - t95 * t17 - t103 * t42 + (t15 * t54 + t144 + (t21 * t94 - t35 * t54 + t145) * qJD(5)) * t87 + t150 * t84, -t14 * t102 - t95 * t16 - t9 * t42 + (-(-qJD(5) * t35 + t15) * t54 - t144 + (-qJD(5) * t21 + t11) * t94 - qJD(5) * t145) * t84 + t150 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, -t133, 0, 0, 0, 0, 0, -t39, -t143, 0, 0, 0, 0, 0, -t105 + t92, t132 - t151; 0, 0, 0, 0, -t86 * t90 * t88, t131 * t90, 0, 0, 0, -t65 * t129, -t65 * t128, t139, t30, 0, t26, 0, t154 * t120 + t27 * t79 + (-t118 + (-pkin(3) * t79 - t46) * t85) * qJD(4) + t97, t28 * t79 + (-t79 * t113 + t57 * t129) * pkin(3) + t91, t2, t1, t4, t3, t140, t77 * t17 + t104 * t84 + t106 * t47 + (-t108 * t87 + t98 * t84) * t54 + t99, t77 * t16 + t104 * t87 - t106 * t102 + (t108 * t84 + t98 * t87) * t54 + t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t30, 0, t26, 0, t97 + (-qJD(4) + t79) * t23, t22 * t79 + t91, t2, t1, t4, t3, t140, -pkin(4) * t17 - (-t84 * t22 + t87 * t38) * t54 - t23 * t47 - t84 * t146 + (-t54 * t125 - t137) * pkin(8) + t99, -pkin(4) * t16 + (t87 * t22 + t84 * t38) * t54 + t23 * t102 - t87 * t146 + (t54 * t126 - t135) * pkin(8) + t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102 * t47, t102 ^ 2 - t47 ^ 2, t16 + t142, -t17 - t141, t37, t20 * t102 - t153 * t9 - t84 * t6 + t10, t153 * t103 - t84 * t11 + t20 * t47 - t87 * t6;];
tauc_reg = t5;
