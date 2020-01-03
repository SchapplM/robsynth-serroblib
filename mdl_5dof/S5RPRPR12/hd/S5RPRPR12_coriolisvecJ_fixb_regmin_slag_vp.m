% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR12_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:30:15
% EndTime: 2019-12-31 18:30:21
% DurationCPUTime: 1.71s
% Computational Cost: add. (2205->240), mult. (5992->346), div. (0->0), fcn. (4642->8), ass. (0->128)
t122 = cos(pkin(8));
t173 = cos(qJ(3));
t143 = qJD(1) * t173;
t111 = t122 * t143;
t120 = sin(pkin(8));
t124 = sin(qJ(3));
t156 = t124 * t120;
t145 = qJD(1) * t156;
t89 = -t111 + t145;
t84 = qJD(5) + t89;
t168 = pkin(6) + qJ(2);
t105 = t168 * t120;
t101 = qJD(1) * t105;
t107 = t168 * t122;
t102 = qJD(1) * t107;
t67 = -t124 * t101 + t173 * t102;
t180 = qJD(3) * t67;
t179 = qJD(5) - t84;
t119 = sin(pkin(9));
t121 = cos(pkin(9));
t155 = t124 * t122;
t100 = t173 * t120 + t155;
t91 = t100 * qJD(1);
t74 = -t121 * qJD(3) + t119 * t91;
t178 = t74 * t89;
t76 = qJD(3) * t119 + t121 * t91;
t177 = t76 * t89;
t123 = sin(qJ(5));
t125 = cos(qJ(5));
t157 = t119 * t125;
t99 = t121 * t123 + t157;
t164 = t84 * t99;
t176 = t125 * t74;
t136 = t123 * t74 - t125 * t76;
t175 = t136 * t84;
t66 = -t173 * t101 - t124 * t102;
t72 = t173 * t105 + t124 * t107;
t154 = t125 * t121;
t158 = t119 * t123;
t97 = -t154 + t158;
t165 = t84 * t97;
t95 = t100 * qJD(3);
t81 = qJD(1) * t95;
t174 = t165 * t84 - t99 * t81;
t110 = qJD(3) * t111;
t80 = -qJD(3) * t145 + t110;
t14 = -t136 * qJD(5) + t99 * t80;
t85 = t89 ^ 2;
t172 = pkin(7) * t121;
t35 = t123 * t76 + t176;
t171 = t35 * t91;
t170 = t136 * t91;
t167 = pkin(7) + qJ(4);
t34 = pkin(3) * t81 - qJ(4) * t80 - qJD(4) * t91;
t132 = t173 * t122 - t156;
t128 = t132 * qJD(2);
t38 = qJD(1) * t128 + (qJD(4) + t66) * qJD(3);
t10 = t119 * t34 + t121 * t38;
t94 = t132 * qJD(3);
t43 = pkin(3) * t95 - qJ(4) * t94 - qJD(4) * t100;
t48 = -t72 * qJD(3) + t128;
t16 = t119 * t43 + t121 * t48;
t115 = -pkin(2) * t122 - pkin(1);
t103 = t115 * qJD(1) + qJD(2);
t47 = t89 * pkin(3) - t91 * qJ(4) + t103;
t62 = qJD(3) * qJ(4) + t67;
t20 = t119 * t47 + t121 * t62;
t63 = pkin(3) * t91 + qJ(4) * t89;
t26 = t119 * t63 + t121 * t66;
t64 = -pkin(3) * t132 - qJ(4) * t100 + t115;
t73 = -t124 * t105 + t173 * t107;
t28 = t119 * t64 + t121 * t73;
t151 = qJD(5) * t125;
t166 = -t74 * t151 + t80 * t154;
t163 = t119 * t80;
t162 = t119 * t89;
t161 = t119 * t94;
t160 = t121 * t80;
t159 = t100 * t119;
t153 = t120 ^ 2 + t122 ^ 2;
t152 = qJD(5) * t100;
t150 = qJD(1) * qJD(2);
t9 = -t119 * t38 + t121 * t34;
t4 = pkin(4) * t81 - pkin(7) * t160 + t9;
t5 = -pkin(7) * t163 + t10;
t149 = -t123 * t5 + t125 * t4;
t148 = t20 * t89 + t9;
t19 = -t119 * t62 + t121 * t47;
t144 = -t19 * t89 + t10;
t15 = -t119 * t48 + t121 * t43;
t25 = -t119 * t66 + t121 * t63;
t27 = -t119 * t73 + t121 * t64;
t142 = t153 * qJD(1) ^ 2;
t41 = t120 * qJD(2) * t143 + t150 * t155 + t180;
t141 = -t164 * t84 - t97 * t81;
t140 = t123 * t4 + t125 * t5;
t11 = -pkin(7) * t74 + t20;
t7 = pkin(4) * t89 - pkin(7) * t76 + t19;
t2 = t11 * t125 + t123 * t7;
t139 = t11 * t123 - t125 * t7;
t17 = -pkin(4) * t132 - t100 * t172 + t27;
t21 = -pkin(7) * t159 + t28;
t138 = -t123 * t21 + t125 * t17;
t137 = t123 * t17 + t125 * t21;
t135 = -qJD(5) * t76 - t163;
t134 = 0.2e1 * t153 * t150;
t104 = t167 * t119;
t131 = pkin(7) * t162 - qJD(4) * t121 + qJD(5) * t104 + t26;
t106 = t167 * t121;
t130 = pkin(4) * t91 + qJD(4) * t119 + qJD(5) * t106 + t89 * t172 + t25;
t60 = -qJD(3) * pkin(3) + qJD(4) - t66;
t129 = t100 * t41 + t60 * t94 + t72 * t80;
t127 = -pkin(3) * t80 - qJ(4) * t81 + (-qJD(4) + t60) * t89;
t13 = t135 * t123 + t166;
t49 = t100 * qJD(2) + t73 * qJD(3);
t114 = -pkin(4) * t121 - pkin(3);
t57 = t97 * t100;
t56 = t99 * t100;
t50 = pkin(4) * t159 + t72;
t42 = -pkin(4) * t162 + t67;
t30 = t74 * pkin(4) + t60;
t29 = pkin(4) * t161 + t49;
t24 = pkin(4) * t163 + t41;
t23 = t94 * t157 - t152 * t158 + (t100 * t151 + t123 * t94) * t121;
t22 = -t99 * t152 - t97 * t94;
t8 = -pkin(7) * t161 + t16;
t6 = pkin(4) * t95 - t94 * t172 + t15;
t1 = [0, 0, 0, 0, 0, t134, qJ(2) * t134, t100 * t80 + t91 * t94, -t100 * t81 + t132 * t80 - t89 * t94 - t91 * t95, t94 * qJD(3), -t95 * qJD(3), 0, -qJD(3) * t49 + t103 * t95 + t115 * t81, -qJD(3) * t48 + t103 * t94 + t115 * t80, t129 * t119 - t132 * t9 + t15 * t89 + t19 * t95 + t27 * t81 + t49 * t74, t10 * t132 + t129 * t121 - t16 * t89 - t20 * t95 - t28 * t81 + t49 * t76, -t15 * t76 - t16 * t74 + (-t100 * t9 - t19 * t94 - t27 * t80) * t121 + (-t10 * t100 - t20 * t94 - t28 * t80) * t119, t10 * t28 + t15 * t19 + t16 * t20 + t27 * t9 + t41 * t72 + t49 * t60, -t13 * t57 - t136 * t22, -t13 * t56 + t136 * t23 + t14 * t57 - t22 * t35, -t13 * t132 - t136 * t95 + t22 * t84 - t57 * t81, t132 * t14 - t23 * t84 - t35 * t95 - t56 * t81, -t132 * t81 + t84 * t95, (-t123 * t8 + t125 * t6) * t84 + t138 * t81 - t149 * t132 - t139 * t95 + t29 * t35 + t50 * t14 + t24 * t56 + t30 * t23 + (t132 * t2 - t137 * t84) * qJD(5), -(t123 * t6 + t125 * t8) * t84 - t137 * t81 + t140 * t132 - t2 * t95 - t29 * t136 + t50 * t13 - t24 * t57 + t30 * t22 + (-t132 * t139 - t138 * t84) * qJD(5); 0, 0, 0, 0, 0, -t142, -qJ(2) * t142, 0, 0, 0, 0, 0, 0.2e1 * t91 * qJD(3), t110 + (-t89 - t145) * qJD(3), -t119 * t85 + t121 * t81 - t74 * t91, -t119 * t81 - t121 * t85 - t76 * t91, (-t160 - t178) * t121 + (-t163 + t177) * t119, t144 * t119 + t148 * t121 - t60 * t91, 0, 0, 0, 0, 0, t141 - t171, t170 + t174; 0, 0, 0, 0, 0, 0, 0, t91 * t89, t91 ^ 2 - t85, t110 + (t89 - t145) * qJD(3), 0, 0, -t103 * t91 + t180 - t41, t103 * t89 - t132 * t150, t127 * t119 - t121 * t41 - t19 * t91 - t25 * t89 - t67 * t74, t119 * t41 + t127 * t121 + t20 * t91 + t26 * t89 - t67 * t76, t25 * t76 + t26 * t74 + (-qJD(4) * t74 + t144) * t121 + (qJD(4) * t76 - t148) * t119, -pkin(3) * t41 - t19 * t25 - t20 * t26 - t60 * t67 + (-t119 * t19 + t121 * t20) * qJD(4) + (t10 * t121 - t9 * t119) * qJ(4), t13 * t99 + t136 * t165, -t13 * t97 + t136 * t164 - t99 * t14 + t165 * t35, t170 - t174, t141 + t171, -t84 * t91, (-t104 * t125 - t106 * t123) * t81 + t114 * t14 + t24 * t97 + t139 * t91 - t42 * t35 + (t131 * t123 - t130 * t125) * t84 + t164 * t30, -(-t104 * t123 + t106 * t125) * t81 + t114 * t13 + t24 * t99 + t2 * t91 + t42 * t136 + (t123 * t130 + t125 * t131) * t84 - t165 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163 + t177, t160 - t178, -t74 ^ 2 - t76 ^ 2, t19 * t76 + t20 * t74 + t41, 0, 0, 0, 0, 0, t14 - t175, -t84 * t176 + (-t76 * t84 + t135) * t123 + t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136 * t35, t136 ^ 2 - t35 ^ 2, t35 * t84 + t13, -t14 - t175, t81, t30 * t136 - t179 * t2 + t149, t179 * t139 + t30 * t35 - t140;];
tauc_reg = t1;
