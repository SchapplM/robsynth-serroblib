% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:17
% EndTime: 2019-12-31 17:41:19
% DurationCPUTime: 1.06s
% Computational Cost: add. (777->228), mult. (1511->250), div. (0->0), fcn. (785->4), ass. (0->132)
t89 = cos(qJ(3));
t145 = qJ(4) * t89;
t158 = pkin(3) + pkin(4);
t88 = sin(qJ(3));
t161 = t158 * t88;
t103 = t145 - t161;
t163 = qJDD(3) * qJ(4) + qJD(3) * qJD(4);
t127 = qJD(2) * qJD(3);
t121 = t89 * t127;
t69 = t88 * qJDD(2);
t166 = t121 + t69;
t143 = qJD(2) * t88;
t66 = pkin(6) * t143;
t28 = t89 * qJD(1) - t66;
t165 = qJD(4) - t28;
t123 = t158 * qJD(3);
t16 = qJ(5) * t143 + t28;
t134 = qJD(4) - t16;
t9 = -t123 + t134;
t119 = t158 * qJDD(3);
t138 = t88 * qJD(1);
t29 = t89 * qJD(2) * pkin(6) + t138;
t83 = qJD(3) * qJ(4);
t24 = t29 + t83;
t80 = pkin(7) + qJ(2);
t67 = sin(t80);
t68 = cos(t80);
t149 = g(1) * t68 + g(2) * t67;
t164 = 0.2e1 * t163;
t151 = pkin(6) - qJ(5);
t36 = t151 * t89;
t18 = qJD(2) * t36 + t138;
t14 = t83 + t18;
t74 = t88 * qJ(4);
t120 = pkin(2) + t74;
t12 = qJD(5) + (t158 * t89 + t120) * qJD(2);
t133 = qJD(5) + t12;
t162 = t133 * t88;
t21 = -qJD(3) * pkin(3) + t165;
t154 = t68 * t88;
t155 = t67 * t88;
t160 = -g(1) * t154 - g(2) * t155 + g(3) * t89;
t144 = pkin(6) * qJDD(3);
t78 = t89 * pkin(3);
t104 = -t120 - t78;
t22 = t104 * qJD(2);
t148 = t78 + t74;
t30 = -pkin(2) - t148;
t159 = (qJD(2) * t30 + t22) * qJD(3) - t144;
t84 = t88 ^ 2;
t85 = t89 ^ 2;
t147 = -t84 + t85;
t71 = t89 * qJDD(2);
t10 = 0.2e1 * t147 * t127 + 0.2e1 * t88 * t71;
t54 = g(1) * t67;
t157 = g(2) * t68;
t77 = t89 * pkin(4);
t156 = t14 * t89;
t153 = t68 * t89;
t92 = qJD(2) ^ 2;
t152 = t88 * t92;
t150 = t68 * pkin(2) + t67 * pkin(6);
t146 = pkin(6) * qJD(3);
t142 = qJD(3) * t36;
t141 = qJD(3) * t88;
t87 = qJDD(2) * pkin(2);
t140 = qJDD(3) * pkin(3);
t139 = t14 * qJD(3);
t137 = t88 * qJD(4);
t136 = t88 * qJD(5);
t135 = t89 * qJD(5);
t132 = qJ(5) * qJD(3);
t131 = t84 * qJDD(2);
t130 = t85 * qJDD(2);
t129 = qJ(5) * qJDD(2);
t128 = qJD(1) * qJD(3);
t126 = pkin(6) * t71 + t88 * qJDD(1) + t89 * t128;
t125 = t77 + t148;
t124 = t54 - t157;
t122 = t88 * t127;
t118 = pkin(3) * t153 + t68 * t74 + t150;
t117 = qJ(5) * t122 + t126;
t116 = t166 * pkin(6) - t89 * qJDD(1) + t88 * t128;
t115 = -t149 + (t130 + t131) * pkin(6);
t23 = pkin(2) + t125;
t114 = qJD(2) * t23 + t12;
t111 = t88 * t121;
t110 = pkin(3) * t71 + t166 * qJ(4) + qJD(2) * t137 + t87;
t91 = qJD(3) ^ 2;
t109 = pkin(6) * t91 + t157;
t107 = pkin(3) * t88 - t145;
t106 = -qJDD(4) - t116;
t105 = g(3) * t88 - t126;
t102 = t116 + t160;
t101 = -0.2e1 * pkin(2) * t127 - t144;
t100 = pkin(4) * t71 + qJDD(5) + t110;
t99 = -t109 + 0.2e1 * t87;
t7 = -pkin(6) * t122 + t126;
t6 = -t106 - t140;
t11 = t103 * qJD(3) + t137;
t3 = -t158 * t122 + t100;
t98 = qJD(2) * t11 + qJDD(2) * t23 - t157 + t3;
t97 = t29 * qJD(3) - t102;
t20 = t107 * qJD(3) - t137;
t4 = pkin(3) * t122 - t110;
t96 = -qJD(2) * t20 - qJDD(2) * t30 - t109 - t4;
t95 = -t88 * t129 + qJDD(4) + t102;
t5 = t7 + t163;
t94 = t5 * t89 + t6 * t88 + (t21 * t89 - t24 * t88) * qJD(3);
t93 = t7 * t89 + t116 * t88 + (-t28 * t89 - t29 * t88) * qJD(3);
t86 = qJDD(1) - g(3);
t51 = t68 * pkin(6);
t48 = t89 * t152;
t45 = -t84 * t92 - t91;
t43 = t89 * t54;
t42 = g(1) * t155;
t39 = t68 * t145;
t37 = t67 * t145;
t35 = t151 * t88;
t34 = -qJDD(3) - t48;
t33 = t147 * t92;
t32 = qJDD(3) * t89 - t91 * t88;
t31 = qJDD(3) * t88 + t91 * t89;
t27 = t107 * qJD(2);
t26 = -0.2e1 * t111 + t130;
t25 = 0.2e1 * t111 + t131;
t19 = -t136 + t142;
t17 = -t151 * t141 - t135;
t15 = t103 * qJD(2);
t2 = -t89 * t129 + (-pkin(6) * t141 - t135) * qJD(2) + t117 + t163;
t1 = -t166 * qJ(5) - qJD(2) * t136 - t106 - t119;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, 0, 0, 0, 0, 0, t32, -t31, 0, t7 * t88 - t116 * t89 - g(3) + (-t28 * t88 + t29 * t89) * qJD(3), 0, 0, 0, 0, 0, 0, t32, 0, t31, t5 * t88 - t6 * t89 - g(3) + (t21 * t88 + t24 * t89) * qJD(3), 0, 0, 0, 0, 0, 0, t32, t31, 0, -t1 * t89 + t2 * t88 - g(3) + (t88 * t9 + t156) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t124, t149, 0, 0, t25, t10, t31, t26, t32, 0, t101 * t88 + t99 * t89 + t43, t101 * t89 - t99 * t88 - t42, t93 + t115, qJDD(2) * pkin(2) ^ 2 - g(1) * (-t67 * pkin(2) + t51) - g(2) * t150 + t93 * pkin(6), t25, t31, -t10, 0, -t32, t26, t159 * t88 + t96 * t89 + t43, t94 + t115, -t159 * t89 + t96 * t88 + t42, t94 * pkin(6) - g(1) * t51 - g(2) * t118 - t104 * t54 + t22 * t20 + t4 * t30, t25, -t10, -t31, t26, t32, 0, -t35 * qJDD(3) + t43 + (-t114 * t88 - t19) * qJD(3) + t98 * t89, t36 * qJDD(3) + t42 + (t114 * t89 + t17) * qJD(3) + t98 * t88, (-qJD(3) * t9 - qJDD(2) * t36 - t2 + (-qJD(3) * t35 - t17) * qJD(2)) * t89 + (t139 - qJDD(2) * t35 - t1 + (-t19 + t142) * qJD(2)) * t88 + t149, t2 * t36 + t14 * t17 + t1 * t35 + t9 * t19 + t3 * t23 + t12 * t11 - g(1) * (-t68 * qJ(5) + t51) - g(2) * (pkin(4) * t153 + t118) + (-g(1) * (t104 - t77) + g(2) * qJ(5)) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t33, t69, t48, t71, qJDD(3), pkin(2) * t152 + t97, (t28 + t66) * qJD(3) + (pkin(2) * t92 + t149) * t89 + t105, 0, 0, -t48, t69, t33, qJDD(3), -t71, t48, 0.2e1 * t140 - qJDD(4) + (-t22 * t88 + t27 * t89) * qJD(2) + t97, -t107 * qJDD(2), -t28 * qJD(3) - t149 * t89 + (t22 * t89 + (t27 - t146) * t88) * qJD(2) - t105 + t164, t5 * qJ(4) - t6 * pkin(3) - t22 * t27 - t21 * t29 - g(1) * (-pkin(3) * t154 + t39) - g(2) * (-pkin(3) * t155 + t37) - g(3) * t148 + t165 * t24, -t48, t33, -t69, t48, t71, qJDD(3), t18 * qJD(3) + 0.2e1 * t119 + ((-t15 + t132) * t89 + t162) * qJD(2) - t95, -t16 * qJD(3) + (-g(3) + (-t15 - t146) * qJD(2)) * t88 + (-t133 * qJD(2) - t129 - t149) * t89 + t117 + t164, -t103 * qJDD(2), -g(1) * t39 - g(2) * t37 - g(3) * t125 + t2 * qJ(4) - t1 * t158 - t12 * t15 + t134 * t14 + t149 * t161 - t9 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t69, t45, -t24 * qJD(3) + t22 * t143 + t160 + t6, 0, 0, 0, 0, 0, 0, t34, t45, -t69, -t139 - t119 + (-t89 * t132 - t162) * qJD(2) + t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 - 0.2e1 * t122, t69 + 0.2e1 * t121, (-t84 - t85) * t92, (t156 + (t9 - t123) * t88) * qJD(2) + t100 + t124;];
tau_reg = t8;
