% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRP4
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:25
% EndTime: 2019-12-31 17:52:26
% DurationCPUTime: 0.95s
% Computational Cost: add. (1121->224), mult. (1934->256), div. (0->0), fcn. (1085->6), ass. (0->128)
t133 = qJ(2) * qJDD(1);
t95 = -pkin(1) - pkin(2);
t51 = qJDD(1) * t95 + qJDD(2);
t89 = sin(pkin(7));
t90 = cos(pkin(7));
t154 = t133 * t90 + t51 * t89;
t131 = qJD(1) * qJD(2);
t58 = t90 * t131;
t22 = t58 + t154;
t18 = -qJDD(1) * pkin(6) + t22;
t169 = qJD(3) * qJD(4) + t18;
t162 = sin(qJ(1));
t163 = cos(qJ(1));
t39 = -t162 * t89 - t163 * t90;
t40 = -t162 * t90 + t163 * t89;
t113 = -g(1) * t39 - g(2) * t40;
t94 = cos(qJ(4));
t134 = t94 * qJDD(1);
t168 = pkin(4) * t134 + qJDD(5);
t130 = qJD(1) * qJD(4);
t124 = t94 * t130;
t129 = qJD(1) * qJD(5);
t132 = qJ(5) * qJDD(1);
t93 = sin(qJ(4));
t167 = qJ(5) * t124 + (t129 + t132) * t93;
t139 = qJ(2) * qJD(1);
t52 = qJD(1) * t95 + qJD(2);
t29 = t139 * t90 + t52 * t89;
t25 = -qJD(1) * pkin(6) + t29;
t15 = qJD(3) * t94 - t25 * t93;
t157 = t94 * t25;
t16 = qJD(3) * t93 + t157;
t141 = qJD(4) * t93;
t3 = t93 * qJDD(3) - t25 * t141 + t169 * t94;
t75 = t94 * qJDD(3);
t4 = -t16 * qJD(4) - t93 * t18 + t75;
t98 = -(t15 * t94 + t16 * t93) * qJD(4) + t3 * t94 - t4 * t93;
t83 = g(3) * t94;
t138 = qJ(5) * qJD(1);
t11 = t138 * t93 + t15;
t148 = qJD(4) * pkin(4);
t8 = t11 + t148;
t164 = -t11 + t8;
t12 = -t138 * t94 + t16;
t161 = t12 * t94;
t160 = t39 * t93;
t159 = t40 * t93;
t97 = qJD(1) ^ 2;
t158 = t90 * t97;
t156 = t94 * t97;
t155 = g(1) * t159 - g(2) * t160;
t45 = qJ(2) * t90 + t89 * t95;
t153 = t163 * pkin(1) + qJ(2) * t162;
t152 = g(1) * t162 - g(2) * t163;
t86 = t93 ^ 2;
t87 = t94 ^ 2;
t151 = t86 - t87;
t150 = t86 + t87;
t96 = qJD(4) ^ 2;
t149 = t96 + t97;
t42 = -pkin(6) + t45;
t147 = qJ(5) - t42;
t146 = pkin(1) * qJDD(1);
t143 = qJD(1) * t94;
t28 = -t139 * t89 + t52 * t90;
t24 = qJD(1) * pkin(3) - t28;
t19 = pkin(4) * t143 + qJD(5) + t24;
t145 = qJD(1) * t19;
t144 = qJD(1) * t24;
t142 = qJD(2) * t90;
t140 = qJDD(4) * pkin(4);
t137 = qJDD(4) * t93;
t136 = qJDD(4) * t94;
t135 = t93 * qJDD(1);
t127 = t163 * pkin(2) + t153;
t126 = 0.2e1 * t131;
t56 = t89 * t131;
t125 = t93 * t130;
t123 = -t133 * t89 + t90 * t51;
t44 = -qJ(2) * t89 + t90 * t95;
t122 = -g(1) * t160 - g(2) * t159 + t75 + t83;
t41 = pkin(3) - t44;
t80 = t94 * pkin(4);
t36 = t41 + t80;
t121 = -qJD(1) * t36 - t19;
t120 = qJD(4) * t147;
t119 = qJDD(1) * t150;
t118 = -qJD(5) + t142;
t117 = qJDD(2) - t146;
t116 = t93 * t124;
t21 = t123 - t56;
t115 = -pkin(1) * t162 + qJ(2) * t163;
t114 = -g(1) * t40 + g(2) * t39;
t111 = -t8 * t93 + t161;
t109 = t15 * t93 - t16 * t94;
t108 = t28 * t89 - t29 * t90;
t88 = qJDD(1) * pkin(3);
t17 = -t21 + t88;
t47 = -pkin(4) * t141 + qJD(2) * t89;
t5 = -pkin(4) * t125 + t168 + t17;
t107 = -qJD(1) * t47 - qJDD(1) * t36 - t5;
t106 = 0.2e1 * t124 + t135;
t105 = 0.2e1 * t125 - t134;
t104 = g(1) * t163 + g(2) * t162;
t103 = -g(3) * t93 + t113 * t94 - t3;
t102 = t114 - t123;
t101 = -pkin(2) * t162 + t115;
t100 = -qJDD(1) * t41 + t42 * t96 - t17 - t56;
t99 = -qJDD(4) * t42 + (-qJD(1) * t41 - t142 - t24) * qJD(4);
t92 = -qJ(5) - pkin(6);
t70 = t80 + pkin(3);
t54 = t93 * t156;
t50 = t151 * t97;
t49 = -t93 * t96 + t136;
t48 = -t94 * t96 - t137;
t38 = qJDD(1) * t87 - 0.2e1 * t116;
t37 = qJDD(1) * t86 + 0.2e1 * t116;
t27 = t147 * t94;
t26 = t147 * t93;
t23 = -0.2e1 * t130 * t151 + 0.2e1 * t134 * t93;
t14 = -t119 * t89 + t150 * t158;
t10 = -t118 * t93 + t120 * t94;
t9 = t118 * t94 + t120 * t93;
t7 = t106 * t90 + (t149 * t93 - t136) * t89;
t6 = t105 * t90 + (-t149 * t94 - t137) * t89;
t2 = -t94 * t129 + (t125 - t134) * qJ(5) + t3;
t1 = t4 + t140 + t167;
t13 = [0, 0, 0, 0, 0, qJDD(1), t152, t104, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -qJDD(2) + 0.2e1 * t146 + t152, 0, -t104 + t126 + 0.2e1 * t133, -t117 * pkin(1) - g(1) * t115 - g(2) * t153 + (t126 + t133) * qJ(2), 0, 0, 0, 0, 0, qJDD(1), -qJDD(1) * t44 + t102 + 0.2e1 * t56, qJDD(1) * t45 - t113 + t154 + 0.2e1 * t58, 0, -g(1) * t101 - g(2) * t127 - qJD(2) * t108 + t21 * t44 + t22 * t45, t37, t23, t48, t38, -t49, 0, t99 * t93 + (-t100 + t114) * t94, t100 * t93 + t94 * t99 + t155, -t119 * t42 - t150 * t58 + t113 - t98, t17 * t41 - g(1) * (pkin(3) * t40 + pkin(6) * t39 + t101) - g(2) * (-pkin(3) * t39 + pkin(6) * t40 + t127) + (-t109 * t90 + t24 * t89) * qJD(2) + t98 * t42, t37, t23, t48, t38, -t49, 0, t26 * qJDD(4) + (t121 * t93 + t10) * qJD(4) + (-t107 + t114) * t94, t27 * qJDD(4) + t107 * t93 + (t121 * t94 - t9) * qJD(4) + t155, (qJD(4) * t8 + qJDD(1) * t27 - t2 + (qJD(4) * t26 - t9) * qJD(1)) * t94 + (t12 * qJD(4) + qJDD(1) * t26 + t1 + (-qJD(4) * t27 + t10) * qJD(1)) * t93 + t113, -t2 * t27 + t12 * t9 + t1 * t26 + t8 * t10 + t5 * t36 + t19 * t47 - g(1) * (-t39 * t92 + t40 * t70 + t101) - g(2) * (-t39 * t70 - t40 * t92 + t127); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t97, -qJ(2) * t97 + t117 - t152, 0, 0, 0, 0, 0, 0, -qJDD(1) * t90 - t89 * t97, qJDD(1) * t89 - t158, 0, qJD(1) * t108 + t21 * t90 + t22 * t89 - t152, 0, 0, 0, 0, 0, 0, t6, t7, t14, (qJD(1) * t109 - t17) * t90 + (t98 - t144) * t89 - t152, 0, 0, 0, 0, 0, 0, t6, t7, t14, (-qJD(1) * t111 - t5) * t90 + (-t145 - t1 * t93 + t2 * t94 + (-t12 * t93 - t8 * t94) * qJD(4)) * t89 - t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) + g(3), 0, 0, 0, 0, 0, 0, t49, t48, 0, -qJD(4) * t109 + t3 * t93 + t4 * t94 + g(3), 0, 0, 0, 0, 0, 0, t49, t48, 0, qJD(4) * t111 + t1 * t94 + t2 * t93 + g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t50, -t135, t54, -t134, qJDD(4), (-t18 + t144) * t93 + t122, qJD(4) * t15 + t143 * t24 + t103, 0, 0, -t54, t50, -t135, t54, -t134, qJDD(4), 0.2e1 * t140 + (t12 - t157) * qJD(4) + (pkin(4) * t156 + t145 - t169) * t93 + t122 + t167, -t86 * t97 * pkin(4) + t94 * t132 + t11 * qJD(4) + (-qJ(5) * t141 + (qJD(5) + t19) * t94) * qJD(1) + t103, pkin(4) * t135 + (t148 - t164) * t143, t164 * t12 + (t83 + t1 + (t113 + t145) * t93) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, -t106, -t150 * t97, t56 + t88 + (t161 + (-t8 - t148) * t93) * qJD(1) + t102 + t168;];
tau_reg = t13;
