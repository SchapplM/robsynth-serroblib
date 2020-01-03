% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:15
% EndTime: 2019-12-31 18:09:18
% DurationCPUTime: 0.82s
% Computational Cost: add. (1362->204), mult. (2801->259), div. (0->0), fcn. (1786->12), ass. (0->111)
t151 = 2 * qJD(3);
t86 = sin(pkin(7));
t68 = t86 * pkin(1) + pkin(6);
t135 = qJ(4) + t68;
t85 = sin(pkin(8));
t87 = cos(pkin(8));
t90 = sin(qJ(3));
t92 = cos(qJ(3));
t52 = t85 * t92 + t87 * t90;
t47 = t52 * qJD(1);
t42 = t47 ^ 2;
t139 = t87 * t92;
t124 = qJD(1) * t139;
t134 = qJD(1) * t90;
t44 = t134 * t85 - t124;
t150 = -t44 ^ 2 - t42;
t88 = cos(pkin(7));
t70 = -t88 * pkin(1) - pkin(2);
t78 = t92 * pkin(3);
t149 = t70 - t78;
t81 = qJ(1) + pkin(7);
t73 = sin(t81);
t75 = cos(t81);
t148 = -g(1) * t73 + g(2) * t75;
t80 = qJ(3) + pkin(8);
t72 = sin(t80);
t74 = cos(t80);
t147 = t74 * pkin(4) + t72 * qJ(5);
t116 = g(1) * t75 + g(2) * t73;
t56 = t68 * qJDD(1);
t101 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + (qJD(3) * qJD(2)) + t56;
t118 = t135 * qJD(1);
t109 = t118 * qJD(3);
t76 = t92 * qJDD(2);
t14 = qJDD(3) * pkin(3) - t101 * t90 - t109 * t92 + t76;
t17 = (qJDD(2) - t109) * t90 + t101 * t92;
t3 = t87 * t14 - t85 * t17;
t123 = -qJDD(5) + t3;
t43 = qJD(1) * t149 + qJD(4);
t23 = t44 * pkin(4) - t47 * qJ(5) + t43;
t146 = -g(3) * t74 + t116 * t72 - t23 * t47 + t123;
t127 = qJD(1) * qJD(3);
t122 = t90 * t127;
t103 = t52 * qJDD(1) - t85 * t122;
t142 = g(3) * t92;
t4 = t85 * t14 + t87 * t17;
t38 = t90 * qJD(2) + t118 * t92;
t140 = t85 * t38;
t33 = t87 * t38;
t37 = t92 * qJD(2) - t118 * t90;
t35 = qJD(3) * pkin(3) + t37;
t16 = t85 * t35 + t33;
t71 = t78 + pkin(2);
t93 = cos(qJ(1));
t138 = t93 * pkin(1) + t75 * t71;
t83 = t90 ^ 2;
t137 = -t92 ^ 2 + t83;
t59 = qJD(1) * t70;
t133 = qJDD(3) * pkin(4);
t132 = t90 * qJD(3);
t20 = t87 * t37 - t140;
t131 = qJD(5) - t20;
t130 = qJDD(2) - g(3);
t129 = t90 * qJDD(1);
t128 = t92 * qJDD(1);
t126 = qJDD(3) * qJ(5) + t4;
t125 = pkin(3) * t132;
t121 = t92 * t127;
t120 = t135 * t90;
t119 = qJD(3) * t135;
t91 = sin(qJ(1));
t114 = g(1) * t91 - g(2) * t93;
t89 = -qJ(4) - pkin(6);
t113 = -t91 * pkin(1) - t75 * t89;
t112 = -t87 * t128 + t129 * t85;
t15 = t87 * t35 - t140;
t46 = t52 * qJD(3);
t29 = qJD(1) * t46 + t112;
t30 = t121 * t87 + t103;
t49 = qJD(3) * t139 - t132 * t85;
t51 = t85 * t90 - t139;
t108 = -t52 * t29 + t51 * t30 - t49 * t44 + t46 * t47;
t105 = -t90 * qJD(4) - t119 * t92;
t104 = -t59 * qJD(1) + t116 - t56;
t102 = -qJDD(3) * t68 + t59 * t151;
t100 = pkin(3) * t122 + qJDD(1) * t149 + qJDD(4);
t94 = qJD(3) ^ 2;
t99 = -0.2e1 * qJDD(1) * t70 - t68 * t94 - t148;
t39 = t92 * qJD(4) - t119 * t90;
t21 = -t105 * t87 + t85 * t39;
t22 = t105 * t85 + t87 * t39;
t50 = t135 * t92;
t27 = t120 * t87 + t85 * t50;
t28 = -t120 * t85 + t87 * t50;
t98 = t21 * t47 - t22 * t44 + t27 * t30 - t28 * t29 - t116;
t97 = t29 * pkin(4) - t30 * qJ(5) + t100;
t95 = qJD(1) ^ 2;
t69 = -t87 * pkin(3) - pkin(4);
t65 = t85 * pkin(3) + qJ(5);
t55 = qJDD(3) * t92 - t94 * t90;
t54 = qJDD(3) * t90 + t94 * t92;
t26 = t51 * pkin(4) - t52 * qJ(5) + t149;
t25 = pkin(3) * t134 + t47 * pkin(4) + t44 * qJ(5);
t19 = t85 * t37 + t33;
t18 = t46 * pkin(4) - t49 * qJ(5) - t52 * qJD(5) + t125;
t11 = qJD(3) * qJ(5) + t16;
t10 = -qJD(3) * pkin(4) + qJD(5) - t15;
t5 = -t47 * qJD(5) + t97;
t2 = -t123 - t133;
t1 = qJD(3) * qJD(5) + t126;
t6 = [qJDD(1), t114, g(1) * t93 + g(2) * t91, (t114 + (t86 ^ 2 + t88 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t83 * qJDD(1) + 0.2e1 * t121 * t90, -0.2e1 * t127 * t137 + 0.2e1 * t128 * t90, t54, t55, 0, t102 * t90 + t92 * t99, t102 * t92 - t90 * t99, -t15 * t49 - t16 * t46 - t3 * t52 - t4 * t51 + t98, t4 * t28 + t16 * t22 - t3 * t27 - t15 * t21 + t100 * t149 + t43 * t125 - g(1) * (-t73 * t71 + t113) - g(2) * (-t73 * t89 + t138), -t21 * qJD(3) - t27 * qJDD(3) - t148 * t74 + t18 * t44 + t23 * t46 + t26 * t29 + t5 * t51, -t1 * t51 + t10 * t49 - t11 * t46 + t2 * t52 + t98, t22 * qJD(3) + t28 * qJDD(3) - t148 * t72 - t18 * t47 - t23 * t49 - t26 * t30 - t5 * t52, t1 * t28 + t11 * t22 + t5 * t26 + t23 * t18 + t2 * t27 + t10 * t21 - g(1) * t113 - g(2) * (t147 * t75 + t138) + (-g(1) * (-t147 - t71) + g(2) * t89) * t73; 0, 0, 0, t130, 0, 0, 0, 0, 0, t55, -t54, t108, -t15 * t46 + t16 * t49 - t3 * t51 + t4 * t52 - g(3), -t46 * qJD(3) - t51 * qJDD(3), t108, t49 * qJD(3) + t52 * qJDD(3), t1 * t52 + t10 * t46 + t11 * t49 + t2 * t51 - g(3); 0, 0, 0, 0, -t90 * t95 * t92, t137 * t95, t129, t128, qJDD(3), t104 * t90 - t142 + t76, t104 * t92 - t130 * t90, (t16 - t19) * t47 + (-t15 + t20) * t44 + (-t29 * t85 - t30 * t87) * pkin(3), t15 * t19 - t16 * t20 + (-t142 + t3 * t87 + t4 * t85 + (-qJD(1) * t43 + t116) * t90) * pkin(3), t19 * qJD(3) - t25 * t44 + (pkin(4) - t69) * qJDD(3) + t146, -t65 * t29 + t69 * t30 + (t11 - t19) * t47 + (t10 - t131) * t44, -g(3) * t72 + t65 * qJDD(3) - t23 * t44 + t25 * t47 - t116 * t74 + (0.2e1 * qJD(5) - t20) * qJD(3) + t126, t1 * t65 + t2 * t69 - t23 * t25 - t10 * t19 - g(3) * (t147 + t78) + t131 * t11 + t116 * (pkin(3) * t90 + pkin(4) * t72 - qJ(5) * t74); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, t15 * t47 + t16 * t44 + t100 + t148, t47 * t151 + t112, t150, (t44 - t124) * qJD(3) - t103, t11 * t44 + (-qJD(5) - t10) * t47 + t97 + t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t44 - qJDD(3), (t44 + t124) * qJD(3) + t103, -t42 - t94, -t11 * qJD(3) - t133 - t146;];
tau_reg = t6;
