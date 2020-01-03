% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:26
% EndTime: 2019-12-31 17:49:28
% DurationCPUTime: 0.74s
% Computational Cost: add. (1104->186), mult. (2197->219), div. (0->0), fcn. (1536->12), ass. (0->110)
t84 = qJ(1) + pkin(7);
t77 = sin(t84);
t79 = cos(t84);
t149 = -g(1) * t77 + g(2) * t79;
t88 = cos(pkin(7));
t142 = t88 * pkin(1);
t71 = -pkin(2) - t142;
t128 = qJDD(1) * t71;
t55 = qJDD(3) + t128;
t150 = -t149 - t55;
t114 = g(1) * t79 + g(2) * t77;
t139 = cos(qJ(4));
t85 = sin(pkin(8));
t87 = cos(pkin(8));
t90 = sin(qJ(4));
t54 = t139 * t85 + t90 * t87;
t45 = t54 * qJD(1);
t86 = sin(pkin(7));
t67 = t86 * pkin(1) + qJ(3);
t140 = pkin(6) + t67;
t49 = t140 * t85;
t50 = t140 * t87;
t17 = t139 * t50 - t90 * t49;
t83 = pkin(8) + qJ(4);
t76 = sin(t83);
t120 = t139 * t87;
t136 = t90 * t85;
t101 = t120 - t136;
t102 = -t139 * t49 - t90 * t50;
t8 = qJD(3) * t101 + t102 * qJD(4);
t148 = -t8 * qJD(4) - t17 * qJDD(4) + t149 * t76;
t118 = qJD(4) * t139;
t51 = qJD(1) * qJD(3) + t67 * qJDD(1);
t73 = t87 * qJDD(2);
t31 = t73 + (-pkin(6) * qJDD(1) - t51) * t85;
t126 = t87 * qJDD(1);
t36 = t85 * qJDD(2) + t87 * t51;
t32 = pkin(6) * t126 + t36;
t133 = pkin(6) * qJD(1);
t60 = t67 * qJD(1);
t75 = t87 * qJD(2);
t33 = t75 + (-t60 - t133) * t85;
t123 = t33 * t118 + t139 * t32 + t90 * t31;
t78 = cos(t83);
t147 = -g(3) * t76 - t114 * t78 + t123;
t146 = t45 ^ 2;
t91 = sin(qJ(1));
t141 = t91 * pkin(1);
t115 = qJD(1) * t120;
t121 = qJD(1) * t136;
t43 = -t115 + t121;
t138 = t45 * t43;
t38 = t85 * qJD(2) + t87 * t60;
t34 = t87 * t133 + t38;
t137 = t90 * t34;
t117 = qJDD(1) * t139;
t127 = t85 * qJDD(1);
t111 = -t87 * t117 + t90 * t127;
t48 = t54 * qJD(4);
t20 = qJD(1) * t48 + t111;
t130 = qJD(4) * t90;
t47 = -t118 * t87 + t85 * t130;
t135 = -t54 * t20 + t47 * t43;
t134 = t85 ^ 2 + t87 ^ 2;
t7 = t139 * t34 + t90 * t33;
t132 = t7 * qJD(4);
t6 = t139 * t33 - t137;
t131 = qJD(5) - t6;
t129 = qJDD(4) * pkin(4);
t124 = qJDD(4) * qJ(5);
t122 = qJD(4) * t115 + t85 * t117 + t90 * t126;
t70 = t87 * pkin(3) + pkin(2);
t119 = t6 + t137;
t116 = t34 * t118 + t33 * t130 - t139 * t31 + t90 * t32;
t92 = cos(qJ(1));
t112 = g(1) * t91 - g(2) * t92;
t110 = t78 * pkin(4) + t76 * qJ(5);
t19 = qJD(4) * t121 - t122;
t108 = t101 * t19 + t48 * t45;
t35 = -t85 * t51 + t73;
t107 = -t35 * t85 + t36 * t87;
t106 = (-t85 * t60 + t75) * t85 - t38 * t87;
t59 = -t70 - t142;
t22 = -t48 * qJD(4) + qJDD(4) * t101;
t103 = t110 + t70;
t100 = -g(3) * t78 + t114 * t76 - t116;
t99 = -t128 + t150;
t41 = qJD(1) * t59 + qJD(3);
t39 = qJDD(1) * t59 + qJDD(3);
t9 = qJD(3) * t54 + t17 * qJD(4);
t98 = -t9 * qJD(4) + t102 * qJDD(4) - t149 * t78;
t12 = t43 * pkin(4) - t45 * qJ(5) + t41;
t97 = t12 * t45 + qJDD(5) - t100;
t96 = t20 * pkin(4) + t19 * qJ(5) + t39;
t95 = 0.2e1 * t45 * qJD(4) + t111;
t89 = -pkin(6) - qJ(3);
t80 = t92 * pkin(1);
t42 = t43 ^ 2;
t21 = -t47 * qJD(4) + t54 * qJDD(4);
t18 = t45 * pkin(4) + t43 * qJ(5);
t15 = -pkin(4) * t101 - t54 * qJ(5) + t59;
t13 = t48 * pkin(4) + t47 * qJ(5) - t54 * qJD(5);
t11 = (t43 - t121) * qJD(4) + t122;
t10 = (t43 + t121) * qJD(4) - t122;
t5 = qJD(4) * qJ(5) + t7;
t4 = -qJD(4) * pkin(4) + t131;
t3 = -t45 * qJD(5) + t96;
t2 = qJDD(5) + t116 - t129;
t1 = t124 + (qJD(5) - t137) * qJD(4) + t123;
t14 = [qJDD(1), t112, g(1) * t92 + g(2) * t91, (t112 + (t86 ^ 2 + t88 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t99 * t87, -t99 * t85, t51 * t134 + t107 - t114, t55 * t71 - g(1) * (-t77 * pkin(2) + t79 * qJ(3) - t141) - g(2) * (t79 * pkin(2) + t77 * qJ(3) + t80) + t107 * t67 - t106 * qJD(3), -t19 * t54 - t45 * t47, -t108 + t135, t21, t22, 0, -t101 * t39 + t59 * t20 + t41 * t48 + t98, -t59 * t19 + t39 * t54 - t41 * t47 + t148, -t101 * t3 + t12 * t48 + t13 * t43 + t15 * t20 + t98, t1 * t101 + t102 * t19 - t17 * t20 + t2 * t54 - t4 * t47 - t8 * t43 + t9 * t45 - t5 * t48 - t114, t12 * t47 - t13 * t45 + t15 * t19 - t3 * t54 - t148, g(1) * t141 - g(2) * t80 + t1 * t17 + t12 * t13 + t3 * t15 - t2 * t102 + t4 * t9 + t5 * t8 + (g(1) * t89 - g(2) * t103) * t79 + (g(1) * t103 + g(2) * t89) * t77; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, t35 * t87 + t36 * t85 - g(3), 0, 0, 0, 0, 0, t22, -t21, t22, t108 + t135, t21, t1 * t54 - t101 * t2 + t4 * t48 - t5 * t47 - g(3); 0, 0, 0, 0, -t126, t127, -t134 * qJD(1) ^ 2, qJD(1) * t106 - t150, 0, 0, 0, 0, 0, t95, -t10, t95, -t42 - t146, t10, t5 * t43 + (-qJD(5) - t4) * t45 + t96 + t149; 0, 0, 0, 0, 0, 0, 0, 0, t138, -t42 + t146, t11, -t111, qJDD(4), -t41 * t45 + t100 + t132, t119 * qJD(4) + t41 * t43 - t147, -t18 * t43 + 0.2e1 * t129 + t132 - t97, pkin(4) * t19 - t20 * qJ(5) + (t5 - t7) * t45 + (t4 - t131) * t43, 0.2e1 * t124 - t12 * t43 + t18 * t45 + (0.2e1 * qJD(5) - t119) * qJD(4) + t147, -t2 * pkin(4) - g(3) * t110 + t1 * qJ(5) - t12 * t18 + t131 * t5 - t4 * t7 + t114 * (pkin(4) * t76 - qJ(5) * t78); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t138, t11, -qJD(4) ^ 2 - t146, -t5 * qJD(4) - t129 + t97;];
tau_reg = t14;
