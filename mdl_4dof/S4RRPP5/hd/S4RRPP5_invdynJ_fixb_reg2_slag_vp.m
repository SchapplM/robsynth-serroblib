% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPP5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:39
% EndTime: 2019-12-31 17:00:41
% DurationCPUTime: 0.90s
% Computational Cost: add. (520->197), mult. (1169->213), div. (0->0), fcn. (577->4), ass. (0->119)
t76 = sin(qJ(2));
t63 = t76 * qJ(3);
t78 = cos(qJ(2));
t126 = t78 * pkin(2) + t63;
t146 = -pkin(1) - t126;
t17 = t146 * qJD(1);
t147 = qJDD(1) * t146;
t140 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t118 = qJD(1) * t78;
t55 = pkin(5) * t118;
t24 = pkin(3) * t118 + t55;
t145 = -qJD(4) - t24;
t128 = pkin(2) + qJ(4);
t144 = t128 * qJD(2);
t143 = t128 * qJDD(2);
t113 = qJD(1) * qJD(2);
t108 = t78 * t113;
t59 = t76 * qJDD(1);
t142 = (-t108 - t59) * pkin(3);
t141 = 0.2e1 * t140;
t114 = qJD(2) * qJ(3);
t16 = t114 - t145;
t77 = sin(qJ(1));
t79 = cos(qJ(1));
t99 = g(1) * t79 + g(2) * t77;
t130 = t76 * t79;
t131 = t76 * t77;
t139 = -g(1) * t130 - g(2) * t131 + g(3) * t78;
t138 = pkin(3) + pkin(5);
t137 = g(1) * t77;
t134 = g(2) * t79;
t133 = g(3) * t76;
t74 = t78 ^ 2;
t81 = qJD(1) ^ 2;
t132 = t74 * t81;
t129 = t76 * t81;
t60 = t78 * qJDD(1);
t127 = qJ(3) * t60 + qJD(3) * t118;
t125 = t79 * pkin(1) + t77 * pkin(5);
t73 = t76 ^ 2;
t124 = t73 + t74;
t123 = qJ(3) * t78;
t122 = qJD(2) * pkin(2);
t121 = t78 * qJ(4);
t120 = pkin(5) * qJDD(2);
t119 = qJD(1) * t76;
t34 = t138 * t78;
t25 = qJD(2) * t34;
t117 = qJDD(2) * pkin(2);
t116 = t76 * qJD(3);
t54 = pkin(5) * t119;
t22 = -pkin(3) * t119 - t54;
t115 = qJD(3) - t22;
t40 = pkin(5) * t108;
t50 = pkin(5) * t59;
t112 = qJDD(3) + t40 + t50;
t52 = pkin(5) * t60;
t111 = pkin(3) * t60 + qJDD(4) + t52;
t110 = t138 * qJD(2);
t109 = t76 * t113;
t102 = t121 + t126;
t18 = -pkin(1) - t102;
t91 = -t128 * t78 - pkin(1) - t63;
t8 = t91 * qJD(1);
t107 = qJD(1) * t18 + t8;
t106 = t126 * t79 + t125;
t105 = t50 + t139;
t23 = t76 * t110;
t103 = t76 * t108;
t101 = t124 * qJDD(1) * pkin(5);
t80 = qJD(2) ^ 2;
t100 = pkin(5) * t80 + t134;
t31 = -t55 - t114;
t97 = (qJD(3) + t54 - t122) * t78 + t31 * t76;
t96 = -qJDD(3) - t105;
t95 = qJ(4) * t76 - t123;
t13 = t112 - t117;
t93 = -0.2e1 * pkin(1) * t113 - t120;
t92 = -t78 * t114 - t116;
t90 = 0.2e1 * qJDD(1) * pkin(1) - t100;
t41 = pkin(2) * t109;
t83 = t95 * qJD(2) - t78 * qJD(4) - t116;
t1 = t83 * qJD(1) + qJDD(1) * t91 + t41;
t57 = t76 * t122;
t5 = t57 + t83;
t89 = -qJD(1) * t5 - qJDD(1) * t18 - t1 - t134;
t88 = -0.2e1 * t17 * qJD(2) + t120;
t87 = t112 - t142 - t143;
t86 = t8 * t118 - t99 * t78 + t111;
t15 = t57 + t92;
t3 = t92 * qJD(1) + t147 + t41;
t85 = qJD(1) * t15 + t100 + t147 + t3;
t9 = pkin(5) * t109 - t140 - t52;
t84 = t97 * qJD(2) + t13 * t76 - t9 * t78;
t68 = t79 * pkin(5);
t61 = t73 * t81;
t58 = pkin(2) * t119;
t46 = t78 * t137;
t45 = g(1) * t131;
t39 = t79 * t123;
t37 = t77 * t123;
t36 = t78 * t129;
t35 = -t61 - t80;
t33 = t138 * t76;
t32 = qJDD(2) + t36;
t30 = -t61 + t132;
t29 = qJDD(2) * t78 - t80 * t76;
t28 = qJDD(2) * t76 + t80 * t78;
t21 = -qJ(3) * t118 + t58;
t20 = t74 * qJDD(1) - 0.2e1 * t103;
t19 = t73 * qJDD(1) + 0.2e1 * t103;
t14 = t95 * qJD(1) + t58;
t12 = t115 - t144;
t11 = t17 * t119;
t10 = t76 * t60 + (-t73 + t74) * t113;
t7 = 0.2e1 * t10;
t4 = -qJD(1) * t23 + t111 + t140;
t2 = -qJD(2) * qJD(4) + t87;
t6 = [0, 0, 0, 0, 0, qJDD(1), -t134 + t137, t99, 0, 0, t19, t7, t28, t20, t29, 0, t76 * t93 + t78 * t90 + t46, -t76 * t90 + t78 * t93 - t45, 0.2e1 * t101 - t99, -g(1) * (-t77 * pkin(1) + t68) - g(2) * t125 + (t124 * pkin(5) ^ 2 + pkin(1) ^ 2) * qJDD(1), 0, -t28, -t29, t19, t7, t20, t101 + t84 - t99, t76 * t88 + t78 * t85 - t46, -t76 * t85 + t78 * t88 + t45, pkin(5) * t84 - g(1) * t68 - g(2) * t106 + t17 * t15 + (-t137 + t3) * t146, 0, -t29, t28, t20, -0.2e1 * t10, t19, (t12 * qJD(2) + qJDD(1) * t34 + t4 + (qJD(2) * t33 - t23) * qJD(1)) * t78 + (-t16 * qJD(2) + qJDD(1) * t33 + t2) * t76 - t99, t34 * qJDD(2) + t45 + (-t107 * t78 - t23) * qJD(2) + t89 * t76, -t33 * qJDD(2) + t46 + (t107 * t76 - t25) * qJD(2) + t89 * t78, t1 * t18 + t8 * t5 + t2 * t33 + t12 * t25 + t4 * t34 - t16 * t23 - g(1) * (t79 * pkin(3) + t68) - g(2) * (t79 * t121 + t106) + (-g(1) * (t146 - t121) - g(2) * pkin(3)) * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t30, t59, t36, t60, qJDD(2), pkin(1) * t129 - t105, t133 - t52 + (pkin(1) * t81 + t99) * t78, 0, 0, qJDD(2), -t59, -t60, -t36, -t30, t36, -pkin(2) * t59 + (-qJD(2) * t126 - t97) * qJD(1) + t127, -t118 * t21 + t11 - 0.2e1 * t117 - t96, t52 + (qJD(1) * t21 - g(3)) * t76 + (qJD(1) * t17 - t99) * t78 + t141, -t9 * qJ(3) - t31 * qJD(3) - t13 * pkin(2) - t17 * t21 - g(1) * (-pkin(2) * t130 + t39) - g(2) * (-pkin(2) * t131 + t37) - g(3) * t126 - t97 * qJD(1) * pkin(5), qJDD(2), -t60, t59, t36, t30, -t36, -t128 * t59 + (-t12 - t22 - t144) * t118 + t127, -t22 * qJD(2) + (-g(3) + (t14 - t110) * qJD(1)) * t76 + t86 + t141, -t40 + (0.2e1 * qJD(4) + t24) * qJD(2) + (t14 * t78 - t76 * t8) * qJD(1) + 0.2e1 * t143 + t96 + t142, -g(1) * t39 - g(2) * t37 - g(3) * t102 + t4 * qJ(3) + t115 * t16 + t145 * t12 - t8 * t14 + (t76 * t99 - t2) * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t32, t35, t31 * qJD(2) + t11 + t13 + t139, 0, 0, 0, 0, 0, 0, t59, t35, -t32, t8 * t119 + (-qJD(4) - t16) * qJD(2) + t87 + t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, qJDD(2) - t36, -t80 - t132, -t133 + (-t138 * t119 + t12) * qJD(2) + t86 + t140;];
tau_reg = t6;
