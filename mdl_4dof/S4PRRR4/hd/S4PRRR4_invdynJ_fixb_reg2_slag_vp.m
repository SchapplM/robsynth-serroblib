% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4PRRR4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4PRRR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:42
% EndTime: 2019-12-31 16:32:43
% DurationCPUTime: 0.69s
% Computational Cost: add. (936->168), mult. (2093->224), div. (0->0), fcn. (1372->8), ass. (0->107)
t111 = qJD(2) * qJD(3);
t79 = cos(qJ(3));
t100 = t79 * t111;
t77 = sin(qJ(3));
t114 = t77 * qJDD(2);
t134 = -t100 - t114;
t70 = pkin(7) + qJ(2);
t62 = sin(t70);
t63 = cos(t70);
t95 = g(1) * t63 + g(2) * t62;
t76 = sin(qJ(4));
t78 = cos(qJ(4));
t35 = t76 * t79 + t77 * t78;
t33 = t35 * qJD(2);
t71 = qJD(3) + qJD(4);
t118 = qJD(4) * t76;
t133 = pkin(6) + pkin(5);
t103 = qJD(2) * t133;
t115 = t77 * qJD(1);
t47 = t133 * t77;
t65 = t79 * qJDD(1);
t18 = qJDD(3) * pkin(3) + t65 - qJDD(2) * t47 + (-t103 * t79 - t115) * qJD(3);
t101 = t77 * t111;
t112 = qJD(1) * qJD(3);
t113 = t79 * qJDD(2);
t109 = pkin(5) * t113 + qJDD(1) * t77 + t112 * t79;
t24 = -pkin(5) * t101 + t109;
t88 = t101 - t113;
t19 = -pkin(6) * t88 + t24;
t122 = qJD(3) * pkin(3);
t66 = t79 * qJD(1);
t28 = -t103 * t77 + t66;
t27 = t28 + t122;
t48 = t133 * t79;
t29 = qJD(2) * t48 + t115;
t1 = (qJD(4) * t27 + t19) * t78 - t29 * t118 + t76 * t18;
t130 = g(3) * t79;
t21 = t71 * t35;
t93 = -t113 * t78 + t114 * t76;
t10 = qJD(2) * t21 + t93;
t117 = qJD(4) * t78;
t124 = t78 * t79;
t126 = t76 * t77;
t91 = t71 * t126;
t20 = -qJD(3) * t124 - t117 * t79 + t91;
t119 = qJD(2) * t79;
t104 = t78 * t119;
t120 = qJD(2) * t77;
t105 = t76 * t120;
t31 = -t104 + t105;
t129 = -t10 * t35 + t20 * t31;
t128 = t33 * t31;
t127 = t76 * t29;
t125 = t78 * t29;
t72 = t77 ^ 2;
t73 = t79 ^ 2;
t123 = t72 - t73;
t121 = pkin(5) * qJDD(2);
t116 = qJDD(2) * pkin(2);
t82 = qJD(2) ^ 2;
t110 = t77 * t82 * t79;
t108 = t77 * t122;
t107 = pkin(5) * t120;
t106 = pkin(5) * t119;
t60 = pkin(3) * t79 + pkin(2);
t102 = qJD(3) * t133;
t98 = -qJD(4) * t104 - t76 * t113 + t134 * t78;
t96 = t77 * t100;
t94 = g(1) * t62 - g(2) * t63;
t34 = -t124 + t126;
t9 = qJD(2) * t91 + t98;
t92 = -t21 * t33 + t34 * t9;
t69 = qJDD(3) + qJDD(4);
t90 = t20 * t71 - t35 * t69;
t15 = t27 * t76 + t125;
t22 = -t47 * t78 - t48 * t76;
t23 = -t47 * t76 + t48 * t78;
t89 = -0.2e1 * pkin(2) * t111 - pkin(5) * qJDD(3);
t87 = pkin(2) * t82 + t95;
t81 = qJD(3) ^ 2;
t86 = -pkin(5) * t81 + 0.2e1 * t116 + t94;
t2 = -qJD(4) * t15 + t18 * t78 - t76 * t19;
t46 = t60 * qJD(2);
t75 = qJ(3) + qJ(4);
t67 = sin(t75);
t68 = cos(t75);
t85 = g(3) * t67 - t46 * t31 + t68 * t95 - t1;
t25 = pkin(5) * t134 - t77 * t112 + t65;
t38 = t66 - t107;
t39 = t106 + t115;
t84 = t24 * t79 - t25 * t77 + (-t38 * t79 - t39 * t77) * qJD(3) - t95;
t83 = -g(3) * t68 + t46 * t33 + t67 * t95 + t2;
t74 = qJDD(1) - g(3);
t45 = qJDD(3) * t79 - t77 * t81;
t44 = qJDD(3) * t77 + t79 * t81;
t37 = t79 * t102;
t36 = t77 * t102;
t30 = pkin(3) * t88 - t116;
t17 = t28 * t78 - t127;
t16 = -t28 * t76 - t125;
t14 = t27 * t78 - t127;
t11 = -t31 ^ 2 + t33 ^ 2;
t8 = -t21 * t71 - t34 * t69;
t7 = -qJD(4) * t23 + t76 * t36 - t78 * t37;
t6 = qJD(4) * t22 - t78 * t36 - t76 * t37;
t3 = -t98 + (-t105 + t31) * t71;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, 0, 0, 0, 0, 0, t45, -t44, 0, t24 * t77 + t25 * t79 - g(3) + (-t38 * t77 + t39 * t79) * qJD(3), 0, 0, 0, 0, 0, 0, t8, t90, -t92 + t129, t1 * t35 - t14 * t21 - t15 * t20 - t2 * t34 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t94, t95, 0, 0, qJDD(2) * t72 + 0.2e1 * t96, -0.2e1 * t111 * t123 + 0.2e1 * t113 * t77, t44, qJDD(2) * t73 - 0.2e1 * t96, t45, 0, t77 * t89 + t79 * t86, -t77 * t86 + t79 * t89, (t72 + t73) * t121 + t84, (t94 + t116) * pkin(2) + t84 * pkin(5), -t20 * t33 - t35 * t9, t92 + t129, -t90, t10 * t34 + t21 * t31, t8, 0, -t60 * t10 + t108 * t31 - t46 * t21 + t22 * t69 + t30 * t34 + t68 * t94 + t7 * t71, t108 * t33 + t46 * t20 - t23 * t69 + t30 * t35 - t6 * t71 + t60 * t9 - t67 * t94, -t1 * t34 - t10 * t23 + t14 * t20 - t15 * t21 - t2 * t35 + t22 * t9 - t31 * t6 - t33 * t7 - t95, t1 * t23 + t15 * t6 + t2 * t22 + t14 * t7 - t30 * t60 - t46 * t108 - g(1) * (t133 * t63 - t60 * t62) - g(2) * (t133 * t62 + t60 * t63); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, t123 * t82, t114, t110, t113, qJDD(3), -t130 + t65 + (t39 - t106) * qJD(3) + (-t112 + t87 - t121) * t77, g(3) * t77 + (t38 + t107) * qJD(3) + t87 * t79 - t109, 0, 0, t128, t11, t3, -t128, -t93, t69, -t16 * t71 + (-t118 * t71 - t120 * t31 + t69 * t78) * pkin(3) + t83, t17 * t71 + (-t117 * t71 - t120 * t33 - t69 * t76) * pkin(3) + t85, (t15 + t16) * t33 + (-t14 + t17) * t31 + (-t10 * t76 + t78 * t9 + (-t31 * t78 + t33 * t76) * qJD(4)) * pkin(3), -t14 * t16 - t15 * t17 + (-t130 + t1 * t76 + t2 * t78 + (-t14 * t76 + t15 * t78) * qJD(4) + (qJD(2) * t46 + t95) * t77) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t11, t3, -t128, -t93, t69, t15 * t71 + t83, t14 * t71 + t85, 0, 0;];
tau_reg = t4;
