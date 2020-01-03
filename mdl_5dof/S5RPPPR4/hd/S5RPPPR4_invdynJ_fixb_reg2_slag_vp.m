% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPPR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:20
% EndTime: 2019-12-31 17:45:21
% DurationCPUTime: 0.62s
% Computational Cost: add. (1060->158), mult. (1783->194), div. (0->0), fcn. (1168->12), ass. (0->101)
t78 = sin(pkin(7));
t54 = t78 * pkin(1) + qJ(3);
t72 = qJ(1) + pkin(7);
t65 = sin(t72);
t67 = cos(t72);
t132 = g(1) * t65 - g(2) * t67;
t118 = pkin(1) * qJDD(1);
t77 = sin(pkin(8));
t79 = cos(pkin(8));
t82 = sin(qJ(5));
t84 = cos(qJ(5));
t131 = -t82 * t77 + t84 * t79;
t80 = cos(pkin(7));
t58 = -t80 * pkin(1) - pkin(2);
t130 = t58 * qJDD(1);
t117 = qJD(1) * t77;
t110 = t82 * t117;
t42 = t84 * t77 + t82 * t79;
t91 = -qJD(5) * t110 + t42 * qJDD(1);
t116 = qJD(1) * t79;
t109 = t84 * t116;
t36 = t109 - t110;
t129 = t36 ^ 2;
t128 = t77 * pkin(4);
t51 = -qJ(4) + t58;
t126 = -pkin(6) + t51;
t12 = qJD(5) * t109 + t91;
t34 = t42 * qJD(1);
t38 = t42 * qJD(5);
t125 = -t12 * t131 + t38 * t34;
t124 = t36 * t34;
t60 = t78 * t118;
t120 = qJDD(1) * qJ(3) + t60;
t46 = t54 * qJD(1);
t113 = qJD(1) * qJD(4);
t27 = t51 * qJDD(1) + qJDD(3) - t113;
t21 = t79 * qJDD(2) + t77 * t27;
t41 = t51 * qJD(1) + qJD(3);
t23 = t79 * qJD(2) + t77 * t41;
t69 = t77 ^ 2;
t70 = t79 ^ 2;
t119 = t69 + t70;
t115 = t77 * qJDD(1);
t114 = t79 * qJDD(1);
t74 = qJD(3) * qJD(1);
t85 = cos(qJ(1));
t112 = t85 * pkin(1) + t67 * pkin(2) + t65 * qJ(3);
t111 = -t74 - t120;
t44 = qJD(4) + t46;
t83 = sin(qJ(1));
t108 = -t83 * pkin(1) + t67 * qJ(3);
t20 = -t77 * qJDD(2) + t79 * t27;
t15 = -pkin(6) * t114 + t20;
t16 = -pkin(6) * t115 + t21;
t107 = t84 * t15 - t82 * t16;
t22 = -t77 * qJD(2) + t79 * t41;
t106 = t119 * qJDD(1);
t40 = qJDD(4) - t111;
t105 = g(1) * t67 + g(2) * t65;
t103 = g(1) * t83 - g(2) * t85;
t102 = t84 * t114 - t82 * t115;
t11 = qJD(1) * t38 - t102;
t101 = -t11 * t131 - t36 * t38;
t39 = t131 * qJD(5);
t100 = -t42 * t11 + t36 * t39;
t99 = t42 * t12 + t39 * t34;
t98 = t82 * t15 + t84 * t16;
t18 = -pkin(6) * t116 + t22;
t19 = -pkin(6) * t117 + t23;
t3 = t84 * t18 - t82 * t19;
t4 = t82 * t18 + t84 * t19;
t97 = t20 * t79 + t21 * t77;
t96 = t22 * t79 + t23 * t77;
t32 = t126 * t77;
t33 = t126 * t79;
t9 = t84 * t32 + t82 * t33;
t8 = -t82 * t32 + t84 * t33;
t13 = -t38 * qJD(5) + qJDD(5) * t131;
t14 = -t39 * qJD(5) - t42 * qJDD(5);
t94 = -t97 + t132;
t93 = qJDD(3) + t130;
t92 = t54 * qJDD(1) - t105;
t90 = -t105 + t40;
t1 = t3 * qJD(5) + t98;
t2 = -t4 * qJD(5) + t107;
t89 = t1 * t42 + t131 * t2 - t3 * t38 + t4 * t39 - t132;
t88 = t40 + t92 + t74;
t86 = qJD(1) ^ 2;
t81 = -pkin(6) - qJ(4);
t76 = qJDD(2) - g(3);
t71 = pkin(8) + qJ(5);
t66 = cos(t71);
t64 = sin(t71);
t59 = pkin(4) * t115;
t45 = t54 + t128;
t31 = t34 ^ 2;
t30 = pkin(4) * t117 + t44;
t26 = t59 + t40;
t6 = -qJD(4) * t131 - t9 * qJD(5);
t5 = -t42 * qJD(4) + t8 * qJD(5);
t7 = [0, 0, 0, 0, 0, qJDD(1), t103, g(1) * t85 + g(2) * t83, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t80 * t118 + t132, t105 - 0.2e1 * t60, 0, (t103 + (t78 ^ 2 + t80 ^ 2) * t118) * pkin(1), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) - t132 + 0.2e1 * t130, 0.2e1 * t74 + t92 + t120, -t111 * t54 + t46 * qJD(3) + t93 * t58 - g(1) * (-t65 * pkin(2) + t108) - g(2) * t112, t70 * qJDD(1), -0.2e1 * t77 * t114, 0, t69 * qJDD(1), 0, 0, t88 * t77, t88 * t79, -t51 * t106 + t119 * t113 + t94, t40 * t54 + t44 * qJD(3) - g(1) * ((-pkin(2) - qJ(4)) * t65 + t108) - g(2) * (t67 * qJ(4) + t112) + t97 * t51 - t96 * qJD(4), t101, -t100 + t125, t13, t99, t14, 0, qJD(3) * t34 + t6 * qJD(5) + t8 * qJDD(5) - t105 * t64 + t45 * t12 + t26 * t42 + t30 * t39, qJD(3) * t36 - t5 * qJD(5) - t9 * qJDD(5) - t105 * t66 - t45 * t11 + t131 * t26 - t30 * t38, t8 * t11 - t9 * t12 - t5 * t34 - t6 * t36 - t89, t1 * t9 + t4 * t5 + t2 * t8 + t3 * t6 + t26 * t45 + t30 * qJD(3) - g(1) * (t67 * t128 + (-pkin(2) + t81) * t65 + t108) - g(2) * (t65 * t128 - t67 * t81 + t112); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20 * t77 + t21 * t79 - g(3), 0, 0, 0, 0, 0, 0, t14, -t13, t100 + t125, t1 * t131 - t2 * t42 - t3 * t39 - t4 * t38 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t86, -t46 * qJD(1) - t132 + t93, 0, 0, 0, 0, 0, 0, -t86 * t77, -t86 * t79, -t106, -t44 * qJD(1) - t94, 0, 0, 0, 0, 0, 0, -qJD(1) * t34 + t13, -qJD(1) * t36 + t14, -t101 - t99, -t30 * qJD(1) + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, t114, -t119 * t86, t96 * qJD(1) + t90, 0, 0, 0, 0, 0, 0, (t36 + t109) * qJD(5) + t91, -0.2e1 * t34 * qJD(5) + t102, -t31 - t129, t3 * t36 + t4 * t34 + t59 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, -t31 + t129, t102, -t124, (t36 - t109) * qJD(5) - t91, qJDD(5), g(3) * t64 - t132 * t66 - t30 * t36 + t107, g(3) * t66 + t132 * t64 + t30 * t34 - t98, 0, 0;];
tau_reg = t7;
