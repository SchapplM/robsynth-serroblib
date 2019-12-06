% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PPRPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:28
% EndTime: 2019-12-05 15:05:30
% DurationCPUTime: 0.84s
% Computational Cost: add. (975->171), mult. (2103->258), div. (0->0), fcn. (1775->12), ass. (0->111)
t68 = cos(pkin(8));
t48 = -t68 * qJDD(1) + qJDD(4);
t65 = sin(pkin(8));
t66 = sin(pkin(7));
t69 = cos(pkin(7));
t92 = g(1) * t69 + g(2) * t66;
t84 = g(3) * t68 - t92 * t65;
t130 = t48 + t84;
t113 = qJD(1) * t65;
t71 = sin(qJ(3));
t73 = cos(qJ(3));
t129 = t73 * qJD(2) - t71 * t113;
t64 = sin(pkin(9));
t67 = cos(pkin(9));
t39 = t64 * t71 - t67 * t73;
t128 = g(3) * t65;
t127 = t68 ^ 2 * qJDD(1) - g(3);
t42 = t71 * qJD(2) + t73 * t113;
t126 = t64 * t42;
t124 = t66 * t68;
t123 = t66 * t71;
t122 = t66 * t73;
t34 = t67 * t42;
t61 = qJ(3) + pkin(9);
t55 = sin(t61);
t120 = t69 * t55;
t56 = cos(t61);
t119 = t69 * t56;
t118 = t69 * t71;
t117 = t69 * t73;
t70 = sin(qJ(5));
t72 = cos(qJ(5));
t116 = t70 * t72;
t109 = qJDD(3) * pkin(3);
t105 = t65 * qJDD(1);
t57 = t73 * qJDD(2);
t23 = -t42 * qJD(3) - t71 * t105 + t57;
t20 = t23 + t109;
t103 = t71 * qJDD(2);
t24 = qJD(3) * t129 + t73 * t105 + t103;
t6 = t64 * t20 + t67 * t24;
t36 = qJD(3) * pkin(3) + t129;
t14 = t64 * t36 + t34;
t62 = t70 ^ 2;
t63 = t72 ^ 2;
t115 = t62 - t63;
t114 = t62 + t63;
t112 = qJD(3) * t39;
t111 = qJD(3) * t65;
t40 = t64 * t73 + t67 * t71;
t32 = t40 * t65;
t110 = qJD(5) * t32;
t18 = t129 * t67 - t126;
t108 = t18 * qJD(3);
t107 = t112 * qJD(3);
t106 = qJDD(1) - g(3);
t104 = t70 * qJDD(3);
t102 = t72 * qJDD(3);
t101 = qJD(3) * qJD(5);
t100 = t68 * t118;
t75 = qJD(3) ^ 2;
t99 = t75 * t116;
t97 = -g(1) * t66 + g(2) * t69;
t96 = -t67 * t20 + t64 * t24;
t95 = t65 * t106;
t94 = qJDD(3) * t114;
t93 = t101 * t116;
t12 = qJD(3) * pkin(6) + t14;
t49 = -t68 * qJD(1) + qJD(4);
t10 = t72 * t12 + t70 * t49;
t9 = -t70 * t12 + t72 * t49;
t91 = t10 * t72 - t9 * t70;
t33 = t39 * t65;
t22 = -t72 * t33 - t68 * t70;
t21 = t70 * t33 - t68 * t72;
t13 = t67 * t36 - t126;
t90 = t73 * qJDD(3) - t75 * t71;
t89 = -qJDD(3) * t71 - t75 * t73;
t37 = t40 * qJD(3);
t88 = -t37 * qJD(3) - t39 * qJDD(3);
t87 = -t68 * t123 - t117;
t28 = t56 * t124 - t120;
t30 = t68 * t119 + t66 * t55;
t83 = g(1) * t30 + g(2) * t28 + t56 * t128;
t74 = qJD(5) ^ 2;
t82 = t40 * t74 - t88;
t81 = 0.2e1 * t112 * qJD(5) - qJDD(5) * t40;
t11 = -qJD(3) * pkin(4) - t13;
t52 = t64 * pkin(3) + pkin(6);
t53 = -t67 * pkin(3) - pkin(4);
t80 = -qJDD(5) * t52 + (qJD(3) * t53 + t11 + t18) * qJD(5);
t17 = t129 * t64 + t34;
t27 = -t55 * t124 - t119;
t29 = -t68 * t120 + t66 * t56;
t79 = -g(1) * t29 - g(2) * t27 + t17 * qJD(3) + t55 * t128;
t4 = qJDD(3) * pkin(6) + t6;
t78 = -t11 * qJD(3) - t4 + t83;
t1 = t9 * qJD(5) + t72 * t4 + t70 * t48;
t43 = t72 * t48;
t2 = -t10 * qJD(5) - t70 * t4 + t43;
t77 = t1 * t72 - t2 * t70 + (-t10 * t70 - t72 * t9) * qJD(5);
t3 = -qJDD(3) * pkin(4) + t96;
t76 = -qJDD(3) * t53 - t52 * t74 - t3 + t79;
t50 = pkin(3) * t122;
t46 = qJDD(5) * t72 - t74 * t70;
t45 = qJDD(5) * t70 + t74 * t72;
t26 = t40 * t111;
t25 = t39 * t111;
t8 = -t22 * qJD(5) + t70 * t26;
t7 = t21 * qJD(5) - t72 * t26;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t106, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65 ^ 2 * qJDD(1) + t127, 0, 0, 0, 0, 0, 0, t89 * t65, -t90 * t65, 0, (-t23 * t71 + t24 * t73 + (-t129 * t73 - t42 * t71) * qJD(3)) * t65 + t127, 0, 0, 0, 0, 0, 0, t25 * qJD(3) - t32 * qJDD(3), t26 * qJD(3) + t33 * qJDD(3), 0, t13 * t25 - t14 * t26 + t32 * t96 - t6 * t33 - t48 * t68 - g(3), 0, 0, 0, 0, 0, 0, -t32 * t102 + t8 * qJD(5) + t21 * qJDD(5) + (t70 * t110 + t25 * t72) * qJD(3), t32 * t104 - t7 * qJD(5) - t22 * qJDD(5) + (t72 * t110 - t25 * t70) * qJD(3), (-t21 * t70 + t22 * t72) * qJDD(3) + (t7 * t72 - t70 * t8 + (-t21 * t72 - t22 * t70) * qJD(5)) * qJD(3), t1 * t22 + t10 * t7 - t11 * t25 + t2 * t21 + t3 * t32 + t9 * t8 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) + t97, 0, 0, 0, 0, 0, 0, t90, t89, 0, t23 * t73 + t24 * t71 + (-t129 * t71 + t42 * t73) * qJD(3) + t97, 0, 0, 0, 0, 0, 0, t88, -t40 * qJDD(3) + t107, 0, -t112 * t14 - t13 * t37 + t39 * t96 + t6 * t40 + t97, 0, 0, 0, 0, 0, 0, t81 * t70 - t82 * t72, t82 * t70 + t81 * t72, -t114 * t107 + t40 * t94, t11 * t37 - t112 * t91 + t3 * t39 + t77 * t40 + t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t57 - g(1) * (-t100 + t122) - g(2) * t87 - t71 * t95, -t103 - g(1) * (-t68 * t117 - t123) - g(2) * (-t68 * t122 + t118) - t73 * t95, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t67 * t109 + t79 - t96, -t64 * t109 + t108 - t6 + t83, 0, -g(1) * t50 + t13 * t17 - t14 * t18 + (g(2) * t117 - t96 * t67 + t6 * t64 + (t92 * t68 + t128) * t71) * pkin(3), t62 * qJDD(3) + 0.2e1 * t93, -0.2e1 * t115 * t101 + 0.2e1 * t70 * t102, t45, t63 * qJDD(3) - 0.2e1 * t93, t46, 0, t80 * t70 + t76 * t72, -t76 * t70 + t80 * t72, -t114 * t108 + t52 * t94 + t77 - t83, t3 * t53 - t11 * t17 - g(1) * (-pkin(3) * t100 + t29 * pkin(4) + t30 * pkin(6) + t50) - g(2) * (t87 * pkin(3) + t27 * pkin(4) + t28 * pkin(6)) - t91 * t18 - (-pkin(3) * t71 - pkin(4) * t55 + pkin(6) * t56) * t128 + t77 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, 0, 0, 0, 0, 0, 0, t46, -t45, 0, t91 * qJD(5) + t1 * t70 + t2 * t72 + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, t115 * t75, t104, t99, t102, qJDD(5), t78 * t70 + t84 * t72 + t43, -t130 * t70 + t78 * t72, 0, 0;];
tau_reg = t5;
