% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PPRRR5
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:49
% EndTime: 2019-12-31 17:35:50
% DurationCPUTime: 0.79s
% Computational Cost: add. (1014->151), mult. (1805->201), div. (0->0), fcn. (1315->10), ass. (0->105)
t56 = qJDD(3) + qJDD(4);
t129 = t56 * pkin(4);
t64 = sin(qJ(3));
t113 = qJD(2) * t64;
t100 = qJD(4) * t113;
t106 = qJD(2) * qJD(3);
t67 = cos(qJ(3));
t52 = t67 * qJDD(2);
t32 = qJDD(3) * pkin(3) - t64 * t106 + t52;
t66 = cos(qJ(4));
t27 = t66 * t32;
t42 = qJD(3) * pkin(3) + t67 * qJD(2);
t63 = sin(qJ(4));
t108 = t64 * qJDD(2);
t82 = -t67 * t106 - t108;
t8 = -(qJD(4) * t42 - t82) * t63 - t66 * t100 + t27;
t6 = -t129 - t8;
t114 = sin(pkin(8));
t61 = cos(pkin(8));
t109 = qJ(3) + qJ(4);
t98 = sin(t109);
t99 = cos(t109);
t28 = -t114 * t98 - t61 * t99;
t29 = -t114 * t99 + t61 * t98;
t90 = g(1) * t29 - g(2) * t28;
t139 = t6 - t90;
t21 = t66 * t113 + t63 * t42;
t57 = qJD(3) + qJD(4);
t17 = t57 * pkin(7) + t21;
t62 = sin(qJ(5));
t65 = cos(qJ(5));
t9 = -t65 * qJD(1) - t62 * t17;
t115 = t9 * qJD(5);
t112 = qJD(4) * t66;
t91 = -t42 * t112 + t82 * t66 + (t100 - t32) * t63;
t5 = t56 * pkin(7) - t91;
t2 = -t62 * qJDD(1) + t65 * t5 + t115;
t1 = t2 * t65;
t10 = -t62 * qJD(1) + t65 * t17;
t107 = qJD(1) * qJD(5);
t48 = t62 * t107;
t94 = -qJD(5) * t17 - qJDD(1);
t3 = -t62 * t5 + t94 * t65 + t48;
t138 = (-t10 * t62 - t65 * t9) * qJD(5) - t3 * t62 + t1;
t36 = t63 * t64 - t66 * t67;
t31 = t36 * qJD(2);
t135 = pkin(3) * t112 + t31;
t137 = t114 * t67 - t61 * t64;
t89 = g(1) * t28 + g(2) * t29;
t70 = t89 + t138;
t136 = t36 * t57;
t127 = t63 * pkin(3);
t37 = t63 * t67 + t66 * t64;
t30 = t37 * qJD(2);
t92 = qJD(4) * t127 - t30;
t88 = t10 * t65 - t9 * t62;
t68 = qJD(5) ^ 2;
t134 = pkin(7) * t68 - t129;
t50 = pkin(7) + t127;
t126 = t66 * pkin(3);
t51 = -pkin(4) - t126;
t133 = t50 * t68 + t51 * t56 + t92 * t57;
t128 = t57 * pkin(4);
t123 = t136 * t57;
t20 = -t63 * t113 + t66 * t42;
t122 = t20 * t57;
t121 = t21 * t57;
t118 = t65 * t56;
t58 = t62 ^ 2;
t59 = t65 ^ 2;
t117 = t58 - t59;
t116 = t58 + t59;
t111 = qJD(5) * t65;
t110 = t10 * qJD(5);
t55 = t57 ^ 2;
t105 = t62 * t55 * t65;
t102 = t116 * t56;
t16 = -t20 - t128;
t96 = t16 * t111 + t139 * t62;
t95 = t62 * t57 * t111;
t93 = t137 * pkin(3);
t87 = -g(1) * t114 + g(2) * t61;
t12 = t57 * t37;
t86 = -t12 * t57 - t36 * t56;
t84 = -g(3) - t94;
t33 = -t114 * t64 - t61 * t67;
t81 = t90 + t121;
t80 = t37 * t68 - t86;
t79 = -t16 * t57 - t5 - t89;
t78 = t33 * pkin(3);
t77 = -pkin(7) * qJDD(5) + (t20 - t128) * qJD(5);
t76 = 0.2e1 * t136 * qJD(5) - qJDD(5) * t37;
t75 = -t89 + t91;
t72 = -qJDD(5) * t50 + (t51 * t57 - t135) * qJD(5);
t69 = qJD(3) ^ 2;
t60 = qJDD(1) - g(3);
t39 = qJDD(5) * t65 - t68 * t62;
t38 = qJDD(5) * t62 + t68 * t65;
t25 = t29 * pkin(4);
t24 = t28 * pkin(4);
t23 = t59 * t56 - 0.2e1 * t95;
t22 = t58 * t56 + 0.2e1 * t95;
t15 = -0.2e1 * t117 * t57 * qJD(5) + 0.2e1 * t62 * t118;
t13 = t16 * qJD(5) * t62;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, 0, 0, 0, -t39, t38, 0, -t88 * qJD(5) - t2 * t62 - t3 * t65 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) + t87, 0, 0, 0, 0, 0, 0, t67 * qJDD(3) - t69 * t64, -qJDD(3) * t64 - t69 * t67, 0, (t64 ^ 2 + t67 ^ 2) * qJDD(2) + t87, 0, 0, 0, 0, 0, 0, t86, -t37 * t56 + t123, 0, -t20 * t12 - t136 * t21 - t8 * t36 - t37 * t91 + t87, 0, 0, 0, 0, 0, 0, t76 * t62 - t80 * t65, t80 * t62 + t76 * t65, t37 * t102 - t116 * t123, t16 * t12 - t136 * t88 + t138 * t37 + t6 * t36 + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -g(1) * t137 - g(2) * t33 + t52, -g(1) * t33 + g(2) * t137 - t108, 0, 0, 0, 0, 0, 0, 0, t56, t30 * t57 + t27 + (pkin(3) * t56 - t100) * t66 + ((-pkin(3) * t57 - t42) * qJD(4) + t82) * t63 + t90, -t31 * t57 + (-t57 * t112 - t56 * t63) * pkin(3) + t75, 0, -g(1) * t93 - g(2) * t78 + t8 * t126 - t127 * t91 + t135 * t21 - t92 * t20, t22, t15, t38, t23, t39, 0, t13 + t72 * t62 + (-t133 - t139) * t65, t133 * t62 + t72 * t65 + t96, t135 * t57 * t116 + t50 * t102 + t70, t6 * t51 - g(1) * (-t28 * pkin(7) - t25 + t93) - g(2) * (-t29 * pkin(7) + t24 + t78) + t92 * t16 + t135 * t88 + (-t9 * t111 + t1 + (-t110 - t3) * t62) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t8 + t81, t75 + t122, 0, 0, t22, t15, t38, t23, t39, 0, t13 + t77 * t62 + (-t6 + t81 - t134) * t65, t77 * t65 + (-t121 + t134) * t62 + t96, pkin(7) * t102 - t116 * t122 + t70, -t6 * pkin(4) + t70 * pkin(7) + g(1) * t25 - g(2) * t24 - t16 * t21 - t88 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, t117 * t55, t62 * t56, t105, t118, qJDD(5), t79 * t62 - t84 * t65 + t110 + t48, t115 + t84 * t62 + (t79 + t107) * t65, 0, 0;];
tau_reg = t4;
