% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PPRRP4
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:45
% EndTime: 2019-12-31 17:34:46
% DurationCPUTime: 0.74s
% Computational Cost: add. (621->179), mult. (1341->221), div. (0->0), fcn. (856->6), ass. (0->111)
t111 = sin(pkin(7));
t112 = cos(pkin(7));
t63 = sin(qJ(3));
t65 = cos(qJ(3));
t29 = -t111 * t63 - t112 * t65;
t30 = -t111 * t65 + t112 * t63;
t77 = g(1) * t29 + g(2) * t30;
t64 = cos(qJ(4));
t40 = qJD(3) * pkin(6) + t63 * qJD(2);
t83 = qJ(5) * qJD(3) + t40;
t131 = t83 * t64;
t53 = t64 * qJDD(3);
t62 = sin(qJ(4));
t92 = qJD(3) * qJD(4);
t88 = t62 * t92;
t130 = -0.2e1 * t88 + t53;
t54 = t62 * qJD(1);
t11 = -t54 + t131;
t103 = t64 * qJD(1);
t10 = -t83 * t62 - t103;
t113 = qJD(4) * pkin(4);
t9 = t10 + t113;
t76 = t11 * t64 - t62 * t9;
t50 = t64 * pkin(4) + pkin(3);
t100 = qJDD(3) * t50;
t43 = pkin(4) * t88;
t93 = qJD(2) * qJD(3);
t47 = t63 * t93;
t97 = t65 * qJDD(2);
t82 = t47 - t97;
t8 = qJDD(5) + t43 + t82 - t100;
t129 = t76 * qJD(3) - t8;
t58 = t62 ^ 2;
t128 = pkin(4) * t58;
t57 = g(3) * t64;
t125 = -t10 + t9;
t124 = t29 * t62;
t123 = t30 * t62;
t122 = t64 * t40;
t67 = qJD(3) ^ 2;
t121 = t64 * t67;
t120 = t67 * t65;
t119 = qJ(5) + pkin(6);
t102 = t65 * qJD(2);
t89 = qJD(4) * t102;
t118 = t64 * t47 + t62 * t89;
t59 = t64 ^ 2;
t117 = t58 - t59;
t116 = t58 + t59;
t66 = qJD(4) ^ 2;
t115 = t66 + t67;
t114 = qJD(3) * pkin(3);
t108 = qJD(3) * t50;
t24 = qJD(5) - t102 - t108;
t110 = qJD(3) * t24;
t41 = -t102 - t114;
t109 = qJD(3) * t41;
t107 = qJD(3) * t64;
t106 = qJD(4) * t62;
t105 = qJDD(3) * pkin(3);
t104 = qJDD(4) * pkin(4);
t101 = -qJD(5) - t24;
t99 = qJDD(4) * t62;
t51 = t62 * qJDD(3);
t98 = t63 * qJDD(2);
t96 = t65 * qJDD(3);
t95 = qJ(5) * qJDD(3);
t94 = qJD(1) * qJD(4);
t91 = qJD(3) * qJD(5);
t90 = -g(1) * t123 + g(2) * t124 + t64 * t89;
t87 = t64 * t92;
t86 = qJD(4) * t119;
t26 = qJDD(3) * pkin(6) + t65 * t93 + t98;
t3 = -t62 * qJDD(1) - t40 * t106 + (t26 - t94) * t64;
t85 = t8 - t100;
t84 = t116 * qJDD(3);
t81 = -t26 - t95;
t80 = t62 * t87;
t79 = t116 * t102;
t78 = g(1) * t30 - g(2) * t29;
t75 = -g(1) * t111 + g(2) * t112;
t20 = -t62 * t40 - t103;
t21 = -t54 + t122;
t74 = t20 * t62 - t21 * t64;
t25 = t82 - t105;
t73 = pkin(6) * t66 - t105 + t25;
t48 = t62 * t94;
t72 = -g(1) * t124 - g(2) * t123 - t64 * qJDD(1) + t48 + t57;
t71 = -g(3) * t62 - t77 * t64 - t3;
t70 = -pkin(6) * qJDD(4) + (t41 - t114) * qJD(4);
t4 = -t62 * t26 + t48 + (-qJD(4) * t40 - qJDD(1)) * t64;
t69 = t3 * t64 + (-t20 * t64 - t21 * t62) * qJD(4) - t4 * t62;
t68 = t69 + t77;
t60 = qJDD(1) - g(3);
t42 = t62 * t121;
t36 = t119 * t64;
t35 = t119 * t62;
t34 = t117 * t67;
t33 = qJDD(4) * t64 - t66 * t62;
t32 = t66 * t64 + t99;
t28 = t59 * qJDD(3) - 0.2e1 * t80;
t27 = t58 * qJDD(3) + 0.2e1 * t80;
t23 = -t62 * qJD(5) - t64 * t86;
t22 = t64 * qJD(5) - t62 * t86;
t12 = -0.2e1 * t117 * t92 + 0.2e1 * t62 * t53;
t7 = t116 * t120 + t63 * t84;
t6 = t130 * t65 + (-t115 * t64 - t99) * t63;
t5 = (-qJDD(4) * t63 - 0.2e1 * t65 * t92) * t64 + (t115 * t63 - t96) * t62;
t2 = t64 * t91 + (-t88 + t53) * qJ(5) + t3;
t1 = t104 + t48 + (-t83 * qJD(4) - qJDD(1)) * t64 + (t81 - t91) * t62;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, 0, 0, 0, -t33, t32, 0, t74 * qJD(4) - t3 * t62 - t4 * t64 - g(3), 0, 0, 0, 0, 0, 0, -t33, t32, 0, -t76 * qJD(4) - t1 * t64 - t2 * t62 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) + t75, 0, 0, 0, 0, 0, 0, -t67 * t63 + t96, -qJDD(3) * t63 - t120, 0, (t63 ^ 2 + t65 ^ 2) * qJDD(2) + t75, 0, 0, 0, 0, 0, 0, t6, t5, t7, (-t74 * qJD(3) - t25) * t65 + (t69 + t109) * t63 + t75, 0, 0, 0, 0, 0, 0, t6, t5, t7, t129 * t65 + (t110 - t1 * t62 + t2 * t64 + (-t11 * t62 - t64 * t9) * qJD(4)) * t63 + t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t78 + t97, -t77 - t98, 0, 0, t27, t12, t32, t28, t33, 0, t70 * t62 + (-t73 + t78) * t64 + t118, t70 * t64 + (t73 - t47) * t62 + t90, pkin(6) * t84 - qJD(3) * t79 + t68, (-t41 * t63 + t74 * t65) * qJD(2) + (-t25 + t78) * pkin(3) + t68 * pkin(6), t27, t12, t32, t28, t33, 0, -t35 * qJDD(4) + (t23 + (t24 - t108) * t62) * qJD(4) + (-t43 + t78 - t85) * t64 + t118, -t36 * qJDD(4) + (t85 - t47) * t62 + (t24 * t64 - t22 + (-t50 * t64 + t128) * qJD(3)) * qJD(4) + t90, (-qJD(4) * t9 + qJDD(3) * t36 + t2) * t64 + (-t11 * qJD(4) + qJDD(3) * t35 - t1) * t62 + (t22 * t64 - t23 * t62 + (t35 * t64 - t36 * t62) * qJD(4) - t79) * qJD(3) + t77, t2 * t36 + t11 * t22 - t1 * t35 + t9 * t23 - t8 * t50 + t24 * pkin(4) * t106 - g(1) * (-t119 * t29 - t30 * t50) - g(2) * (-t119 * t30 + t29 * t50) + (-t24 * t63 - t76 * t65) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t34, t51, t42, t53, qJDD(4), (-t26 - t109) * t62 + (t21 - t122) * qJD(4) + t72, t20 * qJD(4) - t41 * t107 + t71, 0, 0, -t42, t34, t51, t42, t53, qJDD(4), 0.2e1 * t104 + (t11 - t131) * qJD(4) + (pkin(4) * t121 + t101 * qJD(3) + t81) * t62 + t72, -t67 * t128 - t64 * t95 + t10 * qJD(4) + (qJ(5) * t106 + t101 * t64) * qJD(3) + t71, -pkin(4) * t51 + (-t113 + t125) * t107, t125 * t11 + (t57 + t1 + (-t77 - t110) * t62) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t130, t51 + 0.2e1 * t87, -t116 * t67, -t129 - t78;];
tau_reg = t13;
