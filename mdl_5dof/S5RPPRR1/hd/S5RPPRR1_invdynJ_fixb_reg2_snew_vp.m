% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRR1
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:15
% EndTime: 2019-12-05 17:38:19
% DurationCPUTime: 0.77s
% Computational Cost: add. (1658->165), mult. (3301->194), div. (0->0), fcn. (1818->6), ass. (0->111)
t81 = sin(qJ(5));
t82 = sin(qJ(4));
t84 = cos(qJ(5));
t85 = cos(qJ(4));
t52 = (-t81 * t85 - t82 * t84) * qJD(1);
t111 = qJD(1) * t85;
t54 = -qJD(1) * t81 * t82 + t111 * t84;
t128 = t54 * t52;
t75 = qJDD(4) + qJDD(5);
t132 = t75 + t128;
t135 = t132 * t81;
t134 = t132 * t84;
t76 = qJD(4) + qJD(5);
t127 = t76 * t52;
t107 = qJD(1) * qJD(4);
t102 = t85 * t107;
t109 = t82 * qJDD(1);
t59 = -t102 - t109;
t103 = t82 * t107;
t108 = t85 * qJDD(1);
t60 = -t103 + t108;
t27 = qJD(5) * t52 + t59 * t81 + t60 * t84;
t133 = t27 + t127;
t88 = qJD(1) ^ 2;
t115 = t85 * t88;
t114 = pkin(1) + qJ(3);
t77 = qJDD(1) * qJ(2);
t83 = sin(qJ(1));
t86 = cos(qJ(1));
t98 = t86 * g(1) + t83 * g(2);
t95 = 0.2e1 * qJD(2) * qJD(1) - t98;
t94 = qJDD(3) + t95;
t44 = -t114 * t88 + t77 + t94;
t39 = -qJDD(1) * pkin(6) + t44;
t116 = t85 * t39;
t129 = t60 * pkin(7);
t131 = (-pkin(4) * t115 - pkin(7) * t107 + g(3)) * t82 + qJDD(4) * pkin(4) + t116 - t129;
t50 = t52 ^ 2;
t51 = t54 ^ 2;
t74 = t76 ^ 2;
t130 = 2 * qJD(3);
t126 = t76 * t81;
t125 = t76 * t84;
t78 = t82 ^ 2;
t124 = t78 * t88;
t79 = t85 ^ 2;
t123 = t79 * t88;
t65 = qJD(4) * pkin(4) - pkin(7) * t111;
t80 = t88 * pkin(6);
t100 = t114 * qJDD(1);
t104 = t83 * g(1) - t86 * g(2);
t97 = -qJDD(2) + t104;
t93 = t88 * qJ(2) + t97;
t89 = t100 + t93;
t25 = -t59 * pkin(4) - pkin(7) * t124 - t80 + (t65 * t85 + t130) * qJD(1) + t89;
t122 = t81 * t25;
t31 = -t128 + t75;
t121 = t81 * t31;
t106 = t82 * t115;
t120 = t82 * (qJDD(4) + t106);
t34 = g(3) * t85 - t39 * t82;
t24 = -pkin(4) * t124 + pkin(7) * t59 - qJD(4) * t65 - t34;
t119 = t84 * t24;
t118 = t84 * t25;
t117 = t84 * t31;
t113 = qJ(2) - pkin(6);
t112 = t78 + t79;
t110 = qJDD(1) * pkin(1);
t105 = qJD(1) * t130;
t10 = -t131 * t84 + t81 * t24;
t11 = t131 * t81 + t119;
t4 = t10 * t81 + t11 * t84;
t101 = -t59 * t84 + t81 * t60;
t99 = pkin(4) * t82 + t114;
t3 = -t10 * t84 + t11 * t81;
t1 = t85 * t3 + t82 * t4;
t33 = t82 * g(3) + t116;
t96 = qJDD(4) - t106;
t14 = t85 * t33 - t82 * t34;
t91 = (-qJD(5) + t76) * t54 - t101;
t43 = t89 + t105;
t87 = qJD(4) ^ 2;
t70 = 0.2e1 * t77;
t63 = t112 * t88;
t62 = t112 * qJDD(1);
t61 = -0.2e1 * t103 + t108;
t58 = 0.2e1 * t102 + t109;
t57 = t85 * t96;
t49 = t93 + t110;
t46 = -t51 + t74;
t45 = t50 - t74;
t42 = -t51 - t74;
t41 = -t120 + t85 * (-t87 - t123);
t40 = t82 * (-t87 - t124) + t57;
t38 = -t80 + t43;
t35 = t51 - t50;
t29 = -t74 - t50;
t28 = -t50 - t51;
t26 = -qJD(5) * t54 - t101;
t23 = -t42 * t81 - t117;
t22 = t42 * t84 - t121;
t21 = t27 - t127;
t16 = (qJD(5) + t76) * t54 + t101;
t13 = t29 * t84 - t135;
t12 = t29 * t81 + t134;
t9 = t22 * t85 + t23 * t82;
t8 = t21 * t81 + t84 * t91;
t7 = -t21 * t84 + t81 * t91;
t5 = t12 * t85 + t13 * t82;
t2 = t7 * t85 + t8 * t82;
t6 = [0, 0, 0, 0, 0, qJDD(1), t104, t98, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t97 - 0.2e1 * t110, t70 + t95, qJ(2) * (-pkin(1) * t88 + t77 + t95) + pkin(1) * t49, qJDD(1), 0, 0, 0, 0, 0, 0, t70 + t94, 0.2e1 * t100 + t105 + t97, qJ(2) * t44 + t114 * t43, (t60 - t103) * t85, -t58 * t85 - t61 * t82, t57 - t82 * (t87 - t123), (-t59 + t102) * t82, t85 * (-t87 + t124) - t120, 0, t113 * t40 + t114 * t58 + t82 * t38, t113 * t41 + t114 * t61 + t85 * t38, -t113 * t62 - t114 * t63 - t14, t113 * t14 + t114 * t38, t85 * (-t126 * t54 + t27 * t84) - t82 * (t125 * t54 + t27 * t81), t85 * (-t133 * t81 - t16 * t84) - t82 * (t133 * t84 - t16 * t81), t85 * (-t46 * t81 + t134) - t82 * (t46 * t84 + t135), t85 * (-t125 * t52 - t26 * t81) - t82 * (-t126 * t52 + t26 * t84), t85 * (t45 * t84 - t121) - t82 * (t45 * t81 + t117), (t85 * (t52 * t84 + t54 * t81) - t82 * (t52 * t81 - t54 * t84)) * t76, t85 * (-pkin(7) * t12 + t122) - t82 * (pkin(7) * t13 - t118) + t113 * t5 + t99 * t16, t85 * (-pkin(7) * t22 + t118) - t82 * (pkin(7) * t23 + t122) + t113 * t9 + t99 * t133, t85 * (-pkin(7) * t7 - t3) - t82 * (pkin(7) * t8 + t4) + t99 * t28 + t113 * t2, t99 * t25 + (-pkin(7) + t113) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t88, -t49, 0, 0, 0, 0, 0, 0, 0, -t88, -qJDD(1), -t43, 0, 0, 0, 0, 0, 0, -t58, -t61, t63, -t38, 0, 0, 0, 0, 0, 0, -t16, -t133, -t28, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t88, t44, 0, 0, 0, 0, 0, 0, t40, t41, -t62, t14, 0, 0, 0, 0, 0, 0, t5, t9, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, (-t78 + t79) * t88, t108, -t106, -t109, qJDD(4), t33, t34, 0, 0, -t128, t35, t21, t128, t91, t75, pkin(4) * t12 - t10, -t119 - t81 * (-pkin(7) * t103 - t129 + t33) + (-t81 * t96 + t22) * pkin(4), pkin(4) * t7, pkin(4) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, t35, t21, t128, t91, t75, -t10, -t11, 0, 0;];
tauJ_reg = t6;
