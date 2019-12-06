% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PPRRR3
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PPRRR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:17:00
% EndTime: 2019-12-05 15:17:04
% DurationCPUTime: 0.89s
% Computational Cost: add. (2162->176), mult. (4102->259), div. (0->0), fcn. (2947->10), ass. (0->106)
t92 = sin(qJ(4));
t110 = qJD(3) * t92;
t91 = sin(qJ(5));
t94 = cos(qJ(5));
t95 = cos(qJ(4));
t55 = -t94 * t95 * qJD(3) + t91 * t110;
t57 = (t95 * t91 + t92 * t94) * qJD(3);
t41 = t57 * t55;
t82 = qJDD(4) + qJDD(5);
t125 = -t41 + t82;
t127 = t125 * t91;
t126 = t125 * t94;
t107 = qJD(3) * qJD(4);
t106 = t92 * t107;
t78 = t95 * qJDD(3);
t102 = t78 - t106;
t105 = t95 * t107;
t109 = t92 * qJDD(3);
t62 = t105 + t109;
t30 = -t55 * qJD(5) + t91 * t102 + t94 * t62;
t83 = qJD(4) + qJD(5);
t52 = t83 * t55;
t124 = -t52 + t30;
t111 = sin(pkin(8));
t112 = cos(pkin(8));
t100 = -t111 * g(1) + t112 * g(2) + qJDD(2);
t68 = -t112 * g(1) - t111 * g(2);
t86 = -g(3) + qJDD(1);
t88 = sin(pkin(9));
t89 = cos(pkin(9));
t48 = t89 * t68 + t88 * t86;
t93 = sin(qJ(3));
t96 = cos(qJ(3));
t35 = t93 * t100 + t96 * t48;
t98 = qJD(3) ^ 2;
t33 = -t98 * pkin(3) + qJDD(3) * pkin(6) + t35;
t47 = t88 * t68 - t89 * t86;
t19 = t92 * t33 - t95 * t47;
t123 = -t19 + (-t62 + t105) * pkin(7);
t53 = t55 ^ 2;
t54 = t57 ^ 2;
t81 = t83 ^ 2;
t20 = t95 * t33 + t92 * t47;
t71 = qJD(4) * pkin(4) - pkin(7) * t110;
t85 = t95 ^ 2;
t80 = t85 * t98;
t14 = -pkin(4) * t80 + t102 * pkin(7) - qJD(4) * t71 + t20;
t74 = t92 * t98 * t95;
t108 = qJDD(4) + t74;
t99 = t108 * pkin(4) + t123;
t6 = t91 * t14 - t94 * t99;
t116 = t94 * t14;
t7 = t91 * t99 + t116;
t2 = -t94 * t6 + t91 * t7;
t122 = t92 * t2;
t121 = t83 * t91;
t120 = t83 * t94;
t34 = t96 * t100 - t93 * t48;
t32 = -qJDD(3) * pkin(3) - t98 * pkin(6) - t34;
t17 = -t102 * pkin(4) - pkin(7) * t80 + t71 * t110 + t32;
t119 = t91 * t17;
t38 = t41 + t82;
t118 = t91 * t38;
t117 = t92 * t108;
t115 = t94 * t17;
t114 = t94 * t38;
t70 = qJDD(4) - t74;
t113 = t95 * t70;
t3 = t91 * t6 + t94 * t7;
t9 = t92 * t19 + t95 * t20;
t104 = -t94 * t102 + t91 * t62;
t63 = t78 - 0.2e1 * t106;
t101 = (-qJD(5) + t83) * t57 - t104;
t97 = qJD(4) ^ 2;
t84 = t92 ^ 2;
t79 = t84 * t98;
t73 = -t80 - t97;
t72 = -t79 - t97;
t67 = t79 + t80;
t66 = -t93 * qJDD(3) - t96 * t98;
t65 = t96 * qJDD(3) - t93 * t98;
t64 = (t84 + t85) * qJDD(3);
t61 = 0.2e1 * t105 + t109;
t50 = -t54 + t81;
t49 = t53 - t81;
t45 = -t54 - t81;
t44 = -t92 * t72 - t113;
t43 = t95 * t73 - t117;
t42 = t89 * t47;
t40 = t54 - t53;
t36 = -t81 - t53;
t31 = -t53 - t54;
t29 = -t57 * qJD(5) - t104;
t28 = -t91 * t45 - t114;
t27 = t94 * t45 - t118;
t26 = t52 + t30;
t21 = (qJD(5) + t83) * t57 + t104;
t16 = t94 * t36 - t127;
t15 = t91 * t36 + t126;
t12 = -t92 * t27 + t95 * t28;
t11 = t101 * t94 + t91 * t26;
t10 = t101 * t91 - t94 * t26;
t8 = -t92 * t15 + t95 * t16;
t4 = -t92 * t10 + t95 * t11;
t1 = t95 * t3 - t122;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88 * t48 - t42, 0, 0, 0, 0, 0, 0, t88 * t66, -t88 * t65, 0, t88 * (-t93 * t34 + t96 * t35) - t42, 0, 0, 0, 0, 0, 0, t88 * (t96 * t43 - t93 * t63) + t89 * (-t108 * t95 - t92 * t73), t88 * (t96 * t44 + t93 * t61) + t89 * (t92 * t70 - t95 * t72), t88 * (t96 * t64 - t93 * t67), t88 * (t93 * t32 + t96 * t9) + t89 * (t95 * t19 - t92 * t20), 0, 0, 0, 0, 0, 0, t88 * (t93 * t21 + t96 * t8) + t89 * (-t95 * t15 - t92 * t16), t88 * (t96 * t12 + t124 * t93) + t89 * (-t95 * t27 - t92 * t28), t88 * (t93 * t31 + t96 * t4) + t89 * (-t95 * t10 - t92 * t11), t88 * (t96 * t1 + t93 * t17) + t89 * (-t95 * t2 - t92 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, 0, 0, 0, 0, 0, 0, t65, t66, 0, t96 * t34 + t93 * t35, 0, 0, 0, 0, 0, 0, t93 * t43 + t96 * t63, t93 * t44 - t96 * t61, t93 * t64 + t96 * t67, -t96 * t32 + t93 * t9, 0, 0, 0, 0, 0, 0, -t96 * t21 + t93 * t8, t93 * t12 - t124 * t96, -t96 * t31 + t93 * t4, t93 * t1 - t96 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t34, -t35, 0, 0, (t62 + t105) * t92, t95 * t61 + t92 * t63, t117 + t95 * (-t79 + t97), t63 * t95, t92 * (t80 - t97) + t113, 0, pkin(3) * t63 + pkin(6) * t43 - t95 * t32, -pkin(3) * t61 + pkin(6) * t44 + t92 * t32, pkin(3) * t67 + pkin(6) * t64 + t9, -pkin(3) * t32 + pkin(6) * t9, t92 * (-t57 * t121 + t94 * t30) + t95 * (t57 * t120 + t91 * t30), t92 * (-t124 * t91 - t94 * t21) + t95 * (t124 * t94 - t91 * t21), t92 * (-t91 * t50 + t126) + t95 * (t94 * t50 + t127), t92 * (t55 * t120 - t91 * t29) + t95 * (t55 * t121 + t94 * t29), t92 * (t94 * t49 - t118) + t95 * (t91 * t49 + t114), (t92 * (-t55 * t94 + t57 * t91) + t95 * (-t55 * t91 - t57 * t94)) * t83, t92 * (-pkin(7) * t15 + t119) + t95 * (-pkin(4) * t21 + pkin(7) * t16 - t115) - pkin(3) * t21 + pkin(6) * t8, t92 * (-pkin(7) * t27 + t115) + t95 * (-pkin(4) * t124 + pkin(7) * t28 + t119) - pkin(3) * t124 + pkin(6) * t12, t92 * (-pkin(7) * t10 - t2) + t95 * (-pkin(4) * t31 + pkin(7) * t11 + t3) - pkin(3) * t31 + pkin(6) * t4, -pkin(7) * t122 + t95 * (-pkin(4) * t17 + pkin(7) * t3) - pkin(3) * t17 + pkin(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t79 - t80, t109, t74, t78, qJDD(4), -t19, -t20, 0, 0, t41, t40, t26, -t41, t101, t82, pkin(4) * t15 - t6, -t116 - t91 * t123 + (-t108 * t91 + t27) * pkin(4), pkin(4) * t10, pkin(4) * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t40, t26, -t41, t101, t82, -t6, -t7, 0, 0;];
tauJ_reg = t5;
