% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPRR4
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPRR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:35
% EndTime: 2019-12-31 16:50:37
% DurationCPUTime: 0.68s
% Computational Cost: add. (1521->162), mult. (2969->219), div. (0->0), fcn. (1832->8), ass. (0->106)
t96 = qJD(1) ^ 2;
t89 = sin(qJ(4));
t90 = sin(qJ(3));
t118 = qJD(1) * t90;
t92 = cos(qJ(4));
t52 = -t92 * qJD(3) + t89 * t118;
t54 = t89 * qJD(3) + t92 * t118;
t128 = t54 * t52;
t113 = qJD(1) * qJD(3);
t74 = t90 * t113;
t93 = cos(qJ(3));
t76 = t93 * qJDD(1);
t60 = t76 - t74;
t51 = -qJDD(4) + t60;
t98 = -t51 - t128;
t130 = t89 * t98;
t129 = t92 * t98;
t109 = t93 * t113;
t114 = t90 * qJDD(1);
t59 = t109 + t114;
t107 = -t92 * qJDD(3) + t89 * t59;
t70 = t93 * qJD(1) - qJD(4);
t21 = (qJD(4) + t70) * t54 + t107;
t49 = t52 ^ 2;
t50 = t54 ^ 2;
t68 = t70 ^ 2;
t127 = t70 * t89;
t126 = t70 * t92;
t104 = -t93 * pkin(3) - t90 * pkin(6);
t91 = sin(qJ(1));
t94 = cos(qJ(1));
t110 = t91 * g(1) - t94 * g(2);
t55 = qJDD(1) * pkin(1) + t110;
t103 = t94 * g(1) + t91 * g(2);
t56 = -t96 * pkin(1) - t103;
t85 = sin(pkin(7));
t86 = cos(pkin(7));
t119 = t85 * t55 + t86 * t56;
t35 = -t96 * pkin(2) + qJDD(1) * pkin(5) + t119;
t105 = t96 * t104 + t35;
t117 = -g(3) + qJDD(2);
t73 = t93 * t117;
t95 = qJD(3) ^ 2;
t19 = -qJDD(3) * pkin(3) - t95 * pkin(6) + t105 * t90 - t73;
t125 = t89 * t19;
t31 = t51 - t128;
t124 = t89 * t31;
t69 = t93 * t96 * t90;
t64 = qJDD(3) + t69;
t123 = t90 * t64;
t122 = t92 * t19;
t121 = t92 * t31;
t65 = qJDD(3) - t69;
t120 = t93 * t65;
t116 = qJD(4) - t70;
t112 = t93 * t128;
t111 = pkin(1) * t85 + pkin(5);
t100 = -t60 + t74;
t101 = t59 + t109;
t108 = t86 * t55 - t85 * t56;
t34 = -qJDD(1) * pkin(2) - t96 * pkin(5) - t108;
t16 = t100 * pkin(3) - t101 * pkin(6) + t34;
t106 = t90 * t117;
t20 = -t95 * pkin(3) + qJDD(3) * pkin(6) + t105 * t93 + t106;
t6 = -t92 * t16 + t89 * t20;
t7 = t89 * t16 + t92 * t20;
t3 = t89 * t6 + t92 * t7;
t28 = t90 * t35 - t73;
t29 = t93 * t35 + t106;
t12 = t90 * t28 + t93 * t29;
t102 = t92 * t6 - t89 * t7;
t99 = -t89 * qJDD(3) - t92 * t59;
t97 = -pkin(1) * t86 - pkin(2) + t104;
t37 = -t52 * qJD(4) - t99;
t82 = t93 ^ 2;
t81 = t90 ^ 2;
t79 = t82 * t96;
t77 = t81 * t96;
t67 = -t79 - t95;
t66 = -t77 - t95;
t63 = t77 + t79;
t62 = (t81 + t82) * qJDD(1);
t61 = t76 - 0.2e1 * t74;
t58 = 0.2e1 * t109 + t114;
t45 = t52 * t70;
t44 = -t50 + t68;
t43 = t49 - t68;
t42 = -t90 * t66 - t120;
t41 = t93 * t67 - t123;
t40 = t50 - t49;
t39 = -t50 - t68;
t38 = -t68 - t49;
t36 = -t54 * qJD(4) - t107;
t30 = t49 + t50;
t26 = t116 * t52 + t99;
t25 = t37 - t45;
t24 = t37 + t45;
t22 = -t116 * t54 - t107;
t18 = -t89 * t39 + t121;
t17 = t92 * t39 + t124;
t15 = t92 * t38 - t130;
t14 = t89 * t38 + t129;
t11 = -t21 * t92 + t89 * t25;
t9 = t93 * t18 - t90 * t26;
t8 = t93 * t15 - t90 * t22;
t1 = [0, 0, 0, 0, 0, qJDD(1), t110, t103, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t86 * qJDD(1) - t85 * t96) + t108, pkin(1) * (-t85 * qJDD(1) - t86 * t96) - t119, 0, pkin(1) * (t108 * t86 + t119 * t85), t101 * t90, t93 * t58 + t90 * t61, t123 + t93 * (-t77 + t95), -t100 * t93, t90 * (t79 - t95) + t120, 0, -t93 * t34 + pkin(2) * t61 + pkin(5) * t41 + pkin(1) * (t85 * t41 + t86 * t61), t90 * t34 - pkin(2) * t58 + pkin(5) * t42 + pkin(1) * (t85 * t42 - t86 * t58), pkin(2) * t63 + pkin(5) * t62 + pkin(1) * (t85 * t62 + t86 * t63) + t12, -pkin(2) * t34 + pkin(5) * t12 + pkin(1) * (t85 * t12 - t86 * t34), t90 * (t127 * t54 + t92 * t37) - t112, t90 * (t92 * t22 - t89 * t24) - t93 * t40, t90 * (-t89 * t44 + t129) - t93 * t25, t90 * (-t126 * t52 - t89 * t36) + t112, t90 * (t92 * t43 + t124) + t93 * t21, t93 * t51 + t90 * (t52 * t92 - t54 * t89) * t70, t90 * (-pkin(6) * t14 + t125) + t93 * (-pkin(3) * t14 + t6) - pkin(2) * t14 + pkin(5) * t8 + pkin(1) * (-t86 * t14 + t85 * t8), t90 * (-pkin(6) * t17 + t122) + t93 * (-pkin(3) * t17 + t7) - pkin(2) * t17 + pkin(5) * t9 + pkin(1) * (-t86 * t17 + t85 * t9), t90 * t102 + t111 * (t93 * t11 - t90 * t30) + t97 * (-t21 * t89 - t92 * t25), t111 * (t90 * t19 + t93 * t3) - t97 * t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, 0, 0, 0, 0, 0, 0, t93 * t64 + t90 * t67, -t90 * t65 + t93 * t66, 0, -t93 * t28 + t90 * t29, 0, 0, 0, 0, 0, 0, t90 * t15 + t93 * t22, t90 * t18 + t93 * t26, t90 * t11 + t93 * t30, -t93 * t19 + t90 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t77 - t79, t114, t69, t76, qJDD(3), -t28, -t29, 0, 0, -t126 * t54 + t89 * t37, t89 * t22 + t92 * t24, t92 * t44 + t130, -t127 * t52 + t92 * t36, t89 * t43 - t121, (t52 * t89 + t54 * t92) * t70, pkin(3) * t22 + pkin(6) * t15 - t122, pkin(3) * t26 + pkin(6) * t18 + t125, pkin(3) * t30 + pkin(6) * t11 + t3, -pkin(3) * t19 + pkin(6) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t40, t25, -t128, -t21, -t51, -t6, -t7, 0, 0;];
tauJ_reg = t1;
