% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:53
% EndTime: 2019-12-31 17:59:56
% DurationCPUTime: 0.73s
% Computational Cost: add. (2063->158), mult. (3851->210), div. (0->0), fcn. (2256->8), ass. (0->114)
t92 = cos(qJ(4));
t121 = qJD(1) * t92;
t88 = sin(qJ(5));
t91 = cos(qJ(5));
t54 = -t91 * qJD(4) + t88 * t121;
t56 = t88 * qJD(4) + t91 * t121;
t42 = t56 * t54;
t117 = qJD(1) * qJD(4);
t112 = t92 * t117;
t89 = sin(qJ(4));
t118 = t89 * qJDD(1);
t60 = -t112 - t118;
t53 = qJDD(5) - t60;
t135 = -t42 + t53;
t137 = t135 * t88;
t136 = t135 * t91;
t113 = t89 * t117;
t74 = t92 * qJDD(1);
t61 = t74 - t113;
t109 = -t91 * qJDD(4) + t88 * t61;
t72 = t89 * qJD(1) + qJD(5);
t21 = (qJD(5) - t72) * t56 + t109;
t51 = t54 ^ 2;
t52 = t56 ^ 2;
t71 = t72 ^ 2;
t134 = -pkin(2) - pkin(6);
t133 = t72 * t88;
t132 = t72 * t91;
t80 = t89 ^ 2;
t95 = qJD(1) ^ 2;
t131 = t80 * t95;
t81 = t92 ^ 2;
t130 = t81 * t95;
t35 = t42 + t53;
t129 = t88 * t35;
t115 = t89 * t95 * t92;
t67 = qJDD(4) + t115;
t128 = t89 * t67;
t127 = t91 * t35;
t82 = -g(3) + qJDD(2);
t90 = sin(qJ(1));
t93 = cos(qJ(1));
t104 = t93 * g(1) + t90 * g(2);
t57 = -t95 * pkin(1) - t104;
t85 = sin(pkin(8));
t86 = cos(pkin(8));
t114 = t90 * g(1) - t93 * g(2);
t98 = qJDD(1) * pkin(1) + t114;
t110 = -t85 * t57 + t86 * t98;
t84 = qJDD(1) * pkin(2);
t32 = -t95 * qJ(3) + qJDD(3) - t110 - t84;
t96 = -qJDD(1) * pkin(6) + t32;
t27 = t89 * t82 - t92 * t96;
t105 = t89 * pkin(4) - t92 * pkin(7);
t58 = t105 * qJD(1);
t94 = qJD(4) ^ 2;
t17 = -qJDD(4) * pkin(4) - t94 * pkin(7) + t58 * t121 + t27;
t126 = t92 * t17;
t68 = qJDD(4) - t115;
t125 = t92 * t68;
t124 = t92 * t82;
t123 = t86 * t57 + t85 * t98;
t122 = t80 + t81;
t119 = qJD(5) + t72;
t116 = t89 * t42;
t101 = -t61 + t113;
t102 = -t60 + t112;
t76 = 0.2e1 * qJD(3) * qJD(1);
t78 = qJDD(1) * qJ(3);
t108 = t76 + t78 + t123;
t30 = t134 * t95 + t108;
t14 = t102 * pkin(4) + t101 * pkin(7) + t30;
t18 = -t94 * pkin(4) + qJDD(4) * pkin(7) + t124 + (-qJD(1) * t58 + t96) * t89;
t5 = -t91 * t14 + t88 * t18;
t6 = t88 * t14 + t91 * t18;
t3 = t88 * t5 + t91 * t6;
t111 = pkin(1) * t85 + qJ(3);
t107 = pkin(1) * t86 - t134;
t106 = -pkin(1) * (t85 * qJDD(1) + t86 * t95) - t123;
t103 = t91 * t5 - t88 * t6;
t28 = t89 * t96 + t124;
t12 = -t92 * t27 + t89 * t28;
t100 = -t88 * qJDD(4) - t91 * t61;
t99 = -pkin(1) * (-t86 * qJDD(1) + t85 * t95) + t110;
t97 = t105 + t111;
t38 = -t54 * qJD(5) - t100;
t70 = -t94 - t130;
t69 = -t94 - t131;
t65 = t122 * qJDD(1);
t62 = t74 - 0.2e1 * t113;
t59 = 0.2e1 * t112 + t118;
t47 = t72 * t54;
t46 = -t52 + t71;
t45 = t51 - t71;
t44 = t92 * t70 - t128;
t43 = t89 * t69 + t125;
t41 = t52 - t51;
t40 = -t52 - t71;
t39 = -t71 - t51;
t37 = -t56 * qJD(5) - t109;
t33 = t51 + t52;
t31 = -t95 * pkin(2) + t108;
t26 = t119 * t54 + t100;
t25 = t38 + t47;
t24 = t38 - t47;
t22 = -t119 * t56 - t109;
t20 = -t88 * t40 - t127;
t16 = t91 * t39 - t137;
t11 = -t21 * t91 + t88 * t25;
t9 = t89 * t20 + t92 * t26;
t8 = t89 * t16 + t92 * t22;
t7 = t89 * t11 + t92 * t33;
t1 = t89 * t3 - t126;
t2 = [0, 0, 0, 0, 0, qJDD(1), t114, t104, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t99, t106, 0, pkin(1) * (t86 * t110 + t85 * t123), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) - 0.2e1 * t84 - t99, -t106 + t76 + 0.2e1 * t78, pkin(1) * (t85 * t31 - t86 * t32) - pkin(2) * t32 + qJ(3) * t31, -t101 * t92, -t92 * t59 - t89 * t62, t125 - t89 * (t94 - t130), t102 * t89, t92 * (-t94 + t131) - t128, 0, -t107 * t43 + t111 * t59 + t89 * t30, -t107 * t44 + t111 * t62 + t92 * t30, -t111 * t122 * t95 + t107 * t65 - t12, -t107 * t12 + t111 * t30, t92 * (-t56 * t133 + t91 * t38) + t116, t92 * (t91 * t22 - t88 * t24) + t89 * t41, t92 * (-t88 * t46 + t136) + t89 * t25, t92 * (t54 * t132 - t88 * t37) - t116, t92 * (t91 * t45 - t129) - t89 * t21, t89 * t53 + t92 * (-t54 * t91 + t56 * t88) * t72, t88 * t126 - t89 * t5 + t97 * (t88 * t39 + t136) - t107 * t8, t91 * t126 - t89 * t6 + t97 * (t91 * t40 - t129) - t107 * t9, t92 * t103 - t107 * t7 + t97 * (-t21 * t88 - t91 * t25), -t107 * t1 - t97 * t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, 0, 0, 0, 0, 0, -t89 * t68 + t92 * t69, -t92 * t67 - t89 * t70, 0, t89 * t27 + t92 * t28, 0, 0, 0, 0, 0, 0, t92 * t16 - t89 * t22, t92 * t20 - t89 * t26, t92 * t11 - t89 * t33, t89 * t17 + t92 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t95, t32, 0, 0, 0, 0, 0, 0, t43, t44, -t65, t12, 0, 0, 0, 0, 0, 0, t8, t9, t7, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, (-t80 + t81) * t95, t74, -t115, -t118, qJDD(4), -t27, -t28, 0, 0, t56 * t132 + t88 * t38, t88 * t22 + t91 * t24, t91 * t46 + t137, t54 * t133 + t91 * t37, t88 * t45 + t127, (-t54 * t88 - t56 * t91) * t72, pkin(4) * t22 + pkin(7) * t16 - t91 * t17, pkin(4) * t26 + pkin(7) * t20 + t88 * t17, pkin(4) * t33 + pkin(7) * t11 + t3, -pkin(4) * t17 + pkin(7) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t41, t25, -t42, -t21, t53, -t5, -t6, 0, 0;];
tauJ_reg = t2;
