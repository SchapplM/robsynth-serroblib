% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRPRR3
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRPRR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:47:33
% EndTime: 2019-12-05 15:47:40
% DurationCPUTime: 1.01s
% Computational Cost: add. (2728->192), mult. (5088->280), div. (0->0), fcn. (3443->10), ass. (0->116)
t101 = sin(qJ(5));
t104 = cos(qJ(5));
t105 = cos(qJ(4));
t102 = sin(qJ(4));
t122 = qJD(2) * t102;
t66 = -t104 * t105 * qJD(2) + t101 * t122;
t68 = (t105 * t101 + t102 * t104) * qJD(2);
t49 = t68 * t66;
t93 = qJDD(4) + qJDD(5);
t138 = -t49 + t93;
t140 = t101 * t138;
t139 = t104 * t138;
t120 = qJD(2) * qJD(4);
t117 = t102 * t120;
t89 = t105 * qJDD(2);
t114 = t89 - t117;
t116 = t105 * t120;
t119 = t102 * qJDD(2);
t73 = t116 + t119;
t40 = -t66 * qJD(5) + t101 * t114 + t104 * t73;
t94 = qJD(4) + qJD(5);
t63 = t94 * t66;
t137 = -t63 + t40;
t108 = qJD(2) ^ 2;
t103 = sin(qJ(2));
t106 = cos(qJ(2));
t124 = sin(pkin(8));
t125 = cos(pkin(8));
t111 = -t125 * g(1) - t124 * g(2);
t123 = -g(3) + qJDD(1);
t60 = -t103 * t111 + t106 * t123;
t110 = qJDD(2) * pkin(2) + t60;
t61 = t103 * t123 + t106 * t111;
t57 = -t108 * pkin(2) + t61;
t98 = sin(pkin(9));
t99 = cos(pkin(9));
t38 = t98 * t110 + t99 * t57;
t36 = -t108 * pkin(3) + qJDD(2) * pkin(6) + t38;
t75 = -t124 * g(1) + t125 * g(2) + qJDD(3);
t23 = t102 * t36 - t105 * t75;
t136 = -t23 + (-t73 + t116) * pkin(7);
t64 = t66 ^ 2;
t65 = t68 ^ 2;
t92 = t94 ^ 2;
t85 = t102 * t108 * t105;
t121 = qJDD(4) + t85;
t109 = t121 * pkin(4) + t136;
t24 = t102 * t75 + t105 * t36;
t82 = qJD(4) * pkin(4) - pkin(7) * t122;
t96 = t105 ^ 2;
t91 = t96 * t108;
t20 = -pkin(4) * t91 + t114 * pkin(7) - qJD(4) * t82 + t24;
t10 = t101 * t20 - t104 * t109;
t130 = t104 * t20;
t11 = t101 * t109 + t130;
t3 = -t104 * t10 + t101 * t11;
t135 = t102 * t3;
t118 = -t99 * t110 + t98 * t57;
t35 = -qJDD(2) * pkin(3) - t108 * pkin(6) + t118;
t21 = -t114 * pkin(4) - pkin(7) * t91 + t82 * t122 + t35;
t134 = t101 * t21;
t46 = t49 + t93;
t133 = t101 * t46;
t132 = t101 * t94;
t131 = t102 * t121;
t129 = t104 * t21;
t128 = t104 * t46;
t127 = t104 * t94;
t81 = qJDD(4) - t85;
t126 = t105 * t81;
t4 = t101 * t10 + t104 * t11;
t115 = t101 * t73 - t104 * t114;
t13 = t102 * t23 + t105 * t24;
t74 = t89 - 0.2e1 * t117;
t112 = (-qJD(5) + t94) * t68 - t115;
t107 = qJD(4) ^ 2;
t95 = t102 ^ 2;
t90 = t95 * t108;
t84 = -t91 - t107;
t83 = -t90 - t107;
t79 = t90 + t91;
t78 = (t95 + t96) * qJDD(2);
t77 = -t98 * qJDD(2) - t99 * t108;
t76 = t99 * qJDD(2) - t98 * t108;
t72 = 0.2e1 * t116 + t119;
t59 = -t65 + t92;
t58 = t64 - t92;
t56 = -t65 - t92;
t55 = -t102 * t83 - t126;
t54 = t105 * t84 - t131;
t50 = t98 * t78 + t99 * t79;
t48 = t65 - t64;
t44 = -t92 - t64;
t43 = t98 * t55 - t99 * t72;
t42 = t98 * t54 + t99 * t74;
t41 = -t64 - t65;
t39 = -t68 * qJD(5) - t115;
t34 = -t101 * t56 - t128;
t33 = t104 * t56 - t133;
t32 = t63 + t40;
t27 = (qJD(5) + t94) * t68 + t115;
t26 = t104 * t44 - t140;
t25 = t101 * t44 + t139;
t18 = -t118 * t99 + t98 * t38;
t17 = -t102 * t33 + t105 * t34;
t16 = t101 * t32 + t104 * t112;
t15 = t101 * t112 - t104 * t32;
t14 = -t102 * t25 + t105 * t26;
t12 = -t137 * t99 + t98 * t17;
t8 = t98 * t13 - t99 * t35;
t7 = t98 * t14 - t99 * t27;
t6 = -t102 * t15 + t105 * t16;
t5 = -t99 * t41 + t98 * t6;
t2 = t105 * t4 - t135;
t1 = t98 * t2 - t99 * t21;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t123, 0, 0, 0, 0, 0, 0, t106 * qJDD(2) - t103 * t108, -t103 * qJDD(2) - t106 * t108, 0, t103 * t61 + t106 * t60, 0, 0, 0, 0, 0, 0, t103 * t77 + t106 * t76, -t103 * t76 + t106 * t77, 0, t103 * (t118 * t98 + t99 * t38) + t106 * t18, 0, 0, 0, 0, 0, 0, t103 * (t99 * t54 - t98 * t74) + t106 * t42, t103 * (t99 * t55 + t98 * t72) + t106 * t43, t103 * (t99 * t78 - t98 * t79) + t106 * t50, t103 * (t99 * t13 + t98 * t35) + t106 * t8, 0, 0, 0, 0, 0, 0, t103 * (t99 * t14 + t98 * t27) + t106 * t7, t103 * (t137 * t98 + t99 * t17) + t106 * t12, t103 * (t98 * t41 + t99 * t6) + t106 * t5, t103 * (t99 * t2 + t98 * t21) + t106 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t60, -t61, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(2) * t76 - t118, pkin(2) * t77 - t38, 0, pkin(2) * t18, (t73 + t116) * t102, t102 * t74 + t105 * t72, t131 + t105 * (-t90 + t107), t74 * t105, t102 * (t91 - t107) + t126, 0, pkin(2) * t42 + pkin(3) * t74 + pkin(6) * t54 - t105 * t35, pkin(2) * t43 - pkin(3) * t72 + pkin(6) * t55 + t102 * t35, pkin(2) * t50 + pkin(3) * t79 + pkin(6) * t78 + t13, pkin(2) * t8 - pkin(3) * t35 + pkin(6) * t13, t102 * (t104 * t40 - t68 * t132) + t105 * (t101 * t40 + t68 * t127), t102 * (-t101 * t137 - t104 * t27) + t105 * (-t101 * t27 + t104 * t137), t102 * (-t101 * t59 + t139) + t105 * (t104 * t59 + t140), t102 * (-t101 * t39 + t66 * t127) + t105 * (t104 * t39 + t66 * t132), t102 * (t104 * t58 - t133) + t105 * (t101 * t58 + t128), (t102 * (t101 * t68 - t104 * t66) + t105 * (-t101 * t66 - t104 * t68)) * t94, t102 * (-pkin(7) * t25 + t134) + t105 * (-pkin(4) * t27 + pkin(7) * t26 - t129) - pkin(3) * t27 + pkin(6) * t14 + pkin(2) * t7, t102 * (-pkin(7) * t33 + t129) + t105 * (-pkin(4) * t137 + pkin(7) * t34 + t134) - pkin(3) * t137 + pkin(6) * t17 + pkin(2) * t12, t102 * (-pkin(7) * t15 - t3) + t105 * (-pkin(4) * t41 + pkin(7) * t16 + t4) - pkin(3) * t41 + pkin(6) * t6 + pkin(2) * t5, -pkin(7) * t135 + t105 * (-pkin(4) * t21 + pkin(7) * t4) - pkin(3) * t21 + pkin(6) * t2 + pkin(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, 0, 0, 0, 0, 0, t102 * t84 + t105 * t121, -t102 * t81 + t105 * t83, 0, t102 * t24 - t105 * t23, 0, 0, 0, 0, 0, 0, t102 * t26 + t105 * t25, t102 * t34 + t105 * t33, t102 * t16 + t105 * t15, t102 * t4 + t105 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, t90 - t91, t119, t85, t89, qJDD(4), -t23, -t24, 0, 0, t49, t48, t32, -t49, t112, t93, pkin(4) * t25 - t10, -t130 - t101 * t136 + (-t101 * t121 + t33) * pkin(4), pkin(4) * t15, pkin(4) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t48, t32, -t49, t112, t93, -t10, -t11, 0, 0;];
tauJ_reg = t9;
