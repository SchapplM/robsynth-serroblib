% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPRP6
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRP6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:26
% EndTime: 2019-12-05 15:41:29
% DurationCPUTime: 0.89s
% Computational Cost: add. (678->175), mult. (1259->207), div. (0->0), fcn. (721->6), ass. (0->122)
t59 = sin(qJ(2));
t115 = t59 * qJD(1);
t58 = sin(qJ(4));
t60 = cos(qJ(4));
t84 = t58 * pkin(4) - t60 * qJ(5);
t28 = qJ(3) + t84;
t120 = qJD(2) * t28;
t11 = t115 + t120;
t117 = t11 * qJD(2);
t111 = qJD(4) * qJ(5);
t139 = pkin(2) + pkin(6);
t61 = cos(qJ(2));
t114 = t61 * qJD(1);
t92 = qJD(3) - t114;
t26 = -t139 * qJD(2) + t92;
t133 = t58 * t26;
t12 = t111 + t133;
t129 = t60 * t26;
t104 = t61 * qJDD(1);
t102 = qJD(1) * qJD(2);
t47 = t59 * t102;
t82 = qJDD(3) + t47 - t104;
t13 = -t139 * qJDD(2) + t82;
t7 = t58 * t13;
t98 = qJDD(4) * qJ(5);
t2 = t98 + t7 + (qJD(5) + t129) * qJD(4);
t118 = qJDD(4) * pkin(4);
t144 = qJD(4) * t133 - t118;
t8 = t60 * t13;
t3 = qJDD(5) - t8 + t144;
t93 = qJD(4) * pkin(4) - qJD(5);
t9 = -t93 - t129;
t68 = t2 * t58 - t3 * t60 + (t12 * t60 + t9 * t58) * qJD(4);
t148 = -t68 + t117;
t54 = t58 ^ 2;
t55 = t60 ^ 2;
t125 = t54 + t55;
t105 = t59 * qJDD(1);
t110 = qJDD(2) * t28;
t85 = pkin(4) * t60 + qJ(5) * t58;
t15 = t85 * qJD(4) - t60 * qJD(5) + qJD(3);
t1 = t105 + t110 + (t15 + t114) * qJD(2);
t63 = qJD(4) ^ 2;
t56 = sin(pkin(7));
t57 = cos(pkin(7));
t90 = g(1) * t57 + g(2) * t56;
t66 = -(t90 + t102) * t61 - g(3) * t59 + t139 * t63;
t147 = qJD(2) * t15 + t1 + t110 + t66;
t146 = t125 * t13;
t145 = t125 * t26;
t112 = qJD(2) * qJ(3);
t37 = t112 + t115;
t143 = (t112 + t37 - t115) * qJD(4);
t64 = qJD(2) ^ 2;
t142 = t61 * qJDD(2) - t64 * t59;
t113 = qJDD(1) - g(3);
t140 = -t113 * t59 + t90 * t61;
t53 = g(3) * t61;
t136 = t37 * t61;
t135 = t56 * t59;
t134 = t57 * t59;
t132 = t58 * t59;
t131 = t58 * t60;
t130 = t59 * t60;
t127 = t61 * pkin(2) + t59 * qJ(3);
t126 = t54 - t55;
t124 = t63 + t64;
t123 = qJ(3) * t61;
t99 = qJDD(2) * qJ(3);
t14 = t99 + t105 + (qJD(3) + t114) * qJD(2);
t122 = t14 * qJ(3);
t119 = qJDD(2) * pkin(2);
t116 = t37 * qJD(2);
t109 = qJDD(2) * t59;
t108 = qJDD(4) * t58;
t107 = qJDD(4) * t139;
t106 = t58 * qJDD(2);
t48 = t60 * qJDD(2);
t101 = qJD(2) * qJD(3);
t100 = qJD(2) * qJD(4);
t97 = -g(1) * t134 - g(2) * t135 + t53;
t96 = t14 * t59 - g(3);
t29 = t125 * qJDD(2);
t91 = t100 * t131;
t86 = t12 * t58 - t9 * t60;
t83 = (-qJD(2) * pkin(2) + t92) * t59 + t136;
t19 = t57 * t132 + t56 * t60;
t21 = -t56 * t132 + t57 * t60;
t81 = g(1) * t19 - g(2) * t21 - t7;
t80 = t97 - t116;
t79 = -t97 + t104;
t18 = -t57 * t130 + t56 * t58;
t20 = t56 * t130 + t57 * t58;
t78 = g(1) * t18 - g(2) * t20 + t60 * t53 + t8;
t76 = t124 * t61 + t109;
t75 = -qJDD(4) * t61 + 0.2e1 * t59 * t100;
t73 = -qJDD(5) + t78;
t17 = t82 - t119;
t72 = (t11 - t115 + t120) * qJD(4);
t71 = t90 * t139;
t40 = t56 * t123;
t41 = t57 * t123;
t70 = -g(1) * t41 - g(2) * t40 - g(3) * (t61 * pkin(6) + t127);
t69 = t125 * t47 + t139 * t29 - t97;
t65 = t14 + t66 + t99 + t101;
t49 = qJDD(4) * t60;
t43 = t64 * t131;
t42 = t60 * t107;
t34 = t126 * t64;
t32 = t64 * t61 + t109;
t31 = -t63 * t58 + t49;
t30 = t63 * t60 + t108;
t27 = t85 * qJD(2);
t25 = t55 * qJDD(2) - 0.2e1 * t91;
t24 = t54 * qJDD(2) + 0.2e1 * t91;
t23 = -t124 * t58 + t49;
t22 = t124 * t60 + t108;
t10 = t126 * t100 - t58 * t48;
t6 = t142 * t125;
t5 = -t75 * t58 + t76 * t60;
t4 = t76 * t58 + t75 * t60;
t16 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t113, 0, 0, 0, 0, 0, 0, t142, -t32, 0, -g(3) + (t59 ^ 2 + t61 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, 0, -t142, t32, t83 * qJD(2) - t17 * t61 + t96, 0, 0, 0, 0, 0, 0, t4, t5, t6, -t61 * t146 + (t59 * t145 + t136) * qJD(2) + t96, 0, 0, 0, 0, 0, 0, t4, t6, -t5, -g(3) + (t86 * qJD(2) + t1) * t59 + t148 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t79, t140, 0, 0, qJDD(2), 0, 0, 0, 0, 0, 0, qJDD(3) - t79 - 0.2e1 * t119, 0.2e1 * t101 - t140 + 0.2e1 * t99, t122 + t37 * qJD(3) - t17 * pkin(2) - g(1) * (-pkin(2) * t134 + t41) - g(2) * (-pkin(2) * t135 + t40) - g(3) * t127 - t83 * qJD(1), t25, 0.2e1 * t10, t31, t24, -t30, 0, t60 * t143 + t65 * t58 - t42, (t107 - t143) * t58 + t65 * t60, -t146 + t69, t122 + t92 * t37 - t139 * t146 + (-qJD(1) * t145 + t71) * t59 + t70, t25, t31, -0.2e1 * t10, 0, t30, t24, t147 * t58 + t60 * t72 - t42, -t68 + t69, (t72 - t107) * t58 - t147 * t60, t1 * t28 + t11 * t15 + (-t11 * qJD(1) - t90 * t84) * t61 - t68 * t139 + (-g(3) * t84 - t86 * qJD(1) + t71) * t59 + t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t64, t17 + t80, 0, 0, 0, 0, 0, 0, t23, -t22, -t29, t146 + t80, 0, 0, 0, 0, 0, 0, t23, -t29, t22, t97 - t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t34, t48, -t43, -t106, qJDD(4), -t60 * t116 + t78, (t116 - t53) * t58 + t81, 0, 0, t43, t48, t34, qJDD(4), t106, -t43, 0.2e1 * t118 + (-t11 * t60 - t27 * t58) * qJD(2) + t73, -t85 * qJDD(2) + ((t12 - t111) * t60 + (t9 + t93) * t58) * qJD(2), t58 * t53 + 0.2e1 * t98 + 0.2e1 * qJD(4) * qJD(5) + (-t11 * t58 + t27 * t60) * qJD(2) - t81, t2 * qJ(5) - t3 * pkin(4) - t11 * t27 - t9 * t133 - g(1) * (-t18 * pkin(4) + t19 * qJ(5)) - g(2) * (t20 * pkin(4) - t21 * qJ(5)) + t85 * t53 + (qJD(5) - t129) * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t43, t48, -t55 * t64 - t63, -t12 * qJD(4) + t60 * t117 + t144 - t73;];
tau_reg = t16;
