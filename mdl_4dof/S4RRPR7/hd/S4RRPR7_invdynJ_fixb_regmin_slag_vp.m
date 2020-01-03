% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRPR7
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [4x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:45
% EndTime: 2019-12-31 17:06:48
% DurationCPUTime: 0.79s
% Computational Cost: add. (1078->195), mult. (2638->284), div. (0->0), fcn. (1873->10), ass. (0->113)
t114 = qJD(1) * qJD(2);
t85 = cos(qJ(2));
t109 = t85 * t114;
t82 = sin(qJ(2));
t110 = t82 * t114;
t78 = sin(pkin(7));
t79 = cos(pkin(7));
t55 = t78 * t85 + t79 * t82;
t29 = t55 * qJDD(1) + t79 * t109 - t78 * t110;
t48 = t55 * qJD(1);
t81 = sin(qJ(4));
t84 = cos(qJ(4));
t37 = t81 * qJD(2) + t84 * t48;
t10 = t37 * qJD(4) - t84 * qJDD(2) + t81 * t29;
t120 = qJD(1) * t82;
t129 = t79 * t85;
t46 = qJD(1) * t129 - t78 * t120;
t69 = t85 * pkin(2) + pkin(1);
t59 = -t69 * qJD(1) + qJD(3);
t11 = -t46 * pkin(3) - t48 * pkin(6) + t59;
t123 = qJ(3) + pkin(5);
t111 = t123 * t82;
t105 = qJD(2) * t123;
t93 = -t82 * qJD(3) - t85 * t105;
t23 = qJDD(2) * pkin(2) + t93 * qJD(1) - qJDD(1) * t111;
t45 = t85 * qJD(3) - t82 * t105;
t60 = t123 * t85;
t30 = t45 * qJD(1) + qJDD(1) * t60;
t8 = t78 * t23 + t79 * t30;
t107 = qJDD(2) * pkin(6) + qJD(4) * t11 + t8;
t13 = t79 * t45 + t78 * t93;
t58 = qJD(1) * t60;
t130 = t78 * t58;
t121 = qJD(2) * pkin(2);
t57 = qJD(1) * t111;
t53 = -t57 + t121;
t26 = t79 * t53 - t130;
t18 = -qJD(2) * pkin(3) - t26;
t47 = t55 * qJD(2);
t115 = t85 * qJDD(1);
t116 = t82 * qJDD(1);
t99 = t79 * t115 - t78 * t116;
t24 = qJD(1) * t47 + qJDD(4) - t99;
t54 = t78 * t82 - t129;
t25 = t54 * pkin(3) - t55 * pkin(6) - t69;
t34 = -t78 * t111 + t79 * t60;
t40 = qJD(4) - t46;
t7 = t79 * t23 - t78 * t30;
t5 = -qJDD(2) * pkin(3) - t7;
t50 = t54 * qJD(2);
t140 = -(qJD(4) * t25 + t13) * t40 - t107 * t54 - t18 * t50 - t34 * t24 + t5 * t55;
t83 = sin(qJ(1));
t86 = cos(qJ(1));
t101 = g(1) * t86 + g(2) * t83;
t66 = t78 * pkin(2) + pkin(6);
t75 = qJ(2) + pkin(7);
t70 = sin(t75);
t71 = cos(t75);
t139 = t101 * t70 - (pkin(2) * t120 + t48 * pkin(3) - t46 * pkin(6) + qJD(4) * t66) * t40 - g(3) * t71 - t5;
t138 = g(3) * t85;
t117 = t84 * qJD(2);
t119 = qJD(4) * t81;
t9 = qJD(4) * t117 + t81 * qJDD(2) - t48 * t119 + t84 * t29;
t137 = t9 * t81;
t136 = t18 * t55;
t135 = t25 * t24;
t35 = t81 * t48 - t117;
t134 = t35 * t40;
t133 = t37 * t40;
t132 = t37 * t48;
t131 = t48 * t35;
t51 = t79 * t58;
t128 = t81 * t24;
t127 = t83 * t81;
t126 = t83 * t84;
t17 = t84 * t24;
t125 = t86 * t81;
t124 = t86 * t84;
t27 = t78 * t53 + t51;
t76 = t82 ^ 2;
t122 = -t85 ^ 2 + t76;
t118 = qJD(4) * t84;
t113 = t82 * t121;
t19 = qJD(2) * pkin(6) + t27;
t28 = -qJD(2) * t48 + t99;
t92 = pkin(2) * t110 - t69 * qJDD(1) + qJDD(3);
t4 = -t28 * pkin(3) - t29 * pkin(6) + t92;
t106 = qJD(4) * t19 - t4;
t104 = t40 * t84;
t100 = g(1) * t83 - g(2) * t86;
t98 = t17 + (t46 * t81 - t119) * t40;
t97 = g(3) * t70 - t107;
t96 = -t55 * t119 - t84 * t50;
t95 = -0.2e1 * pkin(1) * t114 - pkin(5) * qJDD(2);
t32 = -t79 * t57 - t130;
t91 = -t66 * t24 + (t18 + t32) * t40;
t87 = qJD(2) ^ 2;
t90 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t87 + t100;
t88 = qJD(1) ^ 2;
t89 = pkin(1) * t88 - pkin(5) * qJDD(1) + t101;
t67 = -t79 * pkin(2) - pkin(3);
t44 = t71 * t124 + t127;
t43 = -t71 * t125 + t126;
t42 = -t71 * t126 + t125;
t41 = t71 * t127 + t124;
t33 = t79 * t111 + t78 * t60;
t31 = -t78 * t57 + t51;
t15 = t47 * pkin(3) + t50 * pkin(6) + t113;
t12 = t78 * t45 - t79 * t93;
t3 = t81 * t11 + t84 * t19;
t2 = t84 * t11 - t81 * t19;
t1 = t84 * t4;
t6 = [qJDD(1), t100, t101, t76 * qJDD(1) + 0.2e1 * t82 * t109, -0.2e1 * t122 * t114 + 0.2e1 * t82 * t115, qJDD(2) * t82 + t87 * t85, qJDD(2) * t85 - t87 * t82, 0, t95 * t82 + t90 * t85, -t90 * t82 + t95 * t85, t12 * t48 + t13 * t46 + t26 * t50 - t27 * t47 + t34 * t28 + t33 * t29 - t8 * t54 - t7 * t55 - t101, t8 * t34 + t27 * t13 - t7 * t33 - t26 * t12 - t92 * t69 + t59 * t113 - g(1) * (t123 * t86 - t83 * t69) - g(2) * (t123 * t83 + t86 * t69), t9 * t84 * t55 + t96 * t37, -(-t35 * t84 - t37 * t81) * t50 + (-t10 * t84 - t137 + (t35 * t81 - t37 * t84) * qJD(4)) * t55, t55 * t17 + t37 * t47 + t96 * t40 + t9 * t54, -t55 * t128 - t10 * t54 - t35 * t47 + (-t55 * t118 + t81 * t50) * t40, t24 * t54 + t40 * t47, -g(1) * t42 - g(2) * t44 + t1 * t54 + t33 * t10 + t12 * t35 + t2 * t47 + (t15 * t40 + t135 + (-t19 * t54 - t34 * t40 + t136) * qJD(4)) * t84 + t140 * t81, -g(1) * t41 - g(2) * t43 + t12 * t37 - t3 * t47 + t33 * t9 + (-(-qJD(4) * t34 + t15) * t40 - t135 + t106 * t54 - qJD(4) * t136) * t81 + t140 * t84; 0, 0, 0, -t82 * t88 * t85, t122 * t88, t116, t115, qJDD(2), t89 * t82 - t138, g(3) * t82 + t89 * t85, (t27 - t31) * t48 + (t26 - t32) * t46 + (t28 * t78 - t29 * t79) * pkin(2), t26 * t31 - t27 * t32 + (-t138 + t7 * t79 + t78 * t8 + (-qJD(1) * t59 + t101) * t82) * pkin(2), t37 * t104 + t137, (t9 - t134) * t84 + (-t10 - t133) * t81, t40 * t104 + t128 - t132, t98 + t131, -t40 * t48, t67 * t10 + t139 * t84 - t2 * t48 - t31 * t35 + t91 * t81, -t139 * t81 + t3 * t48 - t31 * t37 + t67 * t9 + t91 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46 ^ 2 - t48 ^ 2, t26 * t48 - t27 * t46 - t100 + t92, 0, 0, 0, 0, 0, t98 - t131, -t40 ^ 2 * t84 - t128 - t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t35, -t35 ^ 2 + t37 ^ 2, t9 + t134, -t10 + t133, t24, -g(1) * t43 + g(2) * t41 - t19 * t118 - t18 * t37 + t3 * t40 + t97 * t81 + t1, g(1) * t44 - g(2) * t42 + t106 * t81 + t18 * t35 + t2 * t40 + t97 * t84;];
tau_reg = t6;
