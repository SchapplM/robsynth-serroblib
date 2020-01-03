% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRRR1
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:16
% EndTime: 2019-12-31 17:22:18
% DurationCPUTime: 0.72s
% Computational Cost: add. (1153->162), mult. (1924->211), div. (0->0), fcn. (1027->12), ass. (0->114)
t67 = sin(qJ(4));
t64 = t67 ^ 2;
t71 = cos(qJ(4));
t65 = t71 ^ 2;
t121 = t64 + t65;
t62 = qJDD(1) + qJDD(2);
t56 = qJDD(3) + t62;
t73 = cos(qJ(2));
t119 = qJD(2) * t73;
t101 = qJD(1) * t119;
t72 = cos(qJ(3));
t117 = qJD(3) * t72;
t120 = pkin(1) * qJD(1);
t69 = sin(qJ(2));
t111 = t69 * t120;
t132 = t73 * pkin(1);
t54 = qJDD(1) * t132;
t21 = t62 * pkin(2) - qJD(2) * t111 + t54;
t63 = qJD(1) + qJD(2);
t29 = t63 * pkin(2) + t73 * t120;
t128 = t69 * t72;
t44 = pkin(1) * t128;
t68 = sin(qJ(3));
t96 = qJD(3) * t111;
t93 = -t72 * pkin(1) * t101 - qJDD(1) * t44 - t29 * t117 + (-t21 + t96) * t68;
t4 = t56 * pkin(7) - t93;
t107 = t121 * t4;
t66 = qJ(1) + qJ(2);
t60 = qJ(3) + t66;
t50 = cos(t60);
t137 = g(2) * t50;
t135 = t56 * pkin(3);
t109 = t69 * t117;
t115 = qJDD(1) * t69;
t118 = qJD(3) * t68;
t18 = t72 * t21;
t7 = -(qJD(1) * (t68 * t119 + t109) + t68 * t115) * pkin(1) - t29 * t118 + t18;
t5 = -t135 - t7;
t144 = t5 + t137;
t49 = sin(t60);
t124 = -g(1) * t50 - g(2) * t49;
t129 = t68 * t69;
t88 = t72 * t73 - t129;
t25 = t88 * t120;
t143 = pkin(2) * t117 - t25;
t42 = g(1) * t49;
t142 = t42 - t137;
t58 = sin(t66);
t59 = cos(t66);
t141 = g(1) * t58 - g(2) * t59;
t51 = t68 * pkin(2) + pkin(7);
t52 = -t72 * pkin(2) - pkin(3);
t57 = qJD(3) + t63;
t75 = qJD(4) ^ 2;
t89 = t68 * t73 + t128;
t24 = t89 * t120;
t94 = pkin(2) * t118 - t24;
t140 = -t51 * t75 - t52 * t56 - t57 * t94;
t138 = pkin(2) * t58;
t134 = t57 * pkin(3);
t70 = sin(qJ(1));
t133 = t70 * pkin(1);
t53 = pkin(2) + t132;
t9 = t53 * t118 + (qJD(2) * t89 + t109) * pkin(1);
t131 = t9 * t57;
t16 = t111 * t72 + t68 * t29;
t130 = t16 * t57;
t127 = t71 * t56;
t15 = -t111 * t68 + t72 * t29;
t13 = -t15 - t134;
t126 = t13 * qJD(4) * t67 + t71 * t42;
t125 = t50 * pkin(3) + t49 * pkin(7);
t27 = t68 * t53 + t44;
t123 = g(1) * t59 + g(2) * t58;
t122 = t64 - t65;
t116 = qJD(4) * t71;
t114 = t13 * t116 + t144 * t67;
t55 = t57 ^ 2;
t113 = t67 * t55 * t71;
t48 = pkin(2) * t59;
t112 = t48 + t125;
t8 = t53 * t117 + (qJD(2) * t88 - t69 * t118) * pkin(1);
t106 = t121 * t8;
t104 = -t49 * pkin(3) + t50 * pkin(7);
t103 = t121 * t15;
t102 = t121 * t56;
t100 = t124 + t107;
t99 = qJD(1) * (-qJD(2) + t63);
t98 = qJD(2) * (-qJD(1) - t63);
t97 = t67 * t57 * t116;
t95 = t54 + t141;
t74 = cos(qJ(1));
t92 = g(1) * t70 - g(2) * t74;
t87 = t104 - t138;
t26 = -pkin(1) * t129 + t72 * t53;
t86 = t93 - t124;
t85 = pkin(7) * t75 - t130 - t135;
t22 = -pkin(3) - t26;
t23 = pkin(7) + t27;
t84 = t22 * t56 + t23 * t75 + t131;
t82 = -t13 * t57 - t124 - t4;
t81 = -pkin(7) * qJDD(4) + (t15 - t134) * qJD(4);
t80 = -qJDD(4) * t23 + (t22 * t57 - t8) * qJD(4);
t79 = t143 * t121;
t78 = -qJDD(4) * t51 + (t52 * t57 - t143) * qJD(4);
t77 = t7 + t142;
t61 = t74 * pkin(1);
t31 = qJDD(4) * t71 - t75 * t67;
t30 = qJDD(4) * t67 + t75 * t71;
t20 = t65 * t56 - 0.2e1 * t97;
t19 = t64 * t56 + 0.2e1 * t97;
t14 = t57 * pkin(7) + t16;
t12 = -0.2e1 * t122 * t57 * qJD(4) + 0.2e1 * t67 * t127;
t1 = [0, 0, 0, 0, 0, qJDD(1), t92, g(1) * t74 + g(2) * t70, 0, 0, 0, 0, 0, 0, 0, t62, (t62 * t73 + t69 * t98) * pkin(1) + t95, ((-qJDD(1) - t62) * t69 + t73 * t98) * pkin(1) + t123, 0, (t92 + (t69 ^ 2 + t73 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), 0, 0, 0, 0, 0, t56, t26 * t56 - t131 + t77, -t27 * t56 - t8 * t57 + t86, 0, -t93 * t27 + t16 * t8 + t7 * t26 - t15 * t9 - g(1) * (-t133 - t138) - g(2) * (t48 + t61), t19, t12, t30, t20, t31, 0, t80 * t67 + (-t144 - t84) * t71 + t126, t80 * t71 + (t84 - t42) * t67 + t114, t102 * t23 + t106 * t57 + t100, t5 * t22 + t13 * t9 - g(1) * (t87 - t133) - g(2) * (t61 + t112) + t23 * t107 + t14 * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, pkin(1) * t69 * t99 + t95, (t73 * t99 - t115) * pkin(1) + t123, 0, 0, 0, 0, 0, 0, 0, t56, t24 * t57 + t18 + (pkin(2) * t56 - t96) * t72 + ((-pkin(2) * t57 - t29) * qJD(3) + (-t101 - t115) * pkin(1)) * t68 + t142, t25 * t57 + (-t117 * t57 - t56 * t68) * pkin(2) + t86, 0, t15 * t24 - t16 * t25 + (-t93 * t68 + t7 * t72 + (-t15 * t68 + t16 * t72) * qJD(3) + t141) * pkin(2), t19, t12, t30, t20, t31, 0, t78 * t67 + (-t144 + t140) * t71 + t126, t78 * t71 + (-t140 - t42) * t67 + t114, t102 * t51 + t57 * t79 + t100, -g(1) * t87 - g(2) * t112 + t107 * t51 + t13 * t94 + t14 * t79 + t5 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t77 + t130, t15 * t57 + t86, 0, 0, t19, t12, t30, t20, t31, 0, t81 * t67 + (-t144 - t85) * t71 + t126, t81 * t71 + (t85 - t42) * t67 + t114, pkin(7) * t102 - t103 * t57 + t100, -t5 * pkin(3) + pkin(7) * t107 - g(1) * t104 - g(2) * t125 - t14 * t103 - t13 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, t122 * t55, t67 * t56, t113, t127, qJDD(4), -g(3) * t71 + t67 * t82, g(3) * t67 + t71 * t82, 0, 0;];
tau_reg = t1;
