% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRPP4
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:59:21
% EndTime: 2019-12-31 16:59:22
% DurationCPUTime: 0.82s
% Computational Cost: add. (519->202), mult. (1177->229), div. (0->0), fcn. (585->4), ass. (0->127)
t76 = cos(qJ(2));
t129 = qJ(3) * t76;
t142 = pkin(2) + pkin(3);
t74 = sin(qJ(2));
t145 = t142 * t74;
t89 = t129 - t145;
t114 = qJD(1) * qJD(2);
t106 = t76 * t114;
t54 = t74 * qJDD(1);
t148 = t106 + t54;
t108 = t142 * qJD(2);
t117 = qJ(4) * qJD(1);
t126 = qJD(1) * t74;
t52 = pkin(5) * t126;
t21 = t74 * t117 - t52;
t118 = qJD(3) - t21;
t10 = -t108 + t118;
t104 = t142 * qJDD(2);
t130 = pkin(5) * qJD(1);
t53 = t76 * t130;
t23 = -t76 * t117 + t53;
t70 = qJD(2) * qJ(3);
t16 = t23 + t70;
t75 = sin(qJ(1));
t77 = cos(qJ(1));
t147 = g(1) * t77 + g(2) * t75;
t56 = t74 * qJ(3);
t105 = pkin(1) + t56;
t7 = qJD(4) + (t142 * t76 + t105) * qJD(1);
t128 = qJD(4) + t7;
t146 = t128 * t74;
t139 = t74 * t77;
t140 = t74 * t75;
t144 = -g(1) * t139 - g(2) * t140 + g(3) * t76;
t127 = pkin(5) * qJDD(2);
t61 = t76 * pkin(2);
t90 = -t105 - t61;
t17 = t90 * qJD(1);
t135 = t61 + t56;
t25 = -pkin(1) - t135;
t143 = (qJD(1) * t25 + t17) * qJD(2) - t127;
t71 = t74 ^ 2;
t72 = t76 ^ 2;
t132 = -t71 + t72;
t55 = t76 * qJDD(1);
t5 = 0.2e1 * t132 * t114 + 0.2e1 * t74 * t55;
t66 = g(1) * t75;
t141 = g(2) * t77;
t60 = t76 * pkin(3);
t80 = qJD(1) ^ 2;
t138 = t74 * t80;
t137 = t76 * t77;
t136 = pkin(5) - qJ(4);
t134 = t77 * pkin(1) + t75 * pkin(5);
t131 = t71 + t72;
t32 = t136 * t76;
t125 = qJD(2) * t32;
t124 = qJD(2) * t74;
t73 = qJDD(1) * pkin(1);
t123 = qJDD(2) * pkin(2);
t122 = t16 * qJD(2);
t121 = t74 * qJD(3);
t120 = t74 * qJD(4);
t119 = t76 * qJD(4);
t116 = qJ(4) * qJD(2);
t115 = qJ(4) * qJDD(1);
t49 = pkin(5) * t55;
t68 = qJDD(2) * qJ(3);
t69 = qJD(2) * qJD(3);
t113 = t49 + 0.2e1 * t68 + 0.2e1 * t69;
t112 = t49 + t68 + t69;
t111 = t60 + t135;
t40 = pkin(5) * t106;
t48 = pkin(5) * t54;
t110 = qJDD(3) + t40 + t48;
t109 = t66 - t141;
t107 = t74 * t114;
t18 = pkin(1) + t111;
t103 = qJD(1) * t18 + t7;
t102 = pkin(2) * t137 + t77 * t56 + t134;
t101 = t48 + t144;
t99 = -qJD(2) * pkin(2) + qJD(3);
t98 = t74 * t106;
t97 = t131 * qJDD(1) * pkin(5);
t96 = pkin(2) * t55 + t148 * qJ(3) + qJD(1) * t121 + t73;
t79 = qJD(2) ^ 2;
t95 = pkin(5) * t79 + t141;
t93 = pkin(2) * t74 - t129;
t24 = t52 + t99;
t29 = t53 + t70;
t92 = t24 * t76 - t29 * t74;
t91 = -qJDD(3) - t101;
t11 = t110 - t123;
t88 = -0.2e1 * pkin(1) * t114 - t127;
t87 = pkin(3) * t55 + qJDD(4) + t96;
t86 = -t95 + 0.2e1 * t73;
t1 = -t142 * t107 + t87;
t6 = t89 * qJD(2) + t121;
t85 = qJD(1) * t6 + qJDD(1) * t18 + t1 - t141;
t84 = -t74 * t115 + t40 - t91;
t15 = t93 * qJD(2) - t121;
t4 = pkin(2) * t107 - t96;
t83 = -qJD(1) * t15 - qJDD(1) * t25 - t4 - t95;
t8 = -pkin(5) * t107 + t112;
t82 = t92 * qJD(2) + t11 * t74 + t8 * t76;
t62 = t77 * pkin(5);
t44 = t76 * t66;
t43 = g(1) * t140;
t39 = t77 * t129;
t37 = t75 * t129;
t36 = t76 * t138;
t34 = qJ(4) * t107;
t33 = -t71 * t80 - t79;
t31 = t136 * t74;
t30 = -qJDD(2) - t36;
t28 = t132 * t80;
t27 = qJDD(2) * t76 - t79 * t74;
t26 = qJDD(2) * t74 + t79 * t76;
t22 = t93 * qJD(1);
t20 = t72 * qJDD(1) - 0.2e1 * t98;
t19 = t71 * qJDD(1) + 0.2e1 * t98;
t14 = -t120 + t125;
t13 = -t136 * t124 - t119;
t12 = t89 * qJD(1);
t3 = -t76 * t115 + t34 + (-pkin(5) * t124 - t119) * qJD(1) + t112;
t2 = -t148 * qJ(4) - qJD(1) * t120 - t104 + t110;
t9 = [0, 0, 0, 0, 0, qJDD(1), t109, t147, 0, 0, t19, t5, t26, t20, t27, 0, t74 * t88 + t76 * t86 + t44, -t74 * t86 + t76 * t88 - t43, -t147 + 0.2e1 * t97, -g(1) * (-t75 * pkin(1) + t62) - g(2) * t134 + (t131 * pkin(5) ^ 2 + pkin(1) ^ 2) * qJDD(1), t19, t26, -t5, 0, -t27, t20, t143 * t74 + t83 * t76 + t44, t97 + t82 - t147, -t143 * t76 + t83 * t74 + t43, pkin(5) * t82 - g(1) * t62 - g(2) * t102 + t17 * t15 + t4 * t25 - t90 * t66, t19, -t5, -t26, t20, t27, 0, -t31 * qJDD(2) + t44 + (-t103 * t74 - t14) * qJD(2) + t85 * t76, t32 * qJDD(2) + t43 + (t103 * t76 + t13) * qJD(2) + t85 * t74, (-qJD(2) * t10 - qJDD(1) * t32 - t3 + (-qJD(2) * t31 - t13) * qJD(1)) * t76 + (t122 - qJDD(1) * t31 - t2 + (-t14 + t125) * qJD(1)) * t74 + t147, t3 * t32 + t16 * t13 + t2 * t31 + t10 * t14 + t1 * t18 + t7 * t6 - g(1) * (-t77 * qJ(4) + t62) - g(2) * (pkin(3) * t137 + t102) + (-g(1) * (t90 - t60) + g(2) * qJ(4)) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t28, t54, t36, t55, qJDD(2), pkin(1) * t138 - t101, g(3) * t74 - t49 + (pkin(1) * t80 + t147) * t76, 0, 0, -t36, t54, t28, qJDD(2), -t55, t36, 0.2e1 * t123 + (-t17 * t74 + t22 * t76) * qJD(1) + t91, -t93 * qJDD(1) + ((t29 - t70) * t74 + (-t24 + t99) * t76) * qJD(1), (qJD(1) * t22 - g(3)) * t74 + (qJD(1) * t17 - t147) * t76 + t113, t8 * qJ(3) + t29 * qJD(3) - t11 * pkin(2) - t17 * t22 - g(1) * (-pkin(2) * t139 + t39) - g(2) * (-pkin(2) * t140 + t37) - g(3) * t135 - t92 * t130, -t36, t28, -t54, t36, t55, qJDD(2), t23 * qJD(2) + 0.2e1 * t104 + ((-t12 + t116) * t76 + t146) * qJD(1) - t84, -t21 * qJD(2) + t34 + (-g(3) + (-pkin(5) * qJD(2) - t12) * qJD(1)) * t74 + (-t128 * qJD(1) - t115 - t147) * t76 + t113, -t89 * qJDD(1), -g(1) * t39 - g(2) * t37 - g(3) * t111 + t3 * qJ(3) - t10 * t23 + t118 * t16 - t7 * t12 - t142 * t2 + t147 * t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t54, t33, -t29 * qJD(2) + t17 * t126 + t11 + t144, 0, 0, 0, 0, 0, 0, t30, t33, -t54, -t122 - t104 + (-t76 * t116 - t146) * qJD(1) + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 - 0.2e1 * t107, t54 + 0.2e1 * t106, -t131 * t80, (t16 * t76 + (t10 - t108) * t74) * qJD(1) + t87 + t109;];
tau_reg = t9;
