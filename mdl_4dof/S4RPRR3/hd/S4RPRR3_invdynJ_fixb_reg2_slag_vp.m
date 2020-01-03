% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RPRR3
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
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RPRR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:21
% EndTime: 2019-12-31 16:49:22
% DurationCPUTime: 0.80s
% Computational Cost: add. (1180->185), mult. (2531->249), div. (0->0), fcn. (1620->12), ass. (0->118)
t83 = sin(pkin(7));
t63 = t83 * pkin(1) + pkin(5);
t51 = t63 * qJDD(1);
t146 = qJD(2) * qJD(3) + t51;
t79 = qJ(1) + pkin(7);
t69 = sin(t79);
t70 = cos(t79);
t107 = g(1) * t70 + g(2) * t69;
t130 = pkin(1) * qJDD(1);
t85 = sin(qJ(4));
t86 = sin(qJ(3));
t88 = cos(qJ(4));
t89 = cos(qJ(3));
t43 = t85 * t89 + t88 * t86;
t37 = t43 * qJD(1);
t78 = qJD(3) + qJD(4);
t53 = t63 * qJD(1);
t111 = pkin(6) * qJD(1) + t53;
t125 = t86 * qJD(2);
t29 = t111 * t89 + t125;
t72 = t89 * qJDD(2);
t11 = qJDD(3) * pkin(3) + t72 + (-pkin(6) * qJDD(1) - t51) * t86 - t29 * qJD(3);
t127 = qJD(4) * t85;
t122 = qJD(1) * qJD(3);
t115 = t86 * t122;
t123 = t89 * qJDD(1);
t119 = t86 * qJDD(2) + t146 * t89;
t128 = qJD(3) * t86;
t20 = -t53 * t128 + t119;
t14 = (-t115 + t123) * pkin(6) + t20;
t73 = t89 * qJD(2);
t28 = -t111 * t86 + t73;
t27 = qJD(3) * pkin(3) + t28;
t1 = (qJD(4) * t27 + t14) * t88 + t85 * t11 - t29 * t127;
t143 = g(3) * t89;
t84 = cos(pkin(7));
t142 = t84 * pkin(1);
t87 = sin(qJ(1));
t141 = t87 * pkin(1);
t140 = pkin(6) + t63;
t124 = t86 * qJDD(1);
t104 = -t88 * t123 + t85 * t124;
t25 = t78 * t43;
t17 = qJD(1) * t25 + t104;
t136 = t85 * t86;
t103 = t78 * t136;
t126 = qJD(4) * t88;
t133 = t88 * t89;
t24 = -qJD(3) * t133 - t89 * t126 + t103;
t116 = qJD(1) * t133;
t129 = qJD(1) * t86;
t117 = t85 * t129;
t35 = -t116 + t117;
t139 = -t43 * t17 + t24 * t35;
t138 = t37 * t35;
t137 = t85 * t29;
t135 = t86 * t53;
t134 = t88 * t29;
t132 = t89 * t53;
t80 = t86 ^ 2;
t81 = t89 ^ 2;
t131 = t80 - t81;
t64 = -pkin(2) - t142;
t54 = qJD(1) * t64;
t52 = qJDD(1) * t64;
t93 = qJD(1) ^ 2;
t120 = t86 * t93 * t89;
t118 = pkin(3) * t128;
t68 = t89 * pkin(3) + pkin(2);
t114 = t89 * t122;
t112 = qJD(3) * t140;
t110 = -qJD(4) * t116 - t85 * t123 + (-t114 - t124) * t88;
t108 = t86 * t114;
t106 = g(1) * t69 - g(2) * t70;
t90 = cos(qJ(1));
t105 = g(1) * t87 - g(2) * t90;
t16 = qJD(1) * t103 + t110;
t42 = -t133 + t136;
t102 = -t42 * t16 + t37 * t25;
t77 = qJDD(3) + qJDD(4);
t101 = t24 * t78 - t43 * t77;
t10 = t85 * t27 + t134;
t39 = t140 * t86;
t40 = t140 * t89;
t22 = -t88 * t39 - t85 * t40;
t23 = -t85 * t39 + t88 * t40;
t32 = t125 + t132;
t48 = -t68 - t142;
t100 = -qJD(1) * t54 + t107;
t99 = 0.2e1 * t54 * qJD(3) - qJDD(3) * t63;
t2 = -qJD(4) * t10 + t88 * t11 - t85 * t14;
t92 = qJD(3) ^ 2;
t98 = -t63 * t92 + t106 - 0.2e1 * t52;
t21 = -t32 * qJD(3) - t86 * t51 + t72;
t31 = t73 - t135;
t97 = t20 * t89 - t21 * t86 + (-t31 * t89 - t32 * t86) * qJD(3);
t38 = t48 * qJD(1);
t82 = qJ(3) + qJ(4);
t74 = sin(t82);
t75 = cos(t82);
t96 = g(3) * t74 + t107 * t75 + t38 * t35 - t1;
t95 = -g(3) * t75 + t107 * t74 - t38 * t37 + t2;
t91 = -pkin(6) - pkin(5);
t76 = t90 * pkin(1);
t50 = qJDD(3) * t89 - t92 * t86;
t49 = qJDD(3) * t86 + t92 * t89;
t34 = t89 * t112;
t33 = t86 * t112;
t30 = pkin(3) * t115 + qJDD(1) * t48;
t18 = -t35 ^ 2 + t37 ^ 2;
t15 = -t25 * t78 - t42 * t77;
t13 = t88 * t28 - t137;
t12 = -t85 * t28 - t134;
t9 = t88 * t27 - t137;
t5 = -t110 + (-t117 + t35) * t78;
t4 = -qJD(4) * t23 + t85 * t33 - t88 * t34;
t3 = qJD(4) * t22 - t88 * t33 - t85 * t34;
t6 = [0, 0, 0, 0, 0, qJDD(1), t105, g(1) * t90 + g(2) * t87, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t84 * t130 + t106, -0.2e1 * t83 * t130 + t107, 0, (t105 + (t83 ^ 2 + t84 ^ 2) * t130) * pkin(1), t80 * qJDD(1) + 0.2e1 * t108, -0.2e1 * t131 * t122 + 0.2e1 * t86 * t123, t49, t81 * qJDD(1) - 0.2e1 * t108, t50, 0, t86 * t99 + t89 * t98, -t86 * t98 + t89 * t99, (t80 + t81) * t51 + t97 - t107, t52 * t64 - g(1) * (-t69 * pkin(2) + t70 * pkin(5) - t141) - g(2) * (t70 * pkin(2) + t69 * pkin(5) + t76) + t97 * t63, -t16 * t43 - t37 * t24, -t102 + t139, -t101, t17 * t42 + t35 * t25, t15, 0, t106 * t75 + t35 * t118 + t48 * t17 + t22 * t77 + t38 * t25 + t30 * t42 + t4 * t78, -t106 * t74 + t37 * t118 - t48 * t16 - t23 * t77 - t38 * t24 - t3 * t78 + t30 * t43, -t1 * t42 - t10 * t25 + t22 * t16 - t23 * t17 - t2 * t43 + t9 * t24 - t3 * t35 - t4 * t37 - t107, t1 * t23 + t10 * t3 + t2 * t22 + t9 * t4 + t30 * t48 + t38 * t118 - g(1) * (-t69 * t68 - t70 * t91 - t141) - g(2) * (t70 * t68 - t69 * t91 + t76); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t50, -t49, 0, t20 * t86 + t21 * t89 - g(3) + (-t31 * t86 + t32 * t89) * qJD(3), 0, 0, 0, 0, 0, 0, t15, t101, t102 + t139, t1 * t43 - t10 * t24 - t2 * t42 - t9 * t25 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, t131 * t93, t124, t120, t123, qJDD(3), -t143 + t72 + (t32 - t132) * qJD(3) + (t100 - t146) * t86, g(3) * t86 + (t31 + t135) * qJD(3) + t100 * t89 - t119, 0, 0, t138, t18, t5, -t138, -t104, t77, -t12 * t78 + (-t78 * t127 - t35 * t129 + t77 * t88) * pkin(3) + t95, t13 * t78 + (-t78 * t126 - t37 * t129 - t77 * t85) * pkin(3) + t96, (t10 + t12) * t37 + (t13 - t9) * t35 + (t16 * t88 - t17 * t85 + (-t35 * t88 + t37 * t85) * qJD(4)) * pkin(3), -t10 * t13 - t9 * t12 + (-t143 + t1 * t85 + t2 * t88 + (t10 * t88 - t85 * t9) * qJD(4) + (-qJD(1) * t38 + t107) * t86) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, t18, t5, -t138, -t104, t77, t10 * t78 + t95, t9 * t78 + t96, 0, 0;];
tau_reg = t6;
