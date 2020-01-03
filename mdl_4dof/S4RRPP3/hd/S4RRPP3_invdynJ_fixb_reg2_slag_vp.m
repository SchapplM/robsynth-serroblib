% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S4RRPP3
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% tau_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRPP3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_invdynJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:54
% EndTime: 2019-12-31 16:57:55
% DurationCPUTime: 0.84s
% Computational Cost: add. (1044->212), mult. (2509->267), div. (0->0), fcn. (1674->8), ass. (0->110)
t81 = sin(qJ(1));
t83 = cos(qJ(1));
t143 = g(1) * t81 - g(2) * t83;
t106 = g(1) * t83 + g(2) * t81;
t77 = sin(pkin(6));
t78 = cos(pkin(6));
t80 = sin(qJ(2));
t82 = cos(qJ(2));
t51 = t77 * t82 + t78 * t80;
t142 = t51 * qJD(1);
t139 = t142 ^ 2;
t130 = t78 * t82;
t114 = qJD(1) * t130;
t125 = qJD(1) * t80;
t41 = t77 * t125 - t114;
t38 = t41 ^ 2;
t145 = -t38 - t139;
t144 = -t38 + t139;
t129 = qJ(3) + pkin(5);
t112 = t129 * t80;
t56 = t129 * t82;
t33 = -t77 * t112 + t78 * t56;
t73 = qJ(2) + pkin(6);
t70 = sin(t73);
t141 = -t33 * qJDD(2) - t143 * t70;
t109 = qJD(2) * t129;
t95 = -t80 * qJD(3) - t82 * t109;
t22 = qJDD(2) * pkin(2) + t95 * qJD(1) - qJDD(1) * t112;
t37 = t82 * qJD(3) - t80 * t109;
t29 = t37 * qJD(1) + qJDD(1) * t56;
t6 = t77 * t22 + t78 * t29;
t71 = cos(t73);
t140 = g(3) * t70 + t106 * t71 - t6;
t138 = pkin(2) * t80;
t134 = g(3) * t82;
t133 = t82 * pkin(2);
t132 = t142 * t41;
t54 = qJD(1) * t56;
t131 = t77 * t54;
t47 = t78 * t54;
t5 = t78 * t22 - t77 * t29;
t53 = qJD(1) * t112;
t49 = qJD(2) * pkin(2) - t53;
t25 = t77 * t49 + t47;
t75 = t80 ^ 2;
t76 = t82 ^ 2;
t128 = t75 - t76;
t127 = t75 + t76;
t126 = pkin(5) * qJDD(1);
t124 = qJD(2) * t80;
t123 = qJDD(2) * pkin(3);
t30 = -t77 * t53 + t47;
t122 = t30 * qJD(2);
t31 = -t78 * t53 - t131;
t121 = qJD(4) - t31;
t119 = t80 * qJDD(1);
t118 = t82 * qJDD(1);
t117 = qJD(1) * qJD(2);
t85 = qJD(1) ^ 2;
t116 = t80 * t85 * t82;
t115 = pkin(2) * t124;
t69 = pkin(1) + t133;
t111 = t80 * t117;
t110 = t82 * t117;
t107 = t80 * t110;
t104 = -t78 * t118 + t77 * t119;
t103 = t71 * pkin(3) + t70 * qJ(4);
t43 = t51 * qJD(2);
t26 = qJD(1) * t43 + t104;
t50 = t77 * t80 - t130;
t102 = t26 * t50 + t41 * t43;
t24 = t78 * t49 - t131;
t101 = qJD(2) * t43 + qJDD(2) * t50;
t100 = -g(3) * t71 + t106 * t70 + t5;
t99 = -0.2e1 * pkin(1) * t117 - pkin(5) * qJDD(2);
t55 = -t69 * qJD(1) + qJD(3);
t32 = t78 * t112 + t77 * t56;
t97 = -t32 * qJDD(2) + t143 * t71;
t36 = pkin(2) * t111 - t69 * qJDD(1) + qJDD(3);
t94 = t51 * qJDD(1) - t77 * t111;
t84 = qJD(2) ^ 2;
t93 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t84 + t143;
t92 = pkin(1) * t85 + t106 - t126;
t27 = t78 * t110 + t94;
t46 = qJD(2) * t130 - t77 * t124;
t91 = t142 * t43 + t51 * t26 + t27 * t50 + t46 * t41;
t12 = t41 * pkin(3) - qJ(4) * t142 + t55;
t90 = -t12 * t142 - qJDD(4) + t100;
t14 = t77 * t37 - t78 * t95;
t15 = t78 * t37 + t77 * t95;
t89 = t14 * t142 - t15 * t41 - t33 * t26 + t32 * t27 - t106;
t88 = t26 * pkin(3) - t27 * qJ(4) + t36;
t87 = 0.2e1 * t142 * qJD(2) + t104;
t74 = qJDD(2) * qJ(4);
t68 = -t78 * pkin(2) - pkin(3);
t66 = t77 * pkin(2) + qJ(4);
t61 = t83 * t69;
t28 = t46 * qJD(2) + t51 * qJDD(2);
t23 = t50 * pkin(3) - t51 * qJ(4) - t69;
t19 = qJD(2) * qJ(4) + t25;
t16 = -qJD(2) * pkin(3) + qJD(4) - t24;
t13 = pkin(2) * t125 + pkin(3) * t142 + t41 * qJ(4);
t11 = (t41 + t114) * qJD(2) + t94;
t10 = (-t41 + t114) * qJD(2) + t94;
t9 = t43 * pkin(3) - t46 * qJ(4) - t51 * qJD(4) + t115;
t4 = qJDD(4) - t123 - t5;
t3 = t142 * t46 + t27 * t51;
t2 = qJD(2) * qJD(4) + t6 + t74;
t1 = -qJD(4) * t142 + t88;
t7 = [0, 0, 0, 0, 0, qJDD(1), t143, t106, 0, 0, t75 * qJDD(1) + 0.2e1 * t107, -0.2e1 * t128 * t117 + 0.2e1 * t80 * t118, qJDD(2) * t80 + t84 * t82, t76 * qJDD(1) - 0.2e1 * t107, qJDD(2) * t82 - t84 * t80, 0, t80 * t99 + t82 * t93, -t80 * t93 + t82 * t99, 0.2e1 * t126 * t127 - t106, -g(1) * (-t81 * pkin(1) + t83 * pkin(5)) - g(2) * (t83 * pkin(1) + t81 * pkin(5)) + (t127 * pkin(5) ^ 2 + pkin(1) ^ 2) * qJDD(1), t3, -t91, t28, t102, -t101, 0, -t69 * t26 + t36 * t50 + t55 * t43 + (t41 * t138 - t14) * qJD(2) + t97, -t69 * t27 + t36 * t51 + t55 * t46 + (t138 * t142 - t15) * qJD(2) + t141, -t24 * t46 - t25 * t43 - t5 * t51 - t6 * t50 + t89, t6 * t33 + t25 * t15 - t5 * t32 - t24 * t14 - t36 * t69 + t55 * t115 - g(1) * (t129 * t83 - t81 * t69) - g(2) * (t129 * t81 + t61), t3, t28, t91, 0, t101, t102, -t14 * qJD(2) + t1 * t50 + t12 * t43 + t23 * t26 + t9 * t41 + t97, t16 * t46 - t19 * t43 - t2 * t50 + t4 * t51 + t89, t15 * qJD(2) - t1 * t51 - t12 * t46 - t142 * t9 - t23 * t27 - t141, -g(2) * t61 + t1 * t23 + t12 * t9 + t16 * t14 + t19 * t15 + t2 * t33 + t4 * t32 + (-g(1) * t129 - g(2) * t103) * t83 + (-g(1) * (-t103 - t69) - g(2) * t129) * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, t128 * t85, t119, t116, t118, qJDD(2), t80 * t92 - t134, g(3) * t80 + t82 * t92, 0, 0, t132, t144, t11, -t132, -t104, qJDD(2), t122 - t55 * t142 + (qJDD(2) * t78 - t41 * t125) * pkin(2) + t100, t31 * qJD(2) + t55 * t41 + (-qJDD(2) * t77 - t125 * t142) * pkin(2) + t140, (t25 - t30) * t142 + (-t24 + t31) * t41 + (-t26 * t77 - t27 * t78) * pkin(2), t24 * t30 - t25 * t31 + (-t134 + t5 * t78 + t6 * t77 + (-qJD(1) * t55 + t106) * t80) * pkin(2), t132, t11, -t144, qJDD(2), t104, -t132, t122 - t13 * t41 + (pkin(3) - t68) * qJDD(2) + t90, -t66 * t26 + t68 * t27 + (t19 - t30) * t142 + (t16 - t121) * t41, t66 * qJDD(2) - t12 * t41 + t13 * t142 + t74 + (0.2e1 * qJD(4) - t31) * qJD(2) - t140, t2 * t66 + t4 * t68 - t12 * t13 - t16 * t30 - g(3) * (t103 + t133) + t121 * t19 + t106 * (pkin(3) * t70 - qJ(4) * t71 + t138); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t10, t145, t142 * t24 + t25 * t41 - t143 + t36, 0, 0, 0, 0, 0, 0, t87, t145, -t10, t19 * t41 + (-qJD(4) - t16) * t142 + t88 - t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2) + t132, t11, -t139 - t84, -t19 * qJD(2) - t123 - t90;];
tau_reg = t7;
