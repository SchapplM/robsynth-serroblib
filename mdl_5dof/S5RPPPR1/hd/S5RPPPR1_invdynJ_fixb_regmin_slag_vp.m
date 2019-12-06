% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:29:09
% EndTime: 2019-12-05 17:29:12
% DurationCPUTime: 0.89s
% Computational Cost: add. (865->181), mult. (1783->271), div. (0->0), fcn. (1301->14), ass. (0->125)
t77 = sin(pkin(8));
t129 = qJD(1) * t77;
t76 = sin(pkin(9));
t116 = t76 * t129;
t82 = sin(qJ(5));
t110 = t82 * t116;
t79 = cos(pkin(9));
t84 = cos(qJ(5));
t133 = t84 * t79;
t115 = qJD(5) * t133;
t97 = t84 * t76 + t82 * t79;
t11 = -qJD(5) * t110 + (qJD(1) * t115 + t97 * qJDD(1)) * t77;
t75 = qJ(1) + pkin(7);
t69 = cos(t75);
t145 = g(2) * t69;
t67 = sin(t75);
t108 = g(3) * t67 + t145;
t81 = cos(pkin(7));
t144 = t81 * pkin(1);
t65 = -pkin(2) - t144;
t123 = qJDD(1) * t65;
t56 = qJDD(3) + t123;
t150 = t108 - t56;
t80 = cos(pkin(8));
t62 = t80 * qJD(1) - qJD(5);
t149 = qJD(5) + t62;
t78 = sin(pkin(7));
t64 = t78 * pkin(1) + qJ(3);
t52 = qJD(1) * qJD(3) + qJDD(1) * t64;
t147 = pkin(6) * t77;
t146 = g(1) * t77;
t85 = cos(qJ(1));
t143 = t85 * pkin(1);
t117 = t77 * t133;
t122 = qJDD(1) * t76;
t90 = t97 * qJD(5);
t10 = qJDD(1) * t117 + (-qJD(1) * t90 - t82 * t122) * t77;
t142 = t10 * t80;
t35 = t97 * t129;
t141 = t35 * t62;
t37 = qJD(1) * t117 - t110;
t140 = t37 * t62;
t139 = t64 * t80;
t138 = t67 * t80;
t137 = t69 * t80;
t136 = t80 * t11;
t86 = qJD(1) ^ 2;
t135 = t80 * t86;
t134 = t82 * t76;
t39 = (-qJD(5) * t134 + t115) * t77;
t43 = t97 * t77;
t120 = t80 * qJDD(1);
t61 = -qJDD(5) + t120;
t132 = t39 * t62 + t43 * t61;
t126 = qJD(4) * t77;
t94 = t80 * pkin(3) + t77 * qJ(4) + pkin(2);
t51 = -t94 - t144;
t19 = -qJD(1) * t126 + qJDD(1) * t51 + qJDD(3);
t32 = t77 * qJDD(2) + t80 * t52;
t6 = t76 * t19 + t79 * t32;
t33 = qJD(1) * t51 + qJD(3);
t60 = t64 * qJD(1);
t42 = t77 * qJD(2) + t80 * t60;
t9 = t76 * t33 + t79 * t42;
t18 = t79 * t139 + t76 * t51;
t131 = -t76 ^ 2 - t79 ^ 2;
t71 = t77 ^ 2;
t130 = t80 ^ 2 + t71;
t128 = qJD(3) * t77;
t127 = qJD(3) * t80;
t125 = qJDD(2) - g(1);
t121 = t77 * qJDD(1);
t118 = t79 * t147;
t5 = t79 * t19 - t76 * t32;
t95 = -t80 * pkin(4) - t118;
t2 = qJDD(1) * t95 + t5;
t113 = t76 * t121;
t3 = -pkin(6) * t113 + t6;
t114 = t84 * t2 - t82 * t3;
t83 = sin(qJ(1));
t112 = -t83 * pkin(1) + t69 * qJ(3);
t111 = t130 * t86;
t8 = t79 * t33 - t76 * t42;
t41 = t80 * qJD(2) - t77 * t60;
t48 = t77 * t52;
t31 = t80 * qJDD(2) - t48;
t107 = g(2) * t67 - g(3) * t69;
t106 = g(2) * t85 + g(3) * t83;
t105 = t82 * t2 + t84 * t3;
t4 = qJD(1) * t95 + t8;
t7 = -pkin(6) * t116 + t9;
t104 = t84 * t4 - t82 * t7;
t103 = -t82 * t4 - t84 * t7;
t46 = t79 * t51;
t12 = -t118 + t46 + (-t64 * t76 - pkin(4)) * t80;
t13 = -t76 * t147 + t18;
t102 = t84 * t12 - t82 * t13;
t101 = t82 * t12 + t84 * t13;
t100 = -t31 * t77 + t32 * t80;
t38 = t77 * t90;
t96 = -t133 + t134;
t44 = t96 * t77;
t99 = t38 * t62 + t44 * t61;
t98 = t41 * t77 - t42 * t80;
t40 = qJD(4) - t41;
t26 = qJDD(4) - t31;
t93 = t96 * t62;
t17 = -t76 * t139 + t46;
t49 = -t79 * t126 - t76 * t127;
t92 = -t49 * qJD(1) - t17 * qJDD(1) - t5;
t50 = -t76 * t126 + t79 * t127;
t91 = t50 * qJD(1) + t18 * qJDD(1) + t6;
t89 = -t123 + t150;
t88 = t26 * t77 + t52 * t71 + t107;
t74 = pkin(9) + qJ(5);
t68 = cos(t74);
t66 = sin(t74);
t47 = (pkin(4) * t76 + t64) * t77;
t30 = -t68 * t137 - t67 * t66;
t29 = t66 * t137 - t67 * t68;
t28 = t68 * t138 - t69 * t66;
t27 = t66 * t138 + t69 * t68;
t21 = pkin(4) * t116 + t40;
t16 = pkin(4) * t113 + t26;
t1 = [qJDD(1), t106, -g(2) * t83 + g(3) * t85, (t106 + (t78 ^ 2 + t81 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t89 * t80, -t89 * t77, t52 * t130 + t100 + t107, t56 * t65 - g(2) * (-t69 * pkin(2) - t67 * qJ(3) - t143) - g(3) * (-t67 * pkin(2) + t112) + t100 * t64 - t98 * qJD(3), (t108 * t79 + t92) * t80 + t88 * t76, (-t108 * t76 + t91) * t80 + t88 * t79, (-t76 * t91 + t79 * t92 + t108) * t77, t6 * t18 + t9 * t50 + t5 * t17 + t8 * t49 + g(2) * t143 - g(3) * t112 + (t40 * qJD(3) + t26 * t64) * t77 + t94 * t145 + (g(2) * qJ(3) + g(3) * t94) * t67, -t10 * t44 - t37 * t38, -t10 * t43 + t44 * t11 + t38 * t35 - t37 * t39, t99 - t142, t132 + t136, t61 * t80, -(t84 * t49 - t82 * t50) * t62 - t102 * t61 - t114 * t80 + t35 * t128 + t47 * t11 + t16 * t43 + t21 * t39 - g(2) * t30 + g(3) * t28 + (t101 * t62 - t103 * t80) * qJD(5), (t82 * t49 + t84 * t50) * t62 + t101 * t61 + t105 * t80 + t37 * t128 + t47 * t10 - t16 * t44 - t21 * t38 - g(2) * t29 - g(3) * t27 + (t102 * t62 + t104 * t80) * qJD(5); 0, 0, 0, t125, 0, 0, 0, t31 * t80 + t32 * t77 - g(1), 0, 0, 0, -t26 * t80 - g(1) + (-t5 * t76 + t6 * t79) * t77, 0, 0, 0, 0, 0, t132 - t136, -t99 - t142; 0, 0, 0, 0, -t120, t121, -t111, t98 * qJD(1) - t150, -t76 * t111 - t79 * t120, -t79 * t111 + t76 * t120, t131 * t121, t5 * t79 + t6 * t76 + (-t40 * t77 + (t76 * t8 - t79 * t9) * t80) * qJD(1) - t108, 0, 0, 0, 0, 0, t96 * t61 + t62 * t90 + (-t97 * t62 * t80 - t77 * t35) * qJD(1), t97 * t61 - qJD(5) * t93 + (-t77 * t37 + t80 * t93) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, (-t79 * t135 + t122) * t77, (qJDD(1) * t79 + t76 * t135) * t77, t131 * t86 * t71, qJDD(4) + t48 - t125 * t80 + ((t76 * t9 + t79 * t8) * qJD(1) + t107) * t77, 0, 0, 0, 0, 0, t11 - t140, t10 + t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t35, -t35 ^ 2 + t37 ^ 2, t10 - t141, -t11 - t140, -t61, -g(2) * t27 + g(3) * t29 + t149 * t103 + t66 * t146 - t21 * t37 + t114, -g(2) * t28 - g(3) * t30 - t149 * t104 + t68 * t146 + t21 * t35 - t105;];
tau_reg = t1;
