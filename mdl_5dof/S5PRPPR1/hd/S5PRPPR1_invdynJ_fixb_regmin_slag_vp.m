% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPPR1
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:12
% EndTime: 2019-12-05 15:22:15
% DurationCPUTime: 0.86s
% Computational Cost: add. (723->173), mult. (1560->257), div. (0->0), fcn. (1185->10), ass. (0->117)
t75 = pkin(7) + qJ(2);
t67 = sin(t75);
t142 = g(1) * t67;
t69 = cos(t75);
t63 = g(2) * t69;
t146 = -t63 + t142;
t89 = qJ(3) * qJDD(2) + qJD(2) * qJD(3);
t77 = sin(pkin(8));
t125 = qJD(2) * t77;
t76 = sin(pkin(9));
t111 = t76 * t125;
t80 = sin(qJ(5));
t105 = t80 * t111;
t78 = cos(pkin(9));
t81 = cos(qJ(5));
t132 = t81 * t78;
t110 = qJD(5) * t132;
t93 = t81 * t76 + t80 * t78;
t8 = -qJD(5) * t105 + (qJD(2) * t110 + qJDD(2) * t93) * t77;
t121 = qJDD(2) * pkin(2);
t147 = t121 + t146;
t79 = cos(pkin(8));
t57 = t79 * qJD(2) - qJD(5);
t145 = qJD(5) + t57;
t143 = pkin(6) * t77;
t141 = g(3) * t77;
t112 = t77 * t132;
t118 = qJDD(2) * t76;
t85 = t93 * qJD(5);
t7 = qJDD(2) * t112 + (-qJD(2) * t85 - t80 * t118) * t77;
t140 = t7 * t79;
t139 = t79 * t8;
t28 = t93 * t125;
t138 = t28 * t57;
t30 = qJD(2) * t112 - t105;
t137 = t30 * t57;
t136 = t67 * t79;
t135 = t69 * t79;
t82 = qJD(2) ^ 2;
t134 = t79 * t82;
t133 = t80 * t76;
t32 = (-qJD(5) * t133 + t110) * t77;
t36 = t93 * t77;
t116 = t79 * qJDD(2);
t56 = -qJDD(5) + t116;
t131 = t32 * t57 + t36 * t56;
t122 = qJD(4) * t77;
t51 = -t79 * pkin(3) - t77 * qJ(4) - pkin(2);
t23 = -qJD(2) * t122 + qJDD(2) * t51 + qJDD(3);
t39 = t77 * qJDD(1) + t79 * t89;
t6 = t76 * t23 + t78 * t39;
t40 = qJD(2) * t51 + qJD(3);
t119 = qJ(3) * qJD(2);
t49 = t77 * qJD(1) + t79 * t119;
t11 = t76 * t40 + t78 * t49;
t126 = qJ(3) * t79;
t25 = t78 * t126 + t76 * t51;
t130 = t89 * t77;
t129 = t69 * pkin(2) + t67 * qJ(3);
t128 = -t76 ^ 2 - t78 ^ 2;
t71 = t77 ^ 2;
t127 = t79 ^ 2 + t71;
t124 = qJD(3) * t77;
t123 = qJD(3) * t79;
t120 = qJDD(1) - g(3);
t117 = t77 * qJDD(2);
t113 = t78 * t143;
t5 = t78 * t23 - t76 * t39;
t90 = -t79 * pkin(4) - t113;
t2 = qJDD(2) * t90 + t5;
t108 = t76 * t117;
t3 = -pkin(6) * t108 + t6;
t109 = t81 * t2 - t80 * t3;
t106 = t127 * t82;
t10 = t78 * t40 - t76 * t49;
t48 = t79 * qJD(1) - t77 * t119;
t103 = g(1) * t69 + g(2) * t67;
t101 = t80 * t2 + t81 * t3;
t4 = qJD(2) * t90 + t10;
t9 = -pkin(6) * t111 + t11;
t100 = t81 * t4 - t80 * t9;
t99 = -t80 * t4 - t81 * t9;
t38 = t79 * qJDD(1) - t130;
t45 = t78 * t51;
t12 = -t113 + t45 + (-qJ(3) * t76 - pkin(4)) * t79;
t13 = -t76 * t143 + t25;
t98 = t81 * t12 - t80 * t13;
t97 = t80 * t12 + t81 * t13;
t31 = t77 * t85;
t92 = -t132 + t133;
t37 = t92 * t77;
t96 = t31 * t57 + t37 * t56;
t95 = -t38 * t77 + t39 * t79;
t94 = t48 * t77 - t49 * t79;
t46 = qJD(4) - t48;
t35 = qJDD(4) - t38;
t88 = t92 * t57;
t24 = -t76 * t126 + t45;
t41 = -t78 * t122 - t76 * t123;
t87 = -t41 * qJD(2) - t24 * qJDD(2) - t5;
t42 = -t76 * t122 + t78 * t123;
t86 = t42 * qJD(2) + t25 * qJDD(2) + t6;
t65 = qJDD(3) - t121;
t84 = -t65 + t147;
t83 = t35 * t77 + t71 * t89 - t103;
t74 = pkin(9) + qJ(5);
t68 = cos(t74);
t66 = sin(t74);
t59 = t69 * qJ(3);
t47 = (pkin(4) * t76 + qJ(3)) * t77;
t27 = pkin(4) * t111 + t46;
t22 = t68 * t135 + t67 * t66;
t21 = -t66 * t135 + t67 * t68;
t20 = -t68 * t136 + t69 * t66;
t19 = t66 * t136 + t69 * t68;
t15 = pkin(4) * t108 + t35;
t1 = [t120, 0, 0, 0, 0, 0, 0, t38 * t79 + t39 * t77 - g(3), 0, 0, 0, -t35 * t79 - g(3) + (-t5 * t76 + t6 * t78) * t77, 0, 0, 0, 0, 0, t131 - t139, -t96 - t140; 0, qJDD(2), t146, t103, t84 * t79, -t84 * t77, t89 * t127 - t103 + t95, -t65 * pkin(2) - g(1) * (-t67 * pkin(2) + t59) - g(2) * t129 - t94 * qJD(3) + t95 * qJ(3), (t146 * t78 + t87) * t79 + t83 * t76, (-t146 * t76 + t86) * t79 + t83 * t78, (-t76 * t86 + t78 * t87 + t146) * t77, t6 * t25 + t11 * t42 + t5 * t24 + t10 * t41 - g(1) * t59 - g(2) * (pkin(3) * t135 + t129) + (t35 * qJ(3) - qJ(4) * t63 + t46 * qJD(3)) * t77 - t51 * t142, -t30 * t31 - t7 * t37, t31 * t28 - t30 * t32 - t7 * t36 + t37 * t8, t96 - t140, t131 + t139, t56 * t79, -(t81 * t41 - t80 * t42) * t57 - t98 * t56 - t109 * t79 + t28 * t124 + t47 * t8 + t15 * t36 + t27 * t32 - g(1) * t20 - g(2) * t22 + (t97 * t57 - t99 * t79) * qJD(5), (t80 * t41 + t81 * t42) * t57 + t97 * t56 + t101 * t79 + t30 * t124 + t47 * t7 - t15 * t37 - t27 * t31 - g(1) * t19 - g(2) * t21 + (t100 * t79 + t98 * t57) * qJD(5); 0, 0, 0, 0, -t116, t117, -t106, qJD(2) * t94 + qJDD(3) - t147, -t76 * t106 - t78 * t116, -t78 * t106 + t76 * t116, t128 * t117, t5 * t78 + t6 * t76 + (-t46 * t77 + (t10 * t76 - t11 * t78) * t79) * qJD(2) - t146, 0, 0, 0, 0, 0, t92 * t56 + t57 * t85 + (-t57 * t79 * t93 - t77 * t28) * qJD(2), t93 * t56 - qJD(5) * t88 + (-t77 * t30 + t79 * t88) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, (-t78 * t134 + t118) * t77, (qJDD(2) * t78 + t76 * t134) * t77, t128 * t82 * t71, qJDD(4) - t120 * t79 + ((t10 * t78 + t11 * t76) * qJD(2) - t103) * t77 + t130, 0, 0, 0, 0, 0, t8 - t137, t7 + t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 * t28, -t28 ^ 2 + t30 ^ 2, t7 - t138, -t137 - t8, -t56, -g(1) * t21 + g(2) * t19 + t66 * t141 + t145 * t99 - t27 * t30 + t109, g(1) * t22 - g(2) * t20 - t145 * t100 + t68 * t141 + t27 * t28 - t101;];
tau_reg = t1;
