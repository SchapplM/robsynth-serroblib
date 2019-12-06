% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PPRPR1
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:28
% EndTime: 2019-12-05 15:01:31
% DurationCPUTime: 1.16s
% Computational Cost: add. (1167->183), mult. (2649->243), div. (0->0), fcn. (2189->14), ass. (0->114)
t134 = qJDD(3) * pkin(3);
t89 = sin(pkin(7));
t92 = cos(pkin(7));
t116 = g(1) * t92 + g(2) * t89;
t86 = pkin(8) + qJ(3);
t79 = sin(t86);
t81 = cos(t86);
t158 = -g(3) * t81 + t116 * t79;
t88 = sin(pkin(8));
t91 = cos(pkin(8));
t95 = sin(qJ(3));
t97 = cos(qJ(3));
t52 = t97 * t88 + t95 * t91;
t44 = t52 * qJD(1);
t155 = t44 * qJD(3) + t158;
t127 = qJD(1) * qJD(3);
t120 = t95 * t127;
t130 = qJDD(1) * t97;
t157 = qJDD(1) * t95 + t97 * t127;
t118 = (t120 - t130) * t91 + t157 * t88;
t107 = qJDD(4) + t118;
t22 = t107 - t134;
t159 = t134 - t22 + t155;
t141 = t97 * t91;
t143 = t95 * t88;
t67 = qJD(1) * t143;
t41 = qJD(1) * t141 - t67;
t105 = qJD(4) - t41;
t87 = sin(pkin(9));
t94 = sin(qJ(5));
t144 = t94 * t87;
t90 = cos(pkin(9));
t96 = cos(qJ(5));
t49 = -t96 * t90 + t144;
t21 = t49 * t52;
t51 = t96 * t87 + t94 * t90;
t42 = t51 * qJD(3);
t150 = g(3) * t79;
t101 = -t116 * t81 - t150;
t153 = t42 ^ 2;
t135 = qJD(3) * t90;
t123 = t96 * t135;
t124 = qJD(3) * t144;
t39 = -t123 + t124;
t149 = t42 * t39;
t85 = pkin(9) + qJ(5);
t78 = sin(t85);
t148 = t89 * t78;
t80 = cos(t85);
t147 = t89 * t80;
t146 = t92 * t78;
t145 = t92 * t80;
t140 = pkin(6) + qJ(4);
t54 = t140 * t87;
t55 = t140 * t90;
t32 = -t96 * t54 - t94 * t55;
t139 = t32 * qJD(5) - t105 * t49;
t33 = -t94 * t54 + t96 * t55;
t138 = -t33 * qJD(5) - t105 * t51;
t125 = t88 * t130 + t157 * t91;
t19 = qJDD(3) * qJ(4) + (qJD(4) - t67) * qJD(3) + t125;
t11 = t87 * qJDD(2) + t90 * t19;
t128 = t90 * qJDD(3);
t129 = t87 * qJDD(3);
t114 = -t96 * t128 + t94 * t129;
t46 = t51 * qJD(5);
t26 = qJD(3) * t46 + t114;
t45 = t49 * qJD(5);
t137 = -t51 * t26 + t45 * t39;
t37 = qJD(3) * qJ(4) + t44;
t31 = t87 * qJD(2) + t90 * t37;
t83 = t87 ^ 2;
t84 = t90 ^ 2;
t136 = t83 + t84;
t50 = -t141 + t143;
t47 = t50 * qJD(3);
t132 = t47 * qJD(3);
t126 = qJD(5) * t123 + t94 * t128 + t96 * t129;
t73 = t90 * pkin(4) + pkin(3);
t75 = t90 * qJDD(2);
t8 = t75 + (-pkin(6) * qJDD(3) - t19) * t87;
t9 = pkin(6) * t128 + t11;
t122 = t96 * t8 - t94 * t9;
t121 = -g(1) * t89 + g(2) * t92;
t117 = t136 * qJDD(3);
t115 = t94 * t8 + t96 * t9;
t10 = -t87 * t19 + t75;
t112 = -t10 * t87 + t11 * t90;
t77 = t90 * qJD(2);
t23 = t77 + (-pkin(6) * qJD(3) - t37) * t87;
t24 = pkin(6) * t135 + t31;
t5 = t96 * t23 - t94 * t24;
t6 = t94 * t23 + t96 * t24;
t25 = qJD(5) * t124 - t126;
t111 = -t49 * t25 + t42 * t46;
t30 = -t87 * t37 + t77;
t110 = t30 * t87 - t31 * t90;
t108 = t73 * qJDD(3);
t48 = t52 * qJD(3);
t106 = t48 * qJD(3) + t50 * qJDD(3);
t20 = t51 * t52;
t99 = t107 - t158;
t53 = qJDD(2) + t121;
t38 = t39 ^ 2;
t36 = -qJD(3) * pkin(3) + t105;
t34 = -t73 * qJD(3) + t105;
t28 = -t46 * qJD(5) - t49 * qJDD(5);
t27 = -t45 * qJD(5) + t51 * qJDD(5);
t15 = -t108 + t107;
t4 = qJD(5) * t21 + t51 * t47;
t3 = -qJD(5) * t20 + t49 * t47;
t2 = -t6 * qJD(5) + t122;
t1 = t5 * qJD(5) + t115;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) + (t88 ^ 2 + t91 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -t106, -t52 * qJDD(3) + t132, 0, (-t88 * t120 + t125) * t52 - t44 * t47 + t118 * t50 - t41 * t48 - g(3), 0, 0, 0, 0, 0, 0, -t106 * t90, t106 * t87, t52 * t117 - t136 * t132, t110 * t47 + t112 * t52 + t22 * t50 + t36 * t48 - g(3), 0, 0, 0, 0, 0, 0, t4 * qJD(5) - t20 * qJDD(5) + t50 * t26 + t48 * t39, -t3 * qJD(5) + t21 * qJDD(5) - t50 * t25 + t48 * t42, -t20 * t25 + t21 * t26 - t3 * t39 - t4 * t42, -t1 * t21 + t15 * t50 - t2 * t20 + t6 * t3 + t34 * t48 + t5 * t4 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t90 + t11 * t87 + t121, 0, 0, 0, 0, 0, 0, t28, -t27, t111 + t137, t1 * t51 - t2 * t49 - t6 * t45 - t5 * t46 + t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t118 + t155, (t41 + t67) * qJD(3) - t125 - t101, 0, 0, t83 * qJDD(3), 0.2e1 * t87 * t128, 0, t84 * qJDD(3), 0, 0, t159 * t90, -t159 * t87, t105 * qJD(3) * t136 + qJ(4) * t117 + t101 + t112, -t22 * pkin(3) - t36 * t44 - g(3) * (t81 * pkin(3) + t79 * qJ(4)) + (t11 * qJ(4) + t105 * t31) * t90 + (-t10 * qJ(4) - t105 * t30) * t87 + t116 * (pkin(3) * t79 - qJ(4) * t81), -t25 * t51 - t42 * t45, -t111 + t137, t27, t26 * t49 + t39 * t46, t28, 0, t138 * qJD(5) + t32 * qJDD(5) + t15 * t49 + t158 * t80 - t73 * t26 + t34 * t46 - t44 * t39, -t139 * qJD(5) - t33 * qJDD(5) + t15 * t51 - t158 * t78 + t73 * t25 - t34 * t45 - t44 * t42, -t1 * t49 - t138 * t42 - t139 * t39 - t2 * t51 + t32 * t25 - t33 * t26 + t5 * t45 - t6 * t46 + t101, t1 * t33 + t2 * t32 - t15 * t73 - t34 * t44 - g(3) * (t140 * t79 + t81 * t73) + t139 * t6 + t138 * t5 + t116 * (-t140 * t81 + t73 * t79); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, t129, -t136 * qJD(3) ^ 2, t110 * qJD(3) - t134 + t99, 0, 0, 0, 0, 0, 0, 0.2e1 * t42 * qJD(5) + t114, (-t39 - t124) * qJD(5) + t126, -t38 - t153, t6 * t39 + t5 * t42 - t108 + t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, -t38 + t153, (t39 - t124) * qJD(5) + t126, -t149, -t114, qJDD(5), -t34 * t42 - g(1) * (-t81 * t146 + t147) - g(2) * (-t81 * t148 - t145) + t78 * t150 + t122, t34 * t39 - g(1) * (-t81 * t145 - t148) - g(2) * (-t81 * t147 + t146) + t80 * t150 - t115, 0, 0;];
tau_reg = t7;
