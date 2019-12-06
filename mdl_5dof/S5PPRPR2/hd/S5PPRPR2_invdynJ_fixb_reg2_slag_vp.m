% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PPRPR2
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:27
% EndTime: 2019-12-05 15:03:28
% DurationCPUTime: 0.56s
% Computational Cost: add. (649->125), mult. (1263->153), div. (0->0), fcn. (954->10), ass. (0->94)
t51 = sin(qJ(3));
t53 = cos(qJ(3));
t90 = qJD(1) * qJD(3);
t126 = qJDD(1) * t51 + t53 * t90;
t47 = sin(pkin(7));
t49 = cos(pkin(7));
t125 = g(1) * t47 - g(2) * t49;
t48 = cos(pkin(8));
t109 = t53 * t48;
t46 = sin(pkin(8));
t110 = t51 * t46;
t85 = qJD(1) * t110;
t12 = -qJD(1) * t109 + t85;
t41 = pkin(8) + qJ(3);
t38 = sin(t41);
t39 = cos(t41);
t77 = g(1) * t49 + g(2) * t47;
t60 = -g(3) * t38 - t77 * t39;
t96 = qJDD(1) * t53;
t87 = t126 * t48 + t46 * t96;
t124 = (-t12 + t85) * qJD(3) - t87 - t60;
t123 = t77 * t38;
t64 = -qJDD(2) + t125;
t16 = -t109 + t110;
t17 = t53 * t46 + t51 * t48;
t13 = t17 * qJD(1);
t98 = qJD(3) * qJ(4);
t11 = t13 + t98;
t101 = t11 * qJD(3);
t36 = g(3) * t39;
t122 = -t101 + t36 - t123;
t118 = pkin(3) + pkin(6);
t120 = (t11 - t13 + t98) * qJD(5) - qJDD(5) * t118;
t50 = sin(qJ(5));
t52 = cos(qJ(5));
t65 = qJD(4) + t12;
t9 = -t118 * qJD(3) + t65;
t7 = -t50 * qJD(2) + t52 * t9;
t103 = t7 * qJD(5);
t84 = t51 * t90;
t82 = (t84 - t96) * t48 + t126 * t46;
t71 = qJDD(4) + t82;
t4 = -t118 * qJDD(3) + t71;
t1 = t52 * qJDD(2) + t50 * t4 + t103;
t3 = t52 * t4;
t8 = t52 * qJD(2) + t50 * t9;
t2 = -t8 * qJD(5) - t50 * qJDD(2) + t3;
t58 = -(t50 * t7 - t52 * t8) * qJD(5) + t1 * t50 + t2 * t52;
t15 = qJD(3) * t17;
t119 = 0.2e1 * qJD(5) * t15 + qJDD(5) * t16;
t113 = t38 * t47;
t112 = t38 * t49;
t111 = t50 * t52;
t108 = t39 * pkin(3) + t38 * qJ(4);
t44 = t50 ^ 2;
t45 = t52 ^ 2;
t107 = t44 - t45;
t106 = t44 + t45;
t55 = qJD(5) ^ 2;
t56 = qJD(3) ^ 2;
t105 = -t55 - t56;
t104 = qJ(4) * t39;
t102 = qJDD(3) * pkin(3);
t100 = t13 * qJD(3);
t99 = t15 * qJD(3);
t94 = qJDD(5) * t50;
t93 = qJDD(5) * t52;
t91 = t52 * qJDD(3);
t89 = qJD(3) * qJD(5);
t42 = qJDD(3) * qJ(4);
t88 = t56 * t111;
t86 = -g(1) * t112 - g(2) * t113 + t36;
t79 = t106 * qJDD(3);
t78 = t89 * t111;
t74 = t8 * t50 + t7 * t52;
t14 = t16 * qJD(3);
t43 = qJD(3) * qJD(4);
t62 = -t46 * t84 + t87;
t5 = -t42 - t43 - t62;
t69 = -t11 * t14 - t5 * t17 - g(3);
t68 = t86 - t101;
t67 = t14 * qJD(3) - t17 * qJDD(3);
t66 = t16 * qJDD(3) + t99;
t63 = -t5 * qJ(4) + t65 * t11;
t6 = t71 - t102;
t61 = -t16 * t55 - t67;
t59 = -t82 - t86 + t100;
t57 = t65 * qJD(3) + t118 * t55 + t42 - t5 + t60;
t22 = t49 * t104;
t21 = t47 * t104;
t20 = -t55 * t50 + t93;
t19 = -t55 * t52 - t94;
t10 = -qJD(3) * pkin(3) + t65;
t18 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) + (t46 ^ 2 + t48 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -t66, t67, 0, t12 * t15 - t13 * t14 + t82 * t16 + t62 * t17 - g(3), 0, 0, 0, 0, 0, 0, 0, t66, -t67, t10 * t15 + t6 * t16 + t69, 0, 0, 0, 0, 0, 0, t119 * t52 + t61 * t50, -t119 * t50 + t61 * t52, -t106 * t99 - t16 * t79, t74 * t15 + t58 * t16 + t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, 0, 0, 0, 0, 0, 0, t19, -t20, 0, -t74 * qJD(5) + t1 * t52 - t2 * t50 - t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t59, t124, 0, 0, qJDD(3), 0, 0, 0, 0, 0, 0, qJDD(4) - t59 - 0.2e1 * t102, 0.2e1 * t42 + 0.2e1 * t43 - t124, -t6 * pkin(3) - t10 * t13 - g(1) * (-pkin(3) * t112 + t22) - g(2) * (-pkin(3) * t113 + t21) - g(3) * t108 + t63, t45 * qJDD(3) - 0.2e1 * t78, 0.2e1 * t107 * t89 - 0.2e1 * t50 * t91, t20, t44 * qJDD(3) + 0.2e1 * t78, t19, 0, t120 * t52 + t57 * t50, -t120 * t50 + t57 * t52, t106 * t100 + t118 * t79 - t58 - t86, -g(1) * t22 - g(2) * t21 - g(3) * (t39 * pkin(6) + t108) - t74 * t13 + t63 + (-t58 + t123) * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t56, t6 + t68, 0, 0, 0, 0, 0, 0, t105 * t50 + t93, t105 * t52 - t94, -t79, t58 + t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, -t107 * t56, t91, -t88, -t50 * qJDD(3), qJDD(5), t122 * t52 + t64 * t50 + t3, t103 + (-qJD(5) * t9 + t64) * t52 + (qJD(5) * qJD(2) - t122 - t4) * t50, 0, 0;];
tau_reg = t18;
