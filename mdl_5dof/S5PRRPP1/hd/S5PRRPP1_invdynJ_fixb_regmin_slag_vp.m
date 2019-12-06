% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRPP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% tau_reg [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:56
% EndTime: 2019-12-05 16:06:59
% DurationCPUTime: 0.69s
% Computational Cost: add. (1128->193), mult. (2377->246), div. (0->0), fcn. (1547->8), ass. (0->97)
t77 = sin(pkin(8));
t78 = cos(pkin(8));
t80 = sin(qJ(3));
t81 = cos(qJ(3));
t49 = t77 * t81 + t78 * t80;
t45 = t49 * qJD(2);
t41 = t45 ^ 2;
t120 = t78 * t81;
t107 = qJD(2) * t120;
t117 = qJD(2) * t80;
t42 = t117 * t77 - t107;
t130 = -t42 ^ 2 - t41;
t72 = pkin(7) + qJ(2);
t66 = sin(t72);
t68 = cos(t72);
t129 = g(1) * t66 - g(2) * t68;
t100 = g(1) * t68 + g(2) * t66;
t119 = qJ(4) + pkin(6);
t127 = qJD(3) * qJD(1) + qJD(2) * qJD(4) + t119 * qJDD(2);
t70 = t81 * qJDD(1);
t101 = t119 * qJD(3);
t97 = qJD(2) * t101;
t16 = qJDD(3) * pkin(3) - t127 * t80 - t81 * t97 + t70;
t20 = (qJDD(1) - t97) * t80 + t127 * t81;
t4 = t78 * t16 - t77 * t20;
t106 = -qJDD(5) + t4;
t122 = t81 * pkin(3);
t65 = pkin(2) + t122;
t52 = -qJD(2) * t65 + qJD(4);
t17 = t42 * pkin(4) - t45 * qJ(5) + t52;
t73 = qJ(3) + pkin(8);
t67 = sin(t73);
t69 = cos(t73);
t128 = -g(3) * t69 + t100 * t67 - t17 * t45 + t106;
t110 = qJD(2) * qJD(3);
t103 = t80 * t110;
t88 = t49 * qJDD(2) - t77 * t103;
t123 = g(3) * t81;
t5 = t77 * t16 + t78 * t20;
t55 = t119 * t81;
t39 = t80 * qJD(1) + qJD(2) * t55;
t121 = t77 * t39;
t33 = t78 * t39;
t104 = t119 * t80;
t37 = t81 * qJD(1) - qJD(2) * t104;
t36 = qJD(3) * pkin(3) + t37;
t19 = t77 * t36 + t33;
t75 = t80 ^ 2;
t118 = -t81 ^ 2 + t75;
t116 = qJDD(3) * pkin(4);
t115 = t80 * qJD(3);
t24 = t78 * t37 - t121;
t114 = qJD(5) - t24;
t113 = qJDD(1) - g(3);
t112 = t80 * qJDD(2);
t111 = t81 * qJDD(2);
t109 = qJDD(3) * qJ(5) + t5;
t108 = pkin(3) * t115;
t102 = t81 * t110;
t98 = -t78 * t111 + t112 * t77;
t96 = t69 * pkin(4) + t67 * qJ(5);
t18 = t78 * t36 - t121;
t94 = -0.2e1 * pkin(2) * t110 - pkin(6) * qJDD(3);
t44 = t49 * qJD(3);
t27 = qJD(2) * t44 + t98;
t28 = t102 * t78 + t88;
t47 = qJD(3) * t120 - t115 * t77;
t48 = t77 * t80 - t120;
t93 = -t49 * t27 + t48 * t28 - t47 * t42 + t44 * t45;
t91 = -t80 * qJD(4) - t101 * t81;
t90 = pkin(3) * t103 - qJDD(2) * t65 + qJDD(4);
t82 = qJD(3) ^ 2;
t87 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t82 + t129;
t83 = qJD(2) ^ 2;
t86 = t83 * pkin(2) - qJDD(2) * pkin(6) + t100;
t38 = t81 * qJD(4) - t101 * t80;
t23 = t77 * t38 - t78 * t91;
t25 = t78 * t38 + t77 * t91;
t30 = t104 * t78 + t77 * t55;
t31 = -t104 * t77 + t78 * t55;
t85 = t23 * t45 - t25 * t42 - t31 * t27 + t30 * t28 - t100;
t84 = t27 * pkin(4) - t28 * qJ(5) + t90;
t64 = -t78 * pkin(3) - pkin(4);
t61 = t77 * pkin(3) + qJ(5);
t54 = qJDD(3) * t81 - t82 * t80;
t53 = qJDD(3) * t80 + t82 * t81;
t51 = t68 * t65;
t26 = t48 * pkin(4) - t49 * qJ(5) - t65;
t22 = t77 * t37 + t33;
t21 = pkin(3) * t117 + t45 * pkin(4) + t42 * qJ(5);
t13 = qJD(3) * qJ(5) + t19;
t11 = -qJD(3) * pkin(4) + qJD(5) - t18;
t8 = t44 * pkin(4) - t47 * qJ(5) - t49 * qJD(5) + t108;
t3 = -t106 - t116;
t2 = qJD(3) * qJD(5) + t109;
t1 = -t45 * qJD(5) + t84;
t6 = [t113, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t53, t93, -t18 * t44 + t19 * t47 - t4 * t48 + t5 * t49 - g(3), -t44 * qJD(3) - t48 * qJDD(3), t93, t47 * qJD(3) + t49 * qJDD(3), t11 * t44 + t13 * t47 + t2 * t49 + t3 * t48 - g(3); 0, qJDD(2), t129, t100, t75 * qJDD(2) + 0.2e1 * t102 * t80, -0.2e1 * t110 * t118 + 0.2e1 * t111 * t80, t53, t54, 0, t80 * t94 + t81 * t87, -t80 * t87 + t81 * t94, -t18 * t47 - t19 * t44 - t4 * t49 - t5 * t48 + t85, t5 * t31 + t19 * t25 - t4 * t30 - t18 * t23 - t90 * t65 + t52 * t108 - g(1) * (t119 * t68 - t66 * t65) - g(2) * (t119 * t66 + t51), -t23 * qJD(3) - t30 * qJDD(3) + t1 * t48 + t129 * t69 + t17 * t44 + t26 * t27 + t8 * t42, t11 * t47 - t13 * t44 - t2 * t48 + t3 * t49 + t85, t25 * qJD(3) + t31 * qJDD(3) - t1 * t49 + t129 * t67 - t17 * t47 - t26 * t28 - t8 * t45, -g(2) * t51 + t1 * t26 + t11 * t23 + t13 * t25 + t17 * t8 + t2 * t31 + t3 * t30 + (-g(1) * t119 - g(2) * t96) * t68 + (-g(1) * (-t65 - t96) - g(2) * t119) * t66; 0, 0, 0, 0, -t80 * t83 * t81, t118 * t83, t112, t111, qJDD(3), t80 * t86 - t123 + t70, -t113 * t80 + t81 * t86, (t19 - t22) * t45 + (-t18 + t24) * t42 + (-t27 * t77 - t28 * t78) * pkin(3), t18 * t22 - t19 * t24 + (-t123 + t4 * t78 + t5 * t77 + (-qJD(2) * t52 + t100) * t80) * pkin(3), t22 * qJD(3) - t21 * t42 + (pkin(4) - t64) * qJDD(3) + t128, -t61 * t27 + t64 * t28 + (t13 - t22) * t45 + (t11 - t114) * t42, -g(3) * t67 + t61 * qJDD(3) - t17 * t42 + t21 * t45 - t100 * t69 + (0.2e1 * qJD(5) - t24) * qJD(3) + t109, t2 * t61 + t3 * t64 - t17 * t21 - t11 * t22 - g(3) * (t96 + t122) + t114 * t13 + t100 * (pkin(3) * t80 + pkin(4) * t67 - qJ(5) * t69); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t18 * t45 + t19 * t42 - t129 + t90, 0.2e1 * t45 * qJD(3) + t98, t130, (t42 - t107) * qJD(3) - t88, t13 * t42 + (-qJD(5) - t11) * t45 + t84 - t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t42 - qJDD(3), (t42 + t107) * qJD(3) + t88, -t41 - t82, -t13 * qJD(3) - t116 - t128;];
tau_reg = t6;
