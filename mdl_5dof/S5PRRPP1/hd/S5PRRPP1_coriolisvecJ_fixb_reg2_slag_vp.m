% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:53
% EndTime: 2019-12-05 16:06:57
% DurationCPUTime: 0.78s
% Computational Cost: add. (929->163), mult. (2474->215), div. (0->0), fcn. (1617->4), ass. (0->105)
t74 = sin(pkin(8));
t76 = cos(qJ(3));
t108 = cos(pkin(8));
t75 = sin(qJ(3));
t94 = t108 * t75;
t56 = t74 * t76 + t94;
t49 = t56 * qJD(3);
t37 = qJD(2) * t49;
t106 = qJD(2) * t75;
t93 = t108 * t76;
t89 = qJD(2) * t93;
t47 = t74 * t106 - t89;
t105 = qJD(3) * t75;
t52 = qJD(3) * t93 - t74 * t105;
t110 = -t56 * t37 - t52 * t47;
t101 = qJD(2) * qJD(3);
t95 = t75 * t101;
t62 = t74 * t95;
t38 = qJD(3) * t89 - t62;
t55 = t74 * t75 - t93;
t83 = qJD(2) * t56;
t111 = t38 * t55 + t49 * t83;
t129 = t110 + t111;
t128 = t111 - t110;
t127 = -0.2e1 * t101;
t126 = 0.2e1 * t83;
t123 = t83 ^ 2;
t44 = t47 ^ 2;
t125 = -t44 - t123;
t124 = -t44 + t123;
t122 = pkin(3) * t75;
t112 = -qJ(4) - pkin(6);
t61 = t112 * t76;
t30 = -t112 * t94 - t74 * t61;
t92 = qJD(3) * t112;
t42 = t76 * qJD(4) + t75 * t92;
t102 = qJD(1) * qJD(3);
t69 = t76 * t102;
t29 = qJD(2) * t42 + t69;
t84 = -t75 * qJD(4) + t76 * t92;
t79 = qJD(2) * t84 - t75 * t102;
t4 = -t108 * t79 + t74 * t29;
t121 = t4 * t30;
t120 = t4 * t55;
t70 = -t76 * pkin(3) - pkin(2);
t107 = qJD(2) * t70;
t60 = qJD(4) + t107;
t11 = t47 * pkin(4) - qJ(5) * t83 + t60;
t119 = t11 * t83;
t118 = t83 * t47;
t104 = t75 * qJD(1);
t43 = -qJD(2) * t61 + t104;
t116 = t74 * t43;
t78 = qJD(2) ^ 2;
t115 = t76 * t78;
t77 = qJD(3) ^ 2;
t114 = t77 * t75;
t113 = t77 * t76;
t5 = t108 * t29 + t74 * t79;
t33 = t108 * t43;
t71 = t76 * qJD(1);
t96 = t112 * t75;
t41 = qJD(2) * t96 + t71;
t36 = qJD(3) * pkin(3) + t41;
t13 = t74 * t36 + t33;
t109 = t75 ^ 2 - t76 ^ 2;
t39 = qJD(3) * t49;
t17 = t108 * t41 - t116;
t103 = qJD(5) - t17;
t100 = t75 * t115;
t99 = pkin(3) * t105;
t98 = pkin(3) * t106;
t97 = pkin(6) * t106;
t91 = pkin(2) * t127;
t90 = t76 * t95;
t88 = t37 * t55 + t47 * t49;
t65 = pkin(3) * t95;
t86 = t37 * pkin(4) - t38 * qJ(5) + t65;
t15 = t74 * t41 + t33;
t85 = t15 * qJD(3) - t4;
t59 = t76 * qJD(2) * pkin(6) + t104;
t12 = t108 * t36 - t116;
t16 = -t108 * t84 + t74 * t42;
t18 = t108 * t42 + t74 * t84;
t31 = -t108 * t61 + t74 * t96;
t82 = t16 * t83 - t18 * t47 + t30 * t38 - t31 * t37 + t4 * t56;
t81 = t126 * qJD(3);
t53 = -pkin(6) * t95 + t69;
t54 = t59 * qJD(3);
t58 = t71 - t97;
t80 = t53 * t76 + t54 * t75 + (-t58 * t76 - t59 * t75) * qJD(3);
t68 = -t108 * pkin(3) - pkin(4);
t66 = t74 * pkin(3) + qJ(5);
t40 = t52 * qJD(3);
t21 = t55 * pkin(4) - t56 * qJ(5) + t70;
t20 = -t62 + (t89 + t47) * qJD(3);
t19 = -t62 + (t89 - t47) * qJD(3);
t14 = pkin(4) * t83 + t47 * qJ(5) + t98;
t9 = qJD(3) * qJ(5) + t13;
t8 = -qJD(3) * pkin(4) + qJD(5) - t12;
t7 = t49 * pkin(4) - t52 * qJ(5) - t56 * qJD(5) + t99;
t3 = qJD(3) * qJD(5) + t5;
t2 = t38 * t56 + t52 * t83;
t1 = -qJD(5) * t83 + t86;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, -t113, 0, t53 * t75 - t54 * t76 + (-t58 * t75 + t59 * t76) * qJD(3), 0, 0, 0, 0, 0, 0, -t39, -t40, t129, -t12 * t49 + t13 * t52 + t5 * t56 + t120, 0, 0, 0, 0, 0, 0, -t39, t129, t40, t3 * t56 + t8 * t49 + t9 * t52 + t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t90, t109 * t127, t113, -0.2e1 * t90, -t114, 0, -pkin(6) * t113 + t75 * t91, pkin(6) * t114 + t76 * t91, t80, t80 * pkin(6), t2, -t128, t40, t88, -t39, 0, t70 * t37 + t60 * t49 + (-t16 + (qJD(2) * t55 + t47) * t122) * qJD(3), t70 * t38 + t60 * t52 + (t126 * t122 - t18) * qJD(3), -t12 * t52 - t13 * t49 - t5 * t55 + t82, -t12 * t16 + t13 * t18 + t121 + t5 * t31 + (t60 + t107) * t99, t2, t40, t128, 0, t39, t88, -t16 * qJD(3) + t1 * t55 + t11 * t49 + t21 * t37 + t7 * t47, -t3 * t55 - t9 * t49 + t8 * t52 + t82, t18 * qJD(3) - t1 * t56 - t11 * t52 - t21 * t38 - t7 * t83, t1 * t21 + t11 * t7 + t8 * t16 + t9 * t18 + t3 * t31 + t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, t109 * t78, 0, t100, 0, 0, t78 * pkin(2) * t75, pkin(2) * t115 - t69 + (t58 + t97) * qJD(3), 0, 0, t118, t124, t20, -t118, 0, 0, -t47 * t98 - t60 * t83 + t85, t17 * qJD(3) + t60 * t47 - t83 * t98 - t5, (t13 - t15) * t83 + (-t12 + t17) * t47 + (-t108 * t38 - t37 * t74) * pkin(3), t12 * t15 - t13 * t17 + (-t60 * t106 - t108 * t4 + t5 * t74) * pkin(3), t118, t20, -t124, 0, 0, -t118, -t14 * t47 - t119 + t85, -t66 * t37 + t68 * t38 + (-t15 + t9) * t83 + (t8 - t103) * t47, -t11 * t47 + t14 * t83 + (0.2e1 * qJD(5) - t17) * qJD(3) + t5, t103 * t9 - t11 * t14 - t8 * t15 + t3 * t66 + t4 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t19, t125, t12 * t83 + t13 * t47 + t65, 0, 0, 0, 0, 0, 0, t81, t125, -t19, t9 * t47 + (-qJD(5) - t8) * t83 + t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t20, -t123 - t77, -t9 * qJD(3) + t119 + t4;];
tauc_reg = t6;
