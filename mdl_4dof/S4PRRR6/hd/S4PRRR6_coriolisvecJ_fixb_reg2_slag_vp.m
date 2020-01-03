% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:06
% EndTime: 2019-12-31 16:35:09
% DurationCPUTime: 0.65s
% Computational Cost: add. (747->132), mult. (1995->203), div. (0->0), fcn. (1337->6), ass. (0->97)
t54 = sin(qJ(4));
t55 = sin(qJ(3));
t57 = cos(qJ(4));
t58 = cos(qJ(3));
t37 = t54 * t58 + t57 * t55;
t51 = qJD(3) + qJD(4);
t113 = t51 * t37;
t15 = t113 * qJD(2);
t85 = qJD(2) * qJD(3);
t112 = -0.2e1 * t85;
t101 = t57 * t58;
t103 = t54 * t55;
t36 = -t101 + t103;
t56 = sin(qJ(2));
t30 = t36 * t56;
t60 = qJD(3) ^ 2;
t61 = qJD(2) ^ 2;
t111 = (t60 + t61) * t56;
t87 = t56 * qJD(1);
t92 = qJD(2) * pkin(5);
t44 = t87 + t92;
t59 = cos(qJ(2));
t86 = t59 * qJD(1);
t89 = qJD(3) * t55;
t19 = -t44 * t89 + (-pkin(6) * t89 + t58 * t86) * qJD(2);
t74 = pkin(6) * qJD(2) + t44;
t27 = t74 * t55;
t24 = qJD(3) * pkin(3) - t27;
t110 = (qJD(4) * t24 + t19) * t57;
t109 = -pkin(6) - pkin(5);
t40 = t109 * t55;
t41 = t109 * t58;
t22 = t54 * t40 - t57 * t41;
t78 = qJD(3) * t109;
t38 = t55 * t78;
t39 = t58 * t78;
t65 = t37 * t59;
t108 = -qJD(1) * t65 + qJD(4) * t22 + t54 * t38 - t57 * t39;
t21 = t57 * t40 + t54 * t41;
t64 = t36 * t59;
t107 = -qJD(1) * t64 - qJD(4) * t21 - t57 * t38 - t54 * t39;
t80 = qJD(2) * t101;
t91 = qJD(2) * t55;
t81 = t54 * t91;
t31 = -t80 + t81;
t33 = t37 * qJD(2);
t106 = t33 * t31;
t50 = -t58 * pkin(3) - pkin(2);
t34 = qJD(2) * t50 - t86;
t105 = t34 * t33;
t28 = t74 * t58;
t104 = t54 * t28;
t102 = t57 * t28;
t100 = t60 * t55;
t99 = t60 * t58;
t98 = t61 * t59;
t77 = t58 * t85;
t97 = -qJD(4) * t80 - t57 * t77;
t52 = t55 ^ 2;
t53 = t58 ^ 2;
t96 = t52 - t53;
t95 = t52 + t53;
t93 = qJD(2) * pkin(2);
t90 = qJD(2) * t56;
t88 = qJD(3) * t58;
t84 = t55 * t61 * t58;
t83 = pkin(3) * t91;
t82 = pkin(3) * t89;
t79 = -pkin(3) * t51 - t24;
t20 = -t44 * t88 + (-pkin(6) * t88 - t55 * t86) * qJD(2);
t76 = -t54 * t19 + t57 * t20;
t75 = -qJD(4) * t104 + t54 * t20;
t45 = -t86 - t93;
t73 = -t45 - t86;
t71 = t59 * t112;
t70 = t55 * t77;
t69 = t51 * t103;
t11 = t54 * t24 + t102;
t68 = qJD(2) * t73;
t67 = t34 * t31 - t75;
t66 = t82 - t87;
t63 = qJD(3) * (-t73 - t93);
t2 = -qJD(4) * t11 + t76;
t35 = (t82 + t87) * qJD(2);
t29 = t37 * t56;
t17 = -qJD(4) * t101 - t57 * t88 + t69;
t14 = qJD(2) * t69 + t97;
t13 = -t57 * t27 - t104;
t12 = t54 * t27 - t102;
t10 = t57 * t24 - t104;
t9 = -t31 ^ 2 + t33 ^ 2;
t8 = t33 * t51 - t15;
t7 = -t97 + (t31 - t81) * t51;
t4 = -qJD(2) * t65 + t51 * t30;
t3 = -qJD(2) * t64 - t113 * t56;
t1 = t75 + t110;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61 * t56, -t98, 0, 0, 0, 0, 0, 0, 0, 0, -t58 * t111 + t55 * t71, t55 * t111 + t58 * t71, t95 * t98, (t45 * t56 + (-t87 + (t44 + t87) * t95) * t59) * qJD(2), 0, 0, 0, 0, 0, 0, -t59 * t15 + t31 * t90 + t4 * t51, t59 * t14 - t3 * t51 + t33 * t90, -t29 * t14 + t30 * t15 - t3 * t31 - t4 * t33, -t1 * t30 + t10 * t4 + t11 * t3 - t2 * t29 + t34 * t90 - t35 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t70, t96 * t112, t99, -0.2e1 * t70, -t100, 0, -pkin(5) * t99 + t55 * t63, pkin(5) * t100 + t58 * t63, 0, ((-t45 - t93) * t56 + (-t44 + t92) * t59 * t95) * qJD(1), -t14 * t37 - t33 * t17, -t113 * t33 + t14 * t36 - t37 * t15 + t17 * t31, -t17 * t51, t113 * t31 + t15 * t36, -t113 * t51, 0, -t108 * t51 + t113 * t34 + t50 * t15 + t66 * t31 + t35 * t36, t107 * t51 - t50 * t14 - t34 * t17 + t66 * t33 + t35 * t37, -t1 * t36 + t10 * t17 + t107 * t31 + t108 * t33 - t11 * t113 + t21 * t14 - t22 * t15 - t2 * t37, t1 * t22 - t108 * t10 - t107 * t11 + t2 * t21 + t66 * t34 + t35 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, t96 * t61, 0, t84, 0, 0, t55 * t68, t58 * t68, 0, 0, t106, t9, t7, -t106, t8, 0, -t31 * t83 - t12 * t51 - t105 + (t54 * t79 - t102) * qJD(4) + t76, -t33 * t83 + t13 * t51 + (qJD(4) * t79 - t19) * t57 + t67, (t11 + t12) * t33 + (-t10 + t13) * t31 + (t14 * t57 - t15 * t54 + (-t31 * t57 + t33 * t54) * qJD(4)) * pkin(3), -t10 * t12 - t11 * t13 + (-t34 * t91 + t1 * t54 + t2 * t57 + (-t10 * t54 + t11 * t57) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t9, t7, -t106, t8, 0, t11 * t51 - t105 + t2, t10 * t51 - t110 + t67, 0, 0;];
tauc_reg = t5;
