% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4PRRR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:48
% EndTime: 2019-12-31 16:36:52
% DurationCPUTime: 1.11s
% Computational Cost: add. (1009->192), mult. (2819->310), div. (0->0), fcn. (2006->8), ass. (0->111)
t51 = sin(qJ(2));
t47 = sin(pkin(4));
t97 = qJD(1) * t47;
t76 = t51 * t97;
t37 = qJD(2) * pkin(6) + t76;
t53 = cos(qJ(3));
t48 = cos(pkin(4));
t54 = cos(qJ(2));
t95 = qJD(2) * t47;
t80 = t54 * t95;
t59 = qJD(1) * (qJD(3) * t48 + t80);
t50 = sin(qJ(3));
t93 = qJD(3) * t50;
t11 = -t37 * t93 + t53 * t59;
t68 = pkin(3) * t50 - pkin(7) * t53;
t36 = t68 * qJD(3);
t24 = (t36 + t76) * qJD(2);
t49 = sin(qJ(4));
t52 = cos(qJ(4));
t96 = qJD(1) * t48;
t26 = t53 * t37 + t50 * t96;
t19 = qJD(3) * pkin(7) + t26;
t39 = -t53 * pkin(3) - t50 * pkin(7) - pkin(2);
t75 = t54 * t97;
t27 = t39 * qJD(2) - t75;
t64 = t49 * t19 - t52 * t27;
t1 = -t64 * qJD(4) + t52 * t11 + t49 * t24;
t86 = t53 * qJD(2);
t41 = -qJD(4) + t86;
t126 = -t64 * t41 + t1;
t6 = t52 * t19 + t49 * t27;
t2 = -qJD(4) * t6 - t49 * t11 + t52 * t24;
t125 = t6 * t41 - t2;
t88 = t49 * qJD(3);
t90 = qJD(4) * t52;
t60 = t50 * t90 + t53 * t88;
t84 = qJD(3) * qJD(4);
t21 = t60 * qJD(2) + t49 * t84;
t102 = t53 * t54;
t89 = qJD(4) * t53;
t91 = qJD(4) * t49;
t122 = (-t49 * t102 + t51 * t52) * t97 + t39 * t91 - t52 * t36 - (t50 * t88 - t52 * t89) * pkin(6);
t87 = t52 * qJD(3);
t121 = (t52 * t102 + t49 * t51) * t97 - t39 * t90 - t49 * t36 - (-t49 * t89 - t50 * t87) * pkin(6);
t92 = qJD(3) * t53;
t12 = t37 * t92 + t50 * t59;
t107 = t47 * t51;
t30 = t50 * t107 - t48 * t53;
t120 = t12 * t30;
t119 = t12 * t49;
t118 = t12 * t50;
t117 = t12 * t52;
t25 = -t50 * t37 + t53 * t96;
t18 = -qJD(3) * pkin(3) - t25;
t116 = t18 * t49;
t115 = t18 * t52;
t94 = qJD(2) * t50;
t74 = t49 * t94;
t79 = t53 * t87;
t20 = -qJD(2) * t79 + qJD(4) * t74 - t52 * t84;
t114 = t20 * t49;
t113 = t21 * t52;
t32 = t74 - t87;
t112 = t32 * t41;
t34 = t52 * t94 + t88;
t111 = t34 * t32;
t110 = t34 * t41;
t109 = t41 * t49;
t108 = t41 * t52;
t106 = t47 * t54;
t56 = qJD(2) ^ 2;
t105 = t47 * t56;
t104 = t49 * t53;
t103 = t52 * t53;
t55 = qJD(3) ^ 2;
t101 = t55 * t50;
t100 = t55 * t53;
t45 = t50 ^ 2;
t46 = t53 ^ 2;
t99 = t45 - t46;
t98 = qJD(2) * pkin(2);
t85 = qJD(2) * qJD(3);
t83 = t51 * t105;
t82 = t50 * t56 * t53;
t81 = t51 * t95;
t78 = t50 * t91;
t77 = t41 * t90;
t72 = t50 * t85;
t71 = t50 * t80;
t70 = t53 * t80;
t69 = t53 * t72;
t38 = -t75 - t98;
t67 = -t38 - t75;
t66 = -t49 * t6 + t52 * t64;
t63 = qJD(2) * t45 - t41 * t53;
t31 = t53 * t107 + t48 * t50;
t15 = -t52 * t106 - t31 * t49;
t62 = t49 * t106 - t31 * t52;
t61 = qJD(2) * t67;
t58 = qJD(3) * (-t67 - t98);
t57 = t11 * t53 + t118 + (-t25 * t53 - t26 * t50) * qJD(3);
t35 = t68 * qJD(2);
t29 = pkin(6) * t103 + t49 * t39;
t28 = -pkin(6) * t104 + t52 * t39;
t14 = t31 * qJD(3) + t71;
t13 = -t30 * qJD(3) + t70;
t10 = t52 * t25 + t49 * t35;
t9 = -t49 * t25 + t52 * t35;
t4 = t15 * qJD(4) + t13 * t52 + t49 * t81;
t3 = t62 * qJD(4) - t13 * t49 + t52 * t81;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -t54 * t105, 0, 0, 0, 0, 0, 0, 0, 0, -t53 * t83 + (-t14 - t71) * qJD(3), t50 * t83 + (-t13 - t70) * qJD(3), (t13 * t53 + t14 * t50 + (t30 * t53 - t31 * t50) * qJD(3)) * qJD(2), t11 * t31 + t120 + t26 * t13 - t25 * t14 + (t38 - t75) * t81, 0, 0, 0, 0, 0, 0, t14 * t32 + t15 * t72 + t30 * t21 - t3 * t41, t14 * t34 - t30 * t20 + t4 * t41 + t62 * t72, t15 * t20 + t21 * t62 - t3 * t34 - t4 * t32, -t1 * t62 + t18 * t14 + t2 * t15 - t3 * t64 + t6 * t4 + t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t69, -0.2e1 * t99 * t85, t100, -0.2e1 * t69, -t101, 0, -pkin(6) * t100 + t50 * t58, pkin(6) * t101 + t53 * t58, (-t45 - t46) * qJD(2) * t75 + t57, ((t25 * t50 - t26 * t53) * t54 + (-t38 - t98) * t51) * t97 + t57 * pkin(6), -t20 * t52 * t50 + (-t78 + t79) * t34, (-t32 * t52 - t34 * t49) * t92 + (t114 - t113 + (t32 * t49 - t34 * t52) * qJD(4)) * t50, t41 * t78 + t20 * t53 + (t34 * t50 + t52 * t63) * qJD(3), t21 * t49 * t50 + t32 * t60, t50 * t77 + t21 * t53 + (-t32 * t50 - t49 * t63) * qJD(3), (-t41 - t86) * t93, t122 * t41 + (-t2 + (pkin(6) * t32 + t116) * qJD(3)) * t53 + (-t32 * t75 + t18 * t90 + pkin(6) * t21 + t119 + (qJD(2) * t28 - t64) * qJD(3)) * t50, -t121 * t41 + (t1 + (pkin(6) * t34 + t115) * qJD(3)) * t53 + (-t34 * t75 - t18 * t91 - pkin(6) * t20 + t117 + (-qJD(2) * t29 - t6) * qJD(3)) * t50, t28 * t20 - t29 * t21 + t122 * t34 + t121 * t32 + t66 * t92 + (-t1 * t49 - t2 * t52 + (-t49 * t64 - t52 * t6) * qJD(4)) * t50, -t18 * t50 * t75 + t1 * t29 + t2 * t28 - t121 * t6 + t122 * t64 + (t18 * t92 + t118) * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82, t99 * t56, 0, t82, 0, 0, t50 * t61, t53 * t61, 0, 0, -t108 * t34 - t114, (-t20 + t112) * t52 + (-t21 + t110) * t49, -t77 + (t41 * t103 + (-t34 + t88) * t50) * qJD(2), -t109 * t32 - t113, t41 * t91 + (-t41 * t104 + (t32 + t87) * t50) * qJD(2), t41 * t94, -pkin(3) * t21 - t117 - t26 * t32 + t9 * t41 + (pkin(7) * t108 + t116) * qJD(4) + (t64 * t50 + (-pkin(7) * t93 - t18 * t53) * t49) * qJD(2), pkin(3) * t20 - t10 * t41 + t119 - t26 * t34 + (-pkin(7) * t109 + t115) * qJD(4) + (-t18 * t103 + (-pkin(7) * t87 + t6) * t50) * qJD(2), t10 * t32 + t9 * t34 + ((qJD(4) * t34 - t21) * pkin(7) + t126) * t52 + ((qJD(4) * t32 - t20) * pkin(7) + t125) * t49, -t12 * pkin(3) - t6 * t10 - t18 * t26 + t64 * t9 + (qJD(4) * t66 + t1 * t52 - t2 * t49) * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, -t32 ^ 2 + t34 ^ 2, -t20 - t112, -t111, -t110 - t21, t72, -t18 * t34 - t125, t18 * t32 - t126, 0, 0;];
tauc_reg = t5;
