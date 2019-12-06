% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:35
% EndTime: 2019-12-05 16:17:38
% DurationCPUTime: 0.58s
% Computational Cost: add. (974->119), mult. (1723->194), div. (0->0), fcn. (968->6), ass. (0->97)
t55 = cos(qJ(3));
t89 = pkin(2) * qJD(2);
t62 = -t55 * t89 + qJD(4);
t47 = qJD(2) + qJD(3);
t88 = pkin(2) * qJD(3);
t72 = qJD(2) * t88;
t34 = t47 * qJD(4) + t55 * t72;
t50 = sin(pkin(9));
t51 = cos(pkin(9));
t37 = -t51 * pkin(4) - t50 * pkin(7) - pkin(3);
t15 = t37 * t47 + t62;
t53 = sin(qJ(3));
t36 = t47 * qJ(4) + t53 * t89;
t22 = t50 * qJD(1) + t51 * t36;
t52 = sin(qJ(5));
t54 = cos(qJ(5));
t5 = t54 * t15 - t52 * t22;
t68 = t53 * t72;
t95 = t51 * t54;
t2 = qJD(5) * t5 + t34 * t95 + t52 * t68;
t109 = t2 * t52;
t108 = t53 * pkin(2);
t107 = t55 * pkin(2);
t45 = t50 ^ 2;
t99 = t45 * t54;
t106 = t2 * t51 + t34 * t99;
t86 = qJ(4) * t51;
t19 = t54 * t37 - t52 * t86;
t83 = qJD(5) * t19;
t84 = qJD(4) * t51;
t94 = t51 * t55;
t105 = t54 * t84 - (t52 * t53 + t54 * t94) * t89 + t83;
t21 = -t51 * qJD(1) + t50 * t36;
t104 = t21 * t50;
t103 = t34 * t51;
t41 = t55 * t88 + qJD(4);
t102 = t41 * t47;
t28 = t45 * t34;
t44 = t47 ^ 2;
t101 = t45 * t44;
t100 = t45 * t47;
t46 = t51 ^ 2;
t29 = t46 * t34;
t98 = t47 * t50;
t97 = t51 * t47;
t96 = t51 * t52;
t20 = t52 * t37 + t54 * t86;
t82 = qJD(5) * t20;
t93 = -t52 * t84 - (-t52 * t94 + t53 * t54) * t89 - t82;
t92 = t28 + t29;
t91 = t45 + t46;
t90 = t52 ^ 2 - t54 ^ 2;
t87 = qJ(4) * t34;
t6 = t52 * t15 + t54 * t22;
t85 = qJD(5) * t6;
t81 = qJD(5) * t52;
t80 = qJD(5) * t54;
t79 = qJ(4) * qJD(5);
t78 = t21 * t98;
t77 = t53 * t88;
t75 = qJD(5) * t100;
t74 = t50 * t81;
t73 = t50 * t80;
t71 = t98 * t108;
t70 = t44 * t52 * t99;
t40 = -qJD(5) + t97;
t69 = t40 * t74;
t3 = -t34 * t96 + t54 * t68 - t85;
t67 = t21 * t73 + t52 * t28 - t3 * t51;
t66 = (-qJD(5) - t40) * t98;
t65 = t5 * t52 - t54 * t6;
t64 = (-qJD(3) + t47) * t89;
t63 = (-qJD(2) - t47) * t88;
t61 = t54 * t52 * t75;
t60 = t22 * t51 + t104;
t59 = t53 * t63;
t58 = t53 * t64;
t30 = t37 - t107;
t42 = qJ(4) + t108;
t14 = t52 * t30 + t42 * t95;
t13 = t54 * t30 - t42 * t96;
t57 = -t109 + (-t3 - t85) * t54;
t56 = -t40 ^ 2 - t101;
t38 = t50 * t68;
t35 = -t47 * pkin(3) + t62;
t33 = t74 * t97;
t32 = -0.2e1 * t61;
t31 = 0.2e1 * t61;
t27 = t45 * t87;
t18 = t42 * t28;
t17 = 0.2e1 * t90 * t75;
t12 = (t40 + t97) * t73;
t11 = t33 + t69;
t8 = -qJD(5) * t14 - t41 * t96 + t54 * t77;
t7 = qJD(5) * t13 + t41 * t95 + t52 * t77;
t4 = t5 * t74;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t40 - t97) * t73, t33 - t69, 0, (t2 * t54 - t3 * t52 - t103 + (-t5 * t54 - t52 * t6) * qJD(5)) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t55 * t63, 0, 0, 0, 0, 0, 0, 0, 0, t51 * t59, qJD(3) * t71 + t38, t91 * t102 + t92, t42 * t29 + t18 + t60 * t41 + (t35 + (-pkin(3) - t107) * qJD(2)) * t77, t32, t17, t11, t31, t12, 0, -t8 * t40 + (t41 * t52 + t42 * t80) * t100 + t67, t99 * t102 + t7 * t40 + (-t42 * t100 - t104) * t81 + t106, t4 + ((-t52 * t7 - t54 * t8 + (t13 * t52 - t14 * t54) * qJD(5)) * t47 + t57) * t50, t41 * t104 + t3 * t13 + t2 * t14 + t5 * t8 + t6 * t7 + t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t55 * t64, 0, 0, 0, 0, 0, 0, 0, 0, t51 * t58, -qJD(2) * t71 + t38, t62 * t47 * t91 + t92, t46 * t87 + t27 + t60 * qJD(4) + (-t60 * t55 + (-pkin(3) * qJD(3) - t35) * t53) * t89, t32, t17, t11, t31, t12, 0, -t93 * t40 + (t52 * t62 + t54 * t79) * t100 + t67, -t21 * t74 + t105 * t40 + (-t52 * t79 + t54 * t62) * t100 + t106, t4 + (((-t82 - t93) * t54 + (t83 - t105) * t52) * t47 + t57) * t50, t62 * t104 + t105 * t6 + t3 * t19 + t2 * t20 + t93 * t5 + t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91 * t44, -t47 * t60 + t68, 0, 0, 0, 0, 0, 0, t56 * t52, t56 * t54, 0, t109 + t3 * t54 - t65 * qJD(5) + (t65 * t51 - t104) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t90 * t101, t52 * t66, -t70, t54 * t66, 0, -t6 * t40 - t54 * t78 + t3, -t5 * t40 + (-qJD(5) * t15 - t103) * t54 + (qJD(5) * t22 - t68 + t78) * t52, 0, 0;];
tauc_reg = t1;
