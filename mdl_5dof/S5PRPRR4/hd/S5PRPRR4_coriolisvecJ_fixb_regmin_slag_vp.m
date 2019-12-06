% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [5x19]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:35
% EndTime: 2019-12-05 15:51:38
% DurationCPUTime: 0.75s
% Computational Cost: add. (657->161), mult. (1796->259), div. (0->0), fcn. (1378->10), ass. (0->105)
t50 = sin(pkin(10));
t51 = sin(pkin(5));
t52 = cos(pkin(10));
t56 = sin(qJ(2));
t59 = cos(qJ(2));
t28 = (t50 * t56 - t52 * t59) * t51;
t58 = cos(qJ(4));
t89 = t58 * qJD(2);
t43 = -qJD(5) + t89;
t117 = -qJD(5) - t43;
t55 = sin(qJ(4));
t48 = t55 ^ 2;
t68 = qJD(2) * t48 - t43 * t58;
t54 = sin(qJ(5));
t93 = qJD(5) * t54;
t85 = t55 * t93;
t57 = cos(qJ(5));
t90 = t57 * qJD(4);
t116 = -t43 * t85 - t68 * t90;
t92 = qJD(5) * t57;
t84 = t55 * t92;
t87 = qJD(4) * qJD(5);
t91 = t54 * qJD(4);
t23 = (t58 * t91 + t84) * qJD(2) + t54 * t87;
t98 = qJD(1) * t51;
t83 = t59 * t98;
t38 = qJD(2) * pkin(2) + t83;
t86 = t56 * t98;
t40 = t52 * t86;
t19 = t50 * t38 + t40;
t17 = qJD(2) * pkin(7) + t19;
t53 = cos(pkin(5));
t42 = t53 * qJD(1) + qJD(3);
t10 = t58 * t17 + t55 * t42;
t26 = qJD(2) * t28;
t21 = qJD(1) * t26;
t4 = t10 * qJD(4) - t55 * t21;
t115 = t4 * t54;
t114 = t4 * t57;
t113 = t52 * pkin(2);
t81 = t58 * t90;
t97 = qJD(2) * t55;
t82 = t54 * t97;
t22 = qJD(2) * t81 - qJD(5) * t82 + t57 * t87;
t112 = t22 * t54;
t33 = t82 - t90;
t111 = t33 * t43;
t35 = t57 * t97 + t91;
t110 = t35 * t43;
t109 = t43 * t54;
t61 = qJD(2) ^ 2;
t108 = t51 * t61;
t107 = t54 * t58;
t106 = t55 * t33;
t105 = t57 * t43;
t104 = t57 * t58;
t103 = t58 * t23;
t60 = qJD(4) ^ 2;
t102 = t60 * t55;
t101 = t60 * t58;
t24 = t50 * t83 + t40;
t73 = pkin(4) * t55 - pkin(8) * t58;
t37 = t73 * qJD(4);
t100 = t24 - t37;
t99 = -t58 ^ 2 + t48;
t44 = t50 * pkin(2) + pkin(7);
t96 = qJD(4) * t44;
t95 = qJD(4) * t55;
t94 = qJD(4) * t58;
t88 = qJD(2) * qJD(4);
t8 = qJD(4) * pkin(8) + t10;
t80 = t43 * t44 + t8;
t78 = t55 * t88;
t77 = -t22 * t58 + t35 * t95;
t39 = t50 * t86;
t18 = t52 * t38 - t39;
t16 = -qJD(2) * pkin(3) - t18;
t76 = -qJD(2) * t16 + t21;
t74 = t43 * t84;
t67 = -t58 * pkin(4) - t55 * pkin(8) - pkin(3);
t12 = t67 * qJD(2) - t18;
t1 = t57 * t12 - t54 * t8;
t2 = t54 * t12 + t57 * t8;
t29 = (t50 * t59 + t52 * t56) * t51;
t15 = t29 * t58 + t53 * t55;
t72 = t15 * t57 + t28 * t54;
t71 = -t15 * t54 + t28 * t57;
t70 = t55 * t17 - t58 * t42;
t14 = t29 * t55 - t53 * t58;
t25 = qJD(2) * t29;
t20 = qJD(1) * t25;
t66 = qJD(2) * t24 - t44 * t60 - t20;
t27 = t52 * t83 - t39;
t65 = qJD(4) * (qJD(2) * (-pkin(3) - t113) + t16 + t27);
t64 = t68 * t54;
t3 = -t70 * qJD(4) - t58 * t21;
t7 = -qJD(4) * pkin(4) + t70;
t62 = qJD(4) * t7 + qJD(5) * t12 - t27 * t43 + t3;
t36 = t73 * qJD(2);
t32 = t67 - t113;
t13 = (qJD(1) * t29 + t37) * qJD(2);
t11 = t57 * t13;
t6 = -t14 * qJD(4) - t26 * t58;
t5 = t15 * qJD(4) - t26 * t55;
t9 = [0, 0, -t56 * t108, -t59 * t108, -t18 * t25 - t19 * t26 + t20 * t28 - t21 * t29, 0, 0, 0, 0, 0, -t5 * qJD(4) + (-t25 * t58 + t28 * t95) * qJD(2), -t6 * qJD(4) + (t25 * t55 + t28 * t94) * qJD(2), 0, 0, 0, 0, 0, -(-t72 * qJD(5) + t25 * t57 - t6 * t54) * t43 + t71 * t78 + t5 * t33 + t14 * t23, (t71 * qJD(5) + t25 * t54 + t6 * t57) * t43 - t72 * t78 + t5 * t35 + t14 * t22; 0, 0, 0, 0, t18 * t24 - t19 * t27 + (-t20 * t52 - t21 * t50) * pkin(2), 0.2e1 * t58 * t78, -0.2e1 * t99 * t88, t101, -t102, 0, t55 * t65 + t66 * t58, -t66 * t55 + t58 * t65, t22 * t57 * t55 + (t81 - t85) * t35, (-t33 * t57 - t35 * t54) * t94 + (-t112 - t23 * t57 + (t33 * t54 - t35 * t57) * qJD(5)) * t55, -t116 + t77, t74 + t103 + (-t64 - t106) * qJD(4), (-t43 - t89) * t95, (t100 * t57 + t32 * t93) * t43 + (t33 * t96 + t62 * t54 + t80 * t92 - t11) * t58 + (t7 * t92 + t44 * t23 - t27 * t33 + t115 + (-t44 * t109 + (-t44 * t107 + t57 * t32) * qJD(2) + t1) * qJD(4)) * t55, (-t100 * t54 + t32 * t92) * t43 + (t35 * t96 + (-t80 * qJD(5) + t13) * t54 + t62 * t57) * t58 + (-t7 * t93 + t44 * t22 - t27 * t35 + t114 + (-t44 * t105 - (t44 * t104 + t54 * t32) * qJD(2) - t2) * qJD(4)) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t101, 0, 0, 0, 0, 0, t74 - t103 + (-t64 + t106) * qJD(4), t116 + t77; 0, 0, 0, 0, 0, -t55 * t61 * t58, t99 * t61, 0, 0, 0, t76 * t55, t76 * t58, -t35 * t105 + t112, (t22 + t111) * t57 + (-t23 + t110) * t54, -t43 * t92 + (t43 * t104 + (-t35 + t91) * t55) * qJD(2), t43 * t93 + (-t43 * t107 + (t33 + t90) * t55) * qJD(2), t43 * t97, -pkin(4) * t23 - t114 + (t57 * t36 + t54 * t70) * t43 - t10 * t33 + (pkin(8) * t105 + t7 * t54) * qJD(5) + (-t1 * t55 + (-pkin(8) * t95 - t58 * t7) * t54) * qJD(2), -pkin(4) * t22 + t115 - (t54 * t36 - t57 * t70) * t43 - t10 * t35 + (-pkin(8) * t109 + t7 * t57) * qJD(5) + (-t7 * t104 + (-pkin(8) * t90 + t2) * t55) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t33, -t33 ^ 2 + t35 ^ 2, t22 - t111, -t110 - t23, t78, t117 * t2 - t54 * t3 - t7 * t35 + t11, t117 * t1 - t54 * t13 - t57 * t3 + t7 * t33;];
tauc_reg = t9;
