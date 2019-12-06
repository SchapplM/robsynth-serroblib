% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRP2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:05
% EndTime: 2019-12-05 15:31:07
% DurationCPUTime: 0.64s
% Computational Cost: add. (694->134), mult. (1784->208), div. (0->0), fcn. (1087->4), ass. (0->104)
t61 = sin(pkin(8));
t62 = cos(pkin(8));
t45 = -t62 * pkin(3) - t61 * pkin(6) - pkin(2);
t64 = cos(qJ(4));
t53 = t64 * t62 * qJ(3);
t63 = sin(qJ(4));
t25 = t63 * t45 + t53;
t31 = qJD(2) * t45 + qJD(3);
t96 = qJ(3) * qJD(2);
t43 = t61 * qJD(1) + t62 * t96;
t11 = t63 * t31 + t64 * t43;
t68 = qJD(4) * t11;
t94 = qJD(2) * qJD(3);
t80 = t62 * t94;
t7 = -t63 * t80 - t68;
t10 = t64 * t31 - t63 * t43;
t71 = t10 * t63 - t11 * t64;
t100 = qJD(4) * t64;
t101 = qJD(4) * t63;
t91 = -t31 * t100 + t43 * t101 - t64 * t80;
t116 = -t71 * qJD(4) - t63 * t91 + t7 * t64;
t104 = qJD(2) * t61;
t81 = qJ(5) * t100;
t99 = qJD(5) * t63;
t1 = (-t81 - t99) * t104 - t91;
t86 = t61 * t101;
t77 = qJD(2) * t86;
t46 = qJ(5) * t77;
t111 = t61 * t64;
t103 = qJD(3) * t63;
t84 = t62 * t103;
t67 = -qJD(5) * t111 - t84;
t2 = qJD(2) * t67 + t46 - t68;
t98 = t62 * qJD(2);
t54 = -qJD(4) + t98;
t82 = qJ(5) * t104;
t8 = -t64 * t82 + t10;
t3 = -t54 * pkin(4) + t8;
t9 = -t63 * t82 + t11;
t73 = t3 * t63 - t64 * t9;
t115 = -t73 * qJD(4) + t1 * t63 + t2 * t64;
t114 = t3 - t8;
t56 = t62 * qJD(1);
t42 = t61 * t96 - t56;
t113 = t42 * t61;
t57 = t61 ^ 2;
t65 = qJD(2) ^ 2;
t112 = t57 * t65;
t102 = qJD(3) * t64;
t110 = t45 * t100 + t62 * t102;
t93 = qJD(2) * qJD(4);
t79 = t64 * t93;
t76 = t61 * t79;
t30 = pkin(4) * t76 + t61 * t94;
t58 = t62 ^ 2;
t109 = t57 + t58;
t59 = t63 ^ 2;
t60 = t64 ^ 2;
t108 = t59 - t60;
t107 = qJ(3) * t63;
t106 = qJ(5) * t61;
t105 = qJD(2) * t57;
t39 = (pkin(4) * t63 + qJ(3)) * t61;
t26 = qJD(2) * t39 + qJD(5) - t56;
t97 = qJD(5) + t26;
t95 = qJ(3) * qJD(4);
t92 = t63 * t112;
t90 = t62 * t107;
t89 = t64 * t106;
t88 = t63 * t104;
t87 = t64 * t104;
t85 = t61 * t100;
t83 = t63 * t95;
t78 = t54 * t86;
t74 = t3 * t64 + t63 * t9;
t70 = -t43 * t62 - t113;
t69 = t57 * t63 * t79;
t66 = -t54 ^ 2 - t112;
t51 = t57 * qJ(3) * t94;
t47 = t64 * t92;
t44 = t62 * t77;
t41 = -0.2e1 * t69;
t40 = 0.2e1 * t69;
t38 = t64 * t45;
t36 = (pkin(4) * t100 + qJD(3)) * t61;
t35 = t54 * t87;
t33 = t108 * t112;
t28 = 0.2e1 * t108 * t57 * t93;
t24 = t38 - t90;
t23 = -t35 - t76;
t22 = (-qJD(4) - t54) * t88;
t21 = (t54 - t98) * t85;
t20 = (t54 + t98) * t85;
t19 = t44 + t78;
t18 = t44 - t78;
t17 = -t63 * t106 + t25;
t16 = -t25 * qJD(4) - t84;
t15 = -t62 * t83 + t110;
t14 = -t89 + t38 + (-pkin(4) - t107) * t62;
t13 = t66 * t64;
t12 = t66 * t63;
t5 = (-t53 + (-t45 + t106) * t63) * qJD(4) + t67;
t4 = -t61 * t99 + (-t89 - t90) * qJD(4) + t110;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t18, 0, (-t80 - t91 * t64 - t63 * t7 + (-t10 * t64 - t11 * t63) * qJD(4)) * t61, 0, 0, 0, 0, 0, 0, t21, t18, 0, -t30 * t62 + (-t74 * qJD(4) + t1 * t64 - t2 * t63) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t109 * t94, t51 + (t58 * t96 - t70) * qJD(3), t41, t28, t19, t40, t20, 0, t42 * t85 - t16 * t54 - t7 * t62 + (t64 * t95 + 0.2e1 * t103) * t105, -t42 * t86 + t15 * t54 - t91 * t62 + (-t83 + 0.2e1 * t102) * t105, ((-t15 * t63 - t16 * t64 + (t24 * t63 - t25 * t64) * qJD(4)) * qJD(2) - t116) * t61, qJD(3) * t113 + t10 * t16 + t11 * t15 + t7 * t24 - t25 * t91 + t51, t41, t28, t19, t40, t20, 0, -t2 * t62 - t5 * t54 + (t26 * t100 + t30 * t63 + (t39 * t100 + t36 * t63) * qJD(2)) * t61, t1 * t62 + t4 * t54 + (-t26 * t101 + t30 * t64 + (-t39 * t101 + t36 * t64) * qJD(2)) * t61, ((-t4 * t63 - t5 * t64 + (t14 * t63 - t17 * t64) * qJD(4)) * qJD(2) - t115) * t61, t1 * t17 + t2 * t14 + t26 * t36 + t3 * t5 + t30 * t39 + t9 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109 * t65, t70 * qJD(2), 0, 0, 0, 0, 0, 0, t12, t13, 0, (t71 * t62 - t113) * qJD(2) + t116, 0, 0, 0, 0, 0, 0, t12, t13, 0, (-t26 * t61 + t73 * t62) * qJD(2) + t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t33, t22, -t47, t23, 0, -t11 * t54 - t68 + (-t42 * t111 - t84) * qJD(2), -t10 * t54 + t42 * t88 + t91, 0, 0, t47, -t33, t22, -t47, t23, 0, -t9 * t54 + t46 + (-qJD(4) * t31 - t80) * t63 + (-pkin(4) * t92 - qJD(4) * t43 - t97 * t104) * t64, -t60 * pkin(4) * t112 - t8 * t54 + (t97 * t63 + t81) * t104 + t91, (pkin(4) * qJD(4) - t114) * t88, t114 * t9 + (-t26 * t87 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35 + t76, (-qJD(4) + t54) * t88, (-t59 - t60) * t112, t74 * t104 + t30;];
tauc_reg = t6;
