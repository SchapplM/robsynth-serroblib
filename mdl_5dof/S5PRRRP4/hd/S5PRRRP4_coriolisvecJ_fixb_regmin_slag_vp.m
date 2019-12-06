% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:46:20
% EndTime: 2019-12-05 16:46:22
% DurationCPUTime: 0.53s
% Computational Cost: add. (864->132), mult. (1698->163), div. (0->0), fcn. (1090->6), ass. (0->100)
t122 = 2 * qJD(4);
t60 = sin(qJ(2));
t86 = qJD(1) * qJD(2);
t94 = qJD(1) * t60;
t121 = qJD(3) * t94 + t60 * t86;
t63 = cos(qJ(2));
t88 = t63 * qJD(1);
t48 = qJD(2) * pkin(2) + t88;
t62 = cos(qJ(3));
t78 = t63 * t86;
t59 = sin(qJ(3));
t97 = t121 * t59;
t13 = (qJD(3) * t48 + t78) * t62 - t97;
t61 = cos(qJ(4));
t10 = t61 * t13;
t31 = t59 * t48 + t62 * t94;
t55 = qJD(2) + qJD(3);
t25 = t55 * pkin(7) + t31;
t58 = sin(qJ(4));
t105 = t58 * t25;
t4 = t10 + (qJD(5) - t105) * qJD(4);
t9 = t58 * t13;
t91 = qJD(4) * t61;
t6 = t25 * t91 + t9;
t120 = t4 * t61 + t6 * t58;
t39 = t59 * t60 - t62 * t63;
t119 = t39 * t55;
t40 = t59 * t63 + t62 * t60;
t35 = t40 * qJD(1);
t93 = qJD(3) * t59;
t83 = pkin(2) * t93;
t118 = (-t35 + t83) * t55;
t49 = t59 * t94;
t36 = t62 * t88 - t49;
t114 = t62 * pkin(2);
t84 = qJD(3) * t114;
t117 = -t36 + t84;
t64 = qJD(4) ^ 2;
t116 = pkin(7) * t64;
t115 = t55 * pkin(3);
t113 = t119 * t55;
t19 = t55 * t40;
t112 = t19 * t55;
t30 = t62 * t48 - t49;
t111 = t30 * t55;
t110 = t31 * t55;
t42 = -t61 * pkin(4) - t58 * qJ(5) - pkin(3);
t109 = t42 * t55;
t50 = t59 * pkin(2) + pkin(7);
t108 = t50 * t64;
t107 = t55 * t58;
t106 = t55 * t61;
t104 = t61 * t25;
t103 = t64 * t58;
t14 = t121 * t62 + t48 * t93 + t59 * t78;
t24 = -t30 - t115;
t102 = t14 * t58 + t24 * t91;
t92 = qJD(4) * t58;
t101 = t31 * t106 + t30 * t92;
t87 = qJD(4) * qJ(5);
t89 = t58 * qJD(5);
t33 = pkin(4) * t92 - t61 * t87 - t89;
t27 = t33 + t83;
t100 = t27 - t35;
t99 = t35 * t106 + t36 * t92;
t98 = t31 - t33;
t56 = t58 ^ 2;
t57 = t61 ^ 2;
t96 = t56 - t57;
t95 = t56 + t57;
t17 = t87 + t104;
t90 = t17 * qJD(4);
t54 = t55 ^ 2;
t85 = t58 * t54 * t61;
t82 = t55 * t92;
t74 = pkin(4) * t58 - qJ(5) * t61;
t7 = (t74 * qJD(4) - t89) * t55 + t14;
t81 = -t7 - t116;
t80 = -t7 - t108;
t76 = -qJD(4) * pkin(4) + qJD(5);
t16 = t76 + t105;
t73 = t16 * t58 + t17 * t61;
t72 = t40 * t64 + t112;
t71 = t119 * t122;
t37 = t42 - t114;
t70 = t37 * t55 - t84;
t69 = (-pkin(3) - t114) * t55 - t84;
t68 = t16 * t91 - t58 * t90 + t120;
t66 = (t16 * t61 - t17 * t58) * qJD(4) + t120;
t65 = qJD(2) ^ 2;
t53 = t64 * t61;
t41 = 0.2e1 * t61 * t82;
t34 = t74 * t55;
t32 = -0.2e1 * t96 * t55 * qJD(4);
t21 = t24 * t92;
t12 = -t30 + t109;
t8 = t12 * t92;
t2 = t58 * t71 - t72 * t61;
t1 = t72 * t58 + t61 * t71;
t3 = [0, 0, -t65 * t60, -t65 * t63, 0, -t112, t113, 0, 0, 0, 0, 0, t2, t1, t2, -t95 * t113, -t1, -t119 * t73 + t12 * t19 + t7 * t39 + t66 * t40; 0, 0, 0, 0, 0, -t14 - t118, t36 * t55 + (-t78 + (-pkin(2) * t55 - t48) * qJD(3)) * t62 + t97, t41, t32, t53, -t103, 0, t21 + t69 * t92 + (-t55 * t83 - t108 - t14) * t61 + t99, (t108 + t118) * t58 + (t36 + t69) * t91 + t102, t8 + t70 * t92 + (-t27 * t55 + t80) * t61 + t99, t117 * t55 * t95 + t68, (-t100 * t55 + t80) * t58 + (-t12 - t36 - t70) * t91, t100 * t12 + t117 * t73 + t7 * t37 + t66 * t50; 0, 0, 0, 0, 0, -t14 + t110, t111 - t13, t41, t32, t53, -t103, 0, -pkin(3) * t82 + t21 + (-t14 - t116) * t61 + t101, (-t110 + t116) * t58 + (t30 - t115) * t91 + t102, t42 * t82 + t8 + (-t33 * t55 + t81) * t61 + t101, -t95 * t111 + t68, (-t12 - t30 - t109) * t91 + (t98 * t55 + t81) * t58, t66 * pkin(7) - t98 * t12 - t73 * t30 + t7 * t42; 0, 0, 0, 0, 0, 0, 0, -t85, t96 * t54, 0, 0, 0, -t24 * t107 - t9, -t24 * t106 - t10, -t9 + (-t12 * t58 + t34 * t61) * t55, ((t17 - t87) * t58 + (-t16 + t76) * t61) * t55, qJD(5) * t122 + t10 + (t12 * t61 + t34 * t58) * t55, -t16 * t104 - t6 * pkin(4) + t4 * qJ(5) - t12 * t34 + (qJD(5) + t105) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, 0, -t56 * t54 - t64, t12 * t107 + t6 - t90;];
tauc_reg = t3;
