% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:33
% EndTime: 2019-12-05 17:51:36
% DurationCPUTime: 0.70s
% Computational Cost: add. (1482->130), mult. (2995->198), div. (0->0), fcn. (1736->8), ass. (0->101)
t56 = cos(pkin(8)) * pkin(1) + pkin(2);
t50 = t56 * qJD(1);
t69 = sin(qJ(3));
t71 = cos(qJ(3));
t120 = sin(pkin(8)) * pkin(1);
t87 = qJD(1) * t120;
t36 = t69 * t50 + t71 * t87;
t121 = t36 * qJD(3);
t101 = t71 * t120 + t69 * t56;
t41 = t101 * qJD(3);
t61 = qJD(1) + qJD(3);
t123 = -t41 * t61 - t121;
t122 = t36 * t61 - t121;
t53 = t69 * t87;
t35 = t71 * t50 - t53;
t77 = qJD(4) - t35;
t66 = cos(pkin(9));
t70 = cos(qJ(5));
t106 = t66 * t70;
t95 = qJD(3) * t71;
t102 = qJD(3) * t53 - t50 * t95;
t27 = t61 * qJD(4) - t102;
t68 = sin(qJ(5));
t64 = sin(pkin(9));
t48 = -t66 * pkin(4) - t64 * pkin(7) - pkin(3);
t14 = t48 * t61 + t77;
t31 = t61 * qJ(4) + t36;
t17 = t64 * qJD(2) + t66 * t31;
t7 = t70 * t14 - t68 * t17;
t2 = t7 * qJD(5) + t27 * t106 + t121 * t68;
t119 = t2 * t68;
t59 = t64 ^ 2;
t110 = t59 * t70;
t118 = t27 * t110 + t2 * t66;
t16 = -t66 * qJD(2) + t64 * t31;
t117 = t16 * t64;
t116 = t27 * t66;
t54 = t69 * t120;
t40 = -qJD(3) * t54 + t56 * t95;
t39 = qJD(4) + t40;
t114 = t39 * t61;
t21 = t59 * t27;
t58 = t61 ^ 2;
t112 = t59 * t58;
t111 = t59 * t61;
t60 = t66 ^ 2;
t22 = t60 * t27;
t109 = t61 * t64;
t108 = t66 * t61;
t107 = t66 * t68;
t97 = qJ(4) * t66;
t38 = t68 * t48 + t70 * t97;
t92 = qJD(5) * t38;
t94 = qJD(4) * t66;
t105 = -t35 * t107 + t70 * t36 + t68 * t94 + t92;
t37 = t70 * t48 - t68 * t97;
t93 = qJD(5) * t37;
t104 = -t35 * t106 - t68 * t36 + t70 * t94 + t93;
t103 = t21 + t22;
t100 = t59 + t60;
t99 = t68 ^ 2 - t70 ^ 2;
t98 = qJ(4) * t27;
t8 = t68 * t14 + t70 * t17;
t96 = qJD(5) * t8;
t91 = qJD(5) * t68;
t90 = qJD(5) * t70;
t89 = qJ(4) * qJD(5);
t88 = t16 * t109;
t86 = qJD(5) * t111;
t85 = t64 * t91;
t84 = t64 * t90;
t83 = t71 * t56 - t54;
t82 = t58 * t68 * t110;
t52 = -qJD(5) + t108;
t81 = t52 * t85;
t3 = -t27 * t107 + t121 * t70 - t96;
t80 = t16 * t84 + t68 * t21 - t3 * t66;
t79 = (-qJD(5) - t52) * t109;
t78 = t68 * t7 - t70 * t8;
t76 = t70 * t68 * t86;
t75 = t17 * t66 + t117;
t32 = t48 - t83;
t42 = qJ(4) + t101;
t10 = t42 * t106 + t68 * t32;
t9 = -t42 * t107 + t70 * t32;
t74 = -t119 + (-t3 - t96) * t70;
t73 = -t52 ^ 2 - t112;
t46 = t85 * t108;
t44 = -0.2e1 * t76;
t43 = 0.2e1 * t76;
t33 = 0.2e1 * t99 * t86;
t29 = t121 * t64;
t28 = -t61 * pkin(3) + t77;
t26 = (t52 + t108) * t84;
t25 = t46 + t81;
t20 = t59 * t98;
t13 = t42 * t21;
t6 = t7 * t85;
t5 = -qJD(5) * t10 - t39 * t107 + t70 * t41;
t4 = qJD(5) * t9 + t39 * t106 + t68 * t41;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, -t40 * t61 + t102, 0, -t102 * t101 - t121 * t83 - t35 * t41 + t36 * t40, 0, 0, 0, 0, 0, 0, t123 * t66, t41 * t109 + t29, t100 * t114 + t103, t42 * t22 + t13 + t121 * (-pkin(3) - t83) + t28 * t41 + t75 * t39, t44, t33, t25, t43, t26, 0, -t5 * t52 + (t39 * t68 + t42 * t90) * t111 + t80, t110 * t114 + t4 * t52 + (-t42 * t111 - t117) * t91 + t118, t6 + ((-t4 * t68 - t5 * t70 + (-t10 * t70 + t68 * t9) * qJD(5)) * t61 + t74) * t64, t2 * t10 + t39 * t117 + t3 * t9 + t8 * t4 + t7 * t5 + t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t52 - t108) * t84, t46 - t81, 0, (t2 * t70 - t116 - t3 * t68 + (-t68 * t8 - t7 * t70) * qJD(5)) * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, t35 * t61 + t102, 0, 0, 0, 0, 0, 0, 0, 0, t122 * t66, -t36 * t109 + t29, t100 * t61 * t77 + t103, -pkin(3) * t121 - t28 * t36 + t60 * t98 + t75 * t77 + t20, t44, t33, t25, t43, t26, 0, t105 * t52 + (t77 * t68 + t70 * t89) * t111 + t80, -t16 * t85 + t104 * t52 + (-t68 * t89 + t77 * t70) * t111 + t118, t6 + (((-t92 + t105) * t70 + (t93 - t104) * t68) * t61 + t74) * t64, t104 * t8 - t105 * t7 + t77 * t117 + t2 * t38 + t3 * t37 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100 * t58, -t75 * t61 + t121, 0, 0, 0, 0, 0, 0, t73 * t68, t73 * t70, 0, t119 + t3 * t70 - t78 * qJD(5) + (t78 * t66 - t117) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t99 * t112, t68 * t79, -t82, t70 * t79, 0, -t8 * t52 - t70 * t88 + t3, -t7 * t52 + (-qJD(5) * t14 - t116) * t70 + (qJD(5) * t17 - t121 + t88) * t68, 0, 0;];
tauc_reg = t1;
