% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:53
% EndTime: 2019-12-31 17:59:57
% DurationCPUTime: 1.02s
% Computational Cost: add. (1351->175), mult. (2738->261), div. (0->0), fcn. (1511->6), ass. (0->104)
t35 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(6);
t27 = t35 * qJD(1) + qJD(3);
t45 = sin(qJ(4));
t47 = cos(qJ(4));
t57 = t45 * qJD(2) - t47 * t27;
t15 = t57 * qJD(4);
t65 = pkin(4) * t47 + pkin(7) * t45;
t28 = t65 * qJD(4) + qJD(3);
t23 = t28 * qJD(1);
t44 = sin(qJ(5));
t46 = cos(qJ(5));
t20 = t47 * qJD(2) + t45 * t27;
t14 = qJD(4) * pkin(7) + t20;
t37 = sin(pkin(8)) * pkin(1) + qJ(3);
t25 = t45 * pkin(4) - t47 * pkin(7) + t37;
t22 = t25 * qJD(1);
t62 = t44 * t14 - t46 * t22;
t1 = -t62 * qJD(5) - t46 * t15 + t44 * t23;
t84 = t45 * qJD(1);
t36 = qJD(5) + t84;
t125 = t62 * t36 + t1;
t85 = t44 * qJD(4);
t78 = t45 * t85;
t87 = qJD(5) * t47;
t124 = t46 * t87 - t78;
t6 = t46 * t14 + t44 * t22;
t2 = -qJD(5) * t6 + t44 * t15 + t46 * t23;
t123 = -t6 * t36 - t2;
t74 = t44 * t87;
t83 = t46 * qJD(4);
t53 = t45 * t83 + t74;
t17 = t53 * qJD(1) - qJD(5) * t83;
t94 = qJD(1) * t47;
t72 = t46 * t94;
t31 = t72 + t85;
t92 = qJD(4) * t47;
t98 = -t17 * t45 + t31 * t92;
t90 = qJD(5) * t31;
t18 = -qJD(1) * t78 + t90;
t41 = t47 ^ 2;
t95 = qJD(1) * t41;
t122 = -(-t36 * t45 + t95) * t83 + t36 * t74;
t16 = t20 * qJD(4);
t63 = -t44 * t62 - t46 * t6;
t121 = qJD(4) * t63 + t16;
t112 = t16 * t47;
t120 = (-t20 * t47 - t45 * t57) * qJD(4) + t15 * t45 + t112;
t119 = 0.2e1 * qJD(3);
t13 = -qJD(4) * pkin(4) + t57;
t116 = t13 * t44;
t115 = t13 * t46;
t114 = t16 * t44;
t113 = t16 * t46;
t29 = t44 * t94 - t83;
t111 = t29 * t36;
t110 = t29 * t47;
t109 = t31 * t29;
t108 = t31 * t36;
t107 = t36 * t44;
t106 = t36 * t46;
t105 = t44 * t45;
t104 = t45 * t18;
t103 = t45 * t46;
t102 = t47 * t17;
t101 = t47 * t18;
t48 = qJD(4) ^ 2;
t100 = t48 * t45;
t99 = t48 * t47;
t97 = t45 ^ 2 - t41;
t49 = qJD(1) ^ 2;
t96 = -t48 - t49;
t33 = qJD(1) * t37;
t93 = qJD(4) * t45;
t91 = qJD(5) * t29;
t89 = qJD(5) * t44;
t88 = qJD(5) * t46;
t86 = t33 * qJD(1);
t82 = qJD(1) * qJD(4);
t81 = t47 * t49 * t45;
t80 = t44 * t95;
t79 = t35 * t92;
t77 = t29 * t93;
t76 = t31 * t93;
t75 = t31 * t87;
t71 = t47 * t82;
t70 = t36 + t84;
t69 = -t17 + t91;
t67 = qJD(5) * t45 + qJD(1);
t66 = t45 * t71;
t64 = t44 * t6 - t46 * t62;
t56 = t33 * t119;
t55 = t124 * t36;
t10 = t35 * t103 + t44 * t25;
t9 = -t35 * t105 + t46 * t25;
t54 = -pkin(7) * t92 + t13 * t45;
t52 = -t64 * qJD(5) + t1 * t46 - t2 * t44;
t50 = qJD(4) * t13 + t52;
t32 = t65 * qJD(1);
t11 = t46 * t101;
t8 = t44 * t32 - t46 * t57;
t7 = t46 * t32 + t44 * t57;
t4 = -t10 * qJD(5) + t46 * t28 - t44 * t79;
t3 = t9 * qJD(5) + t44 * t28 + t46 * t79;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1) * t119, t56, -0.2e1 * t66, 0.2e1 * t97 * t82, -t100, 0.2e1 * t66, -t99, 0, t33 * t92 - t35 * t100 + (t45 * t119 + t37 * t92) * qJD(1), -t33 * t93 - t35 * t99 + (t47 * t119 - t37 * t93) * qJD(1), t120, -t120 * t35 + t56, -t46 * t102 - t53 * t31, -t11 + (-t75 + t77) * t46 + (t76 + (t17 + t91) * t47) * t44, -t122 + t98, t44 * t101 + t124 * t29, -t104 + (-t80 - t110) * qJD(4) - t55, t70 * t92, t4 * t36 + (t2 + (t29 * t35 - t116) * qJD(4)) * t45 + (t13 * t88 + t114 - t18 * t35 + (qJD(1) * t9 - t62) * qJD(4)) * t47, -t3 * t36 + (-t1 + (t31 * t35 - t115) * qJD(4)) * t45 + (-t13 * t89 + t113 + t17 * t35 + (-qJD(1) * t10 - t6) * qJD(4)) * t47, -t10 * t18 + t9 * t17 - t3 * t29 - t4 * t31 + t64 * t93 + (t63 * qJD(5) - t1 * t44 - t2 * t46) * t47, t1 * t10 + t2 * t9 + t6 * t3 - t62 * t4 + (t13 * t93 - t112) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, t100, 0, -t15 * t47 + t16 * t45 + (-t20 * t45 + t47 * t57) * qJD(4), 0, 0, 0, 0, 0, 0, t104 + (-t80 + t110) * qJD(4) - t55, t122 + t98, -t11 + (t75 + t77) * t46 + (t69 * t47 - t76) * t44, t121 * t45 + t50 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t86, 0, 0, 0, 0, 0, 0, t96 * t45, t96 * t47, 0, -t120 - t86, 0, 0, 0, 0, 0, 0, -t101 - t67 * t106 + (-t70 * t47 * t44 + t29 * t45) * qJD(4), t102 + t67 * t107 + (-t47 * t106 + (t31 - t72) * t45) * qJD(4), (-t29 * t92 + t67 * t31 - t104) * t46 + (t67 * t29 + t98) * t44, -t64 * qJD(1) - t121 * t47 + t50 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t97 * t49, 0, -t81, 0, 0, -t47 * t86, t33 * t84, 0, 0, t31 * t106 - t17 * t44, (-t17 - t111) * t46 + (-t18 - t108) * t44, t36 * t88 + (t36 * t103 + (-t31 + t85) * t47) * qJD(1), t29 * t107 - t18 * t46, -t36 * t89 + (-t36 * t105 + (t29 + t83) * t47) * qJD(1), -t36 * t94, -pkin(4) * t18 - t113 - t20 * t29 - t7 * t36 + (-pkin(7) * t106 + t116) * qJD(5) + (t54 * t44 + t47 * t62) * qJD(1), pkin(4) * t17 + t114 - t20 * t31 + t8 * t36 + (pkin(7) * t107 + t115) * qJD(5) + (t54 * t46 + t47 * t6) * qJD(1), t8 * t29 + t7 * t31 + ((-t18 + t90) * pkin(7) + t125) * t46 + (t69 * pkin(7) + t123) * t44, -t16 * pkin(4) + pkin(7) * t52 - t13 * t20 - t6 * t8 + t62 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, -t29 ^ 2 + t31 ^ 2, t111 - t17, -t109, t108 - t18, t71, -t13 * t31 - t123, t13 * t29 - t125, 0, 0;];
tauc_reg = t5;
