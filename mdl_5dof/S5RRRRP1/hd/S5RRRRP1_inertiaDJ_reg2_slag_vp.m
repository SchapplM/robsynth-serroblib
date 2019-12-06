% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:46:00
% EndTime: 2019-12-05 18:46:05
% DurationCPUTime: 1.32s
% Computational Cost: add. (2446->146), mult. (5593->260), div. (0->0), fcn. (5100->6), ass. (0->93)
t120 = -pkin(7) - pkin(6);
t71 = sin(qJ(2));
t118 = cos(qJ(3));
t92 = t118 * t120;
t89 = t71 * t92;
t116 = sin(qJ(3));
t72 = cos(qJ(2));
t99 = t116 * t72;
t122 = t120 * t99 + t89;
t125 = t122 * qJD(3);
t101 = t118 * t72;
t91 = t116 * t120;
t57 = t71 * t91;
t123 = t120 * t101 - t57;
t124 = t123 * qJD(3);
t100 = t116 * t71;
t96 = t118 * qJD(3);
t98 = qJD(3) * t116;
t121 = -(t118 * qJD(2) + t96) * t72 + qJD(2) * t100 + t71 * t98;
t70 = sin(qJ(4));
t119 = pkin(3) * t70;
t117 = cos(qJ(4));
t56 = t118 * t71 + t99;
t75 = (-qJD(2) - qJD(3)) * t56;
t85 = t100 - t101;
t82 = t70 * t85;
t97 = qJD(4) * t117;
t15 = -qJD(4) * t82 - t117 * t75 - t121 * t70 + t56 * t97;
t81 = t117 * t85;
t45 = t56 * t70 + t81;
t93 = pkin(3) * t97;
t115 = -t15 * t119 - t45 * t93;
t77 = -t56 * pkin(8) + t122;
t28 = t70 * t77;
t34 = -t85 * pkin(8) - t123;
t20 = t117 * t34 + t28;
t106 = t118 * pkin(2);
t67 = t106 + pkin(3);
t95 = pkin(2) * t96;
t39 = -t67 * t97 - t117 * t95 + (t116 * qJD(4) + t98) * t70 * pkin(2);
t90 = t117 * t116;
t54 = pkin(2) * t90 + t70 * t67;
t114 = -t39 * t119 + t54 * t93;
t112 = qJD(4) * t70;
t111 = t71 * qJD(2);
t110 = t72 * qJD(2);
t109 = -0.2e1 * pkin(1) * qJD(2);
t76 = (-qJD(4) * t90 + (-t118 * t70 - t90) * qJD(3)) * pkin(2);
t40 = -t67 * t112 + t76;
t46 = t117 * t56 - t82;
t108 = -t54 * t15 + t39 * t45 - t40 * t46;
t69 = pkin(2) * t111;
t107 = pkin(3) * t112;
t105 = t117 * pkin(3);
t104 = t116 * pkin(2);
t103 = t46 * t112;
t102 = t71 * t110;
t68 = -t72 * pkin(2) - pkin(1);
t29 = t117 * t77;
t19 = -t34 * t70 + t29;
t94 = pkin(2) * t98;
t53 = -t70 * t104 + t117 * t67;
t87 = qJD(2) * t92;
t86 = qJD(2) * t91;
t73 = t75 * pkin(8) + t71 * t87 + t72 * t86 + t125;
t74 = pkin(8) * t121 - t71 * t86 + t72 * t87 + t124;
t3 = -qJD(4) * t29 + t34 * t112 - t117 * t73 - t70 * t74;
t83 = t68 * t56;
t50 = t85 * pkin(3) + t68;
t38 = -t75 * pkin(3) + t69;
t4 = -qJD(4) * t28 + t117 * t74 - t34 * t97 - t70 * t73;
t66 = t105 + pkin(4);
t65 = -0.2e1 * t93;
t64 = -0.2e1 * t107;
t51 = pkin(4) + t53;
t36 = 0.2e1 * t40;
t35 = 0.2e1 * t39;
t32 = -t93 + t39;
t31 = (-pkin(3) - t67) * t112 + t76;
t26 = t54 * t39;
t25 = t45 * pkin(4) + t50;
t24 = t124 + (t72 * t92 - t57) * qJD(2);
t23 = -t125 + (-t72 * t91 - t89) * qJD(2);
t14 = qJD(4) * t81 + t56 * t112 + t117 * t121 - t70 * t75;
t10 = -qJ(5) * t45 + t20;
t9 = -qJ(5) * t46 + t19;
t8 = -0.2e1 * t46 * t14;
t7 = 0.2e1 * t45 * t15;
t6 = t15 * pkin(4) + t38;
t5 = 0.2e1 * t14 * t45 - 0.2e1 * t15 * t46;
t2 = t14 * qJ(5) - t46 * qJD(5) + t4;
t1 = t15 * qJ(5) + t45 * qJD(5) + t3;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t102, 0.2e1 * (-t71 ^ 2 + t72 ^ 2) * qJD(2), 0, -0.2e1 * t102, 0, 0, t71 * t109, t72 * t109, 0, 0, -0.2e1 * t56 * t121, 0.2e1 * t121 * t85 + 0.2e1 * t56 * t75, 0, -0.2e1 * t85 * t75, 0, 0, 0.2e1 * qJD(3) * t83 + 0.2e1 * (t71 * pkin(2) * t85 + t83) * qJD(2), -0.2e1 * t121 * t68 + 0.2e1 * t56 * t69, 0.2e1 * t121 * t122 - 0.2e1 * t123 * t75 + 0.2e1 * t23 * t85 - 0.2e1 * t24 * t56, 0.2e1 * t122 * t24 + 0.2e1 * t123 * t23 + 0.2e1 * t68 * t69, t8, t5, 0, t7, 0, 0, 0.2e1 * t15 * t50 + 0.2e1 * t38 * t45, -0.2e1 * t14 * t50 + 0.2e1 * t38 * t46, 0.2e1 * t14 * t19 - 0.2e1 * t15 * t20 + 0.2e1 * t3 * t45 - 0.2e1 * t4 * t46, 0.2e1 * t19 * t4 - 0.2e1 * t20 * t3 + 0.2e1 * t38 * t50, t8, t5, 0, t7, 0, 0, 0.2e1 * t15 * t25 + 0.2e1 * t45 * t6, -0.2e1 * t14 * t25 + 0.2e1 * t46 * t6, 0.2e1 * t1 * t45 - 0.2e1 * t10 * t15 + 0.2e1 * t14 * t9 - 0.2e1 * t2 * t46, -0.2e1 * t1 * t10 + 0.2e1 * t2 * t9 + 0.2e1 * t25 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, 0, -t111, 0, -pkin(6) * t110, pkin(6) * t111, 0, 0, 0, 0, -t121, 0, t75, 0, t24, t23, t75 * t104 + t106 * t121 + t56 * t94 - t85 * t95, (t118 * t24 - t23 * t116 + (-t116 * t122 - t118 * t123) * qJD(3)) * pkin(2), 0, 0, -t14, 0, -t15, 0, t4, t3, t14 * t53 + t108, t19 * t40 - t20 * t39 - t3 * t54 + t4 * t53, 0, 0, -t14, 0, -t15, 0, t2, t1, t14 * t51 + t108, -t1 * t54 - t10 * t39 + t2 * t51 + t40 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t94, -0.2e1 * t95, 0, 0, 0, 0, 0, 0, 0, 0, t36, t35, 0, 0.2e1 * t40 * t53 - 0.2e1 * t26, 0, 0, 0, 0, 0, 0, t36, t35, 0, 0.2e1 * t40 * t51 - 0.2e1 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, 0, t75, 0, t24, t23, 0, 0, 0, 0, -t14, 0, -t15, 0, t4, t3, (t117 * t14 + t103) * pkin(3) + t115, (t117 * t4 - t3 * t70 + (t117 * t20 - t19 * t70) * qJD(4)) * pkin(3), 0, 0, -t14, 0, -t15, 0, t2, t1, pkin(3) * t103 + t14 * t66 + t115, t2 * t66 + (-t1 * t70 + (t10 * t117 - t70 * t9) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t95, 0, 0, 0, 0, 0, 0, 0, 0, t31, t32, 0, (-t112 * t53 + t117 * t40) * pkin(3) + t114, 0, 0, 0, 0, 0, 0, t31, t32, 0, -t107 * t51 + t40 * t66 + t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t65, 0, 0, 0, 0, 0, 0, 0, 0, t64, t65, 0, 0.2e1 * (t105 - t66) * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, -t15, 0, t4, t3, 0, 0, 0, 0, -t14, 0, -t15, 0, t2, t1, t14 * pkin(4), t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t39, 0, 0, 0, 0, 0, 0, 0, 0, t40, t39, 0, t40 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t93, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t93, 0, -pkin(4) * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
