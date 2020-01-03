% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRP8
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
% MMD_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:08
% EndTime: 2019-12-31 22:02:13
% DurationCPUTime: 1.26s
% Computational Cost: add. (1399->176), mult. (3546->333), div. (0->0), fcn. (2967->6), ass. (0->102)
t85 = cos(qJ(2));
t114 = t85 * qJD(2);
t81 = sin(qJ(3));
t101 = t81 * t114;
t84 = cos(qJ(3));
t118 = qJD(3) * t84;
t82 = sin(qJ(2));
t137 = t82 * t118 + t101;
t136 = -0.4e1 * t82;
t129 = t82 * t84;
t133 = pkin(6) * t81;
t89 = -t85 * pkin(2) - t82 * pkin(7);
t55 = -pkin(1) + t89;
t45 = t84 * t55;
t27 = -pkin(8) * t129 + t45 + (-pkin(3) - t133) * t85;
t127 = t84 * t85;
t69 = pkin(6) * t127;
t122 = t81 * t55 + t69;
t130 = t81 * t82;
t32 = -pkin(8) * t130 + t122;
t83 = cos(qJ(4));
t29 = t83 * t32;
t80 = sin(qJ(4));
t126 = t80 * t27 + t29;
t134 = pkin(7) + pkin(8);
t59 = t134 * t81;
t60 = t134 * t84;
t123 = -t80 * t59 + t83 * t60;
t119 = qJD(3) * t81;
t98 = t84 * t114;
t135 = -t82 * t119 + t98;
t77 = t82 ^ 2;
t92 = (-t85 ^ 2 + t77) * qJD(2);
t78 = t84 ^ 2;
t121 = t81 ^ 2 - t78;
t93 = t121 * qJD(3);
t113 = qJD(3) + qJD(4);
t88 = pkin(2) * t82 - pkin(7) * t85;
t53 = t88 * qJD(2);
t117 = qJD(3) * t85;
t105 = t81 * t117;
t74 = t82 * qJD(2);
t87 = t84 * t74 + t105;
t19 = t87 * pkin(6) - t55 * t118 - t81 * t53;
t132 = t83 * pkin(3);
t131 = t80 * t81;
t128 = t83 * t84;
t102 = t81 * t74;
t124 = pkin(6) * t102 + t84 * t53;
t54 = pkin(3) * t130 + t82 * pkin(6);
t116 = qJD(4) * t80;
t115 = qJD(4) * t83;
t112 = -0.2e1 * pkin(1) * qJD(2);
t111 = -0.2e1 * pkin(2) * qJD(3);
t110 = t80 * t130;
t73 = pkin(6) * t114;
t33 = t137 * pkin(3) + t73;
t109 = pkin(3) * t119;
t108 = pkin(3) * t116;
t107 = pkin(3) * t115;
t103 = t84 * t117;
t100 = t81 * t118;
t99 = t82 * t114;
t72 = -t84 * pkin(3) - pkin(2);
t97 = qJD(3) * t134;
t12 = (pkin(3) * t82 - pkin(8) * t127) * qJD(2) + (-t69 + (pkin(8) * t82 - t55) * t81) * qJD(3) + t124;
t14 = -pkin(8) * t137 - t19;
t96 = t83 * t12 - t80 * t14;
t95 = t83 * t27 - t80 * t32;
t94 = -t83 * t59 - t80 * t60;
t91 = 0.2e1 * t99;
t90 = t81 * t98;
t47 = t80 * t84 + t83 * t81;
t3 = -t27 * t115 + t32 * t116 - t80 * t12 - t83 * t14;
t51 = t81 * t97;
t52 = t84 * t97;
t17 = t59 * t115 + t60 * t116 + t83 * t51 + t80 * t52;
t4 = -t126 * qJD(4) + t96;
t18 = -t123 * qJD(4) + t80 * t51 - t83 * t52;
t31 = t113 * t47;
t71 = pkin(4) + t132;
t65 = -0.2e1 * t99;
t46 = -t128 + t131;
t36 = t82 * t128 - t110;
t35 = t47 * t82;
t34 = t46 * pkin(4) + t72;
t30 = t113 * t131 - t84 * t115 - t83 * t118;
t28 = t35 * pkin(4) + t54;
t23 = t31 * pkin(4) + t109;
t22 = -t46 * qJ(5) + t123;
t21 = -t47 * qJ(5) + t94;
t20 = -t122 * qJD(3) + t124;
t16 = -qJD(4) * t110 + (t113 * t129 + t101) * t83 + t135 * t80;
t15 = t80 * t101 + t31 * t82 - t83 * t98;
t9 = t16 * pkin(4) + t33;
t8 = -t35 * qJ(5) + t126;
t7 = -t85 * pkin(4) - t36 * qJ(5) + t95;
t6 = t30 * qJ(5) - t47 * qJD(5) + t18;
t5 = -t31 * qJ(5) - t46 * qJD(5) - t17;
t2 = -t16 * qJ(5) - t35 * qJD(5) - t3;
t1 = pkin(4) * t74 + t15 * qJ(5) - t36 * qJD(5) + t4;
t10 = [0, 0, 0, t91, -0.2e1 * t92, 0, 0, 0, t82 * t112, t85 * t112, -0.2e1 * t77 * t100 + 0.2e1 * t78 * t99, t90 * t136 + 0.2e1 * t77 * t93, 0.2e1 * t82 * t105 + 0.2e1 * t84 * t92, 0.2e1 * t82 * t103 - 0.2e1 * t81 * t92, t65, 0.2e1 * t45 * t74 - 0.2e1 * t20 * t85 + 0.2e1 * (t77 * t118 + t81 * t99) * pkin(6), -0.2e1 * t19 * t85 - 0.2e1 * t122 * t74 + 0.2e1 * (-t77 * t119 + t84 * t91) * pkin(6), -0.2e1 * t36 * t15, 0.2e1 * t15 * t35 - 0.2e1 * t36 * t16, 0.2e1 * t15 * t85 + 0.2e1 * t36 * t74, 0.2e1 * t16 * t85 - 0.2e1 * t35 * t74, t65, 0.2e1 * t54 * t16 + 0.2e1 * t33 * t35 - 0.2e1 * t4 * t85 + 0.2e1 * t95 * t74, -0.2e1 * t126 * t74 - 0.2e1 * t54 * t15 - 0.2e1 * t3 * t85 + 0.2e1 * t33 * t36, -0.2e1 * t1 * t36 + 0.2e1 * t7 * t15 - 0.2e1 * t8 * t16 - 0.2e1 * t2 * t35, 0.2e1 * t7 * t1 + 0.2e1 * t8 * t2 + 0.2e1 * t28 * t9; 0, 0, 0, 0, 0, t114, -t74, 0, -t73, pkin(6) * t74, -t82 * t93 + t90, t100 * t136 - t121 * t114, t102 - t103, t87, 0, (pkin(7) * t127 + (-pkin(2) * t84 + t133) * t82) * qJD(3) + (t89 * t81 - t69) * qJD(2), (pkin(6) * t129 + t88 * t81) * qJD(3) + (t85 * t133 + t89 * t84) * qJD(2), -t15 * t47 - t36 * t30, t15 * t46 - t47 * t16 + t30 * t35 - t36 * t31, t30 * t85 + t47 * t74, t31 * t85 - t46 * t74, 0, t35 * t109 + t72 * t16 - t18 * t85 + t54 * t31 + t33 * t46 + t94 * t74, t36 * t109 - t123 * t74 - t72 * t15 - t17 * t85 - t54 * t30 + t33 * t47, -t1 * t47 + t21 * t15 - t22 * t16 - t2 * t46 + t7 * t30 - t8 * t31 - t5 * t35 - t6 * t36, t1 * t21 + t2 * t22 + t28 * t23 + t9 * t34 + t8 * t5 + t7 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t100, -0.2e1 * t93, 0, 0, 0, t81 * t111, t84 * t111, -0.2e1 * t47 * t30, 0.2e1 * t30 * t46 - 0.2e1 * t47 * t31, 0, 0, 0, 0.2e1 * t46 * t109 + 0.2e1 * t72 * t31, 0.2e1 * t47 * t109 - 0.2e1 * t72 * t30, 0.2e1 * t21 * t30 - 0.2e1 * t22 * t31 - 0.2e1 * t5 * t46 - 0.2e1 * t6 * t47, 0.2e1 * t21 * t6 + 0.2e1 * t22 * t5 + 0.2e1 * t34 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, -t137, t74, t20, t19, 0, 0, -t15, -t16, t74, t74 * t132 + (-t29 + (t85 * pkin(3) - t27) * t80) * qJD(4) + t96, (t85 * t115 - t80 * t74) * pkin(3) + t3, t71 * t15 + (-t16 * t80 + (-t35 * t83 + t36 * t80) * qJD(4)) * pkin(3), t1 * t71 + (t2 * t80 + (-t7 * t80 + t8 * t83) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, -t119, 0, -pkin(7) * t118, pkin(7) * t119, 0, 0, -t30, -t31, 0, t18, t17, t71 * t30 + (-t31 * t80 + (-t46 * t83 + t47 * t80) * qJD(4)) * pkin(3), t6 * t71 + (t5 * t80 + (-t21 * t80 + t22 * t83) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t108, -0.2e1 * t107, 0, 0.2e1 * (-t71 + t132) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, t74, t4, t3, pkin(4) * t15, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t31, 0, t18, t17, pkin(4) * t30, t6 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t107, 0, -pkin(4) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
