% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPP1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:31
% EndTime: 2019-12-31 19:24:35
% DurationCPUTime: 1.29s
% Computational Cost: add. (1211->190), mult. (3673->364), div. (0->0), fcn. (3162->6), ass. (0->109)
t71 = sin(pkin(5));
t126 = -0.2e1 * t71;
t72 = cos(pkin(8));
t74 = sin(qJ(2));
t113 = t74 * t72;
t70 = sin(pkin(8));
t73 = cos(pkin(5));
t75 = cos(qJ(2));
t41 = (t113 * t73 + t70 * t75) * qJD(2);
t108 = qJD(2) * t75;
t109 = qJD(2) * t74;
t98 = t70 * t109;
t42 = t108 * t72 - t73 * t98;
t114 = t73 * t75;
t46 = -t114 * t72 + t70 * t74;
t47 = t114 * t70 + t113;
t77 = 0.2e1 * t41 * t47 + 0.2e1 * t42 * t46;
t110 = qJ(3) * t71;
t54 = -pkin(2) * t75 - t110 * t74 - pkin(1);
t92 = qJ(3) * t73 + pkin(7);
t55 = t92 * t74;
t125 = (t54 * t71 - t55 * t73) * t72;
t80 = (t41 * t70 - t42 * t72) * t71;
t118 = t70 * t73;
t119 = t70 * t71;
t105 = qJD(3) * t74;
t34 = -t71 * t105 + (pkin(2) * t74 - t110 * t75) * qJD(2);
t104 = qJD(3) * t75;
t84 = qJD(2) * t92;
t35 = t104 * t73 - t74 * t84;
t36 = -t105 * t73 - t75 * t84;
t10 = t118 * t36 + t119 * t34 + t35 * t72;
t5 = -(qJ(4) * t109 - qJD(4) * t75) * t71 - t10;
t124 = 0.2e1 * t71;
t123 = t34 * t72;
t106 = qJD(3) * t71;
t53 = qJD(4) * t73 + t106 * t72;
t122 = t53 * t46;
t121 = t53 * t73;
t120 = t53 * t75;
t117 = t71 * t72;
t116 = t71 * t75;
t115 = t72 * t73;
t112 = pkin(3) + qJ(5);
t56 = t92 * t75;
t50 = t70 * t56;
t111 = pkin(3) * t116 + t50;
t49 = pkin(2) * t118 + t72 * t110;
t69 = t71 ^ 2;
t107 = qJD(3) * t69;
t103 = qJD(4) * t69;
t102 = qJD(4) * t70;
t101 = t47 * qJD(4);
t20 = 0.2e1 * t46 * t41;
t21 = 0.2e1 * t47 * t42;
t100 = -0.2e1 * pkin(1) * qJD(2);
t99 = t41 * t117;
t33 = t42 * t119;
t16 = -t55 * t118 + t119 * t54 + t72 * t56;
t97 = t69 * t104;
t96 = t73 * t106;
t65 = t71 * t109;
t95 = t74 * t108;
t94 = -pkin(2) * t72 - pkin(3);
t93 = -qJ(4) * t70 - pkin(2);
t17 = t34 * t73 - t71 * t36;
t24 = t73 * t54 + t71 * t55;
t91 = -0.2e1 * t96;
t90 = -0.2e1 * t95;
t89 = t70 * t97;
t30 = t70 * t35;
t88 = -t115 * t36 + t30;
t37 = -qJ(4) * t73 - t49;
t83 = -t47 * qJ(4) + t24;
t82 = t109 * t46 - t41 * t75;
t81 = t109 * t47 - t42 * t75;
t13 = qJ(4) * t116 - t16;
t78 = t109 * t69 * t72 - t41 * t73;
t22 = t42 * t73 + t69 * t98;
t76 = -t42 * qJ(4) - t101 + t17;
t68 = t70 ^ 2;
t63 = t70 * t110;
t62 = t70 * t106;
t61 = t68 * t107;
t58 = t73 * t65;
t57 = t69 * t90;
t52 = -qJD(5) * t73 + t62;
t48 = pkin(2) * t115 - t63;
t45 = (-qJD(5) * t72 - t102) * t71;
t39 = (-pkin(3) * t72 + t93) * t71;
t38 = t73 * t94 + t63;
t29 = (-t112 * t72 + t93) * t71;
t28 = pkin(4) * t117 - t37;
t25 = pkin(4) * t119 + t63 + (-qJ(5) + t94) * t73;
t19 = t81 * t124;
t18 = t82 * t124;
t15 = -t50 + t125;
t14 = t111 - t125;
t12 = pkin(3) * t46 + t83;
t11 = -pkin(4) * t46 - t13;
t9 = -t30 + (t34 * t71 + t36 * t73) * t72;
t8 = t112 * t46 + t83;
t7 = t55 * t115 + t47 * pkin(4) + (qJ(5) * t75 - t54 * t72) * t71 + t111;
t6 = (-pkin(3) * t109 - t123) * t71 + t88;
t4 = pkin(3) * t41 + t76;
t3 = -t41 * pkin(4) - t5;
t2 = t42 * pkin(4) + (qJD(5) * t75 - t109 * t112 - t123) * t71 + t88;
t1 = t46 * qJD(5) + t112 * t41 + t76;
t23 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t95, 0.2e1 * (-t74 ^ 2 + t75 ^ 2) * qJD(2), 0, t90, 0, 0, t74 * t100, t75 * t100, 0, 0, t21, -t77, t19, t20, t82 * t126, t57, 0.2e1 * t17 * t46 + 0.2e1 * t24 * t41 + 0.2e1 * (t109 * t15 - t75 * t9) * t71, 0.2e1 * t17 * t47 + 0.2e1 * t24 * t42 + 0.2e1 * (t10 * t75 - t109 * t16) * t71, -0.2e1 * t10 * t46 - 0.2e1 * t15 * t42 - 0.2e1 * t16 * t41 - 0.2e1 * t47 * t9, 0.2e1 * t10 * t16 + 0.2e1 * t15 * t9 + 0.2e1 * t17 * t24, t57, t81 * t126, t18, t21, -t77, t20, 0.2e1 * t13 * t41 + 0.2e1 * t14 * t42 + 0.2e1 * t46 * t5 + 0.2e1 * t47 * t6, -0.2e1 * t12 * t41 - 0.2e1 * t4 * t46 + 0.2e1 * (t109 * t14 - t6 * t75) * t71, -0.2e1 * t12 * t42 - 0.2e1 * t4 * t47 + 0.2e1 * (-t109 * t13 + t5 * t75) * t71, 0.2e1 * t12 * t4 + 0.2e1 * t13 * t5 + 0.2e1 * t14 * t6, t57, t18, t19, t20, t77, t21, -0.2e1 * t11 * t41 + 0.2e1 * t2 * t47 - 0.2e1 * t3 * t46 + 0.2e1 * t42 * t7, -0.2e1 * t1 * t47 - 0.2e1 * t8 * t42 + 0.2e1 * (t109 * t11 - t3 * t75) * t71, 0.2e1 * t1 * t46 + 0.2e1 * t8 * t41 + 0.2e1 * (-t109 * t7 + t2 * t75) * t71, 0.2e1 * t1 * t8 + 0.2e1 * t11 * t3 + 0.2e1 * t2 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, 0, -t109, 0, -pkin(7) * t108, pkin(7) * t109, 0, 0, t33, -t80, t22, -t99, t78, t58, t89 + t9 * t73 + (-pkin(2) * t41 + t109 * t48 - t17 * t72) * t71, t72 * t97 - t10 * t73 + (-pkin(2) * t42 - t109 * t49 + t17 * t70) * t71, -t49 * t41 - t48 * t42 + (t10 * t72 - t70 * t9 + (-t46 * t72 + t47 * t70) * qJD(3)) * t71, t10 * t49 + t9 * t48 + (-pkin(2) * t17 + (-t15 * t70 + t16 * t72) * qJD(3)) * t71, t58, -t22, -t78, t33, -t80, -t99, t37 * t41 + t38 * t42 - t122 + (-t5 * t72 + (qJD(3) * t47 + t6) * t70) * t71, -t89 - t39 * t41 + t6 * t73 + (t102 * t46 + t109 * t38 + t4 * t72) * t71, -t39 * t42 - t5 * t73 + (-t37 * t109 - t120 + (-t4 + t101) * t70) * t71, -t13 * t53 + t5 * t37 + t6 * t38 + t4 * t39 + (qJD(3) * t14 - qJD(4) * t12) * t119, t58, -t78, t22, -t99, t80, t33, t25 * t42 - t28 * t41 - t122 + t52 * t47 + (t2 * t70 + t3 * t72) * t71, -t29 * t42 + t3 * t73 - t45 * t47 + (-t1 * t70 + t109 * t28 - t120) * t71, -t2 * t73 + t29 * t41 + t45 * t46 + (-t1 * t72 - t109 * t25 + t52 * t75) * t71, t1 * t29 + t11 * t53 + t2 * t25 + t28 * t3 + t45 * t8 + t52 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70 * t91, t72 * t91, 0.2e1 * t107 * t72 ^ 2 + 0.2e1 * t61, 0.2e1 * (-t48 * t70 + t49 * t72) * t106, 0, 0, 0, 0, 0, 0, 0.2e1 * t117 * t53 + 0.2e1 * t61, 0.2e1 * (-t103 * t72 + t96) * t70, 0.2e1 * t103 * t68 + 0.2e1 * t121, -0.2e1 * t37 * t53 + 0.2e1 * (qJD(3) * t38 - qJD(4) * t39) * t119, 0, 0, 0, 0, 0, 0, (t52 * t70 + t53 * t72) * t124, -0.2e1 * t119 * t45 + 0.2e1 * t121, -0.2e1 * t117 * t45 - 0.2e1 * t52 * t73, 0.2e1 * t25 * t52 + 0.2e1 * t28 * t53 + 0.2e1 * t29 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t42, 0, t17, 0, 0, 0, 0, 0, 0, 0, -t41, -t42, t4, 0, 0, 0, 0, 0, 0, 0, -t42, t41, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71 * t102, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t65, 0, t6, 0, 0, 0, 0, 0, 0, t42, 0, -t65, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t65, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t23;
