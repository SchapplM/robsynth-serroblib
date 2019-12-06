% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRRR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_inertiaDJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:44
% EndTime: 2019-12-05 15:19:49
% DurationCPUTime: 1.38s
% Computational Cost: add. (1107->163), mult. (3572->323), div. (0->0), fcn. (3796->12), ass. (0->110)
t45 = sin(qJ(4));
t124 = -0.4e1 * t45;
t41 = sin(pkin(6));
t46 = sin(qJ(3));
t115 = t41 * t46;
t43 = cos(pkin(6));
t48 = cos(qJ(4));
t25 = t45 * t115 - t48 * t43;
t100 = t25 * qJD(4);
t49 = cos(qJ(3));
t105 = qJD(3) * t49;
t86 = t41 * t105;
t123 = t48 * t86 - t100;
t44 = sin(qJ(5));
t36 = t44 ^ 2;
t47 = cos(qJ(5));
t38 = t47 ^ 2;
t112 = t36 - t38;
t79 = qJD(5) * t112;
t40 = sin(pkin(11));
t42 = sin(pkin(5));
t108 = cos(pkin(11));
t80 = t43 * t108;
t109 = cos(pkin(5));
t81 = t41 * t109;
t16 = -t49 * t81 + (t40 * t46 - t49 * t80) * t42;
t39 = t48 ^ 2;
t17 = t46 * t81 + (t40 * t49 + t46 * t80) * t42;
t55 = -t42 * t108 * t41 + t109 * t43;
t10 = t17 * t48 + t55 * t45;
t13 = t16 * qJD(3);
t3 = t10 * qJD(4) - t13 * t45;
t9 = t17 * t45 - t55 * t48;
t122 = t9 * t3;
t121 = t3 * t45;
t120 = t44 * pkin(8);
t14 = t17 * qJD(3);
t119 = t14 * t49;
t118 = t16 * t14;
t26 = t48 * t115 + t45 * t43;
t18 = t26 * qJD(4) + t45 * t86;
t117 = t18 * t45;
t116 = t25 * t18;
t114 = t41 * t49;
t113 = t47 * t48;
t37 = t45 ^ 2;
t111 = t37 - t39;
t110 = qJD(4) * t9;
t107 = qJD(3) * t46;
t106 = qJD(3) * t48;
t104 = qJD(4) * t47;
t103 = qJD(5) * t44;
t102 = qJD(5) * t47;
t101 = qJD(5) * t48;
t99 = t45 * qJD(4);
t98 = t48 * qJD(4);
t97 = t48 * t120;
t96 = pkin(8) * t113;
t95 = -0.2e1 * pkin(3) * qJD(4);
t94 = -0.2e1 * pkin(4) * qJD(5);
t92 = pkin(8) * t98;
t91 = t46 * t99;
t90 = t43 * t98;
t89 = t44 * t101;
t88 = t47 * t101;
t87 = t41 * t107;
t85 = t44 * t98;
t84 = t44 * t102;
t83 = t45 * t98;
t82 = t47 * t98;
t78 = t111 * qJD(4);
t77 = 0.2e1 * t83;
t75 = t44 * t82;
t74 = t37 * t84;
t72 = -t48 * pkin(4) - t45 * pkin(9);
t71 = pkin(4) * t45 - pkin(9) * t48;
t70 = t9 * t18 + t3 * t25;
t5 = -t10 * t44 + t16 * t47;
t6 = t10 * t47 + t16 * t44;
t69 = -t44 * t6 - t47 * t5;
t20 = -t44 * t114 + t47 * t26;
t62 = t47 * t114 + t44 * t26;
t68 = -t20 * t44 + t47 * t62;
t65 = pkin(3) - t72;
t58 = t47 * t65;
t22 = -t58 - t97;
t23 = -t44 * t65 + t96;
t67 = -t22 * t47 - t23 * t44;
t64 = t9 * t102 + t3 * t44;
t63 = t9 * t103 - t3 * t47;
t61 = t71 * t44;
t60 = t25 * t102 + t18 * t44;
t59 = t25 * t103 - t18 * t47;
t57 = t48 * t105 - t91;
t56 = t47 * t99 + t89;
t4 = -t13 * t48 - t110;
t1 = -t6 * qJD(5) + t14 * t47 - t4 * t44;
t2 = t5 * qJD(5) + t14 * t44 + t4 * t47;
t54 = t69 * qJD(5) - t1 * t44 + t2 * t47;
t53 = t121 + t4 * t48 + (-t10 * t45 + t48 * t9) * qJD(4);
t7 = t62 * qJD(5) - t47 * t123 - t44 * t87;
t8 = -t43 * t85 - t26 * t102 + (t47 * t107 + (t91 + (qJD(5) - t106) * t49) * t44) * t41;
t52 = t68 * qJD(5) - t8 * t44 - t7 * t47;
t11 = t56 * pkin(8) - qJD(4) * t61 + qJD(5) * t58;
t12 = -t23 * qJD(5) + (t45 * t120 + t47 * t71) * qJD(4);
t51 = t67 * qJD(5) - t11 * t47 - t12 * t44;
t50 = -t26 * t99 + t39 * t86 + t117;
t34 = -0.2e1 * t83;
t21 = t45 * t79 - t75;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t17 * t13 + 0.2e1 * t118, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t10 * t4 + 0.2e1 * t118 + 0.2e1 * t122, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t5 * t1 + 0.2e1 * t6 * t2 + 0.2e1 * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t13 * t46 - t119 + (t16 * t46 + t17 * t49) * qJD(3)) * t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t90 + t4 * t26 + (t10 * t57 + t16 * t107 - t119) * t41 + t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t62 + t2 * t20 + t5 * t8 - t6 * t7 + t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t26 * t90 + 0.2e1 * t116 + 0.2e1 * (t26 * t57 - t46 * t86) * t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t20 * t7 - 0.2e1 * t62 * t8 + 0.2e1 * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t13, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * t48 + t16 * t99, t14 * t45 + t16 * t98, t53, -t14 * pkin(3) + t53 * pkin(8), 0, 0, 0, 0, 0, 0, (t44 * t110 - t1) * t48 + (qJD(4) * t5 + t64) * t45, (t9 * t104 + t2) * t48 + (-qJD(4) * t6 - t63) * t45, t69 * t98 + (-t1 * t47 - t2 * t44 + (t44 * t5 - t47 * t6) * qJD(5)) * t45, t1 * t22 - t6 * t11 + t5 * t12 + t2 * t23 + (t9 * t98 + t121) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t86, 0, 0, 0, 0, 0, 0, 0, 0, (-t46 * t106 - t49 * t99) * t41, (t45 * t107 - t49 * t98) * t41, t50, -pkin(3) * t87 + t50 * pkin(8), 0, 0, 0, 0, 0, 0, (t44 * t100 - t8) * t48 + (-qJD(4) * t62 + t60) * t45, (t47 * t100 - t7) * t48 + (-qJD(4) * t20 - t59) * t45, t68 * t98 + (t44 * t7 - t47 * t8 + (-t20 * t47 - t44 * t62) * qJD(5)) * t45, -t20 * t11 - t62 * t12 + t8 * t22 - t7 * t23 + (t25 * t98 + t117) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -0.2e1 * t78, 0, t34, 0, 0, t45 * t95, t48 * t95, 0, 0, 0.2e1 * t38 * t83 - 0.2e1 * t74, t75 * t124 + 0.2e1 * t37 * t79, 0.2e1 * t111 * t104 + 0.2e1 * t45 * t89, 0.2e1 * t36 * t83 + 0.2e1 * t74, -0.2e1 * t44 * t78 + 0.2e1 * t45 * t88, t34, 0.2e1 * t22 * t99 - 0.2e1 * t12 * t48 + 0.2e1 * (t37 * t102 + t44 * t77) * pkin(8), -0.2e1 * t23 * t99 - 0.2e1 * t11 * t48 + 0.2e1 * (-t37 * t103 + t47 * t77) * pkin(8), 0.2e1 * t67 * t98 + 0.2e1 * (t11 * t44 - t12 * t47 + (t22 * t44 - t23 * t47) * qJD(5)) * t45, 0.2e1 * pkin(8) ^ 2 * t83 - 0.2e1 * t23 * t11 + 0.2e1 * t22 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, 0, 0, 0, 0, 0, 0, 0, t63, t64, t54, -t3 * pkin(4) + pkin(9) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t123, 0, 0, 0, 0, 0, 0, 0, 0, t59, t60, t52, -t18 * pkin(4) + pkin(9) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, 0, -t99, 0, -t92, pkin(8) * t99, 0, 0, -t21, -t112 * t98 + t84 * t124, t44 * t99 - t88, t21, t56, 0, (pkin(9) * t113 + (-t47 * pkin(4) + t120) * t45) * qJD(5) + (t72 * t44 - t96) * qJD(4), (pkin(8) * t45 * t47 + t61) * qJD(5) + (t72 * t47 + t97) * qJD(4), t51, -pkin(4) * t92 + pkin(9) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t84, -0.2e1 * t79, 0, -0.2e1 * t84, 0, 0, t44 * t94, t47 * t94, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45 * t103 + t82, 0, -t45 * t102 - t85, t99, t12, t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, 0, -t103, 0, -pkin(9) * t102, pkin(9) * t103, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t15;
