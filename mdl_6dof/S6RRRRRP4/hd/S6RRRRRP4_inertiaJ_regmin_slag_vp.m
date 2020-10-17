% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:51:40
% EndTime: 2019-05-08 04:51:45
% DurationCPUTime: 1.14s
% Computational Cost: add. (1717->159), mult. (3222->249), div. (0->0), fcn. (3778->8), ass. (0->111)
t81 = sin(qJ(5));
t82 = sin(qJ(4));
t106 = t81 * t82;
t112 = cos(qJ(5));
t85 = cos(qJ(4));
t90 = t112 * t85;
t57 = -t90 + t106;
t128 = 0.2e1 * t57;
t113 = cos(qJ(3));
t114 = cos(qJ(2));
t83 = sin(qJ(3));
t84 = sin(qJ(2));
t58 = -t113 * t114 + t83 * t84;
t127 = -0.2e1 * t58;
t126 = 0.2e1 * t58;
t91 = t112 * t82;
t59 = t81 * t85 + t91;
t125 = -0.2e1 * t59;
t118 = t85 * pkin(4);
t93 = t113 * pkin(2);
t74 = -t93 - pkin(3);
t63 = t74 - t118;
t124 = 0.2e1 * t63;
t75 = -pkin(3) - t118;
t123 = 0.2e1 * t75;
t76 = -pkin(2) * t114 - pkin(1);
t122 = 0.2e1 * t76;
t121 = -pkin(10) - pkin(9);
t51 = t58 * qJ(6);
t60 = t113 * t84 + t114 * t83;
t102 = t85 * t60;
t34 = pkin(3) * t58 - pkin(9) * t60 + t76;
t64 = (-pkin(8) - pkin(7)) * t84;
t94 = t114 * pkin(7);
t66 = pkin(8) * t114 + t94;
t44 = t113 * t66 + t64 * t83;
t16 = t34 * t85 - t44 * t82;
t13 = pkin(4) * t58 - pkin(10) * t102 + t16;
t103 = t85 * t44;
t15 = t103 + (-pkin(10) * t60 + t34) * t82;
t6 = t112 * t15 + t13 * t81;
t3 = t51 + t6;
t5 = t112 * t13 - t15 * t81;
t52 = t58 * pkin(5);
t4 = -t5 - t52;
t120 = -t3 * t57 + t4 * t59;
t77 = t81 * pkin(4);
t119 = t83 * pkin(2);
t117 = t85 * pkin(9);
t116 = pkin(3) - t74;
t71 = pkin(9) + t119;
t115 = -pkin(10) - t71;
t101 = t85 * t71;
t78 = t85 * pkin(10);
t53 = t78 + t101;
t31 = -t115 * t91 + t53 * t81;
t111 = t31 * t58;
t32 = t106 * t115 + t112 * t53;
t110 = t32 * t58;
t65 = t78 + t117;
t41 = -t121 * t91 + t65 * t81;
t109 = t41 * t58;
t42 = -t113 * t64 + t66 * t83;
t108 = t42 * t85;
t43 = t106 * t121 + t112 * t65;
t107 = t43 * t58;
t105 = t82 * t60;
t104 = t82 * t85;
t100 = t31 * t59 - t32 * t57;
t99 = t41 * t59 - t43 * t57;
t33 = pkin(5) * t57 - qJ(6) * t59 + t75;
t27 = -t93 + t33;
t98 = t27 + t33;
t97 = t63 + t75;
t96 = 0.2e1 * t114;
t95 = t60 * t127;
t92 = t112 * pkin(4);
t23 = pkin(4) * t105 + t42;
t89 = -pkin(3) * t60 - pkin(9) * t58;
t88 = -t58 * t71 + t60 * t74;
t87 = 0.2e1 * pkin(5);
t86 = 0.2e1 * qJ(6);
t80 = t85 ^ 2;
t79 = t82 ^ 2;
t72 = t92 + pkin(5);
t68 = t77 + qJ(6);
t67 = 0.2e1 * t104;
t56 = t60 ^ 2;
t55 = t59 ^ 2;
t54 = t58 ^ 2;
t50 = t85 * t58;
t49 = t82 * t58;
t46 = t82 * t102;
t40 = t59 * t58;
t39 = t57 * t58;
t38 = t57 * t125;
t37 = t42 * t82;
t36 = -pkin(5) * t59 - qJ(6) * t57;
t35 = (-t79 + t80) * t60;
t29 = -t105 * t81 + t60 * t90;
t28 = t59 * t60;
t24 = -t57 * t68 - t59 * t72;
t20 = t29 * t59;
t19 = t23 * t59;
t18 = t23 * t57;
t17 = t34 * t82 + t103;
t10 = -t28 * t59 - t29 * t57;
t9 = pkin(5) * t28 - qJ(6) * t29 + t23;
t8 = t9 * t59;
t7 = t9 * t57;
t1 = [1, 0, 0, t84 ^ 2, t84 * t96, 0, 0, 0, pkin(1) * t96, -0.2e1 * pkin(1) * t84, t56, t95, 0, 0, 0, t58 * t122, t60 * t122, t80 * t56, -0.2e1 * t56 * t104, t102 * t126, t82 * t95, t54, 0.2e1 * t105 * t42 + 0.2e1 * t16 * t58, 0.2e1 * t102 * t42 - 0.2e1 * t17 * t58, t29 ^ 2, -0.2e1 * t29 * t28, t29 * t126, t28 * t127, t54, 0.2e1 * t23 * t28 + 0.2e1 * t5 * t58, 0.2e1 * t23 * t29 - 0.2e1 * t58 * t6, 0.2e1 * t28 * t9 - 0.2e1 * t4 * t58, -0.2e1 * t28 * t3 + 0.2e1 * t29 * t4, -0.2e1 * t29 * t9 + 0.2e1 * t3 * t58, t3 ^ 2 + t4 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, t84, t114, 0, -t84 * pkin(7), -t94, 0, 0, t60, -t58, 0, -t42, -t44, t46, t35, t49, t50, 0, t82 * t88 - t108, t85 * t88 + t37, t20, t10, t40, -t39, 0, t28 * t63 - t111 + t18, t29 * t63 - t110 + t19, t27 * t28 - t111 + t7, -t28 * t32 + t29 * t31 + t120, -t27 * t29 + t110 - t8, t27 * t9 + t3 * t32 + t31 * t4; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t93, -0.2e1 * t119, t79, t67, 0, 0, 0, -0.2e1 * t74 * t85, 0.2e1 * t74 * t82, t55, t38, 0, 0, 0, t57 * t124, t59 * t124, t27 * t128, 0.2e1 * t100, t27 * t125, t27 ^ 2 + t31 ^ 2 + t32 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t58, 0, -t42, -t44, t46, t35, t49, t50, 0, t82 * t89 - t108, t85 * t89 + t37, t20, t10, t40, -t39, 0, t28 * t75 - t109 + t18, t29 * t75 - t107 + t19, t28 * t33 - t109 + t7, -t28 * t43 + t29 * t41 + t120, -t29 * t33 + t107 - t8, t3 * t43 + t33 * t9 + t4 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t93, -t119, t79, t67, 0, 0, 0, t116 * t85, -t116 * t82, t55, t38, 0, 0, 0, t97 * t57, t97 * t59, t98 * t57, t99 + t100, -t98 * t59, t27 * t33 + t31 * t41 + t32 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t79, t67, 0, 0, 0, 0.2e1 * pkin(3) * t85, -0.2e1 * pkin(3) * t82, t55, t38, 0, 0, 0, t57 * t123, t59 * t123, t33 * t128, 0.2e1 * t99, t33 * t125, t33 ^ 2 + t41 ^ 2 + t43 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, -t105, t58, t16, -t17, 0, 0, t29, -t28, t58, t58 * t92 + t5, -t58 * t77 - t6, t58 * t72 - t4, -t28 * t68 - t29 * t72, t58 * t68 + t3, t3 * t68 - t4 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t85, 0, -t82 * t71, -t101, 0, 0, t59, -t57, 0, -t31, -t32, -t31, t24, t32, -t31 * t72 + t32 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t85, 0, -t82 * pkin(9), -t117, 0, 0, t59, -t57, 0, -t41, -t43, -t41, t24, t43, -t41 * t72 + t43 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t92, -0.2e1 * t77, 0.2e1 * t72, 0, 0.2e1 * t68, t68 ^ 2 + t72 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, t58, t5, -t6, -t4 + t52, -pkin(5) * t29 - qJ(6) * t28, 0.2e1 * t51 + t6, -pkin(5) * t4 + qJ(6) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t57, 0, -t31, -t32, -t31, t36, t32, -pkin(5) * t31 + qJ(6) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t57, 0, -t41, -t43, -t41, t36, t43, -pkin(5) * t41 + qJ(6) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t92, -t77, t87 + t92, 0, t86 + t77, pkin(5) * t72 + qJ(6) * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t87, 0, t86, pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t29, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t1;
