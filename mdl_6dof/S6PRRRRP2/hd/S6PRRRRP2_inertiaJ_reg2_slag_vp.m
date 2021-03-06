% Calculate inertial parameters regressor of joint inertia matrix for
% S6PRRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRRP2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:37:55
% EndTime: 2019-05-05 09:38:02
% DurationCPUTime: 2.15s
% Computational Cost: add. (1151->158), mult. (2468->254), div. (0->0), fcn. (2906->10), ass. (0->111)
t81 = sin(qJ(5));
t75 = t81 ^ 2;
t85 = cos(qJ(5));
t77 = t85 ^ 2;
t145 = t75 + t77;
t82 = sin(qJ(4));
t134 = t82 * pkin(3);
t67 = pkin(10) + t134;
t114 = t145 * t67;
t79 = sin(pkin(6));
t88 = cos(qJ(2));
t122 = t79 * t88;
t84 = sin(qJ(2));
t123 = t79 * t84;
t80 = cos(pkin(6));
t83 = sin(qJ(3));
t87 = cos(qJ(3));
t37 = -t83 * t123 + t80 * t87;
t38 = t87 * t123 + t80 * t83;
t86 = cos(qJ(4));
t22 = t82 * t37 + t86 * t38;
t15 = t85 * t122 + t81 * t22;
t17 = -t81 * t122 + t85 * t22;
t106 = t15 * t81 + t17 * t85;
t20 = -t86 * t37 + t82 * t38;
t19 = t20 ^ 2;
t138 = -pkin(9) - pkin(8);
t107 = t138 * t83;
t57 = t138 * t87;
t31 = -t86 * t107 - t82 * t57;
t144 = t31 ^ 2;
t48 = t82 * t83 - t86 * t87;
t46 = t48 ^ 2;
t50 = t82 * t87 + t86 * t83;
t143 = 0.2e1 * t50;
t69 = -t87 * pkin(3) - pkin(2);
t142 = 0.2e1 * t69;
t141 = -0.2e1 * t81;
t140 = -0.2e1 * t85;
t139 = 0.2e1 * t87;
t137 = pkin(10) * t48;
t136 = t48 * pkin(5);
t135 = t81 * pkin(10);
t133 = t85 * pkin(10);
t132 = t86 * pkin(3);
t68 = -pkin(4) - t132;
t131 = pkin(4) - t68;
t97 = pkin(5) * t81 - t85 * qJ(6);
t10 = t97 * t50 + t31;
t130 = t10 * t81;
t129 = t10 * t85;
t127 = t20 * t31;
t126 = t20 * t85;
t125 = t31 * t85;
t124 = t48 * t67;
t121 = t81 * t50;
t120 = t81 * t67;
t119 = t81 * t85;
t44 = t85 * t50;
t118 = t85 * t67;
t25 = t48 * pkin(4) - t50 * pkin(10) + t69;
t33 = t82 * t107 - t86 * t57;
t8 = t81 * t25 + t85 * t33;
t98 = -t85 * pkin(5) - t81 * qJ(6);
t53 = -pkin(4) + t98;
t45 = t53 - t132;
t117 = -t45 - t53;
t116 = t114 * pkin(10);
t115 = t145 * t67 ^ 2;
t113 = t145 * pkin(10) ^ 2;
t112 = t145 * pkin(10);
t76 = t83 ^ 2;
t78 = t87 ^ 2;
t111 = t76 + t78;
t110 = t48 * qJ(6);
t109 = t48 * t121;
t47 = t50 ^ 2;
t108 = t47 * t119;
t105 = t20 * t121 - t15 * t48;
t104 = -t85 * t25 + t81 * t33;
t103 = t15 ^ 2 + t17 ^ 2 + t19;
t102 = t106 * pkin(10);
t101 = t17 * t118 + t15 * t120;
t100 = -pkin(4) * t50 - t137;
t5 = t110 + t8;
t6 = t104 - t136;
t1 = t5 * t85 + t6 * t81;
t2 = t104 * t81 + t8 * t85;
t99 = -t50 * t53 + t137;
t96 = -t37 * t83 + t38 * t87;
t95 = t45 * t50 - t124;
t94 = t50 * t68 - t124;
t93 = -t17 * t48 + t20 * t44;
t92 = (t15 * t85 - t17 * t81) * t50;
t74 = t79 ^ 2;
t65 = t74 * t88 ^ 2;
t63 = -0.2e1 * t119;
t62 = 0.2e1 * t119;
t52 = 0.2e1 * t112;
t43 = t85 * t48;
t41 = t77 * t47;
t40 = t81 * t48;
t39 = t75 * t47;
t36 = t81 * t44;
t34 = 0.2e1 * t114;
t30 = t112 + t114;
t29 = t31 * t81;
t27 = 0.2e1 * t48 * t44;
t26 = (t75 - t77) * t50;
t18 = t20 * t81;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t84 ^ 2 + t80 ^ 2 + t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 ^ 2 + t38 ^ 2 + t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 ^ 2 + t19 + t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, -t123, 0, 0, 0, 0, 0, 0, 0, 0, t87 * t122, -t83 * t122, t96, pkin(2) * t122 + t96 * pkin(8), 0, 0, 0, 0, 0, 0, -t48 * t122, -t50 * t122, t20 * t50 - t22 * t48, -t69 * t122 + t22 * t33 + t127, 0, 0, 0, 0, 0, 0, t105, t93, t92, t104 * t15 + t17 * t8 + t127, 0, 0, 0, 0, 0, 0, t105, t92, -t93, t20 * t10 + t15 * t6 + t17 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t76, t83 * t139, 0, t78, 0, 0, pkin(2) * t139, -0.2e1 * pkin(2) * t83, 0.2e1 * t111 * pkin(8), t111 * pkin(8) ^ 2 + pkin(2) ^ 2, t47, -0.2e1 * t50 * t48, 0, t46, 0, 0, t48 * t142, t50 * t142, 0.2e1 * t31 * t50 - 0.2e1 * t33 * t48, t33 ^ 2 + t69 ^ 2 + t144, t41, -0.2e1 * t108, t27, t39, -0.2e1 * t109, t46, -0.2e1 * t104 * t48 + 0.2e1 * t31 * t121, 0.2e1 * t31 * t44 - 0.2e1 * t8 * t48 (t104 * t85 - t8 * t81) * t143, t104 ^ 2 + t8 ^ 2 + t144, t41, t27, 0.2e1 * t108, t46, 0.2e1 * t109, t39, 0.2e1 * t10 * t121 - 0.2e1 * t6 * t48 (-t5 * t81 + t6 * t85) * t143, -0.2e1 * t10 * t44 + 0.2e1 * t5 * t48, t10 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t38, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t22, 0 (-t20 * t86 + t22 * t82) * pkin(3), 0, 0, 0, 0, 0, 0, -t126, t18, t106, t20 * t68 + t101, 0, 0, 0, 0, 0, 0, -t126, t106, -t18, t20 * t45 + t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, 0, t87, 0, -t83 * pkin(8), -t87 * pkin(8), 0, 0, 0, 0, t50, 0, -t48, 0, -t31, -t33 (-t48 * t82 - t50 * t86) * pkin(3) (-t31 * t86 + t33 * t82) * pkin(3), t36, -t26, t40, -t36, t43, 0, t94 * t81 - t125, t94 * t85 + t29, t2, t2 * t67 + t31 * t68, t36, t40, t26, 0, -t43, -t36, t95 * t81 - t129, t1, -t95 * t85 - t130, t1 * t67 + t10 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t132, -0.2e1 * t134, 0 (t82 ^ 2 + t86 ^ 2) * pkin(3) ^ 2, t75, t62, 0, t77, 0, 0, t68 * t140, 0.2e1 * t68 * t81, t34, t68 ^ 2 + t115, t75, 0, t63, 0, 0, t77, t45 * t140, t34, t45 * t141, t45 ^ 2 + t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t22, 0, 0, 0, 0, 0, 0, 0, 0, -t126, t18, t106, -t20 * pkin(4) + t102, 0, 0, 0, 0, 0, 0, -t126, t106, -t18, t20 * t53 + t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t48, 0, -t31, -t33, 0, 0, t36, -t26, t40, -t36, t43, 0, t100 * t81 - t125, t100 * t85 + t29, t2, -t31 * pkin(4) + t2 * pkin(10), t36, t40, t26, 0, -t43, -t36, -t99 * t81 - t129, t1, t99 * t85 - t130, t1 * pkin(10) + t10 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t132, -t134, 0, 0, t75, t62, 0, t77, 0, 0, t131 * t85, -t131 * t81, t30, -t68 * pkin(4) + t116, t75, 0, t63, 0, 0, t77, t117 * t85, t30, t117 * t81, t45 * t53 + t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t75, t62, 0, t77, 0, 0, 0.2e1 * pkin(4) * t85, pkin(4) * t141, t52, pkin(4) ^ 2 + t113, t75, 0, t63, 0, 0, t77, t53 * t140, t52, t53 * t141, t53 ^ 2 + t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t17, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, t17, -t15 * pkin(5) + t17 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, -t121, t48, -t104, -t8, 0, 0, 0, t44, 0, t48, t121, 0, -t104 + 0.2e1 * t136, t98 * t50, 0.2e1 * t110 + t8, -t6 * pkin(5) + t5 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0, t85, 0, -t120, -t118, 0, 0, 0, t81, 0, 0, -t85, 0, -t120, -t97, t118, -t97 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0, t85, 0, -t135, -t133, 0, 0, 0, t81, 0, 0, -t85, 0, -t135, -t97, t133, -t97 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t44, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0, t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0, t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
