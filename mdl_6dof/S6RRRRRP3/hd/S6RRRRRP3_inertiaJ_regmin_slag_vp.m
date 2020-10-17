% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRRP3
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
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:42:47
% EndTime: 2019-05-08 04:42:50
% DurationCPUTime: 0.86s
% Computational Cost: add. (1302->129), mult. (2496->218), div. (0->0), fcn. (2948->8), ass. (0->103)
t103 = cos(qJ(3));
t104 = cos(qJ(2));
t82 = sin(qJ(3));
t83 = sin(qJ(2));
t58 = -t103 * t104 + t82 * t83;
t117 = -0.2e1 * t58;
t116 = 0.2e1 * t58;
t85 = cos(qJ(4));
t107 = t85 * pkin(4);
t88 = t103 * pkin(2);
t73 = -t88 - pkin(3);
t63 = t73 - t107;
t115 = 0.2e1 * t63;
t74 = -pkin(3) - t107;
t114 = 0.2e1 * t74;
t75 = -t104 * pkin(2) - pkin(1);
t113 = 0.2e1 * t75;
t60 = t103 * t83 + t82 * t104;
t81 = sin(qJ(4));
t101 = t81 * t60;
t80 = sin(qJ(5));
t84 = cos(qJ(5));
t97 = t85 * t60;
t27 = -t80 * t101 + t84 * t97;
t32 = t58 * pkin(3) - t60 * pkin(9) + t75;
t65 = (-pkin(8) - pkin(7)) * t83;
t89 = t104 * pkin(7);
t67 = t104 * pkin(8) + t89;
t41 = t103 * t67 + t82 * t65;
t98 = t85 * t41;
t10 = t98 + (-pkin(10) * t60 + t32) * t81;
t110 = t58 * pkin(4);
t12 = t85 * t32 - t81 * t41;
t9 = -pkin(10) * t97 + t110 + t12;
t5 = -t80 * t10 + t84 * t9;
t2 = t58 * pkin(5) - t27 * qJ(6) + t5;
t59 = t80 * t85 + t84 * t81;
t26 = t59 * t60;
t99 = t84 * t10;
t6 = t80 * t9 + t99;
t4 = -t26 * qJ(6) + t6;
t57 = t80 * t81 - t84 * t85;
t112 = -t2 * t59 - t4 * t57;
t111 = t57 * pkin(5);
t109 = t80 * pkin(4);
t108 = t82 * pkin(2);
t76 = t84 * pkin(4);
t106 = t85 * pkin(9);
t105 = pkin(3) - t73;
t39 = -t103 * t65 + t82 * t67;
t102 = t39 * t85;
t100 = t81 * t85;
t71 = pkin(9) + t108;
t96 = t85 * t71;
t52 = (-pkin(10) - t71) * t81;
t77 = t85 * pkin(10);
t53 = t77 + t96;
t30 = t84 * t52 - t80 * t53;
t92 = t59 * qJ(6);
t21 = t30 - t92;
t31 = t80 * t52 + t84 * t53;
t50 = t57 * qJ(6);
t22 = t31 - t50;
t95 = -t21 * t59 - t22 * t57;
t64 = (-pkin(9) - pkin(10)) * t81;
t66 = t77 + t106;
t38 = t84 * t64 - t80 * t66;
t24 = t38 - t92;
t40 = t80 * t64 + t84 * t66;
t25 = t40 - t50;
t94 = -t24 * t59 - t25 * t57;
t93 = t63 + t74;
t91 = 0.2e1 * t104;
t90 = t60 * t117;
t23 = pkin(4) * t101 + t39;
t87 = -pkin(3) * t60 - pkin(9) * t58;
t86 = -t58 * t71 + t60 * t73;
t79 = t85 ^ 2;
t78 = t81 ^ 2;
t72 = t76 + pkin(5);
t68 = 0.2e1 * t100;
t56 = t60 ^ 2;
t55 = t59 ^ 2;
t54 = t58 ^ 2;
t51 = pkin(5) * t59;
t49 = t85 * t58;
t48 = t81 * t58;
t45 = t81 * t97;
t43 = t74 + t111;
t42 = t63 + t111;
t37 = t59 * t58;
t36 = t57 * t58;
t35 = -0.2e1 * t59 * t57;
t34 = t39 * t81;
t33 = (-t78 + t79) * t60;
t29 = -t57 * t109 - t72 * t59;
t20 = t27 * t59;
t17 = t23 * t59;
t16 = t23 * t57;
t13 = t81 * t32 + t98;
t11 = t26 * pkin(5) + t23;
t7 = -t59 * t26 - t27 * t57;
t1 = [1, 0, 0, t83 ^ 2, t83 * t91, 0, 0, 0, pkin(1) * t91, -0.2e1 * pkin(1) * t83, t56, t90, 0, 0, 0, t58 * t113, t60 * t113, t79 * t56, -0.2e1 * t56 * t100, t97 * t116, t81 * t90, t54, 0.2e1 * t39 * t101 + 0.2e1 * t12 * t58, -0.2e1 * t13 * t58 + 0.2e1 * t39 * t97, t27 ^ 2, -0.2e1 * t27 * t26, t27 * t116, t26 * t117, t54, 0.2e1 * t23 * t26 + 0.2e1 * t5 * t58, 0.2e1 * t23 * t27 - 0.2e1 * t6 * t58, -0.2e1 * t2 * t27 - 0.2e1 * t4 * t26, t11 ^ 2 + t2 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, t83, t104, 0, -t83 * pkin(7), -t89, 0, 0, t60, -t58, 0, -t39, -t41, t45, t33, t48, t49, 0, t81 * t86 - t102, t85 * t86 + t34, t20, t7, t37, -t36, 0, t63 * t26 + t30 * t58 + t16, t63 * t27 - t31 * t58 + t17, -t21 * t27 - t22 * t26 + t112, t11 * t42 + t2 * t21 + t4 * t22; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t88, -0.2e1 * t108, t78, t68, 0, 0, 0, -0.2e1 * t73 * t85, 0.2e1 * t73 * t81, t55, t35, 0, 0, 0, t57 * t115, t59 * t115, 0.2e1 * t95, t21 ^ 2 + t22 ^ 2 + t42 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t58, 0, -t39, -t41, t45, t33, t48, t49, 0, t81 * t87 - t102, t85 * t87 + t34, t20, t7, t37, -t36, 0, t74 * t26 + t38 * t58 + t16, t74 * t27 - t40 * t58 + t17, -t24 * t27 - t25 * t26 + t112, t11 * t43 + t2 * t24 + t4 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t88, -t108, t78, t68, 0, 0, 0, t105 * t85, -t105 * t81, t55, t35, 0, 0, 0, t93 * t57, t93 * t59, t94 + t95, t21 * t24 + t22 * t25 + t42 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t78, t68, 0, 0, 0, 0.2e1 * pkin(3) * t85, -0.2e1 * pkin(3) * t81, t55, t35, 0, 0, 0, t57 * t114, t59 * t114, 0.2e1 * t94, t24 ^ 2 + t25 ^ 2 + t43 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, -t101, t58, t12, -t13, 0, 0, t27, -t26, t58, t58 * t76 + t5, -t99 + (-t9 - t110) * t80, -t26 * t109 - t72 * t27, t4 * t109 + t2 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t85, 0, -t81 * t71, -t96, 0, 0, t59, -t57, 0, t30, -t31, t29, t22 * t109 + t21 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t85, 0, -t81 * pkin(9), -t106, 0, 0, t59, -t57, 0, t38, -t40, t29, t25 * t109 + t24 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t76, -0.2e1 * t109, 0, t80 ^ 2 * pkin(4) ^ 2 + t72 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, t58, t5, -t6, -pkin(5) * t27, t2 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t57, 0, t30, -t31, -t51, t21 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t57, 0, t38, -t40, -t51, t24 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t76, -t109, 0, t72 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t1;
