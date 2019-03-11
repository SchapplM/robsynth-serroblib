% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRRP6
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
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t62 = sin(qJ(4));
t63 = sin(qJ(3));
t49 = t62 * t63;
t65 = cos(qJ(4));
t66 = cos(qJ(3));
t86 = t65 * t66;
t74 = -t49 + t86;
t87 = t65 * t63;
t39 = t62 * t66 + t87;
t61 = sin(qJ(5));
t94 = cos(qJ(5));
t22 = t94 * t39 + t61 * t74;
t107 = -0.2e1 * t22;
t53 = -t66 * pkin(3) - pkin(2);
t31 = -t74 * pkin(4) + t53;
t106 = 0.2e1 * t31;
t105 = 0.2e1 * t39;
t64 = sin(qJ(2));
t104 = -0.2e1 * t64;
t67 = cos(qJ(2));
t103 = -0.2e1 * t67;
t102 = 0.2e1 * t67;
t101 = pkin(8) + pkin(9);
t41 = -t67 * pkin(2) - t64 * pkin(8) - pkin(1);
t38 = t66 * t41;
t85 = t66 * t64;
t99 = pkin(7) * t63;
t23 = -pkin(9) * t85 + t38 + (-pkin(3) - t99) * t67;
t84 = t66 * t67;
t80 = pkin(7) * t84;
t25 = t80 + (-pkin(9) * t64 + t41) * t63;
t88 = t65 * t25;
t15 = t62 * t23 + t88;
t83 = t62 * t85 + t64 * t87;
t13 = -t83 * pkin(10) + t15;
t14 = t65 * t23 - t62 * t25;
t33 = t74 * t64;
t96 = t67 * pkin(4);
t8 = -t33 * pkin(10) + t14 - t96;
t4 = t94 * t13 + t61 * t8;
t100 = pkin(2) * t66;
t98 = t62 * pkin(3);
t97 = t67 * pkin(3);
t95 = t67 * pkin(5);
t78 = t101 * t63;
t73 = t62 * t78;
t19 = -t49 * pkin(10) + (pkin(10) + t101) * t86 - t73;
t77 = t101 * t66;
t26 = -t62 * t77 - t65 * t78;
t70 = -t39 * pkin(10) + t26;
t11 = t61 * t19 - t94 * t70;
t93 = t11 * t67;
t12 = t94 * t19 + t61 * t70;
t92 = t12 * t67;
t91 = t63 * t64;
t90 = t63 * t66;
t89 = t63 * t67;
t56 = t65 * pkin(3);
t52 = t56 + pkin(4);
t36 = t61 * t52 + t94 * t98;
t55 = t64 * pkin(7);
t40 = pkin(3) * t91 + t55;
t82 = t67 * qJ(6);
t81 = t64 * t102;
t68 = 0.2e1 * qJ(6);
t79 = t68 + t36;
t76 = t94 * pkin(4);
t75 = t61 * t13 - t94 * t8;
t72 = -t94 * t52 + t61 * t98;
t24 = t83 * pkin(4) + t40;
t71 = t76 - t72;
t69 = 0.2e1 * pkin(5);
t60 = t67 ^ 2;
t59 = t66 ^ 2;
t58 = t64 ^ 2;
t57 = t63 ^ 2;
t54 = t61 * pkin(4);
t50 = t76 + pkin(5);
t48 = t54 + qJ(6);
t34 = -pkin(5) + t72;
t32 = qJ(6) + t36;
t30 = t63 * t41 + t80;
t29 = -pkin(7) * t89 + t38;
t27 = t65 * t77 - t73;
t21 = t61 * t39 - t94 * t74;
t17 = t94 * t33 - t61 * t83;
t16 = t61 * t33 + t94 * t83;
t10 = t21 * pkin(5) - t22 * qJ(6) + t31;
t5 = t16 * pkin(5) - t17 * qJ(6) + t24;
t2 = t75 + t95;
t1 = -t82 + t4;
t3 = [1, 0, 0, t58, t81, 0, 0, 0, pkin(1) * t102, pkin(1) * t104, t59 * t58, -0.2e1 * t58 * t90, t84 * t104, t63 * t81, t60, -0.2e1 * t29 * t67 + 0.2e1 * t58 * t99, 0.2e1 * t58 * pkin(7) * t66 + 0.2e1 * t30 * t67, t33 ^ 2, -0.2e1 * t33 * t83, t33 * t103, t83 * t102, t60, -0.2e1 * t14 * t67 + 0.2e1 * t40 * t83, 0.2e1 * t15 * t67 + 0.2e1 * t40 * t33, t17 ^ 2, -0.2e1 * t17 * t16, t17 * t103, t16 * t102, t60, 0.2e1 * t24 * t16 + 0.2e1 * t67 * t75, 0.2e1 * t24 * t17 + 0.2e1 * t4 * t67, 0.2e1 * t5 * t16 + 0.2e1 * t2 * t67, -0.2e1 * t1 * t16 + 0.2e1 * t2 * t17, -0.2e1 * t1 * t67 - 0.2e1 * t5 * t17, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t64, t67, 0, -t55, -t67 * pkin(7), t63 * t85 (-t57 + t59) * t64, -t89, -t84, 0, -pkin(7) * t85 + (-pkin(2) * t64 + pkin(8) * t67) * t63, pkin(8) * t84 + (t99 - t100) * t64, t33 * t39, t33 * t74 - t39 * t83, -t39 * t67, -t74 * t67, 0, -t26 * t67 - t40 * t74 + t53 * t83, t27 * t67 + t53 * t33 + t40 * t39, t17 * t22, -t22 * t16 - t17 * t21, -t22 * t67, t21 * t67, 0, t31 * t16 + t24 * t21 + t93, t31 * t17 + t24 * t22 + t92, t10 * t16 + t5 * t21 + t93, -t1 * t21 + t11 * t17 - t12 * t16 + t2 * t22, -t10 * t17 - t5 * t22 - t92, t1 * t12 + t5 * t10 + t2 * t11; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t57, 0.2e1 * t90, 0, 0, 0, 0.2e1 * t100, -0.2e1 * pkin(2) * t63, t39 ^ 2, t74 * t105, 0, 0, 0, -0.2e1 * t53 * t74, t53 * t105, t22 ^ 2, t21 * t107, 0, 0, 0, t21 * t106, t22 * t106, 0.2e1 * t10 * t21, 0.2e1 * t11 * t22 - 0.2e1 * t12 * t21, t10 * t107, t10 ^ 2 + t11 ^ 2 + t12 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, -t91, -t67, t29, -t30, 0, 0, t33, -t83, -t67, -t65 * t97 + t14, -t88 + (-t23 + t97) * t62, 0, 0, t17, -t16, -t67, t67 * t72 - t75, t36 * t67 - t4 (-pkin(5) + t34) * t67 - t75, -t32 * t16 + t34 * t17 (-qJ(6) - t32) * t67 + t4, t1 * t32 + t2 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, t66, 0, -t63 * pkin(8), -t66 * pkin(8), 0, 0, t39, t74, 0, t26, -t27, 0, 0, t22, -t21, 0, -t11, -t12, -t11, -t32 * t21 + t34 * t22, t12, t11 * t34 + t12 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t56, -0.2e1 * t98, 0, 0, 0, 0, 1, -0.2e1 * t72, -0.2e1 * t36, -0.2e1 * t34, 0, 0.2e1 * t32, t32 ^ 2 + t34 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t83, -t67, t14, -t15, 0, 0, t17, -t16, -t67, -t67 * t76 - t75, t61 * t96 - t4 (-pkin(5) - t50) * t67 - t75, -t48 * t16 - t50 * t17 (-qJ(6) - t48) * t67 + t4, t1 * t48 - t2 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t74, 0, t26, -t27, 0, 0, t22, -t21, 0, -t11, -t12, -t11, -t48 * t21 - t50 * t22, t12, -t11 * t50 + t12 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t56, -t98, 0, 0, 0, 0, 1, t71, -t54 - t36, t69 + t71, 0, t54 + t79, t32 * t48 - t34 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t76, -0.2e1 * t54, 0.2e1 * t50, 0, 0.2e1 * t48, t48 ^ 2 + t50 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, -t67, -t75, -t4, -t75 - 0.2e1 * t95, -pkin(5) * t17 - t16 * qJ(6), -0.2e1 * t82 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, 0, -t11, -t12, -t11, -pkin(5) * t22 - t21 * qJ(6), t12, -t11 * pkin(5) + t12 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t72, -t36, t69 - t72, 0, t79, -t34 * pkin(5) + t32 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t76, -t54, t69 + t76, 0, t68 + t54, t50 * pkin(5) + t48 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t69, 0, t68, pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t17, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
