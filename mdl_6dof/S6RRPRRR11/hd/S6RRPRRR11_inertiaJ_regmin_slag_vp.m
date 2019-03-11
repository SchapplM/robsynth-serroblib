% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRR11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t66 = sin(qJ(2));
t64 = sin(qJ(5));
t65 = sin(qJ(4));
t68 = cos(qJ(5));
t69 = cos(qJ(4));
t38 = t64 * t69 + t68 * t65;
t39 = -t64 * t65 + t68 * t69;
t63 = sin(qJ(6));
t67 = cos(qJ(6));
t74 = t67 * t38 + t63 * t39;
t103 = t74 * t66;
t76 = -t63 * t38 + t67 * t39;
t102 = t76 * t66;
t49 = t65 * pkin(4) + qJ(3);
t27 = t38 * pkin(5) + t49;
t101 = 0.2e1 * t27;
t100 = 0.2e1 * t49;
t99 = -0.2e1 * t66;
t98 = 0.2e1 * t66;
t70 = cos(qJ(2));
t97 = 0.2e1 * t70;
t96 = 0.2e1 * qJ(3);
t71 = -pkin(2) - pkin(8);
t95 = t63 * pkin(5);
t94 = t64 * pkin(4);
t93 = t66 * pkin(4);
t92 = t66 * pkin(5);
t55 = t67 * pkin(5);
t83 = t69 * t70;
t87 = t65 * t70;
t28 = t64 * t87 - t68 * t83;
t54 = t66 * pkin(7);
t44 = t66 * pkin(3) + t54;
t40 = t69 * t44;
t77 = -t66 * qJ(3) - pkin(1);
t36 = t71 * t70 + t77;
t78 = pkin(9) * t70 - t36;
t15 = t78 * t65 + t40 + t93;
t89 = t65 * t44;
t16 = -t78 * t69 + t89;
t85 = t68 * t16;
t9 = t64 * t15 + t85;
t5 = t28 * pkin(10) + t9;
t91 = t67 * t5;
t56 = t68 * pkin(4);
t90 = t38 * t66;
t88 = t65 * t66;
t86 = t66 * t70;
t84 = t69 * t65;
t57 = t70 * pkin(7);
t45 = t70 * pkin(3) + t57;
t60 = t66 ^ 2;
t62 = t70 ^ 2;
t82 = t60 + t62;
t81 = t70 * qJ(3);
t80 = -0.2e1 * t86;
t79 = t67 * t94;
t32 = pkin(4) * t83 + t45;
t29 = t38 * t70;
t8 = t68 * t15 - t64 * t16;
t4 = t29 * pkin(10) + t8 + t92;
t1 = t67 * t4 - t63 * t5;
t41 = (-pkin(9) + t71) * t65;
t52 = t69 * t71;
t42 = -t69 * pkin(9) + t52;
t24 = -t64 * t41 + t68 * t42;
t51 = t56 + pkin(5);
t30 = t67 * t51 - t63 * t94;
t2 = t63 * t4 + t91;
t75 = -t66 * pkin(2) + t81;
t25 = t68 * t41 + t64 * t42;
t73 = t66 * t71 + t81;
t61 = t69 ^ 2;
t59 = t65 ^ 2;
t50 = t69 * t66;
t43 = -t70 * pkin(2) + t77;
t35 = t39 * t66;
t31 = t63 * t51 + t79;
t23 = t69 * t36 + t89;
t22 = -t65 * t36 + t40;
t17 = -t28 * pkin(5) + t32;
t14 = -t38 * pkin(10) + t25;
t13 = -t39 * pkin(10) + t24;
t12 = t63 * t28 - t67 * t29;
t11 = -t67 * t28 - t63 * t29;
t7 = t63 * t13 + t67 * t14;
t6 = t67 * t13 - t63 * t14;
t3 = [1, 0, 0, t60, 0.2e1 * t86, 0, 0, 0, pkin(1) * t97, pkin(1) * t99, 0.2e1 * t82 * pkin(7), t43 * t97, t43 * t99, t82 * pkin(7) ^ 2 + t43 ^ 2, t59 * t62, 0.2e1 * t62 * t84, t65 * t80, t69 * t80, t60, 0.2e1 * t22 * t66 + 0.2e1 * t45 * t83, -0.2e1 * t23 * t66 - 0.2e1 * t45 * t87, t29 ^ 2, -0.2e1 * t29 * t28, -t29 * t98, t28 * t98, t60, -0.2e1 * t32 * t28 + 0.2e1 * t8 * t66, -0.2e1 * t32 * t29 - 0.2e1 * t9 * t66, t12 ^ 2, -0.2e1 * t12 * t11, t12 * t98, t11 * t99, t60, 0.2e1 * t1 * t66 + 0.2e1 * t17 * t11, 0.2e1 * t17 * t12 - 0.2e1 * t2 * t66; 0, 0, 0, 0, 0, t66, t70, 0, -t54, -t57, t75, t54, t57, t75 * pkin(7), -t65 * t83 (t59 - t61) * t70, t50, -t88, 0, t45 * t65 + t73 * t69, t45 * t69 - t73 * t65, -t29 * t39, t39 * t28 + t29 * t38, t35, -t90, 0, t24 * t66 - t49 * t28 + t32 * t38, -t25 * t66 - t49 * t29 + t32 * t39, t12 * t76, -t11 * t76 - t12 * t74, t102, -t103, 0, t27 * t11 + t17 * t74 + t6 * t66, t27 * t12 + t17 * t76 - t7 * t66; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(2), t96, pkin(2) ^ 2 + qJ(3) ^ 2, t61, -0.2e1 * t84, 0, 0, 0, t65 * t96, t69 * t96, t39 ^ 2, -0.2e1 * t39 * t38, 0, 0, 0, t38 * t100, t39 * t100, t76 ^ 2, -0.2e1 * t76 * t74, 0, 0, 0, t74 * t101, t76 * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, t54, 0, 0, 0, 0, 0, t50, -t88, 0, 0, 0, 0, 0, t35, -t90, 0, 0, 0, 0, 0, t102, -t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, -t83, t66, t22, -t23, 0, 0, -t29, t28, t66, t66 * t56 + t8, -t85 + (-t15 - t93) * t64, 0, 0, t12, -t11, t66, t30 * t66 + t1, -t31 * t66 - t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t65, 0, t52, -t65 * t71, 0, 0, t39, -t38, 0, t24, -t25, 0, 0, t76, -t74, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t65, 0, 0, 0, 0, 0, t39, -t38, 0, 0, 0, 0, 0, t76, -t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t56, -0.2e1 * t94, 0, 0, 0, 0, 1, 0.2e1 * t30, -0.2e1 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, t28, t66, t8, -t9, 0, 0, t12, -t11, t66, t66 * t55 + t1, -t91 + (-t4 - t92) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t38, 0, t24, -t25, 0, 0, t76, -t74, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t38, 0, 0, 0, 0, 0, t76, -t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t56, -t94, 0, 0, 0, 0, 1, t30 + t55, -t79 + (-pkin(5) - t51) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t55, -0.2e1 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, t66, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t74, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t30, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t55, -t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
