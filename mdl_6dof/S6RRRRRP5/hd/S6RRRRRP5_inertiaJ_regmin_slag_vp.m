% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRRP5
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
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t61 = sin(qJ(4));
t62 = sin(qJ(3));
t64 = cos(qJ(4));
t65 = cos(qJ(3));
t68 = -t61 * t62 + t64 * t65;
t52 = -t65 * pkin(3) - pkin(2);
t32 = -pkin(4) * t68 + t52;
t93 = 0.2e1 * t32;
t76 = t64 * t62;
t39 = t61 * t65 + t76;
t92 = 0.2e1 * t39;
t63 = sin(qJ(2));
t91 = -0.2e1 * t63;
t66 = cos(qJ(2));
t90 = -0.2e1 * t66;
t89 = 0.2e1 * t66;
t88 = pkin(8) + pkin(9);
t87 = pkin(2) * t65;
t86 = pkin(7) * t62;
t60 = sin(qJ(5));
t85 = t60 * pkin(4);
t84 = t61 * pkin(3);
t83 = t66 * pkin(3);
t82 = t66 * pkin(4);
t81 = cos(qJ(5));
t80 = t62 * t63;
t79 = t62 * t65;
t78 = t62 * t66;
t42 = -pkin(2) * t66 - pkin(8) * t63 - pkin(1);
t73 = t65 * t66;
t70 = pkin(7) * t73;
t26 = t70 + (-pkin(9) * t63 + t42) * t62;
t77 = t64 * t26;
t74 = t65 * t63;
t72 = t61 * t74 + t63 * t76;
t53 = t63 * pkin(7);
t41 = pkin(3) * t80 + t53;
t71 = t63 * t89;
t54 = t81 * pkin(4);
t38 = t65 * t42;
t24 = -pkin(9) * t74 + t38 + (-pkin(3) - t86) * t66;
t14 = t24 * t61 + t77;
t11 = -pkin(10) * t72 + t14;
t69 = t81 * t11;
t13 = t24 * t64 - t26 * t61;
t33 = t68 * t63;
t8 = -pkin(10) * t33 + t13 - t82;
t3 = -t11 * t60 + t8 * t81;
t43 = t88 * t62;
t44 = t88 * t65;
t27 = -t43 * t64 - t44 * t61;
t19 = -pkin(10) * t39 + t27;
t28 = -t61 * t43 + t64 * t44;
t20 = pkin(10) * t68 + t28;
t9 = t19 * t81 - t20 * t60;
t67 = t81 * t84;
t55 = t64 * pkin(3);
t51 = t55 + pkin(4);
t35 = t51 * t81 - t60 * t84;
t25 = pkin(4) * t72 + t41;
t4 = t60 * t8 + t69;
t10 = t19 * t60 + t20 * t81;
t59 = t66 ^ 2;
t58 = t65 ^ 2;
t57 = t63 ^ 2;
t56 = t62 ^ 2;
t50 = t54 + pkin(5);
t36 = t51 * t60 + t67;
t34 = pkin(5) + t35;
t31 = t42 * t62 + t70;
t30 = -pkin(7) * t78 + t38;
t23 = t39 * t81 + t60 * t68;
t22 = t39 * t60 - t68 * t81;
t17 = t33 * t81 - t60 * t72;
t16 = t33 * t60 + t72 * t81;
t15 = t22 * pkin(5) + t32;
t12 = pkin(5) * t16 + t25;
t6 = -qJ(6) * t22 + t10;
t5 = -qJ(6) * t23 + t9;
t2 = -qJ(6) * t16 + t4;
t1 = -pkin(5) * t66 - qJ(6) * t17 + t3;
t7 = [1, 0, 0, t57, t71, 0, 0, 0, pkin(1) * t89, pkin(1) * t91, t58 * t57, -0.2e1 * t57 * t79, t73 * t91, t62 * t71, t59, -0.2e1 * t30 * t66 + 0.2e1 * t57 * t86, 0.2e1 * pkin(7) * t57 * t65 + 0.2e1 * t31 * t66, t33 ^ 2, -0.2e1 * t33 * t72, t33 * t90, t72 * t89, t59, -0.2e1 * t13 * t66 + 0.2e1 * t41 * t72, 0.2e1 * t14 * t66 + 0.2e1 * t33 * t41, t17 ^ 2, -0.2e1 * t17 * t16, t17 * t90, t16 * t89, t59, 0.2e1 * t16 * t25 - 0.2e1 * t3 * t66, 0.2e1 * t17 * t25 + 0.2e1 * t4 * t66, -0.2e1 * t1 * t17 - 0.2e1 * t16 * t2, t1 ^ 2 + t12 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, t63, t66, 0, -t53, -t66 * pkin(7), t62 * t74 (-t56 + t58) * t63, -t78, -t73, 0, -pkin(7) * t74 + (-pkin(2) * t63 + pkin(8) * t66) * t62, pkin(8) * t73 + (t86 - t87) * t63, t33 * t39, t33 * t68 - t39 * t72, -t39 * t66, -t68 * t66, 0, -t27 * t66 - t41 * t68 + t52 * t72, t28 * t66 + t33 * t52 + t39 * t41, t17 * t23, -t16 * t23 - t17 * t22, -t23 * t66, t22 * t66, 0, t16 * t32 + t22 * t25 - t66 * t9, t10 * t66 + t17 * t32 + t23 * t25, -t1 * t23 - t16 * t6 - t17 * t5 - t2 * t22, t1 * t5 + t12 * t15 + t2 * t6; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t56, 0.2e1 * t79, 0, 0, 0, 0.2e1 * t87, -0.2e1 * pkin(2) * t62, t39 ^ 2, t68 * t92, 0, 0, 0, -0.2e1 * t52 * t68, t52 * t92, t23 ^ 2, -0.2e1 * t23 * t22, 0, 0, 0, t22 * t93, t23 * t93, -0.2e1 * t22 * t6 - 0.2e1 * t23 * t5, t15 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t80, -t66, t30, -t31, 0, 0, t33, -t72, -t66, -t64 * t83 + t13, -t77 + (-t24 + t83) * t61, 0, 0, t17, -t16, -t66, -t35 * t66 + t3, t36 * t66 - t4, -t16 * t36 - t17 * t34, t1 * t34 + t2 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t65, 0, -t62 * pkin(8), -t65 * pkin(8), 0, 0, t39, t68, 0, t27, -t28, 0, 0, t23, -t22, 0, t9, -t10, -t22 * t36 - t23 * t34, t34 * t5 + t36 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t55, -0.2e1 * t84, 0, 0, 0, 0, 1, 0.2e1 * t35, -0.2e1 * t36, 0, t34 ^ 2 + t36 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t72, -t66, t13, -t14, 0, 0, t17, -t16, -t66, -t54 * t66 + t3, -t69 + (-t8 + t82) * t60, -t16 * t85 - t17 * t50, t1 * t50 + t2 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t68, 0, t27, -t28, 0, 0, t23, -t22, 0, t9, -t10, -t22 * t85 - t23 * t50, t5 * t50 + t6 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t55, -t84, 0, 0, 0, 0, 1, t35 + t54, -t67 + (-pkin(4) - t51) * t60, 0, t34 * t50 + t36 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t54, -0.2e1 * t85, 0, pkin(4) ^ 2 * t60 ^ 2 + t50 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, -t66, t3, -t4, -pkin(5) * t17, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, 0, t9, -t10, -pkin(5) * t23, t5 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t35, -t36, 0, t34 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t54, -t85, 0, t50 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t7;
