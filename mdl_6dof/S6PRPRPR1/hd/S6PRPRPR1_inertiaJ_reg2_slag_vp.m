% Calculate inertial parameters regressor of joint inertia matrix for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRPR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_inertiaJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t48 = sin(pkin(12));
t51 = cos(pkin(12));
t55 = sin(qJ(4));
t58 = cos(qJ(4));
t32 = t48 * t58 + t51 * t55;
t93 = -0.2e1 * t32;
t49 = sin(pkin(11));
t50 = sin(pkin(6));
t52 = cos(pkin(11));
t56 = sin(qJ(2));
t59 = cos(qJ(2));
t19 = (t49 * t59 + t52 * t56) * t50;
t53 = cos(pkin(6));
t14 = -t19 * t55 + t53 * t58;
t15 = t19 * t58 + t53 * t55;
t5 = -t51 * t14 + t48 * t15;
t92 = t5 ^ 2;
t85 = t49 * pkin(2);
t39 = pkin(8) + t85;
t72 = qJ(5) + t39;
t26 = t72 * t58;
t68 = t72 * t55;
t9 = t48 * t26 + t51 * t68;
t91 = t9 ^ 2;
t78 = t50 * t59;
t79 = t50 * t56;
t17 = t49 * t79 - t52 * t78;
t16 = t17 ^ 2;
t30 = t48 * t55 - t51 * t58;
t90 = t30 ^ 2;
t82 = t52 * pkin(2);
t41 = -pkin(3) - t82;
t33 = -t58 * pkin(4) + t41;
t89 = 0.2e1 * t33;
t88 = 0.2e1 * t55;
t87 = t5 * t9;
t86 = t48 * pkin(4);
t84 = t5 * t30;
t83 = t51 * pkin(4);
t81 = t9 * t30;
t54 = sin(qJ(6));
t44 = t54 ^ 2;
t80 = t44 * t32;
t22 = t54 * t30;
t77 = t54 * t32;
t57 = cos(qJ(6));
t76 = t54 * t57;
t75 = t57 * t32;
t46 = t57 ^ 2;
t74 = t44 + t46;
t45 = t55 ^ 2;
t47 = t58 ^ 2;
t73 = t45 + t47;
t71 = t30 * t93;
t70 = t54 * t75;
t38 = pkin(9) + t86;
t69 = t74 * t38;
t7 = t48 * t14 + t51 * t15;
t1 = t17 * t57 - t54 * t7;
t2 = t17 * t54 + t57 * t7;
t67 = t1 * t57 + t2 * t54;
t66 = -t1 * t54 + t2 * t57;
t11 = t51 * t26 - t48 * t68;
t8 = t30 * pkin(5) - t32 * pkin(9) + t33;
t3 = -t54 * t11 + t57 * t8;
t4 = t57 * t11 + t54 * t8;
t65 = t3 * t57 + t4 * t54;
t64 = -t3 * t54 + t4 * t57;
t63 = -t14 * t55 + t15 * t58;
t40 = -pkin(5) - t83;
t62 = -t30 * t38 + t32 * t40;
t43 = t53 ^ 2;
t29 = t32 ^ 2;
t25 = t57 * t30;
t24 = t46 * t32;
t23 = t46 * t29;
t21 = t44 * t29;
t12 = -t24 - t80;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 + (t56 ^ 2 + t59 ^ 2) * t50 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 ^ 2 + t16 + t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 ^ 2 + t15 ^ 2 + t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t16 + t92, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t79, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t19, 0 (-t17 * t52 + t19 * t49) * pkin(2), 0, 0, 0, 0, 0, 0, -t17 * t58, t17 * t55, t63, t17 * t41 + t39 * t63, 0, 0, 0, 0, 0, 0, t17 * t30, t17 * t32, -t7 * t30 + t5 * t32, t7 * t11 + t17 * t33 + t87, 0, 0, 0, 0, 0, 0, t1 * t30 + t5 * t77, -t2 * t30 + t5 * t75, -t67 * t32, t1 * t3 + t2 * t4 + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t82, -0.2e1 * t85, 0 (t49 ^ 2 + t52 ^ 2) * pkin(2) ^ 2, t45, t58 * t88, 0, t47, 0, 0, -0.2e1 * t41 * t58, t41 * t88, 0.2e1 * t73 * t39, t39 ^ 2 * t73 + t41 ^ 2, t29, t71, 0, t90, 0, 0, t30 * t89, t32 * t89, -0.2e1 * t11 * t30 + 0.2e1 * t9 * t32, t11 ^ 2 + t33 ^ 2 + t91, t23, -0.2e1 * t29 * t76, 0.2e1 * t30 * t75, t21, t54 * t71, t90, 0.2e1 * t3 * t30 + 0.2e1 * t77 * t9, -0.2e1 * t4 * t30 + 0.2e1 * t75 * t9, t65 * t93, t3 ^ 2 + t4 ^ 2 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t58 + t15 * t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t32 + t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t32 + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t32 + t81, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64 * t32 + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 + t90, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 + t21 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t15, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t7, 0 (t48 * t7 - t5 * t51) * pkin(4), 0, 0, 0, 0, 0, 0, -t5 * t57, t5 * t54, t66, t38 * t66 + t5 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, t58, 0, -t55 * t39, -t58 * t39, 0, 0, 0, 0, t32, 0, -t30, 0, -t9, -t11 (-t30 * t48 - t32 * t51) * pkin(4) (t11 * t48 - t51 * t9) * pkin(4), t70, t24 - t80, t22, -t70, t25, 0, t54 * t62 - t9 * t57, t9 * t54 + t57 * t62, t64, t38 * t64 + t9 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t55, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t32, 0 (-t30 * t51 + t32 * t48) * pkin(4), 0, 0, 0, 0, 0, 0, -t25, t22, -t12, t30 * t40 + t32 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t83, -0.2e1 * t86, 0 (t48 ^ 2 + t51 ^ 2) * pkin(4) ^ 2, t44, 0.2e1 * t76, 0, t46, 0, 0, -0.2e1 * t40 * t57, 0.2e1 * t40 * t54, 0.2e1 * t69, t38 ^ 2 * t74 + t40 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t32, 0, t33, 0, 0, 0, 0, 0, 0, t25, -t22, t12, t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, -t77, t30, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t75, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, t57, 0, -t54 * t38, -t57 * t38, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t54, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t6;
