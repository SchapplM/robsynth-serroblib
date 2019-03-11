% Calculate inertial parameters regressor of joint inertia matrix for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PPPRRR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_inertiaJ_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t46 = sin(pkin(14));
t47 = sin(pkin(13));
t50 = sin(pkin(6));
t51 = cos(pkin(14));
t52 = cos(pkin(13));
t54 = cos(pkin(7));
t80 = t52 * t54;
t49 = sin(pkin(7));
t55 = cos(pkin(6));
t81 = t49 * t55;
t17 = t51 * t81 + (-t46 * t47 + t51 * t80) * t50;
t29 = -t50 * t52 * t49 + t55 * t54;
t48 = sin(pkin(8));
t53 = cos(pkin(8));
t12 = -t17 * t48 + t29 * t53;
t57 = sin(qJ(5));
t60 = cos(qJ(5));
t18 = t50 * t47 * t51 + (t50 * t80 + t81) * t46;
t58 = sin(qJ(4));
t61 = cos(qJ(4));
t8 = t18 * t61 + (t17 * t53 + t29 * t48) * t58;
t3 = -t12 * t60 + t8 * t57;
t100 = t3 ^ 2;
t79 = t53 * t61;
t83 = t48 * t61;
t6 = -t17 * t79 + t18 * t58 - t29 * t83;
t99 = t6 ^ 2;
t84 = t48 * t58;
t21 = t54 * t84 + (t51 * t53 * t58 + t46 * t61) * t49;
t82 = t49 * t51;
t28 = -t48 * t82 + t53 * t54;
t13 = t57 * t21 - t60 * t28;
t98 = t13 ^ 2;
t19 = t58 * t49 * t46 - t54 * t83 - t79 * t82;
t97 = t19 ^ 2;
t30 = -t60 * t53 + t57 * t84;
t96 = t30 ^ 2;
t95 = -0.2e1 * t57;
t94 = 0.2e1 * t60;
t59 = cos(qJ(6));
t93 = pkin(5) * t59;
t92 = t3 * t13;
t91 = t3 * t30;
t43 = t57 ^ 2;
t90 = t43 * pkin(10);
t89 = t57 * pkin(10);
t88 = t6 * t19;
t87 = t13 * t30;
t86 = t13 * t57;
t85 = t30 * t57;
t56 = sin(qJ(6));
t78 = t56 * t57;
t77 = t59 * t56;
t76 = t59 * t57;
t75 = t59 * t60;
t74 = t60 * t56;
t42 = t56 ^ 2;
t44 = t59 ^ 2;
t73 = t42 + t44;
t72 = t57 * t94;
t71 = t56 * t76;
t5 = t12 * t57 + t8 * t60;
t1 = -t5 * t56 + t6 * t59;
t2 = t5 * t59 + t6 * t56;
t70 = -t1 * t56 + t2 * t59;
t69 = t3 * t57 + t5 * t60;
t15 = t60 * t21 + t57 * t28;
t10 = t59 * t15 + t56 * t19;
t9 = -t56 * t15 + t59 * t19;
t68 = t10 * t59 - t9 * t56;
t67 = t15 * t60 + t86;
t32 = t57 * t53 + t60 * t84;
t22 = -t56 * t32 - t59 * t83;
t23 = t59 * t32 - t56 * t83;
t66 = -t22 * t56 + t23 * t59;
t34 = -t60 * pkin(5) - t57 * pkin(11) - pkin(4);
t25 = -pkin(10) * t74 + t59 * t34;
t26 = pkin(10) * t75 + t56 * t34;
t65 = -t25 * t56 + t26 * t59;
t64 = t32 * t60 + t85;
t63 = pkin(10) ^ 2;
t45 = t60 ^ 2;
t39 = t48 ^ 2;
t38 = t43 * t63;
t36 = t39 * t61 ^ 2;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 ^ 2 + (t47 ^ 2 + t52 ^ 2) * t50 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 ^ 2 + t18 ^ 2 + t29 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 ^ 2 + t8 ^ 2 + t99, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 + t100 + t99, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t54 + (t17 * t51 + t18 * t46) * t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t28 + t8 * t21 + t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t15 + t88 + t92, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t9 + t2 * t10 + t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 ^ 2 + (t46 ^ 2 + t51 ^ 2) * t49 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 ^ 2 + t28 ^ 2 + t97, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 ^ 2 + t97 + t98, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t9 ^ 2 + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t53 + (t58 * t8 - t6 * t61) * t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t32 - t6 * t83 + t91, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t22 + t2 * t23 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t53 + (-t19 * t61 + t21 * t58) * t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * t32 - t19 * t83 + t87, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t23 + t9 * t22 + t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t58 ^ 2 + t53 ^ 2 + t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 ^ 2 + t36 + t96, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 ^ 2 + t23 ^ 2 + t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t8, 0, 0, 0, 0, 0, 0, 0, 0, -t6 * t60, t6 * t57, t69, -t6 * pkin(4) + t69 * pkin(10), 0, 0, 0, 0, 0, 0, -t1 * t60 + t3 * t78, t2 * t60 + t3 * t76 (-t1 * t59 - t2 * t56) * t57, t1 * t25 + t2 * t26 + t3 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t21, 0, 0, 0, 0, 0, 0, 0, 0, -t19 * t60, t19 * t57, t67, -t19 * pkin(4) + t67 * pkin(10), 0, 0, 0, 0, 0, 0, t13 * t78 - t9 * t60, t10 * t60 + t13 * t76 (-t10 * t56 - t59 * t9) * t57, pkin(10) * t86 + t10 * t26 + t9 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t84, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t83, -t57 * t83, t64, pkin(4) * t83 + t64 * pkin(10), 0, 0, 0, 0, 0, 0, -t22 * t60 + t30 * t78, t23 * t60 + t30 * t76 (-t22 * t59 - t23 * t56) * t57, pkin(10) * t85 + t22 * t25 + t23 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t43, t72, 0, t45, 0, 0, pkin(4) * t94, pkin(4) * t95, 0.2e1 * (t43 + t45) * pkin(10), pkin(4) ^ 2 + t45 * t63 + t38, t44 * t43, -0.2e1 * t43 * t77, t75 * t95, t42 * t43, t56 * t72, t45, -0.2e1 * t25 * t60 + 0.2e1 * t56 * t90, 0.2e1 * t26 * t60 + 0.2e1 * t59 * t90, 0.2e1 * (-t25 * t59 - t26 * t56) * t57, t25 ^ 2 + t26 ^ 2 + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t5, 0, 0, 0, 0, 0, 0, 0, 0, -t3 * t59, t3 * t56, t70, -t3 * pkin(5) + t70 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t59, t13 * t56, t68, -t13 * pkin(5) + t68 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t32, 0, 0, 0, 0, 0, 0, 0, 0, -t30 * t59, t30 * t56, t66, -t30 * pkin(5) + t66 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, t60, 0, -t89, -t60 * pkin(10), 0, 0, t71 (-t42 + t44) * t57, -t74, -t71, -t75, 0, -pkin(10) * t76 + (-pkin(5) * t57 + pkin(11) * t60) * t56, pkin(11) * t75 + (pkin(10) * t56 - t93) * t57, t65, -pkin(5) * t89 + t65 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t42, 0.2e1 * t77, 0, t44, 0, 0, 0.2e1 * t93, -0.2e1 * pkin(5) * t56, 0.2e1 * t73 * pkin(11), t73 * pkin(11) ^ 2 + pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t23, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, -t78, -t60, t25, -t26, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, t59, 0, -t56 * pkin(11), -t59 * pkin(11), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t4;
