% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t61 = cos(qJ(2));
t49 = -t61 * pkin(2) - pkin(1);
t86 = 0.2e1 * t49;
t56 = sin(qJ(5));
t85 = 0.2e1 * t56;
t59 = cos(qJ(5));
t84 = -0.2e1 * t59;
t83 = 0.2e1 * t61;
t82 = pkin(7) + pkin(8);
t81 = t56 * pkin(5);
t57 = sin(qJ(3));
t80 = t57 * pkin(2);
t79 = t59 * pkin(5);
t58 = sin(qJ(2));
t41 = t82 * t58;
t42 = t82 * t61;
t60 = cos(qJ(3));
t26 = t57 * t41 - t60 * t42;
t37 = t57 * t58 - t60 * t61;
t15 = -t37 * qJ(4) - t26;
t54 = sin(pkin(10));
t55 = cos(pkin(10));
t25 = -t60 * t41 - t57 * t42;
t38 = t57 * t61 + t60 * t58;
t63 = -t38 * qJ(4) + t25;
t10 = t54 * t15 - t55 * t63;
t78 = t10 * t59;
t20 = t55 * t37 + t54 * t38;
t17 = t56 * t20;
t21 = -t54 * t37 + t55 * t38;
t77 = t56 * t21;
t76 = t56 * t59;
t12 = t55 * t15 + t54 * t63;
t75 = t59 * t12;
t74 = t59 * t21;
t51 = t60 * pkin(2);
t48 = t51 + pkin(3);
t33 = t54 * t48 + t55 * t80;
t30 = pkin(9) + t33;
t73 = t59 * t30;
t46 = t54 * pkin(3) + pkin(9);
t72 = t59 * t46;
t32 = t55 * t48 - t54 * t80;
t29 = -pkin(4) - t32;
t47 = -t55 * pkin(3) - pkin(4);
t71 = t29 + t47;
t52 = t56 ^ 2;
t53 = t59 ^ 2;
t70 = t52 + t53;
t50 = t59 * qJ(6);
t28 = t37 * pkin(3) + t49;
t9 = t20 * pkin(4) - t21 * pkin(9) + t28;
t4 = -t56 * t12 + t59 * t9;
t1 = t20 * pkin(5) - t21 * t50 + t4;
t3 = t75 + (-qJ(6) * t21 + t9) * t56;
t69 = -t1 * t56 + t3 * t59;
t68 = t1 * t59 + t3 * t56;
t67 = -t20 * t30 + t21 * t29;
t66 = -t20 * t46 + t21 * t47;
t23 = (-qJ(6) - t30) * t56;
t24 = t50 + t73;
t65 = t23 * t59 + t24 * t56;
t35 = (-qJ(6) - t46) * t56;
t36 = t50 + t72;
t64 = t35 * t59 + t36 * t56;
t45 = 0.2e1 * t76;
t39 = t47 - t79;
t31 = t36 * t59;
t27 = t29 - t79;
t22 = t24 * t59;
t19 = t21 ^ 2;
t18 = t59 * t20;
t16 = t56 * t74;
t13 = (-t52 + t53) * t21;
t8 = t10 * t56;
t6 = pkin(5) * t77 + t10;
t5 = t56 * t9 + t75;
t2 = [1, 0, 0, t58 ^ 2, t58 * t83, 0, 0, 0, pkin(1) * t83, -0.2e1 * pkin(1) * t58, t38 ^ 2, -0.2e1 * t38 * t37, 0, 0, 0, t37 * t86, t38 * t86, 0.2e1 * t10 * t21 - 0.2e1 * t12 * t20, t10 ^ 2 + t12 ^ 2 + t28 ^ 2, t53 * t19, -0.2e1 * t19 * t76, 0.2e1 * t20 * t74, -0.2e1 * t20 * t77, t20 ^ 2, 0.2e1 * t10 * t77 + 0.2e1 * t4 * t20, 0.2e1 * t10 * t74 - 0.2e1 * t5 * t20, -0.2e1 * t68 * t21, t1 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, t58, t61, 0, -t58 * pkin(7), -t61 * pkin(7), 0, 0, t38, -t37, 0, t25, t26, -t33 * t20 - t32 * t21, -t10 * t32 + t12 * t33, t16, t13, t17, t18, 0, t56 * t67 - t78, t59 * t67 + t8, -t21 * t65 + t69, t1 * t23 + t3 * t24 + t6 * t27; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t51, -0.2e1 * t80, 0, t32 ^ 2 + t33 ^ 2, t52, t45, 0, 0, 0, t29 * t84, t29 * t85, -0.2e1 * t23 * t56 + 0.2e1 * t22, t23 ^ 2 + t24 ^ 2 + t27 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t37, 0, t25, t26 (-t20 * t54 - t21 * t55) * pkin(3) (-t10 * t55 + t12 * t54) * pkin(3), t16, t13, t17, t18, 0, t56 * t66 - t78, t59 * t66 + t8, -t21 * t64 + t69, t1 * t35 + t3 * t36 + t6 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t51, -t80, 0 (t32 * t55 + t33 * t54) * pkin(3), t52, t45, 0, 0, 0, -t71 * t59, t71 * t56, t22 + t31 + (-t23 - t35) * t56, t23 * t35 + t24 * t36 + t27 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t54 ^ 2 + t55 ^ 2) * pkin(3) ^ 2, t52, t45, 0, 0, 0, t47 * t84, t47 * t85, -0.2e1 * t35 * t56 + 0.2e1 * t31, t35 ^ 2 + t36 ^ 2 + t39 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, t18, -t17, -t70 * t21, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t77, t20, t4, -t5, -pkin(5) * t74, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t59, 0, -t56 * t30, -t73, -t81, t23 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t59, 0, -t56 * t46, -t72, -t81, t35 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t56, 0, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t2;
