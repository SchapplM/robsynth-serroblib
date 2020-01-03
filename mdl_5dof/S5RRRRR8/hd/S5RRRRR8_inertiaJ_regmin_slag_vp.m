% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x31]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t58 = sin(qJ(3));
t59 = sin(qJ(2));
t75 = cos(qJ(3));
t76 = cos(qJ(2));
t37 = t58 * t59 - t75 * t76;
t89 = -0.2e1 * t37;
t88 = 0.2e1 * t37;
t64 = t75 * pkin(2);
t50 = -t64 - pkin(3);
t61 = cos(qJ(4));
t79 = t61 * pkin(4);
t41 = t50 - t79;
t87 = 0.2e1 * t41;
t51 = -pkin(3) - t79;
t86 = 0.2e1 * t51;
t52 = -pkin(2) * t76 - pkin(1);
t85 = 0.2e1 * t52;
t84 = t37 * pkin(4);
t56 = sin(qJ(5));
t83 = t56 * pkin(4);
t82 = t58 * pkin(2);
t60 = cos(qJ(5));
t81 = t60 * pkin(4);
t39 = t58 * t76 + t59 * t75;
t18 = t37 * pkin(3) - t39 * pkin(8) + t52;
t57 = sin(qJ(4));
t43 = (-pkin(6) - pkin(7)) * t59;
t65 = t76 * pkin(6);
t45 = pkin(7) * t76 + t65;
t27 = t58 * t43 + t45 * t75;
t71 = t61 * t27;
t6 = t71 + (-pkin(9) * t39 + t18) * t57;
t80 = t60 * t6;
t78 = t61 * pkin(8);
t77 = pkin(3) - t50;
t25 = -t75 * t43 + t45 * t58;
t74 = t25 * t61;
t73 = t57 * t39;
t72 = t57 * t61;
t70 = t61 * t39;
t49 = pkin(8) + t82;
t69 = t61 * t49;
t68 = t41 + t51;
t67 = 0.2e1 * t76;
t66 = t39 * t89;
t7 = t61 * t18 - t27 * t57;
t5 = -pkin(9) * t70 + t7 + t84;
t1 = t60 * t5 - t56 * t6;
t63 = -pkin(3) * t39 - pkin(8) * t37;
t62 = -t37 * t49 + t39 * t50;
t38 = t56 * t61 + t57 * t60;
t36 = t56 * t57 - t60 * t61;
t55 = t61 ^ 2;
t54 = t57 ^ 2;
t53 = t61 * pkin(9);
t46 = 0.2e1 * t72;
t44 = t53 + t78;
t42 = (-pkin(8) - pkin(9)) * t57;
t35 = t39 ^ 2;
t34 = t38 ^ 2;
t33 = t37 ^ 2;
t32 = t53 + t69;
t31 = (-pkin(9) - t49) * t57;
t30 = t61 * t37;
t29 = t57 * t37;
t28 = t57 * t70;
t26 = t42 * t56 + t44 * t60;
t24 = t42 * t60 - t44 * t56;
t23 = t38 * t37;
t22 = t36 * t37;
t21 = -0.2e1 * t38 * t36;
t20 = t25 * t57;
t19 = (-t54 + t55) * t39;
t17 = t31 * t56 + t32 * t60;
t16 = t31 * t60 - t32 * t56;
t14 = t36 * t39;
t13 = t38 * t39;
t12 = pkin(4) * t73 + t25;
t11 = t14 * t38;
t10 = t12 * t38;
t9 = t12 * t36;
t8 = t18 * t57 + t71;
t3 = -t13 * t38 + t14 * t36;
t2 = t5 * t56 + t80;
t4 = [1, 0, 0, t59 ^ 2, t59 * t67, 0, 0, 0, pkin(1) * t67, -0.2e1 * pkin(1) * t59, t35, t66, 0, 0, 0, t37 * t85, t39 * t85, t55 * t35, -0.2e1 * t35 * t72, t70 * t88, t57 * t66, t33, 0.2e1 * t25 * t73 + 0.2e1 * t37 * t7, 0.2e1 * t25 * t70 - 0.2e1 * t37 * t8, t14 ^ 2, 0.2e1 * t14 * t13, -t14 * t88, t13 * t89, t33, 0.2e1 * t1 * t37 + 0.2e1 * t12 * t13, -0.2e1 * t12 * t14 - 0.2e1 * t2 * t37; 0, 0, 0, 0, 0, t59, t76, 0, -t59 * pkin(6), -t65, 0, 0, t39, -t37, 0, -t25, -t27, t28, t19, t29, t30, 0, t57 * t62 - t74, t61 * t62 + t20, -t11, t3, t23, -t22, 0, t13 * t41 + t16 * t37 + t9, -t14 * t41 - t17 * t37 + t10; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t64, -0.2e1 * t82, t54, t46, 0, 0, 0, -0.2e1 * t50 * t61, 0.2e1 * t50 * t57, t34, t21, 0, 0, 0, t36 * t87, t38 * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t37, 0, -t25, -t27, t28, t19, t29, t30, 0, t57 * t63 - t74, t61 * t63 + t20, -t11, t3, t23, -t22, 0, t13 * t51 + t24 * t37 + t9, -t14 * t51 - t26 * t37 + t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t64, -t82, t54, t46, 0, 0, 0, t77 * t61, -t77 * t57, t34, t21, 0, 0, 0, t68 * t36, t68 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t54, t46, 0, 0, 0, 0.2e1 * pkin(3) * t61, -0.2e1 * pkin(3) * t57, t34, t21, 0, 0, 0, t36 * t86, t38 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t73, t37, t7, -t8, 0, 0, -t14, -t13, t37, t37 * t81 + t1, -t80 + (-t5 - t84) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t61, 0, -t57 * t49, -t69, 0, 0, t38, -t36, 0, t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t61, 0, -t57 * pkin(8), -t78, 0, 0, t38, -t36, 0, t24, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t81, -0.2e1 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t13, t37, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t36, 0, t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t36, 0, t24, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t81, -t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t4;
