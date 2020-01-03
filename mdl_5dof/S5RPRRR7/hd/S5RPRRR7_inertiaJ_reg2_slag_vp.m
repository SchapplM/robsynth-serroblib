% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t44 = cos(qJ(4));
t33 = -t44 * pkin(4) - pkin(3);
t80 = 0.2e1 * t33;
t42 = sin(qJ(3));
t79 = 0.2e1 * t42;
t45 = cos(qJ(3));
t78 = -0.2e1 * t45;
t77 = -pkin(8) - pkin(7);
t76 = pkin(3) * t44;
t38 = sin(pkin(9));
t75 = t38 * pkin(1);
t39 = cos(pkin(9));
t74 = t39 * pkin(1);
t40 = sin(qJ(5));
t73 = t40 * pkin(4);
t43 = cos(qJ(5));
t72 = t43 * pkin(4);
t28 = -pkin(2) - t74;
t70 = t45 * pkin(3);
t18 = -t42 * pkin(7) + t28 - t70;
t41 = sin(qJ(4));
t27 = pkin(6) + t75;
t56 = t45 * t27;
t51 = t44 * t56;
t5 = t51 + (-pkin(8) * t42 + t18) * t41;
t71 = t43 * t5;
t69 = t45 * pkin(4);
t21 = t40 * t44 + t43 * t41;
t13 = t21 * t42;
t68 = t21 * t13;
t67 = t21 * t45;
t66 = t27 * t41;
t34 = t41 ^ 2;
t65 = t34 * t42;
t35 = t42 ^ 2;
t64 = t35 * t27;
t63 = t41 * t42;
t62 = t41 * t44;
t61 = t41 * t45;
t60 = t42 * t27;
t59 = t44 * t42;
t58 = t44 * t45;
t19 = t40 * t41 - t43 * t44;
t57 = t45 * t19;
t36 = t44 ^ 2;
t55 = t34 + t36;
t37 = t45 ^ 2;
t54 = t35 + t37;
t53 = t45 * t79;
t52 = t41 * t59;
t16 = t44 * t18;
t4 = -pkin(8) * t59 + t16 + (-pkin(4) - t66) * t45;
t1 = t43 * t4 - t40 * t5;
t50 = t55 * pkin(7);
t7 = -t41 * t56 + t16;
t8 = t41 * t18 + t51;
t49 = -t7 * t41 + t8 * t44;
t32 = t36 * t42;
t31 = t36 * t35;
t29 = t34 * t35;
t26 = t27 ^ 2;
t24 = t77 * t44;
t23 = t77 * t41;
t22 = t35 * t26;
t17 = (pkin(4) * t41 + t27) * t42;
t15 = -t40 * t63 + t43 * t59;
t12 = t15 ^ 2;
t11 = t13 ^ 2;
t10 = t40 * t23 - t43 * t24;
t9 = t43 * t23 + t40 * t24;
t6 = t15 * t19;
t2 = t40 * t4 + t71;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t74, -0.2e1 * t75, 0, (t38 ^ 2 + t39 ^ 2) * pkin(1) ^ 2, t35, t53, 0, t37, 0, 0, t28 * t78, t28 * t79, 0.2e1 * t54 * t27, t37 * t26 + t28 ^ 2 + t22, t31, -0.2e1 * t35 * t62, -0.2e1 * t42 * t58, t29, t41 * t53, t37, 0.2e1 * t41 * t64 - 0.2e1 * t7 * t45, 0.2e1 * t44 * t64 + 0.2e1 * t8 * t45, (-t41 * t8 - t44 * t7) * t79, t7 ^ 2 + t8 ^ 2 + t22, t12, -0.2e1 * t15 * t13, t15 * t78, t11, -t13 * t78, t37, -0.2e1 * t1 * t45 + 0.2e1 * t17 * t13, 0.2e1 * t17 * t15 + 0.2e1 * t2 * t45, -0.2e1 * t1 * t15 - 0.2e1 * t2 * t13, t1 ^ 2 + t17 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t49 - t56) * t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t13 + t2 * t15 - t17 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 + t29 + t37, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 + t11 + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t45, 0, -t60, -t56, 0, 0, t52, t32 - t65, -t61, -t52, -t58, 0, -t27 * t59 + (-pkin(3) * t42 + pkin(7) * t45) * t41, pkin(7) * t58 + (t66 - t76) * t42, t49, -pkin(3) * t60 + t49 * pkin(7), t15 * t21, -t6 - t68, -t67, t13 * t19, t57, 0, t33 * t13 + t17 * t19 - t9 * t45, t10 * t45 + t33 * t15 + t17 * t21, -t1 * t21 - t10 * t13 - t9 * t15 - t2 * t19, t1 * t9 + t2 * t10 + t17 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t42, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t61, t32 + t65, t42 * t50 + t70, 0, 0, 0, 0, 0, 0, -t57, -t67, -t6 + t68, t15 * t10 - t13 * t9 - t45 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t34, 0.2e1 * t62, 0, t36, 0, 0, 0.2e1 * t76, -0.2e1 * pkin(3) * t41, 0.2e1 * t50, t55 * pkin(7) ^ 2 + pkin(3) ^ 2, t21 ^ 2, -0.2e1 * t21 * t19, 0, t19 ^ 2, 0, 0, t19 * t80, t21 * t80, -0.2e1 * t10 * t19 - 0.2e1 * t9 * t21, t10 ^ 2 + t33 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, -t63, -t45, t7, -t8, 0, 0, 0, 0, t15, 0, -t13, -t45, -t43 * t69 + t1, -t71 + (-t4 + t69) * t40, (-t13 * t40 - t15 * t43) * pkin(4), (t1 * t43 + t2 * t40) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t59, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, (-t13 * t43 + t15 * t40) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, t44, 0, -t41 * pkin(7), -t44 * pkin(7), 0, 0, 0, 0, t21, 0, -t19, 0, t9, -t10, (-t19 * t40 - t21 * t43) * pkin(4), (t10 * t40 + t43 * t9) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t72, -0.2e1 * t73, 0, (t40 ^ 2 + t43 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, -t13, -t45, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, -t19, 0, t9, -t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t72, -t73, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t3;
