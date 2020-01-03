% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP11_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t42 = sin(qJ(4));
t38 = t42 ^ 2;
t44 = cos(qJ(4));
t39 = t44 ^ 2;
t31 = t38 + t39;
t40 = sin(pkin(8));
t41 = cos(pkin(8));
t43 = sin(qJ(3));
t65 = cos(qJ(3));
t27 = t65 * t40 + t43 * t41;
t74 = -0.2e1 * t27;
t62 = pkin(6) + qJ(2);
t30 = t62 * t41;
t55 = t62 * t40;
t13 = t43 * t30 + t65 * t55;
t73 = t13 ^ 2;
t25 = t43 * t40 - t65 * t41;
t22 = t25 ^ 2;
t33 = -t41 * pkin(2) - pkin(1);
t72 = 0.2e1 * t33;
t71 = 0.2e1 * t41;
t70 = -0.2e1 * t42;
t69 = pkin(7) * t25;
t68 = t25 * pkin(4);
t67 = t42 * pkin(7);
t66 = t44 * pkin(7);
t15 = t65 * t30 - t43 * t55;
t8 = t25 * pkin(3) - t27 * pkin(7) + t33;
t4 = t44 * t15 + t42 * t8;
t18 = t42 * t25;
t64 = t42 * t27;
t63 = t42 * t44;
t20 = t44 * t25;
t21 = t44 * t27;
t61 = t31 * pkin(7) ^ 2;
t36 = t40 ^ 2;
t37 = t41 ^ 2;
t60 = t36 + t37;
t59 = t25 * qJ(5);
t58 = t25 * t64;
t23 = t27 ^ 2;
t57 = t23 * t63;
t56 = t42 * t15 - t44 * t8;
t54 = -pkin(3) * t27 - t69;
t1 = t59 + t4;
t2 = t56 - t68;
t53 = t1 * t44 + t2 * t42;
t52 = t1 * t42 - t2 * t44;
t51 = t4 * t42 - t44 * t56;
t50 = t4 * t44 + t42 * t56;
t48 = t44 * pkin(4) + t42 * qJ(5);
t29 = -pkin(3) - t48;
t49 = -t27 * t29 + t69;
t47 = pkin(4) * t42 - t44 * qJ(5);
t28 = 0.2e1 * t31 * pkin(7);
t19 = t39 * t23;
t17 = t38 * t23;
t16 = t42 * t21;
t12 = 0.2e1 * t25 * t21;
t10 = t31 * t27;
t9 = (t38 - t39) * t27;
t5 = t47 * t27 + t13;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t36, t40 * t71, 0, t37, 0, 0, pkin(1) * t71, -0.2e1 * pkin(1) * t40, 0.2e1 * t60 * qJ(2), t60 * qJ(2) ^ 2 + pkin(1) ^ 2, t23, t25 * t74, 0, t22, 0, 0, t25 * t72, t27 * t72, 0.2e1 * t13 * t27 - 0.2e1 * t15 * t25, t15 ^ 2 + t33 ^ 2 + t73, t19, -0.2e1 * t57, t12, t17, -0.2e1 * t58, t22, 0.2e1 * t13 * t64 - 0.2e1 * t25 * t56, 0.2e1 * t13 * t21 - 0.2e1 * t4 * t25, t51 * t74, t4 ^ 2 + t56 ^ 2 + t73, t19, t12, 0.2e1 * t57, t22, 0.2e1 * t58, t17, -0.2e1 * t2 * t25 + 0.2e1 * t5 * t64, t52 * t74, 0.2e1 * t1 * t25 - 0.2e1 * t5 * t21, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t40, 0, -pkin(1), 0, 0, 0, 0, 0, 0, t25, t27, 0, t33, 0, 0, 0, 0, 0, 0, t20, -t18, -t10, t51, 0, 0, 0, 0, 0, 0, t20, -t10, t18, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, -t25, 0, -t13, -t15, 0, 0, t16, -t9, t18, -t16, t20, 0, -t13 * t44 + t54 * t42, t13 * t42 + t54 * t44, t50, -t13 * pkin(3) + t50 * pkin(7), t16, t18, t9, 0, -t20, -t16, -t49 * t42 - t5 * t44, t53, -t5 * t42 + t49 * t44, t53 * pkin(7) + t5 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t38, 0.2e1 * t63, 0, t39, 0, 0, 0.2e1 * pkin(3) * t44, pkin(3) * t70, t28, pkin(3) ^ 2 + t61, t38, 0, -0.2e1 * t63, 0, 0, t39, -0.2e1 * t29 * t44, t28, t29 * t70, t29 ^ 2 + t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, -t64, t25, -t56, -t4, 0, 0, 0, t21, 0, t25, t64, 0, -t56 + 0.2e1 * t68, -t48 * t27, 0.2e1 * t59 + t4, -t2 * pkin(4) + t1 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t42, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, t42, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t44, 0, -t67, -t66, 0, 0, 0, t42, 0, 0, -t44, 0, -t67, -t47, t66, -t47 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t21, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t3;
