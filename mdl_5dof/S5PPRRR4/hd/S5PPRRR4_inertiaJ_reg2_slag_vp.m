% Calculate inertial parameters regressor of joint inertia matrix for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRRR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:44
% EndTime: 2019-12-05 15:19:46
% DurationCPUTime: 0.65s
% Computational Cost: add. (381->95), mult. (1095->178), div. (0->0), fcn. (1369->12), ass. (0->61)
t31 = sin(pkin(6));
t32 = sin(pkin(5));
t33 = cos(pkin(11));
t34 = cos(pkin(6));
t35 = cos(pkin(5));
t14 = -t32 * t33 * t31 + t35 * t34;
t37 = sin(qJ(4));
t40 = cos(qJ(4));
t30 = sin(pkin(11));
t38 = sin(qJ(3));
t41 = cos(qJ(3));
t57 = t33 * t34;
t59 = t31 * t38;
t8 = t35 * t59 + (t30 * t41 + t38 * t57) * t32;
t3 = -t14 * t40 + t8 * t37;
t69 = t3 ^ 2;
t58 = t31 * t41;
t6 = -t35 * t58 + (t30 * t38 - t41 * t57) * t32;
t68 = t6 ^ 2;
t15 = -t40 * t34 + t37 * t59;
t67 = t15 ^ 2;
t66 = -0.2e1 * t37;
t65 = 0.2e1 * t40;
t39 = cos(qJ(5));
t64 = pkin(4) * t39;
t27 = t37 ^ 2;
t63 = t27 * pkin(8);
t62 = t3 * t15;
t61 = t37 * pkin(8);
t60 = t15 * t37;
t36 = sin(qJ(5));
t56 = t36 * t37;
t55 = t36 * t39;
t54 = t36 * t40;
t53 = t39 * t37;
t52 = t39 * t40;
t26 = t36 ^ 2;
t28 = t39 ^ 2;
t51 = t26 + t28;
t50 = t37 * t65;
t49 = t36 * t53;
t5 = t14 * t37 + t8 * t40;
t1 = -t5 * t36 + t6 * t39;
t2 = t6 * t36 + t5 * t39;
t48 = -t1 * t36 + t2 * t39;
t47 = t3 * t37 + t5 * t40;
t17 = t37 * t34 + t40 * t59;
t10 = t39 * t17 - t36 * t58;
t9 = -t36 * t17 - t39 * t58;
t46 = t10 * t39 - t9 * t36;
t19 = -t40 * pkin(4) - t37 * pkin(9) - pkin(3);
t11 = -pkin(8) * t54 + t39 * t19;
t12 = pkin(8) * t52 + t36 * t19;
t45 = -t11 * t36 + t12 * t39;
t44 = t17 * t40 + t60;
t43 = pkin(8) ^ 2;
t29 = t40 ^ 2;
t24 = t31 ^ 2;
t23 = t27 * t43;
t21 = t24 * t41 ^ 2;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 ^ 2 + (t30 ^ 2 + t33 ^ 2) * t32 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 ^ 2 + t8 ^ 2 + t68, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 + t68 + t69, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t34 + (t38 * t8 - t41 * t6) * t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t17 - t6 * t58 + t62, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t9 + t2 * t10 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t38 ^ 2 + t34 ^ 2 + t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 ^ 2 + t21 + t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t9 ^ 2 + t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t8, 0, 0, 0, 0, 0, 0, 0, 0, -t6 * t40, t6 * t37, t47, -t6 * pkin(3) + t47 * pkin(8), 0, 0, 0, 0, 0, 0, -t1 * t40 + t3 * t56, t2 * t40 + t3 * t53, (-t1 * t39 - t2 * t36) * t37, t1 * t11 + t2 * t12 + t3 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t59, 0, 0, 0, 0, 0, 0, 0, 0, t40 * t58, -t37 * t58, t44, pkin(3) * t58 + pkin(8) * t44, 0, 0, 0, 0, 0, 0, t15 * t56 - t9 * t40, t10 * t40 + t15 * t53, (-t10 * t36 - t39 * t9) * t37, pkin(8) * t60 + t10 * t12 + t9 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t27, t50, 0, t29, 0, 0, pkin(3) * t65, pkin(3) * t66, 0.2e1 * (t27 + t29) * pkin(8), pkin(3) ^ 2 + t29 * t43 + t23, t28 * t27, -0.2e1 * t27 * t55, t52 * t66, t26 * t27, t36 * t50, t29, -0.2e1 * t11 * t40 + 0.2e1 * t36 * t63, 0.2e1 * t12 * t40 + 0.2e1 * t39 * t63, 0.2e1 * (-t11 * t39 - t12 * t36) * t37, t11 ^ 2 + t12 ^ 2 + t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t5, 0, 0, 0, 0, 0, 0, 0, 0, -t3 * t39, t3 * t36, t48, -t3 * pkin(4) + t48 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t17, 0, 0, 0, 0, 0, 0, 0, 0, -t15 * t39, t15 * t36, t46, -t15 * pkin(4) + t46 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t40, 0, -t61, -t40 * pkin(8), 0, 0, t49, (-t26 + t28) * t37, -t54, -t49, -t52, 0, -pkin(8) * t53 + (-pkin(4) * t37 + pkin(9) * t40) * t36, pkin(9) * t52 + (pkin(8) * t36 - t64) * t37, t45, -pkin(4) * t61 + pkin(9) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t26, 0.2e1 * t55, 0, t28, 0, 0, 0.2e1 * t64, -0.2e1 * pkin(4) * t36, 0.2e1 * t51 * pkin(9), t51 * pkin(9) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, -t56, -t40, t11, -t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, t39, 0, -t36 * pkin(9), -t39 * pkin(9), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t4;
