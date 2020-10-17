% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:34:02
% EndTime: 2020-01-03 11:34:05
% DurationCPUTime: 0.52s
% Computational Cost: add. (325->49), mult. (608->89), div. (0->0), fcn. (618->8), ass. (0->46)
t38 = sin(pkin(9));
t36 = t38 ^ 2;
t40 = cos(pkin(9));
t37 = t40 ^ 2;
t47 = t36 + t37;
t48 = t47 * qJ(4);
t41 = cos(pkin(8));
t53 = t41 * pkin(1);
t31 = pkin(2) + t53;
t43 = sin(qJ(3));
t44 = cos(qJ(3));
t39 = sin(pkin(8));
t55 = t39 * pkin(1);
t19 = t44 * t31 - t43 * t55;
t18 = -pkin(3) - t19;
t54 = t40 * pkin(4);
t13 = t18 - t54;
t60 = 0.2e1 * t13;
t32 = -pkin(3) - t54;
t59 = 0.2e1 * t32;
t58 = 0.2e1 * t40;
t42 = sin(qJ(5));
t51 = cos(qJ(5));
t23 = t42 * t38 - t51 * t40;
t25 = t51 * t38 + t42 * t40;
t20 = t43 * t31 + t44 * t55;
t17 = qJ(4) + t20;
t11 = (-pkin(7) - t17) * t38;
t35 = t40 * pkin(7);
t12 = t40 * t17 + t35;
t3 = t51 * t11 - t42 * t12;
t4 = t42 * t11 + t51 * t12;
t57 = -t4 * t23 - t3 * t25;
t26 = (-pkin(7) - qJ(4)) * t38;
t27 = t40 * qJ(4) + t35;
t10 = t42 * t26 + t51 * t27;
t9 = t51 * t26 - t42 * t27;
t56 = -t10 * t23 - t9 * t25;
t52 = pkin(3) - t18;
t50 = t13 + t32;
t49 = t47 * t17;
t28 = t38 * t58;
t22 = t25 ^ 2;
t21 = t23 ^ 2;
t8 = -0.2e1 * t25 * t23;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t53, -0.2e1 * t55, 0, (t39 ^ 2 + t41 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t19, -0.2e1 * t20, 0, t19 ^ 2 + t20 ^ 2, t36, t28, 0, t37, 0, 0, -0.2e1 * t18 * t40, 0.2e1 * t18 * t38, 0.2e1 * t49, t47 * t17 ^ 2 + t18 ^ 2, t22, t8, 0, t21, 0, 0, t23 * t60, t25 * t60, 0.2e1 * t57, t13 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3 * t23 + t4 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 + t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t19, -t20, 0, 0, t36, t28, 0, t37, 0, 0, t52 * t40, -t52 * t38, t48 + t49, -t18 * pkin(3) + t17 * t48, t22, t8, 0, t21, 0, 0, t50 * t23, t50 * t25, t56 + t57, t4 * t10 + t13 * t32 + t3 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 * t10 - t23 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t36, t28, 0, t37, 0, 0, pkin(3) * t58, -0.2e1 * pkin(3) * t38, 0.2e1 * t48, t47 * qJ(4) ^ 2 + pkin(3) ^ 2, t22, t8, 0, t21, 0, 0, t23 * t59, t25 * t59, 0.2e1 * t56, t10 ^ 2 + t32 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t38, 0, t18, 0, 0, 0, 0, 0, 0, t23, t25, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t38, 0, -pkin(3), 0, 0, 0, 0, 0, 0, t23, t25, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, -t23, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t25, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, -t23, 0, t9, -t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
