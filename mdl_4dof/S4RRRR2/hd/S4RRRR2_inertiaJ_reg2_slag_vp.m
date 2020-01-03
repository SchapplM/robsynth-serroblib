% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t36 = sin(qJ(2));
t52 = t36 * pkin(1);
t26 = pkin(6) + t52;
t35 = sin(qJ(3));
t32 = t35 ^ 2;
t38 = cos(qJ(3));
t33 = t38 ^ 2;
t43 = t32 + t33;
t45 = t43 * t26;
t28 = -pkin(3) * t38 - pkin(2);
t39 = cos(qJ(2));
t49 = t39 * pkin(1);
t18 = t28 - t49;
t58 = 0.2e1 * t18;
t57 = 0.2e1 * t28;
t56 = 0.2e1 * t38;
t34 = sin(qJ(4));
t37 = cos(qJ(4));
t15 = t34 * t35 - t37 * t38;
t17 = t34 * t38 + t35 * t37;
t11 = (-pkin(7) - t26) * t35;
t31 = t38 * pkin(7);
t47 = t38 * t26;
t12 = t31 + t47;
t5 = t11 * t37 - t12 * t34;
t6 = t11 * t34 + t12 * t37;
t55 = -t6 * t15 - t5 * t17;
t19 = (-pkin(7) - pkin(6)) * t35;
t50 = t38 * pkin(6);
t20 = t31 + t50;
t10 = t19 * t34 + t20 * t37;
t9 = t19 * t37 - t20 * t34;
t54 = -t10 * t15 - t9 * t17;
t53 = t34 * pkin(3);
t51 = t37 * pkin(3);
t27 = -pkin(2) - t49;
t48 = pkin(2) - t27;
t46 = t18 + t28;
t44 = t43 * pkin(6);
t23 = t35 * t56;
t14 = t17 ^ 2;
t13 = t15 ^ 2;
t8 = -0.2e1 * t17 * t15;
t7 = (-t15 * t34 - t17 * t37) * pkin(3);
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t49, -0.2e1 * t52, 0, (t36 ^ 2 + t39 ^ 2) * pkin(1) ^ 2, t32, t23, 0, t33, 0, 0, -0.2e1 * t27 * t38, 0.2e1 * t27 * t35, 0.2e1 * t45, t43 * t26 ^ 2 + t27 ^ 2, t14, t8, 0, t13, 0, 0, t15 * t58, t17 * t58, 0.2e1 * t55, t18 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t49, -t52, 0, 0, t32, t23, 0, t33, 0, 0, t48 * t38, -t48 * t35, t44 + t45, -t27 * pkin(2) + pkin(6) * t45, t14, t8, 0, t13, 0, 0, t46 * t15, t46 * t17, t54 + t55, t10 * t6 + t18 * t28 + t5 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t32, t23, 0, t33, 0, 0, pkin(2) * t56, -0.2e1 * pkin(2) * t35, 0.2e1 * t44, t43 * pkin(6) ^ 2 + pkin(2) ^ 2, t14, t8, 0, t13, 0, 0, t15 * t57, t17 * t57, 0.2e1 * t54, t10 ^ 2 + t28 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, t38, 0, -t35 * t26, -t47, 0, 0, 0, 0, t17, 0, -t15, 0, t5, -t6, t7, (t34 * t6 + t37 * t5) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, t38, 0, -t35 * pkin(6), -t50, 0, 0, 0, 0, t17, 0, -t15, 0, t9, -t10, t7, (t10 * t34 + t37 * t9) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t51, -0.2e1 * t53, 0, (t34 ^ 2 + t37 ^ 2) * pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, -t15, 0, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, -t15, 0, t9, -t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t51, -t53, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
