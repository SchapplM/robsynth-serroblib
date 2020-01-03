% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRRR5
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
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRR5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t33 = cos(qJ(3));
t23 = -t33 * pkin(3) - pkin(2);
t60 = 0.2e1 * t23;
t31 = sin(qJ(2));
t59 = -0.2e1 * t31;
t34 = cos(qJ(2));
t58 = -0.2e1 * t34;
t57 = 0.2e1 * t34;
t56 = -pkin(7) - pkin(6);
t55 = pkin(2) * t33;
t30 = sin(qJ(3));
t54 = pkin(5) * t30;
t26 = t31 ^ 2;
t53 = t26 * pkin(5);
t29 = sin(qJ(4));
t52 = t29 * pkin(3);
t51 = t31 * pkin(5);
t32 = cos(qJ(4));
t50 = t32 * pkin(3);
t18 = -t34 * pkin(2) - t31 * pkin(6) - pkin(1);
t43 = t33 * t34;
t40 = pkin(5) * t43;
t5 = t40 + (-pkin(7) * t31 + t18) * t30;
t49 = t32 * t5;
t48 = t34 * pkin(3);
t47 = t30 * t31;
t46 = t30 * t33;
t45 = t30 * t34;
t44 = t33 * t31;
t25 = t30 ^ 2;
t27 = t33 ^ 2;
t42 = t25 + t27;
t41 = t31 * t57;
t39 = t30 * t44;
t13 = t33 * t18;
t4 = -pkin(7) * t44 + t13 + (-pkin(3) - t54) * t34;
t1 = -t29 * t5 + t32 * t4;
t8 = -pkin(5) * t45 + t13;
t9 = t30 * t18 + t40;
t38 = -t8 * t30 + t9 * t33;
t16 = t29 * t33 + t32 * t30;
t36 = pkin(5) ^ 2;
t28 = t34 ^ 2;
t24 = t26 * t36;
t20 = t56 * t33;
t19 = t56 * t30;
t17 = (pkin(3) * t30 + pkin(5)) * t31;
t14 = t29 * t30 - t32 * t33;
t12 = -t29 * t47 + t32 * t44;
t10 = t16 * t31;
t7 = t29 * t19 - t32 * t20;
t6 = t32 * t19 + t29 * t20;
t2 = t29 * t4 + t49;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t26, t41, 0, t28, 0, 0, pkin(1) * t57, pkin(1) * t59, 0.2e1 * (t26 + t28) * pkin(5), pkin(1) ^ 2 + t28 * t36 + t24, t27 * t26, -0.2e1 * t26 * t46, t43 * t59, t25 * t26, t30 * t41, t28, 0.2e1 * t30 * t53 - 0.2e1 * t8 * t34, 0.2e1 * t33 * t53 + 0.2e1 * t9 * t34, 0.2e1 * (-t30 * t9 - t33 * t8) * t31, t8 ^ 2 + t9 ^ 2 + t24, t12 ^ 2, -0.2e1 * t12 * t10, t12 * t58, t10 ^ 2, -t10 * t58, t28, -0.2e1 * t1 * t34 + 0.2e1 * t17 * t10, 0.2e1 * t17 * t12 + 0.2e1 * t2 * t34, -0.2e1 * t1 * t12 - 0.2e1 * t2 * t10, t1 ^ 2 + t17 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t34, 0, -t51, -t34 * pkin(5), 0, 0, t39, (-t25 + t27) * t31, -t45, -t39, -t43, 0, -pkin(5) * t44 + (-pkin(2) * t31 + pkin(6) * t34) * t30, pkin(6) * t43 + (t54 - t55) * t31, t38, -pkin(2) * t51 + t38 * pkin(6), t12 * t16, -t16 * t10 - t12 * t14, -t16 * t34, t10 * t14, t14 * t34, 0, t23 * t10 + t17 * t14 - t6 * t34, t23 * t12 + t17 * t16 + t7 * t34, -t1 * t16 - t7 * t10 - t6 * t12 - t2 * t14, t1 * t6 + t17 * t23 + t2 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t25, 0.2e1 * t46, 0, t27, 0, 0, 0.2e1 * t55, -0.2e1 * pkin(2) * t30, 0.2e1 * t42 * pkin(6), t42 * pkin(6) ^ 2 + pkin(2) ^ 2, t16 ^ 2, -0.2e1 * t16 * t14, 0, t14 ^ 2, 0, 0, t14 * t60, t16 * t60, -0.2e1 * t7 * t14 - 0.2e1 * t6 * t16, t23 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, -t47, -t34, t8, -t9, 0, 0, 0, 0, t12, 0, -t10, -t34, -t32 * t48 + t1, -t49 + (-t4 + t48) * t29, (-t10 * t29 - t12 * t32) * pkin(3), (t1 * t32 + t2 * t29) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, t33, 0, -t30 * pkin(6), -t33 * pkin(6), 0, 0, 0, 0, t16, 0, -t14, 0, t6, -t7, (-t14 * t29 - t16 * t32) * pkin(3), (t29 * t7 + t32 * t6) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t50, -0.2e1 * t52, 0, (t29 ^ 2 + t32 ^ 2) * pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, -t10, -t34, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t14, 0, t6, -t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t50, -t52, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t3;
