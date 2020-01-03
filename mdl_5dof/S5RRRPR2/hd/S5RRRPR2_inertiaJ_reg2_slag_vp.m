% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t43 = cos(qJ(3));
t33 = t43 * pkin(2);
t30 = t33 + pkin(3);
t37 = sin(pkin(9));
t38 = cos(pkin(9));
t40 = sin(qJ(3));
t59 = t40 * pkin(2);
t53 = t38 * t59;
t18 = t37 * t30 + t53;
t15 = pkin(8) + t18;
t39 = sin(qJ(5));
t35 = t39 ^ 2;
t42 = cos(qJ(5));
t36 = t42 ^ 2;
t54 = t35 + t36;
t57 = t15 * t54;
t60 = t37 * pkin(3);
t28 = pkin(8) + t60;
t63 = t54 * t28;
t62 = -0.2e1 * t42;
t44 = cos(qJ(2));
t34 = t44 * pkin(1);
t31 = t34 + pkin(2);
t41 = sin(qJ(2));
t58 = t41 * pkin(1);
t19 = t43 * t31 - t40 * t58;
t16 = pkin(3) + t19;
t52 = t43 * t58;
t20 = t40 * t31 + t52;
t56 = t38 * t20;
t8 = t37 * t16 + t56;
t6 = pkin(8) + t8;
t61 = t54 * t6;
t51 = -t20 - t59;
t12 = t38 * t16;
t49 = t37 * t20 - t12;
t22 = t38 * t30;
t48 = t37 * t59 - t22;
t32 = t38 * pkin(3);
t29 = -t32 - pkin(4);
t27 = 0.2e1 * t39 * t42;
t24 = t29 * t39;
t14 = -pkin(4) + t48;
t11 = t14 * t39;
t5 = -pkin(4) + t49;
t3 = t5 * t39;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t34, -0.2e1 * t58, 0, (t41 ^ 2 + t44 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t19, -0.2e1 * t20, 0, t19 ^ 2 + t20 ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t49, -0.2e1 * t8, 0, t49 ^ 2 + t8 ^ 2, t35, t27, 0, t36, 0, 0, t5 * t62, 0.2e1 * t3, 0.2e1 * t61, t54 * t6 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t34, -t58, 0, 0, 0, 0, 0, 0, 0, 1, t19 + t33, -t52 + (-pkin(2) - t31) * t40, 0, (t19 * t43 + t20 * t40) * pkin(2), 0, 0, 0, 0, 0, 1, t51 * t37 + t12 + t22, t51 * t38 + (-t16 - t30) * t37, 0, t8 * t18 + t48 * t49, t35, t27, 0, t36, 0, 0, (-t14 - t5) * t42, t11 + t3, t57 + t61, t5 * t14 + t6 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t33, -0.2e1 * t59, 0, (t40 ^ 2 + t43 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t48, -0.2e1 * t18, 0, t18 ^ 2 + t48 ^ 2, t35, t27, 0, t36, 0, 0, t14 * t62, 0.2e1 * t11, 0.2e1 * t57, t54 * t15 ^ 2 + t14 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t19, -t20, 0, 0, 0, 0, 0, 0, 0, 1, t32 - t49, -t56 + (-pkin(3) - t16) * t37, 0, (t37 * t8 - t38 * t49) * pkin(3), t35, t27, 0, t36, 0, 0, (-t29 - t5) * t42, t24 + t3, t63 + t61, t5 * t29 + t6 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t33, -t59, 0, 0, 0, 0, 0, 0, 0, 1, t32 - t48, -t53 + (-pkin(3) - t30) * t37, 0, (t18 * t37 - t38 * t48) * pkin(3), t35, t27, 0, t36, 0, 0, (-t14 - t29) * t42, t24 + t11, t63 + t57, t14 * t29 + t15 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t32, -0.2e1 * t60, 0, (t37 ^ 2 + t38 ^ 2) * pkin(3) ^ 2, t35, t27, 0, t36, 0, 0, t29 * t62, 0.2e1 * t24, 0.2e1 * t63, t54 * t28 ^ 2 + t29 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, t42, 0, -t39 * t6, -t42 * t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, t42, 0, -t39 * t15, -t42 * t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, t42, 0, -t39 * t28, -t42 * t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t39, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
