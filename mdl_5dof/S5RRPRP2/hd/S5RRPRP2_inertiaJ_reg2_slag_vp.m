% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t36 = sin(qJ(4));
t32 = t36 ^ 2;
t38 = cos(qJ(4));
t33 = t38 ^ 2;
t23 = t32 + t33;
t34 = sin(pkin(8));
t58 = t34 * pkin(2);
t28 = pkin(7) + t58;
t46 = t23 * t28;
t63 = -0.2e1 * t36;
t62 = 0.2e1 * t36;
t61 = -0.2e1 * t38;
t39 = cos(qJ(2));
t31 = t39 * pkin(1);
t30 = t31 + pkin(2);
t35 = cos(pkin(8));
t37 = sin(qJ(2));
t56 = t37 * pkin(1);
t44 = t35 * t56;
t15 = t34 * t30 + t44;
t13 = pkin(7) + t15;
t60 = t46 * t13;
t59 = t23 * t13 ^ 2;
t57 = t35 * pkin(2);
t55 = t23 * t13;
t43 = t38 * pkin(4) + t36 * qJ(5);
t42 = -pkin(3) - t43;
t16 = t42 - t57;
t45 = -t35 * t30 + t34 * t56;
t5 = t42 + t45;
t54 = -t16 - t5;
t53 = t36 * t13;
t52 = t36 * t28;
t51 = t36 * t38;
t50 = t38 * t13;
t49 = t38 * t28;
t12 = -pkin(3) + t45;
t29 = -pkin(3) - t57;
t48 = t12 + t29;
t47 = t23 * t28 ^ 2;
t17 = -t36 * pkin(4) + t38 * qJ(5);
t27 = -0.2e1 * t51;
t26 = 0.2e1 * t51;
t6 = 0.2e1 * t46;
t2 = 0.2e1 * t55;
t1 = t46 + t55;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t31, -0.2e1 * t56, 0, (t37 ^ 2 + t39 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t45, -0.2e1 * t15, 0, t15 ^ 2 + t45 ^ 2, t32, t26, 0, t33, 0, 0, t12 * t61, t12 * t62, t2, t12 ^ 2 + t59, t32, 0, t27, 0, 0, t33, t5 * t61, t2, t5 * t63, t5 ^ 2 + t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t31, -t56, 0, 0, 0, 0, 0, 0, 0, 1, -t45 + t57, -t44 + (-pkin(2) - t30) * t34, 0, (t15 * t34 - t35 * t45) * pkin(2), t32, t26, 0, t33, 0, 0, -t48 * t38, t48 * t36, t1, t12 * t29 + t60, t32, 0, t27, 0, 0, t33, t54 * t38, t1, t54 * t36, t5 * t16 + t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t57, -0.2e1 * t58, 0, (t34 ^ 2 + t35 ^ 2) * pkin(2) ^ 2, t32, t26, 0, t33, 0, 0, t29 * t61, t29 * t62, t6, t29 ^ 2 + t47, t32, 0, t27, 0, 0, t33, t16 * t61, t6, t16 * t63, t16 ^ 2 + t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, t38, 0, -t53, -t50, 0, 0, 0, t36, 0, 0, -t38, 0, -t53, t17, t50, t17 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, t38, 0, -t52, -t49, 0, 0, 0, t36, 0, 0, -t38, 0, -t52, t17, t49, t17 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t36, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t36, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t3;
