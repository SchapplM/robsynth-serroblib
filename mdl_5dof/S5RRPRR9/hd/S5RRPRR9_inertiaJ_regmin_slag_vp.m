% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t38 = sin(pkin(9));
t39 = cos(pkin(9));
t42 = sin(qJ(2));
t56 = cos(qJ(2));
t24 = t38 * t42 - t39 * t56;
t64 = -0.2e1 * t24;
t63 = 0.2e1 * t24;
t33 = -t39 * pkin(2) - pkin(3);
t44 = cos(qJ(4));
t29 = -t44 * pkin(4) + t33;
t62 = 0.2e1 * t29;
t61 = t24 * pkin(4);
t40 = sin(qJ(5));
t60 = t40 * pkin(4);
t43 = cos(qJ(5));
t59 = t43 * pkin(4);
t25 = t38 * t56 + t39 * t42;
t35 = -t56 * pkin(2) - pkin(1);
t14 = t24 * pkin(3) - t25 * pkin(7) + t35;
t41 = sin(qJ(4));
t48 = t56 * pkin(6);
t30 = t56 * qJ(3) + t48;
t47 = (-qJ(3) - pkin(6)) * t42;
t17 = t39 * t30 + t38 * t47;
t51 = t44 * t17;
t5 = t51 + (-pkin(8) * t25 + t14) * t41;
t58 = t43 * t5;
t32 = t38 * pkin(2) + pkin(7);
t57 = pkin(8) + t32;
t28 = t40 * t44 + t43 * t41;
t55 = t28 * t24;
t54 = t41 * t24;
t53 = t41 * t25;
t52 = t41 * t44;
t50 = t44 * t25;
t49 = 0.2e1 * t56;
t6 = t44 * t14 - t41 * t17;
t4 = -pkin(8) * t50 + t6 + t61;
t1 = t43 * t4 - t40 * t5;
t15 = t38 * t30 - t39 * t47;
t46 = -t24 * t32 + t25 * t33;
t27 = t40 * t41 - t43 * t44;
t37 = t44 ^ 2;
t36 = t41 ^ 2;
t23 = t25 ^ 2;
t22 = t24 ^ 2;
t21 = t57 * t44;
t20 = t57 * t41;
t19 = t44 * t24;
t18 = t27 * t24;
t13 = -t40 * t20 + t43 * t21;
t12 = -t43 * t20 - t40 * t21;
t10 = t27 * t25;
t9 = t28 * t25;
t8 = pkin(4) * t53 + t15;
t7 = t41 * t14 + t51;
t2 = t40 * t4 + t58;
t3 = [1, 0, 0, t42 ^ 2, t42 * t49, 0, 0, 0, pkin(1) * t49, -0.2e1 * pkin(1) * t42, 0.2e1 * t15 * t25 - 0.2e1 * t17 * t24, t15 ^ 2 + t17 ^ 2 + t35 ^ 2, t37 * t23, -0.2e1 * t23 * t52, t50 * t63, t53 * t64, t22, 0.2e1 * t15 * t53 + 0.2e1 * t6 * t24, 0.2e1 * t15 * t50 - 0.2e1 * t7 * t24, t10 ^ 2, 0.2e1 * t10 * t9, -t10 * t63, t9 * t64, t22, 0.2e1 * t1 * t24 + 0.2e1 * t8 * t9, -0.2e1 * t8 * t10 - 0.2e1 * t2 * t24; 0, 0, 0, 0, 0, t42, t56, 0, -t42 * pkin(6), -t48, (-t24 * t38 - t25 * t39) * pkin(2), (-t15 * t39 + t17 * t38) * pkin(2), t41 * t50, (-t36 + t37) * t25, t54, t19, 0, -t15 * t44 + t46 * t41, t15 * t41 + t46 * t44, -t10 * t28, t10 * t27 - t28 * t9, t55, -t18, 0, t12 * t24 + t8 * t27 + t29 * t9, -t29 * t10 - t13 * t24 + t8 * t28; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t38 ^ 2 + t39 ^ 2) * pkin(2) ^ 2, t36, 0.2e1 * t52, 0, 0, 0, -0.2e1 * t33 * t44, 0.2e1 * t33 * t41, t28 ^ 2, -0.2e1 * t28 * t27, 0, 0, 0, t27 * t62, t28 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, t19, -t54, 0, 0, 0, 0, 0, -t18, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t53, t24, t6, -t7, 0, 0, -t10, -t9, t24, t24 * t59 + t1, -t58 + (-t4 - t61) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t44, 0, -t41 * t32, -t44 * t32, 0, 0, t28, -t27, 0, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t41, 0, 0, 0, 0, 0, -t27, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t59, -0.2e1 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t9, t24, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t27, 0, t12, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t59, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
