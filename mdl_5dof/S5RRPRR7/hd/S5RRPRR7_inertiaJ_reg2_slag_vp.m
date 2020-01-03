% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t39 = sin(qJ(5));
t40 = sin(qJ(4));
t42 = cos(qJ(5));
t43 = cos(qJ(4));
t17 = t39 * t43 + t42 * t40;
t19 = -t39 * t40 + t42 * t43;
t8 = (t17 * t39 + t19 * t42) * pkin(4);
t15 = t17 ^ 2;
t16 = t19 ^ 2;
t66 = t16 + t15;
t41 = sin(qJ(2));
t58 = t41 * pkin(1);
t29 = qJ(3) + t58;
t65 = t29 ^ 2;
t35 = t40 * pkin(4);
t23 = t29 + t35;
t64 = 0.2e1 * t23;
t63 = 0.2e1 * t29;
t31 = qJ(3) + t35;
t62 = 0.2e1 * t31;
t46 = 0.2e1 * qJ(3);
t44 = cos(qJ(2));
t55 = t44 * pkin(1);
t33 = -pkin(2) - t55;
t27 = -pkin(7) + t33;
t13 = (-pkin(8) + t27) * t40;
t24 = t43 * t27;
t56 = t43 * pkin(8);
t14 = t24 - t56;
t6 = -t39 * t13 + t42 * t14;
t7 = t42 * t13 + t39 * t14;
t61 = -t7 * t17 - t6 * t19;
t45 = -pkin(2) - pkin(7);
t21 = (-pkin(8) + t45) * t40;
t34 = t43 * t45;
t22 = t34 - t56;
t10 = t42 * t21 + t39 * t22;
t9 = -t39 * t21 + t42 * t22;
t60 = -t10 * t17 - t9 * t19;
t59 = t39 * pkin(4);
t57 = t42 * pkin(4);
t54 = t23 + t31;
t36 = t40 ^ 2;
t37 = t43 ^ 2;
t25 = t36 + t37;
t53 = t29 * qJ(3);
t52 = qJ(3) + t29;
t20 = t25 * t45;
t48 = -0.2e1 * pkin(2);
t47 = qJ(3) ^ 2;
t28 = -0.2e1 * t43 * t40;
t12 = t25 * t27;
t11 = -0.2e1 * t19 * t17;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t55, -0.2e1 * t58, 0, (t41 ^ 2 + t44 ^ 2) * pkin(1) ^ 2, 1, 0, 0, 0, 0, 0, 0, 0.2e1 * t33, t63, t33 ^ 2 + t65, t37, t28, 0, t36, 0, 0, t40 * t63, t43 * t63, -0.2e1 * t12, t25 * t27 ^ 2 + t65, t16, t11, 0, t15, 0, 0, t17 * t64, t19 * t64, 0.2e1 * t61, t23 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t55, -t58, 0, 0, 1, 0, 0, 0, 0, 0, 0, t48 - t55, t46 + t58, -t33 * pkin(2) + t53, t37, t28, 0, t36, 0, 0, t52 * t40, t52 * t43, (-t27 - t45) * t25, t27 * t20 + t53, t16, t11, 0, t15, 0, 0, t54 * t17, t54 * t19, t60 + t61, t7 * t10 + t23 * t31 + t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, t48, t46, pkin(2) ^ 2 + t47, t37, t28, 0, t36, 0, 0, t40 * t46, t43 * t46, -0.2e1 * t20, t25 * t45 ^ 2 + t47, t16, t11, 0, t15, 0, 0, t17 * t62, t19 * t62, 0.2e1 * t60, t10 ^ 2 + t31 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t12, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t25, t20, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, -t40, 0, t24, -t40 * t27, 0, 0, 0, 0, t19, 0, -t17, 0, t6, -t7, -t8, (t39 * t7 + t42 * t6) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, -t40, 0, t34, -t40 * t45, 0, 0, 0, 0, t19, 0, -t17, 0, t9, -t10, -t8, (t10 * t39 + t42 * t9) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t40, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t17, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t57, -0.2e1 * t59, 0, (t39 ^ 2 + t42 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, -t17, 0, t6, -t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, -t17, 0, t9, -t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t17, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t57, -t59, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
