% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPR7
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPR7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t33 = sin(pkin(9));
t35 = cos(pkin(9));
t38 = sin(qJ(3));
t40 = cos(qJ(3));
t20 = t33 * t40 + t35 * t38;
t66 = -0.2e1 * t20;
t34 = sin(pkin(8));
t60 = t34 * pkin(1);
t26 = pkin(6) + t60;
t50 = qJ(4) + t26;
t14 = t50 * t40;
t46 = t50 * t38;
t4 = t33 * t14 + t35 * t46;
t65 = t4 ^ 2;
t18 = t33 * t38 - t35 * t40;
t64 = t18 ^ 2;
t36 = cos(pkin(8));
t58 = t36 * pkin(1);
t28 = -pkin(2) - t58;
t21 = -t40 * pkin(3) + t28;
t63 = 0.2e1 * t21;
t62 = 0.2e1 * t38;
t61 = t33 * pkin(3);
t59 = t35 * pkin(3);
t57 = t4 * t18;
t37 = sin(qJ(5));
t29 = t37 ^ 2;
t56 = t29 * t20;
t10 = t37 * t18;
t55 = t37 * t20;
t39 = cos(qJ(5));
t54 = t37 * t39;
t53 = t39 * t20;
t31 = t39 ^ 2;
t52 = t29 + t31;
t30 = t38 ^ 2;
t32 = t40 ^ 2;
t51 = t30 + t32;
t49 = t18 * t66;
t48 = t37 * t53;
t25 = pkin(7) + t61;
t47 = t52 * t25;
t3 = t18 * pkin(4) - t20 * pkin(7) + t21;
t6 = t35 * t14 - t33 * t46;
t1 = t39 * t3 - t37 * t6;
t2 = t37 * t3 + t39 * t6;
t45 = t1 * t39 + t2 * t37;
t44 = -t1 * t37 + t2 * t39;
t27 = -pkin(4) - t59;
t43 = -t18 * t25 + t20 * t27;
t17 = t20 ^ 2;
t13 = t39 * t18;
t12 = t31 * t20;
t11 = t31 * t17;
t9 = t29 * t17;
t7 = -t12 - t56;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t58, -0.2e1 * t60, 0, (t34 ^ 2 + t36 ^ 2) * pkin(1) ^ 2, t30, t40 * t62, 0, t32, 0, 0, -0.2e1 * t28 * t40, t28 * t62, 0.2e1 * t51 * t26, t51 * t26 ^ 2 + t28 ^ 2, t17, t49, 0, t64, 0, 0, t18 * t63, t20 * t63, -0.2e1 * t6 * t18 + 0.2e1 * t4 * t20, t21 ^ 2 + t6 ^ 2 + t65, t11, -0.2e1 * t17 * t54, 0.2e1 * t18 * t53, t9, t37 * t49, t64, 0.2e1 * t1 * t18 + 0.2e1 * t4 * t55, -0.2e1 * t2 * t18 + 0.2e1 * t4 * t53, t45 * t66, t1 ^ 2 + t2 ^ 2 + t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 * t20 + t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t20 + t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 + t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 + t9 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t40, 0, -t38 * t26, -t40 * t26, 0, 0, 0, 0, t20, 0, -t18, 0, -t4, -t6, (-t18 * t33 - t20 * t35) * pkin(3), (t33 * t6 - t35 * t4) * pkin(3), t48, t12 - t56, t10, -t48, t13, 0, t43 * t37 - t4 * t39, t4 * t37 + t43 * t39, t44, t44 * t25 + t4 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t38, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t20, 0, (-t18 * t35 + t20 * t33) * pkin(3), 0, 0, 0, 0, 0, 0, -t13, t10, -t7, t18 * t27 + t20 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t59, -0.2e1 * t61, 0, (t33 ^ 2 + t35 ^ 2) * pkin(3) ^ 2, t29, 0.2e1 * t54, 0, t31, 0, 0, -0.2e1 * t27 * t39, 0.2e1 * t27 * t37, 0.2e1 * t47, t52 * t25 ^ 2 + t27 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t20, 0, t21, 0, 0, 0, 0, 0, 0, t13, -t10, t7, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, -t55, t18, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t53, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, t39, 0, -t37 * t25, -t39 * t25, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t37, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t5;
