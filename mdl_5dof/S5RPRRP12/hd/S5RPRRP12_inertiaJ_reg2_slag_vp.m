% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP12_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:20
% EndTime: 2019-12-31 18:57:22
% DurationCPUTime: 0.65s
% Computational Cost: add. (273->83), mult. (536->141), div. (0->0), fcn. (491->4), ass. (0->63)
t63 = 2 * pkin(4);
t36 = cos(qJ(4));
t62 = 0.2e1 * t36;
t37 = cos(qJ(3));
t61 = 0.2e1 * t37;
t60 = 2 * qJ(2);
t34 = sin(qJ(4));
t59 = t34 * pkin(4);
t58 = t37 * pkin(3);
t35 = sin(qJ(3));
t23 = t34 * t35;
t57 = t34 * t36;
t56 = t34 * t37;
t38 = -pkin(1) - pkin(6);
t55 = t34 * t38;
t54 = t35 * t38;
t26 = t36 * t37;
t53 = t36 * t38;
t27 = -t36 * pkin(4) - pkin(3);
t52 = t37 * t27;
t51 = t37 * t35;
t50 = t37 * t38;
t49 = -qJ(5) - pkin(7);
t29 = t34 ^ 2;
t31 = t36 ^ 2;
t48 = t29 + t31;
t30 = t35 ^ 2;
t32 = t37 ^ 2;
t20 = t30 + t32;
t47 = qJ(5) * t37;
t46 = -0.2e1 * t51;
t45 = t35 * t53;
t10 = t48 * t35;
t13 = t35 * pkin(3) - t37 * pkin(7) + qJ(2);
t6 = t36 * t13;
t44 = -t36 * t47 + t6;
t43 = -pkin(7) * t35 - t58;
t3 = -t34 * t54 + t6;
t4 = t34 * t13 + t45;
t42 = -t3 * t34 + t4 * t36;
t14 = t49 * t34;
t15 = t49 * t36;
t41 = -t14 * t34 - t15 * t36;
t39 = qJ(2) ^ 2;
t33 = t38 ^ 2;
t28 = t32 * t33;
t25 = t36 * t35;
t24 = t31 * t32;
t22 = t29 * t32;
t21 = 0.2e1 * t57;
t19 = t34 * t26;
t18 = t51 * t62;
t17 = -0.2e1 * t32 * t57;
t16 = t34 * t46;
t12 = t20 * t38;
t11 = t20 * t36;
t9 = t20 * t34;
t8 = (-t29 + t31) * t37;
t7 = (-t38 + t59) * t37;
t5 = t48 * t30 + t32;
t2 = t45 + (t13 - t47) * t34;
t1 = (pkin(4) - t55) * t35 + t44;
t40 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t60, pkin(1) ^ 2 + t39, t32, t46, 0, t30, 0, 0, t35 * t60, t37 * t60, -0.2e1 * t12, t30 * t33 + t28 + t39, t24, t17, t18, t22, t16, t30, 0.2e1 * t3 * t35 - 0.2e1 * t32 * t55, -0.2e1 * t32 * t53 - 0.2e1 * t4 * t35, (-t3 * t36 - t34 * t4) * t61, t3 ^ 2 + t4 ^ 2 + t28, t24, t17, t18, t22, t16, t30, 0.2e1 * t1 * t35 + 0.2e1 * t7 * t56, -0.2e1 * t2 * t35 + 0.2e1 * t7 * t26, (-t1 * t36 - t2 * t34) * t61, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t20, t12, 0, 0, 0, 0, 0, 0, -t9, -t11, 0, t32 * t38 + t42 * t35, 0, 0, 0, 0, 0, 0, -t9, -t11, 0, -t7 * t37 + (-t1 * t34 + t2 * t36) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, -t35, 0, t50, -t54, 0, 0, t19, t8, t23, -t19, t25, 0, t43 * t34 + t36 * t50, -t34 * t50 + t43 * t36, t42, pkin(3) * t50 + t42 * pkin(7), t19, t8, t23, -t19, t25, 0, t14 * t35 + t34 * t52 - t7 * t36, t15 * t35 + t7 * t34 + t36 * t52, (-t14 * t37 + t2) * t36 + (t15 * t37 - t1) * t34, t1 * t14 - t2 * t15 + t7 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t35, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t56, t10, pkin(7) * t10 + t58, 0, 0, 0, 0, 0, 0, t26, -t56, t10, t41 * t35 - t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t29, t21, 0, t31, 0, 0, pkin(3) * t62, -0.2e1 * pkin(3) * t34, 0.2e1 * t48 * pkin(7), t48 * pkin(7) ^ 2 + pkin(3) ^ 2, t29, t21, 0, t31, 0, 0, -0.2e1 * t27 * t36, 0.2e1 * t27 * t34, 0.2e1 * t41, t14 ^ 2 + t15 ^ 2 + t27 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, -t56, t35, t3, -t4, 0, 0, 0, 0, t26, 0, -t56, t35, (t63 - t55) * t35 + t44, -t2, -pkin(4) * t26, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t25, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t25, 0, -pkin(4) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, t36, 0, -t34 * pkin(7), -t36 * pkin(7), 0, 0, 0, 0, t34, 0, t36, 0, t14, t15, -t59, t14 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t63, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t26, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t34, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t40;
