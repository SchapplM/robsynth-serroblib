% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:16
% EndTime: 2019-12-31 18:43:18
% DurationCPUTime: 0.61s
% Computational Cost: add. (302->86), mult. (602->143), div. (0->0), fcn. (558->6), ass. (0->62)
t35 = sin(qJ(3));
t63 = 0.2e1 * t35;
t36 = cos(qJ(4));
t62 = pkin(3) * t36;
t32 = sin(pkin(8));
t61 = t32 * pkin(1);
t33 = cos(pkin(8));
t60 = t33 * pkin(1);
t34 = sin(qJ(4));
t59 = t34 * pkin(4);
t37 = cos(qJ(3));
t58 = t37 * pkin(3);
t20 = pkin(6) + t61;
t57 = t20 * t34;
t28 = t34 ^ 2;
t56 = t28 * t35;
t29 = t35 ^ 2;
t55 = t29 * t20;
t54 = t34 * t35;
t53 = t34 * t36;
t52 = t34 * t37;
t51 = t35 * t20;
t25 = t36 * t35;
t26 = t36 * t37;
t50 = t37 * t20;
t49 = -qJ(5) - pkin(7);
t30 = t36 ^ 2;
t48 = t28 + t30;
t31 = t37 ^ 2;
t47 = t29 + t31;
t46 = qJ(5) * t35;
t45 = t37 * t63;
t44 = t36 * t50;
t21 = -pkin(2) - t60;
t43 = t48 * pkin(7);
t8 = -pkin(7) * t35 + t21 - t58;
t5 = t36 * t8;
t42 = -t36 * t46 + t5;
t3 = -t34 * t50 + t5;
t4 = t34 * t8 + t44;
t41 = -t3 * t34 + t36 * t4;
t11 = t49 * t34;
t12 = t49 * t36;
t40 = -t11 * t34 - t12 * t36;
t27 = -pkin(4) * t36 - pkin(3);
t24 = t30 * t35;
t23 = t30 * t29;
t22 = t28 * t29;
t19 = 0.2e1 * t53;
t18 = t20 ^ 2;
t17 = t34 * t25;
t16 = -0.2e1 * t35 * t26;
t15 = -0.2e1 * t29 * t53;
t14 = t34 * t45;
t13 = t29 * t18;
t10 = t24 + t56;
t9 = t24 - t56;
t7 = t23 + t22 + t31;
t6 = (t20 + t59) * t35;
t2 = t44 + (t8 - t46) * t34;
t1 = (-pkin(4) - t57) * t37 + t42;
t38 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t60, -0.2e1 * t61, 0, (t32 ^ 2 + t33 ^ 2) * pkin(1) ^ 2, t29, t45, 0, t31, 0, 0, -0.2e1 * t21 * t37, t21 * t63, 0.2e1 * t47 * t20, t18 * t31 + t21 ^ 2 + t13, t23, t15, t16, t22, t14, t31, -0.2e1 * t3 * t37 + 0.2e1 * t34 * t55, 0.2e1 * t36 * t55 + 0.2e1 * t37 * t4, (-t3 * t36 - t34 * t4) * t63, t3 ^ 2 + t4 ^ 2 + t13, t23, t15, t16, t22, t14, t31, -0.2e1 * t1 * t37 + 0.2e1 * t54 * t6, 0.2e1 * t2 * t37 + 0.2e1 * t25 * t6, (-t1 * t36 - t2 * t34) * t63, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t41 - t50) * t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6 * t37 + (-t1 * t34 + t2 * t36) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, t37, 0, -t51, -t50, 0, 0, t17, t9, -t52, -t17, -t26, 0, -t20 * t25 + (-pkin(3) * t35 + pkin(7) * t37) * t34, pkin(7) * t26 + (t57 - t62) * t35, t41, -pkin(3) * t51 + pkin(7) * t41, t17, t9, -t52, -t17, -t26, 0, -t11 * t37 + t27 * t54 - t36 * t6, -t12 * t37 + t25 * t27 + t34 * t6, (-t11 * t35 + t2) * t36 + (t12 * t35 - t1) * t34, t1 * t11 - t12 * t2 + t27 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t35, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t52, t10, t35 * t43 + t58, 0, 0, 0, 0, 0, 0, t26, -t52, t10, -t37 * t27 + t35 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t28, t19, 0, t30, 0, 0, 0.2e1 * t62, -0.2e1 * pkin(3) * t34, 0.2e1 * t43, pkin(7) ^ 2 * t48 + pkin(3) ^ 2, t28, t19, 0, t30, 0, 0, -0.2e1 * t27 * t36, 0.2e1 * t27 * t34, 0.2e1 * t40, t11 ^ 2 + t12 ^ 2 + t27 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, -t54, -t37, t3, -t4, 0, 0, 0, 0, t25, 0, -t54, -t37, (-0.2e1 * pkin(4) - t57) * t37 + t42, -t2, -pkin(4) * t25, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t25, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t25, 0, -pkin(4) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, t36, 0, -t34 * pkin(7), -t36 * pkin(7), 0, 0, 0, 0, t34, 0, t36, 0, t11, t12, -t59, t11 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t25, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, t34, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t38;
