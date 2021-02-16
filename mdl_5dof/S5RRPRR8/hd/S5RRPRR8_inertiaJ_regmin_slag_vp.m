% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRR8
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
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:36
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:35:11
% EndTime: 2021-01-15 21:35:14
% DurationCPUTime: 0.36s
% Computational Cost: add. (445->56), mult. (900->120), div. (0->0), fcn. (1095->8), ass. (0->54)
t34 = sin(pkin(9));
t35 = cos(pkin(9));
t38 = sin(qJ(2));
t40 = cos(qJ(2));
t24 = t34 * t38 - t35 * t40;
t31 = -t40 * pkin(2) - pkin(1);
t18 = t24 * pkin(3) + t31;
t57 = 0.2e1 * t18;
t56 = 0.2e1 * t31;
t55 = 0.2e1 * t40;
t54 = t34 * pkin(2);
t53 = t35 * pkin(2);
t39 = cos(qJ(5));
t37 = sin(qJ(4));
t46 = -qJ(3) - pkin(6);
t27 = t46 * t38;
t28 = t46 * t40;
t16 = t35 * t27 + t34 * t28;
t25 = t34 * t40 + t35 * t38;
t42 = -t25 * pkin(7) + t16;
t50 = cos(qJ(4));
t17 = t34 * t27 - t35 * t28;
t9 = -t24 * pkin(7) + t17;
t5 = t37 * t9 - t50 * t42;
t52 = t5 * t39;
t30 = pkin(3) + t53;
t21 = t50 * t30 - t37 * t54;
t19 = -pkin(4) - t21;
t51 = pkin(4) - t19;
t14 = t50 * t24 + t37 * t25;
t36 = sin(qJ(5));
t11 = t36 * t14;
t15 = -t37 * t24 + t50 * t25;
t49 = t36 * t15;
t48 = t36 * t39;
t47 = t39 * t15;
t45 = -0.2e1 * t15 * t14;
t44 = -pkin(4) * t15 - pkin(8) * t14;
t22 = -t37 * t30 - t50 * t54;
t20 = pkin(8) - t22;
t43 = -t14 * t20 + t15 * t19;
t33 = t39 ^ 2;
t32 = t36 ^ 2;
t29 = 0.2e1 * t48;
t13 = t15 ^ 2;
t12 = t39 * t14;
t10 = t36 * t47;
t7 = (-t32 + t33) * t15;
t6 = t37 * t42 + t50 * t9;
t4 = t14 * pkin(4) - t15 * pkin(8) + t18;
t3 = t5 * t36;
t2 = t36 * t4 + t39 * t6;
t1 = -t36 * t6 + t39 * t4;
t8 = [1, 0, 0, t38 ^ 2, t38 * t55, 0, 0, 0, pkin(1) * t55, -0.2e1 * pkin(1) * t38, t24 * t56, t25 * t56, -0.2e1 * t16 * t25 - 0.2e1 * t17 * t24, t16 ^ 2 + t17 ^ 2 + t31 ^ 2, t13, t45, 0, 0, 0, t14 * t57, t15 * t57, t33 * t13, -0.2e1 * t13 * t48, 0.2e1 * t14 * t47, t36 * t45, t14 ^ 2, 0.2e1 * t1 * t14 + 0.2e1 * t5 * t49, -0.2e1 * t2 * t14 + 0.2e1 * t5 * t47; 0, 0, 0, 0, 0, t38, t40, 0, -t38 * pkin(6), -t40 * pkin(6), t16, -t17, (-t24 * t34 - t25 * t35) * pkin(2), (t16 * t35 + t17 * t34) * pkin(2), 0, 0, t15, -t14, 0, -t5, -t6, t10, t7, t11, t12, 0, t43 * t36 - t52, t43 * t39 + t3; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t53, -0.2e1 * t54, 0, (t34 ^ 2 + t35 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t21, 0.2e1 * t22, t32, t29, 0, 0, 0, -0.2e1 * t19 * t39, 0.2e1 * t19 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t25, 0, t31, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, t12, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, -t5, -t6, t10, t7, t11, t12, 0, t44 * t36 - t52, t44 * t39 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t21, t22, t32, t29, 0, 0, 0, t51 * t39, -t51 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t32, t29, 0, 0, 0, 0.2e1 * pkin(4) * t39, -0.2e1 * pkin(4) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t49, t14, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t39, 0, -t36 * t20, -t39 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t39, 0, -t36 * pkin(8), -t39 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t8;
