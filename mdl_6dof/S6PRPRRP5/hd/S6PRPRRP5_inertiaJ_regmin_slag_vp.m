% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:51:50
% EndTime: 2021-01-16 01:51:52
% DurationCPUTime: 0.50s
% Computational Cost: add. (293->92), mult. (620->162), div. (0->0), fcn. (691->8), ass. (0->61)
t61 = 2 * pkin(5);
t31 = cos(qJ(5));
t60 = 0.2e1 * t31;
t59 = 2 * qJ(3);
t28 = sin(qJ(5));
t58 = t28 * pkin(5);
t27 = cos(pkin(6));
t29 = sin(qJ(4));
t32 = cos(qJ(4));
t26 = sin(pkin(6));
t33 = cos(qJ(2));
t55 = t26 * t33;
t10 = t27 * t29 + t32 * t55;
t57 = t10 * t31;
t30 = sin(qJ(2));
t56 = t26 * t30;
t54 = t28 * t29;
t53 = t28 * t31;
t52 = t28 * t32;
t34 = -pkin(2) - pkin(8);
t51 = t28 * t34;
t50 = t29 * t34;
t49 = t31 * t29;
t19 = t31 * t32;
t48 = t31 * t34;
t20 = -t31 * pkin(5) - pkin(4);
t47 = t32 * t20;
t46 = t32 * t29;
t45 = t32 * t34;
t44 = -qJ(6) - pkin(9);
t22 = t28 ^ 2;
t24 = t31 ^ 2;
t43 = t22 + t24;
t23 = t29 ^ 2;
t25 = t32 ^ 2;
t42 = -t23 - t25;
t41 = qJ(6) * t32;
t40 = -0.2e1 * t46;
t39 = t29 * t48;
t16 = t29 * pkin(4) - t32 * pkin(9) + qJ(3);
t12 = t31 * t16;
t38 = -t31 * t41 + t12;
t37 = -pkin(4) * t32 - pkin(9) * t29;
t11 = t27 * t32 - t29 * t55;
t5 = -t11 * t28 + t31 * t56;
t6 = t11 * t31 + t28 * t56;
t36 = -t5 * t28 + t6 * t31;
t17 = t44 * t28;
t18 = t44 * t31;
t35 = -t17 * t28 - t18 * t31;
t15 = t42 * t31;
t14 = t42 * t28;
t13 = (-t34 + t58) * t32;
t9 = t10 * t28;
t8 = t28 * t16 + t39;
t7 = -t28 * t50 + t12;
t4 = t39 + (t16 - t41) * t28;
t3 = (pkin(5) - t51) * t29 + t38;
t2 = t10 * t19 - t6 * t29;
t1 = t10 * t52 + t5 * t29;
t21 = [1, 0, 0, 0, 0, 0, t27 ^ 2 + (t30 ^ 2 + t33 ^ 2) * t26 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, t55, -t56, -t55, t56, (pkin(2) * t33 + qJ(3) * t30) * t26, 0, 0, 0, 0, 0, t29 * t56, t32 * t56, 0, 0, 0, 0, 0, t1, t2, t1, t2, (-t28 * t6 - t31 * t5) * t32, t10 * t13 + t5 * t3 + t6 * t4; 0, 1, 0, 0, -0.2e1 * pkin(2), t59, pkin(2) ^ 2 + (qJ(3) ^ 2), t25, t40, 0, 0, 0, t29 * t59, t32 * t59, t24 * t25, -0.2e1 * t25 * t53, t46 * t60, t28 * t40, t23, -0.2e1 * t25 * t51 + 0.2e1 * t7 * t29, -0.2e1 * t25 * t48 - 0.2e1 * t8 * t29, 0.2e1 * t13 * t52 + 0.2e1 * t3 * t29, 0.2e1 * t13 * t19 - 0.2e1 * t4 * t29, 0.2e1 * (-t28 * t4 - t3 * t31) * t32, t13 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, -t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t32 + t36 * t29; 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t15, t14, t15, 0, -t13 * t32 + (-t3 * t28 + t4 * t31) * t29; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 * t23 + t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t57, t9, -t57, t9, t36, t10 * t20 + t5 * t17 - t6 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t29, 0, t45, -t50, t28 * t19, (-t22 + t24) * t32, t54, t49, 0, t37 * t28 + t31 * t45, -t28 * t45 + t37 * t31, -t13 * t31 + t17 * t29 + t28 * t47, t13 * t28 + t18 * t29 + t31 * t47, (-t17 * t32 + t4) * t31 + (t18 * t32 - t3) * t28, t13 * t20 + t3 * t17 - t4 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t29, 0, 0, 0, 0, 0, t19, -t52, t19, -t52, t43 * t29, t35 * t29 - t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t22, 0.2e1 * t53, 0, 0, 0, pkin(4) * t60, -0.2e1 * pkin(4) * t28, -0.2e1 * t20 * t31, 0.2e1 * t20 * t28, 0.2e1 * t35, t17 ^ 2 + t18 ^ 2 + t20 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, t5, -t6, 0, t5 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t52, t29, t7, -t8, (t61 - t51) * t29 + t38, -t4, -pkin(5) * t19, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t49, -t54, -t49, 0, -pkin(5) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t31, 0, -t28 * pkin(9), -t31 * pkin(9), t17, t18, -t58, t17 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t61, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t19, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t28, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t21;
