% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t32 = sin(pkin(10));
t34 = cos(pkin(10));
t36 = sin(qJ(5));
t39 = cos(qJ(5));
t19 = t36 * t32 - t39 * t34;
t37 = sin(qJ(3));
t14 = t19 * t37;
t56 = 0.2e1 * t14;
t27 = -t34 * pkin(4) - pkin(3);
t55 = 0.2e1 * t27;
t40 = cos(qJ(3));
t54 = 0.2e1 * t40;
t53 = pkin(7) * t32;
t28 = t37 * pkin(7);
t52 = t40 * pkin(7);
t51 = t32 * t37;
t33 = sin(pkin(5));
t50 = t33 * sin(qJ(2));
t49 = t33 * cos(qJ(2));
t48 = t34 * t37;
t47 = pkin(8) + qJ(4);
t22 = -t40 * pkin(3) - t37 * qJ(4) - pkin(2);
t12 = t32 * t22 + t34 * t52;
t46 = t32 ^ 2 + t34 ^ 2;
t35 = cos(pkin(5));
t16 = t35 * t37 + t40 * t50;
t6 = -t16 * t32 - t34 * t49;
t7 = t16 * t34 - t32 * t49;
t45 = -t6 * t32 + t7 * t34;
t44 = -pkin(3) * t37 + qJ(4) * t40;
t18 = t34 * t22;
t11 = -t32 * t52 + t18;
t43 = -t11 * t32 + t12 * t34;
t20 = t39 * t32 + t36 * t34;
t31 = t37 ^ 2;
t24 = t47 * t34;
t23 = t47 * t32;
t21 = pkin(4) * t51 + t28;
t15 = -t35 * t40 + t37 * t50;
t13 = t20 * t37;
t10 = -t36 * t23 + t39 * t24;
t9 = -t39 * t23 - t36 * t24;
t8 = -pkin(8) * t51 + t12;
t5 = -pkin(8) * t48 + t18 + (-pkin(4) - t53) * t40;
t4 = t36 * t6 + t39 * t7;
t3 = -t36 * t7 + t39 * t6;
t2 = t36 * t5 + t39 * t8;
t1 = -t36 * t8 + t39 * t5;
t17 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 ^ 2 + t6 ^ 2 + t7 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t49, -t50, 0, 0, 0, 0, 0, t40 * t49, -t37 * t49, t15 * t51 - t6 * t40, t15 * t48 + t7 * t40, (-t32 * t7 - t34 * t6) * t37, t6 * t11 + t7 * t12 + t15 * t28, 0, 0, 0, 0, 0, t15 * t13 - t3 * t40, -t15 * t14 + t4 * t40; 0, 1, 0, 0, t31, t37 * t54, 0, 0, 0, pkin(2) * t54, -0.2e1 * pkin(2) * t37, -0.2e1 * t11 * t40 + 0.2e1 * t31 * t53, 0.2e1 * t31 * pkin(7) * t34 + 0.2e1 * t12 * t40, 0.2e1 * (-t11 * t34 - t12 * t32) * t37, t31 * pkin(7) ^ 2 + t11 ^ 2 + t12 ^ 2, t14 ^ 2, t13 * t56, t40 * t56, t13 * t54, t40 ^ 2, -0.2e1 * t1 * t40 + 0.2e1 * t21 * t13, -0.2e1 * t21 * t14 + 0.2e1 * t2 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, -t15 * t34, t15 * t32, t45, -t15 * pkin(3) + t45 * qJ(4), 0, 0, 0, 0, 0, t15 * t19, t15 * t20; 0, 0, 0, 0, 0, 0, t37, t40, 0, -t28, -t52, -pkin(7) * t48 + t44 * t32, pkin(7) * t51 + t44 * t34, t43, -pkin(3) * t28 + t43 * qJ(4), -t14 * t20, -t20 * t13 + t14 * t19, -t20 * t40, t19 * t40, 0, t27 * t13 + t21 * t19 - t9 * t40, t10 * t40 - t27 * t14 + t21 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t34, -0.2e1 * pkin(3) * t32, 0.2e1 * t46 * qJ(4), t46 * qJ(4) ^ 2 + pkin(3) ^ 2, t20 ^ 2, -0.2e1 * t20 * t19, 0, 0, 0, t19 * t55, t20 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t48, 0, t28, 0, 0, 0, 0, 0, t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t32, 0, -pkin(3), 0, 0, 0, 0, 0, t19, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t13, -t40, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, t9, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t17;
