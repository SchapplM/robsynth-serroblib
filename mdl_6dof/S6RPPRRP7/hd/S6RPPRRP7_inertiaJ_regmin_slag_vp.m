% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRP7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t31 = sin(pkin(9));
t32 = cos(pkin(9));
t35 = sin(qJ(4));
t51 = cos(qJ(4));
t15 = t51 * t31 + t35 * t32;
t11 = t15 ^ 2;
t14 = t35 * t31 - t51 * t32;
t12 = t14 ^ 2;
t58 = -t11 - t12;
t57 = -0.2e1 * t14;
t56 = 0.2e1 * t15;
t55 = 2 * qJ(2);
t36 = cos(qJ(5));
t54 = t36 * pkin(5);
t33 = -pkin(1) - qJ(3);
t52 = -pkin(7) + t33;
t17 = t52 * t31;
t18 = t52 * t32;
t9 = t51 * t17 + t35 * t18;
t53 = t36 * t9;
t34 = sin(qJ(5));
t50 = t14 * t34;
t49 = t14 * t36;
t48 = t34 * t15;
t47 = t34 * t36;
t10 = t36 * t15;
t46 = -qJ(6) - pkin(8);
t21 = t31 ^ 2 + t32 ^ 2;
t29 = t34 ^ 2;
t30 = t36 ^ 2;
t45 = t29 + t30;
t44 = qJ(6) * t14;
t22 = t31 * pkin(3) + qJ(2);
t43 = t14 * t56;
t7 = t15 * pkin(4) + t14 * pkin(8) + t22;
t3 = -t34 * t9 + t36 * t7;
t42 = pkin(4) * t14 - pkin(8) * t15;
t1 = t15 * pkin(5) + t36 * t44 + t3;
t2 = t53 + (t7 + t44) * t34;
t41 = t1 * t36 + t2 * t34;
t40 = -t1 * t34 + t2 * t36;
t19 = t46 * t34;
t20 = t46 * t36;
t39 = t36 * t19 - t34 * t20;
t38 = -t19 * t34 - t20 * t36;
t8 = t35 * t17 - t51 * t18;
t37 = qJ(2) ^ 2;
t24 = -pkin(4) - t54;
t13 = t21 * t33;
t5 = -pkin(5) * t50 + t8;
t4 = t34 * t7 + t53;
t6 = [1, 0, 0, -2 * pkin(1), t55, pkin(1) ^ 2 + t37, t31 * t55, t32 * t55, -0.2e1 * t13, t21 * t33 ^ 2 + t37, t12, t43, 0, 0, 0, t22 * t56, t22 * t57, t30 * t12, -0.2e1 * t12 * t47, t10 * t57, t34 * t43, t11, 0.2e1 * t3 * t15 - 0.2e1 * t8 * t50, -0.2e1 * t4 * t15 - 0.2e1 * t8 * t49, 0.2e1 * t41 * t14, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, -t21, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t34, t58 * t36, 0, t5 * t14 + t15 * t40; 0, 0, 0, 0, 0, 1, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t11 + t12; 0, 0, 0, 0, 0, 0, t31, t32, 0, qJ(2), 0, 0, 0, 0, 0, t15, -t14, 0, 0, 0, 0, 0, t10, -t48, t45 * t14, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, 0, -t8, -t9, -t14 * t47 (t29 - t30) * t14, t48, t10, 0, t42 * t34 - t8 * t36, t8 * t34 + t42 * t36, t14 * t39 + t40, t1 * t19 - t2 * t20 + t5 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, 0, 0, 0, 0, 0, -t49, t50, t45 * t15, t14 * t24 + t15 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t29, 0.2e1 * t47, 0, 0, 0, 0.2e1 * pkin(4) * t36, -0.2e1 * pkin(4) * t34, 0.2e1 * t38, t19 ^ 2 + t20 ^ 2 + t24 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, t50, t15, t3, -t4, pkin(5) * t49, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t10, 0, -pkin(5) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t34, 0, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t36, 0, -t34 * pkin(8), -t36 * pkin(8), -t34 * pkin(5), t19 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
