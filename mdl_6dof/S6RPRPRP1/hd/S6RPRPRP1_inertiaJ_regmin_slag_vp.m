% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t34 = sin(pkin(10));
t36 = cos(pkin(10));
t39 = sin(qJ(3));
t58 = cos(qJ(3));
t22 = t34 * t39 - t36 * t58;
t20 = t22 ^ 2;
t60 = 0.2e1 * t39;
t40 = cos(qJ(5));
t59 = t40 * pkin(5);
t24 = t34 * t58 + t36 * t39;
t38 = sin(qJ(5));
t32 = t38 ^ 2;
t57 = t32 * t24;
t13 = t38 * t22;
t56 = t38 * t24;
t55 = t38 * t40;
t35 = sin(pkin(9));
t29 = t35 * pkin(1) + pkin(7);
t49 = t58 * t29;
t19 = t58 * qJ(4) + t49;
t48 = (-qJ(4) - t29) * t39;
t10 = t36 * t19 + t34 * t48;
t54 = t40 * t10;
t53 = t40 * t24;
t52 = qJ(6) * t24;
t28 = t34 * pkin(3) + pkin(8);
t51 = qJ(6) + t28;
t50 = pkin(5) * t56;
t37 = cos(pkin(9));
t31 = -t37 * pkin(1) - pkin(2);
t30 = -t36 * pkin(3) - pkin(4);
t26 = -t58 * pkin(3) + t31;
t7 = t22 * pkin(4) - t24 * pkin(8) + t26;
t3 = -t38 * t10 + t40 * t7;
t8 = t34 * t19 - t36 * t48;
t1 = t22 * pkin(5) - t40 * t52 + t3;
t2 = t54 + (t7 - t52) * t38;
t47 = t1 * t40 + t2 * t38;
t46 = -t1 * t38 + t2 * t40;
t17 = t51 * t38;
t18 = t51 * t40;
t45 = -t17 * t40 + t18 * t38;
t44 = t17 * t38 + t18 * t40;
t43 = -t22 * t28 + t24 * t30;
t33 = t40 ^ 2;
t25 = t30 - t59;
t21 = t24 ^ 2;
t16 = t40 * t22;
t15 = t33 * t24;
t14 = t33 * t21;
t11 = -t15 - t57;
t5 = t8 + t50;
t4 = t38 * t7 + t54;
t6 = [1, 0, 0 (t35 ^ 2 + t37 ^ 2) * pkin(1) ^ 2, t39 ^ 2, t58 * t60, 0, 0, 0, -0.2e1 * t31 * t58, t31 * t60, -0.2e1 * t10 * t22 + 0.2e1 * t8 * t24, t10 ^ 2 + t26 ^ 2 + t8 ^ 2, t14, -0.2e1 * t21 * t55, 0.2e1 * t22 * t53, -0.2e1 * t22 * t56, t20, 0.2e1 * t3 * t22 + 0.2e1 * t8 * t56, -0.2e1 * t4 * t22 + 0.2e1 * t8 * t53, -0.2e1 * t47 * t24, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 * t24 + t8 * t22, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t22 + t46 * t24; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t21 + t20, 0, 0, 0, 0, 0, 0, 0, 0, t32 * t21 + t14 + t20; 0, 0, 0, 0, 0, 0, t39, t58, 0, -t39 * t29, -t49 (-t22 * t34 - t24 * t36) * pkin(3) (t10 * t34 - t36 * t8) * pkin(3), t38 * t53, t15 - t57, t13, t16, 0, t43 * t38 - t8 * t40, t8 * t38 + t43 * t40, -t45 * t24 + t46, -t1 * t17 + t2 * t18 + t5 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t39, 0 (-t22 * t36 + t24 * t34) * pkin(3), 0, 0, 0, 0, 0, -t16, t13, -t11, t22 * t25 + t44 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t34 ^ 2 + t36 ^ 2) * pkin(3) ^ 2, t32, 0.2e1 * t55, 0, 0, 0, -0.2e1 * t30 * t40, 0.2e1 * t30 * t38, 0.2e1 * t44, t17 ^ 2 + t18 ^ 2 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, 0, 0, 0, t16, -t13, t11, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t32 + t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t56, t22, t3, -t4, -pkin(5) * t53, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t53, 0, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t40, 0, -t38 * t28, -t40 * t28, -t38 * pkin(5), -t17 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t38, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
