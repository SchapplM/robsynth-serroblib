% Calculate minimal parameter regressor of joint inertia matrix for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PPRRRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t33 = sin(qJ(4));
t54 = -0.2e1 * t33;
t36 = cos(qJ(4));
t53 = 0.2e1 * t36;
t35 = cos(qJ(5));
t52 = pkin(4) * t35;
t32 = sin(qJ(5));
t51 = pkin(9) * t32;
t50 = t32 * pkin(5);
t27 = sin(pkin(7));
t34 = sin(qJ(3));
t49 = t27 * t34;
t37 = cos(qJ(3));
t48 = t27 * t37;
t29 = cos(pkin(12));
t30 = cos(pkin(7));
t47 = t29 * t30;
t46 = t32 * t33;
t45 = t32 * t35;
t44 = t32 * t36;
t43 = t35 * t33;
t42 = t35 * t36;
t41 = -qJ(6) - pkin(10);
t40 = qJ(6) * t33;
t39 = t33 * t53;
t38 = pkin(9) * t42;
t31 = cos(pkin(6));
t28 = sin(pkin(6));
t26 = sin(pkin(12));
t25 = t35 ^ 2;
t24 = t33 ^ 2;
t23 = t32 ^ 2;
t21 = -t35 * pkin(5) - pkin(4);
t20 = t41 * t35;
t19 = t41 * t32;
t18 = -t36 * pkin(4) - t33 * pkin(10) - pkin(3);
t17 = (pkin(9) + t50) * t33;
t16 = t35 * t18;
t15 = t33 * t30 + t36 * t49;
t14 = -t36 * t30 + t33 * t49;
t13 = -t28 * t29 * t27 + t31 * t30;
t12 = t32 * t18 + t38;
t11 = -pkin(9) * t44 + t16;
t10 = t38 + (t18 - t40) * t32;
t9 = t35 * t15 - t32 * t48;
t8 = -t32 * t15 - t35 * t48;
t7 = t31 * t49 + (t26 * t37 + t34 * t47) * t28;
t6 = -t31 * t48 + (t26 * t34 - t37 * t47) * t28;
t5 = -t35 * t40 + t16 + (-pkin(5) - t51) * t36;
t4 = t13 * t33 + t7 * t36;
t3 = -t13 * t36 + t7 * t33;
t2 = t6 * t32 + t4 * t35;
t1 = -t4 * t32 + t6 * t35;
t22 = [1, t31 ^ 2 + (t26 ^ 2 + t29 ^ 2) * t28 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t8 + t3 * t14 + t2 * t9; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, -t6, -t7, 0, 0, 0, 0, 0, -t6 * t36, t6 * t33, 0, 0, 0, 0, 0, -t1 * t36 + t3 * t46, t2 * t36 + t3 * t43 (-t1 * t35 - t2 * t32) * t33, t1 * t5 + t2 * t10 + t3 * t17; 0, 0, 0, t48, -t49, 0, 0, 0, 0, 0, t36 * t48, -t33 * t48, 0, 0, 0, 0, 0, t14 * t46 - t8 * t36, t14 * t43 + t9 * t36 (-t32 * t9 - t35 * t8) * t33, t9 * t10 + t14 * t17 + t8 * t5; 0, 0, 1, 0, 0, t24, t39, 0, 0, 0, pkin(3) * t53, pkin(3) * t54, t25 * t24, -0.2e1 * t24 * t45, t42 * t54, t32 * t39, t36 ^ 2, -0.2e1 * t11 * t36 + 0.2e1 * t24 * t51, 0.2e1 * t24 * pkin(9) * t35 + 0.2e1 * t12 * t36, 0.2e1 * (-t10 * t32 - t35 * t5) * t33, t10 ^ 2 + t17 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, 0, 0, 0, 0, -t3 * t35, t3 * t32, -t1 * t32 + t2 * t35, t1 * t19 - t2 * t20 + t3 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, 0, 0, 0, 0, 0, -t14 * t35, t14 * t32, -t8 * t32 + t9 * t35, t14 * t21 + t8 * t19 - t9 * t20; 0, 0, 0, 0, 0, 0, 0, t33, t36, 0, -t33 * pkin(9), -t36 * pkin(9), t32 * t43 (-t23 + t25) * t33, -t44, -t42, 0, -pkin(9) * t43 + (-pkin(4) * t33 + pkin(10) * t36) * t32, pkin(10) * t42 + (t51 - t52) * t33 (-t19 * t33 + t10) * t35 + (t20 * t33 - t5) * t32, -t10 * t20 + t17 * t21 + t5 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t23, 0.2e1 * t45, 0, 0, 0, 0.2e1 * t52, -0.2e1 * pkin(4) * t32, -0.2e1 * t19 * t32 - 0.2e1 * t20 * t35, t19 ^ 2 + t20 ^ 2 + t21 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9, 0, t8 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t46, -t36, t11, -t12, -pkin(5) * t43, t5 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t35, 0, -t32 * pkin(10), -t35 * pkin(10), -t50, t19 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t22;
