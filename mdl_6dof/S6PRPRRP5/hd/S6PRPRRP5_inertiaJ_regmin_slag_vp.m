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
% MM_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
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
t26 = cos(qJ(5));
t52 = 0.2e1 * t26;
t51 = 2 * qJ(3);
t23 = sin(qJ(5));
t50 = t23 * pkin(5);
t21 = sin(pkin(6));
t25 = sin(qJ(2));
t49 = t21 * t25;
t28 = cos(qJ(2));
t48 = t21 * t28;
t24 = sin(qJ(4));
t47 = t23 * t24;
t46 = t23 * t26;
t27 = cos(qJ(4));
t45 = t23 * t27;
t29 = -pkin(2) - pkin(8);
t44 = t23 * t29;
t43 = t24 * t29;
t42 = t26 * t24;
t14 = t26 * t27;
t41 = t26 * t29;
t40 = t27 * t24;
t39 = t27 * t29;
t38 = -qJ(6) - pkin(9);
t17 = t23 ^ 2;
t19 = t26 ^ 2;
t37 = t17 + t19;
t18 = t24 ^ 2;
t20 = t27 ^ 2;
t36 = -t18 - t20;
t35 = qJ(6) * t27;
t34 = -0.2e1 * t40;
t33 = t24 * t41;
t32 = -pkin(4) * t27 - pkin(9) * t24;
t22 = cos(pkin(6));
t8 = t22 * t27 - t24 * t48;
t3 = -t8 * t23 + t26 * t49;
t4 = t23 * t49 + t8 * t26;
t31 = -t3 * t23 + t4 * t26;
t12 = t38 * t23;
t13 = t38 * t26;
t30 = -t12 * t23 - t13 * t26;
t15 = -t26 * pkin(5) - pkin(4);
t11 = t24 * pkin(4) - t27 * pkin(9) + qJ(3);
t10 = (-t29 + t50) * t27;
t9 = t26 * t11;
t7 = t22 * t24 + t27 * t48;
t6 = t23 * t11 + t33;
t5 = -t23 * t43 + t9;
t2 = t33 + (t11 - t35) * t23;
t1 = -t26 * t35 + t9 + (pkin(5) - t44) * t24;
t16 = [1, 0, 0, 0, 0, 0, t22 ^ 2 + (t25 ^ 2 + t28 ^ 2) * t21 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t7 ^ 2; 0, 0, t48, -t49, -t48, t49 (pkin(2) * t28 + qJ(3) * t25) * t21, 0, 0, 0, 0, 0, t24 * t49, t27 * t49, 0, 0, 0, 0, 0, t3 * t24 + t7 * t45, t7 * t14 - t4 * t24 (-t23 * t4 - t26 * t3) * t27, t3 * t1 + t7 * t10 + t4 * t2; 0, 1, 0, 0, -0.2e1 * pkin(2), t51, pkin(2) ^ 2 + (qJ(3) ^ 2) t20, t34, 0, 0, 0, t24 * t51, t27 * t51, t19 * t20, -0.2e1 * t20 * t46, t40 * t52, t23 * t34, t18, -0.2e1 * t20 * t44 + 0.2e1 * t5 * t24, -0.2e1 * t20 * t41 - 0.2e1 * t6 * t24, 0.2e1 * (-t1 * t26 - t2 * t23) * t27, t1 ^ 2 + t10 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, -t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 * t24 - t7 * t27; 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 * t23, t36 * t26, 0, -t10 * t27 + (-t1 * t23 + t2 * t26) * t24; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t18 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, 0, 0, 0, 0, -t7 * t26, t7 * t23, t31, t3 * t12 - t4 * t13 + t7 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t24, 0, t39, -t43, t23 * t14 (-t17 + t19) * t27, t47, t42, 0, t32 * t23 + t26 * t39, -t23 * t39 + t32 * t26 (-t12 * t27 + t2) * t26 + (t13 * t27 - t1) * t23, t1 * t12 + t10 * t15 - t2 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t24, 0, 0, 0, 0, 0, t14, -t45, t37 * t24, -t27 * t15 + t30 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t17, 0.2e1 * t46, 0, 0, 0, pkin(4) * t52, -0.2e1 * pkin(4) * t23, 0.2e1 * t30, t12 ^ 2 + t13 ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t45, t24, t5, -t6, -pkin(5) * t14, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t42, 0, -pkin(5) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t26, 0, -t23 * pkin(9), -t26 * pkin(9), -t50, t12 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t16;
