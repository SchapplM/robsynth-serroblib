% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t28 = sin(qJ(2));
t29 = sin(qJ(1));
t31 = cos(qJ(2));
t41 = cos(pkin(5));
t53 = cos(qJ(1));
t36 = t41 * t53;
t12 = t29 * t28 - t31 * t36;
t23 = pkin(10) + qJ(5);
t21 = sin(t23);
t22 = cos(t23);
t13 = t28 * t36 + t29 * t31;
t27 = sin(qJ(3));
t30 = cos(qJ(3));
t25 = sin(pkin(5));
t40 = t25 * t53;
t5 = t13 * t30 - t27 * t40;
t59 = -t12 * t22 + t5 * t21;
t58 = t12 * t21 + t5 * t22;
t39 = t29 * t41;
t14 = t53 * t28 + t31 * t39;
t57 = -g(1) * t14 - g(2) * t12;
t54 = g(3) * t25;
t50 = t21 * t30;
t49 = t22 * t30;
t24 = sin(pkin(10));
t48 = t24 * t30;
t47 = t25 * t28;
t46 = t25 * t29;
t45 = t25 * t30;
t44 = t25 * t31;
t26 = cos(pkin(10));
t43 = t26 * t30;
t42 = t30 * t31;
t4 = t13 * t27 + t30 * t40;
t15 = -t28 * t39 + t53 * t31;
t8 = t15 * t27 - t29 * t45;
t38 = -g(1) * t4 + g(2) * t8;
t37 = -g(1) * t15 - g(2) * t13;
t10 = t27 * t47 - t41 * t30;
t34 = g(1) * t8 + g(2) * t4 + g(3) * t10;
t11 = t41 * t27 + t28 * t45;
t9 = t15 * t30 + t27 * t46;
t33 = g(1) * t9 + g(2) * t5 + g(3) * t11;
t32 = g(3) * t44 + t57;
t3 = t32 * t27;
t2 = t14 * t21 + t9 * t22;
t1 = t14 * t22 - t9 * t21;
t6 = [0, g(1) * t29 - g(2) * t53, g(1) * t53 + g(2) * t29, 0, 0, 0, 0, 0, g(1) * t13 - g(2) * t15, -g(1) * t12 + g(2) * t14, 0, 0, 0, 0, 0, g(1) * t5 - g(2) * t9, t38, -g(1) * (-t12 * t24 - t26 * t5) - g(2) * (t14 * t24 + t9 * t26), -g(1) * (-t12 * t26 + t24 * t5) - g(2) * (t14 * t26 - t9 * t24), -t38, -g(1) * (-t29 * pkin(1) - t13 * pkin(2) - pkin(3) * t5 + pkin(7) * t40 - t12 * pkin(8) - qJ(4) * t4) - g(2) * (t53 * pkin(1) + t15 * pkin(2) + t9 * pkin(3) + pkin(7) * t46 + t14 * pkin(8) + t8 * qJ(4)), 0, 0, 0, 0, 0, g(1) * t58 - g(2) * t2, -g(1) * t59 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -t32, g(3) * t47 - t37, 0, 0, 0, 0, 0, -t32 * t30, t3, -g(1) * (-t14 * t43 + t15 * t24) - g(2) * (-t12 * t43 + t13 * t24) - (t24 * t28 + t26 * t42) * t54, -g(1) * (t14 * t48 + t15 * t26) - g(2) * (t12 * t48 + t13 * t26) - (-t24 * t42 + t26 * t28) * t54, -t3, (-t28 * t54 + t37) * pkin(8) + (-t31 * t54 - t57) * (pkin(3) * t30 + qJ(4) * t27 + pkin(2)), 0, 0, 0, 0, 0, -g(1) * (-t14 * t49 + t15 * t21) - g(2) * (-t12 * t49 + t13 * t21) - (t21 * t28 + t22 * t42) * t54, -g(1) * (t14 * t50 + t15 * t22) - g(2) * (t12 * t50 + t13 * t22) - (-t21 * t42 + t22 * t28) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t33, t34 * t26, -t34 * t24, -t33, -g(1) * (-t8 * pkin(3) + t9 * qJ(4)) - g(2) * (-t4 * pkin(3) + t5 * qJ(4)) - g(3) * (-t10 * pkin(3) + t11 * qJ(4)), 0, 0, 0, 0, 0, t34 * t22, -t34 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t59 - g(3) * (-t11 * t21 - t22 * t44), g(1) * t2 + g(2) * t58 - g(3) * (-t11 * t22 + t21 * t44);];
taug_reg = t6;
