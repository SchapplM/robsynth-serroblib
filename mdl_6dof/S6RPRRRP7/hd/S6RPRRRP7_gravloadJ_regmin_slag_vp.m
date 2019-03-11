% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t37 = sin(qJ(1));
t39 = cos(qJ(1));
t21 = g(1) * t39 + g(2) * t37;
t31 = pkin(10) + qJ(3);
t27 = cos(t31);
t38 = cos(qJ(4));
t48 = t39 * t38;
t36 = sin(qJ(4));
t53 = t37 * t36;
t15 = t27 * t53 + t48;
t49 = t39 * t36;
t52 = t37 * t38;
t17 = -t27 * t49 + t52;
t26 = sin(t31);
t59 = g(3) * t26;
t64 = -g(1) * t17 + g(2) * t15 + t36 * t59;
t42 = -g(3) * t27 + t21 * t26;
t32 = qJ(4) + qJ(5);
t28 = sin(t32);
t57 = t26 * t28;
t29 = cos(t32);
t56 = t26 * t29;
t55 = t37 * t28;
t54 = t37 * t29;
t51 = t39 * t28;
t50 = t39 * t29;
t46 = pkin(4) * t36 + pkin(7) + qJ(2);
t10 = t27 * t55 + t50;
t12 = t27 * t51 - t54;
t45 = g(1) * t10 - g(2) * t12;
t20 = g(1) * t37 - g(2) * t39;
t25 = t38 * pkin(4) + pkin(3);
t44 = pkin(5) * t29 + qJ(6) * t28 + t25;
t34 = cos(pkin(10));
t40 = -pkin(9) - pkin(8);
t43 = t34 * pkin(2) + t27 * t25 - t26 * t40 + pkin(1);
t1 = g(1) * t12 + g(2) * t10 + g(3) * t57;
t11 = t27 * t54 - t51;
t13 = t27 * t50 + t55;
t3 = g(1) * t13 + g(2) * t11 + g(3) * t56;
t41 = -g(1) * (-t12 * pkin(5) + t13 * qJ(6)) - g(2) * (-t10 * pkin(5) + t11 * qJ(6)) - g(3) * (-pkin(5) * t57 + qJ(6) * t56);
t18 = t27 * t48 + t53;
t16 = -t27 * t52 + t49;
t14 = t20 * t26;
t7 = t21 * t27 + t59;
t6 = t42 * t29;
t5 = t42 * t28;
t4 = g(1) * t11 - g(2) * t13;
t2 = [0, t20, t21, t20 * t34, -t20 * sin(pkin(10)) -t21, -g(1) * (-t37 * pkin(1) + t39 * qJ(2)) - g(2) * (t39 * pkin(1) + t37 * qJ(2)) 0, 0, 0, 0, 0, t20 * t27, -t14, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, -g(1) * t15 - g(2) * t17, 0, 0, 0, 0, 0, t4, -t45, t4, t14, t45, -g(1) * (-t11 * pkin(5) - t10 * qJ(6)) - g(2) * (t13 * pkin(5) + t12 * qJ(6)) + (-g(1) * t46 - g(2) * t43) * t39 + (g(1) * t43 - g(2) * t46) * t37; 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t7, 0, 0, 0, 0, 0, t42 * t38, -t42 * t36, 0, 0, 0, 0, 0, t6, -t5, t6, -t7, t5 (-g(3) * t44 + t21 * t40) * t27 + (g(3) * t40 + t21 * t44) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, g(1) * t18 - g(2) * t16 + t38 * t59, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t64 * pkin(4) + t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
