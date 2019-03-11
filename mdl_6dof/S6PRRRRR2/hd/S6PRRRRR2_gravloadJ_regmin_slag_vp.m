% Calculate minimal parameter regressor of gravitation load for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t27 = sin(pkin(6));
t52 = g(3) * t27;
t24 = qJ(5) + qJ(6);
t20 = sin(t24);
t25 = qJ(3) + qJ(4);
t23 = cos(t25);
t51 = t20 * t23;
t22 = cos(t24);
t50 = t22 * t23;
t28 = sin(qJ(5));
t49 = t23 * t28;
t31 = cos(qJ(5));
t48 = t23 * t31;
t33 = cos(qJ(2));
t47 = t23 * t33;
t26 = sin(pkin(12));
t46 = t26 * t27;
t30 = sin(qJ(2));
t45 = t27 * t30;
t32 = cos(qJ(3));
t44 = t27 * t32;
t43 = t27 * t33;
t42 = t28 * t33;
t41 = t31 * t33;
t40 = cos(pkin(6));
t39 = cos(pkin(12));
t38 = t26 * t40;
t37 = t27 * t39;
t36 = t40 * t39;
t16 = t26 * t33 + t30 * t36;
t18 = -t30 * t38 + t39 * t33;
t21 = sin(t25);
t35 = g(1) * (-t18 * t21 + t23 * t46) + g(2) * (-t16 * t21 - t23 * t37) + g(3) * (-t21 * t45 + t40 * t23);
t15 = t26 * t30 - t33 * t36;
t17 = t39 * t30 + t33 * t38;
t34 = -g(1) * t17 - g(2) * t15 + g(3) * t43;
t29 = sin(qJ(3));
t14 = t40 * t21 + t23 * t45;
t12 = t18 * t23 + t21 * t46;
t10 = t16 * t23 - t21 * t37;
t8 = g(1) * t12 + g(2) * t10 + g(3) * t14;
t6 = t35 * t31;
t5 = t35 * t28;
t4 = t35 * t22;
t3 = t35 * t20;
t2 = -g(1) * (-t12 * t22 - t17 * t20) - g(2) * (-t10 * t22 - t15 * t20) - g(3) * (-t14 * t22 + t20 * t43);
t1 = -g(1) * (-t12 * t20 + t17 * t22) - g(2) * (-t10 * t20 + t15 * t22) - g(3) * (-t14 * t20 - t22 * t43);
t7 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t34, g(1) * t18 + g(2) * t16 + g(3) * t45, 0, 0, 0, 0, 0, -t34 * t32, t34 * t29, 0, 0, 0, 0, 0, -t34 * t23, t34 * t21, 0, 0, 0, 0, 0, -g(1) * (-t17 * t48 + t18 * t28) - g(2) * (-t15 * t48 + t16 * t28) - (t23 * t41 + t28 * t30) * t52, -g(1) * (t17 * t49 + t18 * t31) - g(2) * (t15 * t49 + t16 * t31) - (-t23 * t42 + t30 * t31) * t52, 0, 0, 0, 0, 0, -g(1) * (-t17 * t50 + t18 * t20) - g(2) * (-t15 * t50 + t16 * t20) - (t20 * t30 + t22 * t47) * t52, -g(1) * (t17 * t51 + t18 * t22) - g(2) * (t15 * t51 + t16 * t22) - (-t20 * t47 + t22 * t30) * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t18 * t29 + t26 * t44) - g(2) * (-t16 * t29 - t32 * t37) - g(3) * (-t29 * t45 + t40 * t32) -g(1) * (-t18 * t32 - t29 * t46) - g(2) * (-t16 * t32 + t29 * t37) - g(3) * (-t40 * t29 - t30 * t44) 0, 0, 0, 0, 0, -t35, t8, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t8, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t12 * t28 + t17 * t31) - g(2) * (-t10 * t28 + t15 * t31) - g(3) * (-t14 * t28 - t27 * t41) -g(1) * (-t12 * t31 - t17 * t28) - g(2) * (-t10 * t31 - t15 * t28) - g(3) * (-t14 * t31 + t27 * t42) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t7;
