% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t26 = qJ(3) + pkin(10);
t20 = sin(t26);
t22 = cos(t26);
t60 = t22 * pkin(4) + t20 * pkin(8);
t27 = qJ(1) + pkin(9);
t21 = sin(t27);
t23 = cos(t27);
t15 = g(1) * t23 + g(2) * t21;
t33 = cos(qJ(5));
t46 = t23 * t33;
t30 = sin(qJ(5));
t49 = t21 * t30;
t10 = t22 * t49 + t46;
t47 = t23 * t30;
t48 = t21 * t33;
t12 = -t22 * t47 + t48;
t54 = g(3) * t20;
t1 = -g(1) * t12 + g(2) * t10 + t30 * t54;
t7 = -g(3) * t22 + t15 * t20;
t31 = sin(qJ(3));
t59 = pkin(3) * t31;
t32 = sin(qJ(1));
t50 = t32 * pkin(1);
t34 = cos(qJ(3));
t24 = t34 * pkin(3);
t19 = t24 + pkin(2);
t35 = cos(qJ(1));
t25 = t35 * pkin(1);
t45 = t23 * t19 + t25;
t29 = -qJ(4) - pkin(7);
t43 = pkin(5) * t30 - t29;
t14 = g(1) * t21 - g(2) * t23;
t41 = g(1) * t32 - g(2) * t35;
t40 = -t23 * t29 - t50;
t18 = t33 * pkin(5) + pkin(4);
t28 = -qJ(6) - pkin(8);
t39 = t22 * t18 - t20 * t28;
t36 = -g(3) * t34 + t15 * t31;
t13 = t22 * t46 + t49;
t11 = -t22 * t48 + t47;
t9 = t14 * t20;
t8 = t15 * t22 + t54;
t6 = t7 * t33;
t5 = t7 * t30;
t4 = -g(1) * t11 - g(2) * t13;
t3 = -g(1) * t10 - g(2) * t12;
t2 = g(1) * t13 - g(2) * t11 + t33 * t54;
t16 = [0, 0, 0, 0, 0, 0, t41, g(1) * t35 + g(2) * t32, 0, 0, 0, 0, 0, 0, 0, 0, t14, t15, 0, t41 * pkin(1), 0, 0, 0, 0, 0, 0, t14 * t34, -t14 * t31, -t15, -g(1) * (-t21 * pkin(2) + t23 * pkin(7) - t50) - g(2) * (t23 * pkin(2) + t21 * pkin(7) + t25) 0, 0, 0, 0, 0, 0, t14 * t22, -t9, -t15, -g(1) * (-t21 * t19 + t40) - g(2) * (-t21 * t29 + t45) 0, 0, 0, 0, 0, 0, t4, t3, t9, -g(1) * t40 - g(2) * (t60 * t23 + t45) + (-g(1) * (-t19 - t60) + g(2) * t29) * t21, 0, 0, 0, 0, 0, 0, t4, t3, t9, g(1) * t50 - g(2) * t45 + (-g(1) * t43 - g(2) * t39) * t23 + (-g(1) * (-t19 - t39) - g(2) * t43) * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, g(3) * t31 + t15 * t34, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t36 * pkin(3), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * (t24 + t60) + t15 * (pkin(4) * t20 - pkin(8) * t22 + t59) 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * (t24 + t39) + t15 * (t18 * t20 + t22 * t28 + t59); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg  = t16;
