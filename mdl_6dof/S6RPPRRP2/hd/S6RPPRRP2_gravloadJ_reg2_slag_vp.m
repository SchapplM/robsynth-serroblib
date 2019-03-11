% Calculate inertial parameters regressor of gravitation load for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t30 = qJ(1) + pkin(9);
t25 = sin(t30);
t27 = cos(t30);
t12 = g(1) * t27 + g(2) * t25;
t29 = pkin(10) + qJ(4);
t24 = sin(t29);
t59 = t12 * t24;
t26 = cos(t29);
t47 = t26 * pkin(4) + t24 * pkin(8);
t56 = g(3) * t24;
t35 = sin(qJ(1));
t55 = t35 * pkin(1);
t54 = t24 * t27;
t34 = sin(qJ(5));
t53 = t25 * t34;
t36 = cos(qJ(5));
t52 = t25 * t36;
t51 = t26 * t27;
t50 = t27 * t34;
t49 = t27 * t36;
t32 = cos(pkin(10));
t23 = t32 * pkin(3) + pkin(2);
t37 = cos(qJ(1));
t28 = t37 * pkin(1);
t48 = t27 * t23 + t28;
t46 = pkin(4) * t51 + pkin(8) * t54 + t48;
t7 = t26 * t53 + t49;
t9 = t26 * t50 - t52;
t45 = g(1) * t7 - g(2) * t9;
t11 = g(1) * t25 - g(2) * t27;
t44 = g(1) * t35 - g(2) * t37;
t33 = -pkin(7) - qJ(3);
t43 = -t27 * t33 - t55;
t42 = pkin(5) * t36 + qJ(6) * t34;
t1 = g(1) * t9 + g(2) * t7 + t34 * t56;
t10 = t26 * t49 + t53;
t8 = t26 * t52 - t50;
t40 = g(1) * t10 + g(2) * t8 + t36 * t56;
t39 = -g(3) * t26 + t59;
t38 = (-g(1) * (-t23 - t47) + g(2) * t33) * t25;
t16 = pkin(8) * t51;
t14 = t25 * t26 * pkin(8);
t6 = t11 * t24;
t5 = t12 * t26 + t56;
t4 = t39 * t36;
t3 = t39 * t34;
t2 = g(1) * t8 - g(2) * t10;
t13 = [0, 0, 0, 0, 0, 0, t44, g(1) * t37 + g(2) * t35, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, t44 * pkin(1), 0, 0, 0, 0, 0, 0, t11 * t32, -t11 * sin(pkin(10)) -t12, -g(1) * (-t25 * pkin(2) + t27 * qJ(3) - t55) - g(2) * (t27 * pkin(2) + t25 * qJ(3) + t28) 0, 0, 0, 0, 0, 0, t11 * t26, -t6, -t12, -g(1) * (-t25 * t23 + t43) - g(2) * (-t25 * t33 + t48) 0, 0, 0, 0, 0, 0, t2, -t45, t6, -g(1) * t43 - g(2) * t46 + t38, 0, 0, 0, 0, 0, 0, t2, t6, t45, -g(1) * (-t8 * pkin(5) - t7 * qJ(6) + t43) - g(2) * (t10 * pkin(5) + t9 * qJ(6) + t46) + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t5, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t5, -g(1) * (-pkin(4) * t54 + t16) - g(2) * (-t25 * t24 * pkin(4) + t14) - g(3) * t47, 0, 0, 0, 0, 0, 0, t4, -t5, t3, -g(1) * t16 - g(2) * t14 - g(3) * (t42 * t26 + t47) + (pkin(4) + t42) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t40, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t40, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - (-pkin(5) * t34 + qJ(6) * t36) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t13;
