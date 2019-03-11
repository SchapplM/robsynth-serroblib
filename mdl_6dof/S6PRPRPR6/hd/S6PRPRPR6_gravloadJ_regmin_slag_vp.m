% Calculate minimal parameter regressor of gravitation load for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% taug_reg [6x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRPR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t26 = sin(pkin(10));
t30 = sin(qJ(2));
t32 = cos(qJ(2));
t42 = cos(pkin(10));
t43 = cos(pkin(6));
t37 = t43 * t42;
t12 = t26 * t32 + t30 * t37;
t40 = t26 * t43;
t14 = -t30 * t40 + t42 * t32;
t55 = -g(1) * t14 - g(2) * t12;
t27 = sin(pkin(6));
t52 = g(3) * t27;
t24 = pkin(11) + qJ(6);
t22 = sin(t24);
t29 = sin(qJ(4));
t51 = t22 * t29;
t23 = cos(t24);
t50 = t23 * t29;
t25 = sin(pkin(11));
t49 = t25 * t29;
t48 = t26 * t27;
t47 = t27 * t30;
t46 = t27 * t32;
t28 = cos(pkin(11));
t45 = t28 * t29;
t44 = t29 * t30;
t41 = g(3) * (pkin(2) * t46 + qJ(3) * t47);
t39 = t27 * t42;
t31 = cos(qJ(4));
t38 = pkin(4) * t29 - qJ(5) * t31;
t15 = t43 * t29 + t31 * t46;
t13 = t42 * t30 + t32 * t40;
t3 = -t13 * t31 + t29 * t48;
t11 = t26 * t30 - t32 * t37;
t5 = t11 * t31 + t29 * t39;
t35 = g(1) * t3 - g(2) * t5 + g(3) * t15;
t16 = -t29 * t46 + t43 * t31;
t4 = t13 * t29 + t31 * t48;
t6 = -t11 * t29 + t31 * t39;
t34 = g(1) * t4 - g(2) * t6 + g(3) * t16;
t2 = -g(1) * t13 - g(2) * t11 + g(3) * t46;
t33 = g(3) * t47 - t55;
t10 = t13 * pkin(2);
t9 = t11 * pkin(2);
t1 = t33 * t31;
t7 = [-g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t2, t33, t2, -t33, -g(1) * (t14 * qJ(3) - t10) - g(2) * (t12 * qJ(3) - t9) - t41, 0, 0, 0, 0, 0, -t33 * t29, -t1, -g(1) * (-t13 * t25 + t14 * t45) - g(2) * (-t11 * t25 + t12 * t45) - (t25 * t32 + t28 * t44) * t52, -g(1) * (-t13 * t28 - t14 * t49) - g(2) * (-t11 * t28 - t12 * t49) - (-t25 * t44 + t28 * t32) * t52, t1, -g(1) * (-t13 * pkin(8) - t10) - g(2) * (-t11 * pkin(8) - t9) - t41 - (pkin(8) * t32 + t38 * t30) * t52 + t55 * (qJ(3) + t38) 0, 0, 0, 0, 0, -g(1) * (-t13 * t22 + t14 * t50) - g(2) * (-t11 * t22 + t12 * t50) - (t22 * t32 + t23 * t44) * t52, -g(1) * (-t13 * t23 - t14 * t51) - g(2) * (-t11 * t23 - t12 * t51) - (-t22 * t44 + t23 * t32) * t52; 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t34, t35 * t28, -t35 * t25, -t34, -g(1) * (-t3 * pkin(4) + t4 * qJ(5)) - g(2) * (t5 * pkin(4) - t6 * qJ(5)) - g(3) * (-t15 * pkin(4) + t16 * qJ(5)) 0, 0, 0, 0, 0, t35 * t23, -t35 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t14 * t23 - t4 * t22) - g(2) * (t12 * t23 + t6 * t22) - g(3) * (-t16 * t22 + t23 * t47) -g(1) * (-t14 * t22 - t4 * t23) - g(2) * (-t12 * t22 + t6 * t23) - g(3) * (-t16 * t23 - t22 * t47);];
taug_reg  = t7;
