% Calculate minimal parameter regressor of gravitation load for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t23 = sin(pkin(6));
t47 = g(3) * t23;
t19 = pkin(12) + qJ(4);
t16 = cos(t19);
t20 = qJ(5) + qJ(6);
t17 = sin(t20);
t46 = t16 * t17;
t18 = cos(t20);
t45 = t16 * t18;
t25 = sin(qJ(5));
t44 = t16 * t25;
t27 = cos(qJ(5));
t43 = t16 * t27;
t28 = cos(qJ(2));
t42 = t16 * t28;
t22 = sin(pkin(11));
t41 = t22 * t23;
t26 = sin(qJ(2));
t40 = t23 * t26;
t39 = t23 * t28;
t38 = t25 * t28;
t37 = t27 * t28;
t36 = cos(pkin(6));
t35 = cos(pkin(11));
t34 = t22 * t36;
t33 = t23 * t35;
t32 = t36 * t35;
t10 = t22 * t28 + t26 * t32;
t12 = -t26 * t34 + t35 * t28;
t15 = sin(t19);
t31 = g(1) * (-t12 * t15 + t16 * t41) + g(2) * (-t10 * t15 - t16 * t33) + g(3) * (-t15 * t40 + t36 * t16);
t11 = t35 * t26 + t28 * t34;
t9 = t22 * t26 - t28 * t32;
t30 = -g(1) * t11 - g(2) * t9 + g(3) * t39;
t29 = g(1) * t12 + g(2) * t10 + g(3) * t40;
t8 = t36 * t15 + t16 * t40;
t6 = t12 * t16 + t15 * t41;
t4 = t10 * t16 - t15 * t33;
t2 = -g(1) * (-t11 * t17 - t6 * t18) - g(2) * (-t9 * t17 - t4 * t18) - g(3) * (t17 * t39 - t8 * t18);
t1 = -g(1) * (t11 * t18 - t6 * t17) - g(2) * (-t4 * t17 + t9 * t18) - g(3) * (-t8 * t17 - t18 * t39);
t3 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t30, t29, -t30 * cos(pkin(12)) t30 * sin(pkin(12)) -t29, -g(1) * (-t11 * pkin(2) + t12 * qJ(3)) - g(2) * (-t9 * pkin(2) + t10 * qJ(3)) - (pkin(2) * t28 + qJ(3) * t26) * t47, 0, 0, 0, 0, 0, -t30 * t16, t30 * t15, 0, 0, 0, 0, 0, -g(1) * (-t11 * t43 + t12 * t25) - g(2) * (t10 * t25 - t9 * t43) - (t16 * t37 + t25 * t26) * t47, -g(1) * (t11 * t44 + t12 * t27) - g(2) * (t10 * t27 + t9 * t44) - (-t16 * t38 + t26 * t27) * t47, 0, 0, 0, 0, 0, -g(1) * (-t11 * t45 + t12 * t17) - g(2) * (t10 * t17 - t9 * t45) - (t17 * t26 + t18 * t42) * t47, -g(1) * (t11 * t46 + t12 * t18) - g(2) * (t10 * t18 + t9 * t46) - (-t17 * t42 + t18 * t26) * t47; 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, g(1) * t6 + g(2) * t4 + g(3) * t8, 0, 0, 0, 0, 0, -t31 * t27, t31 * t25, 0, 0, 0, 0, 0, -t31 * t18, t31 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t11 * t27 - t6 * t25) - g(2) * (-t4 * t25 + t9 * t27) - g(3) * (-t23 * t37 - t8 * t25) -g(1) * (-t11 * t25 - t6 * t27) - g(2) * (-t9 * t25 - t4 * t27) - g(3) * (t23 * t38 - t8 * t27) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t3;
