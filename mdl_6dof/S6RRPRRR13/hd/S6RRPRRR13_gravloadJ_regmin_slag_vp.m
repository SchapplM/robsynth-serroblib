% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR13_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR13_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t28 = sin(qJ(2));
t29 = sin(qJ(1));
t32 = cos(qJ(2));
t33 = cos(qJ(1));
t42 = cos(pkin(6));
t40 = t33 * t42;
t17 = t28 * t40 + t29 * t32;
t24 = qJ(5) + qJ(6);
t22 = sin(t24);
t23 = cos(t24);
t16 = t29 * t28 - t32 * t40;
t27 = sin(qJ(4));
t31 = cos(qJ(4));
t25 = sin(pkin(6));
t47 = t25 * t33;
t36 = -t16 * t27 + t31 * t47;
t61 = t17 * t23 + t22 * t36;
t60 = -t17 * t22 + t23 * t36;
t26 = sin(qJ(5));
t30 = cos(qJ(5));
t59 = t17 * t30 + t26 * t36;
t58 = -t17 * t26 + t30 * t36;
t57 = g(3) * t25;
t52 = t22 * t27;
t51 = t23 * t27;
t50 = t25 * t28;
t49 = t25 * t29;
t48 = t25 * t32;
t46 = t26 * t27;
t45 = t27 * t28;
t44 = t27 * t30;
t43 = t28 * t30;
t41 = t29 * t42;
t18 = t33 * t28 + t32 * t41;
t39 = g(1) * t16 - g(2) * t18;
t19 = -t28 * t41 + t33 * t32;
t38 = g(1) * t17 - g(2) * t19;
t37 = g(1) * t33 + g(2) * t29;
t10 = t16 * t31 + t27 * t47;
t8 = t18 * t31 - t27 * t49;
t35 = g(1) * t8 + g(2) * t10 + g(3) * (-t42 * t27 - t31 * t48);
t7 = -g(1) * t18 - g(2) * t16 + g(3) * t48;
t34 = g(1) * t19 + g(2) * t17 + g(3) * t50;
t15 = -t27 * t48 + t42 * t31;
t9 = t18 * t27 + t31 * t49;
t6 = t19 * t26 + t9 * t30;
t5 = t19 * t30 - t9 * t26;
t4 = t19 * t22 + t9 * t23;
t3 = t19 * t23 - t9 * t22;
t2 = g(1) * t4 - g(2) * t60 - g(3) * (-t15 * t23 - t22 * t50);
t1 = -g(1) * t3 - g(2) * t61 - g(3) * (-t15 * t22 + t23 * t50);
t11 = [0, g(1) * t29 - g(2) * t33, t37, 0, 0, 0, 0, 0, t38, -t39, -t37 * t25, -t38, t39, -g(1) * (-t29 * pkin(1) - t17 * pkin(2) + pkin(8) * t47 - t16 * qJ(3)) - g(2) * (t33 * pkin(1) + t19 * pkin(2) + pkin(8) * t49 + t18 * qJ(3)) 0, 0, 0, 0, 0, -g(1) * t36 - g(2) * t9, g(1) * t10 - g(2) * t8, 0, 0, 0, 0, 0, -g(1) * t58 - g(2) * t6, g(1) * t59 - g(2) * t5, 0, 0, 0, 0, 0, -g(1) * t60 - g(2) * t4, g(1) * t61 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, -t7, t34, 0, t7, -t34, -g(1) * (-t18 * pkin(2) + t19 * qJ(3)) - g(2) * (-t16 * pkin(2) + t17 * qJ(3)) - (pkin(2) * t32 + qJ(3) * t28) * t57, 0, 0, 0, 0, 0, -t34 * t27, -t34 * t31, 0, 0, 0, 0, 0, -g(1) * (-t18 * t26 + t19 * t44) - g(2) * (-t16 * t26 + t17 * t44) - (t26 * t32 + t27 * t43) * t57, -g(1) * (-t18 * t30 - t19 * t46) - g(2) * (-t16 * t30 - t17 * t46) - (-t26 * t45 + t30 * t32) * t57, 0, 0, 0, 0, 0, -g(1) * (-t18 * t22 + t19 * t51) - g(2) * (-t16 * t22 + t17 * t51) - (t22 * t32 + t23 * t45) * t57, -g(1) * (-t18 * t23 - t19 * t52) - g(2) * (-t16 * t23 - t17 * t52) - (-t22 * t45 + t23 * t32) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, g(1) * t9 - g(2) * t36 + g(3) * t15, 0, 0, 0, 0, 0, -t35 * t30, t35 * t26, 0, 0, 0, 0, 0, -t35 * t23, t35 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t59 - g(3) * (-t15 * t26 + t25 * t43) g(1) * t6 - g(2) * t58 - g(3) * (-t15 * t30 - t26 * t50) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t11;
