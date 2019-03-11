% Calculate minimal parameter regressor of gravitation load for
% S6PRPRRR2
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
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t26 = sin(pkin(12));
t33 = sin(qJ(2));
t36 = cos(qJ(2));
t45 = cos(pkin(12));
t41 = -t33 * t26 + t36 * t45;
t28 = sin(pkin(6));
t55 = g(3) * t28;
t25 = qJ(5) + qJ(6);
t23 = sin(t25);
t35 = cos(qJ(4));
t54 = t23 * t35;
t24 = cos(t25);
t53 = t24 * t35;
t32 = sin(qJ(4));
t52 = t28 * t32;
t51 = t28 * t35;
t30 = cos(pkin(6));
t50 = t30 * t33;
t49 = t30 * t36;
t31 = sin(qJ(5));
t48 = t31 * t35;
t34 = cos(qJ(5));
t46 = t34 * t35;
t20 = -t36 * t26 - t33 * t45;
t18 = t20 * t30;
t27 = sin(pkin(11));
t29 = cos(pkin(11));
t43 = -t18 * t29 + t27 * t41;
t42 = t18 * t27 + t29 * t41;
t17 = t20 * t28;
t40 = g(1) * (t27 * t51 - t32 * t42) + g(2) * (-t29 * t51 - t32 * t43) + g(3) * (t17 * t32 + t30 * t35);
t38 = t30 * t41;
t11 = t20 * t29 - t27 * t38;
t16 = t41 * t28;
t8 = t27 * t20 + t29 * t38;
t39 = g(1) * t11 + g(2) * t8 + g(3) * t16;
t37 = -g(1) * (-t27 * t49 - t29 * t33) - g(2) * (-t27 * t33 + t29 * t49) - t36 * t55;
t14 = -t17 * t35 + t30 * t32;
t6 = t27 * t52 + t35 * t42;
t4 = -t29 * t52 + t35 * t43;
t2 = -g(1) * (t11 * t23 - t24 * t6) - g(2) * (t23 * t8 - t24 * t4) - g(3) * (-t14 * t24 + t16 * t23);
t1 = -g(1) * (-t11 * t24 - t23 * t6) - g(2) * (-t23 * t4 - t24 * t8) - g(3) * (-t14 * t23 - t16 * t24);
t3 = [-g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t37, -g(1) * (t27 * t50 - t29 * t36) - g(2) * (-t27 * t36 - t29 * t50) + t33 * t55, t37 * pkin(2), 0, 0, 0, 0, 0, -t39 * t35, t39 * t32, 0, 0, 0, 0, 0, -g(1) * (t11 * t46 + t31 * t42) - g(2) * (t31 * t43 + t8 * t46) - g(3) * (t16 * t46 - t17 * t31) -g(1) * (-t11 * t48 + t34 * t42) - g(2) * (t34 * t43 - t8 * t48) - g(3) * (-t16 * t48 - t17 * t34) 0, 0, 0, 0, 0, -g(1) * (t11 * t53 + t23 * t42) - g(2) * (t23 * t43 + t8 * t53) - g(3) * (t16 * t53 - t17 * t23) -g(1) * (-t11 * t54 + t24 * t42) - g(2) * (t24 * t43 - t8 * t54) - g(3) * (-t16 * t54 - t17 * t24); 0, 0, 0, 0, -g(3) * t30 + (-g(1) * t27 + g(2) * t29) * t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, g(1) * t6 + g(2) * t4 + g(3) * t14, 0, 0, 0, 0, 0, -t40 * t34, t40 * t31, 0, 0, 0, 0, 0, -t40 * t24, t40 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t11 * t34 - t31 * t6) - g(2) * (-t31 * t4 - t34 * t8) - g(3) * (-t14 * t31 - t16 * t34) -g(1) * (t11 * t31 - t34 * t6) - g(2) * (t31 * t8 - t34 * t4) - g(3) * (-t14 * t34 + t16 * t31) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t3;
