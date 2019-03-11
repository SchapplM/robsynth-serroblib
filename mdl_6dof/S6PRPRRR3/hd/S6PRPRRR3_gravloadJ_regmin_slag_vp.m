% Calculate minimal parameter regressor of gravitation load for
% S6PRPRRR3
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
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t25 = sin(pkin(6));
t45 = g(3) * t25;
t22 = pkin(12) + qJ(4);
t21 = qJ(5) + t22;
t18 = cos(t21);
t27 = sin(qJ(6));
t44 = t18 * t27;
t29 = cos(qJ(6));
t43 = t18 * t29;
t24 = sin(pkin(11));
t42 = t24 * t25;
t28 = sin(qJ(2));
t41 = t25 * t28;
t30 = cos(qJ(2));
t40 = t27 * t30;
t39 = t29 * t30;
t38 = cos(pkin(6));
t37 = cos(pkin(11));
t36 = t24 * t38;
t35 = t25 * t37;
t34 = t38 * t37;
t12 = t24 * t30 + t28 * t34;
t14 = -t28 * t36 + t37 * t30;
t17 = sin(t21);
t33 = g(1) * (-t14 * t17 + t18 * t42) + g(2) * (-t12 * t17 - t18 * t35) + g(3) * (-t17 * t41 + t38 * t18);
t11 = t24 * t28 - t30 * t34;
t13 = t37 * t28 + t30 * t36;
t32 = -g(1) * t13 - g(2) * t11 + t30 * t45;
t31 = g(1) * t14 + g(2) * t12 + g(3) * t41;
t20 = cos(t22);
t19 = sin(t22);
t10 = t38 * t17 + t18 * t41;
t8 = t14 * t18 + t17 * t42;
t6 = t12 * t18 - t17 * t35;
t4 = g(1) * t8 + g(2) * t6 + g(3) * t10;
t2 = t33 * t29;
t1 = t33 * t27;
t3 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t32, t31, -t32 * cos(pkin(12)) t32 * sin(pkin(12)) -t31, -g(1) * (-t13 * pkin(2) + t14 * qJ(3)) - g(2) * (-t11 * pkin(2) + t12 * qJ(3)) - (pkin(2) * t30 + qJ(3) * t28) * t45, 0, 0, 0, 0, 0, -t32 * t20, t32 * t19, 0, 0, 0, 0, 0, -t32 * t18, t32 * t17, 0, 0, 0, 0, 0, -g(1) * (-t13 * t43 + t14 * t27) - g(2) * (-t11 * t43 + t12 * t27) - (t18 * t39 + t27 * t28) * t45, -g(1) * (t13 * t44 + t14 * t29) - g(2) * (t11 * t44 + t12 * t29) - (-t18 * t40 + t28 * t29) * t45; 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t19 + t20 * t42) - g(2) * (-t12 * t19 - t20 * t35) - g(3) * (-t19 * t41 + t38 * t20) -g(1) * (-t14 * t20 - t19 * t42) - g(2) * (-t12 * t20 + t19 * t35) - g(3) * (-t38 * t19 - t20 * t41) 0, 0, 0, 0, 0, -t33, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t13 * t29 - t8 * t27) - g(2) * (t11 * t29 - t6 * t27) - g(3) * (-t10 * t27 - t25 * t39) -g(1) * (-t13 * t27 - t8 * t29) - g(2) * (-t11 * t27 - t6 * t29) - g(3) * (-t10 * t29 + t25 * t40);];
taug_reg  = t3;
