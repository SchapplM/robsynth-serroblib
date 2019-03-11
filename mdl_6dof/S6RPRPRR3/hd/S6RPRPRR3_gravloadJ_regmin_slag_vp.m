% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t22 = qJ(1) + pkin(10);
t17 = sin(t22);
t19 = cos(t22);
t36 = g(1) * t19 + g(2) * t17;
t25 = sin(qJ(3));
t27 = cos(qJ(3));
t11 = -g(3) * t27 + t36 * t25;
t42 = g(3) * t25;
t40 = t17 * t27;
t39 = t19 * t27;
t23 = sin(pkin(11));
t38 = t23 * t27;
t24 = cos(pkin(11));
t37 = t24 * t27;
t21 = pkin(11) + qJ(5);
t35 = g(1) * t17 - g(2) * t19;
t26 = sin(qJ(1));
t28 = cos(qJ(1));
t34 = g(1) * t26 - g(2) * t28;
t33 = t27 * pkin(3) + t25 * qJ(4);
t31 = pkin(2) + t33;
t30 = t34 * pkin(1);
t20 = qJ(6) + t21;
t18 = cos(t21);
t16 = sin(t21);
t15 = cos(t20);
t14 = sin(t20);
t13 = t35 * t25;
t12 = t36 * t27 + t42;
t10 = t17 * t16 + t18 * t39;
t9 = -t16 * t39 + t17 * t18;
t8 = t19 * t16 - t18 * t40;
t7 = t16 * t40 + t19 * t18;
t6 = t17 * t14 + t15 * t39;
t5 = -t14 * t39 + t17 * t15;
t4 = t19 * t14 - t15 * t40;
t3 = t14 * t40 + t19 * t15;
t2 = g(1) * t6 - g(2) * t4 + t15 * t42;
t1 = -g(1) * t5 + g(2) * t3 + t14 * t42;
t29 = [0, t34, g(1) * t28 + g(2) * t26, t30, 0, 0, 0, 0, 0, t35 * t27, -t13, -g(1) * (-t17 * t37 + t19 * t23) - g(2) * (t17 * t23 + t19 * t37) -g(1) * (t17 * t38 + t19 * t24) - g(2) * (t17 * t24 - t19 * t38) t13, t30 + (-g(1) * pkin(7) - g(2) * t31) * t19 + (-g(2) * pkin(7) + g(1) * t31) * t17, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, t11 * t24, -t11 * t23, -t12, -g(3) * t33 + t36 * (pkin(3) * t25 - qJ(4) * t27) 0, 0, 0, 0, 0, t11 * t18, -t11 * t16, 0, 0, 0, 0, 0, t11 * t15, -t11 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t16 * t42, g(1) * t10 - g(2) * t8 + t18 * t42, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t29;
