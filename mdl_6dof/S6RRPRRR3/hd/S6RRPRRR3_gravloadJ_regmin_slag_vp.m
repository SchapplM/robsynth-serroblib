% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t32 = sin(qJ(1));
t35 = cos(qJ(1));
t18 = g(1) * t35 + g(2) * t32;
t27 = qJ(2) + pkin(11);
t22 = sin(t27);
t23 = cos(t27);
t38 = -g(3) * t23 + t18 * t22;
t52 = g(3) * t22;
t28 = qJ(4) + qJ(5);
t26 = qJ(6) + t28;
t19 = sin(t26);
t50 = t32 * t19;
t20 = cos(t26);
t49 = t32 * t20;
t24 = sin(t28);
t48 = t32 * t24;
t25 = cos(t28);
t47 = t32 * t25;
t30 = sin(qJ(4));
t46 = t32 * t30;
t33 = cos(qJ(4));
t45 = t32 * t33;
t44 = t35 * t19;
t43 = t35 * t20;
t42 = t35 * t24;
t41 = t35 * t25;
t40 = t35 * t30;
t39 = t35 * t33;
t17 = g(1) * t32 - g(2) * t35;
t31 = sin(qJ(2));
t34 = cos(qJ(2));
t36 = -g(3) * t34 + t18 * t31;
t29 = -qJ(3) - pkin(7);
t21 = t34 * pkin(2) + pkin(1);
t16 = t23 * t39 + t46;
t15 = -t23 * t40 + t45;
t14 = -t23 * t45 + t40;
t13 = t23 * t46 + t39;
t12 = t23 * t41 + t48;
t11 = -t23 * t42 + t47;
t10 = -t23 * t47 + t42;
t9 = t23 * t48 + t41;
t8 = t23 * t43 + t50;
t7 = -t23 * t44 + t49;
t6 = -t23 * t49 + t44;
t5 = t23 * t50 + t43;
t4 = g(1) * t12 - g(2) * t10 + t25 * t52;
t3 = -g(1) * t11 + g(2) * t9 + t24 * t52;
t2 = g(1) * t8 - g(2) * t6 + t20 * t52;
t1 = -g(1) * t7 + g(2) * t5 + t19 * t52;
t37 = [0, t17, t18, 0, 0, 0, 0, 0, t17 * t34, -t17 * t31, -t18, -g(1) * (-t32 * t21 - t35 * t29) - g(2) * (t35 * t21 - t32 * t29) 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, t36, g(3) * t31 + t18 * t34, 0, t36 * pkin(2), 0, 0, 0, 0, 0, t38 * t33, -t38 * t30, 0, 0, 0, 0, 0, t38 * t25, -t38 * t24, 0, 0, 0, 0, 0, t38 * t20, -t38 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 + g(2) * t13 + t30 * t52, g(1) * t16 - g(2) * t14 + t33 * t52, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t37;
