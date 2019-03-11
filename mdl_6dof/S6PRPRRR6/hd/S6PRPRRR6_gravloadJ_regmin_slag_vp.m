% Calculate minimal parameter regressor of gravitation load for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t21 = sin(pkin(6));
t44 = g(3) * t21;
t19 = qJ(5) + qJ(6);
t17 = sin(t19);
t23 = sin(qJ(4));
t43 = t17 * t23;
t18 = cos(t19);
t42 = t18 * t23;
t20 = sin(pkin(11));
t41 = t20 * t21;
t24 = sin(qJ(2));
t40 = t21 * t24;
t27 = cos(qJ(2));
t39 = t21 * t27;
t22 = sin(qJ(5));
t38 = t22 * t23;
t37 = t23 * t24;
t25 = cos(qJ(5));
t36 = t23 * t25;
t35 = t24 * t25;
t34 = cos(pkin(6));
t33 = cos(pkin(11));
t32 = t20 * t34;
t31 = t21 * t33;
t30 = t34 * t33;
t11 = t33 * t24 + t27 * t32;
t26 = cos(qJ(4));
t9 = t20 * t24 - t27 * t30;
t29 = g(1) * (t11 * t26 - t23 * t41) + g(2) * (t23 * t31 + t9 * t26) + g(3) * (-t34 * t23 - t26 * t39);
t3 = -g(1) * t11 - g(2) * t9 + g(3) * t39;
t10 = t20 * t27 + t24 * t30;
t12 = -t24 * t32 + t33 * t27;
t28 = g(1) * t12 + g(2) * t10 + g(3) * t40;
t14 = -t23 * t39 + t34 * t26;
t7 = -t9 * t23 + t26 * t31;
t5 = t11 * t23 + t26 * t41;
t2 = -g(1) * (-t12 * t17 - t5 * t18) - g(2) * (-t10 * t17 + t7 * t18) - g(3) * (-t14 * t18 - t17 * t40);
t1 = -g(1) * (t12 * t18 - t5 * t17) - g(2) * (t10 * t18 + t7 * t17) - g(3) * (-t14 * t17 + t18 * t40);
t4 = [-g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t3, t28, t3, -t28, -g(1) * (-t11 * pkin(2) + t12 * qJ(3)) - g(2) * (-t9 * pkin(2) + t10 * qJ(3)) - (pkin(2) * t27 + qJ(3) * t24) * t44, 0, 0, 0, 0, 0, -t28 * t23, -t28 * t26, 0, 0, 0, 0, 0, -g(1) * (-t11 * t22 + t12 * t36) - g(2) * (t10 * t36 - t9 * t22) - (t22 * t27 + t23 * t35) * t44, -g(1) * (-t11 * t25 - t12 * t38) - g(2) * (-t10 * t38 - t9 * t25) - (-t22 * t37 + t25 * t27) * t44, 0, 0, 0, 0, 0, -g(1) * (-t11 * t17 + t12 * t42) - g(2) * (t10 * t42 - t9 * t17) - (t17 * t27 + t18 * t37) * t44, -g(1) * (-t11 * t18 - t12 * t43) - g(2) * (-t10 * t43 - t9 * t18) - (-t17 * t37 + t18 * t27) * t44; 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, g(1) * t5 - g(2) * t7 + g(3) * t14, 0, 0, 0, 0, 0, -t29 * t25, t29 * t22, 0, 0, 0, 0, 0, -t29 * t18, t29 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t25 - t5 * t22) - g(2) * (t10 * t25 + t7 * t22) - g(3) * (-t14 * t22 + t21 * t35) -g(1) * (-t12 * t22 - t5 * t25) - g(2) * (-t10 * t22 + t7 * t25) - g(3) * (-t14 * t25 - t22 * t40) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t4;
