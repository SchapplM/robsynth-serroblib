% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRR2
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
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t29 = sin(qJ(1));
t32 = cos(qJ(1));
t18 = g(1) * t32 + g(2) * t29;
t22 = qJ(2) + pkin(11) + qJ(4);
t19 = sin(t22);
t20 = cos(t22);
t7 = -g(3) * t20 + t18 * t19;
t44 = g(3) * t19;
t25 = qJ(5) + qJ(6);
t23 = sin(t25);
t42 = t29 * t23;
t24 = cos(t25);
t41 = t29 * t24;
t27 = sin(qJ(5));
t40 = t29 * t27;
t30 = cos(qJ(5));
t39 = t29 * t30;
t38 = t32 * t23;
t37 = t32 * t24;
t36 = t32 * t27;
t35 = t32 * t30;
t17 = g(1) * t29 - g(2) * t32;
t28 = sin(qJ(2));
t31 = cos(qJ(2));
t33 = -g(3) * t31 + t18 * t28;
t26 = -qJ(3) - pkin(7);
t21 = t31 * pkin(2) + pkin(1);
t16 = t20 * t35 + t40;
t15 = -t20 * t36 + t39;
t14 = -t20 * t39 + t36;
t13 = t20 * t40 + t35;
t12 = t20 * t37 + t42;
t11 = -t20 * t38 + t41;
t10 = -t20 * t41 + t38;
t9 = t20 * t42 + t37;
t8 = t18 * t20 + t44;
t6 = t7 * t30;
t5 = t7 * t27;
t4 = t7 * t24;
t3 = t7 * t23;
t2 = g(1) * t12 - g(2) * t10 + t24 * t44;
t1 = -g(1) * t11 + g(2) * t9 + t23 * t44;
t34 = [0, t17, t18, 0, 0, 0, 0, 0, t17 * t31, -t17 * t28, -t18, -g(1) * (-t29 * t21 - t32 * t26) - g(2) * (t32 * t21 - t29 * t26) 0, 0, 0, 0, 0, t17 * t20, -t17 * t19, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11; 0, 0, 0, 0, 0, 0, 0, 0, t33, g(3) * t28 + t18 * t31, 0, t33 * pkin(2), 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 + g(2) * t13 + t27 * t44, g(1) * t16 - g(2) * t14 + t30 * t44, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t34;
