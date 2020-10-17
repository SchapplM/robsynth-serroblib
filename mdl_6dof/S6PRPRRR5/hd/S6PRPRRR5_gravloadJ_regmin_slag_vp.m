% Calculate minimal parameter regressor of gravitation load for
% S6PRPRRR5
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
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:18:29
% EndTime: 2019-05-05 01:18:29
% DurationCPUTime: 0.26s
% Computational Cost: add. (231->66), mult. (431->113), div. (0->0), fcn. (534->12), ass. (0->39)
t23 = sin(pkin(6));
t45 = g(3) * t23;
t21 = qJ(4) + qJ(5);
t19 = sin(t21);
t26 = sin(qJ(6));
t44 = t19 * t26;
t29 = cos(qJ(6));
t43 = t19 * t29;
t22 = sin(pkin(11));
t42 = t22 * t23;
t24 = cos(pkin(11));
t41 = t23 * t24;
t27 = sin(qJ(4));
t40 = t23 * t27;
t30 = cos(qJ(4));
t39 = t23 * t30;
t31 = cos(qJ(2));
t38 = t23 * t31;
t25 = cos(pkin(6));
t28 = sin(qJ(2));
t37 = t25 * t28;
t36 = t25 * t31;
t35 = t26 * t28;
t34 = t28 * t29;
t13 = t22 * t28 - t24 * t36;
t15 = t22 * t36 + t24 * t28;
t20 = cos(t21);
t33 = g(1) * (t15 * t20 - t19 * t42) + g(2) * (t13 * t20 + t19 * t41) + g(3) * (-t25 * t19 - t20 * t38);
t5 = -g(1) * t15 - g(2) * t13 + g(3) * t38;
t14 = t22 * t31 + t24 * t37;
t16 = -t22 * t37 + t24 * t31;
t32 = g(1) * t16 + g(2) * t14 + t28 * t45;
t12 = -t19 * t38 + t25 * t20;
t9 = -t13 * t19 + t20 * t41;
t7 = t15 * t19 + t20 * t42;
t4 = g(1) * t7 - g(2) * t9 + g(3) * t12;
t2 = t33 * t29;
t1 = t33 * t26;
t3 = [-g(3), 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t5, t32, t5, -t32, -g(1) * (-t15 * pkin(2) + t16 * qJ(3)) - g(2) * (-t13 * pkin(2) + t14 * qJ(3)) - (pkin(2) * t31 + qJ(3) * t28) * t45, 0, 0, 0, 0, 0, -t32 * t27, -t32 * t30, 0, 0, 0, 0, 0, -t32 * t19, -t32 * t20, 0, 0, 0, 0, 0, -g(1) * (-t15 * t26 + t16 * t43) - g(2) * (-t13 * t26 + t14 * t43) - (t19 * t34 + t26 * t31) * t45, -g(1) * (-t15 * t29 - t16 * t44) - g(2) * (-t13 * t29 - t14 * t44) - (-t19 * t35 + t29 * t31) * t45; 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t15 * t30 - t22 * t40) - g(2) * (t13 * t30 + t24 * t40) - g(3) * (-t25 * t27 - t30 * t38) -g(1) * (-t15 * t27 - t22 * t39) - g(2) * (-t13 * t27 + t24 * t39) - g(3) * (-t25 * t30 + t27 * t38) 0, 0, 0, 0, 0, -t33, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t16 * t29 - t7 * t26) - g(2) * (t14 * t29 + t9 * t26) - g(3) * (-t12 * t26 + t23 * t34) -g(1) * (-t16 * t26 - t7 * t29) - g(2) * (-t14 * t26 + t9 * t29) - g(3) * (-t12 * t29 - t23 * t35);];
taug_reg  = t3;
