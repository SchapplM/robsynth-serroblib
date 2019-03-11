% Calculate minimal parameter regressor of gravitation load for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t25 = sin(pkin(6));
t47 = g(3) * t25;
t23 = qJ(3) + qJ(4);
t22 = qJ(5) + t23;
t19 = cos(t22);
t28 = sin(qJ(6));
t46 = t19 * t28;
t31 = cos(qJ(6));
t45 = t19 * t31;
t24 = sin(pkin(12));
t44 = t24 * t25;
t26 = cos(pkin(12));
t43 = t25 * t26;
t29 = sin(qJ(3));
t42 = t25 * t29;
t30 = sin(qJ(2));
t41 = t25 * t30;
t32 = cos(qJ(3));
t40 = t25 * t32;
t27 = cos(pkin(6));
t39 = t27 * t30;
t33 = cos(qJ(2));
t38 = t27 * t33;
t37 = t28 * t33;
t36 = t31 * t33;
t14 = t24 * t33 + t26 * t39;
t16 = -t24 * t39 + t26 * t33;
t18 = sin(t22);
t35 = g(1) * (-t16 * t18 + t19 * t44) + g(2) * (-t14 * t18 - t19 * t43) + g(3) * (-t18 * t41 + t19 * t27);
t13 = t24 * t30 - t26 * t38;
t15 = t24 * t38 + t26 * t30;
t34 = -g(1) * t15 - g(2) * t13 + t33 * t47;
t21 = cos(t23);
t20 = sin(t23);
t12 = t18 * t27 + t19 * t41;
t10 = t16 * t19 + t18 * t44;
t8 = t14 * t19 - t18 * t43;
t6 = -g(1) * (-t16 * t21 - t20 * t44) - g(2) * (-t14 * t21 + t20 * t43) - g(3) * (-t20 * t27 - t21 * t41);
t5 = -g(1) * (-t16 * t20 + t21 * t44) - g(2) * (-t14 * t20 - t21 * t43) - g(3) * (-t20 * t41 + t21 * t27);
t4 = g(1) * t10 + g(2) * t8 + g(3) * t12;
t2 = t35 * t31;
t1 = t35 * t28;
t3 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t34, g(1) * t16 + g(2) * t14 + g(3) * t41, 0, 0, 0, 0, 0, -t34 * t32, t34 * t29, 0, 0, 0, 0, 0, -t34 * t21, t34 * t20, 0, 0, 0, 0, 0, -t34 * t19, t34 * t18, 0, 0, 0, 0, 0, -g(1) * (-t15 * t45 + t16 * t28) - g(2) * (-t13 * t45 + t14 * t28) - (t19 * t36 + t28 * t30) * t47, -g(1) * (t15 * t46 + t16 * t31) - g(2) * (t13 * t46 + t14 * t31) - (-t19 * t37 + t30 * t31) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t16 * t29 + t24 * t40) - g(2) * (-t14 * t29 - t26 * t40) - g(3) * (t27 * t32 - t29 * t41) -g(1) * (-t16 * t32 - t24 * t42) - g(2) * (-t14 * t32 + t26 * t42) - g(3) * (-t27 * t29 - t30 * t40) 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, -t35, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, -t35, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t10 * t28 + t15 * t31) - g(2) * (t13 * t31 - t28 * t8) - g(3) * (-t12 * t28 - t25 * t36) -g(1) * (-t10 * t31 - t15 * t28) - g(2) * (-t13 * t28 - t31 * t8) - g(3) * (-t12 * t31 + t25 * t37);];
taug_reg  = t3;
