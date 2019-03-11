% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t30 = cos(qJ(2));
t24 = qJ(2) + qJ(3);
t20 = cos(t24);
t21 = qJ(4) + t24;
t16 = sin(t21);
t17 = cos(t21);
t29 = cos(qJ(5));
t18 = t29 * pkin(5) + pkin(4);
t25 = -qJ(6) - pkin(10);
t39 = -t16 * t25 + t17 * t18;
t38 = pkin(3) * t20 + t39;
t56 = t30 * pkin(2) + t38;
t19 = sin(t24);
t35 = t16 * t18 + t17 * t25;
t55 = pkin(3) * t19 + t35;
t28 = sin(qJ(1));
t31 = cos(qJ(1));
t37 = g(1) * t31 + g(2) * t28;
t26 = sin(qJ(5));
t43 = t31 * t26;
t44 = t28 * t29;
t10 = -t17 * t43 + t44;
t48 = g(3) * t16;
t42 = t31 * t29;
t45 = t28 * t26;
t8 = t17 * t45 + t42;
t54 = -g(1) * t10 + g(2) * t8 + t26 * t48;
t3 = -g(3) * t17 + t37 * t16;
t40 = pkin(5) * t26 + pkin(7) + pkin(8) + pkin(9);
t36 = g(1) * t28 - g(2) * t31;
t34 = -pkin(1) - t56;
t27 = sin(qJ(2));
t11 = t17 * t42 + t45;
t9 = -t17 * t44 + t43;
t7 = t36 * t16;
t6 = g(3) * t19 + t37 * t20;
t5 = -g(3) * t20 + t37 * t19;
t4 = t37 * t17 + t48;
t2 = t3 * t29;
t1 = t3 * t26;
t12 = [0, t36, t37, 0, 0, 0, 0, 0, t36 * t30, -t36 * t27, 0, 0, 0, 0, 0, t36 * t20, -t36 * t19, 0, 0, 0, 0, 0, t36 * t17, -t7, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11, -g(1) * t8 - g(2) * t10, t7 (-g(1) * t40 + g(2) * t34) * t31 + (-g(1) * t34 - g(2) * t40) * t28; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t30 + t37 * t27, g(3) * t27 + t37 * t30, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(3) * t56 + t37 * (t27 * pkin(2) + t55); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(3) * t38 + t37 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(3) * t39 + t37 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, g(1) * t11 - g(2) * t9 + t29 * t48, 0, t54 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3;];
taug_reg  = t12;
