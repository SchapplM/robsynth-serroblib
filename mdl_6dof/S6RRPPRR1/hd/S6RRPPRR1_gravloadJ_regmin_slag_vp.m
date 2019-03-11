% Calculate minimal parameter regressor of gravitation load for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t24 = sin(qJ(6));
t22 = qJ(2) + pkin(10);
t18 = sin(t22);
t19 = cos(t22);
t25 = sin(qJ(5));
t41 = cos(qJ(5));
t35 = t18 * t25 + t19 * t41;
t27 = sin(qJ(1));
t45 = t18 * t41 - t19 * t25;
t4 = t45 * t27;
t30 = cos(qJ(1));
t6 = t45 * t30;
t34 = -g(1) * t6 - g(2) * t4 + g(3) * t35;
t47 = t34 * t24;
t28 = cos(qJ(6));
t46 = t34 * t28;
t13 = g(1) * t30 + g(2) * t27;
t44 = g(3) * t45;
t12 = g(1) * t27 - g(2) * t30;
t5 = t35 * t27;
t38 = t5 * t24 - t30 * t28;
t37 = t30 * t24 + t5 * t28;
t36 = t19 * pkin(3) + t18 * qJ(4);
t7 = t35 * t30;
t33 = g(1) * t7 + g(2) * t5 + t44;
t26 = sin(qJ(2));
t29 = cos(qJ(2));
t31 = -g(3) * t29 + t13 * t26;
t23 = -qJ(3) - pkin(7);
t20 = t29 * pkin(2);
t17 = t20 + pkin(1);
t14 = t30 * t17;
t3 = -g(3) * t19 + t13 * t18;
t2 = -t27 * t24 + t7 * t28;
t1 = -t7 * t24 - t27 * t28;
t8 = [0, t12, t13, 0, 0, 0, 0, 0, t12 * t29, -t12 * t26, -t13, -g(1) * (-t27 * t17 - t30 * t23) - g(2) * (-t27 * t23 + t14) t12 * t19, -t13, t12 * t18, -g(2) * t14 + (g(1) * t23 - g(2) * t36) * t30 + (-g(1) * (-t17 - t36) + g(2) * t23) * t27, 0, 0, 0, 0, 0, g(1) * t5 - g(2) * t7, g(1) * t4 - g(2) * t6, 0, 0, 0, 0, 0, g(1) * t37 - g(2) * t2, -g(1) * t38 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, t31, g(3) * t26 + t13 * t29, 0, t31 * pkin(2), t3, 0, -g(3) * t18 - t13 * t19, -g(3) * (t20 + t36) + t13 * (pkin(2) * t26 + pkin(3) * t18 - qJ(4) * t19) 0, 0, 0, 0, 0, -t34, -t33, 0, 0, 0, 0, 0, -t46, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, t33, 0, 0, 0, 0, 0, t46, -t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t38 + t24 * t44, g(1) * t2 + g(2) * t37 + t28 * t44;];
taug_reg  = t8;
