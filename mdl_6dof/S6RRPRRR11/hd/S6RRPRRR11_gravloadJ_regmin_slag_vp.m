% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR11_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR11_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t31 = sin(qJ(1));
t34 = cos(qJ(1));
t21 = g(1) * t34 + g(2) * t31;
t30 = sin(qJ(2));
t33 = cos(qJ(2));
t14 = g(3) * t30 + t21 * t33;
t52 = g(3) * t33;
t28 = qJ(4) + qJ(5);
t26 = qJ(6) + t28;
t22 = sin(t26);
t51 = t31 * t22;
t23 = cos(t26);
t50 = t31 * t23;
t24 = sin(t28);
t49 = t31 * t24;
t25 = cos(t28);
t48 = t31 * t25;
t29 = sin(qJ(4));
t47 = t31 * t29;
t32 = cos(qJ(4));
t46 = t31 * t32;
t45 = t34 * t22;
t44 = t34 * t23;
t43 = t34 * t24;
t42 = t34 * t25;
t41 = t34 * t29;
t40 = t34 * t32;
t39 = g(1) * t31 - g(2) * t34;
t38 = t33 * pkin(2) + t30 * qJ(3);
t36 = pkin(1) + t38;
t20 = t39 * t33;
t19 = t39 * t30;
t18 = -t30 * t47 + t40;
t17 = t30 * t46 + t41;
t16 = t30 * t41 + t46;
t15 = t30 * t40 - t47;
t13 = t21 * t30 - t52;
t12 = -t30 * t49 + t42;
t11 = t30 * t48 + t43;
t10 = t30 * t43 + t48;
t9 = t30 * t42 - t49;
t8 = -t30 * t51 + t44;
t7 = t30 * t50 + t45;
t6 = t30 * t45 + t50;
t5 = t30 * t44 - t51;
t4 = g(1) * t10 - g(2) * t12 - t24 * t52;
t3 = -g(1) * t9 - g(2) * t11 + t25 * t52;
t2 = g(1) * t6 - g(2) * t8 - t22 * t52;
t1 = -g(1) * t5 - g(2) * t7 + t23 * t52;
t27 = [0, t39, t21, 0, 0, 0, 0, 0, t20, -t19, -t21, -t20, t19 (-g(1) * pkin(7) - g(2) * t36) * t34 + (-g(2) * pkin(7) + g(1) * t36) * t31, 0, 0, 0, 0, 0, -g(1) * t18 - g(2) * t16, g(1) * t17 - g(2) * t15, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t6, g(1) * t7 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, -t13, -t14, -g(3) * t38 + t21 * (pkin(2) * t30 - qJ(3) * t33) 0, 0, 0, 0, 0, -t14 * t29, -t14 * t32, 0, 0, 0, 0, 0, -t14 * t24, -t14 * t25, 0, 0, 0, 0, 0, -t14 * t22, -t14 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17 + t32 * t52, g(1) * t16 - g(2) * t18 - t29 * t52, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t27;
