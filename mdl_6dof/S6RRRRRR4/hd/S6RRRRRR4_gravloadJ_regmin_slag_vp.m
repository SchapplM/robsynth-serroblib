% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x38]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t33 = sin(qJ(2));
t36 = cos(qJ(2));
t34 = sin(qJ(1));
t37 = cos(qJ(1));
t41 = g(1) * t37 + g(2) * t34;
t39 = -g(3) * t36 + t41 * t33;
t52 = g(3) * t33;
t50 = t34 * t36;
t31 = qJ(3) + qJ(4);
t30 = qJ(5) + t31;
t27 = qJ(6) + t30;
t23 = sin(t27);
t49 = t37 * t23;
t24 = cos(t27);
t48 = t37 * t24;
t25 = sin(t30);
t47 = t37 * t25;
t26 = cos(t30);
t46 = t37 * t26;
t28 = sin(t31);
t45 = t37 * t28;
t29 = cos(t31);
t44 = t37 * t29;
t32 = sin(qJ(3));
t43 = t37 * t32;
t35 = cos(qJ(3));
t42 = t37 * t35;
t40 = g(1) * t34 - g(2) * t37;
t22 = t34 * t32 + t36 * t42;
t21 = t34 * t35 - t36 * t43;
t20 = -t35 * t50 + t43;
t19 = t32 * t50 + t42;
t18 = t34 * t28 + t36 * t44;
t17 = t34 * t29 - t36 * t45;
t16 = -t29 * t50 + t45;
t15 = t28 * t50 + t44;
t14 = t34 * t25 + t36 * t46;
t13 = t34 * t26 - t36 * t47;
t12 = -t26 * t50 + t47;
t11 = t25 * t50 + t46;
t10 = t34 * t23 + t36 * t48;
t9 = t34 * t24 - t36 * t49;
t8 = -t24 * t50 + t49;
t7 = t23 * t50 + t48;
t6 = g(1) * t18 - g(2) * t16 + t29 * t52;
t5 = -g(1) * t17 + g(2) * t15 + t28 * t52;
t4 = g(1) * t14 - g(2) * t12 + t26 * t52;
t3 = -g(1) * t13 + g(2) * t11 + t25 * t52;
t2 = g(1) * t10 - g(2) * t8 + t24 * t52;
t1 = -g(1) * t9 + g(2) * t7 + t23 * t52;
t38 = [0, t40, t41, 0, 0, 0, 0, 0, t40 * t36, -t40 * t33, 0, 0, 0, 0, 0, -g(1) * t20 - g(2) * t22, -g(1) * t19 - g(2) * t21, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, -g(1) * t15 - g(2) * t17, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t14, -g(1) * t11 - g(2) * t13, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9; 0, 0, 0, 0, 0, 0, 0, 0, t39, t41 * t36 + t52, 0, 0, 0, 0, 0, t39 * t35, -t39 * t32, 0, 0, 0, 0, 0, t39 * t29, -t39 * t28, 0, 0, 0, 0, 0, t39 * t26, -t39 * t25, 0, 0, 0, 0, 0, t39 * t24, -t39 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t21 + g(2) * t19 + t32 * t52, g(1) * t22 - g(2) * t20 + t35 * t52, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t38;
