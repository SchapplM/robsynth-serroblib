% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRR3
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
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t33 = qJ(2) + qJ(3);
t28 = sin(t33);
t30 = cos(t33);
t36 = sin(qJ(1));
t39 = cos(qJ(1));
t42 = g(1) * t39 + g(2) * t36;
t11 = -g(3) * t30 + t42 * t28;
t56 = g(3) * t28;
t32 = qJ(4) + qJ(5);
t31 = qJ(6) + t32;
t25 = sin(t31);
t54 = t36 * t25;
t26 = cos(t31);
t53 = t36 * t26;
t27 = sin(t32);
t52 = t36 * t27;
t29 = cos(t32);
t51 = t36 * t29;
t34 = sin(qJ(4));
t50 = t36 * t34;
t37 = cos(qJ(4));
t49 = t36 * t37;
t48 = t39 * t25;
t47 = t39 * t26;
t46 = t39 * t27;
t45 = t39 * t29;
t44 = t39 * t34;
t43 = t39 * t37;
t41 = g(1) * t36 - g(2) * t39;
t38 = cos(qJ(2));
t35 = sin(qJ(2));
t24 = t30 * t43 + t50;
t23 = -t30 * t44 + t49;
t22 = -t30 * t49 + t44;
t21 = t30 * t50 + t43;
t20 = t30 * t45 + t52;
t19 = -t30 * t46 + t51;
t18 = -t30 * t51 + t46;
t17 = t30 * t52 + t45;
t16 = t30 * t47 + t54;
t15 = -t30 * t48 + t53;
t14 = -t30 * t53 + t48;
t13 = t30 * t54 + t47;
t12 = t42 * t30 + t56;
t10 = t11 * t37;
t9 = t11 * t34;
t8 = t11 * t29;
t7 = t11 * t27;
t6 = t11 * t26;
t5 = t11 * t25;
t4 = g(1) * t20 - g(2) * t18 + t29 * t56;
t3 = -g(1) * t19 + g(2) * t17 + t27 * t56;
t2 = g(1) * t16 - g(2) * t14 + t26 * t56;
t1 = -g(1) * t15 + g(2) * t13 + t25 * t56;
t40 = [0, t41, t42, 0, 0, 0, 0, 0, t41 * t38, -t41 * t35, 0, 0, 0, 0, 0, t41 * t30, -t41 * t28, 0, 0, 0, 0, 0, -g(1) * t22 - g(2) * t24, -g(1) * t21 - g(2) * t23, 0, 0, 0, 0, 0, -g(1) * t18 - g(2) * t20, -g(1) * t17 - g(2) * t19, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t38 + t42 * t35, g(3) * t35 + t42 * t38, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, t10, -t9, 0, 0, 0, 0, 0, t8, -t7, 0, 0, 0, 0, 0, t6, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, t10, -t9, 0, 0, 0, 0, 0, t8, -t7, 0, 0, 0, 0, 0, t6, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t23 + g(2) * t21 + t34 * t56, g(1) * t24 - g(2) * t22 + t37 * t56, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t40;
