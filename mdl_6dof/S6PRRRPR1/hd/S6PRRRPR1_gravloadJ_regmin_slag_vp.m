% Calculate minimal parameter regressor of gravitation load for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t30 = sin(pkin(6));
t53 = g(3) * t30;
t28 = qJ(3) + qJ(4);
t23 = pkin(12) + t28;
t21 = cos(t23);
t33 = sin(qJ(6));
t52 = t21 * t33;
t36 = cos(qJ(6));
t51 = t21 * t36;
t29 = sin(pkin(11));
t50 = t29 * t30;
t31 = cos(pkin(11));
t49 = t30 * t31;
t34 = sin(qJ(3));
t48 = t30 * t34;
t35 = sin(qJ(2));
t47 = t30 * t35;
t37 = cos(qJ(3));
t46 = t30 * t37;
t32 = cos(pkin(6));
t45 = t32 * t35;
t38 = cos(qJ(2));
t44 = t32 * t38;
t43 = t33 * t38;
t42 = t36 * t38;
t25 = cos(t28);
t17 = t37 * pkin(3) + pkin(4) * t25;
t12 = t29 * t38 + t31 * t45;
t14 = -t29 * t45 + t31 * t38;
t20 = sin(t23);
t41 = g(1) * (-t14 * t20 + t21 * t50) + g(2) * (-t12 * t20 - t21 * t49) + g(3) * (-t20 * t47 + t32 * t21);
t11 = t29 * t35 - t31 * t44;
t13 = t29 * t44 + t31 * t35;
t40 = -g(1) * t13 - g(2) * t11 + t38 * t53;
t39 = g(1) * t14 + g(2) * t12 + g(3) * t47;
t24 = sin(t28);
t3 = -g(1) * (-t14 * t24 + t25 * t50) - g(2) * (-t12 * t24 - t25 * t49) - g(3) * (-t24 * t47 + t32 * t25);
t27 = -qJ(5) - pkin(9) - pkin(8);
t16 = -t34 * pkin(3) - pkin(4) * t24;
t15 = pkin(2) + t17;
t10 = t32 * t20 + t21 * t47;
t8 = t14 * t21 + t20 * t50;
t6 = t12 * t21 - t20 * t49;
t4 = -g(1) * (-t14 * t25 - t24 * t50) - g(2) * (-t12 * t25 + t24 * t49) - g(3) * (-t32 * t24 - t25 * t47);
t2 = t41 * t36;
t1 = t41 * t33;
t5 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t40, t39, 0, 0, 0, 0, 0, -t40 * t37, t40 * t34, 0, 0, 0, 0, 0, -t40 * t25, t40 * t24, -t39, -g(1) * (-t13 * t15 - t14 * t27) - g(2) * (-t11 * t15 - t12 * t27) - (t15 * t38 - t27 * t35) * t53, 0, 0, 0, 0, 0, -g(1) * (-t13 * t51 + t14 * t33) - g(2) * (-t11 * t51 + t12 * t33) - (t21 * t42 + t33 * t35) * t53, -g(1) * (t13 * t52 + t14 * t36) - g(2) * (t11 * t52 + t12 * t36) - (-t21 * t43 + t35 * t36) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t34 + t29 * t46) - g(2) * (-t12 * t34 - t31 * t46) - g(3) * (t32 * t37 - t34 * t47) -g(1) * (-t14 * t37 - t29 * t48) - g(2) * (-t12 * t37 + t31 * t48) - g(3) * (-t32 * t34 - t35 * t46) 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * (t14 * t16 + t17 * t50) - g(2) * (t12 * t16 - t17 * t49) - g(3) * (t16 * t47 + t32 * t17) 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t3 * pkin(4), 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t13 * t36 - t8 * t33) - g(2) * (t11 * t36 - t6 * t33) - g(3) * (-t10 * t33 - t30 * t42) -g(1) * (-t13 * t33 - t8 * t36) - g(2) * (-t11 * t33 - t6 * t36) - g(3) * (-t10 * t36 + t30 * t43);];
taug_reg  = t5;
