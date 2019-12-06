% Calculate minimal parameter regressor of gravitation load for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [5x20]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t18 = sin(pkin(9));
t23 = sin(qJ(2));
t26 = cos(qJ(2));
t36 = cos(pkin(9));
t37 = cos(pkin(5));
t32 = t37 * t36;
t6 = t18 * t23 - t26 * t32;
t46 = g(2) * t6;
t7 = t18 * t26 + t23 * t32;
t45 = g(2) * t7;
t19 = sin(pkin(5));
t44 = g(3) * t19;
t43 = t19 * t23;
t25 = cos(qJ(3));
t42 = t19 * t25;
t41 = t19 * t26;
t21 = sin(qJ(4));
t40 = t21 * t25;
t24 = cos(qJ(4));
t39 = t24 * t25;
t38 = t25 * t26;
t35 = pkin(4) * t21 + pkin(7);
t34 = t18 * t37;
t33 = t19 * t36;
t17 = t24 * pkin(4) + pkin(3);
t20 = -qJ(5) - pkin(8);
t22 = sin(qJ(3));
t31 = t17 * t25 - t20 * t22 + pkin(2);
t10 = t22 * t43 - t37 * t25;
t2 = t7 * t22 + t25 * t33;
t9 = -t23 * t34 + t36 * t26;
t4 = -t18 * t42 + t9 * t22;
t30 = g(1) * t4 + g(2) * t2 + g(3) * t10;
t11 = t37 * t22 + t23 * t42;
t3 = -t22 * t33 + t7 * t25;
t5 = t18 * t19 * t22 + t9 * t25;
t29 = g(1) * t5 + g(2) * t3 + g(3) * t11;
t8 = t36 * t23 + t26 * t34;
t28 = -g(1) * t8 + g(3) * t41 - t46;
t27 = -g(1) * (-t5 * t21 + t8 * t24) - g(2) * (-t3 * t21 + t6 * t24) - g(3) * (-t11 * t21 - t24 * t41);
t1 = t28 * t22;
t12 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -t28, g(1) * t9 + g(3) * t43 + t45, 0, 0, 0, 0, 0, -t28 * t25, t1, 0, 0, 0, 0, 0, -g(1) * (t9 * t21 - t8 * t39) - g(2) * (t7 * t21 - t6 * t39) - (t21 * t23 + t24 * t38) * t44, -g(1) * (t9 * t24 + t8 * t40) - g(2) * (t7 * t24 + t6 * t40) - (-t21 * t38 + t23 * t24) * t44, -t1, -g(1) * (-t31 * t8 + t35 * t9) - t35 * t45 + t31 * t46 - (t35 * t23 + t31 * t26) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t29, 0, 0, 0, 0, 0, t30 * t24, -t30 * t21, -t29, -g(1) * (-t4 * t17 - t5 * t20) - g(2) * (-t2 * t17 - t3 * t20) - g(3) * (-t10 * t17 - t11 * t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -g(1) * (-t8 * t21 - t5 * t24) - g(2) * (-t6 * t21 - t3 * t24) - g(3) * (-t11 * t24 + t21 * t41), 0, t27 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30;];
taug_reg = t12;
