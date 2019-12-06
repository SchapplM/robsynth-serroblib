% Calculate minimal parameter regressor of gravitation load for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t22 = sin(pkin(8));
t41 = g(1) * t22;
t26 = cos(qJ(1));
t40 = g(2) * t26;
t20 = pkin(9) + qJ(4);
t18 = qJ(5) + t20;
t14 = sin(t18);
t25 = sin(qJ(1));
t39 = t25 * t14;
t15 = cos(t18);
t38 = t25 * t15;
t16 = sin(t20);
t37 = t25 * t16;
t17 = cos(t20);
t36 = t25 * t17;
t21 = sin(pkin(9));
t35 = t25 * t21;
t23 = cos(pkin(9));
t34 = t25 * t23;
t33 = t26 * t14;
t32 = t26 * t15;
t31 = t26 * t16;
t30 = t26 * t17;
t29 = t26 * t21;
t28 = t26 * t23;
t13 = g(3) * t25 + t40;
t12 = g(2) * t25 - g(3) * t26;
t24 = cos(pkin(8));
t27 = pkin(2) * t24 + qJ(3) * t22 + pkin(1);
t19 = t26 * qJ(2);
t11 = t13 * t22;
t10 = -t24 * t30 - t37;
t9 = t24 * t31 - t36;
t8 = t24 * t36 - t31;
t7 = t24 * t37 + t30;
t6 = -t24 * t32 - t39;
t5 = t24 * t33 - t38;
t4 = t24 * t38 - t33;
t3 = t24 * t39 + t32;
t2 = -g(2) * t4 - g(3) * t6 + t15 * t41;
t1 = -g(2) * t3 + g(3) * t5 + t14 * t41;
t42 = [0, t13, -t12, t13 * t24, -t11, t12, -g(2) * (-t26 * pkin(1) - t25 * qJ(2)) - g(3) * (-t25 * pkin(1) + t19), -g(2) * (-t24 * t28 - t35) - g(3) * (-t24 * t34 + t29), -g(2) * (t24 * t29 - t34) - g(3) * (t24 * t35 + t28), t11, -g(3) * t19 + t27 * t40 + (g(2) * qJ(2) + g(3) * t27) * t25, 0, 0, 0, 0, 0, -g(2) * t10 + g(3) * t8, -g(2) * t9 - g(3) * t7, 0, 0, 0, 0, 0, -g(2) * t6 + g(3) * t4, -g(2) * t5 - g(3) * t3; 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t24 + t12 * t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t7 + g(3) * t9 + t16 * t41, -g(2) * t8 - g(3) * t10 + t17 * t41, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t42;
