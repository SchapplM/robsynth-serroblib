% Calculate minimal parameter regressor of gravitation load for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t22 = cos(pkin(7));
t18 = sin(pkin(8));
t26 = cos(qJ(1));
t31 = t26 * t18;
t21 = cos(pkin(8));
t24 = sin(qJ(1));
t32 = t24 * t21;
t10 = t22 * t31 - t32;
t23 = sin(qJ(5));
t25 = cos(qJ(5));
t30 = t26 * t21;
t33 = t24 * t18;
t11 = t22 * t30 + t33;
t17 = sin(pkin(9));
t20 = cos(pkin(9));
t19 = sin(pkin(7));
t34 = t19 * t26;
t4 = t11 * t20 + t17 * t34;
t41 = -t10 * t25 + t4 * t23;
t40 = t10 * t23 + t4 * t25;
t39 = g(2) * t26;
t36 = t18 * t19;
t35 = t19 * t24;
t8 = t22 * t33 + t30;
t29 = g(2) * t10 + g(3) * t8;
t14 = g(3) * t24 + t39;
t13 = g(2) * t24 - g(3) * t26;
t28 = pkin(2) * t22 + qJ(3) * t19 + pkin(1);
t27 = (g(2) * qJ(2) + g(3) * t28) * t24 + t28 * t39;
t16 = t26 * qJ(2);
t12 = t14 * t19;
t9 = -t22 * t32 + t31;
t7 = t19 * t21 * t20 - t22 * t17;
t6 = g(1) * t22 + t13 * t19;
t3 = -t17 * t35 + t9 * t20;
t2 = -t8 * t23 + t3 * t25;
t1 = -t3 * t23 - t8 * t25;
t5 = [0, t14, -t13, t14 * t22, -t12, t13, -g(2) * (-t26 * pkin(1) - t24 * qJ(2)) - g(3) * (-t24 * pkin(1) + t16), g(2) * t11 - g(3) * t9, -t29, t12, -g(3) * t16 + t27, g(2) * t4 - g(3) * t3, -g(2) * (t11 * t17 - t20 * t34) - g(3) * (-t9 * t17 - t20 * t35), t29, -g(2) * (-t11 * pkin(3) - t10 * qJ(4)) - g(3) * (t9 * pkin(3) - t8 * qJ(4) + t16) + t27, 0, 0, 0, 0, 0, g(2) * t40 - g(3) * t2, -g(2) * t41 - g(3) * t1; 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, -t14, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t36 + g(2) * t8 - g(3) * t10, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t7 * t23 + t25 * t36) - g(2) * t1 + g(3) * t41, -g(1) * (-t23 * t36 - t7 * t25) + g(2) * t2 + g(3) * t40;];
taug_reg = t5;
