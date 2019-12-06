% Calculate inertial parameters regressor of gravitation load for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t20 = sin(qJ(4));
t17 = sin(pkin(8));
t34 = g(1) * t17;
t16 = qJ(1) + pkin(7);
t14 = sin(t16);
t15 = cos(t16);
t22 = cos(qJ(4));
t18 = cos(pkin(8));
t31 = t18 * t20;
t5 = t14 * t31 + t15 * t22;
t7 = -t14 * t22 + t15 * t31;
t1 = -g(2) * t5 + g(3) * t7 + t20 * t34;
t35 = pkin(4) * t20;
t33 = g(2) * t15;
t23 = cos(qJ(1));
t32 = t23 * pkin(1);
t30 = t18 * t22;
t21 = sin(qJ(1));
t28 = -t21 * pkin(1) + t15 * qJ(3);
t11 = g(3) * t14 + t33;
t10 = g(2) * t14 - g(3) * t15;
t27 = g(2) * t23 + g(3) * t21;
t26 = pkin(3) * t18 + pkin(6) * t17 + pkin(2);
t25 = (t22 * pkin(4) + pkin(3)) * t18 - t17 * (-qJ(5) - pkin(6)) + pkin(2);
t24 = g(2) * t32 - g(3) * t28;
t9 = t11 * t17;
t8 = -t14 * t20 - t15 * t30;
t6 = t14 * t30 - t15 * t20;
t4 = -g(2) * t8 + g(3) * t6;
t3 = -g(2) * t7 - g(3) * t5;
t2 = -g(2) * t6 - g(3) * t8 + t22 * t34;
t12 = [0, 0, 0, 0, 0, 0, t27, -g(2) * t21 + g(3) * t23, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, 0, t27 * pkin(1), 0, 0, 0, 0, 0, 0, t11 * t18, -t9, t10, -g(2) * (-t15 * pkin(2) - t14 * qJ(3) - t32) - g(3) * (-t14 * pkin(2) + t28), 0, 0, 0, 0, 0, 0, t4, t3, t9, t26 * t33 + (g(2) * qJ(3) + g(3) * t26) * t14 + t24, 0, 0, 0, 0, 0, 0, t4, t3, t9, (g(2) * t25 - g(3) * t35) * t15 + (-g(2) * (-qJ(3) - t35) + g(3) * t25) * t14 + t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t18 + t10 * t17;];
taug_reg = t12;
