% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t23 = qJ(1) + qJ(2);
t19 = sin(t23);
t38 = pkin(2) * t19;
t21 = cos(t23);
t37 = pkin(2) * t21;
t25 = sin(qJ(1));
t36 = t25 * pkin(1);
t27 = cos(qJ(1));
t35 = t27 * pkin(1);
t17 = pkin(9) + t23;
t14 = sin(t17);
t15 = cos(t17);
t8 = g(2) * t15 + g(3) * t14;
t7 = g(2) * t14 - g(3) * t15;
t10 = g(2) * t21 + g(3) * t19;
t34 = g(2) * t27 + g(3) * t25;
t33 = -t14 * pkin(3) + t15 * pkin(7) - t38;
t26 = cos(qJ(4));
t16 = t26 * pkin(4) + pkin(3);
t28 = -pkin(8) - pkin(7);
t32 = t14 * t28 - t15 * t16 - t37;
t31 = -t15 * pkin(3) - t14 * pkin(7) - t37;
t30 = -t14 * t16 - t15 * t28 - t38;
t24 = sin(qJ(4));
t29 = -g(1) * t26 - t7 * t24;
t22 = qJ(4) + qJ(5);
t20 = cos(t22);
t18 = sin(t22);
t9 = -g(2) * t19 + g(3) * t21;
t6 = t8 * t26;
t5 = t8 * t24;
t4 = t8 * t20;
t3 = t8 * t18;
t2 = -g(1) * t20 - t7 * t18;
t1 = g(1) * t18 - t7 * t20;
t11 = [0, 0, 0, 0, 0, 0, t34, -g(2) * t25 + g(3) * t27, 0, 0, 0, 0, 0, 0, 0, 0, t10, t9, 0, t34 * pkin(1), 0, 0, 0, 0, 0, 0, t8, -t7, 0, -g(2) * (-t35 - t37) - g(3) * (-t36 - t38), 0, 0, 0, 0, 0, 0, t6, -t5, t7, -g(2) * (t31 - t35) - g(3) * (t33 - t36), 0, 0, 0, 0, 0, 0, t4, -t3, t7, -g(2) * (t32 - t35) - g(3) * (t30 - t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t9, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t10 * pkin(2), 0, 0, 0, 0, 0, 0, t6, -t5, t7, -g(2) * t31 - g(3) * t33, 0, 0, 0, 0, 0, 0, t4, -t3, t7, -g(2) * t32 - g(3) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, g(1) * t24 - t7 * t26, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t29 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t11;
