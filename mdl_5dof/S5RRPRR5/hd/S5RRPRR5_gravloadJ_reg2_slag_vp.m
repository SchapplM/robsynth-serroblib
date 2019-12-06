% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR5
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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t31 = sin(qJ(1));
t42 = t31 * pkin(1);
t32 = cos(qJ(1));
t41 = t32 * pkin(1);
t29 = cos(pkin(9));
t17 = t29 * pkin(3) + pkin(2);
t30 = -pkin(7) - qJ(3);
t26 = pkin(9) + qJ(4);
t27 = qJ(1) + qJ(2);
t22 = sin(t27);
t23 = cos(t27);
t40 = -t22 * pkin(2) + t23 * qJ(3);
t20 = cos(t26);
t11 = pkin(4) * t20 + t17;
t25 = -pkin(8) + t30;
t39 = -t23 * t11 + t22 * t25;
t38 = -t23 * t17 + t22 * t30;
t10 = g(2) * t23 + g(3) * t22;
t9 = g(2) * t22 - g(3) * t23;
t37 = g(2) * t32 + g(3) * t31;
t36 = -t23 * pkin(2) - t22 * qJ(3);
t35 = -t22 * t11 - t23 * t25;
t34 = -t22 * t17 - t23 * t30;
t19 = sin(t26);
t33 = -g(1) * t20 - t19 * t9;
t21 = qJ(5) + t26;
t16 = cos(t21);
t15 = sin(t21);
t8 = t10 * t29;
t7 = t10 * sin(pkin(9));
t6 = t10 * t20;
t5 = t10 * t19;
t4 = t10 * t16;
t3 = t10 * t15;
t2 = -g(1) * t16 - t15 * t9;
t1 = g(1) * t15 - t16 * t9;
t12 = [0, 0, 0, 0, 0, 0, t37, -g(2) * t31 + g(3) * t32, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, t37 * pkin(1), 0, 0, 0, 0, 0, 0, t8, -t7, t9, -g(2) * (t36 - t41) - g(3) * (t40 - t42), 0, 0, 0, 0, 0, 0, t6, -t5, t9, -g(2) * (t38 - t41) - g(3) * (t34 - t42), 0, 0, 0, 0, 0, 0, t4, -t3, t9, -g(2) * (t39 - t41) - g(3) * (t35 - t42); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, t9, -g(2) * t36 - g(3) * t40, 0, 0, 0, 0, 0, 0, t6, -t5, t9, -g(2) * t38 - g(3) * t34, 0, 0, 0, 0, 0, 0, t4, -t3, t9, -g(2) * t39 - g(3) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, g(1) * t19 - t20 * t9, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t33 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t12;
