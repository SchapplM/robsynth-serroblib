% Calculate inertial parameters regressor of gravitation load for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:04
% EndTime: 2019-12-05 15:43:04
% DurationCPUTime: 0.12s
% Computational Cost: add. (155->36), mult. (102->39), div. (0->0), fcn. (93->8), ass. (0->20)
t20 = cos(pkin(9));
t9 = t20 * pkin(3) + pkin(2);
t21 = -pkin(6) - qJ(3);
t17 = pkin(9) + qJ(4);
t18 = pkin(8) + qJ(2);
t11 = sin(t18);
t13 = cos(t18);
t4 = g(1) * t13 + g(2) * t11;
t3 = g(1) * t11 - g(2) * t13;
t10 = sin(t17);
t12 = cos(t17);
t22 = -g(3) * t12 + t4 * t10;
t16 = -pkin(7) + t21;
t14 = qJ(5) + t17;
t8 = cos(t14);
t7 = sin(t14);
t5 = pkin(4) * t12 + t9;
t2 = g(3) * t7 + t4 * t8;
t1 = -g(3) * t8 + t4 * t7;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t3 * t20, -t3 * sin(pkin(9)), -t4, -g(1) * (-t11 * pkin(2) + t13 * qJ(3)) - g(2) * (t13 * pkin(2) + t11 * qJ(3)), 0, 0, 0, 0, 0, 0, t3 * t12, -t3 * t10, -t4, -g(1) * (-t11 * t9 - t13 * t21) - g(2) * (-t11 * t21 + t13 * t9), 0, 0, 0, 0, 0, 0, t3 * t8, -t3 * t7, -t4, -g(1) * (-t11 * t5 - t13 * t16) - g(2) * (-t11 * t16 + t13 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, g(3) * t10 + t4 * t12, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t22 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t6;
