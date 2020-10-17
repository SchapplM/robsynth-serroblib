% Calculate inertial parameters regressor of gravitation load for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:15
% EndTime: 2019-12-05 17:38:16
% DurationCPUTime: 0.13s
% Computational Cost: add. (75->34), mult. (122->39), div. (0->0), fcn. (111->6), ass. (0->20)
t13 = sin(qJ(4));
t23 = pkin(4) * t13;
t14 = sin(qJ(1));
t16 = cos(qJ(1));
t22 = t16 * pkin(1) + t14 * qJ(2);
t21 = -pkin(1) - qJ(3);
t20 = t16 * qJ(3) + t22;
t4 = g(1) * t16 + g(2) * t14;
t3 = g(1) * t14 - g(2) * t16;
t9 = t16 * qJ(2);
t19 = t21 * t14 + t9;
t15 = cos(qJ(4));
t18 = g(3) * t13 - t4 * t15;
t17 = -pkin(7) - pkin(6);
t12 = qJ(4) + qJ(5);
t6 = cos(t12);
t5 = sin(t12);
t2 = g(3) * t5 - t4 * t6;
t1 = g(3) * t6 + t4 * t5;
t7 = [0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * (-t14 * pkin(1) + t9) - g(2) * t22, 0, 0, 0, 0, 0, 0, 0, -t4, t3, -g(1) * t19 - g(2) * t20, 0, 0, 0, 0, 0, 0, t3 * t13, t3 * t15, t4, -g(1) * (-t16 * pkin(6) + t19) - g(2) * (-t14 * pkin(6) + t20), 0, 0, 0, 0, 0, 0, t3 * t5, t3 * t6, t4, -g(1) * (t16 * t17 + t9) - g(2) * (t16 * t23 + t20) + (-g(1) * (t21 - t23) - g(2) * t17) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, g(3) * t15 + t4 * t13, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t18 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t7;
