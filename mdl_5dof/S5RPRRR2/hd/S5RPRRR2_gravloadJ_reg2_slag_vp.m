% Calculate inertial parameters regressor of gravitation load for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:14
% EndTime: 2019-12-05 18:12:14
% DurationCPUTime: 0.15s
% Computational Cost: add. (221->44), mult. (175->57), div. (0->0), fcn. (158->10), ass. (0->28)
t27 = cos(pkin(9));
t16 = t27 * pkin(2) + pkin(1);
t28 = -pkin(6) - qJ(2);
t25 = pkin(9) + qJ(3);
t24 = -pkin(7) + t28;
t19 = cos(t25);
t13 = pkin(3) * t19;
t7 = t13 + t16;
t20 = qJ(4) + t25;
t29 = sin(qJ(1));
t30 = cos(qJ(1));
t9 = g(1) * t30 + g(2) * t29;
t8 = g(1) * t29 - g(2) * t30;
t14 = sin(t20);
t15 = cos(t20);
t3 = -g(3) * t15 + t9 * t14;
t18 = sin(t25);
t31 = -g(3) * t19 + t9 * t18;
t21 = -pkin(8) + t24;
t17 = qJ(5) + t20;
t12 = cos(t17);
t11 = sin(t17);
t10 = pkin(4) * t15;
t5 = t10 + t7;
t4 = g(3) * t14 + t9 * t15;
t2 = g(3) * t11 + t9 * t12;
t1 = -g(3) * t12 + t9 * t11;
t6 = [0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, t8 * t27, -t8 * sin(pkin(9)), -t9, -g(1) * (-t29 * pkin(1) + t30 * qJ(2)) - g(2) * (t30 * pkin(1) + t29 * qJ(2)), 0, 0, 0, 0, 0, 0, t8 * t19, -t8 * t18, -t9, -g(1) * (-t29 * t16 - t30 * t28) - g(2) * (t30 * t16 - t29 * t28), 0, 0, 0, 0, 0, 0, t8 * t15, -t8 * t14, -t9, -g(1) * (-t30 * t24 - t29 * t7) - g(2) * (-t29 * t24 + t30 * t7), 0, 0, 0, 0, 0, 0, t8 * t12, -t8 * t11, -t9, -g(1) * (-t30 * t21 - t29 * t5) - g(2) * (-t29 * t21 + t30 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, g(3) * t18 + t9 * t19, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t31 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * (t10 + t13) - t9 * (-pkin(3) * t18 - pkin(4) * t14); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t6;
