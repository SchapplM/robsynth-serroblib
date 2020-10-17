% Calculate inertial parameters regressor of gravitation load for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:19
% EndTime: 2019-12-05 17:06:19
% DurationCPUTime: 0.14s
% Computational Cost: add. (230->34), mult. (102->36), div. (0->0), fcn. (90->8), ass. (0->27)
t20 = pkin(9) + qJ(2);
t19 = qJ(3) + t20;
t16 = qJ(4) + t19;
t11 = sin(t16);
t12 = cos(t16);
t29 = t12 * pkin(4) + t11 * pkin(8);
t17 = sin(t20);
t28 = pkin(2) * t17;
t14 = sin(t19);
t27 = pkin(3) * t14;
t15 = cos(t19);
t10 = pkin(3) * t15;
t26 = t10 + t29;
t25 = -t11 * pkin(4) + t12 * pkin(8);
t4 = g(1) * t12 + g(2) * t11;
t3 = g(1) * t11 - g(2) * t12;
t5 = g(1) * t14 - g(2) * t15;
t18 = cos(t20);
t24 = g(1) * t17 - g(2) * t18;
t23 = t25 - t27;
t22 = cos(qJ(5));
t21 = sin(qJ(5));
t13 = pkin(2) * t18;
t6 = g(1) * t15 + g(2) * t14;
t2 = t3 * t22;
t1 = t3 * t21;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, g(1) * t18 + g(2) * t17, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t24 * pkin(2), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * (-t27 - t28) - g(2) * (t10 + t13), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t23 - t28) - g(2) * (t13 + t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(3), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * t23 - g(2) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * t25 - g(2) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t22 + t4 * t21, g(3) * t21 + t4 * t22, 0, 0;];
taug_reg = t7;
