% Calculate inertial parameters regressor of gravitation load for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:51
% EndTime: 2019-12-05 15:12:51
% DurationCPUTime: 0.14s
% Computational Cost: add. (164->36), mult. (139->47), div. (0->0), fcn. (134->8), ass. (0->29)
t15 = pkin(9) + qJ(3);
t14 = qJ(4) + t15;
t10 = sin(t14);
t11 = cos(t14);
t16 = sin(pkin(8));
t17 = cos(pkin(8));
t21 = g(1) * t17 + g(2) * t16;
t3 = -g(3) * t11 + t21 * t10;
t31 = t11 * pkin(4) + t10 * pkin(7);
t30 = pkin(4) * t10;
t29 = pkin(7) * t11;
t28 = g(3) * t10;
t18 = sin(qJ(5));
t26 = t16 * t18;
t19 = cos(qJ(5));
t25 = t16 * t19;
t24 = t17 * t18;
t23 = t17 * t19;
t12 = sin(t15);
t22 = -pkin(3) * t12 - t30;
t13 = cos(t15);
t20 = -g(3) * t13 + t21 * t12;
t7 = -g(1) * t16 + g(2) * t17;
t6 = t17 * t29;
t5 = t16 * t29;
t4 = t21 * t11 + t28;
t2 = t3 * t19;
t1 = t3 * t18;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, g(3) * t12 + t21 * t13, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t20 * pkin(3), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t22 * t17 + t6) - g(2) * (t22 * t16 + t5) - g(3) * (pkin(3) * t13 + t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-t17 * t30 + t6) - g(2) * (-t16 * t30 + t5) - g(3) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t11 * t24 + t25) - g(2) * (-t11 * t26 - t23) + t18 * t28, -g(1) * (-t11 * t23 - t26) - g(2) * (-t11 * t25 + t24) + t19 * t28, 0, 0;];
taug_reg = t8;
