% Calculate inertial parameters regressor of gravitation load for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:50
% EndTime: 2019-12-05 16:15:51
% DurationCPUTime: 0.13s
% Computational Cost: add. (190->38), mult. (102->36), div. (0->0), fcn. (92->8), ass. (0->26)
t22 = pkin(8) + qJ(2);
t17 = sin(t22);
t31 = pkin(2) * t17;
t20 = qJ(3) + t22;
t13 = sin(t20);
t14 = cos(t20);
t30 = t14 * pkin(3) + t13 * qJ(4);
t29 = -t13 * pkin(3) + t14 * qJ(4);
t24 = cos(pkin(9));
t15 = t24 * pkin(4) + pkin(3);
t25 = -pkin(7) - qJ(4);
t28 = -t13 * t25 + t14 * t15;
t6 = g(1) * t14 + g(2) * t13;
t5 = g(1) * t13 - g(2) * t14;
t19 = cos(t22);
t27 = g(1) * t17 - g(2) * t19;
t26 = -t13 * t15 - t14 * t25;
t21 = pkin(9) + qJ(5);
t18 = cos(t21);
t16 = sin(t21);
t12 = pkin(2) * t19;
t4 = t5 * t24;
t3 = t5 * sin(pkin(9));
t2 = t5 * t18;
t1 = t5 * t16;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, g(1) * t19 + g(2) * t17, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t27 * pkin(2), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t29 - t31) - g(2) * (t12 + t30), 0, 0, 0, 0, 0, 0, t2, -t1, -t6, -g(1) * (t26 - t31) - g(2) * (t12 + t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t29 - g(2) * t30, 0, 0, 0, 0, 0, 0, t2, -t1, -t6, -g(1) * t26 - g(2) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t18 + t6 * t16, g(3) * t16 + t6 * t18, 0, 0;];
taug_reg = t7;
