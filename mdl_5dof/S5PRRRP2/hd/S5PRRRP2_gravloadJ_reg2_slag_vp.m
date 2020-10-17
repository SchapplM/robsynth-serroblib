% Calculate inertial parameters regressor of gravitation load for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:58
% EndTime: 2019-12-05 16:41:59
% DurationCPUTime: 0.16s
% Computational Cost: add. (213->38), mult. (138->37), div. (0->0), fcn. (125->6), ass. (0->25)
t20 = sin(qJ(4));
t21 = cos(qJ(4));
t24 = t21 * pkin(4) + t20 * qJ(5);
t19 = pkin(8) + qJ(2);
t18 = qJ(3) + t19;
t15 = cos(t18);
t14 = sin(t18);
t34 = g(1) * t14;
t5 = -g(2) * t15 + t34;
t6 = g(1) * t15 + g(2) * t14;
t16 = sin(t19);
t35 = pkin(2) * t16;
t29 = t15 * pkin(3) + t14 * pkin(7);
t11 = t15 * pkin(7);
t27 = -t14 * pkin(3) + t11;
t26 = t24 * t15 + t29;
t17 = cos(t19);
t25 = g(1) * t16 - g(2) * t17;
t22 = (-pkin(3) - t24) * t34;
t13 = pkin(2) * t17;
t4 = t5 * t21;
t3 = t5 * t20;
t2 = g(3) * t20 + t6 * t21;
t1 = -g(3) * t21 + t6 * t20;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, g(1) * t17 + g(2) * t16, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t25 * pkin(2), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t27 - t35) - g(2) * (t13 + t29), 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * (t11 - t35) - g(2) * (t13 + t26) - t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t27 - g(2) * t29, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * t11 - g(2) * t26 - t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(3) * t24 + t6 * (pkin(4) * t20 - qJ(5) * t21); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t7;
