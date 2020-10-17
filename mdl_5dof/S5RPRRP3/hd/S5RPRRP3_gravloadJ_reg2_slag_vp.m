% Calculate inertial parameters regressor of gravitation load for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:47:43
% EndTime: 2020-01-03 11:47:44
% DurationCPUTime: 0.18s
% Computational Cost: add. (177->43), mult. (151->49), div. (0->0), fcn. (134->8), ass. (0->27)
t26 = -pkin(7) - pkin(6);
t21 = qJ(3) + qJ(4);
t15 = cos(t21);
t24 = cos(qJ(3));
t17 = t24 * pkin(3);
t29 = pkin(4) * t15 + t17;
t20 = qJ(1) + pkin(8);
t12 = sin(t20);
t13 = cos(t20);
t6 = g(2) * t13 + g(3) * t12;
t5 = g(2) * t12 - g(3) * t13;
t23 = sin(qJ(1));
t25 = cos(qJ(1));
t28 = -g(2) * t25 - g(3) * t23;
t14 = sin(t21);
t2 = -g(1) * t15 + t5 * t14;
t22 = sin(qJ(3));
t27 = -g(1) * t24 + t5 * t22;
t19 = -qJ(5) + t26;
t18 = t25 * pkin(1);
t16 = t23 * pkin(1);
t11 = t17 + pkin(2);
t7 = pkin(2) + t29;
t4 = t6 * t15;
t3 = t6 * t14;
t1 = g(1) * t14 + t5 * t15;
t8 = [0, 0, 0, 0, 0, 0, t28, g(2) * t23 - g(3) * t25, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, t28 * pkin(1), 0, 0, 0, 0, 0, 0, -t6 * t24, t6 * t22, -t5, -g(2) * (t13 * pkin(2) + t12 * pkin(6) + t18) - g(3) * (t12 * pkin(2) - t13 * pkin(6) + t16), 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(2) * (t13 * t11 - t12 * t26 + t18) - g(3) * (t12 * t11 + t13 * t26 + t16), 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(2) * (-t12 * t19 + t13 * t7 + t18) - g(3) * (t12 * t7 + t13 * t19 + t16); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, g(1) * t22 + t5 * t24, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t27 * pkin(3), 0, 0, 0, 0, 0, 0, t2, t1, 0, -g(1) * t29 - t5 * (-t22 * pkin(3) - pkin(4) * t14); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6;];
taug_reg = t8;
