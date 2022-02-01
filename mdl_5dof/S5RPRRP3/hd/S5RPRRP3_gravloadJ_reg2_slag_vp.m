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
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:30:19
% EndTime: 2022-01-23 09:30:19
% DurationCPUTime: 0.15s
% Computational Cost: add. (177->43), mult. (151->49), div. (0->0), fcn. (134->8), ass. (0->27)
t24 = -pkin(7) - pkin(6);
t21 = sin(qJ(1));
t28 = t21 * pkin(1);
t19 = qJ(3) + qJ(4);
t14 = cos(t19);
t22 = cos(qJ(3));
t15 = t22 * pkin(3);
t27 = pkin(4) * t14 + t15;
t18 = qJ(1) + pkin(8);
t11 = sin(t18);
t12 = cos(t18);
t6 = g(1) * t12 + g(2) * t11;
t5 = g(1) * t11 - g(2) * t12;
t23 = cos(qJ(1));
t26 = g(1) * t21 - g(2) * t23;
t13 = sin(t19);
t1 = -g(3) * t14 + t6 * t13;
t20 = sin(qJ(3));
t25 = -g(3) * t22 + t6 * t20;
t17 = -qJ(5) + t24;
t16 = t23 * pkin(1);
t10 = t15 + pkin(2);
t7 = pkin(2) + t27;
t4 = t5 * t14;
t3 = t5 * t13;
t2 = g(3) * t13 + t6 * t14;
t8 = [0, 0, 0, 0, 0, 0, t26, g(1) * t23 + g(2) * t21, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t26 * pkin(1), 0, 0, 0, 0, 0, 0, t5 * t22, -t5 * t20, -t6, -g(1) * (-t11 * pkin(2) + t12 * pkin(6) - t28) - g(2) * (t12 * pkin(2) + t11 * pkin(6) + t16), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-t11 * t10 - t12 * t24 - t28) - g(2) * (t12 * t10 - t11 * t24 + t16), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-t11 * t7 - t12 * t17 - t28) - g(2) * (-t11 * t17 + t12 * t7 + t16); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, g(3) * t20 + t6 * t22, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t25 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t27 - t6 * (-t20 * pkin(3) - pkin(4) * t13); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t8;
