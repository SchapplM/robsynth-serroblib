% Calculate inertial parameters regressor of gravitation load for
% S5RPRRP2
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
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:28:12
% EndTime: 2022-01-23 09:28:13
% DurationCPUTime: 0.14s
% Computational Cost: add. (197->40), mult. (132->43), div. (0->0), fcn. (117->8), ass. (0->26)
t19 = qJ(1) + pkin(8);
t17 = qJ(3) + t19;
t12 = sin(t17);
t13 = cos(t17);
t31 = t13 * pkin(3) + t12 * pkin(7);
t16 = cos(t19);
t24 = cos(qJ(1));
t30 = t24 * pkin(1) + pkin(2) * t16;
t29 = -t12 * pkin(3) + t13 * pkin(7);
t23 = cos(qJ(4));
t14 = t23 * pkin(4) + pkin(3);
t20 = -qJ(5) - pkin(7);
t28 = -t12 * t20 + t13 * t14;
t15 = sin(t19);
t22 = sin(qJ(1));
t27 = -t22 * pkin(1) - pkin(2) * t15;
t6 = g(1) * t13 + g(2) * t12;
t5 = g(1) * t12 - g(2) * t13;
t26 = g(1) * t22 - g(2) * t24;
t25 = -t12 * t14 - t13 * t20;
t21 = sin(qJ(4));
t1 = -g(3) * t23 + t6 * t21;
t4 = t5 * t23;
t3 = t5 * t21;
t2 = g(3) * t21 + t6 * t23;
t7 = [0, 0, 0, 0, 0, 0, t26, g(1) * t24 + g(2) * t22, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t15 - g(2) * t16, g(1) * t16 + g(2) * t15, 0, t26 * pkin(1), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(1) * t27 - g(2) * t30, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t27 + t29) - g(2) * (t30 + t31), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t25 + t27) - g(2) * (t28 + t30); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t29 - g(2) * t31, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * t25 - g(2) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t7;
