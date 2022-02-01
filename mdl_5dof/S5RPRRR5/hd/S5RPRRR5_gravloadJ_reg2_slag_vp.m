% Calculate inertial parameters regressor of gravitation load for
% S5RPRRR5
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
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:49:03
% EndTime: 2022-01-20 09:49:03
% DurationCPUTime: 0.18s
% Computational Cost: add. (224->42), mult. (140->49), div. (0->0), fcn. (125->10), ass. (0->32)
t23 = qJ(1) + pkin(9);
t19 = qJ(3) + t23;
t14 = sin(t19);
t15 = cos(t19);
t37 = t15 * pkin(3) + t14 * pkin(7);
t18 = cos(t23);
t28 = cos(qJ(1));
t36 = t28 * pkin(1) + pkin(2) * t18;
t35 = -t14 * pkin(3) + t15 * pkin(7);
t27 = cos(qJ(4));
t16 = t27 * pkin(4) + pkin(3);
t29 = -pkin(8) - pkin(7);
t34 = -t14 * t29 + t15 * t16;
t17 = sin(t23);
t26 = sin(qJ(1));
t33 = -t26 * pkin(1) - pkin(2) * t17;
t8 = g(1) * t15 + g(2) * t14;
t7 = g(1) * t14 - g(2) * t15;
t32 = g(1) * t26 - g(2) * t28;
t31 = -t14 * t16 - t15 * t29;
t25 = sin(qJ(4));
t30 = -g(3) * t27 + t8 * t25;
t24 = qJ(4) + qJ(5);
t21 = cos(t24);
t20 = sin(t24);
t6 = t7 * t27;
t5 = t7 * t25;
t4 = t7 * t21;
t3 = t7 * t20;
t2 = g(3) * t20 + t8 * t21;
t1 = -g(3) * t21 + t8 * t20;
t9 = [0, 0, 0, 0, 0, 0, t32, g(1) * t28 + g(2) * t26, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t17 - g(2) * t18, g(1) * t18 + g(2) * t17, 0, t32 * pkin(1), 0, 0, 0, 0, 0, 0, t7, t8, 0, -g(1) * t33 - g(2) * t36, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t33 + t35) - g(2) * (t36 + t37), 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * (t31 + t33) - g(2) * (t34 + t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * t35 - g(2) * t37, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * t31 - g(2) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, g(3) * t25 + t27 * t8, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t30 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t9;
