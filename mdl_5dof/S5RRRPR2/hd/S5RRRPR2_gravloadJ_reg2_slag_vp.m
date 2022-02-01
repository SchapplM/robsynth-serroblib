% Calculate inertial parameters regressor of gravitation load for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:53
% EndTime: 2022-01-20 11:30:53
% DurationCPUTime: 0.15s
% Computational Cost: add. (274->41), mult. (138->47), div. (0->0), fcn. (120->10), ass. (0->35)
t23 = qJ(1) + qJ(2);
t19 = sin(t23);
t37 = pkin(2) * t19;
t21 = qJ(3) + t23;
t17 = sin(t21);
t36 = pkin(3) * t17;
t25 = sin(qJ(1));
t35 = t25 * pkin(1);
t20 = cos(t23);
t15 = pkin(2) * t20;
t27 = cos(qJ(1));
t22 = t27 * pkin(1);
t34 = t15 + t22;
t16 = pkin(9) + t21;
t12 = sin(t16);
t13 = cos(t16);
t18 = cos(t21);
t14 = pkin(3) * t18;
t33 = t13 * pkin(4) + t12 * pkin(8) + t14;
t32 = t15 + t33;
t31 = -t36 - t37;
t4 = g(1) * t13 + g(2) * t12;
t3 = g(1) * t12 - g(2) * t13;
t5 = g(1) * t17 - g(2) * t18;
t7 = g(1) * t19 - g(2) * t20;
t30 = g(1) * t25 - g(2) * t27;
t29 = -t12 * pkin(4) + t13 * pkin(8) - t36;
t28 = t29 - t37;
t26 = cos(qJ(5));
t24 = sin(qJ(5));
t8 = g(1) * t20 + g(2) * t19;
t6 = g(1) * t18 + g(2) * t17;
t2 = t3 * t26;
t1 = t3 * t24;
t9 = [0, 0, 0, 0, 0, 0, t30, g(1) * t27 + g(2) * t25, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t30 * pkin(1), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(1) * (-t35 - t37) - g(2) * t34, 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * (t31 - t35) - g(2) * (t14 + t34), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t28 - t35) - g(2) * (t22 + t32); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t7 * pkin(2), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * t31 - g(2) * (t14 + t15), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * t28 - g(2) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(3), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * t29 - g(2) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t26 + t4 * t24, g(3) * t24 + t4 * t26, 0, 0;];
taug_reg = t9;
