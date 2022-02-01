% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:34:20
% EndTime: 2022-01-20 10:34:21
% DurationCPUTime: 0.17s
% Computational Cost: add. (254->40), mult. (128->46), div. (0->0), fcn. (112->10), ass. (0->35)
t23 = qJ(1) + qJ(2);
t20 = sin(t23);
t37 = pkin(2) * t20;
t25 = sin(qJ(1));
t36 = t25 * pkin(1);
t19 = pkin(9) + t23;
t18 = qJ(4) + t19;
t13 = sin(t18);
t14 = cos(t18);
t35 = t14 * pkin(4) + t13 * pkin(8);
t16 = cos(t19);
t12 = pkin(3) * t16;
t21 = cos(t23);
t17 = pkin(2) * t21;
t34 = t12 + t17;
t27 = cos(qJ(1));
t22 = t27 * pkin(1);
t33 = t17 + t22;
t32 = -t13 * pkin(4) + t14 * pkin(8);
t31 = t34 + t35;
t15 = sin(t19);
t30 = -pkin(3) * t15 - t37;
t4 = g(1) * t14 + g(2) * t13;
t3 = g(1) * t13 - g(2) * t14;
t7 = g(1) * t20 - g(2) * t21;
t29 = g(1) * t25 - g(2) * t27;
t28 = t30 + t32;
t26 = cos(qJ(5));
t24 = sin(qJ(5));
t8 = g(1) * t21 + g(2) * t20;
t6 = g(1) * t16 + g(2) * t15;
t5 = g(1) * t15 - g(2) * t16;
t2 = t3 * t26;
t1 = t3 * t24;
t9 = [0, 0, 0, 0, 0, 0, t29, g(1) * t27 + g(2) * t25, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t29 * pkin(1), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(1) * (-t36 - t37) - g(2) * t33, 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * (t30 - t36) - g(2) * (t12 + t33), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t28 - t36) - g(2) * (t22 + t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t7 * pkin(2), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * t30 - g(2) * t34, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * t28 - g(2) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * t32 - g(2) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t26 + t4 * t24, g(3) * t24 + t4 * t26, 0, 0;];
taug_reg = t9;
