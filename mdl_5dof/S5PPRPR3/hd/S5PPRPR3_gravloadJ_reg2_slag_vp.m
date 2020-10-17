% Calculate inertial parameters regressor of gravitation load for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:22
% EndTime: 2019-12-05 15:05:23
% DurationCPUTime: 0.17s
% Computational Cost: add. (118->49), mult. (180->76), div. (0->0), fcn. (202->10), ass. (0->34)
t12 = sin(pkin(8));
t34 = g(3) * t12;
t15 = cos(pkin(7));
t11 = qJ(3) + pkin(9);
t9 = sin(t11);
t33 = t15 * t9;
t16 = sin(qJ(5));
t32 = t12 * t16;
t18 = cos(qJ(5));
t31 = t12 * t18;
t13 = sin(pkin(7));
t14 = cos(pkin(8));
t30 = t13 * t14;
t17 = sin(qJ(3));
t29 = t13 * t17;
t19 = cos(qJ(3));
t28 = t13 * t19;
t10 = cos(t11);
t27 = t15 * t10;
t26 = t15 * t17;
t25 = t15 * t19;
t24 = t14 * t26;
t23 = g(1) * t15 + g(2) * t13;
t22 = -t14 * t29 - t25;
t1 = -t9 * t30 - t27;
t3 = t13 * t10 - t14 * t33;
t21 = -g(1) * t3 - g(2) * t1 + t9 * t34;
t2 = t10 * t30 - t33;
t4 = t13 * t9 + t14 * t27;
t20 = g(1) * t4 + g(2) * t2 + t10 * t34;
t8 = pkin(3) * t28;
t6 = -g(1) * t13 + g(2) * t15;
t5 = g(3) * t14 - t23 * t12;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t24 + t28) - g(2) * t22 + t17 * t34, -g(1) * (-t14 * t25 - t29) - g(2) * (-t14 * t28 + t26) + t19 * t34, 0, 0, 0, 0, 0, 0, 0, 0, t21, t20, 0, -g(1) * t8 + (g(2) * t25 + (t23 * t14 + t34) * t17) * pkin(3), 0, 0, 0, 0, 0, 0, t21 * t18, -t21 * t16, -t20, -g(1) * (-pkin(3) * t24 + t3 * pkin(4) + t4 * pkin(6) + t8) - g(2) * (t22 * pkin(3) + t1 * pkin(4) + t2 * pkin(6)) - (-pkin(3) * t17 - pkin(4) * t9 + pkin(6) * t10) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t15 * t31 - t4 * t16) - g(2) * (t13 * t31 - t2 * t16) - g(3) * (-t10 * t32 - t14 * t18), -g(1) * (-t15 * t32 - t4 * t18) - g(2) * (-t13 * t32 - t2 * t18) - g(3) * (-t10 * t31 + t14 * t16), 0, 0;];
taug_reg = t7;
