% Calculate inertial parameters regressor of gravitation load for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:57
% EndTime: 2019-12-31 20:49:58
% DurationCPUTime: 0.22s
% Computational Cost: add. (246->52), mult. (207->56), div. (0->0), fcn. (187->8), ass. (0->35)
t24 = qJ(3) + pkin(8);
t18 = sin(t24);
t19 = cos(t24);
t49 = t19 * pkin(4) + t18 * qJ(5);
t25 = qJ(1) + qJ(2);
t20 = sin(t25);
t21 = cos(t25);
t7 = g(1) * t20 - g(2) * t21;
t8 = g(1) * t21 + g(2) * t20;
t28 = sin(qJ(1));
t43 = t28 * pkin(1);
t26 = -qJ(4) - pkin(7);
t42 = t21 * t26;
t41 = t21 * pkin(2) + t20 * pkin(7);
t29 = cos(qJ(3));
t22 = t29 * pkin(3);
t17 = t22 + pkin(2);
t12 = t21 * t17;
t39 = t49 * t21 + t12;
t38 = -t20 * pkin(2) + t21 * pkin(7);
t37 = -t20 * t26 + t12;
t30 = cos(qJ(1));
t36 = g(1) * t28 - g(2) * t30;
t34 = -t20 * t17 - t42;
t27 = sin(qJ(3));
t32 = -g(3) * t29 + t8 * t27;
t31 = (-g(1) * (-t17 - t49) + g(2) * t26) * t20;
t23 = t30 * pkin(1);
t6 = t7 * t29;
t5 = t7 * t27;
t4 = t7 * t19;
t3 = t7 * t18;
t2 = g(3) * t18 + t8 * t19;
t1 = -g(3) * t19 + t8 * t18;
t9 = [0, 0, 0, 0, 0, 0, t36, g(1) * t30 + g(2) * t28, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t36 * pkin(1), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t38 - t43) - g(2) * (t23 + t41), 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * (t34 - t43) - g(2) * (t23 + t37), 0, 0, 0, 0, 0, 0, t4, -t8, t3, -g(1) * (-t42 - t43) - g(2) * (t23 + t39) + t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * t38 - g(2) * t41, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * t34 - g(2) * t37, 0, 0, 0, 0, 0, 0, t4, -t8, t3, g(1) * t42 - g(2) * t39 + t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, g(3) * t27 + t8 * t29, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t32 * pkin(3), 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(3) * (t22 + t49) + t8 * (pkin(3) * t27 + pkin(4) * t18 - qJ(5) * t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t9;
