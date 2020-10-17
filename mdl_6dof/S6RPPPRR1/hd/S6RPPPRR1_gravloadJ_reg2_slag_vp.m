% Calculate inertial parameters regressor of gravitation load for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPPRR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:32:24
% EndTime: 2019-05-05 13:32:25
% DurationCPUTime: 0.23s
% Computational Cost: add. (182->57), mult. (178->65), div. (0->0), fcn. (173->8), ass. (0->30)
t19 = sin(qJ(5));
t22 = cos(qJ(5));
t27 = -t19 * pkin(5) + t22 * pkin(8);
t17 = qJ(1) + pkin(9);
t14 = sin(t17);
t15 = cos(t17);
t8 = g(1) * t15 + g(2) * t14;
t42 = -g(3) * t19 + t8 * t22;
t38 = g(3) * t22;
t37 = t15 * pkin(7);
t18 = sin(qJ(6));
t34 = t18 * t19;
t21 = cos(qJ(6));
t33 = t19 * t21;
t32 = -pkin(2) - qJ(4);
t23 = cos(qJ(1));
t31 = t23 * pkin(1) + t15 * pkin(2) + t14 * qJ(3);
t20 = sin(qJ(1));
t30 = -t20 * pkin(1) + t15 * qJ(3);
t29 = t15 * qJ(4) + t31;
t7 = g(1) * t14 - g(2) * t15;
t26 = g(1) * t20 - g(2) * t23;
t25 = t32 * t14 + t30;
t6 = t7 * t22;
t5 = -t14 * t18 + t15 * t33;
t4 = -t14 * t21 - t15 * t34;
t3 = -t14 * t33 - t15 * t18;
t2 = t14 * t34 - t15 * t21;
t1 = t8 * t19 + t38;
t9 = [0, 0, 0, 0, 0, 0, t26, g(1) * t23 + g(2) * t20, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t26 * pkin(1), 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -g(1) * (-t14 * pkin(2) + t30) - g(2) * t31, 0, 0, 0, 0, 0, 0, 0, -t8, t7, -g(1) * t25 - g(2) * t29, 0, 0, 0, 0, 0, 0, t7 * t19, t6, t8, -g(1) * (t25 - t37) - g(2) * (-t14 * pkin(7) + t29) 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5, -g(1) * t2 - g(2) * t4, -t6, -g(1) * (t30 - t37) - g(2) * (-t27 * t15 + t29) + (-g(1) * (t27 + t32) + g(2) * pkin(7)) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t1, 0, 0, 0, 0, 0, 0, 0, 0, -t42 * t21, t42 * t18, -t1, -g(3) * t27 - t8 * (pkin(5) * t22 + pkin(8) * t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 + g(2) * t2 + t18 * t38, g(1) * t5 - g(2) * t3 + t21 * t38, 0, 0;];
taug_reg  = t9;
