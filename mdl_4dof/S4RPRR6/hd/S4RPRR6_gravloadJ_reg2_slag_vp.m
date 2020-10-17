% Calculate inertial parameters regressor of gravitation load for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:41
% EndTime: 2019-12-31 16:52:42
% DurationCPUTime: 0.10s
% Computational Cost: add. (96->30), mult. (102->39), div. (0->0), fcn. (93->8), ass. (0->19)
t17 = cos(pkin(7));
t8 = t17 * pkin(2) + pkin(1);
t18 = -pkin(5) - qJ(2);
t15 = pkin(7) + qJ(3);
t19 = sin(qJ(1));
t20 = cos(qJ(1));
t5 = g(1) * t20 + g(2) * t19;
t4 = g(1) * t19 - g(2) * t20;
t10 = cos(t15);
t9 = sin(t15);
t21 = -g(3) * t10 + t5 * t9;
t14 = -pkin(6) + t18;
t11 = qJ(4) + t15;
t7 = cos(t11);
t6 = sin(t11);
t3 = pkin(3) * t10 + t8;
t2 = g(3) * t6 + t5 * t7;
t1 = -g(3) * t7 + t5 * t6;
t12 = [0, 0, 0, 0, 0, 0, t4, t5, 0, 0, 0, 0, 0, 0, 0, 0, t4 * t17, -t4 * sin(pkin(7)), -t5, -g(1) * (-t19 * pkin(1) + t20 * qJ(2)) - g(2) * (t20 * pkin(1) + t19 * qJ(2)), 0, 0, 0, 0, 0, 0, t4 * t10, -t4 * t9, -t5, -g(1) * (-t20 * t18 - t19 * t8) - g(2) * (-t19 * t18 + t20 * t8), 0, 0, 0, 0, 0, 0, t4 * t7, -t4 * t6, -t5, -g(1) * (-t20 * t14 - t19 * t3) - g(2) * (-t19 * t14 + t20 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, g(3) * t9 + t5 * t10, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t21 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t12;
