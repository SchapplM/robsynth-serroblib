% Calculate inertial parameters regressor of gravitation load for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:05
% EndTime: 2019-12-31 17:17:06
% DurationCPUTime: 0.13s
% Computational Cost: add. (120->39), mult. (158->54), div. (0->0), fcn. (143->6), ass. (0->25)
t15 = qJ(2) + qJ(3);
t12 = sin(t15);
t13 = cos(t15);
t28 = t13 * pkin(3) + t12 * qJ(4);
t27 = pkin(3) * t12;
t25 = qJ(4) * t13;
t16 = sin(qJ(2));
t24 = -pkin(2) * t16 - t27;
t17 = sin(qJ(1));
t19 = cos(qJ(1));
t5 = g(1) * t19 + g(2) * t17;
t23 = g(1) * t17 - g(2) * t19;
t18 = cos(qJ(2));
t21 = -g(3) * t18 + t5 * t16;
t20 = -pkin(6) - pkin(5);
t14 = t18 * pkin(2);
t11 = t14 + pkin(1);
t8 = t19 * t11;
t7 = t19 * t25;
t6 = t17 * t25;
t4 = t23 * t13;
t3 = t23 * t12;
t2 = g(3) * t12 + t5 * t13;
t1 = -g(3) * t13 + t5 * t12;
t9 = [0, 0, 0, 0, 0, 0, t23, t5, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t18, -t23 * t16, -t5, -g(1) * (-t17 * pkin(1) + t19 * pkin(5)) - g(2) * (t19 * pkin(1) + t17 * pkin(5)), 0, 0, 0, 0, 0, 0, t4, -t3, -t5, -g(1) * (-t17 * t11 - t19 * t20) - g(2) * (-t17 * t20 + t8), 0, 0, 0, 0, 0, 0, t4, -t5, t3, -g(2) * t8 + (g(1) * t20 - g(2) * t28) * t19 + (-g(1) * (-t11 - t28) + g(2) * t20) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, g(3) * t16 + t5 * t18, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t21 * pkin(2), 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (t24 * t19 + t7) - g(2) * (t24 * t17 + t6) - g(3) * (t14 + t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (-t19 * t27 + t7) - g(2) * (-t17 * t27 + t6) - g(3) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t9;
