% Calculate inertial parameters regressor of gravitation load for
% S4RPRR7
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:04
% EndTime: 2019-12-31 16:54:05
% DurationCPUTime: 0.14s
% Computational Cost: add. (103->40), mult. (144->58), div. (0->0), fcn. (143->8), ass. (0->27)
t19 = sin(qJ(1));
t21 = cos(qJ(1));
t9 = g(1) * t21 + g(2) * t19;
t14 = pkin(7) + qJ(3);
t11 = sin(t14);
t12 = cos(t14);
t22 = -g(3) * t12 + t9 * t11;
t30 = g(3) * t11;
t18 = sin(qJ(4));
t28 = t19 * t18;
t20 = cos(qJ(4));
t27 = t19 * t20;
t26 = t21 * t18;
t25 = t21 * t20;
t24 = t12 * pkin(3) + t11 * pkin(6);
t8 = g(1) * t19 - g(2) * t21;
t17 = -pkin(5) - qJ(2);
t16 = cos(pkin(7));
t10 = t16 * pkin(2) + pkin(1);
t7 = t21 * t10;
t6 = t12 * t25 + t28;
t5 = -t12 * t26 + t27;
t4 = -t12 * t27 + t26;
t3 = t12 * t28 + t25;
t2 = t8 * t11;
t1 = t9 * t12 + t30;
t13 = [0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, t8 * t16, -t8 * sin(pkin(7)), -t9, -g(1) * (-t19 * pkin(1) + t21 * qJ(2)) - g(2) * (t21 * pkin(1) + t19 * qJ(2)), 0, 0, 0, 0, 0, 0, t8 * t12, -t2, -t9, -g(1) * (-t19 * t10 - t21 * t17) - g(2) * (-t19 * t17 + t7), 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t2, -g(2) * t7 + (g(1) * t17 - g(2) * t24) * t21 + (-g(1) * (-t10 - t24) + g(2) * t17) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, t1, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t20, -t22 * t18, -t1, -g(3) * t24 + t9 * (pkin(3) * t11 - pkin(6) * t12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 + t18 * t30, g(1) * t6 - g(2) * t4 + t20 * t30, 0, 0;];
taug_reg = t13;
