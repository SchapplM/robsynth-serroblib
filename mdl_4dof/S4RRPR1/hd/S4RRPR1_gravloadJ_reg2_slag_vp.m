% Calculate inertial parameters regressor of gravitation load for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:21:41
% EndTime: 2019-05-04 19:21:41
% DurationCPUTime: 0.09s
% Computational Cost: add. (107->22), mult. (60->30), div. (0->0), fcn. (50->8), ass. (0->25)
t18 = qJ(1) + qJ(2);
t15 = sin(t18);
t25 = pkin(2) * t15;
t19 = sin(qJ(1));
t24 = t19 * pkin(1);
t16 = cos(t18);
t12 = pkin(2) * t16;
t20 = cos(qJ(1));
t23 = t20 * pkin(1) + t12;
t14 = pkin(7) + t18;
t10 = sin(t14);
t22 = -pkin(3) * t10 - t25;
t5 = g(1) * t15 - g(2) * t16;
t21 = g(1) * t19 - g(2) * t20;
t13 = qJ(4) + t14;
t11 = cos(t14);
t9 = cos(t13);
t8 = sin(t13);
t7 = pkin(3) * t11;
t6 = g(1) * t16 + g(2) * t15;
t4 = g(1) * t11 + g(2) * t10;
t3 = g(1) * t10 - g(2) * t11;
t2 = g(1) * t9 + g(2) * t8;
t1 = g(1) * t8 - g(2) * t9;
t17 = [0, 0, 0, 0, 0, 0, t21, g(1) * t20 + g(2) * t19, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t21 * pkin(1), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * (-t24 - t25) - g(2) * t23, 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t22 - t24) - g(2) * (t7 + t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(2), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * t22 - g(2) * (t7 + t12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t17;
