% Calculate minimal parameter regressor of gravitation load for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:10:46
% EndTime: 2019-05-05 15:10:47
% DurationCPUTime: 0.12s
% Computational Cost: add. (194->35), mult. (142->51), div. (0->0), fcn. (146->12), ass. (0->30)
t16 = pkin(11) + qJ(4);
t15 = qJ(5) + t16;
t10 = cos(t15);
t17 = qJ(1) + pkin(10);
t12 = sin(t17);
t14 = cos(t17);
t26 = g(1) * t14 + g(2) * t12;
t9 = sin(t15);
t3 = -g(3) * t10 + t26 * t9;
t32 = g(3) * t9;
t20 = sin(qJ(6));
t30 = t12 * t20;
t22 = cos(qJ(6));
t29 = t12 * t22;
t28 = t14 * t20;
t27 = t14 * t22;
t25 = g(1) * t12 - g(2) * t14;
t21 = sin(qJ(1));
t23 = cos(qJ(1));
t24 = g(1) * t21 - g(2) * t23;
t13 = cos(t16);
t11 = sin(t16);
t8 = t10 * t27 + t30;
t7 = -t10 * t28 + t29;
t6 = -t10 * t29 + t28;
t5 = t10 * t30 + t27;
t4 = t26 * t10 + t32;
t2 = t3 * t22;
t1 = t3 * t20;
t18 = [0, t24, g(1) * t23 + g(2) * t21, t24 * pkin(1), t25 * cos(pkin(11)) -t25 * sin(pkin(11)) -t26, -g(1) * (-t21 * pkin(1) - t12 * pkin(2) + t14 * qJ(3)) - g(2) * (t23 * pkin(1) + t14 * pkin(2) + t12 * qJ(3)) 0, 0, 0, 0, 0, t25 * t13, -t25 * t11, 0, 0, 0, 0, 0, t25 * t10, -t25 * t9, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t13 + t26 * t11, g(3) * t11 + t26 * t13, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t5 + t20 * t32, g(1) * t8 - g(2) * t6 + t22 * t32;];
taug_reg  = t18;
