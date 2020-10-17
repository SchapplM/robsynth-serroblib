% Calculate minimal parameter regressor of gravitation load for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x31]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:06:01
% EndTime: 2019-05-05 16:06:01
% DurationCPUTime: 0.14s
% Computational Cost: add. (142->40), mult. (152->48), div. (0->0), fcn. (156->10), ass. (0->28)
t20 = pkin(10) + qJ(4);
t15 = qJ(5) + t20;
t11 = sin(t15);
t12 = cos(t15);
t24 = sin(qJ(1));
t26 = cos(qJ(1));
t9 = g(1) * t24 - g(2) * t26;
t34 = -g(3) * t11 + t9 * t12;
t32 = g(3) * t12;
t23 = sin(qJ(6));
t31 = t24 * t23;
t25 = cos(qJ(6));
t30 = t24 * t25;
t29 = t26 * t23;
t28 = t26 * t25;
t27 = t26 * pkin(1) + t24 * qJ(2);
t10 = g(1) * t26 + g(2) * t24;
t17 = t26 * qJ(2);
t14 = cos(t20);
t13 = sin(t20);
t8 = t11 * t28 - t31;
t7 = t11 * t29 + t30;
t6 = t11 * t30 + t29;
t5 = -t11 * t31 + t28;
t3 = t9 * t11 + t32;
t2 = t34 * t25;
t1 = t34 * t23;
t4 = [0, t9, t10, -t9, -t10, -g(1) * (-t24 * pkin(1) + t17) - g(2) * t27, -t10 * sin(pkin(10)) -t10 * cos(pkin(10)) t9, -g(1) * (t17 + (-pkin(1) - qJ(3)) * t24) - g(2) * (t26 * qJ(3) + t27) 0, 0, 0, 0, 0, -t10 * t13, -t10 * t14, 0, 0, 0, 0, 0, -t10 * t11, -t10 * t12, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t6, g(1) * t7 - g(2) * t5; 0, 0, 0, 0, 0, -t9, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t13 - t9 * t14, g(3) * t14 + t9 * t13, 0, 0, 0, 0, 0, -t34, t3, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t3, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7 + t23 * t32, g(1) * t6 - g(2) * t8 + t25 * t32;];
taug_reg  = t4;
