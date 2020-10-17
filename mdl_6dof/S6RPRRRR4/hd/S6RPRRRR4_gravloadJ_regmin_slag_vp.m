% Calculate minimal parameter regressor of gravitation load for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 03:14:53
% EndTime: 2019-05-06 03:14:53
% DurationCPUTime: 0.15s
% Computational Cost: add. (257->33), mult. (190->50), div. (0->0), fcn. (196->12), ass. (0->31)
t25 = sin(qJ(1));
t27 = cos(qJ(1));
t12 = g(1) * t27 + g(2) * t25;
t21 = pkin(11) + qJ(3);
t20 = qJ(4) + t21;
t17 = qJ(5) + t20;
t13 = sin(t17);
t14 = cos(t17);
t3 = -g(3) * t14 + t12 * t13;
t33 = g(3) * t13;
t24 = sin(qJ(6));
t31 = t25 * t24;
t26 = cos(qJ(6));
t30 = t25 * t26;
t29 = t27 * t24;
t28 = t27 * t26;
t11 = g(1) * t25 - g(2) * t27;
t19 = cos(t21);
t18 = sin(t21);
t16 = cos(t20);
t15 = sin(t20);
t10 = t14 * t28 + t31;
t9 = -t14 * t29 + t30;
t8 = -t14 * t30 + t29;
t7 = t14 * t31 + t28;
t6 = g(3) * t15 + t12 * t16;
t5 = -g(3) * t16 + t12 * t15;
t4 = t12 * t14 + t33;
t2 = t3 * t26;
t1 = t3 * t24;
t22 = [0, t11, t12, t11 * cos(pkin(11)) -t11 * sin(pkin(11)) -t12, -g(1) * (-t25 * pkin(1) + t27 * qJ(2)) - g(2) * (t27 * pkin(1) + t25 * qJ(2)) 0, 0, 0, 0, 0, t11 * t19, -t11 * t18, 0, 0, 0, 0, 0, t11 * t16, -t11 * t15, 0, 0, 0, 0, 0, t11 * t14, -t11 * t13, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9; 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t19 + t12 * t18, g(3) * t18 + t12 * t19, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t24 * t33, g(1) * t10 - g(2) * t8 + t26 * t33;];
taug_reg  = t22;
