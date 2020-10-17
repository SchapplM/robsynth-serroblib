% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:14:52
% EndTime: 2019-05-05 18:14:52
% DurationCPUTime: 0.13s
% Computational Cost: add. (186->35), mult. (144->51), div. (0->0), fcn. (145->10), ass. (0->32)
t14 = qJ(3) + pkin(11) + qJ(5);
t10 = cos(t14);
t15 = qJ(1) + pkin(10);
t12 = sin(t15);
t13 = cos(t15);
t26 = g(1) * t13 + g(2) * t12;
t9 = sin(t14);
t3 = -g(3) * t10 + t26 * t9;
t32 = g(3) * t9;
t17 = sin(qJ(6));
t30 = t12 * t17;
t20 = cos(qJ(6));
t29 = t12 * t20;
t28 = t13 * t17;
t27 = t13 * t20;
t25 = g(1) * t12 - g(2) * t13;
t19 = sin(qJ(1));
t22 = cos(qJ(1));
t24 = g(1) * t19 - g(2) * t22;
t18 = sin(qJ(3));
t21 = cos(qJ(3));
t23 = -g(3) * t21 + t26 * t18;
t16 = -qJ(4) - pkin(7);
t11 = t21 * pkin(3) + pkin(2);
t8 = t10 * t27 + t30;
t7 = -t10 * t28 + t29;
t6 = -t10 * t29 + t28;
t5 = t10 * t30 + t27;
t4 = t26 * t10 + t32;
t2 = t3 * t20;
t1 = t3 * t17;
t31 = [0, t24, g(1) * t22 + g(2) * t19, t24 * pkin(1), 0, 0, 0, 0, 0, t25 * t21, -t25 * t18, -t26, -g(1) * (-t19 * pkin(1) - t12 * t11 - t13 * t16) - g(2) * (t22 * pkin(1) + t13 * t11 - t12 * t16) 0, 0, 0, 0, 0, t25 * t10, -t25 * t9, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, g(3) * t18 + t26 * t21, 0, t23 * pkin(3), 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t5 + t17 * t32, g(1) * t8 - g(2) * t6 + t20 * t32;];
taug_reg  = t31;
