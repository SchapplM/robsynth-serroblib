% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:32:01
% EndTime: 2019-05-06 19:32:01
% DurationCPUTime: 0.14s
% Computational Cost: add. (251->33), mult. (192->50), div. (0->0), fcn. (195->10), ass. (0->33)
t23 = sin(qJ(1));
t26 = cos(qJ(1));
t12 = g(1) * t26 + g(2) * t23;
t19 = qJ(2) + pkin(11) + qJ(4);
t17 = qJ(5) + t19;
t13 = sin(t17);
t14 = cos(t17);
t3 = -g(3) * t14 + t12 * t13;
t33 = g(3) * t13;
t21 = sin(qJ(6));
t31 = t23 * t21;
t24 = cos(qJ(6));
t30 = t23 * t24;
t29 = t26 * t21;
t28 = t26 * t24;
t11 = g(1) * t23 - g(2) * t26;
t22 = sin(qJ(2));
t25 = cos(qJ(2));
t27 = -g(3) * t25 + t12 * t22;
t20 = -qJ(3) - pkin(7);
t18 = t25 * pkin(2) + pkin(1);
t16 = cos(t19);
t15 = sin(t19);
t10 = t14 * t28 + t31;
t9 = -t14 * t29 + t30;
t8 = -t14 * t30 + t29;
t7 = t14 * t31 + t28;
t6 = g(3) * t15 + t12 * t16;
t5 = -g(3) * t16 + t12 * t15;
t4 = t12 * t14 + t33;
t2 = t3 * t24;
t1 = t3 * t21;
t32 = [0, t11, t12, 0, 0, 0, 0, 0, t11 * t25, -t11 * t22, -t12, -g(1) * (-t23 * t18 - t26 * t20) - g(2) * (t26 * t18 - t23 * t20) 0, 0, 0, 0, 0, t11 * t16, -t11 * t15, 0, 0, 0, 0, 0, t11 * t14, -t11 * t13, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9; 0, 0, 0, 0, 0, 0, 0, 0, t27, g(3) * t22 + t12 * t25, 0, t27 * pkin(2), 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t21 * t33, g(1) * t10 - g(2) * t8 + t24 * t33;];
taug_reg  = t32;
