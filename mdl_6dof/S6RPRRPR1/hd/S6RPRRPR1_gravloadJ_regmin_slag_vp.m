% Calculate minimal parameter regressor of gravitation load for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:56:39
% EndTime: 2019-05-05 21:56:40
% DurationCPUTime: 0.15s
% Computational Cost: add. (189->40), mult. (157->59), div. (0->0), fcn. (155->12), ass. (0->36)
t22 = qJ(3) + qJ(4);
t16 = pkin(11) + t22;
t11 = sin(t16);
t12 = cos(t16);
t21 = qJ(1) + pkin(10);
t14 = sin(t21);
t15 = cos(t21);
t31 = g(1) * t15 + g(2) * t14;
t39 = -g(3) * t12 + t31 * t11;
t38 = g(3) * t11;
t23 = sin(qJ(6));
t36 = t14 * t23;
t26 = cos(qJ(6));
t35 = t14 * t26;
t34 = t15 * t23;
t33 = t15 * t26;
t18 = cos(t22);
t27 = cos(qJ(3));
t32 = t27 * pkin(3) + pkin(4) * t18;
t30 = g(1) * t14 - g(2) * t15;
t25 = sin(qJ(1));
t28 = cos(qJ(1));
t29 = g(1) * t25 - g(2) * t28;
t17 = sin(t22);
t3 = -g(3) * t18 + t31 * t17;
t24 = sin(qJ(3));
t20 = -qJ(5) - pkin(8) - pkin(7);
t9 = pkin(2) + t32;
t8 = t12 * t33 + t36;
t7 = -t12 * t34 + t35;
t6 = -t12 * t35 + t34;
t5 = t12 * t36 + t33;
t4 = g(3) * t17 + t31 * t18;
t2 = t39 * t26;
t1 = t39 * t23;
t10 = [0, t29, g(1) * t28 + g(2) * t25, t29 * pkin(1), 0, 0, 0, 0, 0, t30 * t27, -t30 * t24, 0, 0, 0, 0, 0, t30 * t18, -t30 * t17, -t31, -g(1) * (-t25 * pkin(1) - t14 * t9 - t15 * t20) - g(2) * (t28 * pkin(1) - t14 * t20 + t15 * t9) 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t27 + t31 * t24, g(3) * t24 + t31 * t27, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t32 - t31 * (-t24 * pkin(3) - pkin(4) * t17) 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t3 * pkin(4), 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t5 + t23 * t38, g(1) * t8 - g(2) * t6 + t26 * t38;];
taug_reg  = t10;
