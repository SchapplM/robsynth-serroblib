% Calculate minimal parameter regressor of gravitation load for
% S6RPPRRR2
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
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:21:25
% EndTime: 2019-05-05 15:21:25
% DurationCPUTime: 0.18s
% Computational Cost: add. (214->42), mult. (170->65), div. (0->0), fcn. (186->12), ass. (0->37)
t17 = pkin(11) + qJ(4);
t11 = sin(t17);
t13 = cos(t17);
t18 = qJ(1) + pkin(10);
t12 = sin(t18);
t14 = cos(t18);
t30 = g(1) * t14 + g(2) * t12;
t27 = -g(3) * t13 + t30 * t11;
t40 = g(3) * t11;
t19 = qJ(5) + qJ(6);
t15 = sin(t19);
t38 = t12 * t15;
t16 = cos(t19);
t37 = t12 * t16;
t22 = sin(qJ(5));
t36 = t12 * t22;
t24 = cos(qJ(5));
t35 = t12 * t24;
t34 = t14 * t15;
t33 = t14 * t16;
t32 = t14 * t22;
t31 = t14 * t24;
t29 = g(1) * t12 - g(2) * t14;
t23 = sin(qJ(1));
t25 = cos(qJ(1));
t28 = g(1) * t23 - g(2) * t25;
t10 = t13 * t31 + t36;
t9 = -t13 * t32 + t35;
t8 = -t13 * t35 + t32;
t7 = t13 * t36 + t31;
t6 = t13 * t33 + t38;
t5 = -t13 * t34 + t37;
t4 = -t13 * t37 + t34;
t3 = t13 * t38 + t33;
t2 = g(1) * t6 - g(2) * t4 + t16 * t40;
t1 = -g(1) * t5 + g(2) * t3 + t15 * t40;
t20 = [0, t28, g(1) * t25 + g(2) * t23, t28 * pkin(1), t29 * cos(pkin(11)) -t29 * sin(pkin(11)) -t30, -g(1) * (-t23 * pkin(1) - t12 * pkin(2) + t14 * qJ(3)) - g(2) * (t25 * pkin(1) + t14 * pkin(2) + t12 * qJ(3)) 0, 0, 0, 0, 0, t29 * t13, -t29 * t11, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t30 * t13 + t40, 0, 0, 0, 0, 0, t27 * t24, -t27 * t22, 0, 0, 0, 0, 0, t27 * t16, -t27 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t22 * t40, g(1) * t10 - g(2) * t8 + t24 * t40, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t20;
