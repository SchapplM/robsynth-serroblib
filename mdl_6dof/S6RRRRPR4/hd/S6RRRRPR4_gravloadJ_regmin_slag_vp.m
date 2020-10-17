% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:11:16
% EndTime: 2019-05-07 20:11:17
% DurationCPUTime: 0.24s
% Computational Cost: add. (299->51), mult. (299->75), div. (0->0), fcn. (308->10), ass. (0->49)
t33 = cos(qJ(2));
t32 = cos(qJ(4));
t21 = t32 * pkin(4) + pkin(3);
t27 = qJ(2) + qJ(3);
t24 = sin(t27);
t25 = cos(t27);
t28 = -qJ(5) - pkin(9);
t42 = t25 * t21 - t24 * t28;
t61 = t33 * pkin(2) + t42;
t31 = sin(qJ(1));
t34 = cos(qJ(1));
t41 = g(1) * t34 + g(2) * t31;
t45 = t34 * t32;
t29 = sin(qJ(4));
t50 = t31 * t29;
t14 = t25 * t50 + t45;
t46 = t34 * t29;
t49 = t31 * t32;
t16 = -t25 * t46 + t49;
t55 = g(3) * t24;
t60 = -g(1) * t16 + g(2) * t14 + t29 * t55;
t7 = -g(3) * t25 + t41 * t24;
t23 = qJ(4) + pkin(11) + qJ(6);
t19 = sin(t23);
t52 = t31 * t19;
t20 = cos(t23);
t51 = t31 * t20;
t48 = t34 * t19;
t47 = t34 * t20;
t43 = pkin(4) * t29 + pkin(7) + pkin(8);
t40 = g(1) * t31 - g(2) * t34;
t39 = t21 * t24 + t25 * t28;
t38 = pkin(1) + t61;
t30 = sin(qJ(2));
t17 = t25 * t45 + t50;
t15 = -t25 * t49 + t46;
t13 = t40 * t24;
t12 = t25 * t47 + t52;
t11 = -t25 * t48 + t51;
t10 = -t25 * t51 + t48;
t9 = t25 * t52 + t47;
t8 = t41 * t25 + t55;
t6 = t7 * t32;
t5 = t7 * t29;
t4 = t7 * t20;
t3 = t7 * t19;
t2 = g(1) * t12 - g(2) * t10 + t20 * t55;
t1 = -g(1) * t11 + g(2) * t9 + t19 * t55;
t18 = [0, t40, t41, 0, 0, 0, 0, 0, t40 * t33, -t40 * t30, 0, 0, 0, 0, 0, t40 * t25, -t13, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t13 (-g(1) * t43 - g(2) * t38) * t34 + (g(1) * t38 - g(2) * t43) * t31, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t33 + t41 * t30, g(3) * t30 + t41 * t33, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * t61 + t41 * (pkin(2) * t30 + t39) 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * t42 + t41 * t39, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, g(1) * t17 - g(2) * t15 + t32 * t55, 0, t60 * pkin(4), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t18;
