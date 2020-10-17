% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPR6
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
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:51:31
% EndTime: 2019-05-07 20:51:32
% DurationCPUTime: 0.26s
% Computational Cost: add. (323->58), mult. (330->96), div. (0->0), fcn. (355->10), ass. (0->49)
t34 = sin(qJ(1));
t37 = cos(qJ(1));
t43 = g(1) * t37 + g(2) * t34;
t31 = qJ(3) + qJ(4);
t28 = cos(t31);
t36 = cos(qJ(2));
t27 = sin(t31);
t48 = t37 * t27;
t11 = t34 * t28 - t36 * t48;
t33 = sin(qJ(2));
t54 = g(3) * t33;
t47 = t37 * t28;
t51 = t34 * t36;
t9 = t27 * t51 + t47;
t3 = -g(1) * t11 + g(2) * t9 + t27 * t54;
t13 = -g(3) * t36 + t43 * t33;
t32 = sin(qJ(3));
t21 = t32 * pkin(3) + pkin(4) * t27;
t52 = pkin(7) + t21;
t26 = pkin(11) + qJ(6) + t31;
t23 = sin(t26);
t50 = t37 * t23;
t24 = cos(t26);
t49 = t37 * t24;
t46 = t37 * t32;
t35 = cos(qJ(3));
t45 = t37 * t35;
t22 = t35 * pkin(3) + pkin(4) * t28;
t42 = g(1) * t34 - g(2) * t37;
t20 = pkin(2) + t22;
t30 = -qJ(5) - pkin(9) - pkin(8);
t41 = t36 * t20 - t33 * t30;
t39 = pkin(1) + t41;
t19 = t42 * t33;
t18 = t34 * t32 + t36 * t45;
t17 = t34 * t35 - t36 * t46;
t16 = -t35 * t51 + t46;
t15 = t32 * t51 + t45;
t14 = t43 * t36 + t54;
t12 = t34 * t27 + t36 * t47;
t10 = -t28 * t51 + t48;
t8 = t34 * t23 + t36 * t49;
t7 = t34 * t24 - t36 * t50;
t6 = -t24 * t51 + t50;
t5 = t23 * t51 + t49;
t4 = g(1) * t12 - g(2) * t10 + t28 * t54;
t2 = g(1) * t8 - g(2) * t6 + t24 * t54;
t1 = -g(1) * t7 + g(2) * t5 + t23 * t54;
t25 = [0, t42, t43, 0, 0, 0, 0, 0, t42 * t36, -t19, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, -g(1) * t15 - g(2) * t17, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t19 (-g(1) * t52 - g(2) * t39) * t37 + (g(1) * t39 - g(2) * t52) * t34, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, t13 * t35, -t13 * t32, 0, 0, 0, 0, 0, t13 * t28, -t13 * t27, -t14, -g(3) * t41 + t43 * (t20 * t33 + t30 * t36) 0, 0, 0, 0, 0, t13 * t24, -t13 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t17 + g(2) * t15 + t32 * t54, g(1) * t18 - g(2) * t16 + t35 * t54, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * (-t37 * t36 * t21 + t34 * t22) - g(2) * (-t21 * t51 - t37 * t22) + t21 * t54, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t3 * pkin(4), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t25;
