% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRR8
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
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 22:40:41
% EndTime: 2019-05-06 22:40:42
% DurationCPUTime: 0.26s
% Computational Cost: add. (337->57), mult. (311->97), div. (0->0), fcn. (346->12), ass. (0->49)
t32 = sin(qJ(1));
t34 = cos(qJ(1));
t40 = g(1) * t34 + g(2) * t32;
t31 = sin(qJ(2));
t33 = cos(qJ(2));
t17 = -g(3) * t33 + t40 * t31;
t51 = g(3) * t31;
t49 = t32 * t33;
t28 = pkin(11) + qJ(4);
t27 = qJ(5) + t28;
t24 = qJ(6) + t27;
t20 = sin(t24);
t48 = t34 * t20;
t21 = cos(t24);
t47 = t34 * t21;
t22 = sin(t27);
t46 = t34 * t22;
t23 = cos(t27);
t45 = t34 * t23;
t25 = sin(t28);
t44 = t34 * t25;
t26 = cos(t28);
t43 = t34 * t26;
t29 = sin(pkin(11));
t42 = t34 * t29;
t30 = cos(pkin(11));
t41 = t34 * t30;
t39 = g(1) * t32 - g(2) * t34;
t38 = t33 * pkin(2) + t31 * qJ(3);
t36 = pkin(1) + t38;
t19 = t39 * t31;
t18 = t40 * t33 + t51;
t16 = t32 * t25 + t33 * t43;
t15 = t32 * t26 - t33 * t44;
t14 = -t26 * t49 + t44;
t13 = t25 * t49 + t43;
t12 = t32 * t22 + t33 * t45;
t11 = t32 * t23 - t33 * t46;
t10 = -t23 * t49 + t46;
t9 = t22 * t49 + t45;
t8 = t32 * t20 + t33 * t47;
t7 = t32 * t21 - t33 * t48;
t6 = -t21 * t49 + t48;
t5 = t20 * t49 + t47;
t4 = g(1) * t12 - g(2) * t10 + t23 * t51;
t3 = -g(1) * t11 + g(2) * t9 + t22 * t51;
t2 = g(1) * t8 - g(2) * t6 + t21 * t51;
t1 = -g(1) * t7 + g(2) * t5 + t20 * t51;
t35 = [0, t39, t40, 0, 0, 0, 0, 0, t39 * t33, -t19, -g(1) * (-t30 * t49 + t42) - g(2) * (t32 * t29 + t33 * t41) -g(1) * (t29 * t49 + t41) - g(2) * (t32 * t30 - t33 * t42) t19 (-g(1) * pkin(7) - g(2) * t36) * t34 + (-g(2) * pkin(7) + g(1) * t36) * t32, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, t17, t18, t17 * t30, -t17 * t29, -t18, -g(3) * t38 + t40 * (pkin(2) * t31 - qJ(3) * t33) 0, 0, 0, 0, 0, t17 * t26, -t17 * t25, 0, 0, 0, 0, 0, t17 * t23, -t17 * t22, 0, 0, 0, 0, 0, t17 * t21, -t17 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 + g(2) * t13 + t25 * t51, g(1) * t16 - g(2) * t14 + t26 * t51, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t35;
