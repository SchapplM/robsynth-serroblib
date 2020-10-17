% Calculate minimal parameter regressor of gravitation load for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:47:03
% EndTime: 2019-05-06 09:47:04
% DurationCPUTime: 0.25s
% Computational Cost: add. (256->58), mult. (248->91), div. (0->0), fcn. (262->12), ass. (0->48)
t30 = sin(qJ(1));
t32 = cos(qJ(1));
t12 = g(1) * t32 + g(2) * t30;
t25 = qJ(2) + pkin(10);
t18 = sin(t25);
t20 = cos(t25);
t35 = -g(3) * t20 + t12 * t18;
t51 = g(3) * t18;
t24 = pkin(11) + qJ(5);
t21 = qJ(6) + t24;
t14 = sin(t21);
t49 = t30 * t14;
t15 = cos(t21);
t48 = t30 * t15;
t17 = sin(t24);
t47 = t30 * t17;
t19 = cos(t24);
t46 = t30 * t19;
t26 = sin(pkin(11));
t45 = t30 * t26;
t27 = cos(pkin(11));
t44 = t30 * t27;
t43 = t32 * t14;
t42 = t32 * t15;
t41 = t32 * t17;
t40 = t32 * t19;
t39 = t32 * t26;
t38 = t32 * t27;
t11 = g(1) * t30 - g(2) * t32;
t37 = t20 * pkin(3) + t18 * qJ(4);
t29 = sin(qJ(2));
t31 = cos(qJ(2));
t33 = -g(3) * t31 + t12 * t29;
t28 = -qJ(3) - pkin(7);
t22 = t31 * pkin(2);
t16 = t22 + pkin(1);
t13 = t32 * t16;
t10 = t20 * t40 + t47;
t9 = -t20 * t41 + t46;
t8 = -t20 * t46 + t41;
t7 = t20 * t47 + t40;
t6 = t20 * t42 + t49;
t5 = -t20 * t43 + t48;
t4 = -t20 * t48 + t43;
t3 = t20 * t49 + t42;
t2 = g(1) * t6 - g(2) * t4 + t15 * t51;
t1 = -g(1) * t5 + g(2) * t3 + t14 * t51;
t23 = [0, t11, t12, 0, 0, 0, 0, 0, t11 * t31, -t11 * t29, -t12, -g(1) * (-t30 * t16 - t32 * t28) - g(2) * (-t30 * t28 + t13) -g(1) * (-t20 * t44 + t39) - g(2) * (t20 * t38 + t45) -g(1) * (t20 * t45 + t38) - g(2) * (-t20 * t39 + t44) t11 * t18, -g(2) * t13 + (g(1) * t28 - g(2) * t37) * t32 + (-g(1) * (-t16 - t37) + g(2) * t28) * t30, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, t33, g(3) * t29 + t12 * t31, 0, t33 * pkin(2), t35 * t27, -t35 * t26, -t12 * t20 - t51, -g(3) * (t22 + t37) + t12 * (pkin(2) * t29 + pkin(3) * t18 - qJ(4) * t20) 0, 0, 0, 0, 0, t35 * t19, -t35 * t17, 0, 0, 0, 0, 0, t35 * t15, -t35 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t17 * t51, g(1) * t10 - g(2) * t8 + t19 * t51, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t23;
