% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x33]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 11:05:18
% EndTime: 2019-05-07 11:05:19
% DurationCPUTime: 0.25s
% Computational Cost: add. (319->52), mult. (303->86), div. (0->0), fcn. (331->10), ass. (0->48)
t30 = sin(qJ(1));
t33 = cos(qJ(1));
t39 = g(1) * t33 + g(2) * t30;
t28 = sin(qJ(3));
t31 = cos(qJ(3));
t42 = t33 * t31;
t32 = cos(qJ(2));
t48 = t30 * t32;
t15 = t28 * t48 + t42;
t43 = t33 * t28;
t17 = t30 * t31 - t32 * t43;
t29 = sin(qJ(2));
t50 = g(3) * t29;
t55 = -g(1) * t17 + g(2) * t15 + t28 * t50;
t13 = -g(3) * t32 + t39 * t29;
t26 = qJ(3) + pkin(11) + qJ(5);
t24 = qJ(6) + t26;
t20 = sin(t24);
t47 = t33 * t20;
t21 = cos(t24);
t46 = t33 * t21;
t22 = sin(t26);
t45 = t33 * t22;
t23 = cos(t26);
t44 = t33 * t23;
t40 = pkin(3) * t28 + pkin(7);
t38 = g(1) * t30 - g(2) * t33;
t25 = t31 * pkin(3) + pkin(2);
t27 = -qJ(4) - pkin(8);
t37 = t32 * t25 - t29 * t27;
t35 = pkin(1) + t37;
t19 = t38 * t29;
t18 = t30 * t28 + t32 * t42;
t16 = -t31 * t48 + t43;
t14 = t39 * t32 + t50;
t12 = t30 * t22 + t32 * t44;
t11 = t30 * t23 - t32 * t45;
t10 = -t23 * t48 + t45;
t9 = t22 * t48 + t44;
t8 = t30 * t20 + t32 * t46;
t7 = t30 * t21 - t32 * t47;
t6 = -t21 * t48 + t47;
t5 = t20 * t48 + t46;
t4 = g(1) * t12 - g(2) * t10 + t23 * t50;
t3 = -g(1) * t11 + g(2) * t9 + t22 * t50;
t2 = g(1) * t8 - g(2) * t6 + t21 * t50;
t1 = -g(1) * t7 + g(2) * t5 + t20 * t50;
t34 = [0, t38, t39, 0, 0, 0, 0, 0, t38 * t32, -t19, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, -g(1) * t15 - g(2) * t17, t19 (-g(1) * t40 - g(2) * t35) * t33 + (g(1) * t35 - g(2) * t40) * t30, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, t13 * t31, -t13 * t28, -t14, -g(3) * t37 + t39 * (t25 * t29 + t27 * t32) 0, 0, 0, 0, 0, t13 * t23, -t13 * t22, 0, 0, 0, 0, 0, t13 * t21, -t13 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, g(1) * t18 - g(2) * t16 + t31 * t50, 0, t55 * pkin(3), 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t34;
