% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:22:47
% EndTime: 2019-05-06 18:22:48
% DurationCPUTime: 0.33s
% Computational Cost: add. (437->78), mult. (448->115), div. (0->0), fcn. (482->10), ass. (0->53)
t35 = pkin(10) + qJ(4);
t28 = cos(t35);
t37 = cos(pkin(10));
t20 = t37 * pkin(3) + pkin(4) * t28 + pkin(2);
t34 = -pkin(9) - pkin(8) - qJ(3);
t38 = sin(qJ(2));
t40 = cos(qJ(2));
t71 = t40 * t20 - t38 * t34;
t39 = sin(qJ(1));
t41 = cos(qJ(1));
t48 = g(1) * t41 + g(2) * t39;
t27 = sin(t35);
t54 = t41 * t28;
t59 = t39 * t40;
t13 = t27 * t59 + t54;
t55 = t41 * t27;
t15 = t39 * t28 - t40 * t55;
t64 = g(3) * t38;
t70 = -g(1) * t15 + g(2) * t13 + t27 * t64;
t17 = -g(3) * t40 + t48 * t38;
t68 = g(1) * t39;
t29 = qJ(5) + t35;
t25 = sin(t29);
t62 = t25 * t38;
t26 = cos(t29);
t61 = t26 * t38;
t57 = t41 * t25;
t56 = t41 * t26;
t36 = sin(pkin(10));
t53 = t41 * t36;
t52 = t41 * t37;
t51 = t41 * pkin(1) + t39 * pkin(7);
t11 = -t39 * t26 + t40 * t57;
t9 = t25 * t59 + t56;
t49 = g(1) * t9 - g(2) * t11;
t47 = -g(2) * t41 + t68;
t46 = t40 * pkin(2) + t38 * qJ(3);
t44 = pkin(5) * t26 + qJ(6) * t25 + t20;
t1 = g(1) * t11 + g(2) * t9 + g(3) * t62;
t10 = t26 * t59 - t57;
t12 = t39 * t25 + t40 * t56;
t3 = g(1) * t12 + g(2) * t10 + g(3) * t61;
t42 = -g(1) * (-t11 * pkin(5) + t12 * qJ(6)) - g(2) * (-t9 * pkin(5) + t10 * qJ(6)) - g(3) * (-pkin(5) * t62 + qJ(6) * t61);
t31 = t41 * pkin(7);
t21 = t36 * pkin(3) + pkin(4) * t27;
t19 = t47 * t38;
t18 = t48 * t40 + t64;
t16 = t39 * t27 + t40 * t54;
t14 = -t28 * t59 + t55;
t6 = t17 * t26;
t5 = t17 * t25;
t4 = g(1) * t10 - g(2) * t12;
t2 = [0, t47, t48, 0, 0, 0, 0, 0, t47 * t40, -t19, -g(1) * (-t37 * t59 + t53) - g(2) * (t39 * t36 + t40 * t52) -g(1) * (t36 * t59 + t52) - g(2) * (t39 * t37 - t40 * t53) t19, -g(1) * t31 - g(2) * (t46 * t41 + t51) - (-pkin(1) - t46) * t68, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, 0, 0, 0, 0, 0, t4, -t49, t4, t19, t49, -g(1) * (-t10 * pkin(5) - t9 * qJ(6) + t41 * t21 + t31) - g(2) * (t12 * pkin(5) + t11 * qJ(6) + t71 * t41 + t51) + (-g(1) * (-pkin(1) - t71) - g(2) * t21) * t39; 0, 0, 0, 0, 0, 0, 0, 0, t17, t18, t17 * t37, -t17 * t36, -t18, -g(3) * t46 + t48 * (pkin(2) * t38 - qJ(3) * t40) 0, 0, 0, 0, 0, t17 * t28, -t17 * t27, 0, 0, 0, 0, 0, t6, -t5, t6, -t18, t5 (-g(3) * t44 + t48 * t34) * t40 + (g(3) * t34 + t48 * t44) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, g(1) * t16 - g(2) * t14 + t28 * t64, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t70 * pkin(4) + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
