% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:53:34
% EndTime: 2019-05-06 01:53:34
% DurationCPUTime: 0.40s
% Computational Cost: add. (285->83), mult. (442->117), div. (0->0), fcn. (446->8), ass. (0->60)
t38 = sin(qJ(1));
t41 = cos(qJ(1));
t79 = -g(1) * t38 + g(2) * t41;
t37 = sin(qJ(3));
t39 = cos(qJ(4));
t55 = t41 * t39;
t36 = sin(qJ(4));
t63 = t38 * t36;
t13 = -t37 * t63 + t55;
t56 = t41 * t36;
t62 = t38 * t39;
t15 = t37 * t56 + t62;
t40 = cos(qJ(3));
t69 = g(3) * t40;
t78 = -g(1) * t13 - g(2) * t15 + t36 * t69;
t35 = qJ(4) + qJ(5);
t26 = sin(t35);
t27 = cos(t35);
t57 = t41 * t27;
t65 = t38 * t26;
t7 = -t37 * t65 + t57;
t58 = t41 * t26;
t64 = t38 * t27;
t9 = t37 * t58 + t64;
t1 = -g(1) * t7 - g(2) * t9 + t26 * t69;
t12 = -g(3) * t37 - t40 * t79;
t75 = -pkin(1) - pkin(7);
t42 = -pkin(9) - pkin(8);
t68 = t36 * pkin(4);
t67 = t37 * t38;
t66 = t37 * t41;
t61 = t40 * t41;
t60 = t40 * t42;
t19 = pkin(5) * t26 + t68;
t59 = t41 * t19;
t30 = t39 * pkin(4);
t20 = pkin(5) * t27 + t30;
t54 = t41 * pkin(1) + t38 * qJ(2);
t51 = t41 * pkin(7) + t54;
t50 = g(2) * t51;
t48 = t37 * pkin(3) - t40 * pkin(8);
t22 = g(1) * t41 + g(2) * t38;
t18 = pkin(3) + t20;
t34 = -qJ(6) + t42;
t46 = t37 * t18 + t40 * t34;
t25 = t30 + pkin(3);
t44 = t37 * t25 + t60;
t29 = t41 * qJ(2);
t17 = t22 * t40;
t16 = t37 * t55 - t63;
t14 = t37 * t62 + t56;
t11 = g(1) * t67 - g(2) * t66 + t69;
t10 = t37 * t57 - t65;
t8 = t37 * t64 + t58;
t6 = t12 * t27;
t5 = t12 * t26;
t4 = -g(1) * t10 - g(2) * t8;
t3 = g(1) * t9 - g(2) * t7;
t2 = g(1) * t8 - g(2) * t10 + t27 * t69;
t21 = [0, 0, 0, 0, 0, 0, -t79, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, -t22, -g(1) * (-t38 * pkin(1) + t29) - g(2) * t54, 0, 0, 0, 0, 0, 0, -t22 * t37, -t17, -t79, -g(1) * (t75 * t38 + t29) - t50, 0, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t14, g(1) * t15 - g(2) * t13, t17, -g(1) * (pkin(3) * t66 - pkin(8) * t61 + t29) - t50 + (-g(1) * t75 - g(2) * t48) * t38, 0, 0, 0, 0, 0, 0, t4, t3, t17, -g(1) * (t25 * t66 + t41 * t60 + t29) - g(2) * (pkin(4) * t56 + t51) + (-g(1) * (-t68 + t75) - g(2) * t44) * t38, 0, 0, 0, 0, 0, 0, t4, t3, t17, -g(1) * (t18 * t66 + t34 * t61 + t29) - g(2) * (t51 + t59) + (-g(1) * (-t19 + t75) - g(2) * t46) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, 0, 0, 0, 0, 0, 0, 0, 0, -t12 * t39, t12 * t36, -t11, g(3) * t48 + t79 * (pkin(3) * t40 + pkin(8) * t37) 0, 0, 0, 0, 0, 0, -t6, t5, -t11, g(3) * t44 + t79 * (t25 * t40 - t37 * t42) 0, 0, 0, 0, 0, 0, -t6, t5, -t11, g(3) * t46 + t79 * (t18 * t40 - t34 * t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, g(1) * t14 - g(2) * t16 + t39 * t69, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t78 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t19 * t67 + t41 * t20) - g(2) * (t38 * t20 + t37 * t59) + t19 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12;];
taug_reg  = t21;
