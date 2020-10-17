% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 22:20:50
% EndTime: 2019-05-06 22:20:52
% DurationCPUTime: 0.50s
% Computational Cost: add. (378->114), mult. (834->146), div. (0->0), fcn. (940->10), ass. (0->73)
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t20 = g(1) * t46 + g(2) * t43;
t42 = sin(qJ(2));
t89 = t20 * t42;
t39 = qJ(5) + qJ(6);
t30 = sin(t39);
t41 = sin(qJ(4));
t45 = cos(qJ(2));
t76 = cos(qJ(4));
t68 = t42 * t76;
t19 = -t41 * t45 + t68;
t12 = t19 * t43;
t74 = t45 * t46;
t14 = t41 * t74 - t46 * t68;
t18 = t41 * t42 + t45 * t76;
t52 = g(1) * t14 - g(2) * t12 + g(3) * t18;
t88 = t52 * t30;
t31 = cos(t39);
t87 = t52 * t31;
t40 = sin(qJ(5));
t86 = t52 * t40;
t44 = cos(qJ(5));
t85 = t52 * t44;
t32 = t42 * qJ(3);
t72 = pkin(2) * t45 + t32;
t13 = t18 * t43;
t56 = t13 * t40 - t44 * t46;
t77 = g(3) * t19;
t15 = t18 * t46;
t8 = -t15 * t40 - t43 * t44;
t84 = -g(1) * t8 + g(2) * t56 + t40 * t77;
t81 = pkin(2) * t42;
t80 = g(1) * t43;
t34 = t45 * pkin(3);
t73 = t46 * t40;
t71 = pkin(1) * t46 + pkin(7) * t43;
t70 = qJ(3) * t45;
t69 = t34 + t72;
t36 = t46 * pkin(7);
t67 = -t46 * pkin(8) + t36;
t66 = pkin(2) * t74 + t32 * t46 + t71;
t65 = pkin(3) * t74 + t66;
t64 = t12 * pkin(4) + t13 * pkin(9);
t63 = -t14 * pkin(4) + t15 * pkin(9);
t62 = -pkin(4) * t18 + pkin(9) * t19;
t5 = -g(1) * t12 - g(2) * t14;
t61 = -g(2) * t46 + t80;
t29 = pkin(5) * t44 + pkin(4);
t47 = -pkin(10) - pkin(9);
t60 = t12 * t29 - t13 * t47;
t59 = t13 * t31 + t30 * t46;
t58 = t13 * t30 - t31 * t46;
t57 = t13 * t44 + t73;
t55 = -t14 * t29 - t15 * t47;
t54 = -t18 * t29 - t19 * t47;
t53 = -pkin(1) - t72;
t4 = g(1) * t15 + g(2) * t13 + t77;
t51 = g(1) * (t53 - t34);
t49 = (pkin(2) + pkin(3)) * t89;
t48 = (g(2) * pkin(8) - t51) * t43;
t25 = t46 * t70;
t23 = t43 * t70;
t17 = t61 * t45;
t16 = t61 * t42;
t11 = g(3) * t42 + t20 * t45;
t10 = -g(3) * t45 + t89;
t9 = t15 * t44 - t40 * t43;
t7 = t15 * t31 - t30 * t43;
t6 = -t15 * t30 - t31 * t43;
t2 = g(1) * t7 + g(2) * t59 + t31 * t77;
t1 = -g(1) * t6 + g(2) * t58 + t30 * t77;
t3 = [0, 0, 0, 0, 0, 0, t61, t20, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, -t20, -g(1) * (-t43 * pkin(1) + t36) - g(2) * t71, 0, 0, 0, 0, 0, 0, t17, -t20, t16, -g(1) * t36 - g(2) * t66 - t53 * t80, 0, 0, 0, 0, 0, 0, g(1) * t13 - g(2) * t15, -t5, t20, -g(1) * t67 - g(2) * t65 + t48, 0, 0, 0, 0, 0, 0, g(1) * t57 - g(2) * t9, -g(1) * t56 - g(2) * t8, t5, -g(1) * (-t13 * pkin(4) + t12 * pkin(9) + t67) - g(2) * (pkin(4) * t15 + pkin(9) * t14 + t65) + t48, 0, 0, 0, 0, 0, 0, g(1) * t59 - g(2) * t7, -g(1) * t58 - g(2) * t6, t5, -g(1) * (-pkin(5) * t73 - t12 * t47 - t13 * t29 + t67) - g(2) * (-t14 * t47 + t15 * t29 + t65) + (-t51 - g(2) * (-pkin(5) * t40 - pkin(8))) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, -t11, -g(1) * (-t46 * t81 + t25) - g(2) * (-t43 * t81 + t23) - g(3) * t72, 0, 0, 0, 0, 0, 0, -t52, -t4, 0, -g(1) * t25 - g(2) * t23 - g(3) * t69 + t49, 0, 0, 0, 0, 0, 0, -t85, t86, t4, -g(1) * (t25 - t63) - g(2) * (t23 - t64) - g(3) * (-t62 + t69) + t49, 0, 0, 0, 0, 0, 0, -t87, t88, t4, -g(1) * (t25 - t55) - g(2) * (t23 - t60) - g(3) * (-t54 + t69) + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t4, 0, 0, 0, 0, 0, 0, 0, 0, t85, -t86, -t4, -g(1) * t63 - g(2) * t64 - g(3) * t62, 0, 0, 0, 0, 0, 0, t87, -t88, -t4, -g(1) * t55 - g(2) * t60 - g(3) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, g(1) * t9 + g(2) * t57 + t44 * t77, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t84 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t3;
