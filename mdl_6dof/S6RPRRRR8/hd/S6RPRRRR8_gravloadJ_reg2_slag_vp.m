% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:15:55
% EndTime: 2019-05-06 04:15:56
% DurationCPUTime: 0.42s
% Computational Cost: add. (339->97), mult. (415->121), div. (0->0), fcn. (407->10), ass. (0->68)
t43 = sin(qJ(1));
t46 = cos(qJ(1));
t48 = -pkin(8) - pkin(7);
t42 = sin(qJ(3));
t75 = t42 * pkin(3);
t85 = t43 * t48 + t46 * t75;
t40 = qJ(3) + qJ(4);
t32 = sin(t40);
t34 = cos(t40);
t57 = -t32 * pkin(4) + t34 * pkin(9);
t44 = cos(qJ(5));
t63 = t46 * t44;
t41 = sin(qJ(5));
t68 = t43 * t41;
t14 = -t32 * t68 + t63;
t64 = t46 * t41;
t67 = t43 * t44;
t16 = t32 * t64 + t67;
t77 = g(3) * t34;
t84 = -g(1) * t14 - g(2) * t16 + t41 * t77;
t79 = g(2) * t46;
t20 = g(1) * t43 - t79;
t50 = -g(3) * t32 + t20 * t34;
t45 = cos(qJ(3));
t83 = pkin(3) * t45;
t82 = pkin(5) * t41;
t29 = t44 * pkin(5) + pkin(4);
t74 = t29 * t34;
t73 = t32 * t43;
t47 = -pkin(10) - pkin(9);
t72 = t32 * t47;
t71 = t34 * t43;
t39 = qJ(5) + qJ(6);
t31 = sin(t39);
t70 = t43 * t31;
t33 = cos(t39);
t69 = t43 * t33;
t66 = t46 * t31;
t65 = t46 * t33;
t62 = pkin(4) * t71 + pkin(9) * t73;
t61 = t46 * pkin(1) + t43 * qJ(2);
t59 = t43 * t75 + t61;
t36 = t46 * qJ(2);
t58 = -t43 * pkin(1) + t36;
t56 = t29 * t71 - t43 * t72;
t55 = -pkin(4) * t34 - pkin(9) * t32;
t21 = g(1) * t46 + g(2) * t43;
t53 = t32 * t29 + t34 * t47;
t52 = t58 + t85;
t51 = -t46 * t48 + t59;
t49 = g(3) * t42 - t20 * t45;
t26 = t43 * t83;
t19 = t46 * t72;
t17 = t32 * t63 - t68;
t15 = t32 * t67 + t64;
t13 = t21 * t34;
t12 = t32 * t65 - t70;
t11 = t32 * t66 + t69;
t10 = t32 * t69 + t66;
t9 = -t32 * t70 + t65;
t7 = g(1) * t73 - t32 * t79 + t77;
t6 = t50 * t44;
t5 = t50 * t41;
t4 = t50 * t33;
t3 = t50 * t31;
t2 = g(1) * t10 - g(2) * t12 + t33 * t77;
t1 = -g(1) * t9 - g(2) * t11 + t31 * t77;
t8 = [0, 0, 0, 0, 0, 0, t20, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, -t21, -g(1) * t58 - g(2) * t61, 0, 0, 0, 0, 0, 0, -t21 * t42, -t21 * t45, t20, -g(1) * (t36 + (-pkin(1) - pkin(7)) * t43) - g(2) * (t46 * pkin(7) + t61) 0, 0, 0, 0, 0, 0, -t21 * t32, -t13, t20, -g(1) * t52 - g(2) * t51, 0, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t15, g(1) * t16 - g(2) * t14, t13, -g(1) * (-t46 * t57 + t52) - g(2) * (-t43 * t57 + t51) 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, t13, -g(1) * (t36 + t85) - g(2) * t59 + (-g(1) * t53 - g(2) * (-t48 + t82)) * t46 + (-g(1) * (-pkin(1) - t82) - g(2) * t53) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, g(3) * t45 + t20 * t42, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t7, 0, t49 * pkin(3), 0, 0, 0, 0, 0, 0, -t6, t5, -t7, -g(1) * (t26 + t62) - g(3) * (t57 - t75) - (t55 - t83) * t79, 0, 0, 0, 0, 0, 0, -t4, t3, -t7, -g(1) * (t26 + t56) - g(2) * (t19 + (-t74 - t83) * t46) - g(3) * (-t53 - t75); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, -t7, -g(1) * t62 - g(3) * t57 - t55 * t79, 0, 0, 0, 0, 0, 0, -t4, t3, -t7, -g(1) * t56 - g(2) * (-t46 * t74 + t19) + g(3) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, g(1) * t15 - g(2) * t17 + t44 * t77, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t84 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t8;
