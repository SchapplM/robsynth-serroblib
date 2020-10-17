% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP13_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP13_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:13:51
% EndTime: 2019-05-06 19:13:52
% DurationCPUTime: 0.46s
% Computational Cost: add. (293->97), mult. (734->153), div. (0->0), fcn. (889->10), ass. (0->58)
t37 = sin(qJ(2));
t41 = cos(qJ(2));
t42 = cos(qJ(1));
t38 = sin(qJ(1));
t61 = cos(pkin(6));
t58 = t38 * t61;
t23 = -t37 * t58 + t42 * t41;
t35 = sin(qJ(5));
t39 = cos(qJ(5));
t22 = t42 * t37 + t41 * t58;
t36 = sin(qJ(4));
t40 = cos(qJ(4));
t33 = sin(pkin(6));
t68 = t33 * t38;
t6 = t22 * t36 + t40 * t68;
t1 = t23 * t39 - t35 * t6;
t67 = t33 * t41;
t19 = -t36 * t67 + t61 * t40;
t62 = t37 * t39;
t57 = t42 * t61;
t21 = t37 * t57 + t38 * t41;
t20 = t38 * t37 - t41 * t57;
t66 = t33 * t42;
t50 = -t20 * t36 + t40 * t66;
t79 = t21 * t39 + t35 * t50;
t81 = -g(2) * t79 - g(3) * (-t19 * t35 + t33 * t62) - g(1) * t1;
t80 = -g(1) * t23 - g(2) * t21;
t78 = -t21 * t35 + t39 * t50;
t74 = g(3) * t33;
t73 = t20 * t35;
t70 = t22 * t35;
t69 = t33 * t37;
t65 = t35 * t36;
t64 = t35 * t37;
t63 = t36 * t39;
t60 = pkin(5) * t35 + pkin(9);
t59 = g(3) * (pkin(2) * t67 + qJ(3) * t69);
t49 = t20 * t40 + t36 * t66;
t5 = -t22 * t40 + t36 * t68;
t56 = g(1) * t49 + g(2) * t5;
t55 = g(1) * t20 - g(2) * t22;
t54 = g(1) * t21 - g(2) * t23;
t53 = g(1) * t42 + g(2) * t38;
t31 = pkin(5) * t39 + pkin(4);
t34 = -qJ(6) - pkin(10);
t52 = t31 * t36 + t34 * t40;
t51 = t42 * pkin(1) + t23 * pkin(2) + pkin(8) * t68 + t22 * qJ(3);
t18 = t61 * t36 + t40 * t67;
t47 = g(1) * t5 - g(2) * t49 + g(3) * t18;
t46 = g(1) * t6 - g(2) * t50 + g(3) * t19;
t45 = -t38 * pkin(1) - t21 * pkin(2) + pkin(8) * t66 - t20 * qJ(3);
t4 = -g(1) * t22 - g(2) * t20 + g(3) * t67;
t44 = g(3) * t69 - t80;
t16 = t22 * pkin(2);
t14 = t20 * pkin(2);
t3 = t44 * t40;
t2 = t23 * t35 + t39 * t6;
t7 = [0, g(1) * t38 - g(2) * t42, t53, 0, 0, 0, 0, 0, t54, -t55, -t53 * t33, -t54, t55, -g(1) * t45 - g(2) * t51, 0, 0, 0, 0, 0, -g(1) * t50 - g(2) * t6, t56, 0, 0, 0, 0, 0, -g(1) * t78 - g(2) * t2, g(1) * t79 - g(2) * t1, -t56, -g(1) * (pkin(3) * t66 - t60 * t21 + t31 * t50 - t34 * t49 + t45) - g(2) * (pkin(3) * t68 + t60 * t23 + t6 * t31 - t5 * t34 + t51); 0, 0, 0, 0, 0, 0, 0, 0, -t4, t44, 0, t4, -t44, -g(1) * (t23 * qJ(3) - t16) - g(2) * (t21 * qJ(3) - t14) - t59, 0, 0, 0, 0, 0, -t44 * t36, -t3, 0, 0, 0, 0, 0, -g(1) * (t23 * t63 - t70) - g(2) * (t21 * t63 - t73) - (t35 * t41 + t36 * t62) * t74, -g(1) * (-t22 * t39 - t23 * t65) - g(2) * (-t20 * t39 - t21 * t65) - (-t36 * t64 + t39 * t41) * t74, t3, -g(1) * (-pkin(5) * t70 - t22 * pkin(9) - t16) - g(2) * (-pkin(5) * t73 - t20 * pkin(9) - t14) - t59 - (t52 * t37 + t60 * t41) * t74 + t80 * (qJ(3) + t52); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t46, 0, 0, 0, 0, 0, t47 * t39, -t47 * t35, -t46, -g(1) * (-t31 * t5 - t34 * t6) - g(2) * (t31 * t49 + t34 * t50) - g(3) * (-t18 * t31 - t19 * t34); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, g(1) * t2 - g(2) * t78 - g(3) * (-t19 * t39 - t33 * t64) 0, t81 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47;];
taug_reg  = t7;
