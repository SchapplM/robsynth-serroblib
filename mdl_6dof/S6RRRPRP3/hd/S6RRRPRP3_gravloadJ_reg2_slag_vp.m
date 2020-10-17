% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:39:45
% EndTime: 2019-05-07 07:39:47
% DurationCPUTime: 0.52s
% Computational Cost: add. (529->104), mult. (563->133), div. (0->0), fcn. (559->10), ass. (0->67)
t38 = pkin(10) + qJ(5);
t33 = sin(t38);
t34 = cos(t38);
t88 = pkin(5) * t34 + qJ(6) * t33;
t44 = sin(qJ(1));
t46 = cos(qJ(1));
t20 = g(1) * t46 + g(2) * t44;
t39 = qJ(2) + qJ(3);
t35 = sin(t39);
t36 = cos(t39);
t7 = -g(3) * t36 + t20 * t35;
t87 = t88 * t36;
t86 = pkin(3) * t36 + qJ(4) * t35;
t43 = sin(qJ(2));
t85 = pkin(2) * t43;
t81 = g(3) * t35;
t79 = t35 * t44;
t78 = t35 * t46;
t41 = cos(pkin(10));
t29 = pkin(4) * t41 + pkin(3);
t16 = t36 * t29;
t42 = -pkin(9) - qJ(4);
t77 = t36 * t42;
t76 = t44 * t33;
t75 = t44 * t34;
t40 = sin(pkin(10));
t74 = t44 * t40;
t73 = t44 * t41;
t72 = t46 * t33;
t71 = t46 * t34;
t70 = t46 * t40;
t69 = t46 * t41;
t47 = -pkin(8) - pkin(7);
t68 = t46 * t47;
t66 = qJ(4) * t36;
t64 = -t35 * t42 + t16;
t45 = cos(qJ(2));
t37 = t45 * pkin(2);
t32 = t37 + pkin(1);
t23 = t46 * t32;
t63 = -t44 * t47 + t23;
t62 = t37 + t64;
t11 = t36 * t72 - t75;
t9 = t36 * t76 + t71;
t61 = g(1) * t9 - g(2) * t11;
t60 = -pkin(3) * t35 - t85;
t59 = g(1) * t44 - g(2) * t46;
t57 = t29 * t35 + t77;
t56 = t29 + t88;
t54 = t20 * t36;
t1 = g(1) * t11 + g(2) * t9 + t33 * t81;
t10 = t36 * t75 - t72;
t12 = t36 * t71 + t76;
t52 = g(1) * t12 + g(2) * t10 + t34 * t81;
t51 = pkin(4) * t74 + t16 * t46 - t42 * t78 + t63;
t50 = -g(3) * t45 + t20 * t43;
t49 = -t68 + t42 * t79 + pkin(4) * t70 + (-t32 - t16) * t44;
t22 = t46 * t66;
t21 = t44 * t66;
t13 = t59 * t35;
t8 = t54 + t81;
t6 = t7 * t41;
t5 = t7 * t40;
t4 = t7 * t34;
t3 = t7 * t33;
t2 = g(1) * t10 - g(2) * t12;
t14 = [0, 0, 0, 0, 0, 0, t59, t20, 0, 0, 0, 0, 0, 0, 0, 0, t59 * t45, -t59 * t43, -t20, -g(1) * (-pkin(1) * t44 + pkin(7) * t46) - g(2) * (pkin(1) * t46 + pkin(7) * t44) 0, 0, 0, 0, 0, 0, t59 * t36, -t13, -t20, -g(1) * (-t32 * t44 - t68) - g(2) * t63, 0, 0, 0, 0, 0, 0, -g(1) * (-t36 * t73 + t70) - g(2) * (t36 * t69 + t74) -g(1) * (t36 * t74 + t69) - g(2) * (-t36 * t70 + t73) t13, -g(2) * t23 + (g(1) * t47 - g(2) * t86) * t46 + (-g(1) * (-t32 - t86) + g(2) * t47) * t44, 0, 0, 0, 0, 0, 0, t2, -t61, t13, -g(1) * t49 - g(2) * t51, 0, 0, 0, 0, 0, 0, t2, t13, t61, -g(1) * (-t10 * pkin(5) - t9 * qJ(6) + t49) - g(2) * (pkin(5) * t12 + qJ(6) * t11 + t51); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, g(3) * t43 + t20 * t45, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t50 * pkin(2), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t46 * t60 + t22) - g(2) * (t44 * t60 + t21) - g(3) * (t37 + t86) 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(3) * t62 + t20 * (t57 + t85) 0, 0, 0, 0, 0, 0, t4, -t8, t3, -g(3) * (t62 + t87) + t20 * (t35 * t56 + t77 + t85); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (-pkin(3) * t78 + t22) - g(2) * (-pkin(3) * t79 + t21) - g(3) * t86, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(3) * t64 + t20 * t57, 0, 0, 0, 0, 0, 0, t4, -t8, t3, -g(3) * (t16 + t87) + t42 * t54 + (g(3) * t42 + t20 * t56) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t52, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t52, -g(1) * (-pkin(5) * t11 + qJ(6) * t12) - g(2) * (-pkin(5) * t9 + qJ(6) * t10) - (-pkin(5) * t33 + qJ(6) * t34) * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t14;
