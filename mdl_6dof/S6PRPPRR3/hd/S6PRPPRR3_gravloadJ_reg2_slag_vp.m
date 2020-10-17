% Calculate inertial parameters regressor of gravitation load for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPPRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:01:35
% EndTime: 2019-05-04 22:01:37
% DurationCPUTime: 0.38s
% Computational Cost: add. (376->95), mult. (952->144), div. (0->0), fcn. (1203->12), ass. (0->55)
t46 = sin(pkin(6));
t50 = sin(qJ(5));
t80 = t46 * t50;
t51 = sin(qJ(2));
t79 = t46 * t51;
t53 = cos(qJ(5));
t78 = t46 * t53;
t54 = cos(qJ(2));
t77 = t46 * t54;
t49 = sin(qJ(6));
t76 = t49 * t53;
t52 = cos(qJ(6));
t75 = t52 * t53;
t74 = pkin(2) * t77 + qJ(3) * t79;
t73 = cos(pkin(6));
t72 = pkin(3) * t77 + t74;
t45 = sin(pkin(10));
t48 = cos(pkin(10));
t68 = t54 * t73;
t32 = t45 * t51 - t48 * t68;
t69 = t51 * t73;
t33 = t45 * t54 + t48 * t69;
t44 = sin(pkin(11));
t47 = cos(pkin(11));
t71 = -t32 * t47 + t33 * t44;
t13 = t32 * t44 + t33 * t47;
t34 = t45 * t68 + t48 * t51;
t35 = -t45 * t69 + t48 * t54;
t70 = -t34 * t47 + t35 * t44;
t17 = t34 * t44 + t35 * t47;
t67 = -t32 * pkin(2) + t33 * qJ(3);
t66 = -t34 * pkin(2) + t35 * qJ(3);
t65 = -t32 * pkin(3) + t67;
t64 = -t34 * pkin(3) + t66;
t63 = pkin(5) * t53 + pkin(9) * t50;
t26 = (t44 * t51 + t47 * t54) * t46;
t27 = -t44 * t77 + t47 * t79;
t62 = t26 * pkin(4) - t27 * pkin(8) + t72;
t18 = -t27 * t50 - t73 * t53;
t2 = -t13 * t50 + t48 * t78;
t4 = -t17 * t50 - t45 * t78;
t61 = g(1) * t4 + g(2) * t2 + g(3) * t18;
t19 = t27 * t53 - t73 * t50;
t3 = t13 * t53 + t48 * t80;
t5 = t17 * t53 - t45 * t80;
t60 = g(1) * t5 + g(2) * t3 + g(3) * t19;
t59 = -g(1) * t17 - g(2) * t13 - g(3) * t27;
t58 = g(1) * t70 + g(2) * t71 + g(3) * t26;
t57 = pkin(4) * t71 - pkin(8) * t13 + t65;
t56 = pkin(4) * t70 - pkin(8) * t17 + t64;
t6 = -g(1) * t34 - g(2) * t32 + g(3) * t77;
t55 = g(1) * t35 + g(2) * t33 + g(3) * t79;
t21 = g(3) * t73 + (g(1) * t45 - g(2) * t48) * t46;
t1 = t58 * t50;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t55, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, -t55, -g(1) * t66 - g(2) * t67 - g(3) * t74, 0, 0, 0, 0, 0, 0, -t58, t59, 0, -g(1) * t64 - g(2) * t65 - g(3) * t72, 0, 0, 0, 0, 0, 0, -t58 * t53, t1, -t59, -g(1) * t56 - g(2) * t57 - g(3) * t62, 0, 0, 0, 0, 0, 0, -g(1) * (-t17 * t49 + t70 * t75) - g(2) * (-t13 * t49 + t71 * t75) - g(3) * (t26 * t75 - t27 * t49) -g(1) * (-t17 * t52 - t70 * t76) - g(2) * (-t13 * t52 - t71 * t76) - g(3) * (-t26 * t76 - t27 * t52) -t1, -g(1) * (t63 * t70 + t56) - g(2) * (t63 * t71 + t57) - g(3) * (t63 * t26 + t62); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t60, 0, 0, 0, 0, 0, 0, 0, 0, -t61 * t52, t61 * t49, -t60, -g(1) * (t4 * pkin(5) + t5 * pkin(9)) - g(2) * (t2 * pkin(5) + t3 * pkin(9)) - g(3) * (t18 * pkin(5) + t19 * pkin(9)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t5 * t49 + t52 * t70) - g(2) * (-t3 * t49 + t52 * t71) - g(3) * (-t19 * t49 + t26 * t52) -g(1) * (-t49 * t70 - t5 * t52) - g(2) * (-t3 * t52 - t49 * t71) - g(3) * (-t19 * t52 - t26 * t49) 0, 0;];
taug_reg  = t7;
