% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [6x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:47:32
% EndTime: 2019-05-05 03:47:33
% DurationCPUTime: 0.36s
% Computational Cost: add. (434->100), mult. (807->157), div. (0->0), fcn. (991->12), ass. (0->69)
t47 = sin(pkin(10));
t49 = cos(pkin(10));
t56 = cos(qJ(2));
t53 = sin(qJ(2));
t70 = cos(pkin(6));
t65 = t53 * t70;
t33 = t47 * t56 + t49 * t65;
t52 = sin(qJ(3));
t48 = sin(pkin(6));
t55 = cos(qJ(3));
t76 = t48 * t55;
t66 = t49 * t76;
t88 = -t33 * t52 - t66;
t35 = -t47 * t65 + t49 * t56;
t46 = qJ(3) + pkin(11);
t44 = sin(t46);
t45 = cos(t46);
t77 = t48 * t53;
t79 = t48 * t49;
t80 = t47 * t48;
t58 = g(3) * (-t44 * t77 + t70 * t45) + g(2) * (-t33 * t44 - t45 * t79) + g(1) * (-t35 * t44 + t45 * t80);
t83 = t35 * t52;
t51 = sin(qJ(5));
t82 = t45 * t51;
t54 = cos(qJ(5));
t81 = t45 * t54;
t78 = t48 * t52;
t75 = t48 * t56;
t50 = -qJ(4) - pkin(8);
t74 = t50 * t53;
t73 = t54 * t56;
t64 = t56 * t70;
t32 = t47 * t53 - t49 * t64;
t43 = t55 * pkin(3) + pkin(2);
t72 = -t32 * t43 - t33 * t50;
t34 = t47 * t64 + t49 * t53;
t71 = -t34 * t43 - t35 * t50;
t69 = t47 * t76;
t68 = t52 * t77;
t67 = t51 * t75;
t63 = t70 * t55;
t62 = pkin(4) * t45 + pkin(9) * t44;
t25 = t70 * t44 + t45 * t77;
t18 = t25 * t51 + t48 * t73;
t15 = t33 * t45 - t44 * t79;
t5 = t15 * t51 - t32 * t54;
t17 = t35 * t45 + t44 * t80;
t7 = t17 * t51 - t34 * t54;
t1 = g(1) * t7 + g(2) * t5 + g(3) * t18;
t19 = t25 * t54 - t67;
t6 = t15 * t54 + t32 * t51;
t8 = t17 * t54 + t34 * t51;
t60 = g(1) * t8 + g(2) * t6 + g(3) * t19;
t11 = -t34 * t82 - t35 * t54;
t20 = t45 * t67 - t54 * t77;
t9 = -t32 * t82 - t33 * t54;
t59 = g(1) * t11 + g(2) * t9 + g(3) * t20;
t13 = -g(1) * t34 - g(2) * t32 + g(3) * t75;
t57 = g(1) * t35 + g(2) * t33 + g(3) * t77;
t42 = pkin(3) * t63;
t37 = pkin(3) * t69;
t36 = t43 * t75;
t21 = (t45 * t73 + t51 * t53) * t48;
t12 = -t34 * t81 + t35 * t51;
t10 = -t32 * t81 + t33 * t51;
t4 = t58 * t54;
t3 = t58 * t51;
t2 = -g(1) * t12 - g(2) * t10 - g(3) * t21;
t14 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -t13, t57, 0, 0, 0, 0, 0, -t13 * t55, t13 * t52, -t57, -g(1) * t71 - g(2) * t72 - g(3) * (-t48 * t74 + t36) 0, 0, 0, 0, 0, t2, t59, t2, -t13 * t44, -t59, -g(1) * (t12 * pkin(5) + t11 * qJ(6) - t62 * t34 + t71) - g(2) * (t10 * pkin(5) + t9 * qJ(6) - t62 * t32 + t72) + (-t21 * pkin(5) - t20 * qJ(6) - t36 - (t62 * t56 - t74) * t48) * g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t69 - t83) - g(2) * t88 - g(3) * (t63 - t68) -g(1) * (-t35 * t55 - t47 * t78) - g(2) * (-t33 * t55 + t49 * t78) - g(3) * (-t70 * t52 - t53 * t76) 0, -g(1) * t37 - g(3) * t42 + (g(2) * t66 + t57 * t52) * pkin(3), 0, 0, 0, 0, 0, -t4, t3, -t4, -g(1) * t17 - g(2) * t15 - g(3) * t25, -t3, -g(1) * (-pkin(3) * t83 + t17 * pkin(9) + t37) - g(2) * (t88 * pkin(3) + t15 * pkin(9)) - g(3) * (-pkin(3) * t68 + t25 * pkin(9) + t42) - t58 * (pkin(5) * t54 + qJ(6) * t51 + pkin(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t60, t1, 0, -t60, -g(1) * (-t7 * pkin(5) + t8 * qJ(6)) - g(2) * (-t5 * pkin(5) + t6 * qJ(6)) - g(3) * (-t18 * pkin(5) + t19 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t14;
