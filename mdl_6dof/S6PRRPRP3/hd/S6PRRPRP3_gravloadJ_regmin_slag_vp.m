% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRP3
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
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t54 = sin(qJ(2));
t56 = cos(qJ(2));
t72 = cos(pkin(10));
t73 = cos(pkin(6));
t64 = t73 * t72;
t71 = sin(pkin(10));
t31 = t54 * t64 + t71 * t56;
t53 = sin(qJ(3));
t55 = cos(qJ(3));
t50 = sin(pkin(6));
t68 = t50 * t72;
t18 = t31 * t53 + t55 * t68;
t63 = t73 * t71;
t33 = -t54 * t63 + t72 * t56;
t67 = t50 * t71;
t20 = t33 * t53 - t55 * t67;
t78 = t50 * t54;
t34 = t53 * t78 - t73 * t55;
t59 = g(1) * t20 + g(2) * t18 + g(3) * t34;
t83 = g(3) * t50;
t48 = pkin(11) + qJ(5);
t46 = sin(t48);
t82 = t46 * t55;
t47 = cos(t48);
t81 = t47 * t55;
t49 = sin(pkin(11));
t80 = t49 * t54;
t79 = t49 * t55;
t77 = t50 * t56;
t51 = cos(pkin(11));
t76 = t51 * t55;
t75 = t55 * t56;
t74 = pkin(2) * t77 + pkin(8) * t78;
t70 = t46 * t77;
t69 = pkin(4) * t49 + pkin(8);
t66 = pkin(3) * t55 + qJ(4) * t53;
t45 = t51 * pkin(4) + pkin(3);
t52 = -pkin(9) - qJ(4);
t65 = t45 * t55 - t52 * t53;
t35 = t73 * t53 + t55 * t78;
t16 = t35 * t46 + t47 * t77;
t19 = t31 * t55 - t53 * t68;
t30 = t71 * t54 - t56 * t64;
t7 = t19 * t46 - t30 * t47;
t21 = t33 * t55 + t53 * t67;
t32 = t72 * t54 + t56 * t63;
t9 = t21 * t46 - t32 * t47;
t1 = g(1) * t9 + g(2) * t7 + g(3) * t16;
t10 = t21 * t47 + t32 * t46;
t17 = t35 * t47 - t70;
t8 = t19 * t47 + t30 * t46;
t61 = g(1) * t10 + g(2) * t8 + g(3) * t17;
t12 = -t30 * t82 - t31 * t47;
t14 = -t32 * t82 - t33 * t47;
t22 = -t47 * t78 + t55 * t70;
t60 = g(1) * t14 + g(2) * t12 + g(3) * t22;
t58 = g(1) * t21 + g(2) * t19 + g(3) * t35;
t57 = -g(1) * t32 - g(2) * t30 + g(3) * t77;
t29 = t32 * pkin(2);
t28 = t30 * pkin(2);
t23 = (t46 * t54 + t47 * t75) * t50;
t15 = -t32 * t81 + t33 * t46;
t13 = -t30 * t81 + t31 * t46;
t11 = t57 * t53;
t4 = t59 * t47;
t3 = t59 * t46;
t2 = -g(1) * t15 - g(2) * t13 - g(3) * t23;
t5 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -t57, g(1) * t33 + g(2) * t31 + g(3) * t78, 0, 0, 0, 0, 0, -t57 * t55, t11, -g(1) * (-t32 * t76 + t33 * t49) - g(2) * (-t30 * t76 + t31 * t49) - (t51 * t75 + t80) * t83, -g(1) * (t32 * t79 + t33 * t51) - g(2) * (t30 * t79 + t31 * t51) - (-t49 * t75 + t51 * t54) * t83, -t11, -g(1) * (t33 * pkin(8) - t66 * t32 - t29) - g(2) * (t31 * pkin(8) - t66 * t30 - t28) - g(3) * (t66 * t77 + t74) 0, 0, 0, 0, 0, t2, t60, t2, -t11, -t60, -g(1) * (t15 * pkin(5) + t14 * qJ(6) - t65 * t32 + t69 * t33 - t29) - g(2) * (t13 * pkin(5) + t12 * qJ(6) - t65 * t30 + t69 * t31 - t28) - g(3) * (t23 * pkin(5) + t22 * qJ(6) + t74) - (pkin(4) * t80 + t65 * t56) * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t58, t59 * t51, -t59 * t49, -t58, -g(1) * (-t20 * pkin(3) + t21 * qJ(4)) - g(2) * (-t18 * pkin(3) + t19 * qJ(4)) - g(3) * (-t34 * pkin(3) + t35 * qJ(4)) 0, 0, 0, 0, 0, t4, -t3, t4, -t58, t3, t58 * t52 + t59 * (pkin(5) * t47 + qJ(6) * t46 + t45); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t61, t1, 0, -t61, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - g(3) * (-t16 * pkin(5) + t17 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;
