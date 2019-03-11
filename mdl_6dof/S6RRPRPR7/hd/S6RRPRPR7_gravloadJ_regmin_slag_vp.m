% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t41 = sin(qJ(6));
t39 = qJ(4) + pkin(10);
t31 = sin(t39);
t32 = cos(t39);
t44 = sin(qJ(1));
t43 = sin(qJ(2));
t47 = cos(qJ(2));
t54 = t43 * t31 + t47 * t32;
t48 = cos(qJ(1));
t62 = t47 * t48;
t64 = t47 * t31;
t66 = t43 * t48;
t67 = t43 * t44;
t51 = g(1) * (t31 * t62 - t32 * t66) + g(2) * (-t32 * t67 + t44 * t64) + g(3) * t54;
t76 = t51 * t41;
t45 = cos(qJ(6));
t75 = t51 * t45;
t33 = t43 * qJ(3);
t61 = t47 * pkin(2) + t33;
t22 = g(1) * t48 + g(2) * t44;
t74 = g(1) * t44;
t71 = g(3) * (t43 * t32 - t64);
t46 = cos(qJ(4));
t42 = sin(qJ(4));
t68 = t43 * t42;
t17 = t47 * t46 + t68;
t70 = g(3) * t17;
t30 = t46 * pkin(4) + pkin(3);
t65 = t47 * t30;
t63 = t47 * t42;
t60 = qJ(3) * t47;
t59 = pkin(4) * t68;
t58 = t42 * t62;
t57 = pkin(2) * t62 + t44 * pkin(7) + (pkin(1) + t33) * t48;
t21 = -g(2) * t48 + t74;
t4 = t54 * t44;
t56 = t4 * t45 + t48 * t41;
t55 = t4 * t41 - t48 * t45;
t53 = -t43 * t46 + t63;
t52 = -pkin(1) - t61;
t11 = -t46 * t66 + t58;
t9 = t53 * t44;
t50 = g(1) * t11 + g(2) * t9 + t70;
t10 = t17 * t44;
t12 = t17 * t48;
t49 = g(1) * t12 + g(2) * t10 - g(3) * t53;
t40 = -qJ(5) - pkin(8);
t36 = t48 * pkin(7);
t28 = t48 * t60;
t26 = t44 * t60;
t16 = t21 * t47;
t15 = t21 * t43;
t8 = g(3) * t43 + t22 * t47;
t7 = -g(3) * t47 + t22 * t43;
t6 = t54 * t48;
t2 = -t44 * t41 + t6 * t45;
t1 = -t6 * t41 - t44 * t45;
t3 = [0, t21, t22, 0, 0, 0, 0, 0, t16, -t15, t16, -t22, t15, -g(1) * t36 - g(2) * t57 - t52 * t74, 0, 0, 0, 0, 0, g(1) * t10 - g(2) * t12, -g(1) * t9 + g(2) * t11, t22, -g(1) * (t48 * t40 + t36) - g(2) * (t30 * t62 + t48 * t59 + t57) + (-g(1) * (t52 - t59 - t65) - g(2) * t40) * t44, 0, 0, 0, 0, 0, g(1) * t56 - g(2) * t2, -g(1) * t55 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, t7, 0, -t8, -g(1) * (-pkin(2) * t66 + t28) - g(2) * (-pkin(2) * t67 + t26) - g(3) * t61, 0, 0, 0, 0, 0, -t50, -t49, 0, -g(1) * (pkin(4) * t58 + t28) - g(2) * (t44 * pkin(4) * t63 + t26) - g(3) * (t61 + t65) + (-g(3) * pkin(4) * t42 + t22 * (pkin(2) + t30)) * t43, 0, 0, 0, 0, 0, -t75, t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t49, 0 (t22 * t53 + t70) * pkin(4), 0, 0, 0, 0, 0, t75, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t55 + t41 * t71, g(1) * t2 + g(2) * t56 + t45 * t71;];
taug_reg  = t3;
