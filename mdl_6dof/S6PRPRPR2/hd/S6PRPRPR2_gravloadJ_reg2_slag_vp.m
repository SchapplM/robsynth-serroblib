% Calculate inertial parameters regressor of gravitation load for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:23:34
% EndTime: 2019-05-04 22:23:36
% DurationCPUTime: 0.62s
% Computational Cost: add. (543->119), mult. (1321->186), div. (0->0), fcn. (1681->14), ass. (0->65)
t40 = sin(pkin(11));
t47 = sin(qJ(2));
t49 = cos(qJ(2));
t74 = cos(pkin(11));
t25 = -t49 * t40 - t47 * t74;
t41 = sin(pkin(10));
t43 = cos(pkin(10));
t59 = -t47 * t40 + t49 * t74;
t44 = cos(pkin(6));
t76 = t25 * t44;
t15 = t41 * t76 + t43 * t59;
t10 = -t41 * t59 + t43 * t76;
t38 = pkin(12) + qJ(6);
t36 = sin(t38);
t48 = cos(qJ(4));
t86 = t36 * t48;
t37 = cos(t38);
t85 = t37 * t48;
t39 = sin(pkin(12));
t84 = t39 * t48;
t82 = t41 * t47;
t42 = cos(pkin(12));
t81 = t42 * t48;
t79 = t44 * t47;
t78 = t44 * t49;
t73 = sin(pkin(6));
t61 = t74 * t73;
t68 = t47 * t73;
t22 = t40 * t68 - t49 * t61;
t66 = t49 * t73;
t75 = pkin(2) * t66 - t22 * pkin(3);
t72 = t43 * t78;
t71 = pkin(5) * t39 + pkin(8);
t46 = sin(qJ(4));
t70 = t46 * t73;
t67 = t48 * t73;
t23 = t40 * t66 + t47 * t61;
t65 = t23 * pkin(8) + t75;
t64 = pkin(4) * t48 + qJ(5) * t46;
t35 = t42 * pkin(5) + pkin(4);
t45 = -pkin(9) - qJ(5);
t63 = t35 * t48 - t45 * t46;
t53 = t59 * t44;
t11 = t41 * t25 + t43 * t53;
t26 = pkin(2) * t72;
t62 = -pkin(2) * t82 + t11 * pkin(3) + t26;
t60 = -t41 * t78 - t43 * t47;
t16 = t23 * t46 - t44 * t48;
t4 = -t10 * t46 + t43 * t67;
t6 = t15 * t46 - t41 * t67;
t58 = g(1) * t6 + g(2) * t4 + g(3) * t16;
t17 = t23 * t48 + t44 * t46;
t5 = -t10 * t48 - t43 * t70;
t7 = t15 * t48 + t41 * t70;
t57 = g(1) * t7 + g(2) * t5 + g(3) * t17;
t56 = -g(1) * t15 + g(2) * t10 - g(3) * t23;
t14 = t43 * t25 - t41 * t53;
t55 = g(1) * t14 + g(2) * t11 - g(3) * t22;
t54 = -t10 * pkin(8) + t62;
t52 = pkin(2) * t60 + t14 * pkin(3);
t51 = pkin(8) * t15 + t52;
t50 = -g(1) * t60 - g(3) * t66;
t21 = -g(3) * t44 + (-g(1) * t41 + g(2) * t43) * t73;
t3 = t55 * t46;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * (t72 - t82) + t50, -g(1) * (t41 * t79 - t43 * t49) - g(2) * (-t41 * t49 - t43 * t79) + g(3) * t68, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t56, 0, -g(2) * t26 + (g(2) * t82 + t50) * pkin(2), 0, 0, 0, 0, 0, 0, -t55 * t48, t3, t56, -g(1) * t51 - g(2) * t54 - g(3) * t65, 0, 0, 0, 0, 0, 0, -g(1) * (t14 * t81 + t15 * t39) - g(2) * (-t10 * t39 + t11 * t81) - g(3) * (-t22 * t81 + t23 * t39) -g(1) * (-t14 * t84 + t15 * t42) - g(2) * (-t10 * t42 - t11 * t84) - g(3) * (t22 * t84 + t23 * t42) -t3, -g(1) * (t64 * t14 + t51) - g(2) * (t64 * t11 + t54) - g(3) * (-t64 * t22 + t65) 0, 0, 0, 0, 0, 0, -g(1) * (t14 * t85 + t15 * t36) - g(2) * (-t10 * t36 + t11 * t85) - g(3) * (-t22 * t85 + t23 * t36) -g(1) * (-t14 * t86 + t15 * t37) - g(2) * (-t10 * t37 - t11 * t86) - g(3) * (t22 * t86 + t23 * t37) -t3, -g(1) * (t63 * t14 + t15 * t71 + t52) - g(2) * (-t71 * t10 + t63 * t11 + t62) - g(3) * (-t63 * t22 + t71 * t23 + t75); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t57, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t42, -t58 * t39, -t57, -g(1) * (-t6 * pkin(4) + t7 * qJ(5)) - g(2) * (-t4 * pkin(4) + t5 * qJ(5)) - g(3) * (-t16 * pkin(4) + t17 * qJ(5)) 0, 0, 0, 0, 0, 0, t58 * t37, -t58 * t36, -t57, -g(1) * (-t6 * t35 - t7 * t45) - g(2) * (-t4 * t35 - t5 * t45) - g(3) * (-t16 * t35 - t17 * t45); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t37 - t7 * t36) - g(2) * (-t11 * t37 - t5 * t36) - g(3) * (-t17 * t36 + t22 * t37) -g(1) * (t14 * t36 - t7 * t37) - g(2) * (t11 * t36 - t5 * t37) - g(3) * (-t17 * t37 - t22 * t36) 0, 0;];
taug_reg  = t1;
