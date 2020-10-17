% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:43:52
% EndTime: 2019-05-05 21:43:53
% DurationCPUTime: 0.43s
% Computational Cost: add. (270->96), mult. (442->126), div. (0->0), fcn. (451->8), ass. (0->63)
t41 = sin(qJ(3));
t40 = sin(qJ(4));
t45 = cos(qJ(1));
t65 = t45 * t40;
t42 = sin(qJ(1));
t43 = cos(qJ(4));
t70 = t42 * t43;
t13 = t41 * t65 + t70;
t64 = t45 * t43;
t71 = t42 * t40;
t11 = -t41 * t71 + t64;
t78 = g(2) * t45;
t79 = g(1) * t42;
t82 = t78 - t79;
t44 = cos(qJ(3));
t10 = -g(3) * t41 - t44 * t82;
t81 = -pkin(1) - pkin(7);
t80 = pkin(4) * t40;
t76 = g(3) * t44;
t75 = t41 * t42;
t74 = t41 * t45;
t38 = qJ(4) + pkin(9);
t31 = sin(t38);
t73 = t42 * t31;
t32 = cos(t38);
t72 = t42 * t32;
t39 = -qJ(5) - pkin(8);
t69 = t44 * t39;
t68 = t44 * t45;
t67 = t45 * t31;
t66 = t45 * t32;
t63 = t13 * pkin(4);
t62 = t45 * pkin(1) + t42 * qJ(2);
t61 = t40 * t76;
t30 = t43 * pkin(4) + pkin(3);
t34 = t45 * qJ(2);
t58 = t30 * t74 + t39 * t68 + t34;
t57 = t45 * pkin(7) + t62;
t56 = g(2) * t57;
t5 = t41 * t73 - t66;
t7 = t41 * t67 + t72;
t55 = g(1) * t7 + g(2) * t5;
t53 = t41 * pkin(3) - t44 * pkin(8);
t21 = g(1) * t45 + g(2) * t42;
t52 = pkin(5) * t32 + qJ(6) * t31;
t51 = t11 * pkin(4);
t50 = pkin(4) * t65 + t30 * t75 + t42 * t69 + t57;
t49 = t30 + t52;
t48 = (-t80 + t81) * t79;
t1 = g(1) * t5 - g(2) * t7 + t31 * t76;
t6 = t41 * t72 + t67;
t8 = t41 * t66 - t73;
t47 = g(1) * t6 - g(2) * t8 + t32 * t76;
t23 = t39 * t74;
t17 = t42 * t44 * t30;
t15 = t21 * t44;
t14 = t41 * t64 - t71;
t12 = t41 * t70 + t65;
t9 = g(1) * t75 - g(2) * t74 + t76;
t4 = t10 * t32;
t3 = t10 * t31;
t2 = -g(1) * t8 - g(2) * t6;
t16 = [0, 0, 0, 0, 0, 0, -t82, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t21, -g(1) * (-t42 * pkin(1) + t34) - g(2) * t62, 0, 0, 0, 0, 0, 0, -t21 * t41, -t15, -t82, -g(1) * (t42 * t81 + t34) - t56, 0, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t12, g(1) * t13 - g(2) * t11, t15, -g(1) * (pkin(3) * t74 - pkin(8) * t68 + t34) - t56 + (-g(1) * t81 - g(2) * t53) * t42, 0, 0, 0, 0, 0, 0, t2, t55, t15, -g(1) * t58 - g(2) * t50 - t48, 0, 0, 0, 0, 0, 0, t2, t15, -t55, -g(1) * (t8 * pkin(5) + t7 * qJ(6) + t58) - g(2) * (t6 * pkin(5) + t5 * qJ(6) + t50) - t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0, 0, 0, 0, 0, 0, 0, 0, -t10 * t43, t10 * t40, -t9, g(3) * t53 + t82 * (pkin(3) * t44 + pkin(8) * t41) 0, 0, 0, 0, 0, 0, -t4, t3, -t9, -g(1) * (-t39 * t75 + t17) - g(2) * (-t30 * t68 + t23) - g(3) * (-t41 * t30 - t69) 0, 0, 0, 0, 0, 0, -t4, -t9, -t3, -g(1) * t17 - g(2) * t23 + (g(3) * t49 + t39 * t79) * t41 + (g(3) * t39 + t49 * t78 - t52 * t79) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t13 + t61, g(1) * t12 - g(2) * t14 + t43 * t76, 0, 0, 0, 0, 0, 0, 0, 0, t1, t47, 0, pkin(4) * t61 - g(1) * t51 - g(2) * t63, 0, 0, 0, 0, 0, 0, t1, 0, -t47, -g(1) * (-t5 * pkin(5) + t6 * qJ(6) + t51) - g(2) * (t7 * pkin(5) - t8 * qJ(6) + t63) - (-pkin(5) * t31 + qJ(6) * t32 - t80) * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t16;
