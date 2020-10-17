% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:19:27
% EndTime: 2019-05-07 18:19:28
% DurationCPUTime: 0.36s
% Computational Cost: add. (438->89), mult. (484->132), div. (0->0), fcn. (508->10), ass. (0->64)
t49 = qJ(3) + qJ(4);
t42 = cos(t49);
t53 = cos(qJ(3));
t31 = t53 * pkin(3) + pkin(4) * t42;
t29 = pkin(2) + t31;
t54 = cos(qJ(2));
t24 = t54 * t29;
t48 = -qJ(5) - pkin(9) - pkin(8);
t51 = sin(qJ(2));
t86 = -t51 * t48 + t24;
t85 = -pkin(1) - t86;
t52 = sin(qJ(1));
t55 = cos(qJ(1));
t63 = g(1) * t55 + g(2) * t52;
t15 = -g(3) * t54 + t63 * t51;
t41 = sin(t49);
t84 = pkin(4) * t41;
t40 = pkin(10) + t49;
t37 = sin(t40);
t83 = pkin(5) * t37;
t80 = g(3) * t51;
t38 = cos(t40);
t78 = t38 * t51;
t76 = t52 * t42;
t75 = t52 * t54;
t50 = sin(qJ(3));
t30 = t50 * pkin(3) + t84;
t28 = t55 * t30;
t74 = t55 * t37;
t73 = t55 * t38;
t72 = t55 * t41;
t71 = t55 * t42;
t70 = t55 * t50;
t69 = t55 * t53;
t68 = -t54 * t28 + t52 * t31;
t67 = t54 * t72;
t7 = t37 * t75 + t73;
t8 = t38 * t75 - t74;
t66 = -t7 * pkin(5) + t8 * qJ(6);
t10 = t52 * t37 + t54 * t73;
t9 = -t52 * t38 + t54 * t74;
t65 = -t9 * pkin(5) + t10 * qJ(6);
t64 = -t30 * t75 - t55 * t31;
t62 = g(1) * t52 - g(2) * t55;
t61 = pkin(5) * t38 + qJ(6) * t37;
t11 = t41 * t75 + t71;
t58 = (pkin(7) + t30) * t52 - t85 * t55;
t1 = g(1) * t9 + g(2) * t7 + t37 * t80;
t57 = t55 * pkin(7) + t85 * t52 + t28;
t16 = t63 * t54 + t80;
t36 = pkin(4) * t76;
t32 = qJ(6) * t78;
t25 = t62 * t51;
t23 = t52 * t50 + t54 * t69;
t22 = t52 * t53 - t54 * t70;
t21 = -t53 * t75 + t70;
t20 = t50 * t75 + t69;
t14 = t52 * t41 + t54 * t71;
t13 = -t67 + t76;
t12 = -t42 * t75 + t72;
t4 = g(1) * t14 - g(2) * t12 + t42 * t80;
t3 = -g(1) * t13 + g(2) * t11 + t41 * t80;
t2 = -g(1) * t10 - g(2) * t8 - g(3) * t78;
t5 = [0, t62, t63, 0, 0, 0, 0, 0, t62 * t54, -t25, 0, 0, 0, 0, 0, -g(1) * t21 - g(2) * t23, -g(1) * t20 - g(2) * t22, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t14, -g(1) * t11 - g(2) * t13, t25, -g(1) * t57 - g(2) * t58, g(1) * t8 - g(2) * t10, t25, g(1) * t7 - g(2) * t9, -g(1) * (-t8 * pkin(5) - t7 * qJ(6) + t57) - g(2) * (t10 * pkin(5) + t9 * qJ(6) + t58); 0, 0, 0, 0, 0, 0, 0, 0, t15, t16, 0, 0, 0, 0, 0, t15 * t53, -t15 * t50, 0, 0, 0, 0, 0, t15 * t42, -t15 * t41, -t16, -g(3) * t86 + t63 * (t29 * t51 + t48 * t54) t15 * t38, -t16, t15 * t37, -g(3) * t24 + (-g(3) * t61 + t63 * t48) * t54 + (g(3) * t48 + t63 * (t29 + t61)) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t22 + g(2) * t20 + t50 * t80, g(1) * t23 - g(2) * t21 + t53 * t80, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * t68 - g(2) * t64 + t30 * t80, t1, 0, t2, -g(1) * (t65 + t68) - g(2) * (t64 + t66) - g(3) * (t32 + (-t30 - t83) * t51); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * t36 + (g(2) * t71 + t16 * t41) * pkin(4), t1, 0, t2, -g(1) * (-pkin(4) * t67 + t36 + t65) - g(2) * (-t11 * pkin(4) + t66) - g(3) * (t32 + (-t83 - t84) * t51); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, 0, 0, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;
