% Calculate inertial parameters regressor of gravitation load for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:44:15
% EndTime: 2019-05-06 08:44:17
% DurationCPUTime: 0.46s
% Computational Cost: add. (221->99), mult. (546->126), div. (0->0), fcn. (578->8), ass. (0->68)
t41 = sin(qJ(2));
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t19 = g(1) * t45 + g(2) * t42;
t44 = cos(qJ(2));
t86 = t19 * t44;
t8 = g(3) * t41 + t86;
t85 = pkin(2) * t41;
t38 = sin(pkin(9));
t84 = pkin(4) * t38;
t83 = g(1) * t42;
t79 = g(3) * t44;
t33 = t44 * pkin(2);
t78 = t42 * t38;
t39 = cos(pkin(9));
t77 = t42 * t39;
t76 = t44 * t45;
t75 = t45 * t38;
t74 = t45 * t39;
t73 = pkin(2) + qJ(4);
t67 = qJ(3) * t44;
t23 = t42 * t67;
t65 = t44 * t84;
t72 = t42 * t65 + t23;
t26 = t45 * t67;
t71 = t45 * t65 + t26;
t29 = t41 * qJ(3);
t70 = t33 + t29;
t34 = t45 * pkin(7);
t69 = t45 * pkin(3) + t34;
t68 = t45 * pkin(1) + t42 * pkin(7);
t66 = qJ(5) * t39;
t30 = t44 * qJ(4);
t64 = -pkin(8) + t73;
t63 = t30 + t70;
t62 = t44 * t66;
t61 = -pkin(1) - t29;
t60 = pkin(2) * t76 + t45 * t29 + t68;
t11 = -t41 * t74 + t78;
t13 = t41 * t77 + t75;
t59 = g(1) * t13 + g(2) * t11;
t58 = -g(2) * t45 + t83;
t57 = pkin(5) * t38 - t66;
t14 = -t41 * t78 + t74;
t40 = sin(qJ(6));
t43 = cos(qJ(6));
t56 = t13 * t43 - t14 * t40;
t55 = t13 * t40 + t14 * t43;
t54 = t38 * t43 - t39 * t40;
t53 = t38 * t40 + t39 * t43;
t52 = g(3) * (t41 * t84 + t63);
t51 = t14 * pkin(4) + t13 * qJ(5) + t69;
t50 = t42 * pkin(3) + t45 * t30 + t60;
t49 = g(3) * t54;
t12 = t41 * t75 + t77;
t48 = t12 * pkin(4) + t11 * qJ(5) + t50;
t47 = t19 * t73;
t46 = (-t44 * t73 + t61) * t83;
t16 = -g(2) * t76 + t44 * t83;
t15 = t58 * t41;
t7 = t19 * t41 - t79;
t6 = t8 * t39;
t5 = t8 * t38;
t4 = -g(1) * t14 - g(2) * t12;
t3 = t11 * t40 + t12 * t43;
t2 = t11 * t43 - t12 * t40;
t1 = -g(1) * t11 + g(2) * t13 - t39 * t79;
t9 = [0, 0, 0, 0, 0, 0, t58, t19, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, -t19, -g(1) * (-t42 * pkin(1) + t34) - g(2) * t68, 0, 0, 0, 0, 0, 0, -t19, -t16, t15, -g(1) * t34 - g(2) * t60 - (t61 - t33) * t83, 0, 0, 0, 0, 0, 0, t4, t59, t16, -g(1) * t69 - g(2) * t50 - t46, 0, 0, 0, 0, 0, 0, t4, t16, -t59, -g(1) * t51 - g(2) * t48 - t46, 0, 0, 0, 0, 0, 0, -g(1) * t55 - g(2) * t3, -g(1) * t56 - g(2) * t2, -t16, -g(1) * (t14 * pkin(5) + t51) - g(2) * (t12 * pkin(5) - pkin(8) * t76 + t48) - (-t44 * t64 + t61) * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -g(1) * (-t45 * t85 + t26) - g(2) * (-t42 * t85 + t23) - g(3) * t70, 0, 0, 0, 0, 0, 0, -t5, -t6, t7, -g(1) * t26 - g(2) * t23 - g(3) * t63 + t41 * t47, 0, 0, 0, 0, 0, 0, -t5, t7, t6, -g(1) * (-t45 * t62 + t71) - g(2) * (-t42 * t62 + t72) - t52 + (g(3) * t66 + t47) * t41, 0, 0, 0, 0, 0, 0, -t41 * t49 - t54 * t86, t8 * t53, -t7, -g(1) * t71 - g(2) * t72 - t52 + (g(3) * pkin(8) - t19 * t57) * t44 + (-g(3) * t57 + t19 * t64) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 + g(2) * t56 - t53 * t79, g(1) * t3 - g(2) * t55 - t44 * t49, 0, 0;];
taug_reg  = t9;
