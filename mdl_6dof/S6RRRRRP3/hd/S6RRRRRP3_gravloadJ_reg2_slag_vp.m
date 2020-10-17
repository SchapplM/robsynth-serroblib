% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:42:42
% EndTime: 2019-05-08 04:42:44
% DurationCPUTime: 0.55s
% Computational Cost: add. (552->103), mult. (592->136), div. (0->0), fcn. (584->10), ass. (0->72)
t41 = qJ(4) + qJ(5);
t36 = cos(t41);
t46 = cos(qJ(4));
t38 = t46 * pkin(4);
t23 = pkin(5) * t36 + t38;
t21 = pkin(3) + t23;
t42 = qJ(2) + qJ(3);
t35 = sin(t42);
t37 = cos(t42);
t49 = -pkin(10) - pkin(9);
t40 = -qJ(6) + t49;
t96 = t37 * t21 - t35 * t40;
t32 = t38 + pkin(3);
t95 = t37 * t32 - t35 * t49;
t94 = t37 * pkin(3) + t35 * pkin(9);
t45 = sin(qJ(1));
t48 = cos(qJ(1));
t25 = g(1) * t48 + g(2) * t45;
t69 = t48 * t46;
t43 = sin(qJ(4));
t74 = t45 * t43;
t16 = t37 * t74 + t69;
t70 = t48 * t43;
t73 = t45 * t46;
t18 = -t37 * t70 + t73;
t83 = g(3) * t35;
t93 = -g(1) * t18 + g(2) * t16 + t43 * t83;
t71 = t48 * t36;
t34 = sin(t41);
t76 = t45 * t34;
t11 = t37 * t76 + t71;
t72 = t48 * t34;
t75 = t45 * t36;
t13 = -t37 * t72 + t75;
t1 = -g(1) * t13 + g(2) * t11 + t34 * t83;
t9 = -g(3) * t37 + t25 * t35;
t44 = sin(qJ(2));
t92 = pkin(2) * t44;
t91 = pkin(3) * t35;
t47 = cos(qJ(2));
t39 = t47 * pkin(2);
t33 = t39 + pkin(1);
t26 = t48 * t33;
t85 = g(2) * t26;
t81 = t43 * pkin(4);
t78 = t37 * t45;
t77 = t37 * t48;
t22 = pkin(5) * t34 + t81;
t50 = -pkin(8) - pkin(7);
t68 = t22 - t50;
t64 = -t50 + t81;
t61 = -t91 - t92;
t59 = g(1) * t45 - g(2) * t48;
t57 = t21 * t35 + t37 * t40;
t55 = t32 * t35 + t37 * t49;
t51 = -g(3) * t47 + t25 * t44;
t28 = pkin(9) * t77;
t27 = pkin(9) * t78;
t19 = t37 * t69 + t74;
t17 = -t37 * t73 + t70;
t15 = t59 * t35;
t14 = t37 * t71 + t76;
t12 = -t37 * t75 + t72;
t10 = t25 * t37 + t83;
t8 = t9 * t46;
t7 = t9 * t43;
t6 = t9 * t36;
t5 = t9 * t34;
t4 = -g(1) * t12 - g(2) * t14;
t3 = -g(1) * t11 - g(2) * t13;
t2 = g(1) * t14 - g(2) * t12 + t36 * t83;
t20 = [0, 0, 0, 0, 0, 0, t59, t25, 0, 0, 0, 0, 0, 0, 0, 0, t59 * t47, -t59 * t44, -t25, -g(1) * (-t45 * pkin(1) + t48 * pkin(7)) - g(2) * (t48 * pkin(1) + t45 * pkin(7)) 0, 0, 0, 0, 0, 0, t59 * t37, -t15, -t25, -g(1) * (-t45 * t33 - t48 * t50) - g(2) * (-t45 * t50 + t26) 0, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t19, -g(1) * t16 - g(2) * t18, t15, -t85 + (g(1) * t50 - g(2) * t94) * t48 + (-g(1) * (-t33 - t94) + g(2) * t50) * t45, 0, 0, 0, 0, 0, 0, t4, t3, t15, -t85 + (-g(1) * t64 - g(2) * t95) * t48 + (-g(1) * (-t33 - t95) - g(2) * t64) * t45, 0, 0, 0, 0, 0, 0, t4, t3, t15, -t85 + (-g(1) * t68 - g(2) * t96) * t48 + (-g(1) * (-t33 - t96) - g(2) * t68) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, g(3) * t44 + t25 * t47, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t51 * pkin(2), 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * (t61 * t48 + t28) - g(2) * (t61 * t45 + t27) - g(3) * (t39 + t94) 0, 0, 0, 0, 0, 0, t6, -t5, -t10, -g(3) * (t39 + t95) + t25 * (t55 + t92) 0, 0, 0, 0, 0, 0, t6, -t5, -t10, -g(3) * (t39 + t96) + t25 * (t57 + t92); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * (-t48 * t91 + t28) - g(2) * (-t45 * t91 + t27) - g(3) * t94, 0, 0, 0, 0, 0, 0, t6, -t5, -t10, -g(3) * t95 + t25 * t55, 0, 0, 0, 0, 0, 0, t6, -t5, -t10, -g(3) * t96 + t25 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, g(1) * t19 - g(2) * t17 + t46 * t83, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t93 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t22 * t77 + t45 * t23) - g(2) * (-t22 * t78 - t48 * t23) + t22 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9;];
taug_reg  = t20;
