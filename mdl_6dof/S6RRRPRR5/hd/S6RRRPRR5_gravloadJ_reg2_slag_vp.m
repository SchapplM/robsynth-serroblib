% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:43:18
% EndTime: 2019-05-07 10:43:21
% DurationCPUTime: 0.67s
% Computational Cost: add. (447->109), mult. (517->131), div. (0->0), fcn. (503->10), ass. (0->75)
t42 = qJ(2) + qJ(3);
t37 = sin(t42);
t49 = -pkin(10) - pkin(9);
t39 = cos(t42);
t43 = sin(qJ(5));
t80 = t39 * t43;
t96 = pkin(5) * t80 + t37 * t49;
t69 = t39 * pkin(3) + t37 * qJ(4);
t50 = -pkin(8) - pkin(7);
t82 = pkin(4) - t50;
t45 = sin(qJ(1));
t48 = cos(qJ(1));
t23 = g(1) * t48 + g(2) * t45;
t95 = t23 * t37;
t94 = t96 * t45;
t93 = t96 * t48;
t57 = t37 * t43 * pkin(5) - t39 * t49;
t46 = cos(qJ(5));
t71 = t48 * t46;
t76 = t45 * t43;
t15 = t37 * t71 - t76;
t72 = t48 * t43;
t75 = t45 * t46;
t17 = t37 * t75 + t72;
t83 = g(3) * t39;
t92 = -g(1) * t15 - g(2) * t17 + t46 * t83;
t8 = g(3) * t37 + t23 * t39;
t44 = sin(qJ(2));
t90 = pkin(2) * t44;
t89 = pkin(3) * t37;
t31 = t39 * pkin(9);
t41 = qJ(5) + qJ(6);
t36 = sin(t41);
t78 = t45 * t36;
t38 = cos(t41);
t77 = t45 * t38;
t74 = t48 * t36;
t73 = t48 * t38;
t70 = t48 * t50;
t68 = t46 * pkin(5) + t82;
t67 = qJ(4) * t39;
t47 = cos(qJ(2));
t40 = t47 * pkin(2);
t64 = t40 + t69;
t35 = t40 + pkin(1);
t27 = t48 * t35;
t63 = g(2) * (t69 * t48 + t27);
t24 = t45 * t67;
t62 = -t45 * t89 + t24;
t26 = t48 * t67;
t61 = -t48 * t89 + t26;
t60 = -t89 - t90;
t59 = g(1) * t45 - g(2) * t48;
t58 = t69 + t57;
t56 = -t35 - t69;
t55 = t60 * t45 + t24;
t54 = t60 * t48 + t26;
t52 = -g(3) * t47 + t23 * t44;
t51 = (pkin(3) + pkin(9)) * t95;
t18 = -t37 * t76 + t71;
t16 = t37 * t72 + t75;
t14 = t59 * t39;
t13 = t59 * t37;
t12 = -t37 * t78 + t73;
t11 = t37 * t77 + t74;
t10 = t37 * t74 + t77;
t9 = t37 * t73 - t78;
t7 = -t83 + t95;
t6 = t8 * t46;
t5 = t8 * t43;
t4 = t8 * t38;
t3 = t8 * t36;
t2 = g(1) * t10 - g(2) * t12 - t36 * t83;
t1 = -g(1) * t9 - g(2) * t11 + t38 * t83;
t19 = [0, 0, 0, 0, 0, 0, t59, t23, 0, 0, 0, 0, 0, 0, 0, 0, t59 * t47, -t59 * t44, -t23, -g(1) * (-t45 * pkin(1) + t48 * pkin(7)) - g(2) * (t48 * pkin(1) + t45 * pkin(7)) 0, 0, 0, 0, 0, 0, t14, -t13, -t23, -g(1) * (-t45 * t35 - t70) - g(2) * (-t45 * t50 + t27) 0, 0, 0, 0, 0, 0, -t23, -t14, t13, g(1) * t70 - t63 + (-g(1) * t56 + g(2) * t50) * t45, 0, 0, 0, 0, 0, 0, -g(1) * t18 - g(2) * t16, g(1) * t17 - g(2) * t15, t14, -t63 + (-g(1) * t82 - g(2) * t31) * t48 + (-g(1) * (t56 - t31) - g(2) * t82) * t45, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t10, g(1) * t11 - g(2) * t9, t14, -t63 + (-g(1) * t68 - g(2) * t57) * t48 + (-g(1) * (t56 - t57) - g(2) * t68) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, g(3) * t44 + t23 * t47, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t52 * pkin(2), 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -g(1) * t54 - g(2) * t55 - g(3) * t64, 0, 0, 0, 0, 0, 0, -t5, -t6, t7, -g(1) * (-t48 * t90 + t26) - g(2) * (-t45 * t90 + t24) - g(3) * (t31 + t64) + t51, 0, 0, 0, 0, 0, 0, -t3, -t4, t7, -g(1) * (t54 + t93) - g(2) * (t55 + t94) - g(3) * (t40 + t58); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, -g(1) * t61 - g(2) * t62 - g(3) * t69, 0, 0, 0, 0, 0, 0, -t5, -t6, t7, -g(1) * t26 - g(2) * t24 - g(3) * (t31 + t69) + t51, 0, 0, 0, 0, 0, 0, -t3, -t4, t7, -g(1) * (t61 + t93) - g(2) * (t62 + t94) - g(3) * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, g(1) * t16 - g(2) * t18 - g(3) * t80, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t92 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t19;
