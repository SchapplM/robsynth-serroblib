% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:56:39
% EndTime: 2019-05-06 18:56:41
% DurationCPUTime: 0.51s
% Computational Cost: add. (315->108), mult. (528->137), div. (0->0), fcn. (529->8), ass. (0->70)
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t22 = g(1) * t45 + g(2) * t42;
t41 = sin(qJ(2));
t89 = t22 * t41;
t31 = t41 * qJ(3);
t44 = cos(qJ(2));
t57 = t44 * pkin(2) + t31;
t43 = cos(qJ(4));
t58 = t45 * t43;
t40 = sin(qJ(4));
t67 = t42 * t40;
t13 = t41 * t58 - t67;
t59 = t45 * t40;
t66 = t42 * t43;
t15 = t41 * t66 + t59;
t78 = g(3) * t44;
t88 = -g(1) * t13 - g(2) * t15 + t43 * t78;
t39 = qJ(4) + qJ(5);
t30 = cos(t39);
t60 = t45 * t30;
t29 = sin(t39);
t69 = t42 * t29;
t7 = t41 * t60 - t69;
t61 = t45 * t29;
t68 = t42 * t30;
t9 = t41 * t68 + t61;
t1 = -g(1) * t7 - g(2) * t9 + t30 * t78;
t12 = g(3) * t41 + t22 * t44;
t46 = -pkin(9) - pkin(8);
t83 = g(1) * t42;
t77 = t40 * pkin(4);
t76 = t44 * pkin(8);
t73 = t40 * t44;
t72 = t41 * t42;
t71 = t41 * t45;
t20 = pkin(5) * t29 + t77;
t70 = t42 * t20;
t38 = -qJ(6) + t46;
t65 = t44 * t38;
t64 = t44 * t45;
t63 = t44 * t46;
t62 = t45 * t20;
t33 = t43 * pkin(4);
t21 = pkin(5) * t30 + t33;
t56 = pkin(1) * t45 + pkin(7) * t42;
t55 = qJ(3) * t44;
t54 = pkin(4) * t73;
t51 = t41 * t59;
t50 = pkin(2) * t64 + t31 * t45 + t56;
t49 = -g(2) * t45 + t83;
t48 = -pkin(1) - t57;
t35 = t45 * pkin(7);
t28 = t33 + pkin(3);
t25 = t45 * t55;
t23 = t42 * t55;
t19 = pkin(3) + t21;
t18 = t49 * t44;
t17 = t49 * t41;
t16 = -t41 * t67 + t58;
t14 = t51 + t66;
t11 = -t78 + t89;
t10 = -t41 * t69 + t60;
t8 = t41 * t61 + t68;
t6 = t12 * t30;
t5 = t12 * t29;
t4 = -g(1) * t10 - g(2) * t8;
t3 = g(1) * t9 - g(2) * t7;
t2 = g(1) * t8 - g(2) * t10 - t29 * t78;
t24 = [0, 0, 0, 0, 0, 0, t49, t22, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, -t22, -g(1) * (-t42 * pkin(1) + t35) - g(2) * t56, 0, 0, 0, 0, 0, 0, -t22, -t18, t17, -g(1) * t35 - g(2) * t50 - t48 * t83, 0, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t14, g(1) * t15 - g(2) * t13, t18, -g(1) * (t45 * pkin(3) + t35) - g(2) * (pkin(8) * t64 + t50) + (-g(1) * (t48 - t76) - g(2) * pkin(3)) * t42, 0, 0, 0, 0, 0, 0, t4, t3, t18, -g(1) * (t45 * t28 + t35) - g(2) * (pkin(4) * t51 - t45 * t63 + t50) + (-g(1) * (-t41 * t77 + t48 + t63) - g(2) * t28) * t42, 0, 0, 0, 0, 0, 0, t4, t3, t18, -g(1) * (t45 * t19 + t35) - g(2) * (-t38 * t64 + t41 * t62 + t50) + (-g(1) * (-t41 * t20 + t48 + t65) - g(2) * t19) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, -g(1) * (-pkin(2) * t71 + t25) - g(2) * (-pkin(2) * t72 + t23) - g(3) * t57, 0, 0, 0, 0, 0, 0, -t12 * t40, -t12 * t43, t11, -g(1) * t25 - g(2) * t23 - g(3) * (t57 + t76) + (pkin(2) + pkin(8)) * t89, 0, 0, 0, 0, 0, 0, -t5, -t6, t11, -g(1) * (t45 * t54 + t25) - g(2) * (t42 * t54 + t23) - g(3) * (t57 - t63) + (-g(3) * t77 + t22 * (pkin(2) - t46)) * t41, 0, 0, 0, 0, 0, 0, -t5, -t6, t11, -g(1) * (t44 * t62 + t25) - g(2) * (t44 * t70 + t23) - g(3) * (t57 - t65) + (-g(3) * t20 + t22 * (pkin(2) - t38)) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, g(1) * t14 - g(2) * t16 - g(3) * t73, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t88 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t21 * t71 - t70) - g(2) * (t21 * t72 + t62) + t21 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12;];
taug_reg  = t24;
