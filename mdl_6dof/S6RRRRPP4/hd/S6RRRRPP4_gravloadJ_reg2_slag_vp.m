% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:19:22
% EndTime: 2019-05-07 18:19:24
% DurationCPUTime: 0.59s
% Computational Cost: add. (582->119), mult. (660->164), div. (0->0), fcn. (681->10), ass. (0->80)
t55 = qJ(3) + qJ(4);
t48 = cos(t55);
t59 = cos(qJ(3));
t50 = t59 * pkin(3);
t35 = pkin(4) * t48 + t50;
t33 = pkin(2) + t35;
t60 = cos(qJ(2));
t28 = t60 * t33;
t62 = -pkin(9) - pkin(8);
t54 = -qJ(5) + t62;
t57 = sin(qJ(2));
t92 = t57 * t54;
t105 = t28 - t92;
t58 = sin(qJ(1));
t61 = cos(qJ(1));
t40 = g(1) * t61 + g(2) * t58;
t56 = sin(qJ(3));
t82 = t61 * t59;
t89 = t58 * t60;
t24 = t56 * t89 + t82;
t83 = t61 * t56;
t26 = t58 * t59 - t60 * t83;
t96 = g(3) * t57;
t104 = -g(1) * t26 + g(2) * t24 + t56 * t96;
t19 = -g(3) * t60 + t40 * t57;
t47 = sin(t55);
t103 = pkin(4) * t47;
t46 = pkin(10) + t55;
t42 = sin(t46);
t102 = pkin(5) * t42;
t100 = g(1) * t58;
t94 = t56 * pkin(3);
t43 = cos(t46);
t93 = t43 * t57;
t91 = t57 * t62;
t90 = t58 * t48;
t88 = t60 * t61;
t34 = t94 + t103;
t32 = t61 * t34;
t87 = t61 * t42;
t86 = t61 * t43;
t85 = t61 * t47;
t84 = t61 * t48;
t81 = -t60 * t32 + t58 * t35;
t80 = t61 * pkin(1) + t58 * pkin(7);
t78 = t60 * t85;
t11 = t42 * t89 + t86;
t12 = t43 * t89 - t87;
t77 = -t11 * pkin(5) + t12 * qJ(6);
t76 = -t34 * t89 - t61 * t35;
t13 = -t58 * t43 + t60 * t87;
t14 = t58 * t42 + t60 * t86;
t75 = -t13 * pkin(5) + t14 * qJ(6);
t74 = t60 * pkin(2) + t57 * pkin(8);
t72 = g(1) * t11 - g(2) * t13;
t71 = -g(2) * t61 + t100;
t70 = pkin(5) * t43 + qJ(6) * t42;
t45 = t50 + pkin(2);
t68 = t60 * t45 - t91;
t15 = t47 * t89 + t84;
t65 = t33 * t88 + t58 * t34 - t61 * t92 + t80;
t1 = g(1) * t13 + g(2) * t11 + t42 * t96;
t3 = g(1) * t14 + g(2) * t12 + g(3) * t93;
t51 = t61 * pkin(7);
t64 = t32 + t51 + (-pkin(1) - t105) * t58;
t20 = t40 * t60 + t96;
t41 = pkin(4) * t90;
t36 = qJ(6) * t93;
t29 = t71 * t57;
t27 = t58 * t56 + t60 * t82;
t25 = -t59 * t89 + t83;
t18 = t58 * t47 + t60 * t84;
t17 = -t78 + t90;
t16 = -t48 * t89 + t85;
t8 = t19 * t43;
t7 = t19 * t42;
t6 = g(1) * t12 - g(2) * t14;
t5 = g(1) * t18 - g(2) * t16 + t48 * t96;
t4 = -g(1) * t17 + g(2) * t15 + t47 * t96;
t2 = [0, 0, 0, 0, 0, 0, t71, t40, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t60, -t29, -t40, -g(1) * (-t58 * pkin(1) + t51) - g(2) * t80, 0, 0, 0, 0, 0, 0, -g(1) * t25 - g(2) * t27, -g(1) * t24 - g(2) * t26, t29, -g(1) * t51 - g(2) * (t74 * t61 + t80) - (-pkin(1) - t74) * t100, 0, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, -g(1) * t15 - g(2) * t17, t29, -g(1) * (pkin(3) * t83 + t51) - g(2) * (t45 * t88 - t61 * t91 + t80) + (-g(1) * (-pkin(1) - t68) - g(2) * t94) * t58, 0, 0, 0, 0, 0, 0, t6, -t72, t29, -g(1) * t64 - g(2) * t65, 0, 0, 0, 0, 0, 0, t6, t29, t72, -g(1) * (-t12 * pkin(5) - t11 * qJ(6) + t64) - g(2) * (t14 * pkin(5) + t13 * qJ(6) + t65); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, 0, 0, 0, t19 * t59, -t19 * t56, -t20, -g(3) * t74 + t40 * (pkin(2) * t57 - pkin(8) * t60) 0, 0, 0, 0, 0, 0, t19 * t48, -t19 * t47, -t20, -g(3) * t68 + t40 * (t45 * t57 + t60 * t62) 0, 0, 0, 0, 0, 0, t8, -t7, -t20, -g(3) * t105 + t40 * (t33 * t57 + t54 * t60) 0, 0, 0, 0, 0, 0, t8, -t20, t7, -g(3) * t28 + (-g(3) * t70 + t40 * t54) * t60 + (g(3) * t54 + t40 * (t33 + t70)) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, g(1) * t27 - g(2) * t25 + t59 * t96, 0, 0, 0, 0, 0, 0, 0, 0, t4, t5, 0, t104 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t3, 0, -g(1) * t81 - g(2) * t76 + t34 * t96, 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * (t75 + t81) - g(2) * (t76 + t77) - g(3) * (t36 + (-t34 - t102) * t57); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t5, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, -g(1) * t41 + (g(2) * t84 + t20 * t47) * pkin(4), 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * (-pkin(4) * t78 + t41 + t75) - g(2) * (-t15 * pkin(4) + t77) - g(3) * (t36 + (-t102 - t103) * t57); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
