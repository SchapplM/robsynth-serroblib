% Calculate inertial parameters regressor of gravitation load for
% S6PRRPRP2
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
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:47:32
% EndTime: 2019-05-05 03:47:34
% DurationCPUTime: 0.58s
% Computational Cost: add. (580->118), mult. (1055->172), div. (0->0), fcn. (1279->12), ass. (0->77)
t61 = sin(pkin(6));
t65 = sin(qJ(2));
t100 = t61 * t65;
t64 = sin(qJ(3));
t67 = cos(qJ(3));
t94 = cos(pkin(6));
t109 = -t64 * t100 + t94 * t67;
t68 = cos(qJ(2));
t60 = sin(pkin(10));
t89 = t60 * t94;
t93 = cos(pkin(10));
t45 = -t65 * t89 + t68 * t93;
t99 = t61 * t67;
t108 = -t45 * t64 + t60 * t99;
t59 = qJ(3) + pkin(11);
t57 = sin(t59);
t58 = cos(t59);
t107 = pkin(4) * t58 + pkin(9) * t57;
t63 = sin(qJ(5));
t103 = t58 * t63;
t66 = cos(qJ(5));
t102 = t58 * t66;
t101 = t60 * t61;
t98 = t61 * t68;
t97 = t66 * t68;
t80 = t94 * t93;
t42 = t60 * t65 - t68 * t80;
t43 = t60 * t68 + t65 * t80;
t56 = t67 * pkin(3) + pkin(2);
t62 = -qJ(4) - pkin(8);
t96 = -t42 * t56 - t43 * t62;
t44 = t65 * t93 + t68 * t89;
t95 = -t44 * t56 - t45 * t62;
t90 = t63 * t98;
t88 = t61 * t93;
t86 = -t107 * t42 + t96;
t85 = -t107 * t44 + t95;
t84 = t108 * pkin(3);
t83 = -t62 * t100 + t56 * t98;
t82 = pkin(5) * t66 + qJ(6) * t63;
t81 = t109 * pkin(3);
t79 = t107 * t98 + t83;
t35 = t58 * t100 + t57 * t94;
t23 = t35 * t63 + t61 * t97;
t20 = t43 * t58 - t57 * t88;
t7 = t20 * t63 - t42 * t66;
t22 = t57 * t101 + t45 * t58;
t9 = t22 * t63 - t44 * t66;
t1 = g(1) * t9 + g(2) * t7 + g(3) * t23;
t10 = t22 * t66 + t44 * t63;
t24 = t35 * t66 - t90;
t8 = t20 * t66 + t42 * t63;
t78 = g(1) * t10 + g(2) * t8 + g(3) * t24;
t11 = -t42 * t103 - t43 * t66;
t13 = -t44 * t103 - t45 * t66;
t25 = -t66 * t100 + t58 * t90;
t77 = g(1) * t13 + g(2) * t11 + g(3) * t25;
t19 = -t43 * t57 - t58 * t88;
t21 = t58 * t101 - t45 * t57;
t34 = -t57 * t100 + t58 * t94;
t76 = g(1) * t21 + g(2) * t19 + g(3) * t34;
t75 = g(1) * t22 + g(2) * t20 + g(3) * t35;
t74 = t21 * pkin(4) + t22 * pkin(9) + t84;
t73 = -t43 * t64 - t67 * t88;
t15 = -g(1) * t44 - g(2) * t42 + g(3) * t98;
t72 = g(1) * t45 + g(2) * t43 + g(3) * t100;
t71 = t34 * pkin(4) + t35 * pkin(9) + t81;
t70 = t73 * pkin(3);
t69 = t19 * pkin(4) + t20 * pkin(9) + t70;
t26 = (t58 * t97 + t63 * t65) * t61;
t14 = -t44 * t102 + t45 * t63;
t12 = -t42 * t102 + t43 * t63;
t6 = t15 * t57;
t4 = t76 * t66;
t3 = t76 * t63;
t2 = -g(1) * t14 - g(2) * t12 - g(3) * t26;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, t72, 0, 0, 0, 0, 0, 0, 0, 0, -t15 * t67, t15 * t64, -t72, -g(1) * (-t44 * pkin(2) + t45 * pkin(8)) - g(2) * (-t42 * pkin(2) + t43 * pkin(8)) - g(3) * (pkin(2) * t68 + pkin(8) * t65) * t61, 0, 0, 0, 0, 0, 0, -t15 * t58, t6, -t72, -g(1) * t95 - g(2) * t96 - g(3) * t83, 0, 0, 0, 0, 0, 0, t2, t77, -t6, -g(1) * t85 - g(2) * t86 - g(3) * t79, 0, 0, 0, 0, 0, 0, t2, -t6, -t77, -g(1) * (t14 * pkin(5) + t13 * qJ(6) + t85) - g(2) * (t12 * pkin(5) + t11 * qJ(6) + t86) - g(3) * (t26 * pkin(5) + t25 * qJ(6) + t79); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t108 - g(2) * t73 - g(3) * t109, -g(1) * (-t64 * t101 - t45 * t67) - g(2) * (-t43 * t67 + t64 * t88) - g(3) * (-t94 * t64 - t65 * t99) 0, 0, 0, 0, 0, 0, 0, 0, -t76, t75, 0, -g(1) * t84 - g(2) * t70 - g(3) * t81, 0, 0, 0, 0, 0, 0, -t4, t3, -t75, -g(1) * t74 - g(2) * t69 - g(3) * t71, 0, 0, 0, 0, 0, 0, -t4, -t75, -t3, -g(1) * (t82 * t21 + t74) - g(2) * (t82 * t19 + t69) - g(3) * (t82 * t34 + t71); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t78, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t78, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - g(3) * (-t23 * pkin(5) + t24 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t5;
