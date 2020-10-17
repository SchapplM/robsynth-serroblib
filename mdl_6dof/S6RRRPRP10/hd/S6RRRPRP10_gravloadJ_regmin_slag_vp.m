% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:59:58
% EndTime: 2019-05-07 09:00:00
% DurationCPUTime: 0.71s
% Computational Cost: add. (644->141), mult. (1355->215), div. (0->0), fcn. (1692->12), ass. (0->73)
t66 = sin(qJ(2));
t68 = cos(qJ(2));
t100 = cos(qJ(1));
t87 = cos(pkin(6));
t79 = t87 * t100;
t99 = sin(qJ(1));
t42 = t66 * t79 + t99 * t68;
t65 = sin(qJ(3));
t67 = cos(qJ(3));
t62 = sin(pkin(6));
t83 = t62 * t100;
t23 = t42 * t67 - t65 * t83;
t41 = t99 * t66 - t68 * t79;
t60 = pkin(11) + qJ(5);
t57 = sin(t60);
t58 = cos(t60);
t9 = t23 * t57 - t41 * t58;
t10 = t23 * t58 + t41 * t57;
t22 = t42 * t65 + t67 * t83;
t78 = t87 * t99;
t44 = t100 * t68 - t66 * t78;
t82 = t62 * t99;
t26 = t44 * t65 - t67 * t82;
t92 = t62 * t66;
t39 = t65 * t92 - t87 * t67;
t71 = g(1) * t26 + g(2) * t22 + g(3) * t39;
t101 = g(3) * t62;
t96 = t57 * t67;
t95 = t58 * t67;
t61 = sin(pkin(11));
t94 = t61 * t66;
t93 = t61 * t67;
t91 = t62 * t68;
t63 = cos(pkin(11));
t90 = t63 * t67;
t89 = t67 * t68;
t88 = pkin(2) * t91 + pkin(9) * t92;
t86 = t57 * t91;
t85 = t100 * pkin(1) + t44 * pkin(2) + pkin(8) * t82;
t84 = pkin(4) * t61 + pkin(9);
t27 = t44 * t67 + t65 * t82;
t43 = t100 * t66 + t68 * t78;
t13 = t27 * t57 - t43 * t58;
t81 = -g(1) * t9 + g(2) * t13;
t80 = -g(1) * t22 + g(2) * t26;
t77 = pkin(3) * t67 + qJ(4) * t65;
t56 = t63 * pkin(4) + pkin(3);
t64 = -pkin(10) - qJ(4);
t76 = t56 * t67 - t64 * t65;
t75 = -t99 * pkin(1) - t42 * pkin(2) + pkin(8) * t83;
t40 = t87 * t65 + t67 * t92;
t20 = t40 * t57 + t58 * t91;
t1 = g(1) * t13 + g(2) * t9 + g(3) * t20;
t14 = t27 * t58 + t43 * t57;
t21 = t40 * t58 - t86;
t73 = g(1) * t14 + g(2) * t10 + g(3) * t21;
t16 = -t41 * t96 - t42 * t58;
t18 = -t43 * t96 - t44 * t58;
t28 = -t58 * t92 + t67 * t86;
t72 = g(1) * t18 + g(2) * t16 + g(3) * t28;
t70 = g(1) * t27 + g(2) * t23 + g(3) * t40;
t69 = -g(1) * t43 - g(2) * t41 + g(3) * t91;
t37 = t43 * pkin(2);
t35 = t41 * pkin(2);
t29 = (t57 * t66 + t58 * t89) * t62;
t19 = -t43 * t95 + t44 * t57;
t17 = -t41 * t95 + t42 * t57;
t15 = t69 * t65;
t5 = t71 * t58;
t4 = t71 * t57;
t3 = g(1) * t10 - g(2) * t14;
t2 = -g(1) * t19 - g(2) * t17 - g(3) * t29;
t6 = [0, g(1) * t99 - g(2) * t100, g(1) * t100 + g(2) * t99, 0, 0, 0, 0, 0, g(1) * t42 - g(2) * t44, -g(1) * t41 + g(2) * t43, 0, 0, 0, 0, 0, g(1) * t23 - g(2) * t27, t80, -g(1) * (-t23 * t63 - t41 * t61) - g(2) * (t27 * t63 + t43 * t61) -g(1) * (t23 * t61 - t41 * t63) - g(2) * (-t27 * t61 + t43 * t63) -t80, -g(1) * (-pkin(3) * t23 - t41 * pkin(9) - qJ(4) * t22 + t75) - g(2) * (t27 * pkin(3) + t43 * pkin(9) + t26 * qJ(4) + t85) 0, 0, 0, 0, 0, t3, t81, t3, -t80, -t81, -g(1) * (-pkin(5) * t10 - qJ(6) * t9 + t22 * t64 - t23 * t56 - t84 * t41 + t75) - g(2) * (t14 * pkin(5) + t13 * qJ(6) - t26 * t64 + t27 * t56 + t84 * t43 + t85); 0, 0, 0, 0, 0, 0, 0, 0, -t69, g(1) * t44 + g(2) * t42 + g(3) * t92, 0, 0, 0, 0, 0, -t69 * t67, t15, -g(1) * (-t43 * t90 + t44 * t61) - g(2) * (-t41 * t90 + t42 * t61) - (t63 * t89 + t94) * t101, -g(1) * (t43 * t93 + t44 * t63) - g(2) * (t41 * t93 + t42 * t63) - (-t61 * t89 + t63 * t66) * t101, -t15, -g(1) * (t44 * pkin(9) - t77 * t43 - t37) - g(2) * (t42 * pkin(9) - t77 * t41 - t35) - g(3) * (t77 * t91 + t88) 0, 0, 0, 0, 0, t2, t72, t2, -t15, -t72, -g(1) * (t19 * pkin(5) + t18 * qJ(6) - t76 * t43 + t84 * t44 - t37) - g(2) * (t17 * pkin(5) + t16 * qJ(6) - t76 * t41 + t84 * t42 - t35) - g(3) * (t29 * pkin(5) + t28 * qJ(6) + t88) - (pkin(4) * t94 + t76 * t68) * t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t70, t71 * t63, -t71 * t61, -t70, -g(1) * (-t26 * pkin(3) + t27 * qJ(4)) - g(2) * (-t22 * pkin(3) + t23 * qJ(4)) - g(3) * (-t39 * pkin(3) + t40 * qJ(4)) 0, 0, 0, 0, 0, t5, -t4, t5, -t70, t4, t70 * t64 + t71 * (pkin(5) * t58 + qJ(6) * t57 + t56); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t73, t1, 0, -t73, -g(1) * (-t13 * pkin(5) + t14 * qJ(6)) - g(2) * (-t9 * pkin(5) + t10 * qJ(6)) - g(3) * (-t20 * pkin(5) + t21 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t6;
