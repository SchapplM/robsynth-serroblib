% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP14_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP14_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t59 = sin(qJ(2));
t60 = sin(qJ(1));
t63 = cos(qJ(2));
t64 = cos(qJ(1));
t85 = cos(pkin(6));
t82 = t64 * t85;
t43 = t59 * t82 + t60 * t63;
t83 = t60 * t85;
t45 = -t59 * t83 + t64 * t63;
t101 = -g(1) * t45 - g(2) * t43;
t44 = t64 * t59 + t63 * t83;
t58 = sin(qJ(4));
t62 = cos(qJ(4));
t56 = sin(pkin(6));
t92 = t56 * t60;
t20 = -t44 * t62 + t58 * t92;
t42 = t60 * t59 - t63 * t82;
t90 = t56 * t64;
t71 = t42 * t62 + t58 * t90;
t91 = t56 * t63;
t67 = -g(3) * (-t58 * t85 - t62 * t91) - g(2) * t71 + g(1) * t20;
t57 = sin(qJ(5));
t61 = cos(qJ(5));
t72 = -t42 * t58 + t62 * t90;
t8 = t43 * t61 + t57 * t72;
t9 = -t43 * t57 + t61 * t72;
t93 = t56 * t59;
t89 = t57 * t58;
t88 = t58 * t61;
t87 = t59 * t61;
t86 = pkin(2) * t91 + qJ(3) * t93;
t84 = t57 * t93;
t21 = t44 * t58 + t62 * t92;
t6 = t21 * t57 - t45 * t61;
t81 = g(1) * t8 + g(2) * t6;
t80 = pkin(4) * t58 - pkin(10) * t62;
t79 = g(1) * t71 + g(2) * t20;
t78 = g(1) * t42 - g(2) * t44;
t77 = g(1) * t43 - g(2) * t45;
t76 = g(1) * t64 + g(2) * t60;
t75 = t64 * pkin(1) + t45 * pkin(2) + pkin(8) * t92 + t44 * qJ(3);
t41 = -t58 * t91 + t62 * t85;
t18 = t41 * t57 - t56 * t87;
t1 = g(1) * t6 - g(2) * t8 + g(3) * t18;
t19 = t41 * t61 + t84;
t7 = t21 * t61 + t45 * t57;
t70 = g(1) * t7 - g(2) * t9 + g(3) * t19;
t69 = -t60 * pkin(1) - t43 * pkin(2) + pkin(8) * t90 - t42 * qJ(3);
t14 = t42 * t61 + t43 * t89;
t16 = t44 * t61 + t45 * t89;
t26 = t58 * t84 - t61 * t91;
t68 = g(1) * t16 + g(2) * t14 + g(3) * t26;
t66 = g(1) * t21 - g(2) * t72 + g(3) * t41;
t13 = -g(1) * t44 - g(2) * t42 + g(3) * t91;
t65 = g(3) * t93 - t101;
t38 = t44 * pkin(2);
t36 = t42 * pkin(2);
t27 = (t57 * t63 + t58 * t87) * t56;
t17 = -t44 * t57 + t45 * t88;
t15 = -t42 * t57 + t43 * t88;
t12 = t65 * t62;
t5 = t67 * t61;
t4 = t67 * t57;
t3 = -g(1) * t9 - g(2) * t7;
t2 = -g(1) * t17 - g(2) * t15 - g(3) * t27;
t10 = [0, g(1) * t60 - g(2) * t64, t76, 0, 0, 0, 0, 0, t77, -t78, -t76 * t56, -t77, t78, -g(1) * t69 - g(2) * t75, 0, 0, 0, 0, 0, -g(1) * t72 - g(2) * t21, t79, 0, 0, 0, 0, 0, t3, t81, t3, -t79, -t81, -g(1) * (pkin(3) * t90 + pkin(4) * t72 + t9 * pkin(5) - t43 * pkin(9) + pkin(10) * t71 + t8 * qJ(6) + t69) - g(2) * (pkin(3) * t92 + t21 * pkin(4) + t7 * pkin(5) + t45 * pkin(9) + t20 * pkin(10) + t6 * qJ(6) + t75); 0, 0, 0, 0, 0, 0, 0, 0, -t13, t65, 0, t13, -t65, -g(1) * (t45 * qJ(3) - t38) - g(2) * (t43 * qJ(3) - t36) - g(3) * t86, 0, 0, 0, 0, 0, -t65 * t58, -t12, 0, 0, 0, 0, 0, t2, t68, t2, t12, -t68, -g(1) * (t17 * pkin(5) - t44 * pkin(9) + t16 * qJ(6) - t38) - g(2) * (t15 * pkin(5) - t42 * pkin(9) + t14 * qJ(6) - t36) + t101 * (qJ(3) + t80) + (-t27 * pkin(5) - t26 * qJ(6) - t86 - (pkin(9) * t63 + t80 * t59) * t56) * g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t66, 0, 0, 0, 0, 0, t5, -t4, t5, -t66, t4, -pkin(10) * t66 + t67 * (pkin(5) * t61 + qJ(6) * t57 + pkin(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t70, t1, 0, -t70, -g(1) * (-t6 * pkin(5) + t7 * qJ(6)) - g(2) * (pkin(5) * t8 - qJ(6) * t9) - g(3) * (-t18 * pkin(5) + t19 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t10;
