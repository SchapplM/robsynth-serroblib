% Calculate minimal parameter regressor of gravitation load for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRPP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t72 = cos(qJ(3));
t110 = -pkin(3) * t72 - pkin(2);
t68 = sin(qJ(4));
t71 = cos(qJ(4));
t109 = pkin(4) * t71 + qJ(5) * t68 + pkin(3);
t108 = pkin(5) + pkin(9);
t66 = sin(pkin(10));
t70 = sin(qJ(2));
t73 = cos(qJ(2));
t94 = cos(pkin(10));
t95 = cos(pkin(6));
t83 = t95 * t94;
t51 = t66 * t70 - t73 * t83;
t69 = sin(qJ(3));
t105 = t51 * t69;
t87 = t66 * t95;
t53 = t94 * t70 + t73 * t87;
t104 = t53 * t69;
t67 = sin(pkin(6));
t103 = t67 * t70;
t102 = t67 * t72;
t101 = t67 * t73;
t100 = t68 * t72;
t99 = t71 * t72;
t98 = t72 * t73;
t96 = qJ(6) * t71;
t93 = t69 * t101;
t92 = t68 * t101;
t52 = t66 * t73 + t70 * t83;
t86 = t67 * t94;
t29 = -t52 * t69 - t72 * t86;
t91 = t109 * t29;
t54 = -t70 * t87 + t94 * t73;
t31 = t66 * t102 - t54 * t69;
t90 = t109 * t31;
t55 = -t69 * t103 + t95 * t72;
t89 = t109 * t55;
t30 = t52 * t72 - t69 * t86;
t11 = t30 * t68 - t51 * t71;
t12 = t30 * t71 + t51 * t68;
t88 = -t11 * pkin(4) + t12 * qJ(5);
t32 = t66 * t67 * t69 + t54 * t72;
t13 = t32 * t68 - t53 * t71;
t14 = t32 * t71 + t53 * t68;
t85 = -t13 * pkin(4) + t14 * qJ(5);
t56 = t70 * t102 + t95 * t69;
t33 = t71 * t101 + t56 * t68;
t34 = t56 * t71 - t92;
t84 = -t33 * pkin(4) + t34 * qJ(5);
t2 = g(1) * t13 + g(2) * t11 + g(3) * t33;
t82 = g(1) * t14 + g(2) * t12 + g(3) * t34;
t18 = -t51 * t100 - t52 * t71;
t20 = -t53 * t100 - t54 * t71;
t36 = -t71 * t103 + t72 * t92;
t81 = g(1) * t20 + g(2) * t18 + g(3) * t36;
t19 = -t51 * t99 + t52 * t68;
t21 = -t53 * t99 + t54 * t68;
t37 = (t68 * t70 + t71 * t98) * t67;
t80 = g(1) * t21 + g(2) * t19 + g(3) * t37;
t79 = g(1) * t31 + g(2) * t29 + g(3) * t55;
t78 = g(1) * t32 + g(2) * t30 + g(3) * t56;
t77 = t67 * pkin(3) * t98 + pkin(2) * t101 + t37 * pkin(4) + pkin(8) * t103 + pkin(9) * t93 + t36 * qJ(5);
t76 = -g(1) * t53 - g(2) * t51 + g(3) * t101;
t75 = t19 * pkin(4) + t52 * pkin(8) - pkin(9) * t105 + t18 * qJ(5) + t110 * t51;
t74 = t21 * pkin(4) + t54 * pkin(8) - pkin(9) * t104 + t20 * qJ(5) + t110 * t53;
t15 = t76 * t69;
t7 = t79 * t71;
t6 = t79 * t68;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3); 0, 0, -t76, g(1) * t54 + g(2) * t52 + g(3) * t103, 0, 0, 0, 0, 0, -t76 * t72, t15, 0, 0, 0, 0, 0, -t80, t81, -t15, t80, -t81, -g(1) * t74 - g(2) * t75 - g(3) * t77, -t15, -t81, -t80, -g(1) * (-pkin(5) * t104 + t21 * qJ(6) + t74) - g(2) * (-pkin(5) * t105 + t19 * qJ(6) + t75) - g(3) * (pkin(5) * t93 + t37 * qJ(6) + t77); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t78, 0, 0, 0, 0, 0, -t7, t6, -t78, t7, -t6, -g(1) * (t32 * pkin(9) + t90) - g(2) * (t30 * pkin(9) + t91) - g(3) * (t56 * pkin(9) + t89) -t78, -t6, -t7, -g(1) * (t108 * t32 + t31 * t96 + t90) - g(2) * (t108 * t30 + t29 * t96 + t91) - g(3) * (t108 * t56 + t55 * t96 + t89); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t82, 0, -t2, -t82, -g(1) * t85 - g(2) * t88 - g(3) * t84, 0, -t82, t2, -g(1) * (-t13 * qJ(6) + t85) - g(2) * (-t11 * qJ(6) + t88) - g(3) * (-t33 * qJ(6) + t84); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t82;];
taug_reg  = t1;
