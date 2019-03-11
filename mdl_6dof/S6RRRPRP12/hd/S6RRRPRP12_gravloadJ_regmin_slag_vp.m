% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t71 = sin(qJ(3));
t115 = qJ(4) * t71 + pkin(2);
t72 = sin(qJ(2));
t73 = sin(qJ(1));
t76 = cos(qJ(2));
t110 = cos(qJ(1));
t96 = cos(pkin(6));
t85 = t96 * t110;
t52 = t72 * t85 + t73 * t76;
t75 = cos(qJ(3));
t69 = sin(pkin(6));
t92 = t69 * t110;
t30 = t52 * t71 + t75 * t92;
t51 = t73 * t72 - t76 * t85;
t70 = sin(qJ(5));
t74 = cos(qJ(5));
t10 = t30 * t70 + t51 * t74;
t11 = t30 * t74 - t51 * t70;
t31 = t52 * t75 - t71 * t92;
t104 = t69 * t73;
t91 = t73 * t96;
t54 = t110 * t76 - t72 * t91;
t35 = t104 * t71 + t54 * t75;
t103 = t69 * t75;
t50 = t103 * t72 + t71 * t96;
t80 = g(1) * t35 + g(2) * t31 + g(3) * t50;
t114 = pkin(4) + pkin(9);
t107 = t51 * t75;
t53 = t110 * t72 + t76 * t91;
t106 = t53 * t75;
t105 = t69 * t72;
t102 = t69 * t76;
t101 = t70 * t71;
t100 = t70 * t76;
t99 = t71 * t74;
t98 = t75 * t76;
t95 = t74 * t102;
t94 = -pkin(3) * t107 - t115 * t51;
t93 = -pkin(3) * t106 - t115 * t53;
t90 = t69 * pkin(3) * t98 + pkin(9) * t105 + t115 * t102;
t34 = -t103 * t73 + t54 * t71;
t13 = -t34 * t74 + t53 * t70;
t89 = g(1) * t11 + g(2) * t13;
t88 = -g(1) * t30 + g(2) * t34;
t87 = -g(1) * t31 + g(2) * t35;
t86 = g(1) * t51 - g(2) * t53;
t83 = t110 * pkin(1) + t54 * pkin(2) + t35 * pkin(3) + pkin(8) * t104 + t34 * qJ(4);
t49 = t105 * t71 - t75 * t96;
t28 = t100 * t69 + t49 * t74;
t1 = g(1) * t13 - g(2) * t11 - g(3) * t28;
t14 = t34 * t70 + t53 * t74;
t29 = -t49 * t70 + t95;
t82 = g(1) * t14 + g(2) * t10 - g(3) * t29;
t17 = t51 * t99 + t52 * t70;
t19 = t53 * t99 + t54 * t70;
t38 = t105 * t70 - t71 * t95;
t81 = g(1) * t19 + g(2) * t17 + g(3) * t38;
t7 = g(1) * t34 + g(2) * t30 + g(3) * t49;
t79 = -t73 * pkin(1) - t52 * pkin(2) - pkin(3) * t31 + pkin(8) * t92 - qJ(4) * t30;
t78 = -g(1) * t53 - g(2) * t51 + g(3) * t102;
t77 = g(1) * t54 + g(2) * t52 + g(3) * t105;
t44 = t49 * pkin(3);
t39 = (t100 * t71 + t72 * t74) * t69;
t26 = t34 * pkin(3);
t24 = t30 * pkin(3);
t20 = -t101 * t53 + t54 * t74;
t18 = -t101 * t51 + t52 * t74;
t16 = t78 * t75;
t15 = t78 * t71;
t5 = t80 * t74;
t4 = t80 * t70;
t3 = g(1) * t10 - g(2) * t14;
t2 = -g(1) * t20 - g(2) * t18 - g(3) * t39;
t6 = [0, g(1) * t73 - g(2) * t110, g(1) * t110 + g(2) * t73, 0, 0, 0, 0, 0, g(1) * t52 - g(2) * t54, -t86, 0, 0, 0, 0, 0, -t87, t88, t86, t87, -t88, -g(1) * (-t51 * pkin(9) + t79) - g(2) * (t53 * pkin(9) + t83) 0, 0, 0, 0, 0, t3, t89, t3, -t87, -t89, -g(1) * (-pkin(5) * t10 - pkin(10) * t31 + t11 * qJ(6) - t114 * t51 + t79) - g(2) * (t14 * pkin(5) + t35 * pkin(10) + t13 * qJ(6) + t114 * t53 + t83); 0, 0, 0, 0, 0, 0, 0, 0, -t78, t77, 0, 0, 0, 0, 0, -t16, t15, -t77, t16, -t15, -g(1) * (t54 * pkin(9) + t93) - g(2) * (t52 * pkin(9) + t94) - g(3) * t90, 0, 0, 0, 0, 0, t2, t81, t2, -t16, -t81, -g(1) * (t20 * pkin(5) - pkin(10) * t106 + t19 * qJ(6) + t114 * t54 + t93) - g(2) * (t18 * pkin(5) - pkin(10) * t107 + t17 * qJ(6) + t114 * t52 + t94) - g(3) * (t39 * pkin(5) + t38 * qJ(6) + (pkin(4) * t72 + pkin(10) * t98) * t69 + t90); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t80, 0, -t7, -t80, -g(1) * (t35 * qJ(4) - t26) - g(2) * (t31 * qJ(4) - t24) - g(3) * (t50 * qJ(4) - t44) 0, 0, 0, 0, 0, -t4, -t5, -t4, t7, t5, -g(1) * (-t34 * pkin(10) - t26) - g(2) * (-t30 * pkin(10) - t24) - g(3) * (-t49 * pkin(10) - t44) - t80 * (pkin(5) * t70 - qJ(6) * t74 + qJ(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t82, t1, 0, -t82, -g(1) * (-t13 * pkin(5) + t14 * qJ(6)) - g(2) * (pkin(5) * t11 + t10 * qJ(6)) - g(3) * (t28 * pkin(5) - t29 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t6;
