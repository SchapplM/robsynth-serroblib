% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 00:13:35
% EndTime: 2019-05-06 00:13:40
% DurationCPUTime: 1.48s
% Computational Cost: add. (1126->150), mult. (3001->246), div. (0->0), fcn. (3874->16), ass. (0->82)
t106 = sin(pkin(7));
t118 = cos(qJ(1));
t105 = sin(pkin(12));
t116 = sin(qJ(1));
t108 = cos(pkin(12));
t110 = cos(pkin(6));
t92 = t110 * t108;
t77 = t105 * t116 - t118 * t92;
t107 = sin(pkin(6));
t109 = cos(pkin(7));
t89 = t109 * t107;
t121 = t106 * t77 - t118 * t89;
t117 = cos(qJ(3));
t88 = t107 * t106;
t126 = t109 * t77 + t118 * t88;
t90 = t110 * t105;
t38 = t108 * t116 + t118 * t90;
t58 = sin(qJ(3));
t22 = -t38 * t117 + t126 * t58;
t57 = sin(qJ(4));
t59 = cos(qJ(4));
t10 = -t121 * t57 + t22 * t59;
t19 = t117 * t126 + t38 * t58;
t53 = pkin(13) + qJ(6);
t50 = sin(t53);
t51 = cos(t53);
t130 = t10 * t50 + t19 * t51;
t129 = t10 * t51 - t19 * t50;
t9 = t121 * t59 + t22 * t57;
t72 = t105 * t118 + t116 * t92;
t125 = t106 * t72 + t116 * t89;
t122 = t109 * t72 - t116 * t88;
t120 = t106 * t110 + t108 * t89;
t115 = t50 * t59;
t114 = t51 * t59;
t54 = sin(pkin(13));
t113 = t54 * t59;
t55 = cos(pkin(13));
t112 = t55 * t59;
t95 = t107 * t116;
t111 = pkin(1) * t118 + qJ(2) * t95;
t104 = pkin(5) * t54 + pkin(10);
t13 = t19 * pkin(3);
t103 = -pkin(10) * t22 - t13;
t39 = t108 * t118 - t116 * t90;
t23 = t117 * t122 + t39 * t58;
t15 = t23 * pkin(3);
t24 = t117 * t39 - t122 * t58;
t102 = t24 * pkin(10) - t15;
t87 = t107 * t105;
t29 = -t117 * t120 + t58 * t87;
t28 = t29 * pkin(3);
t30 = t117 * t87 + t120 * t58;
t101 = t30 * pkin(10) - t28;
t11 = -t125 * t59 + t24 * t57;
t100 = g(1) * t9 + g(2) * t11;
t96 = t118 * t107;
t99 = -pkin(1) * t116 + qJ(2) * t96;
t98 = -g(1) * t19 + g(2) * t23;
t94 = -pkin(4) * t59 - qJ(5) * t57;
t49 = pkin(5) * t55 + pkin(4);
t56 = -pkin(11) - qJ(5);
t93 = -t49 * t59 + t56 * t57;
t71 = -t108 * t88 + t109 * t110;
t17 = t30 * t57 - t59 * t71;
t86 = g(1) * t11 - g(2) * t9 + g(3) * t17;
t12 = t125 * t57 + t24 * t59;
t18 = t30 * t59 + t57 * t71;
t85 = g(1) * t12 - g(2) * t10 + g(3) * t18;
t84 = g(1) * t23 + g(2) * t19 + g(3) * t29;
t83 = g(1) * t24 - g(2) * t22 + g(3) * t30;
t68 = -t38 * pkin(2) - pkin(9) * t121 + t99;
t65 = pkin(3) * t22 + t68;
t64 = t39 * pkin(2) + pkin(9) * t125 + t111;
t63 = pkin(3) * t24 + t64;
t62 = -pkin(10) * t19 + t65;
t61 = t23 * pkin(10) + t63;
t35 = -g(1) * t95 + g(2) * t96 - g(3) * t110;
t6 = t84 * t57;
t5 = t12 * t51 + t23 * t50;
t4 = -t12 * t50 + t23 * t51;
t1 = [0, 0, 0, 0, 0, 0, g(1) * t116 - g(2) * t118, g(1) * t118 + g(2) * t116, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t38 - g(2) * t39, -g(1) * t77 + g(2) * t72, -g(1) * t96 - g(2) * t95, -g(1) * t99 - g(2) * t111, 0, 0, 0, 0, 0, 0, -g(1) * t22 - g(2) * t24, t98, g(1) * t121 - g(2) * t125, -g(1) * t68 - g(2) * t64, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, t100, -t98, -g(1) * t62 - g(2) * t61, 0, 0, 0, 0, 0, 0, -g(1) * (t10 * t55 - t19 * t54) - g(2) * (t12 * t55 + t23 * t54) -g(1) * (-t10 * t54 - t19 * t55) - g(2) * (-t12 * t54 + t23 * t55) -t100, -g(1) * (t10 * pkin(4) + t9 * qJ(5) + t62) - g(2) * (t12 * pkin(4) + t11 * qJ(5) + t61) 0, 0, 0, 0, 0, 0, -g(1) * t129 - g(2) * t5, g(1) * t130 - g(2) * t4, -t100, -g(1) * (t10 * t49 - t104 * t19 - t9 * t56 + t65) - g(2) * (t104 * t23 - t11 * t56 + t12 * t49 + t63); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t83, 0, 0, 0, 0, 0, 0, 0, 0, t84 * t59, -t6, -t83, -g(1) * t102 - g(2) * t103 - g(3) * t101, 0, 0, 0, 0, 0, 0, -g(1) * (-t112 * t23 + t24 * t54) - g(2) * (-t112 * t19 - t22 * t54) - g(3) * (-t112 * t29 + t30 * t54) -g(1) * (t113 * t23 + t24 * t55) - g(2) * (t113 * t19 - t22 * t55) - g(3) * (t113 * t29 + t30 * t55) t6, -g(1) * (t23 * t94 + t102) - g(2) * (t19 * t94 + t103) - g(3) * (t29 * t94 + t101) 0, 0, 0, 0, 0, 0, -g(1) * (-t114 * t23 + t24 * t50) - g(2) * (-t114 * t19 - t22 * t50) - g(3) * (-t114 * t29 + t30 * t50) -g(1) * (t115 * t23 + t24 * t51) - g(2) * (t115 * t19 - t22 * t51) - g(3) * (t115 * t29 + t30 * t51) t6, -g(1) * (t104 * t24 + t23 * t93 - t15) - g(2) * (-t104 * t22 + t19 * t93 - t13) - g(3) * (t104 * t30 + t29 * t93 - t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t85, 0, 0, 0, 0, 0, 0, 0, 0, t86 * t55, -t86 * t54, -t85, -g(1) * (-pkin(4) * t11 + qJ(5) * t12) - g(2) * (pkin(4) * t9 - qJ(5) * t10) - g(3) * (-pkin(4) * t17 + qJ(5) * t18) 0, 0, 0, 0, 0, 0, t86 * t51, -t86 * t50, -t85, -g(1) * (-t11 * t49 - t12 * t56) - g(2) * (t10 * t56 + t49 * t9) - g(3) * (-t17 * t49 - t18 * t56); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t130 - g(3) * (-t18 * t50 + t29 * t51) g(1) * t5 - g(2) * t129 - g(3) * (-t18 * t51 - t29 * t50) 0, 0;];
taug_reg  = t1;
