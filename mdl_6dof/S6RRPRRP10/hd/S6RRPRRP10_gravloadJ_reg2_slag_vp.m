% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:46:13
% EndTime: 2019-05-06 18:46:17
% DurationCPUTime: 0.80s
% Computational Cost: add. (787->138), mult. (1415->208), div. (0->0), fcn. (1730->12), ass. (0->83)
t78 = pkin(11) + qJ(4);
t75 = sin(t78);
t76 = cos(t78);
t131 = pkin(4) * t76 + pkin(10) * t75;
t128 = cos(qJ(1));
t80 = sin(pkin(6));
t113 = t80 * t128;
t84 = sin(qJ(2));
t85 = sin(qJ(1));
t87 = cos(qJ(2));
t116 = cos(pkin(6));
t99 = t116 * t128;
t57 = t84 * t99 + t85 * t87;
t29 = -t75 * t113 + t57 * t76;
t56 = t85 * t84 - t87 * t99;
t83 = sin(qJ(5));
t86 = cos(qJ(5));
t9 = t29 * t83 - t56 * t86;
t10 = t29 * t86 + t56 * t83;
t125 = t76 * t83;
t124 = t76 * t86;
t123 = t80 * t84;
t122 = t80 * t85;
t121 = t80 * t87;
t120 = t86 * t87;
t81 = cos(pkin(11));
t74 = t81 * pkin(3) + pkin(2);
t82 = -pkin(9) - qJ(3);
t119 = -t56 * t74 - t57 * t82;
t107 = t85 * t116;
t58 = t87 * t107 + t128 * t84;
t59 = -t84 * t107 + t128 * t87;
t118 = -t58 * t74 - t59 * t82;
t117 = t128 * pkin(1) + pkin(8) * t122;
t115 = t83 * t121;
t79 = sin(pkin(11));
t114 = t79 * t122;
t112 = -t85 * pkin(1) + pkin(8) * t113;
t108 = -t76 * t113 - t57 * t75;
t111 = pkin(4) * t108 + t29 * pkin(10);
t32 = -t76 * t122 + t59 * t75;
t33 = t75 * t122 + t59 * t76;
t110 = -t32 * pkin(4) + t33 * pkin(10);
t45 = t116 * t76 - t75 * t123;
t46 = t116 * t75 + t76 * t123;
t109 = t45 * pkin(4) + t46 * pkin(10);
t106 = -t131 * t56 + t119;
t105 = -t131 * t58 + t118;
t104 = t79 * t113;
t103 = t74 * t121 - t82 * t123;
t102 = pkin(3) * t114 - t58 * t82 + t59 * t74 + t117;
t13 = t33 * t83 - t58 * t86;
t101 = -g(1) * t9 + g(2) * t13;
t100 = g(1) * t108 + g(2) * t32;
t21 = g(1) * t56 - g(2) * t58;
t98 = pkin(5) * t86 + qJ(6) * t83;
t97 = g(1) * t128 + g(2) * t85;
t96 = pkin(3) * t104 + t56 * t82 - t57 * t74 + t112;
t95 = t131 * t121 + t103;
t26 = t80 * t120 + t46 * t83;
t1 = g(1) * t13 + g(2) * t9 + g(3) * t26;
t14 = t33 * t86 + t58 * t83;
t27 = t46 * t86 - t115;
t94 = g(1) * t14 + g(2) * t10 + g(3) * t27;
t15 = -t56 * t125 - t57 * t86;
t17 = -t58 * t125 - t59 * t86;
t34 = t76 * t115 - t86 * t123;
t93 = g(1) * t17 + g(2) * t15 + g(3) * t34;
t92 = g(1) * t32 - g(2) * t108 - g(3) * t45;
t91 = g(1) * t33 + g(2) * t29 + g(3) * t46;
t19 = -g(1) * t58 - g(2) * t56 + g(3) * t121;
t90 = g(1) * t59 + g(2) * t57 + g(3) * t123;
t89 = t33 * pkin(4) + t32 * pkin(10) + t102;
t88 = -pkin(4) * t29 + pkin(10) * t108 + t96;
t35 = (t76 * t120 + t83 * t84) * t80;
t18 = -t58 * t124 + t59 * t83;
t16 = -t56 * t124 + t57 * t83;
t8 = t19 * t75;
t5 = t92 * t86;
t4 = t92 * t83;
t3 = -g(1) * t18 - g(2) * t16 - g(3) * t35;
t2 = g(1) * t10 - g(2) * t14;
t6 = [0, 0, 0, 0, 0, 0, g(1) * t85 - g(2) * t128, t97, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t57 - g(2) * t59, -t21, -t97 * t80, -g(1) * t112 - g(2) * t117, 0, 0, 0, 0, 0, 0, -g(1) * (-t57 * t81 + t104) - g(2) * (t59 * t81 + t114) -g(1) * (t81 * t113 + t57 * t79) - g(2) * (t81 * t122 - t59 * t79) t21, -g(1) * (-t57 * pkin(2) - t56 * qJ(3) + t112) - g(2) * (t59 * pkin(2) + t58 * qJ(3) + t117) 0, 0, 0, 0, 0, 0, g(1) * t29 - g(2) * t33, t100, t21, -g(1) * t96 - g(2) * t102, 0, 0, 0, 0, 0, 0, t2, t101, -t100, -g(1) * t88 - g(2) * t89, 0, 0, 0, 0, 0, 0, t2, -t100, -t101, -g(1) * (-pkin(5) * t10 - qJ(6) * t9 + t88) - g(2) * (t14 * pkin(5) + t13 * qJ(6) + t89); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t90, 0, 0, 0, 0, 0, 0, 0, 0, -t19 * t81, t19 * t79, -t90, -g(1) * (-t58 * pkin(2) + t59 * qJ(3)) - g(2) * (-t56 * pkin(2) + t57 * qJ(3)) - g(3) * (pkin(2) * t87 + qJ(3) * t84) * t80, 0, 0, 0, 0, 0, 0, -t19 * t76, t8, -t90, -g(1) * t118 - g(2) * t119 - g(3) * t103, 0, 0, 0, 0, 0, 0, t3, t93, -t8, -g(1) * t105 - g(2) * t106 - g(3) * t95, 0, 0, 0, 0, 0, 0, t3, -t8, -t93, -g(1) * (t18 * pkin(5) + t17 * qJ(6) + t105) - g(2) * (t16 * pkin(5) + t15 * qJ(6) + t106) - g(3) * (t35 * pkin(5) + t34 * qJ(6) + t95); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, t91, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t4, -t91, -g(1) * t110 - g(2) * t111 - g(3) * t109, 0, 0, 0, 0, 0, 0, t5, -t91, t4, -g(1) * (-t98 * t32 + t110) - g(2) * (t108 * t98 + t111) - g(3) * (t98 * t45 + t109); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t94, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t94, -g(1) * (-t13 * pkin(5) + t14 * qJ(6)) - g(2) * (-t9 * pkin(5) + t10 * qJ(6)) - g(3) * (-t26 * pkin(5) + t27 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t6;
