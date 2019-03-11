% Calculate inertial parameters regressor of potential energy for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRPR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:49:33
% EndTime: 2019-03-08 19:49:33
% DurationCPUTime: 0.25s
% Computational Cost: add. (248->97), mult. (524->136), div. (0->0), fcn. (622->12), ass. (0->56)
t91 = sin(pkin(6));
t97 = sin(qJ(2));
t122 = t91 * t97;
t94 = cos(pkin(6));
t127 = t94 * pkin(3) + pkin(8) * t122;
t89 = sin(pkin(11));
t126 = pkin(5) * t89;
t119 = t94 * t97;
t90 = sin(pkin(10));
t93 = cos(pkin(10));
t99 = cos(qJ(2));
t70 = t119 * t93 + t90 * t99;
t125 = t70 * t89;
t124 = t90 * t91;
t96 = sin(qJ(4));
t123 = t91 * t96;
t98 = cos(qJ(4));
t121 = t91 * t98;
t120 = t91 * t99;
t118 = t94 * t99;
t117 = t93 * pkin(1) + pkin(7) * t124;
t116 = qJ(3) * t99;
t115 = t94 * pkin(7) + qJ(1);
t114 = t93 * t91 * pkin(7);
t113 = pkin(2) * t122 + t115;
t112 = (-pkin(3) - pkin(7)) * t93;
t69 = -t118 * t93 + t90 * t97;
t84 = t90 * pkin(1);
t111 = t70 * pkin(2) + t69 * qJ(3) + t84;
t110 = g(1) * t90 - g(2) * t93;
t71 = t118 * t90 + t93 * t97;
t72 = -t119 * t90 + t93 * t99;
t109 = t72 * pkin(2) + t71 * qJ(3) + t117;
t108 = pkin(3) * t124 + t109;
t107 = -t116 * t91 + t113;
t106 = t70 * pkin(8) + t111;
t61 = t123 * t90 - t71 * t98;
t63 = t123 * t93 + t69 * t98;
t73 = t120 * t98 + t94 * t96;
t105 = g(1) * t61 - g(2) * t63 + g(3) * t73;
t104 = -g(1) * t71 - g(2) * t69 + g(3) * t120;
t103 = g(1) * t72 + g(2) * t70 + g(3) * t122;
t102 = t107 + t127;
t101 = t72 * pkin(8) + t108;
t100 = t112 * t91 + t106;
t95 = -pkin(9) - qJ(5);
t92 = cos(pkin(11));
t88 = pkin(11) + qJ(6);
t83 = cos(t88);
t82 = sin(t88);
t81 = t92 * pkin(5) + pkin(4);
t74 = -t120 * t96 + t94 * t98;
t65 = -g(3) * t94 - t110 * t91;
t64 = -t121 * t93 + t69 * t96;
t62 = t121 * t90 + t71 * t96;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t93 - g(2) * t90, t110, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t103, -t104, t65, -g(1) * t117 - g(2) * (t84 - t114) - g(3) * t115, 0, 0, 0, 0, 0, 0, t65, t103, t104, -g(1) * t109 - g(2) * (t111 - t114) - g(3) * t107, 0, 0, 0, 0, 0, 0, -g(1) * t62 - g(2) * t64 - g(3) * t74, t105, -t103, -g(1) * t101 - g(2) * t100 - g(3) * t102, 0, 0, 0, 0, 0, 0, -g(1) * (t62 * t92 + t72 * t89) - g(2) * (t64 * t92 + t125) - g(3) * (t122 * t89 + t74 * t92) -g(1) * (-t62 * t89 + t72 * t92) - g(2) * (-t64 * t89 + t70 * t92) - g(3) * (t122 * t92 - t74 * t89) -t105, -g(1) * (t62 * pkin(4) + t61 * qJ(5) + t101) - g(2) * (t64 * pkin(4) - t63 * qJ(5) + t100) - g(3) * (t74 * pkin(4) + t73 * qJ(5) + t102) 0, 0, 0, 0, 0, 0, -g(1) * (t62 * t83 + t72 * t82) - g(2) * (t64 * t83 + t70 * t82) - g(3) * (t122 * t82 + t74 * t83) -g(1) * (-t62 * t82 + t72 * t83) - g(2) * (-t64 * t82 + t70 * t83) - g(3) * (t122 * t83 - t74 * t82) -t105, -g(1) * (-t61 * t95 + t62 * t81 + (pkin(8) + t126) * t72 + t108) - g(2) * (pkin(5) * t125 + t63 * t95 + t64 * t81 + t106) - g(3) * (-t73 * t95 + t74 * t81 + t113 + t127) + (-g(3) * (t126 * t97 - t116) - g(2) * t112) * t91;];
U_reg  = t1;
