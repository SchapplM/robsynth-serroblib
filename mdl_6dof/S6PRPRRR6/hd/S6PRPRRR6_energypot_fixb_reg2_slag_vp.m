% Calculate inertial parameters regressor of potential energy for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:48:45
% EndTime: 2019-03-08 20:48:45
% DurationCPUTime: 0.25s
% Computational Cost: add. (248->97), mult. (524->137), div. (0->0), fcn. (622->12), ass. (0->56)
t90 = sin(pkin(6));
t95 = sin(qJ(2));
t123 = t90 * t95;
t92 = cos(pkin(6));
t127 = t92 * pkin(3) + pkin(8) * t123;
t120 = t92 * t95;
t89 = sin(pkin(11));
t91 = cos(pkin(11));
t98 = cos(qJ(2));
t70 = t91 * t120 + t89 * t98;
t93 = sin(qJ(5));
t126 = t70 * t93;
t125 = t89 * t90;
t94 = sin(qJ(4));
t124 = t90 * t94;
t97 = cos(qJ(4));
t122 = t90 * t97;
t121 = t90 * t98;
t119 = t92 * t98;
t118 = t93 * t95;
t117 = t91 * pkin(1) + pkin(7) * t125;
t116 = qJ(3) * t98;
t115 = t92 * pkin(7) + qJ(1);
t114 = t91 * t90 * pkin(7);
t113 = pkin(2) * t123 + t115;
t112 = (-pkin(3) - pkin(7)) * t91;
t69 = -t91 * t119 + t89 * t95;
t84 = t89 * pkin(1);
t111 = t70 * pkin(2) + t69 * qJ(3) + t84;
t110 = g(1) * t89 - g(2) * t91;
t71 = t89 * t119 + t91 * t95;
t72 = -t89 * t120 + t91 * t98;
t109 = t72 * pkin(2) + t71 * qJ(3) + t117;
t108 = pkin(3) * t125 + t109;
t107 = -t90 * t116 + t113;
t106 = t70 * pkin(8) + t111;
t61 = t89 * t124 - t71 * t97;
t63 = t91 * t124 + t69 * t97;
t73 = t97 * t121 + t92 * t94;
t105 = g(1) * t61 - g(2) * t63 + g(3) * t73;
t104 = -g(1) * t71 - g(2) * t69 + g(3) * t121;
t103 = g(1) * t72 + g(2) * t70 + g(3) * t123;
t102 = t107 + t127;
t101 = t72 * pkin(8) + t108;
t100 = t90 * t112 + t106;
t99 = -pkin(10) - pkin(9);
t96 = cos(qJ(5));
t88 = qJ(5) + qJ(6);
t83 = cos(t88);
t82 = sin(t88);
t81 = pkin(5) * t96 + pkin(4);
t74 = -t94 * t121 + t92 * t97;
t65 = -g(3) * t92 - t110 * t90;
t64 = -t91 * t122 + t69 * t94;
t62 = t89 * t122 + t71 * t94;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t91 - g(2) * t89, t110, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t103, -t104, t65, -g(1) * t117 - g(2) * (t84 - t114) - g(3) * t115, 0, 0, 0, 0, 0, 0, t65, t103, t104, -g(1) * t109 - g(2) * (t111 - t114) - g(3) * t107, 0, 0, 0, 0, 0, 0, -g(1) * t62 - g(2) * t64 - g(3) * t74, t105, -t103, -g(1) * t101 - g(2) * t100 - g(3) * t102, 0, 0, 0, 0, 0, 0, -g(1) * (t62 * t96 + t72 * t93) - g(2) * (t64 * t96 + t126) - g(3) * (t90 * t118 + t74 * t96) -g(1) * (-t62 * t93 + t72 * t96) - g(2) * (-t64 * t93 + t70 * t96) - g(3) * (t96 * t123 - t74 * t93) -t105, -g(1) * (pkin(4) * t62 + pkin(9) * t61 + t101) - g(2) * (t64 * pkin(4) - t63 * pkin(9) + t100) - g(3) * (t74 * pkin(4) + t73 * pkin(9) + t102) 0, 0, 0, 0, 0, 0, -g(1) * (t62 * t83 + t72 * t82) - g(2) * (t64 * t83 + t70 * t82) - g(3) * (t82 * t123 + t74 * t83) -g(1) * (-t62 * t82 + t72 * t83) - g(2) * (-t64 * t82 + t70 * t83) - g(3) * (t83 * t123 - t74 * t82) -t105, -g(1) * (-t61 * t99 + t62 * t81 + (pkin(5) * t93 + pkin(8)) * t72 + t108) - g(2) * (pkin(5) * t126 + t63 * t99 + t64 * t81 + t106) - g(3) * (-t73 * t99 + t74 * t81 + t113 + t127) + (-g(3) * (pkin(5) * t118 - t116) - g(2) * t112) * t90;];
U_reg  = t1;
