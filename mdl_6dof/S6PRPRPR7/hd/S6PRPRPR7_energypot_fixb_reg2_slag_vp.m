% Calculate inertial parameters regressor of potential energy for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRPR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:53:46
% EndTime: 2019-03-08 19:53:46
% DurationCPUTime: 0.22s
% Computational Cost: add. (225->86), mult. (502->113), div. (0->0), fcn. (592->10), ass. (0->53)
t88 = sin(pkin(6));
t93 = sin(qJ(2));
t122 = t88 * t93;
t90 = cos(pkin(6));
t127 = t90 * pkin(3) + pkin(8) * t122;
t92 = sin(qJ(4));
t123 = t88 * t92;
t96 = cos(qJ(2));
t118 = t90 * t96;
t87 = sin(pkin(10));
t89 = cos(pkin(10));
t70 = -t89 * t118 + t87 * t93;
t95 = cos(qJ(4));
t62 = t89 * t123 + t70 * t95;
t121 = t88 * t95;
t63 = t89 * t121 - t70 * t92;
t126 = -t63 * pkin(4) - t62 * qJ(5);
t119 = t90 * t93;
t73 = -t87 * t119 + t89 * t96;
t125 = t73 * pkin(8);
t124 = t87 * t88;
t120 = t88 * t96;
t117 = t89 * pkin(1) + pkin(7) * t124;
t116 = qJ(3) * t96;
t114 = t90 * pkin(7) + qJ(1);
t113 = t89 * t88 * pkin(7);
t112 = t88 * t116;
t111 = pkin(2) * t122 + t114;
t110 = (-pkin(3) - pkin(7)) * t89;
t71 = t89 * t119 + t87 * t96;
t83 = t87 * pkin(1);
t109 = t71 * pkin(2) + t70 * qJ(3) + t83;
t108 = g(1) * t87 - g(2) * t89;
t72 = t87 * t118 + t89 * t93;
t107 = t73 * pkin(2) + t72 * qJ(3) + t117;
t106 = pkin(3) * t124 + t107;
t105 = t111 - t112;
t104 = t71 * pkin(8) + t109;
t60 = t87 * t123 - t72 * t95;
t74 = t95 * t120 + t90 * t92;
t103 = g(1) * t60 - g(2) * t62 + g(3) * t74;
t61 = t87 * t121 + t72 * t92;
t75 = -t92 * t120 + t90 * t95;
t102 = g(1) * t61 - g(2) * t63 + g(3) * t75;
t101 = -g(1) * t72 - g(2) * t70 + g(3) * t120;
t100 = g(1) * t73 + g(2) * t71 + g(3) * t122;
t99 = t75 * pkin(4) + t74 * qJ(5) + t111 + t127;
t98 = t61 * pkin(4) + t60 * qJ(5) + t106;
t97 = t88 * t110 + t104;
t94 = cos(qJ(6));
t91 = sin(qJ(6));
t64 = -g(3) * t90 - t108 * t88;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t89 - g(2) * t87, t108, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t100, -t101, t64, -g(1) * t117 - g(2) * (t83 - t113) - g(3) * t114, 0, 0, 0, 0, 0, 0, t64, t100, t101, -g(1) * t107 - g(2) * (t109 - t113) - g(3) * t105, 0, 0, 0, 0, 0, 0, -t102, t103, -t100, -g(1) * (t106 + t125) - g(2) * t97 - g(3) * (t105 + t127) 0, 0, 0, 0, 0, 0, -t100, t102, -t103, -g(1) * (t98 + t125) - g(2) * (t97 + t126) - g(3) * (t99 - t112) 0, 0, 0, 0, 0, 0, -g(1) * (t60 * t91 + t73 * t94) - g(2) * (-t62 * t91 + t71 * t94) - g(3) * (t94 * t122 + t74 * t91) -g(1) * (t60 * t94 - t73 * t91) - g(2) * (-t62 * t94 - t71 * t91) - g(3) * (-t91 * t122 + t74 * t94) -t102, -g(1) * (t61 * pkin(9) + (pkin(5) + pkin(8)) * t73 + t98) - g(2) * (t71 * pkin(5) - t63 * pkin(9) + t104 + t126) - g(3) * (t75 * pkin(9) + t99) + (-g(3) * (pkin(5) * t93 - t116) - g(2) * t110) * t88;];
U_reg  = t1;
