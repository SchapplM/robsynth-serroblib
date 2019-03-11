% Calculate inertial parameters regressor of potential energy for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:06:16
% EndTime: 2019-03-09 10:06:16
% DurationCPUTime: 0.13s
% Computational Cost: add. (131->66), mult. (255->71), div. (0->0), fcn. (261->6), ass. (0->42)
t87 = sin(qJ(2));
t103 = qJ(3) * t87;
t88 = sin(qJ(1));
t90 = cos(qJ(2));
t109 = t88 * t90;
t116 = pkin(2) * t109 + t88 * t103;
t86 = sin(qJ(4));
t91 = cos(qJ(1));
t106 = t91 * t86;
t89 = cos(qJ(4));
t110 = t88 * t89;
t68 = t87 * t110 + t106;
t105 = t91 * t89;
t111 = t88 * t86;
t69 = t87 * t111 - t105;
t115 = t69 * pkin(4) - t68 * qJ(5);
t114 = g(3) * pkin(6);
t113 = g(3) * t90;
t112 = t87 * pkin(2) + pkin(6);
t108 = t89 * t90;
t107 = t90 * t91;
t104 = t91 * pkin(1) + t88 * pkin(7);
t83 = t88 * pkin(1);
t101 = -t91 * pkin(7) + t83;
t79 = t87 * pkin(8);
t100 = qJ(5) * t108 + t112 + t79;
t99 = pkin(2) * t107 + t91 * t103 + t104;
t98 = -t90 * qJ(3) + t112;
t97 = g(1) * t91 + g(2) * t88;
t96 = t101 + t116;
t95 = t88 * pkin(3) + pkin(8) * t107 + t99;
t75 = pkin(8) * t109;
t94 = t75 + t83 + (-pkin(3) - pkin(7)) * t91 + t116;
t66 = -t87 * t105 + t111;
t93 = g(1) * t66 - g(2) * t68 + g(3) * t108;
t67 = t87 * t106 + t110;
t92 = t67 * pkin(4) + t66 * qJ(5) + t95;
t70 = g(1) * t88 - g(2) * t91;
t63 = g(3) * t87 + t97 * t90;
t62 = t97 * t87 - t113;
t61 = -g(1) * t67 - g(2) * t69 + t86 * t113;
t1 = [0, 0, 0, 0, 0, 0, -t97, t70, -g(3), -t114, 0, 0, 0, 0, 0, 0, -t63, t62, -t70, -g(1) * t104 - g(2) * t101 - t114, 0, 0, 0, 0, 0, 0, -t70, t63, -t62, -g(1) * t99 - g(2) * t96 - g(3) * t98, 0, 0, 0, 0, 0, 0, t61, t93, -t63, -g(1) * t95 - g(2) * t94 - g(3) * (t79 + t98) 0, 0, 0, 0, 0, 0, t61, -t63, -t93, -g(1) * t92 - g(2) * (t94 + t115) - g(3) * ((-pkin(4) * t86 - qJ(3)) * t90 + t100) 0, 0, 0, 0, 0, 0, t61, -t93, t63, -g(1) * (t67 * pkin(5) + t92) - g(2) * (-t91 * pkin(3) + t69 * pkin(5) + t115 + t75 + t96) - g(3) * (-t87 * qJ(6) + t100) + (-g(3) * (-qJ(3) + (-pkin(4) - pkin(5)) * t86) + t97 * qJ(6)) * t90;];
U_reg  = t1;
