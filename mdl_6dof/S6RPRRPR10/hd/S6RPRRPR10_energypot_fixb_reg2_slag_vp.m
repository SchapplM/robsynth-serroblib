% Calculate inertial parameters regressor of potential energy for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:37:09
% EndTime: 2019-03-09 05:37:09
% DurationCPUTime: 0.15s
% Computational Cost: add. (130->71), mult. (250->86), div. (0->0), fcn. (266->8), ass. (0->45)
t81 = sin(qJ(4));
t87 = cos(qJ(1));
t102 = t87 * t81;
t83 = sin(qJ(1));
t85 = cos(qJ(4));
t106 = t83 * t85;
t82 = sin(qJ(3));
t62 = t82 * t102 + t106;
t101 = t87 * t85;
t107 = t83 * t81;
t63 = -t82 * t101 + t107;
t114 = t63 * pkin(4) - t62 * qJ(5);
t113 = g(3) * pkin(6);
t112 = pkin(2) + pkin(6);
t111 = pkin(3) * t82;
t110 = g(2) * t87;
t86 = cos(qJ(3));
t109 = g(3) * t86;
t108 = t81 * t86;
t105 = t83 * t86;
t104 = t85 * t86;
t103 = t86 * t87;
t100 = t87 * pkin(1) + t83 * qJ(2);
t98 = pkin(8) * t105;
t74 = t83 * pkin(7);
t75 = t83 * pkin(1);
t97 = pkin(8) * t103 + t74 + t75;
t96 = t87 * pkin(7) + t100;
t95 = t86 * pkin(3) + t82 * pkin(8) + t112;
t94 = -qJ(2) - t111;
t93 = -t87 * qJ(2) + t75;
t92 = t83 * t111 + t96;
t64 = g(1) * t83 - t110;
t91 = pkin(4) * t104 + qJ(5) * t108 + t95;
t60 = t82 * t107 - t101;
t61 = t82 * t106 + t102;
t90 = t61 * pkin(4) + t60 * qJ(5) + t92;
t89 = g(1) * t60 - g(2) * t62 + g(3) * t108;
t88 = t94 * t87 + t97;
t84 = cos(qJ(6));
t80 = sin(qJ(6));
t65 = g(1) * t87 + g(2) * t83;
t57 = g(1) * t105 - g(2) * t103 - g(3) * t82;
t56 = -g(1) * t61 - g(2) * t63 - g(3) * t104;
t1 = [0, 0, 0, 0, 0, 0, -t65, t64, -g(3), -t113, 0, 0, 0, 0, 0, 0, -g(3), t65, -t64, -g(1) * t100 - g(2) * t93 - t113, 0, 0, 0, 0, 0, 0, -t64 * t82 - t109, -t57, -t65, -g(1) * t96 - g(2) * (t74 + t93) - g(3) * t112, 0, 0, 0, 0, 0, 0, t56, t89, t57, -g(1) * (t92 - t98) - g(2) * t88 - g(3) * t95, 0, 0, 0, 0, 0, 0, t56, t57, -t89, -g(1) * (t90 - t98) - g(2) * (t88 + t114) - g(3) * t91, 0, 0, 0, 0, 0, 0, -g(1) * (t60 * t80 + t61 * t84) - g(2) * (-t62 * t80 + t63 * t84) - (t80 * t81 + t84 * t85) * t109, -g(1) * (t60 * t84 - t61 * t80) - g(2) * (-t62 * t84 - t63 * t80) - (-t80 * t85 + t81 * t84) * t109, -t57, -g(1) * (t61 * pkin(5) + (-pkin(8) + pkin(9)) * t105 + t90) - g(2) * (t63 * pkin(5) + t114 + t97) - g(3) * (pkin(5) * t104 - t82 * pkin(9) + t91) - (-pkin(9) * t86 + t94) * t110;];
U_reg  = t1;
