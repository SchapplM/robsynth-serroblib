% Calculate inertial parameters regressor of potential energy for
% S6RRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:56:28
% EndTime: 2019-03-09 20:56:28
% DurationCPUTime: 0.15s
% Computational Cost: add. (199->64), mult. (236->74), div. (0->0), fcn. (242->8), ass. (0->41)
t89 = qJ(2) + qJ(3);
t85 = sin(t89);
t93 = cos(qJ(4));
t113 = t85 * t93;
t90 = sin(qJ(4));
t115 = t85 * t90;
t119 = pkin(4) * t113 + qJ(5) * t115;
t118 = g(3) * pkin(6);
t86 = cos(t89);
t117 = pkin(3) * t86;
t91 = sin(qJ(2));
t116 = t91 * pkin(2) + pkin(6);
t92 = sin(qJ(1));
t114 = t85 * t92;
t95 = cos(qJ(1));
t112 = t85 * t95;
t111 = t90 * t92;
t110 = t92 * t93;
t109 = t93 * t95;
t108 = t95 * t90;
t94 = cos(qJ(2));
t83 = pkin(2) * t94 + pkin(1);
t96 = -pkin(8) - pkin(7);
t107 = t92 * t83 + t95 * t96;
t106 = t85 * pkin(3) + t116;
t105 = t95 * t83 - t92 * t96;
t104 = pkin(9) * t114 + t92 * t117 + t107;
t103 = g(1) * t95 + g(2) * t92;
t102 = -pkin(9) * t86 + t106;
t101 = pkin(9) * t112 + t95 * t117 + t105;
t66 = t86 * t111 + t109;
t67 = t86 * t110 - t108;
t100 = t67 * pkin(4) + qJ(5) * t66 + t104;
t68 = t86 * t108 - t110;
t99 = g(1) * t68 + g(2) * t66 + g(3) * t115;
t69 = t86 * t109 + t111;
t98 = g(1) * t69 + g(2) * t67 + g(3) * t113;
t97 = t69 * pkin(4) + t68 * qJ(5) + t101;
t70 = g(1) * t92 - g(2) * t95;
t63 = -g(3) * t86 + t103 * t85;
t1 = [0, 0, 0, 0, 0, 0, -t103, t70, -g(3), -t118, 0, 0, 0, 0, 0, 0, -g(3) * t91 - t103 * t94, -g(3) * t94 + t103 * t91, -t70, -g(1) * (pkin(1) * t95 + pkin(7) * t92) - g(2) * (pkin(1) * t92 - pkin(7) * t95) - t118, 0, 0, 0, 0, 0, 0, -g(3) * t85 - t103 * t86, t63, -t70, -g(1) * t105 - g(2) * t107 - g(3) * t116, 0, 0, 0, 0, 0, 0, -t98, t99, -t63, -g(1) * t101 - g(2) * t104 - g(3) * t102, 0, 0, 0, 0, 0, 0, -t63, t98, -t99, -g(1) * t97 - g(2) * t100 - g(3) * (t102 + t119) 0, 0, 0, 0, 0, 0, -t63, -t99, -t98, -g(1) * (pkin(5) * t112 + t69 * qJ(6) + t97) - g(2) * (pkin(5) * t114 + qJ(6) * t67 + t100) - g(3) * (qJ(6) * t113 + (-pkin(5) - pkin(9)) * t86 + t106 + t119);];
U_reg  = t1;
