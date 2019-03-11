% Calculate inertial parameters regressor of potential energy for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:03:59
% EndTime: 2019-03-09 22:03:59
% DurationCPUTime: 0.11s
% Computational Cost: add. (193->63), mult. (161->75), div. (0->0), fcn. (149->10), ass. (0->38)
t93 = qJ(2) + qJ(3);
t88 = qJ(4) + t93;
t83 = sin(t88);
t106 = qJ(5) * t83;
t84 = cos(t88);
t99 = cos(qJ(1));
t112 = t84 * t99;
t117 = pkin(4) * t112 + t99 * t106;
t116 = g(3) * pkin(6);
t100 = -pkin(8) - pkin(7);
t115 = g(3) * t84;
t95 = sin(qJ(2));
t114 = t95 * pkin(2) + pkin(6);
t98 = cos(qJ(2));
t85 = t98 * pkin(2) + pkin(1);
t96 = sin(qJ(1));
t113 = t84 * t96;
t94 = sin(qJ(6));
t111 = t96 * t94;
t97 = cos(qJ(6));
t110 = t96 * t97;
t109 = t99 * t94;
t108 = t99 * t97;
t87 = cos(t93);
t74 = pkin(3) * t87 + t85;
t92 = -pkin(9) + t100;
t107 = t96 * t74 + t99 * t92;
t86 = sin(t93);
t105 = pkin(3) * t86 + t114;
t73 = t99 * t74;
t104 = -t96 * t92 + t73;
t103 = pkin(4) * t113 + t96 * t106 + t107;
t102 = g(1) * t99 + g(2) * t96;
t101 = t83 * pkin(4) - t84 * qJ(5) + t105;
t79 = g(1) * t96 - g(2) * t99;
t71 = g(3) * t83 + t102 * t84;
t70 = t102 * t83 - t115;
t1 = [0, 0, 0, 0, 0, 0, -t102, t79, -g(3), -t116, 0, 0, 0, 0, 0, 0, -g(3) * t95 - t102 * t98, -g(3) * t98 + t102 * t95, -t79, -g(1) * (t99 * pkin(1) + t96 * pkin(7)) - g(2) * (t96 * pkin(1) - t99 * pkin(7)) - t116, 0, 0, 0, 0, 0, 0, -g(3) * t86 - t102 * t87, -g(3) * t87 + t102 * t86, -t79, -g(1) * (-t96 * t100 + t99 * t85) - g(2) * (t99 * t100 + t96 * t85) - g(3) * t114, 0, 0, 0, 0, 0, 0, -t71, t70, -t79, -g(1) * t104 - g(2) * t107 - g(3) * t105, 0, 0, 0, 0, 0, 0, -t79, t71, -t70, -g(1) * (t104 + t117) - g(2) * t103 - g(3) * t101, 0, 0, 0, 0, 0, 0, -g(1) * (t83 * t109 + t110) - g(2) * (t83 * t111 - t108) + t94 * t115, -g(1) * (t83 * t108 - t111) - g(2) * (t83 * t110 + t109) + t97 * t115, -t71, -g(1) * (pkin(10) * t112 + t73 + (pkin(5) - t92) * t96 + t117) - g(2) * (-t99 * pkin(5) + pkin(10) * t113 + t103) - g(3) * (t83 * pkin(10) + t101);];
U_reg  = t1;
