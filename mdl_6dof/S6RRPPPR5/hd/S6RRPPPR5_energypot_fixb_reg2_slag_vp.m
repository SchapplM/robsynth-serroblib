% Calculate inertial parameters regressor of potential energy for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:23:43
% EndTime: 2019-03-09 08:23:43
% DurationCPUTime: 0.15s
% Computational Cost: add. (157->75), mult. (323->89), div. (0->0), fcn. (351->8), ass. (0->45)
t91 = cos(pkin(9));
t93 = sin(qJ(2));
t121 = t91 * t93;
t90 = sin(pkin(9));
t122 = t90 * t93;
t127 = pkin(3) * t121 + qJ(4) * t122;
t126 = g(3) * pkin(6);
t125 = g(3) * t93;
t96 = cos(qJ(2));
t124 = g(3) * t96;
t123 = t93 * pkin(2) + pkin(6);
t94 = sin(qJ(1));
t120 = t93 * t94;
t97 = cos(qJ(1));
t119 = t93 * t97;
t118 = t94 * t96;
t117 = t97 * t90;
t116 = t97 * t91;
t115 = -pkin(4) - qJ(3);
t114 = pkin(5) + qJ(4);
t113 = t97 * pkin(1) + t94 * pkin(7);
t112 = qJ(3) * t93;
t69 = t90 * t118 + t116;
t111 = t69 * qJ(4);
t71 = t96 * t117 - t94 * t91;
t110 = t71 * qJ(4);
t109 = t94 * pkin(1) - t97 * pkin(7);
t108 = t113 + (pkin(2) * t96 + t112) * t97;
t107 = -t96 * qJ(3) + t123;
t106 = qJ(5) * t121 + t123 + t127;
t72 = t96 * t116 + t94 * t90;
t105 = t72 * pkin(3) + t108;
t104 = g(1) * t97 + g(2) * t94;
t103 = pkin(2) * t118 + t94 * t112 + t109;
t70 = t91 * t118 - t117;
t102 = t70 * pkin(3) + t103;
t101 = g(1) * t71 + g(2) * t69 + g(3) * t122;
t100 = g(1) * t72 + g(2) * t70 + g(3) * t121;
t99 = pkin(4) * t119 + t72 * qJ(5) + t105;
t98 = pkin(4) * t120 + t70 * qJ(5) + t102;
t95 = cos(qJ(6));
t92 = sin(qJ(6));
t73 = g(1) * t94 - g(2) * t97;
t66 = t104 * t93 - t124;
t1 = [0, 0, 0, 0, 0, 0, -t104, t73, -g(3), -t126, 0, 0, 0, 0, 0, 0, -t104 * t96 - t125, t66, -t73, -g(1) * t113 - g(2) * t109 - t126, 0, 0, 0, 0, 0, 0, -t100, t101, -t66, -g(1) * t108 - g(2) * t103 - g(3) * t107, 0, 0, 0, 0, 0, 0, -t66, t100, -t101, -g(1) * (t105 + t110) - g(2) * (t102 + t111) - g(3) * (t107 + t127) 0, 0, 0, 0, 0, 0, -t101, t66, -t100, -g(1) * (t99 + t110) - g(2) * (t98 + t111) - g(3) * (t115 * t96 + t106) 0, 0, 0, 0, 0, 0, -g(1) * (t71 * t95 + t72 * t92) - g(2) * (t69 * t95 + t70 * t92) - (t90 * t95 + t91 * t92) * t125, -g(1) * (-t71 * t92 + t72 * t95) - g(2) * (-t69 * t92 + t70 * t95) - (-t90 * t92 + t91 * t95) * t125, -t66, -g(1) * (pkin(8) * t119 + t114 * t71 + t99) - g(2) * (pkin(8) * t120 + t114 * t69 + t98) - g(3) * (pkin(5) * t122 + t106) - (-pkin(8) + t115) * t124;];
U_reg  = t1;
