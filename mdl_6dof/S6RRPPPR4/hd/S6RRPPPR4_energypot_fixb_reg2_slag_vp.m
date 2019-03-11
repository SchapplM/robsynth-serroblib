% Calculate inertial parameters regressor of potential energy for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:19:43
% EndTime: 2019-03-09 08:19:43
% DurationCPUTime: 0.21s
% Computational Cost: add. (141->74), mult. (281->88), div. (0->0), fcn. (297->8), ass. (0->42)
t95 = sin(qJ(2));
t98 = cos(qJ(2));
t124 = pkin(2) * t98 + qJ(3) * t95;
t96 = sin(qJ(1));
t123 = t124 * t96;
t93 = cos(pkin(9));
t114 = t96 * t93;
t92 = sin(pkin(9));
t99 = cos(qJ(1));
t74 = t95 * t114 + t92 * t99;
t117 = t92 * t96;
t75 = t95 * t117 - t93 * t99;
t122 = t75 * pkin(4) - t74 * qJ(5);
t121 = g(3) * pkin(6);
t119 = g(3) * t98;
t118 = t95 * pkin(2) + pkin(6);
t116 = t93 * t98;
t115 = t95 * t99;
t113 = t99 * pkin(1) + t96 * pkin(7);
t111 = qJ(4) * t98;
t89 = t96 * pkin(1);
t109 = -pkin(7) * t99 + t89;
t85 = t95 * qJ(4);
t108 = qJ(5) * t116 + t118 + t85;
t107 = t124 * t99 + t113;
t106 = -qJ(3) * t98 + t118;
t105 = g(1) * t99 + g(2) * t96;
t104 = t109 + t123;
t103 = t96 * pkin(3) + t99 * t111 + t107;
t80 = t96 * t111;
t102 = t80 + t89 + (-pkin(3) - pkin(7)) * t99 + t123;
t72 = -t93 * t115 + t117;
t101 = g(1) * t72 - g(2) * t74 + g(3) * t116;
t73 = t92 * t115 + t114;
t100 = t73 * pkin(4) + t72 * qJ(5) + t103;
t97 = cos(qJ(6));
t94 = sin(qJ(6));
t76 = g(1) * t96 - g(2) * t99;
t69 = g(3) * t95 + t105 * t98;
t68 = t105 * t95 - t119;
t67 = -g(1) * t73 - g(2) * t75 + t92 * t119;
t1 = [0, 0, 0, 0, 0, 0, -t105, t76, -g(3), -t121, 0, 0, 0, 0, 0, 0, -t69, t68, -t76, -g(1) * t113 - g(2) * t109 - t121, 0, 0, 0, 0, 0, 0, -t76, t69, -t68, -g(1) * t107 - g(2) * t104 - g(3) * t106, 0, 0, 0, 0, 0, 0, t67, t101, -t69, -g(1) * t103 - g(2) * t102 - g(3) * (t106 + t85) 0, 0, 0, 0, 0, 0, t67, -t69, -t101, -g(1) * t100 - g(2) * (t102 + t122) - g(3) * ((-pkin(4) * t92 - qJ(3)) * t98 + t108) 0, 0, 0, 0, 0, 0, -g(1) * (t72 * t94 + t73 * t97) - g(2) * (-t74 * t94 + t75 * t97) - (-t92 * t97 + t93 * t94) * t119, -g(1) * (t72 * t97 - t73 * t94) - g(2) * (-t74 * t97 - t75 * t94) - (t92 * t94 + t93 * t97) * t119, t69, -g(1) * (t73 * pkin(5) + t100) - g(2) * (-pkin(3) * t99 + t75 * pkin(5) + t104 + t122 + t80) - g(3) * (-pkin(8) * t95 + t108) + (-g(3) * (-qJ(3) + (-pkin(4) - pkin(5)) * t92) + t105 * pkin(8)) * t98;];
U_reg  = t1;
