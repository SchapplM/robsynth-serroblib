% Calculate inertial parameters regressor of potential energy for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP12_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP12_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:53:31
% EndTime: 2019-03-09 12:53:31
% DurationCPUTime: 0.15s
% Computational Cost: add. (159->79), mult. (240->92), div. (0->0), fcn. (240->8), ass. (0->47)
t91 = sin(qJ(2));
t109 = qJ(3) * t91;
t92 = sin(qJ(1));
t94 = cos(qJ(2));
t114 = t92 * t94;
t121 = pkin(2) * t114 + t92 * t109;
t120 = g(3) * pkin(6);
t119 = g(3) * t94;
t118 = t91 * pkin(2) + pkin(6);
t117 = t91 * t92;
t95 = cos(qJ(1));
t116 = t91 * t95;
t93 = cos(qJ(4));
t115 = t92 * t93;
t113 = t93 * t95;
t112 = t94 * t95;
t96 = -pkin(9) - pkin(8);
t111 = t94 * t96;
t110 = t95 * pkin(1) + t92 * pkin(7);
t90 = sin(qJ(4));
t108 = t90 * t117;
t107 = t90 * t116;
t86 = t92 * pkin(1);
t106 = t86 + t121;
t105 = -t95 * pkin(7) + t86;
t104 = -pkin(4) * t90 - qJ(3);
t103 = pkin(2) * t112 + t95 * t109 + t110;
t102 = -t91 * t96 + t118;
t101 = -t94 * qJ(3) + t118;
t100 = g(1) * t95 + g(2) * t92;
t99 = t105 + t121;
t81 = pkin(4) * t93 + pkin(3);
t98 = pkin(4) * t107 + t92 * t81 + t103;
t89 = qJ(4) + qJ(5);
t82 = sin(t89);
t83 = cos(t89);
t66 = -t83 * t116 + t82 * t92;
t68 = t83 * t117 + t82 * t95;
t97 = g(1) * t66 - g(2) * t68 + t83 * t119;
t75 = g(1) * t92 - g(2) * t95;
t73 = pkin(4) * t108;
t71 = g(3) * t91 + t100 * t94;
t70 = t100 * t91 - t119;
t69 = t82 * t117 - t83 * t95;
t67 = t82 * t116 + t83 * t92;
t65 = -g(1) * t67 - g(2) * t69 + t82 * t119;
t1 = [0, 0, 0, 0, 0, 0, -t100, t75, -g(3), -t120, 0, 0, 0, 0, 0, 0, -t71, t70, -t75, -g(1) * t110 - g(2) * t105 - t120, 0, 0, 0, 0, 0, 0, -t75, t71, -t70, -g(1) * t103 - g(2) * t99 - g(3) * t101, 0, 0, 0, 0, 0, 0, -g(1) * (t107 + t115) - g(2) * (t108 - t113) + t90 * t119, -g(1) * (t91 * t113 - t90 * t92) - g(2) * (t91 * t115 + t90 * t95) + t93 * t119, -t71, -g(1) * (pkin(3) * t92 + pkin(8) * t112 + t103) - g(2) * (pkin(8) * t114 + (-pkin(3) - pkin(7)) * t95 + t106) - g(3) * (pkin(8) * t91 + t101) 0, 0, 0, 0, 0, 0, t65, t97, -t71, -g(1) * (-t95 * t111 + t98) - g(2) * (-t92 * t111 + t73 + (-pkin(7) - t81) * t95 + t106) - g(3) * (t104 * t94 + t102) 0, 0, 0, 0, 0, 0, t65, -t71, -t97, -g(1) * (t67 * pkin(5) + t66 * qJ(6) + t98) - g(2) * (t69 * pkin(5) - t68 * qJ(6) - t95 * t81 + t73 + t99) - g(3) * t102 + (-g(3) * (-pkin(5) * t82 + qJ(6) * t83 + t104) + t100 * t96) * t94;];
U_reg  = t1;
