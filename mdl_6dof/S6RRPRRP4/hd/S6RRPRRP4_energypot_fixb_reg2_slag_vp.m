% Calculate inertial parameters regressor of potential energy for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:55:38
% EndTime: 2019-03-09 11:55:38
% DurationCPUTime: 0.15s
% Computational Cost: add. (215->70), mult. (215->87), div. (0->0), fcn. (215->10), ass. (0->44)
t102 = -pkin(9) - pkin(8);
t99 = cos(qJ(4));
t85 = t99 * pkin(4) + pkin(3);
t93 = qJ(2) + pkin(10);
t87 = sin(t93);
t88 = cos(t93);
t124 = -t102 * t87 + t85 * t88;
t123 = g(3) * pkin(6);
t122 = g(3) * t87;
t97 = sin(qJ(2));
t121 = t97 * pkin(2) + pkin(6);
t94 = qJ(4) + qJ(5);
t89 = sin(t94);
t98 = sin(qJ(1));
t119 = t98 * t89;
t90 = cos(t94);
t118 = t98 * t90;
t96 = sin(qJ(4));
t117 = t98 * t96;
t116 = t98 * t99;
t101 = cos(qJ(1));
t100 = cos(qJ(2));
t86 = t100 * pkin(2) + pkin(1);
t95 = -pkin(7) - qJ(3);
t115 = t101 * t95 + t98 * t86;
t114 = t101 * t89;
t113 = t101 * t90;
t112 = t101 * t96;
t111 = t101 * t99;
t109 = t101 * t86 - t98 * t95;
t108 = t88 * t102 + t87 * t85 + t121;
t107 = pkin(3) * t88 + pkin(8) * t87;
t106 = g(1) * t101 + g(2) * t98;
t70 = t88 * t119 + t113;
t72 = t88 * t114 - t118;
t105 = g(1) * t72 + g(2) * t70 + t89 * t122;
t104 = pkin(4) * t117 + t124 * t101 + t109;
t103 = -pkin(4) * t112 + t124 * t98 + t115;
t77 = g(1) * t98 - g(2) * t101;
t73 = t88 * t113 + t119;
t71 = t88 * t118 - t114;
t69 = -g(3) * t88 + t106 * t87;
t68 = -g(1) * t73 - g(2) * t71 - t90 * t122;
t1 = [0, 0, 0, 0, 0, 0, -t106, t77, -g(3), -t123, 0, 0, 0, 0, 0, 0, -g(3) * t97 - t106 * t100, -g(3) * t100 + t106 * t97, -t77, -g(1) * (t101 * pkin(1) + t98 * pkin(7)) - g(2) * (t98 * pkin(1) - t101 * pkin(7)) - t123, 0, 0, 0, 0, 0, 0, -t106 * t88 - t122, t69, -t77, -g(1) * t109 - g(2) * t115 - g(3) * t121, 0, 0, 0, 0, 0, 0, -g(1) * (t88 * t111 + t117) - g(2) * (t88 * t116 - t112) - t99 * t122, -g(1) * (-t88 * t112 + t116) - g(2) * (-t88 * t117 - t111) + t96 * t122, -t69, -g(1) * (t107 * t101 + t109) - g(2) * (t107 * t98 + t115) - g(3) * (t87 * pkin(3) - t88 * pkin(8) + t121) 0, 0, 0, 0, 0, 0, t68, t105, -t69, -g(1) * t104 - g(2) * t103 - g(3) * t108, 0, 0, 0, 0, 0, 0, t68, -t69, -t105, -g(1) * (t73 * pkin(5) + t72 * qJ(6) + t104) - g(2) * (t71 * pkin(5) + t70 * qJ(6) + t103) - g(3) * ((pkin(5) * t90 + qJ(6) * t89) * t87 + t108);];
U_reg  = t1;
