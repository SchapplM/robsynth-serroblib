% Calculate inertial parameters regressor of potential energy for
% S6RRPRRP8
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
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:25:14
% EndTime: 2019-03-09 12:25:14
% DurationCPUTime: 0.17s
% Computational Cost: add. (223->85), mult. (243->108), div. (0->0), fcn. (247->10), ass. (0->41)
t119 = g(3) * pkin(6);
t98 = sin(pkin(10));
t118 = t98 * pkin(3);
t99 = cos(pkin(10));
t87 = t99 * pkin(3) + pkin(2);
t101 = sin(qJ(2));
t117 = g(3) * t101;
t100 = -pkin(8) - qJ(3);
t102 = sin(qJ(1));
t104 = cos(qJ(1));
t116 = t104 * pkin(1) + t102 * pkin(7);
t115 = t102 * t98;
t114 = t101 * t102;
t103 = cos(qJ(2));
t113 = t102 * t103;
t112 = t103 * t104;
t97 = pkin(10) + qJ(4);
t89 = cos(t97);
t80 = pkin(4) * t89 + t87;
t96 = -pkin(9) + t100;
t111 = t101 * t80 + t103 * t96 + pkin(6);
t93 = t102 * pkin(1);
t110 = -t104 * pkin(7) + t93;
t109 = g(1) * t104 + g(2) * t102;
t108 = pkin(2) * t103 + qJ(3) * t101;
t88 = sin(t97);
t81 = pkin(4) * t88 + t118;
t107 = -t104 * t101 * t96 + t102 * t81 + t80 * t112 + t116;
t90 = qJ(5) + t97;
t85 = sin(t90);
t86 = cos(t90);
t71 = t104 * t86 + t85 * t113;
t73 = -t102 * t86 + t85 * t112;
t106 = g(1) * t73 + g(2) * t71 + t85 * t117;
t105 = -t96 * t114 + t80 * t113 + t93 + (-pkin(7) - t81) * t104;
t83 = g(1) * t102 - g(2) * t104;
t75 = -g(3) * t103 + t109 * t101;
t74 = t102 * t85 + t86 * t112;
t72 = -t104 * t85 + t86 * t113;
t70 = -g(1) * t74 - g(2) * t72 - t86 * t117;
t1 = [0, 0, 0, 0, 0, 0, -t109, t83, -g(3), -t119, 0, 0, 0, 0, 0, 0, -t109 * t103 - t117, t75, -t83, -g(1) * t116 - g(2) * t110 - t119, 0, 0, 0, 0, 0, 0, -g(1) * (t99 * t112 + t115) - g(2) * (-t104 * t98 + t99 * t113) - t99 * t117, -g(1) * (t102 * t99 - t98 * t112) - g(2) * (-t104 * t99 - t98 * t113) + t98 * t117, -t75, -g(1) * (t108 * t104 + t116) - g(2) * (t108 * t102 + t110) - g(3) * (t101 * pkin(2) - t103 * qJ(3) + pkin(6)) 0, 0, 0, 0, 0, 0, -g(1) * (t102 * t88 + t89 * t112) - g(2) * (-t104 * t88 + t89 * t113) - t89 * t117, -g(1) * (t102 * t89 - t88 * t112) - g(2) * (-t104 * t89 - t88 * t113) + t88 * t117, -t75, -g(1) * (pkin(3) * t115 + t116) - g(2) * (-t100 * t114 + t87 * t113 + t93) - g(3) * (t103 * t100 + t101 * t87 + pkin(6)) + (-g(1) * (-t100 * t101 + t103 * t87) - g(2) * (-pkin(7) - t118)) * t104, 0, 0, 0, 0, 0, 0, t70, t106, -t75, -g(1) * t107 - g(2) * t105 - g(3) * t111, 0, 0, 0, 0, 0, 0, t70, -t75, -t106, -g(1) * (t74 * pkin(5) + t73 * qJ(6) + t107) - g(2) * (t72 * pkin(5) + t71 * qJ(6) + t105) - g(3) * ((pkin(5) * t86 + qJ(6) * t85) * t101 + t111);];
U_reg  = t1;
