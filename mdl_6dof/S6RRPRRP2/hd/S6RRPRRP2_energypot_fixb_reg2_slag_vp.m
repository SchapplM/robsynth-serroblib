% Calculate inertial parameters regressor of potential energy for
% S6RRPRRP2
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
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:45:11
% EndTime: 2019-03-09 11:45:11
% DurationCPUTime: 0.12s
% Computational Cost: add. (215->60), mult. (187->74), div. (0->0), fcn. (183->10), ass. (0->40)
t97 = qJ(2) + pkin(10);
t92 = qJ(4) + t97;
t86 = sin(t92);
t87 = cos(t92);
t122 = pkin(4) * t87 + pkin(9) * t86;
t121 = g(3) * pkin(6);
t118 = g(3) * t86;
t100 = sin(qJ(2));
t117 = t100 * pkin(2) + pkin(6);
t103 = cos(qJ(2));
t89 = t103 * pkin(2) + pkin(1);
t98 = -pkin(7) - qJ(3);
t101 = sin(qJ(1));
t104 = cos(qJ(1));
t91 = cos(t97);
t76 = pkin(3) * t91 + t89;
t96 = -pkin(8) + t98;
t116 = t101 * t76 + t104 * t96;
t99 = sin(qJ(5));
t115 = t101 * t99;
t114 = t104 * t99;
t102 = cos(qJ(5));
t113 = t101 * t102;
t112 = t104 * t102;
t90 = sin(t97);
t111 = pkin(3) * t90 + t117;
t110 = -t101 * t96 + t104 * t76;
t109 = t122 * t101 + t116;
t108 = g(1) * t104 + g(2) * t101;
t107 = t122 * t104 + t110;
t106 = t86 * pkin(4) - t87 * pkin(9) + t111;
t70 = t87 * t115 + t112;
t72 = t87 * t114 - t113;
t105 = g(1) * t72 + g(2) * t70 + t99 * t118;
t81 = g(1) * t101 - g(2) * t104;
t73 = t87 * t112 + t115;
t71 = t87 * t113 - t114;
t69 = -g(3) * t87 + t108 * t86;
t68 = -g(1) * t73 - g(2) * t71 - t102 * t118;
t1 = [0, 0, 0, 0, 0, 0, -t108, t81, -g(3), -t121, 0, 0, 0, 0, 0, 0, -g(3) * t100 - t108 * t103, -g(3) * t103 + t108 * t100, -t81, -g(1) * (t104 * pkin(1) + t101 * pkin(7)) - g(2) * (t101 * pkin(1) - t104 * pkin(7)) - t121, 0, 0, 0, 0, 0, 0, -g(3) * t90 - t108 * t91, -g(3) * t91 + t108 * t90, -t81, -g(1) * (-t101 * t98 + t104 * t89) - g(2) * (t101 * t89 + t104 * t98) - g(3) * t117, 0, 0, 0, 0, 0, 0, -t108 * t87 - t118, t69, -t81, -g(1) * t110 - g(2) * t116 - g(3) * t111, 0, 0, 0, 0, 0, 0, t68, t105, -t69, -g(1) * t107 - g(2) * t109 - g(3) * t106, 0, 0, 0, 0, 0, 0, t68, -t69, -t105, -g(1) * (t73 * pkin(5) + t72 * qJ(6) + t107) - g(2) * (t71 * pkin(5) + t70 * qJ(6) + t109) - g(3) * ((pkin(5) * t102 + qJ(6) * t99) * t86 + t106);];
U_reg  = t1;
