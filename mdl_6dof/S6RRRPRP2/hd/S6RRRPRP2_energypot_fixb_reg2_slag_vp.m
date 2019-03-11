% Calculate inertial parameters regressor of potential energy for
% S6RRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:37:17
% EndTime: 2019-03-09 16:37:17
% DurationCPUTime: 0.12s
% Computational Cost: add. (215->60), mult. (187->74), div. (0->0), fcn. (183->10), ass. (0->40)
t97 = qJ(2) + qJ(3);
t90 = pkin(10) + t97;
t85 = sin(t90);
t86 = cos(t90);
t122 = pkin(4) * t86 + pkin(9) * t85;
t121 = g(3) * pkin(6);
t104 = -pkin(8) - pkin(7);
t118 = g(3) * t85;
t99 = sin(qJ(2));
t117 = t99 * pkin(2) + pkin(6);
t102 = cos(qJ(2));
t89 = t102 * pkin(2) + pkin(1);
t100 = sin(qJ(1));
t103 = cos(qJ(1));
t92 = cos(t97);
t76 = pkin(3) * t92 + t89;
t96 = -qJ(4) + t104;
t116 = t100 * t76 + t103 * t96;
t98 = sin(qJ(5));
t115 = t100 * t98;
t114 = t103 * t98;
t101 = cos(qJ(5));
t113 = t100 * t101;
t112 = t103 * t101;
t91 = sin(t97);
t111 = pkin(3) * t91 + t117;
t110 = -t100 * t96 + t103 * t76;
t109 = t122 * t100 + t116;
t108 = g(1) * t103 + g(2) * t100;
t107 = t122 * t103 + t110;
t106 = t85 * pkin(4) - t86 * pkin(9) + t111;
t70 = t86 * t115 + t112;
t72 = t86 * t114 - t113;
t105 = g(1) * t72 + g(2) * t70 + t98 * t118;
t81 = g(1) * t100 - g(2) * t103;
t73 = t86 * t112 + t115;
t71 = t86 * t113 - t114;
t69 = -g(3) * t86 + t108 * t85;
t68 = -g(1) * t73 - g(2) * t71 - t101 * t118;
t1 = [0, 0, 0, 0, 0, 0, -t108, t81, -g(3), -t121, 0, 0, 0, 0, 0, 0, -g(3) * t99 - t108 * t102, -g(3) * t102 + t108 * t99, -t81, -g(1) * (t103 * pkin(1) + t100 * pkin(7)) - g(2) * (t100 * pkin(1) - t103 * pkin(7)) - t121, 0, 0, 0, 0, 0, 0, -g(3) * t91 - t108 * t92, -g(3) * t92 + t108 * t91, -t81, -g(1) * (-t100 * t104 + t103 * t89) - g(2) * (t100 * t89 + t103 * t104) - g(3) * t117, 0, 0, 0, 0, 0, 0, -t108 * t86 - t118, t69, -t81, -g(1) * t110 - g(2) * t116 - g(3) * t111, 0, 0, 0, 0, 0, 0, t68, t105, -t69, -g(1) * t107 - g(2) * t109 - g(3) * t106, 0, 0, 0, 0, 0, 0, t68, -t69, -t105, -g(1) * (t73 * pkin(5) + t72 * qJ(6) + t107) - g(2) * (t71 * pkin(5) + t70 * qJ(6) + t109) - g(3) * ((pkin(5) * t101 + qJ(6) * t98) * t85 + t106);];
U_reg  = t1;
