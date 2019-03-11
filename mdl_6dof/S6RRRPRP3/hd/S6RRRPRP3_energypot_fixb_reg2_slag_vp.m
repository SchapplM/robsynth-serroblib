% Calculate inertial parameters regressor of potential energy for
% S6RRRPRP3
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
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:41:30
% EndTime: 2019-03-09 16:41:30
% DurationCPUTime: 0.13s
% Computational Cost: add. (215->71), mult. (215->89), div. (0->0), fcn. (215->10), ass. (0->43)
t121 = g(3) * pkin(6);
t94 = qJ(2) + qJ(3);
t89 = sin(t94);
t120 = g(3) * t89;
t98 = sin(qJ(2));
t119 = t98 * pkin(2) + pkin(6);
t97 = -pkin(9) - qJ(4);
t118 = t89 * t97;
t93 = pkin(10) + qJ(5);
t87 = sin(t93);
t99 = sin(qJ(1));
t117 = t99 * t87;
t88 = cos(t93);
t116 = t99 * t88;
t95 = sin(pkin(10));
t115 = t99 * t95;
t96 = cos(pkin(10));
t114 = t99 * t96;
t101 = cos(qJ(1));
t102 = -pkin(8) - pkin(7);
t100 = cos(qJ(2));
t85 = t100 * pkin(2) + pkin(1);
t113 = t101 * t102 + t99 * t85;
t90 = cos(t94);
t112 = t101 * t90;
t111 = t101 * t95;
t110 = t101 * t96;
t83 = t96 * pkin(4) + pkin(3);
t109 = t89 * t83 + t90 * t97 + t119;
t108 = t101 * t85 - t99 * t102;
t107 = g(1) * t101 + g(2) * t99;
t106 = pkin(3) * t90 + qJ(4) * t89;
t70 = t101 * t88 + t90 * t117;
t72 = t87 * t112 - t116;
t105 = g(1) * t72 + g(2) * t70 + t87 * t120;
t104 = pkin(4) * t115 - t101 * t118 + t83 * t112 + t108;
t103 = -pkin(4) * t111 + t113 + (t83 * t90 - t118) * t99;
t77 = g(1) * t99 - g(2) * t101;
t73 = t88 * t112 + t117;
t71 = -t101 * t87 + t90 * t116;
t69 = -g(3) * t90 + t107 * t89;
t68 = -g(1) * t73 - g(2) * t71 - t88 * t120;
t1 = [0, 0, 0, 0, 0, 0, -t107, t77, -g(3), -t121, 0, 0, 0, 0, 0, 0, -g(3) * t98 - t107 * t100, -g(3) * t100 + t107 * t98, -t77, -g(1) * (t101 * pkin(1) + t99 * pkin(7)) - g(2) * (t99 * pkin(1) - t101 * pkin(7)) - t121, 0, 0, 0, 0, 0, 0, -t107 * t90 - t120, t69, -t77, -g(1) * t108 - g(2) * t113 - g(3) * t119, 0, 0, 0, 0, 0, 0, -g(1) * (t90 * t110 + t115) - g(2) * (t90 * t114 - t111) - t96 * t120, -g(1) * (-t90 * t111 + t114) - g(2) * (-t90 * t115 - t110) + t95 * t120, -t69, -g(1) * (t106 * t101 + t108) - g(2) * (t106 * t99 + t113) - g(3) * (t89 * pkin(3) - t90 * qJ(4) + t119) 0, 0, 0, 0, 0, 0, t68, t105, -t69, -g(1) * t104 - g(2) * t103 - g(3) * t109, 0, 0, 0, 0, 0, 0, t68, -t69, -t105, -g(1) * (t73 * pkin(5) + t72 * qJ(6) + t104) - g(2) * (t71 * pkin(5) + t70 * qJ(6) + t103) - g(3) * ((pkin(5) * t88 + qJ(6) * t87) * t89 + t109);];
U_reg  = t1;
