% Calculate inertial parameters regressor of potential energy for
% S6RRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:15:39
% EndTime: 2019-03-10 01:15:40
% DurationCPUTime: 0.13s
% Computational Cost: add. (215->71), mult. (215->89), div. (0->0), fcn. (215->10), ass. (0->43)
t121 = g(3) * pkin(6);
t94 = qJ(2) + qJ(3);
t88 = sin(t94);
t120 = g(3) * t88;
t96 = sin(qJ(2));
t119 = t96 * pkin(2) + pkin(6);
t93 = qJ(4) + qJ(5);
t87 = sin(t93);
t97 = sin(qJ(1));
t118 = t97 * t87;
t89 = cos(t93);
t117 = t97 * t89;
t95 = sin(qJ(4));
t116 = t97 * t95;
t98 = cos(qJ(4));
t115 = t97 * t98;
t100 = cos(qJ(1));
t102 = -pkin(8) - pkin(7);
t99 = cos(qJ(2));
t85 = t99 * pkin(2) + pkin(1);
t114 = t100 * t102 + t97 * t85;
t90 = cos(t94);
t113 = t100 * t90;
t112 = t100 * t95;
t111 = t100 * t98;
t101 = -pkin(10) - pkin(9);
t110 = t101 * t88;
t84 = t98 * pkin(4) + pkin(3);
t109 = t90 * t101 + t88 * t84 + t119;
t108 = t100 * t85 - t97 * t102;
t107 = pkin(3) * t90 + pkin(9) * t88;
t106 = g(1) * t100 + g(2) * t97;
t70 = t100 * t89 + t90 * t118;
t72 = t87 * t113 - t117;
t105 = g(1) * t72 + g(2) * t70 + t87 * t120;
t104 = pkin(4) * t116 - t100 * t110 + t84 * t113 + t108;
t103 = -pkin(4) * t112 + t114 + (t84 * t90 - t110) * t97;
t77 = g(1) * t97 - g(2) * t100;
t73 = t89 * t113 + t118;
t71 = -t100 * t87 + t90 * t117;
t69 = -g(3) * t90 + t106 * t88;
t68 = -g(1) * t73 - g(2) * t71 - t89 * t120;
t1 = [0, 0, 0, 0, 0, 0, -t106, t77, -g(3), -t121, 0, 0, 0, 0, 0, 0, -g(3) * t96 - t106 * t99, -g(3) * t99 + t106 * t96, -t77, -g(1) * (t100 * pkin(1) + t97 * pkin(7)) - g(2) * (t97 * pkin(1) - t100 * pkin(7)) - t121, 0, 0, 0, 0, 0, 0, -t106 * t90 - t120, t69, -t77, -g(1) * t108 - g(2) * t114 - g(3) * t119, 0, 0, 0, 0, 0, 0, -g(1) * (t90 * t111 + t116) - g(2) * (t90 * t115 - t112) - t98 * t120, -g(1) * (-t90 * t112 + t115) - g(2) * (-t90 * t116 - t111) + t95 * t120, -t69, -g(1) * (t107 * t100 + t108) - g(2) * (t107 * t97 + t114) - g(3) * (t88 * pkin(3) - t90 * pkin(9) + t119) 0, 0, 0, 0, 0, 0, t68, t105, -t69, -g(1) * t104 - g(2) * t103 - g(3) * t109, 0, 0, 0, 0, 0, 0, t68, -t69, -t105, -g(1) * (t73 * pkin(5) + t72 * qJ(6) + t104) - g(2) * (t71 * pkin(5) + t70 * qJ(6) + t103) - g(3) * ((pkin(5) * t89 + qJ(6) * t87) * t88 + t109);];
U_reg  = t1;
