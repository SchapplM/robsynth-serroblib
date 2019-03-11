% Calculate inertial parameters regressor of potential energy for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:13:54
% EndTime: 2019-03-09 10:13:54
% DurationCPUTime: 0.11s
% Computational Cost: add. (193->63), mult. (161->75), div. (0->0), fcn. (149->10), ass. (0->38)
t100 = cos(qJ(1));
t93 = qJ(2) + pkin(10);
t88 = qJ(4) + t93;
t84 = cos(t88);
t108 = t100 * t84;
t83 = sin(t88);
t109 = qJ(5) * t83;
t117 = pkin(4) * t108 + t100 * t109;
t116 = g(3) * pkin(6);
t115 = g(3) * t84;
t96 = sin(qJ(2));
t114 = t96 * pkin(2) + pkin(6);
t99 = cos(qJ(2));
t85 = t99 * pkin(2) + pkin(1);
t97 = sin(qJ(1));
t113 = t84 * t97;
t95 = sin(qJ(6));
t112 = t97 * t95;
t98 = cos(qJ(6));
t111 = t97 * t98;
t94 = -pkin(7) - qJ(3);
t87 = cos(t93);
t74 = pkin(3) * t87 + t85;
t92 = -pkin(8) + t94;
t110 = t100 * t92 + t97 * t74;
t107 = t100 * t95;
t106 = t100 * t98;
t86 = sin(t93);
t105 = pkin(3) * t86 + t114;
t73 = t100 * t74;
t104 = -t97 * t92 + t73;
t103 = pkin(4) * t113 + t97 * t109 + t110;
t102 = g(1) * t100 + g(2) * t97;
t101 = t83 * pkin(4) - t84 * qJ(5) + t105;
t79 = g(1) * t97 - g(2) * t100;
t71 = g(3) * t83 + t102 * t84;
t70 = t102 * t83 - t115;
t1 = [0, 0, 0, 0, 0, 0, -t102, t79, -g(3), -t116, 0, 0, 0, 0, 0, 0, -g(3) * t96 - t102 * t99, -g(3) * t99 + t102 * t96, -t79, -g(1) * (t100 * pkin(1) + t97 * pkin(7)) - g(2) * (t97 * pkin(1) - t100 * pkin(7)) - t116, 0, 0, 0, 0, 0, 0, -g(3) * t86 - t102 * t87, -g(3) * t87 + t102 * t86, -t79, -g(1) * (t100 * t85 - t97 * t94) - g(2) * (t100 * t94 + t97 * t85) - g(3) * t114, 0, 0, 0, 0, 0, 0, -t71, t70, -t79, -g(1) * t104 - g(2) * t110 - g(3) * t105, 0, 0, 0, 0, 0, 0, -t79, t71, -t70, -g(1) * (t104 + t117) - g(2) * t103 - g(3) * t101, 0, 0, 0, 0, 0, 0, -g(1) * (t83 * t107 + t111) - g(2) * (t83 * t112 - t106) + t95 * t115, -g(1) * (t83 * t106 - t112) - g(2) * (t83 * t111 + t107) + t98 * t115, -t71, -g(1) * (pkin(9) * t108 + t73 + (pkin(5) - t92) * t97 + t117) - g(2) * (-t100 * pkin(5) + pkin(9) * t113 + t103) - g(3) * (t83 * pkin(9) + t101);];
U_reg  = t1;
