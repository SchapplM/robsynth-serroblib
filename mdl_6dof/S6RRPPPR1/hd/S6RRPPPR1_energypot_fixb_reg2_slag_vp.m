% Calculate inertial parameters regressor of potential energy for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:08:17
% EndTime: 2019-03-09 08:08:17
% DurationCPUTime: 0.16s
% Computational Cost: add. (215->70), mult. (262->92), div. (0->0), fcn. (278->10), ass. (0->44)
t100 = cos(pkin(10));
t98 = qJ(2) + pkin(9);
t94 = sin(t98);
t123 = t100 * t94;
t99 = sin(pkin(10));
t125 = t94 * t99;
t130 = pkin(4) * t123 + qJ(5) * t125;
t129 = g(3) * pkin(6);
t95 = cos(t98);
t128 = pkin(3) * t95;
t127 = g(3) * t94;
t103 = sin(qJ(2));
t126 = t103 * pkin(2) + pkin(6);
t101 = -pkin(7) - qJ(3);
t104 = sin(qJ(1));
t107 = cos(qJ(1));
t106 = cos(qJ(2));
t93 = t106 * pkin(2) + pkin(1);
t124 = t107 * t101 + t104 * t93;
t122 = t104 * t94;
t121 = t104 * t99;
t120 = t107 * t94;
t119 = t107 * t99;
t118 = t104 * t100;
t117 = t107 * t100;
t116 = t94 * pkin(3) + t126;
t115 = qJ(4) * t122 + t104 * t128 + t124;
t114 = -t104 * t101 + t107 * t93;
t113 = g(1) * t107 + g(2) * t104;
t112 = -t95 * qJ(4) + t116;
t111 = qJ(4) * t120 + t107 * t128 + t114;
t76 = t95 * t121 + t117;
t77 = t95 * t118 - t119;
t110 = t77 * pkin(4) + t76 * qJ(5) + t115;
t78 = t95 * t119 - t118;
t109 = g(1) * t78 + g(2) * t76 + g(3) * t125;
t79 = t95 * t117 + t121;
t108 = t79 * pkin(4) + t78 * qJ(5) + t111;
t105 = cos(qJ(6));
t102 = sin(qJ(6));
t84 = g(1) * t104 - g(2) * t107;
t73 = -g(3) * t95 + t113 * t94;
t72 = -g(1) * t79 - g(2) * t77 - g(3) * t123;
t1 = [0, 0, 0, 0, 0, 0, -t113, t84, -g(3), -t129, 0, 0, 0, 0, 0, 0, -g(3) * t103 - t113 * t106, -g(3) * t106 + t113 * t103, -t84, -g(1) * (t107 * pkin(1) + t104 * pkin(7)) - g(2) * (t104 * pkin(1) - t107 * pkin(7)) - t129, 0, 0, 0, 0, 0, 0, -t113 * t95 - t127, t73, -t84, -g(1) * t114 - g(2) * t124 - g(3) * t126, 0, 0, 0, 0, 0, 0, t72, t109, -t73, -g(1) * t111 - g(2) * t115 - g(3) * t112, 0, 0, 0, 0, 0, 0, t72, -t73, -t109, -g(1) * t108 - g(2) * t110 - g(3) * (t112 + t130) 0, 0, 0, 0, 0, 0, -g(1) * (t78 * t102 + t79 * t105) - g(2) * (t76 * t102 + t77 * t105) - (t100 * t105 + t102 * t99) * t127, -g(1) * (-t79 * t102 + t78 * t105) - g(2) * (-t77 * t102 + t76 * t105) - (-t100 * t102 + t105 * t99) * t127, t73, -g(1) * (t79 * pkin(5) - pkin(8) * t120 + t108) - g(2) * (t77 * pkin(5) - pkin(8) * t122 + t110) - g(3) * (pkin(5) * t123 + (pkin(8) - qJ(4)) * t95 + t116 + t130);];
U_reg  = t1;
