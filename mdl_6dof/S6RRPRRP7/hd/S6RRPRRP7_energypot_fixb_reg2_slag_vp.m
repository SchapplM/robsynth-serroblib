% Calculate inertial parameters regressor of potential energy for
% S6RRPRRP7
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
% Datum: 2019-03-09 12:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:20:05
% EndTime: 2019-03-09 12:20:05
% DurationCPUTime: 0.18s
% Computational Cost: add. (162->63), mult. (339->84), div. (0->0), fcn. (377->8), ass. (0->42)
t103 = sin(qJ(4));
t104 = sin(qJ(2));
t107 = cos(qJ(4));
t108 = cos(qJ(2));
t131 = t108 * t103 - t104 * t107;
t130 = g(3) * pkin(6);
t129 = g(3) * t131;
t105 = sin(qJ(1));
t109 = cos(qJ(1));
t128 = t109 * pkin(1) + t105 * pkin(7);
t126 = t104 * t109;
t125 = t105 * t108;
t123 = t108 * t109;
t122 = t105 * pkin(1) - t109 * pkin(7);
t121 = pkin(2) * t123 + qJ(3) * t126 + t128;
t120 = t104 * pkin(2) - t108 * qJ(3) + pkin(6);
t119 = g(1) * t109 + g(2) * t105;
t118 = t105 * t104 * qJ(3) + pkin(2) * t125 + t122;
t117 = t104 * pkin(3) + t120;
t83 = t104 * t103 + t108 * t107;
t78 = t131 * t105;
t80 = t103 * t123 - t107 * t126;
t116 = g(1) * t80 + g(2) * t78 + g(3) * t83;
t115 = pkin(3) * t125 + t109 * pkin(8) + t118;
t114 = pkin(3) * t123 - t105 * pkin(8) + t121;
t102 = sin(qJ(5));
t106 = cos(qJ(5));
t79 = t83 * t105;
t70 = t79 * t102 - t109 * t106;
t81 = t83 * t109;
t72 = t81 * t102 + t105 * t106;
t113 = g(1) * t72 + g(2) * t70 - t102 * t129;
t112 = -pkin(4) * t131 + t83 * pkin(9) + t117;
t111 = t79 * pkin(4) + t78 * pkin(9) + t115;
t110 = t81 * pkin(4) + t80 * pkin(9) + t114;
t85 = g(1) * t105 - g(2) * t109;
t77 = -g(3) * t104 - t119 * t108;
t76 = -g(3) * t108 + t119 * t104;
t73 = -t105 * t102 + t81 * t106;
t71 = t109 * t102 + t79 * t106;
t68 = -g(1) * t73 - g(2) * t71 + t106 * t129;
t1 = [0, 0, 0, 0, 0, 0, -t119, t85, -g(3), -t130, 0, 0, 0, 0, 0, 0, t77, t76, -t85, -g(1) * t128 - g(2) * t122 - t130, 0, 0, 0, 0, 0, 0, t77, -t85, -t76, -g(1) * t121 - g(2) * t118 - g(3) * t120, 0, 0, 0, 0, 0, 0, -g(1) * t81 - g(2) * t79 + t129, t116, t85, -g(1) * t114 - g(2) * t115 - g(3) * t117, 0, 0, 0, 0, 0, 0, t68, t113, -t116, -g(1) * t110 - g(2) * t111 - g(3) * t112, 0, 0, 0, 0, 0, 0, t68, -t116, -t113, -g(1) * (t73 * pkin(5) + t72 * qJ(6) + t110) - g(2) * (t71 * pkin(5) + t70 * qJ(6) + t111) - g(3) * (-(pkin(5) * t106 + qJ(6) * t102) * t131 + t112);];
U_reg  = t1;
