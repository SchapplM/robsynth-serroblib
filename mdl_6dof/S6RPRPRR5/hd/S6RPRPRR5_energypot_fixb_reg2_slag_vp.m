% Calculate inertial parameters regressor of potential energy for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:49:34
% EndTime: 2019-03-09 03:49:34
% DurationCPUTime: 0.18s
% Computational Cost: add. (207->65), mult. (235->85), div. (0->0), fcn. (244->10), ass. (0->40)
t105 = cos(qJ(1));
t96 = pkin(10) + qJ(3);
t92 = cos(t96);
t116 = t105 * t92;
t91 = sin(t96);
t118 = qJ(4) * t91;
t124 = pkin(3) * t116 + t105 * t118;
t101 = sin(qJ(5));
t104 = cos(qJ(5));
t123 = t92 * t101 - t91 * t104;
t122 = g(3) * pkin(6);
t121 = g(3) * t123;
t97 = sin(pkin(10));
t120 = t97 * pkin(2) + pkin(6);
t102 = sin(qJ(1));
t98 = cos(pkin(10));
t89 = t98 * pkin(2) + pkin(1);
t99 = -pkin(7) - qJ(2);
t119 = t102 * t89 + t105 * t99;
t117 = t102 * t92;
t81 = t105 * t89;
t113 = -t102 * t99 + t81;
t112 = pkin(3) * t117 + t102 * t118 + t119;
t111 = g(1) * t105 + g(2) * t102;
t74 = t91 * t101 + t92 * t104;
t110 = t91 * pkin(3) - t92 * qJ(4) + t120;
t109 = pkin(4) * t117 + t105 * pkin(8) + t112;
t108 = t91 * pkin(4) + t110;
t70 = t123 * t102;
t72 = t123 * t105;
t107 = g(1) * t72 + g(2) * t70 + g(3) * t74;
t106 = t81 + pkin(4) * t116 + (-pkin(8) - t99) * t102 + t124;
t103 = cos(qJ(6));
t100 = sin(qJ(6));
t82 = g(1) * t102 - g(2) * t105;
t73 = t74 * t105;
t71 = t74 * t102;
t69 = -g(3) * t91 - t111 * t92;
t68 = -g(3) * t92 + t111 * t91;
t1 = [0, 0, 0, 0, 0, 0, -t111, t82, -g(3), -t122, 0, 0, 0, 0, 0, 0, -g(3) * t97 - t111 * t98, -g(3) * t98 + t111 * t97, -t82, -g(1) * (t105 * pkin(1) + t102 * qJ(2)) - g(2) * (t102 * pkin(1) - t105 * qJ(2)) - t122, 0, 0, 0, 0, 0, 0, t69, t68, -t82, -g(1) * t113 - g(2) * t119 - g(3) * t120, 0, 0, 0, 0, 0, 0, t69, -t82, -t68, -g(1) * (t113 + t124) - g(2) * t112 - g(3) * t110, 0, 0, 0, 0, 0, 0, -g(1) * t73 - g(2) * t71 + t121, t107, t82, -g(1) * t106 - g(2) * t109 - g(3) * t108, 0, 0, 0, 0, 0, 0, -g(1) * (-t102 * t100 + t73 * t103) - g(2) * (t105 * t100 + t71 * t103) + t103 * t121, -g(1) * (-t73 * t100 - t102 * t103) - g(2) * (-t71 * t100 + t105 * t103) - t100 * t121, -t107, -g(1) * (t73 * pkin(5) + t72 * pkin(9) + t106) - g(2) * (t71 * pkin(5) + t70 * pkin(9) + t109) - g(3) * (-pkin(5) * t123 + t74 * pkin(9) + t108);];
U_reg  = t1;
