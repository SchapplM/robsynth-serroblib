% Calculate inertial parameters regressor of potential energy for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:13:21
% EndTime: 2019-03-09 18:13:21
% DurationCPUTime: 0.16s
% Computational Cost: add. (207->65), mult. (235->85), div. (0->0), fcn. (244->10), ass. (0->40)
t108 = cos(qJ(1));
t100 = qJ(2) + qJ(3);
t96 = cos(t100);
t120 = t108 * t96;
t95 = sin(t100);
t122 = qJ(4) * t95;
t128 = pkin(3) * t120 + t108 * t122;
t102 = sin(qJ(5));
t106 = cos(qJ(5));
t127 = t96 * t102 - t95 * t106;
t126 = g(3) * pkin(6);
t125 = g(3) * t127;
t103 = sin(qJ(2));
t124 = t103 * pkin(2) + pkin(6);
t104 = sin(qJ(1));
t109 = -pkin(8) - pkin(7);
t107 = cos(qJ(2));
t93 = t107 * pkin(2) + pkin(1);
t123 = t104 * t93 + t108 * t109;
t121 = t104 * t96;
t117 = pkin(3) * t121 + t104 * t122 + t123;
t86 = t108 * t93;
t116 = -t104 * t109 + t86;
t115 = g(1) * t108 + g(2) * t104;
t78 = t95 * t102 + t96 * t106;
t114 = t95 * pkin(3) - t96 * qJ(4) + t124;
t113 = pkin(4) * t121 + t108 * pkin(9) + t117;
t112 = t95 * pkin(4) + t114;
t74 = t127 * t104;
t76 = t127 * t108;
t111 = g(1) * t76 + g(2) * t74 + g(3) * t78;
t110 = t86 + pkin(4) * t120 + (-pkin(9) - t109) * t104 + t128;
t105 = cos(qJ(6));
t101 = sin(qJ(6));
t82 = g(1) * t104 - g(2) * t108;
t77 = t78 * t108;
t75 = t78 * t104;
t73 = -g(3) * t95 - t115 * t96;
t72 = -g(3) * t96 + t115 * t95;
t1 = [0, 0, 0, 0, 0, 0, -t115, t82, -g(3), -t126, 0, 0, 0, 0, 0, 0, -g(3) * t103 - t115 * t107, -g(3) * t107 + t115 * t103, -t82, -g(1) * (t108 * pkin(1) + t104 * pkin(7)) - g(2) * (t104 * pkin(1) - t108 * pkin(7)) - t126, 0, 0, 0, 0, 0, 0, t73, t72, -t82, -g(1) * t116 - g(2) * t123 - g(3) * t124, 0, 0, 0, 0, 0, 0, t73, -t82, -t72, -g(1) * (t116 + t128) - g(2) * t117 - g(3) * t114, 0, 0, 0, 0, 0, 0, -g(1) * t77 - g(2) * t75 + t125, t111, t82, -g(1) * t110 - g(2) * t113 - g(3) * t112, 0, 0, 0, 0, 0, 0, -g(1) * (-t104 * t101 + t77 * t105) - g(2) * (t108 * t101 + t75 * t105) + t105 * t125, -g(1) * (-t77 * t101 - t104 * t105) - g(2) * (-t75 * t101 + t108 * t105) - t101 * t125, -t111, -g(1) * (t77 * pkin(5) + t76 * pkin(10) + t110) - g(2) * (t75 * pkin(5) + t74 * pkin(10) + t113) - g(3) * (-pkin(5) * t127 + t78 * pkin(10) + t112);];
U_reg  = t1;
