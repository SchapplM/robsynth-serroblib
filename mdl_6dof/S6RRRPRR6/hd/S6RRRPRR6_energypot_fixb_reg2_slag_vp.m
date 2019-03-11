% Calculate inertial parameters regressor of potential energy for
% S6RRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:29:29
% EndTime: 2019-03-09 18:29:29
% DurationCPUTime: 0.21s
% Computational Cost: add. (221->99), mult. (228->128), div. (0->0), fcn. (228->12), ass. (0->40)
t119 = g(3) * pkin(6);
t105 = cos(qJ(3));
t88 = t105 * pkin(3) + pkin(2);
t103 = sin(qJ(2));
t118 = g(3) * t103;
t102 = sin(qJ(3));
t93 = t102 * pkin(3);
t101 = -qJ(4) - pkin(8);
t100 = qJ(3) + pkin(11);
t89 = sin(t100);
t79 = pkin(4) * t89 + t93;
t104 = sin(qJ(1));
t107 = cos(qJ(1));
t117 = t107 * pkin(1) + t104 * pkin(7);
t116 = t103 * t104;
t115 = t104 * t102;
t106 = cos(qJ(2));
t114 = t104 * t106;
t113 = t106 * t107;
t112 = t107 * t102;
t111 = t107 * t105;
t99 = -pkin(9) + t101;
t90 = cos(t100);
t78 = pkin(4) * t90 + t88;
t95 = t104 * pkin(1);
t110 = -t107 * pkin(7) + t95;
t91 = qJ(5) + t100;
t109 = pkin(2) * t106 + pkin(8) * t103;
t108 = g(1) * t107 + g(2) * t104;
t92 = -pkin(10) + t99;
t87 = qJ(6) + t91;
t86 = cos(t91);
t85 = sin(t91);
t82 = cos(t87);
t81 = sin(t87);
t80 = g(1) * t104 - g(2) * t107;
t77 = -g(3) * t106 + t108 * t103;
t76 = pkin(5) * t85 + t79;
t75 = pkin(5) * t86 + t78;
t1 = [0, 0, 0, 0, 0, 0, -t108, t80, -g(3), -t119, 0, 0, 0, 0, 0, 0, -t108 * t106 - t118, t77, -t80, -g(1) * t117 - g(2) * t110 - t119, 0, 0, 0, 0, 0, 0, -g(1) * (t106 * t111 + t115) - g(2) * (t105 * t114 - t112) - t105 * t118, -g(1) * (t104 * t105 - t106 * t112) - g(2) * (-t102 * t114 - t111) + t102 * t118, -t77, -g(1) * (t109 * t107 + t117) - g(2) * (t109 * t104 + t110) - g(3) * (t103 * pkin(2) - t106 * pkin(8) + pkin(6)) 0, 0, 0, 0, 0, 0, -g(1) * (t104 * t89 + t90 * t113) - g(2) * (-t107 * t89 + t90 * t114) - t90 * t118, -g(1) * (t104 * t90 - t89 * t113) - g(2) * (-t107 * t90 - t89 * t114) + t89 * t118, -t77, -g(1) * (pkin(3) * t115 + t117) - g(2) * (-t101 * t116 + t88 * t114 + t95) - g(3) * (t106 * t101 + t103 * t88 + pkin(6)) + (-g(1) * (-t101 * t103 + t106 * t88) - g(2) * (-pkin(7) - t93)) * t107, 0, 0, 0, 0, 0, 0, -g(1) * (t104 * t85 + t86 * t113) - g(2) * (-t107 * t85 + t86 * t114) - t86 * t118, -g(1) * (t104 * t86 - t85 * t113) - g(2) * (-t107 * t86 - t85 * t114) + t85 * t118, -t77, -g(1) * (t104 * t79 + t117) - g(2) * (t78 * t114 - t99 * t116 + t95) - g(3) * (t103 * t78 + t106 * t99 + pkin(6)) + (-g(1) * (-t103 * t99 + t106 * t78) - g(2) * (-pkin(7) - t79)) * t107, 0, 0, 0, 0, 0, 0, -g(1) * (t104 * t81 + t82 * t113) - g(2) * (-t107 * t81 + t82 * t114) - t82 * t118, -g(1) * (t104 * t82 - t81 * t113) - g(2) * (-t107 * t82 - t81 * t114) + t81 * t118, -t77, -g(1) * (t104 * t76 + t117) - g(2) * (t75 * t114 - t92 * t116 + t95) - g(3) * (t103 * t75 + t106 * t92 + pkin(6)) + (-g(1) * (-t103 * t92 + t106 * t75) - g(2) * (-pkin(7) - t76)) * t107;];
U_reg  = t1;
