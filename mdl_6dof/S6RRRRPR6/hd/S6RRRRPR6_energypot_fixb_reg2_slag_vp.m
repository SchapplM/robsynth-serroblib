% Calculate inertial parameters regressor of potential energy for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:23:11
% EndTime: 2019-03-09 22:23:11
% DurationCPUTime: 0.21s
% Computational Cost: add. (221->99), mult. (228->128), div. (0->0), fcn. (228->12), ass. (0->41)
t120 = g(3) * pkin(6);
t107 = -pkin(9) - pkin(8);
t104 = cos(qJ(3));
t88 = t104 * pkin(3) + pkin(2);
t102 = sin(qJ(2));
t119 = g(3) * t102;
t101 = sin(qJ(3));
t93 = t101 * pkin(3);
t100 = qJ(3) + qJ(4);
t90 = sin(t100);
t79 = pkin(4) * t90 + t93;
t103 = sin(qJ(1));
t106 = cos(qJ(1));
t118 = t106 * pkin(1) + t103 * pkin(7);
t117 = t102 * t103;
t116 = t102 * t107;
t115 = t103 * t101;
t105 = cos(qJ(2));
t114 = t103 * t105;
t113 = t105 * t106;
t112 = t106 * t101;
t111 = t106 * t104;
t99 = -qJ(5) + t107;
t91 = cos(t100);
t78 = pkin(4) * t91 + t88;
t95 = t103 * pkin(1);
t110 = -t106 * pkin(7) + t95;
t89 = pkin(11) + t100;
t109 = pkin(2) * t105 + pkin(8) * t102;
t108 = g(1) * t106 + g(2) * t103;
t92 = -pkin(10) + t99;
t87 = qJ(6) + t89;
t84 = cos(t89);
t83 = sin(t89);
t82 = cos(t87);
t81 = sin(t87);
t80 = g(1) * t103 - g(2) * t106;
t77 = -g(3) * t105 + t108 * t102;
t76 = pkin(5) * t83 + t79;
t75 = pkin(5) * t84 + t78;
t1 = [0, 0, 0, 0, 0, 0, -t108, t80, -g(3), -t120, 0, 0, 0, 0, 0, 0, -t108 * t105 - t119, t77, -t80, -g(1) * t118 - g(2) * t110 - t120, 0, 0, 0, 0, 0, 0, -g(1) * (t105 * t111 + t115) - g(2) * (t104 * t114 - t112) - t104 * t119, -g(1) * (t103 * t104 - t105 * t112) - g(2) * (-t101 * t114 - t111) + t101 * t119, -t77, -g(1) * (t109 * t106 + t118) - g(2) * (t109 * t103 + t110) - g(3) * (t102 * pkin(2) - t105 * pkin(8) + pkin(6)) 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t90 + t91 * t113) - g(2) * (-t106 * t90 + t91 * t114) - t91 * t119, -g(1) * (t103 * t91 - t90 * t113) - g(2) * (-t106 * t91 - t90 * t114) + t90 * t119, -t77, -g(1) * (pkin(3) * t115 + t118) - g(2) * (-t103 * t116 + t88 * t114 + t95) - g(3) * (t102 * t88 + t105 * t107 + pkin(6)) + (-g(1) * (t105 * t88 - t116) - g(2) * (-pkin(7) - t93)) * t106, 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t83 + t84 * t113) - g(2) * (-t106 * t83 + t84 * t114) - t84 * t119, -g(1) * (t103 * t84 - t83 * t113) - g(2) * (-t106 * t84 - t83 * t114) + t83 * t119, -t77, -g(1) * (t103 * t79 + t118) - g(2) * (t78 * t114 - t99 * t117 + t95) - g(3) * (t102 * t78 + t105 * t99 + pkin(6)) + (-g(1) * (-t102 * t99 + t105 * t78) - g(2) * (-pkin(7) - t79)) * t106, 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t81 + t82 * t113) - g(2) * (-t106 * t81 + t82 * t114) - t82 * t119, -g(1) * (t103 * t82 - t81 * t113) - g(2) * (-t106 * t82 - t81 * t114) + t81 * t119, -t77, -g(1) * (t103 * t76 + t118) - g(2) * (t75 * t114 - t92 * t117 + t95) - g(3) * (t102 * t75 + t105 * t92 + pkin(6)) + (-g(1) * (-t102 * t92 + t105 * t75) - g(2) * (-pkin(7) - t76)) * t106;];
U_reg  = t1;
