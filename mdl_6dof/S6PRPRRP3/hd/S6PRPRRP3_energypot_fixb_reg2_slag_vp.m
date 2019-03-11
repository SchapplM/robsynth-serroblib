% Calculate inertial parameters regressor of potential energy for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:07:44
% EndTime: 2019-03-08 20:07:44
% DurationCPUTime: 0.24s
% Computational Cost: add. (319->95), mult. (533->135), div. (0->0), fcn. (633->12), ass. (0->51)
t97 = sin(pkin(10));
t98 = sin(pkin(6));
t126 = t97 * t98;
t100 = cos(pkin(10));
t125 = t100 * pkin(1) + pkin(7) * t126;
t124 = t100 * t98;
t101 = cos(pkin(6));
t96 = sin(pkin(11));
t123 = t101 * t96;
t104 = sin(qJ(5));
t105 = sin(qJ(2));
t107 = cos(qJ(2));
t116 = t101 * t107;
t76 = -t100 * t116 + t105 * t97;
t122 = t104 * t76;
t78 = t100 * t105 + t97 * t116;
t121 = t104 * t78;
t120 = t105 * t98;
t119 = t107 * t98;
t118 = t101 * pkin(7) + qJ(1);
t117 = t101 * t105;
t115 = t96 * t126;
t114 = t104 * t119;
t92 = t97 * pkin(1);
t113 = -pkin(7) * t124 + t92;
t103 = -pkin(8) - qJ(3);
t79 = t100 * t107 - t97 * t117;
t99 = cos(pkin(11));
t88 = pkin(3) * t99 + pkin(2);
t112 = pkin(3) * t115 - t78 * t103 + t79 * t88 + t125;
t111 = pkin(3) * t123 + t103 * t119 + t88 * t120 + t118;
t110 = g(1) * t97 - g(2) * t100;
t77 = t100 * t117 + t107 * t97;
t95 = pkin(11) + qJ(4);
t90 = sin(t95);
t91 = cos(t95);
t66 = t91 * t124 + t77 * t90;
t68 = -t91 * t126 + t79 * t90;
t72 = -t101 * t91 + t90 * t120;
t109 = g(1) * t68 + g(2) * t66 + g(3) * t72;
t65 = -g(1) * t78 - g(2) * t76 + g(3) * t119;
t108 = t77 * t88 - t76 * t103 + t92 + (-pkin(3) * t96 - pkin(7)) * t124;
t106 = cos(qJ(5));
t102 = -qJ(6) - pkin(9);
t89 = pkin(5) * t106 + pkin(4);
t73 = t101 * t90 + t91 * t120;
t69 = t90 * t126 + t79 * t91;
t67 = -t90 * t124 + t77 * t91;
t63 = -g(1) * (t106 * t69 + t121) - g(2) * (t106 * t67 + t122) - g(3) * (t73 * t106 - t114);
t62 = -g(1) * (-t104 * t69 + t106 * t78) - g(2) * (-t104 * t67 + t106 * t76) - g(3) * (-t73 * t104 - t106 * t119);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t100 - g(2) * t97, t110, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t79 - g(2) * t77 - g(3) * t120, -t65, -g(3) * t101 - t110 * t98, -g(1) * t125 - g(2) * t113 - g(3) * t118, 0, 0, 0, 0, 0, 0, -g(1) * (t79 * t99 + t115) - g(2) * (-t96 * t124 + t77 * t99) - g(3) * (t99 * t120 + t123) -g(1) * (t99 * t126 - t79 * t96) - g(2) * (-t99 * t124 - t77 * t96) - g(3) * (t101 * t99 - t96 * t120) t65, -g(1) * (pkin(2) * t79 + qJ(3) * t78 + t125) - g(2) * (pkin(2) * t77 + qJ(3) * t76 + t113) - g(3) * ((pkin(2) * t105 - qJ(3) * t107) * t98 + t118) 0, 0, 0, 0, 0, 0, -g(1) * t69 - g(2) * t67 - g(3) * t73, t109, t65, -g(1) * t112 - g(2) * t108 - g(3) * t111, 0, 0, 0, 0, 0, 0, t63, t62, -t109, -g(1) * (pkin(4) * t69 + pkin(9) * t68 + t112) - g(2) * (pkin(4) * t67 + pkin(9) * t66 + t108) - g(3) * (pkin(4) * t73 + pkin(9) * t72 + t111) 0, 0, 0, 0, 0, 0, t63, t62, -t109, -g(1) * (pkin(5) * t121 - t102 * t68 + t69 * t89 + t112) - g(2) * (pkin(5) * t122 - t102 * t66 + t67 * t89 + t108) - g(3) * (-pkin(5) * t114 - t72 * t102 + t73 * t89 + t111);];
U_reg  = t1;
