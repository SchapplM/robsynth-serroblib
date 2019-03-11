% Calculate inertial parameters regressor of potential energy for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:05:15
% EndTime: 2019-03-09 06:05:15
% DurationCPUTime: 0.16s
% Computational Cost: add. (227->68), mult. (201->87), div. (0->0), fcn. (201->10), ass. (0->41)
t90 = cos(qJ(4));
t75 = t90 * pkin(4) + pkin(3);
t88 = sin(qJ(3));
t91 = cos(qJ(3));
t93 = -pkin(9) - pkin(8);
t113 = t75 * t91 - t88 * t93;
t86 = qJ(2) + pkin(6);
t112 = g(3) * t86;
t111 = g(3) * t88;
t84 = qJ(1) + pkin(10);
t77 = sin(t84);
t87 = sin(qJ(4));
t109 = t77 * t87;
t85 = qJ(4) + qJ(5);
t79 = sin(t85);
t108 = t79 * t91;
t80 = cos(t85);
t107 = t80 * t91;
t106 = t87 * t91;
t104 = t90 * t91;
t89 = sin(qJ(1));
t103 = t89 * pkin(1) + t77 * pkin(2);
t78 = cos(t84);
t92 = cos(qJ(1));
t102 = t92 * pkin(1) + t78 * pkin(2) + t77 * pkin(7);
t101 = t88 * t75 + t91 * t93 + t86;
t100 = -t78 * pkin(7) + t103;
t99 = pkin(3) * t91 + pkin(8) * t88;
t98 = g(1) * t78 + g(2) * t77;
t97 = -g(1) * t92 - g(2) * t89;
t61 = t77 * t108 + t78 * t80;
t63 = t78 * t108 - t77 * t80;
t96 = g(1) * t63 + g(2) * t61 + t79 * t111;
t95 = pkin(4) * t109 + t113 * t78 + t102;
t94 = (-pkin(4) * t87 - pkin(7)) * t78 + t103 + t113 * t77;
t68 = g(1) * t77 - g(2) * t78;
t65 = -g(3) * t91 + t98 * t88;
t64 = t78 * t107 + t77 * t79;
t62 = t77 * t107 - t78 * t79;
t60 = -g(1) * t64 - g(2) * t62 - t80 * t111;
t1 = [0, 0, 0, 0, 0, 0, t97, g(1) * t89 - g(2) * t92, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t98, t68, -g(3), t97 * pkin(1) - t112, 0, 0, 0, 0, 0, 0, -t98 * t91 - t111, t65, -t68, -g(1) * t102 - g(2) * t100 - t112, 0, 0, 0, 0, 0, 0, -g(1) * (t78 * t104 + t109) - g(2) * (t77 * t104 - t78 * t87) - t90 * t111, -g(1) * (-t78 * t106 + t77 * t90) - g(2) * (-t77 * t106 - t78 * t90) + t87 * t111, -t65, -g(1) * (t99 * t78 + t102) - g(2) * (t99 * t77 + t100) - g(3) * (t88 * pkin(3) - t91 * pkin(8) + t86) 0, 0, 0, 0, 0, 0, t60, t96, -t65, -g(1) * t95 - g(2) * t94 - g(3) * t101, 0, 0, 0, 0, 0, 0, t60, -t65, -t96, -g(1) * (t64 * pkin(5) + t63 * qJ(6) + t95) - g(2) * (t62 * pkin(5) + t61 * qJ(6) + t94) - g(3) * ((pkin(5) * t80 + qJ(6) * t79) * t88 + t101);];
U_reg  = t1;
