% Calculate inertial parameters regressor of potential energy for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:44:40
% EndTime: 2019-03-09 04:44:41
% DurationCPUTime: 0.16s
% Computational Cost: add. (199->60), mult. (236->72), div. (0->0), fcn. (242->8), ass. (0->40)
t87 = pkin(9) + qJ(3);
t83 = sin(t87);
t84 = cos(t87);
t117 = pkin(3) * t84 + pkin(8) * t83;
t93 = cos(qJ(4));
t110 = t83 * t93;
t91 = sin(qJ(4));
t111 = t83 * t91;
t116 = pkin(4) * t110 + qJ(5) * t111;
t115 = g(3) * pkin(6);
t88 = sin(pkin(9));
t112 = t88 * pkin(2) + pkin(6);
t92 = sin(qJ(1));
t109 = t92 * t91;
t108 = t92 * t93;
t94 = cos(qJ(1));
t107 = t94 * t91;
t106 = t94 * t93;
t89 = cos(pkin(9));
t80 = t89 * pkin(2) + pkin(1);
t90 = -pkin(7) - qJ(2);
t105 = t92 * t80 + t94 * t90;
t104 = qJ(6) * t83;
t103 = t83 * pkin(3) + t112;
t102 = t94 * t80 - t92 * t90;
t101 = t117 * t92 + t105;
t100 = g(1) * t94 + g(2) * t92;
t99 = -pkin(8) * t84 + t103;
t98 = t117 * t94 + t102;
t65 = t84 * t109 + t106;
t66 = t84 * t108 - t107;
t97 = t66 * pkin(4) + qJ(5) * t65 + t101;
t67 = t84 * t107 - t108;
t96 = g(1) * t67 + g(2) * t65 + g(3) * t111;
t68 = t84 * t106 + t109;
t95 = t68 * pkin(4) + t67 * qJ(5) + t98;
t72 = g(1) * t92 - g(2) * t94;
t62 = -g(3) * t84 + t100 * t83;
t61 = -g(1) * t68 - g(2) * t66 - g(3) * t110;
t1 = [0, 0, 0, 0, 0, 0, -t100, t72, -g(3), -t115, 0, 0, 0, 0, 0, 0, -g(3) * t88 - t100 * t89, -g(3) * t89 + t100 * t88, -t72, -g(1) * (pkin(1) * t94 + t92 * qJ(2)) - g(2) * (t92 * pkin(1) - qJ(2) * t94) - t115, 0, 0, 0, 0, 0, 0, -g(3) * t83 - t100 * t84, t62, -t72, -g(1) * t102 - g(2) * t105 - g(3) * t112, 0, 0, 0, 0, 0, 0, t61, t96, -t62, -g(1) * t98 - g(2) * t101 - g(3) * t99, 0, 0, 0, 0, 0, 0, t61, -t62, -t96, -g(1) * t95 - g(2) * t97 - g(3) * (t99 + t116) 0, 0, 0, 0, 0, 0, t61, -t96, t62, -g(1) * (t68 * pkin(5) - t94 * t104 + t95) - g(2) * (pkin(5) * t66 - t92 * t104 + t97) - g(3) * (pkin(5) * t110 + (-pkin(8) + qJ(6)) * t84 + t103 + t116);];
U_reg  = t1;
