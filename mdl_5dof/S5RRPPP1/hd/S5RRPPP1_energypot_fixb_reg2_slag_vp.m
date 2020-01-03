% Calculate inertial parameters regressor of potential energy for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:29
% EndTime: 2019-12-31 19:24:29
% DurationCPUTime: 0.13s
% Computational Cost: add. (148->62), mult. (356->82), div. (0->0), fcn. (406->8), ass. (0->46)
t89 = cos(pkin(5));
t91 = sin(qJ(1));
t111 = t91 * t89;
t87 = sin(pkin(5));
t93 = cos(qJ(1));
t115 = t87 * t93;
t90 = sin(qJ(2));
t120 = t90 * t111 + t115;
t119 = g(3) * pkin(6);
t118 = t90 * pkin(2) + pkin(6);
t117 = t87 * t91;
t92 = cos(qJ(2));
t116 = t87 * t92;
t114 = t89 * t92;
t86 = sin(pkin(8));
t113 = t90 * t86;
t88 = cos(pkin(8));
t112 = t90 * t88;
t110 = t91 * t92;
t109 = t92 * t93;
t108 = t93 * t89;
t107 = t93 * pkin(1) + t91 * pkin(7);
t106 = qJ(3) * t87;
t105 = qJ(3) * t89;
t103 = t90 * t106;
t102 = t92 * t106;
t101 = pkin(2) * t109 + t93 * t103 + t91 * t105 + t107;
t100 = g(1) * t93 + g(2) * t91;
t67 = -t88 * t114 + t113;
t68 = t86 * t114 + t112;
t99 = t68 * pkin(3) + t67 * qJ(4) + t118;
t62 = t86 * t110 + t120 * t88;
t64 = -t88 * t117 + (t89 * t112 + t86 * t92) * t93;
t98 = g(1) * t64 + g(2) * t62 + g(3) * t67;
t63 = t88 * t110 - t120 * t86;
t65 = -t108 * t113 + t88 * t109 + t86 * t117;
t97 = g(1) * t65 + g(2) * t63 + g(3) * t68;
t96 = t65 * pkin(3) + t64 * qJ(4) + t101;
t84 = t91 * pkin(1);
t95 = t91 * t103 + pkin(2) * t110 + t84 + (-pkin(7) - t105) * t93;
t94 = t63 * pkin(3) + t62 * qJ(4) + t95;
t78 = g(1) * t91 - g(2) * t93;
t70 = t90 * t115 + t111;
t69 = t90 * t117 - t108;
t59 = -g(1) * t70 - g(2) * t69 + g(3) * t116;
t1 = [0, 0, 0, 0, 0, 0, -t100, t78, -g(3), -t119, 0, 0, 0, 0, 0, 0, -g(3) * t90 - t100 * t92, -g(3) * t92 + t100 * t90, -t78, -g(1) * t107 - g(2) * (-t93 * pkin(7) + t84) - t119, 0, 0, 0, 0, 0, 0, -t97, t98, t59, -g(1) * t101 - g(2) * t95 - g(3) * (-t102 + t118), 0, 0, 0, 0, 0, 0, t59, t97, -t98, -g(1) * t96 - g(2) * t94 - g(3) * (t99 - t102), 0, 0, 0, 0, 0, 0, t59, -t98, -t97, -g(1) * (t70 * pkin(4) + t65 * qJ(5) + t96) - g(2) * (t69 * pkin(4) + t63 * qJ(5) + t94) - g(3) * (t68 * qJ(5) + (-pkin(4) - qJ(3)) * t116 + t99);];
U_reg = t1;
