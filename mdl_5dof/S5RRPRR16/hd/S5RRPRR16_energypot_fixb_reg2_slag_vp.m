% Calculate inertial parameters regressor of potential energy for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR16_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:47:25
% EndTime: 2019-12-31 20:47:25
% DurationCPUTime: 0.16s
% Computational Cost: add. (155->68), mult. (347->96), div. (0->0), fcn. (405->10), ass. (0->44)
t85 = cos(pkin(5));
t114 = t85 * pkin(7) + pkin(6);
t84 = sin(pkin(5));
t88 = sin(qJ(2));
t113 = t84 * t88;
t89 = sin(qJ(1));
t112 = t84 * t89;
t92 = cos(qJ(2));
t111 = t84 * t92;
t93 = cos(qJ(1));
t110 = t84 * t93;
t109 = t89 * t88;
t108 = t89 * t92;
t107 = t93 * t88;
t106 = t93 * t92;
t105 = t93 * pkin(1) + pkin(7) * t112;
t104 = pkin(7) * t110;
t70 = -t85 * t106 + t109;
t71 = t85 * t107 + t108;
t82 = t89 * pkin(1);
t103 = t71 * pkin(2) + t70 * qJ(3) + t82;
t102 = g(1) * t89 - g(2) * t93;
t72 = t85 * t108 + t107;
t73 = -t85 * t109 + t106;
t101 = t73 * pkin(2) + t72 * qJ(3) + t105;
t100 = pkin(2) * t113 - qJ(3) * t111 + t114;
t87 = sin(qJ(4));
t91 = cos(qJ(4));
t60 = t87 * t112 - t72 * t91;
t62 = t87 * t110 + t70 * t91;
t68 = t91 * t111 + t85 * t87;
t99 = g(1) * t60 - g(2) * t62 + g(3) * t68;
t98 = -g(1) * t72 - g(2) * t70 + g(3) * t111;
t97 = g(1) * t73 + g(2) * t71 + g(3) * t113;
t96 = t85 * pkin(3) + pkin(8) * t113 + t100;
t95 = pkin(3) * t112 + t73 * pkin(8) + t101;
t94 = t71 * pkin(8) + (-pkin(3) - pkin(7)) * t110 + t103;
t90 = cos(qJ(5));
t86 = sin(qJ(5));
t69 = -t87 * t111 + t85 * t91;
t64 = -g(3) * t85 - t102 * t84;
t63 = -t91 * t110 + t70 * t87;
t61 = t91 * t112 + t72 * t87;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t93 - g(2) * t89, t102, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t97, -t98, t64, -g(1) * t105 - g(2) * (t82 - t104) - g(3) * t114, 0, 0, 0, 0, 0, 0, t64, t97, t98, -g(1) * t101 - g(2) * (t103 - t104) - g(3) * t100, 0, 0, 0, 0, 0, 0, -g(1) * t61 - g(2) * t63 - g(3) * t69, t99, -t97, -g(1) * t95 - g(2) * t94 - g(3) * t96, 0, 0, 0, 0, 0, 0, -g(1) * (t61 * t90 + t73 * t86) - g(2) * (t63 * t90 + t71 * t86) - g(3) * (t86 * t113 + t69 * t90), -g(1) * (-t61 * t86 + t73 * t90) - g(2) * (-t63 * t86 + t71 * t90) - g(3) * (t90 * t113 - t69 * t86), -t99, -g(1) * (t61 * pkin(4) + t60 * pkin(9) + t95) - g(2) * (t63 * pkin(4) - t62 * pkin(9) + t94) - g(3) * (t69 * pkin(4) + t68 * pkin(9) + t96);];
U_reg = t1;
