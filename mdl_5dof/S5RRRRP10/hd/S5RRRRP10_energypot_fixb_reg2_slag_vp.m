% Calculate inertial parameters regressor of potential energy for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:12:20
% EndTime: 2019-12-31 22:12:20
% DurationCPUTime: 0.17s
% Computational Cost: add. (187->71), mult. (424->102), div. (0->0), fcn. (511->10), ass. (0->47)
t86 = cos(pkin(5));
t117 = t86 * pkin(7) + pkin(6);
t85 = sin(pkin(5));
t90 = sin(qJ(2));
t116 = t85 * t90;
t91 = sin(qJ(1));
t115 = t85 * t91;
t93 = cos(qJ(3));
t114 = t85 * t93;
t94 = cos(qJ(2));
t113 = t85 * t94;
t95 = cos(qJ(1));
t112 = t85 * t95;
t111 = t91 * t90;
t110 = t91 * t94;
t109 = t95 * t90;
t108 = t95 * t94;
t107 = t95 * pkin(1) + pkin(7) * t115;
t106 = pkin(2) * t116 + t117;
t75 = -t86 * t111 + t108;
t105 = t75 * pkin(2) + t107;
t88 = sin(qJ(4));
t104 = pkin(4) * t88 + pkin(8);
t103 = t91 * pkin(1) - pkin(7) * t112;
t102 = g(1) * t91 - g(2) * t95;
t73 = t86 * t109 + t110;
t101 = t73 * pkin(2) + t103;
t74 = t86 * t110 + t109;
t100 = t74 * pkin(8) + t105;
t99 = -pkin(8) * t113 + t106;
t89 = sin(qJ(3));
t64 = t93 * t112 + t73 * t89;
t66 = -t91 * t114 + t75 * t89;
t70 = t89 * t116 - t86 * t93;
t98 = g(1) * t66 + g(2) * t64 + g(3) * t70;
t72 = -t86 * t108 + t111;
t97 = t72 * pkin(8) + t101;
t96 = -g(1) * t74 - g(2) * t72 + g(3) * t113;
t92 = cos(qJ(4));
t87 = -qJ(5) - pkin(9);
t81 = t92 * pkin(4) + pkin(3);
t71 = t90 * t114 + t86 * t89;
t67 = t89 * t115 + t75 * t93;
t65 = -t89 * t112 + t73 * t93;
t62 = -g(1) * (t67 * t92 + t74 * t88) - g(2) * (t65 * t92 + t72 * t88) - g(3) * (-t88 * t113 + t71 * t92);
t61 = -g(1) * (-t67 * t88 + t74 * t92) - g(2) * (-t65 * t88 + t72 * t92) - g(3) * (-t92 * t113 - t71 * t88);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t91, t102, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -g(1) * t75 - g(2) * t73 - g(3) * t116, -t96, -g(3) * t86 - t102 * t85, -g(1) * t107 - g(2) * t103 - g(3) * t117, 0, 0, 0, 0, 0, 0, -g(1) * t67 - g(2) * t65 - g(3) * t71, t98, t96, -g(1) * t100 - g(2) * t97 - g(3) * t99, 0, 0, 0, 0, 0, 0, t62, t61, -t98, -g(1) * (t67 * pkin(3) + t66 * pkin(9) + t100) - g(2) * (t65 * pkin(3) + t64 * pkin(9) + t97) - g(3) * (t71 * pkin(3) + t70 * pkin(9) + t99), 0, 0, 0, 0, 0, 0, t62, t61, -t98, -g(1) * (t104 * t74 - t66 * t87 + t67 * t81 + t105) - g(2) * (t104 * t72 - t64 * t87 + t65 * t81 + t101) - g(3) * (-t104 * t113 - t70 * t87 + t71 * t81 + t106);];
U_reg = t1;
