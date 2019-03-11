% Calculate inertial parameters regressor of potential energy for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPPRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:23:02
% EndTime: 2019-03-08 19:23:02
% DurationCPUTime: 0.17s
% Computational Cost: add. (288->86), mult. (665->131), div. (0->0), fcn. (818->12), ass. (0->55)
t87 = sin(pkin(10));
t88 = sin(pkin(6));
t121 = t87 * t88;
t90 = cos(pkin(10));
t120 = t88 * t90;
t93 = sin(qJ(5));
t119 = t88 * t93;
t94 = sin(qJ(2));
t118 = t88 * t94;
t96 = cos(qJ(5));
t117 = t88 * t96;
t97 = cos(qJ(2));
t116 = t88 * t97;
t91 = cos(pkin(6));
t115 = t91 * t94;
t114 = t91 * t97;
t113 = t90 * pkin(1) + pkin(7) * t121;
t112 = qJ(4) * t88;
t111 = t91 * pkin(7) + qJ(1);
t110 = t87 * pkin(1) - pkin(7) * t120;
t74 = t87 * t114 + t90 * t94;
t75 = -t87 * t115 + t90 * t97;
t109 = t75 * pkin(2) + t74 * qJ(3) + t113;
t108 = pkin(2) * t118 - qJ(3) * t116 + t111;
t72 = -t90 * t114 + t87 * t94;
t73 = t90 * t115 + t87 * t97;
t86 = sin(pkin(11));
t89 = cos(pkin(11));
t57 = t72 * t86 + t73 * t89;
t49 = -t90 * t117 + t57 * t93;
t59 = t74 * t86 + t75 * t89;
t51 = t87 * t117 + t59 * t93;
t67 = (-t86 * t97 + t89 * t94) * t88;
t60 = t67 * t93 + t91 * t96;
t107 = g(1) * t51 + g(2) * t49 + g(3) * t60;
t56 = -t72 * t89 + t73 * t86;
t58 = -t74 * t89 + t75 * t86;
t66 = (t86 * t94 + t89 * t97) * t88;
t106 = g(1) * t58 + g(2) * t56 + g(3) * t66;
t105 = t73 * pkin(2) + t72 * qJ(3) + t110;
t104 = -g(1) * t74 - g(2) * t72 + g(3) * t116;
t103 = t73 * pkin(3) + t90 * t112 + t105;
t102 = t75 * pkin(3) - t87 * t112 + t109;
t101 = pkin(3) * t118 - t91 * qJ(4) + t108;
t100 = t57 * pkin(4) + t56 * pkin(8) + t103;
t99 = t59 * pkin(4) + t58 * pkin(8) + t102;
t98 = t67 * pkin(4) + t66 * pkin(8) + t101;
t95 = cos(qJ(6));
t92 = sin(qJ(6));
t63 = g(1) * t121 - g(2) * t120 + g(3) * t91;
t61 = t67 * t96 - t91 * t93;
t53 = -g(1) * t75 - g(2) * t73 - g(3) * t118;
t52 = -t87 * t119 + t59 * t96;
t50 = t90 * t119 + t57 * t96;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t90 - g(2) * t87, g(1) * t87 - g(2) * t90, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, t53, -t104, -t63, -g(1) * t113 - g(2) * t110 - g(3) * t111, 0, 0, 0, 0, 0, 0, t53, -t63, t104, -g(1) * t109 - g(2) * t105 - g(3) * t108, 0, 0, 0, 0, 0, 0, -g(1) * t59 - g(2) * t57 - g(3) * t67, t106, t63, -g(1) * t102 - g(2) * t103 - g(3) * t101, 0, 0, 0, 0, 0, 0, -g(1) * t52 - g(2) * t50 - g(3) * t61, t107, -t106, -g(1) * t99 - g(2) * t100 - g(3) * t98, 0, 0, 0, 0, 0, 0, -g(1) * (t52 * t95 + t58 * t92) - g(2) * (t50 * t95 + t56 * t92) - g(3) * (t61 * t95 + t66 * t92) -g(1) * (-t52 * t92 + t58 * t95) - g(2) * (-t50 * t92 + t56 * t95) - g(3) * (-t61 * t92 + t66 * t95) -t107, -g(1) * (t52 * pkin(5) + t51 * pkin(9) + t99) - g(2) * (t50 * pkin(5) + t49 * pkin(9) + t100) - g(3) * (t61 * pkin(5) + t60 * pkin(9) + t98);];
U_reg  = t1;
