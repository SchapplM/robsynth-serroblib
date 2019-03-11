% Calculate inertial parameters regressor of potential energy for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:43:34
% EndTime: 2019-03-08 20:43:35
% DurationCPUTime: 0.25s
% Computational Cost: add. (260->97), mult. (470->137), div. (0->0), fcn. (547->12), ass. (0->55)
t89 = sin(pkin(11));
t126 = g(1) * t89;
t91 = cos(pkin(11));
t125 = g(2) * t91;
t92 = cos(pkin(6));
t98 = cos(qJ(2));
t115 = t92 * t98;
t95 = sin(qJ(2));
t72 = -t115 * t91 + t89 * t95;
t94 = sin(qJ(4));
t124 = t72 * t94;
t74 = t115 * t89 + t91 * t95;
t123 = t74 * t94;
t90 = sin(pkin(6));
t122 = t89 * t90;
t121 = t90 * t91;
t120 = t90 * t94;
t119 = t90 * t95;
t97 = cos(qJ(4));
t118 = t90 * t97;
t117 = t90 * t98;
t116 = t92 * t95;
t114 = pkin(1) * t91 + pkin(7) * t122;
t113 = qJ(3) * t98;
t112 = pkin(7) * t92 + qJ(1);
t111 = pkin(7) * t121;
t110 = pkin(2) * t119 + t112;
t82 = pkin(4) * t97 + pkin(3);
t109 = t82 * t92 + t110;
t73 = t116 * t91 + t89 * t98;
t85 = t89 * pkin(1);
t108 = pkin(2) * t73 + t72 * qJ(3) + t85;
t107 = -t125 + t126;
t75 = -t116 * t89 + t91 * t98;
t106 = pkin(2) * t75 + qJ(3) * t74 + t114;
t99 = -pkin(9) - pkin(8);
t105 = pkin(4) * t124 - t73 * t99 + t108;
t88 = qJ(4) + qJ(5);
t83 = sin(t88);
t84 = cos(t88);
t58 = t122 * t83 - t74 * t84;
t60 = t121 * t83 + t72 * t84;
t66 = t117 * t84 + t83 * t92;
t104 = g(1) * t58 - g(2) * t60 + g(3) * t66;
t103 = -g(1) * t74 - g(2) * t72 + g(3) * t117;
t102 = g(1) * t75 + g(2) * t73 + g(3) * t119;
t101 = pkin(4) * t123 + t122 * t82 - t75 * t99 + t106;
t100 = (-g(3) * (-t95 * t99 + (-pkin(4) * t94 - qJ(3)) * t98) - (-pkin(7) - t82) * t125) * t90;
t96 = cos(qJ(6));
t93 = sin(qJ(6));
t67 = -t117 * t83 + t84 * t92;
t65 = -g(3) * t92 - t107 * t90;
t61 = -t121 * t84 + t72 * t83;
t59 = t122 * t84 + t74 * t83;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t91 - g(2) * t89, t107, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t102, -t103, t65, -g(1) * t114 - g(2) * (t85 - t111) - g(3) * t112, 0, 0, 0, 0, 0, 0, t65, t102, t103, -g(1) * t106 - g(2) * (t108 - t111) - g(3) * (-t113 * t90 + t110) 0, 0, 0, 0, 0, 0, -g(1) * (t118 * t89 + t123) - g(2) * (-t118 * t91 + t124) - g(3) * (-t117 * t94 + t92 * t97) -g(1) * (-t120 * t89 + t74 * t97) - g(2) * (t120 * t91 + t72 * t97) - g(3) * (-t117 * t97 - t92 * t94) -t102, -g(1) * (pkin(8) * t75 + t106) - g(2) * (t73 * pkin(8) + t108) - g(3) * (t92 * pkin(3) + t110) + (-pkin(3) * t126 - g(3) * (pkin(8) * t95 - t113) - (-pkin(3) - pkin(7)) * t125) * t90, 0, 0, 0, 0, 0, 0, -g(1) * t59 - g(2) * t61 - g(3) * t67, t104, -t102, -g(1) * t101 - g(2) * t105 - g(3) * t109 + t100, 0, 0, 0, 0, 0, 0, -g(1) * (t59 * t96 + t75 * t93) - g(2) * (t61 * t96 + t73 * t93) - g(3) * (t119 * t93 + t67 * t96) -g(1) * (-t59 * t93 + t75 * t96) - g(2) * (-t61 * t93 + t73 * t96) - g(3) * (t119 * t96 - t67 * t93) -t104, -g(1) * (pkin(5) * t59 + pkin(10) * t58 + t101) - g(2) * (t61 * pkin(5) - t60 * pkin(10) + t105) - g(3) * (t67 * pkin(5) + t66 * pkin(10) + t109) + t100;];
U_reg  = t1;
