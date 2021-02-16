% Calculate minimal parameter regressor of potential energy for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:57
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:54:47
% EndTime: 2021-01-16 02:54:47
% DurationCPUTime: 0.26s
% Computational Cost: add. (252->73), mult. (425->110), div. (0->0), fcn. (527->12), ass. (0->49)
t110 = sin(pkin(10));
t111 = sin(pkin(6));
t136 = t110 * t111;
t112 = cos(pkin(10));
t135 = t111 * t112;
t116 = sin(qJ(3));
t134 = t111 * t116;
t117 = sin(qJ(2));
t133 = t111 * t117;
t119 = cos(qJ(3));
t132 = t111 * t119;
t120 = cos(qJ(2));
t131 = t111 * t120;
t113 = cos(pkin(6));
t130 = t113 * t116;
t129 = t113 * t117;
t128 = t113 * t120;
t127 = t110 * t134;
t103 = t119 * pkin(3) + pkin(2);
t114 = qJ(4) + pkin(8);
t92 = t110 * t129 - t112 * t120;
t95 = t110 * t128 + t112 * t117;
t126 = t112 * pkin(1) + pkin(3) * t127 + pkin(7) * t136 - t92 * t103 + t95 * t114;
t115 = sin(qJ(5));
t118 = cos(qJ(5));
t109 = qJ(3) + pkin(11);
t104 = sin(t109);
t105 = cos(t109);
t94 = t110 * t120 + t112 * t129;
t81 = -t104 * t135 + t94 * t105;
t93 = t110 * t117 - t112 * t128;
t75 = t81 * t115 - t93 * t118;
t82 = t104 * t136 - t92 * t105;
t77 = t82 * t115 - t95 * t118;
t88 = t113 * t104 + t105 * t133;
t83 = t88 * t115 + t118 * t131;
t125 = g(1) * t77 + g(2) * t75 + g(3) * t83;
t79 = t92 * t104 + t105 * t136;
t80 = t94 * t104 + t105 * t135;
t87 = t104 * t133 - t113 * t105;
t124 = g(1) * t79 - g(2) * t80 - g(3) * t87;
t123 = pkin(3) * t130 + t113 * pkin(7) + t103 * t133 - t114 * t131 + qJ(1);
t122 = -g(1) * t95 - g(2) * t93 + g(3) * t131;
t121 = t93 * t114 + t110 * pkin(1) + t94 * t103 + (-pkin(3) * t116 - pkin(7)) * t135;
t84 = -t115 * t131 + t88 * t118;
t78 = t95 * t115 + t82 * t118;
t76 = t93 * t115 + t81 * t118;
t74 = -g(1) * t78 - g(2) * t76 - g(3) * t84;
t1 = [-g(3) * qJ(1), 0, g(1) * t92 - g(2) * t94 - g(3) * t133, -t122, 0, 0, 0, 0, 0, -g(1) * (-t92 * t119 + t127) - g(2) * (-t112 * t134 + t94 * t119) - g(3) * (t117 * t132 + t130), -g(1) * (t110 * t132 + t92 * t116) - g(2) * (-t112 * t132 - t94 * t116) - g(3) * (t113 * t119 - t116 * t133), -g(1) * t82 - g(2) * t81 - g(3) * t88, -t124, t122, -g(1) * t126 - g(2) * t121 - g(3) * t123, 0, 0, 0, 0, 0, t74, t125, t74, t124, -t125, -g(1) * (t82 * pkin(4) + t78 * pkin(5) - t79 * pkin(9) + t77 * qJ(6) + t126) - g(2) * (t81 * pkin(4) + t76 * pkin(5) + t80 * pkin(9) + t75 * qJ(6) + t121) - g(3) * (t88 * pkin(4) + t84 * pkin(5) + t87 * pkin(9) + t83 * qJ(6) + t123);];
U_reg = t1;
