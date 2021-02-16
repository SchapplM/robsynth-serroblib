% Calculate minimal parameter regressor of potential energy for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:24
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_energypot_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:21:16
% EndTime: 2021-01-16 02:21:16
% DurationCPUTime: 0.24s
% Computational Cost: add. (195->59), mult. (333->88), div. (0->0), fcn. (402->12), ass. (0->44)
t96 = sin(pkin(10));
t97 = sin(pkin(6));
t124 = t96 * t97;
t106 = cos(qJ(2));
t98 = cos(pkin(10));
t115 = t98 * t106;
t103 = sin(qJ(2));
t118 = t96 * t103;
t99 = cos(pkin(6));
t78 = t99 * t118 - t115;
t95 = qJ(3) + pkin(11);
t90 = sin(t95);
t91 = cos(t95);
t69 = t91 * t124 + t78 * t90;
t123 = t97 * t98;
t116 = t98 * t103;
t117 = t96 * t106;
t80 = t99 * t116 + t117;
t70 = t91 * t123 + t80 * t90;
t121 = t103 * t97;
t75 = t90 * t121 - t99 * t91;
t128 = g(1) * t69 - g(2) * t70 - g(3) * t75;
t102 = sin(qJ(3));
t122 = t102 * t97;
t105 = cos(qJ(3));
t120 = t105 * t97;
t119 = t106 * t97;
t114 = t99 * t102;
t113 = t96 * t122;
t100 = qJ(4) + pkin(8);
t81 = t99 * t117 + t116;
t89 = t105 * pkin(3) + pkin(2);
t112 = t98 * pkin(1) + pkin(3) * t113 + pkin(7) * t124 + t81 * t100 - t78 * t89;
t71 = -t90 * t123 + t80 * t91;
t72 = t90 * t124 - t78 * t91;
t76 = t91 * t121 + t99 * t90;
t110 = g(1) * t72 + g(2) * t71 + g(3) * t76;
t109 = pkin(3) * t114 + t99 * pkin(7) - t100 * t119 + t89 * t121 + qJ(1);
t79 = -t99 * t115 + t118;
t68 = -g(1) * t81 - g(2) * t79 + g(3) * t119;
t108 = t79 * t100 + t80 * t89 + t96 * pkin(1) + (-pkin(3) * t102 - pkin(7)) * t123;
t104 = cos(qJ(6));
t101 = sin(qJ(6));
t1 = [-g(3) * qJ(1), 0, g(1) * t78 - g(2) * t80 - g(3) * t121, -t68, 0, 0, 0, 0, 0, -g(1) * (-t78 * t105 + t113) - g(2) * (t80 * t105 - t98 * t122) - g(3) * (t103 * t120 + t114), -g(1) * (t78 * t102 + t96 * t120) - g(2) * (-t80 * t102 - t98 * t120) - g(3) * (-t102 * t121 + t99 * t105), -t110, -t128, t68, -g(1) * t112 - g(2) * t108 - g(3) * t109, t68, t110, t128, -g(1) * (t72 * pkin(4) - t69 * qJ(5) + t112) - g(2) * (t71 * pkin(4) + t70 * qJ(5) + t108) - g(3) * (t76 * pkin(4) + t75 * qJ(5) + t109), 0, 0, 0, 0, 0, t101 * t128 + t68 * t104, -t68 * t101 + t104 * t128;];
U_reg = t1;
