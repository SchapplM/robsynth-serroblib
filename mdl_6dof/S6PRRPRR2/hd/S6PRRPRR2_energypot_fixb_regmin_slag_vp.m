% Calculate minimal parameter regressor of potential energy for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% U_reg [1x29]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:49
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:46:36
% EndTime: 2021-01-16 03:46:36
% DurationCPUTime: 0.15s
% Computational Cost: add. (161->67), mult. (268->115), div. (0->0), fcn. (334->14), ass. (0->36)
t117 = sin(pkin(11));
t118 = sin(pkin(6));
t138 = t117 * t118;
t119 = cos(pkin(11));
t137 = t118 * t119;
t123 = sin(qJ(3));
t136 = t118 * t123;
t124 = sin(qJ(2));
t135 = t118 * t124;
t126 = cos(qJ(3));
t134 = t118 * t126;
t127 = cos(qJ(2));
t133 = t118 * t127;
t120 = cos(pkin(6));
t132 = t120 * t123;
t131 = t120 * t124;
t130 = t120 * t127;
t105 = t117 * t124 - t119 * t130;
t107 = t117 * t130 + t119 * t124;
t128 = -g(1) * t107 - g(2) * t105 + g(3) * t133;
t125 = cos(qJ(5));
t122 = sin(qJ(5));
t121 = qJ(4) + pkin(8);
t116 = qJ(5) + qJ(6);
t115 = qJ(3) + pkin(12);
t114 = cos(t116);
t113 = sin(t116);
t112 = cos(t115);
t111 = sin(t115);
t110 = t126 * pkin(3) + pkin(2);
t106 = t117 * t127 + t119 * t131;
t104 = t117 * t131 - t119 * t127;
t103 = t120 * t111 + t112 * t135;
t102 = -t104 * t112 + t111 * t138;
t101 = t106 * t112 - t111 * t137;
t1 = [-g(3) * qJ(1), 0, g(1) * t104 - g(2) * t106 - g(3) * t135, -t128, 0, 0, 0, 0, 0, -g(1) * (-t104 * t126 + t117 * t136) - g(2) * (t106 * t126 - t119 * t136) - g(3) * (t124 * t134 + t132), -g(1) * (t104 * t123 + t117 * t134) - g(2) * (-t106 * t123 - t119 * t134) - g(3) * (t120 * t126 - t123 * t135), -g(1) * t102 - g(2) * t101 - g(3) * t103, -g(1) * (t104 * t111 + t112 * t138) - g(2) * (-t106 * t111 - t112 * t137) - g(3) * (-t111 * t135 + t120 * t112), t128, -g(1) * (t119 * pkin(1) - t104 * t110 + t107 * t121) - g(2) * (t117 * pkin(1) + t105 * t121 + t106 * t110) - g(3) * (pkin(3) * t132 + t120 * pkin(7) + qJ(1)) + (-g(3) * (t110 * t124 - t121 * t127) + (-g(1) * t117 + g(2) * t119) * (pkin(3) * t123 + pkin(7))) * t118, 0, 0, 0, 0, 0, -g(1) * (t102 * t125 + t107 * t122) - g(2) * (t101 * t125 + t105 * t122) - g(3) * (t103 * t125 - t122 * t133), -g(1) * (-t102 * t122 + t107 * t125) - g(2) * (-t101 * t122 + t105 * t125) - g(3) * (-t103 * t122 - t125 * t133), 0, 0, 0, 0, 0, -g(1) * (t102 * t114 + t107 * t113) - g(2) * (t101 * t114 + t105 * t113) - g(3) * (t103 * t114 - t113 * t133), -g(1) * (-t102 * t113 + t107 * t114) - g(2) * (-t101 * t113 + t105 * t114) - g(3) * (-t103 * t113 - t114 * t133);];
U_reg = t1;
