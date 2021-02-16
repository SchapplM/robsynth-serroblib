% Calculate minimal parameter regressor of potential energy for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:04
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:02:42
% EndTime: 2021-01-15 22:02:42
% DurationCPUTime: 0.17s
% Computational Cost: add. (127->69), mult. (175->120), div. (0->0), fcn. (208->16), ass. (0->52)
t113 = qJ(2) + pkin(10);
t112 = cos(t113);
t114 = sin(pkin(5));
t146 = t112 * t114;
t118 = sin(qJ(4));
t145 = t114 * t118;
t119 = sin(qJ(2));
t144 = t114 * t119;
t117 = sin(qJ(5));
t120 = sin(qJ(1));
t143 = t117 * t120;
t124 = cos(qJ(1));
t142 = t117 * t124;
t111 = sin(t113);
t141 = t120 * t111;
t140 = t120 * t119;
t122 = cos(qJ(4));
t139 = t120 * t122;
t123 = cos(qJ(2));
t138 = t120 * t123;
t121 = cos(qJ(5));
t137 = t121 * t120;
t136 = t121 * t124;
t135 = t122 * t124;
t134 = t124 * t111;
t133 = t124 * t119;
t132 = t124 * t123;
t131 = t120 * t145;
t130 = t122 * t137;
t129 = t117 * t139;
t128 = t124 * t145;
t127 = t117 * t135;
t126 = t121 * t135;
t125 = g(1) * t120 - g(2) * t124;
t116 = pkin(7) + qJ(3);
t115 = cos(pkin(5));
t110 = pkin(5) - t113;
t109 = pkin(5) + t113;
t108 = t123 * pkin(2) + pkin(1);
t107 = cos(t110);
t106 = cos(t109);
t105 = sin(t110);
t104 = sin(t109);
t103 = t124 * t112;
t102 = t120 * t112;
t101 = t106 + t107;
t100 = -t104 + t105;
t99 = t115 * t119 * pkin(2) - t114 * t116;
t98 = t114 * t111 * t122 + t115 * t118;
t97 = t115 * t134 + t102;
t96 = t115 * t141 - t103;
t1 = [0, -g(1) * t124 - g(2) * t120, t125, 0, 0, 0, 0, 0, -g(1) * (-t115 * t140 + t132) - g(2) * (t115 * t133 + t138) - g(3) * t144, -g(1) * (-t115 * t138 - t133) - g(2) * (t115 * t132 - t140) - g(3) * t114 * t123, -g(1) * (t103 + t120 * t100 / 0.2e1) - g(2) * (t102 - t124 * t100 / 0.2e1) - g(3) * (t107 / 0.2e1 - t106 / 0.2e1), -g(1) * (-t134 - t120 * t101 / 0.2e1) - g(2) * (-t141 + t124 * t101 / 0.2e1) - g(3) * (t105 / 0.2e1 + t104 / 0.2e1), -g(3) * t115 - t125 * t114, -g(1) * (t124 * t108 - t99 * t120) - g(2) * (t120 * t108 + t124 * t99) - g(3) * (pkin(2) * t144 + t115 * t116 + pkin(6)), 0, 0, 0, 0, 0, -g(1) * (-t96 * t122 + t131) - g(2) * (t97 * t122 - t128) - g(3) * t98, -g(1) * (t114 * t139 + t96 * t118) - g(2) * (-t114 * t135 - t97 * t118) - g(3) * (-t111 * t145 + t115 * t122), 0, 0, 0, 0, 0, -g(1) * ((t115 * t143 + t126) * t112 + (-t115 * t130 + t142) * t111 + t121 * t131) - g(2) * ((-t115 * t142 + t130) * t112 + (t115 * t126 + t143) * t111 - t121 * t128) - g(3) * (-t117 * t146 + t98 * t121), -g(1) * ((t115 * t137 - t127) * t112 + (t115 * t129 + t136) * t111 - t117 * t131) - g(2) * ((-t115 * t136 - t129) * t112 + (-t115 * t127 + t137) * t111 + t117 * t128) - g(3) * (-t98 * t117 - t121 * t146);];
U_reg = t1;
