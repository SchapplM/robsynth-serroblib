% Calculate minimal parameter regressor of potential energy for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRPRPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:06:17
% EndTime: 2021-01-16 01:06:18
% DurationCPUTime: 0.23s
% Computational Cost: add. (249->89), mult. (255->147), div. (0->0), fcn. (283->18), ass. (0->69)
t112 = qJ(2) + pkin(11);
t105 = pkin(6) + t112;
t150 = sin(t105) / 0.2e1;
t106 = pkin(6) - t112;
t149 = sin(t106);
t124 = cos(qJ(2));
t104 = t124 * pkin(2) + pkin(1);
t113 = sin(pkin(10));
t115 = cos(pkin(10));
t114 = sin(pkin(6));
t118 = pkin(7) + qJ(3);
t116 = cos(pkin(6));
t121 = sin(qJ(2));
t136 = t116 * t121;
t91 = pkin(2) * t136 - t114 * t118;
t148 = t113 * t104 + t115 * t91;
t110 = cos(t112);
t147 = t110 * t116;
t108 = sin(t112);
t146 = t113 * t108;
t145 = t113 * t114;
t144 = t114 * t115;
t119 = sin(qJ(6));
t143 = t114 * t119;
t120 = sin(qJ(4));
t142 = t114 * t120;
t141 = t114 * t121;
t122 = cos(qJ(6));
t140 = t114 * t122;
t123 = cos(qJ(4));
t139 = t114 * t123;
t138 = t115 * t108;
t137 = t116 * t120;
t135 = t116 * t124;
t134 = t119 * t113;
t133 = t119 * t115;
t132 = t122 * t113;
t131 = t122 * t115;
t130 = pkin(2) * t141 + t116 * t118 + qJ(1);
t129 = t115 * t142;
t111 = qJ(4) + pkin(12);
t109 = cos(t111);
t128 = t109 * t134;
t127 = t109 * t133;
t126 = t109 * t132;
t125 = t109 * t131;
t117 = -qJ(5) - pkin(8);
t107 = sin(t111);
t103 = t123 * pkin(4) + pkin(3);
t101 = cos(t105);
t98 = cos(t106) / 0.2e1;
t96 = t116 * t107;
t95 = t115 * t110;
t94 = t113 * t110;
t93 = t115 * t104;
t90 = t98 - t101 / 0.2e1;
t89 = t101 / 0.2e1 + t98;
t88 = t149 / 0.2e1 + t150;
t87 = t150 - t149 / 0.2e1;
t85 = t116 * t138 + t94;
t84 = t116 * t146 - t95;
t83 = t114 * t108 * t109 + t96;
t82 = t107 * t140 + t119 * t147;
t81 = t107 * t143 - t122 * t147;
t80 = -t113 * t87 + t95;
t79 = t113 * t89 + t138;
t78 = t115 * t87 + t94;
t77 = -t115 * t89 + t146;
t1 = [-g(3) * qJ(1), 0, -g(1) * (-t113 * t136 + t115 * t124) - g(2) * (t113 * t124 + t115 * t136) - g(3) * t141, -g(1) * (-t113 * t135 - t115 * t121) - g(2) * (-t113 * t121 + t115 * t135) - g(3) * t114 * t124, -g(1) * (-t113 * t91 + t93) - g(2) * t148 - g(3) * t130, 0, 0, 0, 0, 0, -g(1) * (t113 * t142 - t84 * t123) - g(2) * (t85 * t123 - t129) - g(3) * (t108 * t139 + t137), -g(1) * (t113 * t139 + t84 * t120) - g(2) * (-t115 * t139 - t85 * t120) - g(3) * (-t108 * t142 + t116 * t123), -g(1) * (t107 * t145 + t80 * t109) - g(2) * (-t107 * t144 + t78 * t109) - g(3) * (t90 * t109 + t96), -g(1) * (-t80 * t107 + t109 * t145) - g(2) * (-t78 * t107 - t109 * t144) - g(3) * (-t90 * t107 + t116 * t109), -g(1) * t79 - g(2) * t77 + g(3) * t88, -g(1) * (t80 * t103 - t79 * t117 + t93 + (pkin(4) * t142 - t91) * t113) - g(2) * (-pkin(4) * t129 + t78 * t103 - t77 * t117 + t148) - g(3) * (pkin(4) * t137 + t90 * t103 + t88 * t117 + t130), 0, 0, 0, 0, 0, -g(1) * ((-t116 * t126 + t133) * t108 + t110 * t125 + t113 * t82) - g(2) * ((t116 * t125 + t134) * t108 + t110 * t126 - t115 * t82) - g(3) * (-t110 * t143 + t83 * t122), -g(1) * ((t116 * t128 + t131) * t108 - t110 * t127 - t113 * t81) - g(2) * ((-t116 * t127 + t132) * t108 - t110 * t128 + t115 * t81) - g(3) * (-t110 * t140 - t83 * t119);];
U_reg = t1;
