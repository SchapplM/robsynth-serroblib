% Calculate minimal parameter regressor of potential energy for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% U_reg [1x26]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PRRPPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:06:03
% EndTime: 2021-01-16 02:06:03
% DurationCPUTime: 0.16s
% Computational Cost: add. (226->77), mult. (363->122), div. (0->0), fcn. (445->14), ass. (0->44)
t128 = sin(pkin(10));
t129 = sin(pkin(6));
t152 = t128 * t129;
t131 = cos(pkin(10));
t151 = t129 * t131;
t134 = sin(qJ(3));
t150 = t129 * t134;
t135 = sin(qJ(2));
t149 = t129 * t135;
t136 = cos(qJ(3));
t148 = t129 * t136;
t137 = cos(qJ(2));
t147 = t129 * t137;
t132 = cos(pkin(6));
t146 = t132 * t134;
t145 = t132 * t135;
t144 = t132 * t137;
t143 = t128 * t150;
t106 = t128 * t145 - t131 * t137;
t109 = t128 * t144 + t131 * t135;
t117 = t136 * pkin(3) + pkin(2);
t133 = qJ(4) + pkin(8);
t142 = t131 * pkin(1) + pkin(3) * t143 + pkin(7) * t152 - t106 * t117 + t109 * t133;
t126 = qJ(3) + pkin(11);
t119 = sin(t126);
t121 = cos(t126);
t103 = t119 * t149 - t132 * t121;
t97 = t106 * t119 + t121 * t152;
t108 = t128 * t137 + t131 * t145;
t98 = t108 * t119 + t121 * t151;
t141 = g(1) * t97 - g(2) * t98 - g(3) * t103;
t140 = pkin(3) * t146 + t132 * pkin(7) + t117 * t149 - t133 * t147 + qJ(1);
t107 = t128 * t135 - t131 * t144;
t139 = -g(1) * t109 - g(2) * t107 + g(3) * t147;
t138 = t107 * t133 + t108 * t117 + t128 * pkin(1) + (-pkin(3) * t134 - pkin(7)) * t151;
t130 = cos(pkin(12));
t127 = sin(pkin(12));
t125 = pkin(12) + qJ(6);
t120 = cos(t125);
t118 = sin(t125);
t104 = t132 * t119 + t121 * t149;
t100 = -t106 * t121 + t119 * t152;
t99 = t108 * t121 - t119 * t151;
t1 = [-g(3) * qJ(1), 0, g(1) * t106 - g(2) * t108 - g(3) * t149, -t139, 0, 0, 0, 0, 0, -g(1) * (-t106 * t136 + t143) - g(2) * (t108 * t136 - t131 * t150) - g(3) * (t135 * t148 + t146), -g(1) * (t106 * t134 + t128 * t148) - g(2) * (-t108 * t134 - t131 * t148) - g(3) * (t132 * t136 - t134 * t149), -g(1) * t100 - g(2) * t99 - g(3) * t104, -t141, t139, -g(1) * t142 - g(2) * t138 - g(3) * t140, -g(1) * (t100 * t130 + t109 * t127) - g(2) * (t107 * t127 + t99 * t130) - g(3) * (t104 * t130 - t127 * t147), -g(1) * (-t100 * t127 + t109 * t130) - g(2) * (t107 * t130 - t99 * t127) - g(3) * (-t104 * t127 - t130 * t147), t141, -g(1) * (t100 * pkin(4) - t97 * qJ(5) + t142) - g(2) * (t99 * pkin(4) + t98 * qJ(5) + t138) - g(3) * (t104 * pkin(4) + t103 * qJ(5) + t140), 0, 0, 0, 0, 0, -g(1) * (t100 * t120 + t109 * t118) - g(2) * (t107 * t118 + t99 * t120) - g(3) * (t104 * t120 - t118 * t147), -g(1) * (-t100 * t118 + t109 * t120) - g(2) * (t107 * t120 - t99 * t118) - g(3) * (-t104 * t118 - t120 * t147);];
U_reg = t1;
