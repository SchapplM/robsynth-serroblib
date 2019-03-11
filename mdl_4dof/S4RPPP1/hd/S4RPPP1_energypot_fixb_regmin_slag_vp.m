% Calculate minimal parameter regressor of potential energy for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% U_reg [1x15]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPPP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:25
% EndTime: 2019-03-08 18:26:25
% DurationCPUTime: 0.07s
% Computational Cost: add. (75->38), mult. (164->52), div. (0->0), fcn. (177->6), ass. (0->31)
t127 = sin(qJ(1));
t145 = g(1) * t127;
t128 = cos(qJ(1));
t144 = g(2) * t128;
t126 = cos(pkin(4));
t143 = t126 * qJ(2) + pkin(5);
t124 = sin(pkin(4));
t142 = qJ(2) * t124;
t123 = sin(pkin(6));
t141 = t123 * t124;
t125 = cos(pkin(6));
t140 = t124 * t125;
t139 = t127 * t123;
t138 = t127 * t125;
t137 = t128 * t123;
t136 = t128 * t125;
t135 = t128 * pkin(1) + t127 * t142;
t134 = pkin(2) * t141 + t143;
t133 = t128 * t142;
t132 = -t144 + t145;
t112 = -t126 * t136 + t139;
t113 = t126 * t137 + t138;
t121 = t127 * pkin(1);
t131 = t113 * pkin(2) + t112 * qJ(3) + t121;
t114 = t126 * t138 + t137;
t115 = -t126 * t139 + t136;
t130 = t115 * pkin(2) + t114 * qJ(3) + t135;
t107 = -g(1) * t114 - g(2) * t112 + g(3) * t140;
t129 = g(1) * t115 + g(2) * t113 + g(3) * t141;
t109 = -g(3) * t126 - t132 * t124;
t1 = [0, -g(1) * t128 - g(2) * t127, t132, -t129, -t107, t109, -g(1) * t135 - g(2) * (t121 - t133) - g(3) * t143, t109, t129, t107, -g(1) * t130 - g(2) * (t131 - t133) - g(3) * (-qJ(3) * t140 + t134) t109, t107, -t129, -g(1) * (t115 * qJ(4) + t130) - g(2) * (t113 * qJ(4) + t131) - g(3) * (t126 * pkin(3) + t134) + (-pkin(3) * t145 - g(3) * (-qJ(3) * t125 + qJ(4) * t123) - (-pkin(3) - qJ(2)) * t144) * t124;];
U_reg  = t1;
