% Calculate minimal parameter regressor of potential energy for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% U_reg [1x13]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPPRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:37
% EndTime: 2019-12-05 14:59:37
% DurationCPUTime: 0.08s
% Computational Cost: add. (53->36), mult. (123->64), div. (0->0), fcn. (146->10), ass. (0->28)
t138 = g(3) * qJ(1);
t119 = sin(pkin(9));
t120 = sin(pkin(8));
t137 = t119 * t120;
t126 = sin(qJ(4));
t136 = t120 * t126;
t128 = cos(qJ(4));
t135 = t120 * t128;
t121 = sin(pkin(7));
t123 = cos(pkin(8));
t134 = t121 * t123;
t124 = cos(pkin(7));
t133 = t124 * t119;
t122 = cos(pkin(9));
t132 = t124 * t122;
t131 = t124 * pkin(1) + t121 * qJ(2);
t130 = t121 * pkin(1) - t124 * qJ(2);
t129 = pkin(2) * t123 + qJ(3) * t120;
t127 = cos(qJ(5));
t125 = sin(qJ(5));
t115 = t122 * t135 - t123 * t126;
t114 = t121 * t119 + t123 * t132;
t113 = -t121 * t122 + t123 * t133;
t112 = t122 * t134 - t133;
t111 = t119 * t134 + t132;
t110 = t114 * t128 + t124 * t136;
t109 = t112 * t128 + t121 * t136;
t1 = [-t138, -g(1) * t131 - g(2) * t130 - t138, -g(1) * (t129 * t124 + t131) - g(2) * (t129 * t121 + t130) - g(3) * (t120 * pkin(2) - t123 * qJ(3) + qJ(1)), 0, -g(1) * t110 - g(2) * t109 - g(3) * t115, -g(1) * (-t114 * t126 + t124 * t135) - g(2) * (-t112 * t126 + t121 * t135) - g(3) * (-t122 * t136 - t123 * t128), 0, 0, 0, 0, 0, -g(1) * (t110 * t127 + t113 * t125) - g(2) * (t109 * t127 + t111 * t125) - g(3) * (t115 * t127 + t125 * t137), -g(1) * (-t110 * t125 + t113 * t127) - g(2) * (-t109 * t125 + t111 * t127) - g(3) * (-t115 * t125 + t127 * t137);];
U_reg = t1;
