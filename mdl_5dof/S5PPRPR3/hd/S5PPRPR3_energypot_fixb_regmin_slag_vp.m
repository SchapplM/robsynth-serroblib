% Calculate minimal parameter regressor of potential energy for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% U_reg [1x13]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:21
% EndTime: 2019-12-05 15:05:21
% DurationCPUTime: 0.08s
% Computational Cost: add. (55->39), mult. (88->63), div. (0->0), fcn. (95->10), ass. (0->28)
t137 = g(3) * qJ(1);
t118 = sin(pkin(8));
t136 = g(3) * t118;
t122 = -qJ(4) - pkin(5);
t135 = t118 * t122;
t123 = sin(qJ(5));
t134 = t118 * t123;
t125 = cos(qJ(5));
t133 = t118 * t125;
t119 = sin(pkin(7));
t120 = cos(pkin(8));
t132 = t119 * t120;
t124 = sin(qJ(3));
t131 = t119 * t124;
t126 = cos(qJ(3));
t130 = t119 * t126;
t121 = cos(pkin(7));
t129 = t121 * t124;
t128 = t121 * t126;
t127 = t121 * pkin(1) + t119 * qJ(2);
t117 = qJ(3) + pkin(9);
t115 = t119 * pkin(1);
t113 = cos(t117);
t112 = sin(t117);
t111 = t126 * pkin(3) + pkin(2);
t110 = t121 * t120 * t113 + t119 * t112;
t109 = -t121 * t112 + t113 * t132;
t1 = [-t137, -g(1) * t127 - g(2) * (-t121 * qJ(2) + t115) - t137, 0, -g(1) * (t120 * t128 + t131) - g(2) * (t120 * t130 - t129) - t126 * t136, -g(1) * (-t120 * t129 + t130) - g(2) * (-t120 * t131 - t128) + t124 * t136, -g(1) * (pkin(3) * t131 + t127) - g(2) * (t111 * t132 - t119 * t135 + t115) - g(3) * (t118 * t111 + t120 * t122 + qJ(1)) + (-g(1) * (t111 * t120 - t135) - g(2) * (-pkin(3) * t124 - qJ(2))) * t121, 0, 0, 0, 0, 0, -g(1) * (t110 * t125 + t121 * t134) - g(2) * (t109 * t125 + t119 * t134) - g(3) * (t113 * t133 - t120 * t123), -g(1) * (-t110 * t123 + t121 * t133) - g(2) * (-t109 * t123 + t119 * t133) - g(3) * (-t113 * t134 - t120 * t125);];
U_reg = t1;
