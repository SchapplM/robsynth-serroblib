% Calculate minimal parameter regressor of potential energy for
% S5PPRPR1
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
% U_reg [1x16]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:19
% EndTime: 2019-12-05 15:01:19
% DurationCPUTime: 0.08s
% Computational Cost: add. (79->38), mult. (83->54), div. (0->0), fcn. (84->10), ass. (0->25)
t142 = g(3) * qJ(1);
t125 = pkin(8) + qJ(3);
t121 = sin(t125);
t141 = g(3) * t121;
t124 = pkin(9) + qJ(5);
t120 = sin(t124);
t127 = sin(pkin(7));
t140 = t127 * t120;
t122 = cos(t124);
t139 = t127 * t122;
t126 = sin(pkin(9));
t138 = t127 * t126;
t128 = cos(pkin(9));
t137 = t127 * t128;
t129 = cos(pkin(7));
t136 = t129 * t120;
t135 = t129 * t122;
t134 = t129 * t126;
t133 = t129 * t128;
t132 = g(1) * t129 + g(2) * t127;
t123 = cos(t125);
t131 = pkin(3) * t123 + qJ(4) * t121 + cos(pkin(8)) * pkin(2) + pkin(1);
t130 = -pkin(5) - qJ(2);
t118 = -g(3) * t123 + t132 * t121;
t1 = [-t142, -g(1) * (t129 * pkin(1) + t127 * qJ(2)) - g(2) * (t127 * pkin(1) - t129 * qJ(2)) - t142, 0, -t132 * t123 - t141, t118, -g(1) * (t123 * t133 + t138) - g(2) * (t123 * t137 - t134) - t128 * t141, -g(1) * (-t123 * t134 + t137) - g(2) * (-t123 * t138 - t133) + t126 * t141, -t118, -g(3) * (t121 * pkin(3) - t123 * qJ(4) + sin(pkin(8)) * pkin(2) + qJ(1)) + (-g(1) * t131 - g(2) * t130) * t129 + (g(1) * t130 - g(2) * t131) * t127, 0, 0, 0, 0, 0, -g(1) * (t123 * t135 + t140) - g(2) * (t123 * t139 - t136) - t122 * t141, -g(1) * (-t123 * t136 + t139) - g(2) * (-t123 * t140 - t135) + t120 * t141;];
U_reg = t1;
