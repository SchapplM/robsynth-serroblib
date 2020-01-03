% Calculate minimal parameter regressor of potential energy for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:36:20
% EndTime: 2020-01-03 11:36:20
% DurationCPUTime: 0.06s
% Computational Cost: add. (76->29), mult. (56->41), div. (0->0), fcn. (54->10), ass. (0->18)
t124 = sin(pkin(9));
t135 = g(1) * t124;
t134 = qJ(2) + pkin(5);
t125 = cos(pkin(9));
t126 = sin(qJ(5));
t133 = t125 * t126;
t128 = cos(qJ(5));
t132 = t125 * t128;
t123 = qJ(1) + pkin(8);
t122 = qJ(3) + t123;
t120 = sin(t122);
t121 = cos(t122);
t131 = g(2) * t120 - g(3) * t121;
t127 = sin(qJ(1));
t129 = cos(qJ(1));
t130 = -g(2) * t127 + g(3) * t129;
t119 = g(2) * t121 + g(3) * t120;
t1 = [0, t130, -g(2) * t129 - g(3) * t127, t130 * pkin(1) - g(1) * t134, 0, -t131, -t119, -t131 * t125 - t135, -g(1) * t125 + t131 * t124, t119, -g(1) * (pkin(6) + t134) - g(2) * (t120 * pkin(3) - t121 * qJ(4) + pkin(2) * sin(t123) + t127 * pkin(1)) - g(3) * (-t121 * pkin(3) - t120 * qJ(4) - pkin(2) * cos(t123) - t129 * pkin(1)), 0, 0, 0, 0, 0, -t128 * t135 - g(2) * (t120 * t132 - t121 * t126) - g(3) * (-t120 * t126 - t121 * t132), t126 * t135 - g(2) * (-t120 * t133 - t121 * t128) - g(3) * (-t120 * t128 + t121 * t133);];
U_reg = t1;
