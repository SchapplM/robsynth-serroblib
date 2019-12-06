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
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:51:32
% EndTime: 2019-12-05 17:51:32
% DurationCPUTime: 0.06s
% Computational Cost: add. (76->28), mult. (56->41), div. (0->0), fcn. (54->10), ass. (0->18)
t125 = sin(pkin(9));
t136 = g(1) * t125;
t135 = qJ(2) + pkin(5);
t126 = cos(pkin(9));
t127 = sin(qJ(5));
t134 = t126 * t127;
t129 = cos(qJ(5));
t133 = t126 * t129;
t124 = qJ(1) + pkin(8);
t123 = qJ(3) + t124;
t121 = sin(t123);
t122 = cos(t123);
t132 = g(2) * t121 - g(3) * t122;
t128 = sin(qJ(1));
t130 = cos(qJ(1));
t131 = g(2) * t128 - g(3) * t130;
t120 = g(2) * t122 + g(3) * t121;
t1 = [0, t131, g(2) * t130 + g(3) * t128, t131 * pkin(1) - g(1) * t135, 0, t132, t120, t132 * t126 - t136, -g(1) * t126 - t132 * t125, -t120, -g(1) * (pkin(6) + t135) - g(2) * (-t121 * pkin(3) + t122 * qJ(4) - pkin(2) * sin(t124) - t128 * pkin(1)) - g(3) * (t122 * pkin(3) + t121 * qJ(4) + pkin(2) * cos(t124) + t130 * pkin(1)), 0, 0, 0, 0, 0, -t129 * t136 - g(2) * (-t121 * t133 + t122 * t127) - g(3) * (t121 * t127 + t122 * t133), t127 * t136 - g(2) * (t121 * t134 + t122 * t129) - g(3) * (t121 * t129 - t122 * t134);];
U_reg = t1;
