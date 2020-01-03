% Calculate minimal parameter regressor of potential energy for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:35
% EndTime: 2020-01-03 11:57:35
% DurationCPUTime: 0.06s
% Computational Cost: add. (77->30), mult. (58->45), div. (0->0), fcn. (56->10), ass. (0->20)
t138 = g(1) * (qJ(3) + pkin(6) + pkin(5));
t126 = sin(pkin(9));
t137 = g(1) * t126;
t127 = cos(pkin(9));
t128 = sin(qJ(5));
t136 = t127 * t128;
t130 = cos(qJ(5));
t135 = t127 * t130;
t125 = qJ(1) + qJ(2);
t121 = sin(t125);
t129 = sin(qJ(1));
t134 = t129 * pkin(1) + pkin(2) * t121;
t122 = cos(t125);
t131 = cos(qJ(1));
t133 = -t131 * pkin(1) - pkin(2) * t122;
t120 = pkin(8) + t125;
t117 = sin(t120);
t118 = cos(t120);
t132 = g(2) * t117 - g(3) * t118;
t1 = [0, -g(2) * t129 + g(3) * t131, -g(2) * t131 - g(3) * t129, 0, -g(2) * t121 + g(3) * t122, -g(2) * t122 - g(3) * t121, -g(2) * t134 - g(3) * t133 - t138, -t132 * t127 - t137, -g(1) * t127 + t132 * t126, g(2) * t118 + g(3) * t117, -t138 - g(2) * (t117 * pkin(3) - t118 * qJ(4) + t134) - g(3) * (-t118 * pkin(3) - t117 * qJ(4) + t133), 0, 0, 0, 0, 0, -t130 * t137 - g(2) * (t117 * t135 - t118 * t128) - g(3) * (-t117 * t128 - t118 * t135), t128 * t137 - g(2) * (-t117 * t136 - t118 * t130) - g(3) * (-t117 * t130 + t118 * t136);];
U_reg = t1;
