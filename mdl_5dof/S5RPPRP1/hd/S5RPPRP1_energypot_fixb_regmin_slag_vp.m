% Calculate minimal parameter regressor of potential energy for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:48
% EndTime: 2020-01-03 11:25:48
% DurationCPUTime: 0.07s
% Computational Cost: add. (78->37), mult. (81->53), div. (0->0), fcn. (76->8), ass. (0->24)
t126 = sin(pkin(8));
t143 = g(1) * t126;
t129 = qJ(2) + pkin(5);
t142 = g(1) * t129;
t133 = cos(qJ(1));
t141 = t133 * pkin(1);
t127 = cos(pkin(8));
t130 = sin(qJ(4));
t140 = t127 * t130;
t132 = cos(qJ(4));
t139 = t127 * t132;
t125 = qJ(1) + pkin(7);
t122 = sin(t125);
t131 = sin(qJ(1));
t138 = t131 * pkin(1) + t122 * pkin(2);
t137 = pkin(4) * t130 + qJ(3);
t123 = cos(t125);
t136 = -g(2) * t122 + g(3) * t123;
t135 = -g(2) * t131 + g(3) * t133;
t121 = t132 * pkin(4) + pkin(3);
t128 = -qJ(5) - pkin(6);
t134 = t121 * t127 - t126 * t128;
t119 = g(1) * t127 + t136 * t126;
t1 = [0, t135, -g(2) * t133 - g(3) * t131, t135 * pkin(1) - t142, t136 * t127 - t143, -t119, g(2) * t123 + g(3) * t122, -t142 - g(2) * (-t123 * qJ(3) + t138) - g(3) * (-t123 * pkin(2) - t122 * qJ(3) - t141), 0, 0, 0, 0, 0, -t132 * t143 - g(2) * (t122 * t139 - t123 * t130) - g(3) * (-t122 * t130 - t123 * t139), t130 * t143 - g(2) * (-t122 * t140 - t123 * t132) - g(3) * (-t122 * t132 + t123 * t140), t119, -g(1) * (t126 * t121 + t127 * t128 + t129) - g(2) * t138 + g(3) * t141 + (-g(2) * t134 + g(3) * t137) * t122 + (g(2) * t137 - g(3) * (-pkin(2) - t134)) * t123;];
U_reg = t1;
