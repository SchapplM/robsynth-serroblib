% Calculate minimal parameter regressor of potential energy for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:43
% EndTime: 2019-12-31 17:06:44
% DurationCPUTime: 0.07s
% Computational Cost: add. (30->22), mult. (46->32), div. (0->0), fcn. (47->8), ass. (0->18)
t125 = qJ(2) + pkin(7);
t138 = g(3) * sin(t125);
t127 = sin(qJ(4));
t129 = sin(qJ(1));
t137 = t129 * t127;
t130 = cos(qJ(4));
t136 = t129 * t130;
t132 = cos(qJ(1));
t135 = t132 * t127;
t134 = t132 * t130;
t133 = g(1) * t132 + g(2) * t129;
t131 = cos(qJ(2));
t128 = sin(qJ(2));
t126 = -pkin(5) - qJ(3);
t124 = cos(t125);
t122 = t131 * pkin(2) + pkin(1);
t121 = g(1) * t129 - g(2) * t132;
t1 = [0, -t133, t121, 0, 0, 0, 0, 0, -g(3) * t128 - t133 * t131, -g(3) * t131 + t133 * t128, -t121, -g(1) * (t132 * t122 - t129 * t126) - g(2) * (t129 * t122 + t132 * t126) - g(3) * (t128 * pkin(2) + pkin(4)), 0, 0, 0, 0, 0, -g(1) * (t124 * t134 + t137) - g(2) * (t124 * t136 - t135) - t130 * t138, -g(1) * (-t124 * t135 + t136) - g(2) * (-t124 * t137 - t134) + t127 * t138;];
U_reg = t1;
