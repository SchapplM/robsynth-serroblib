% Calculate minimal parameter regressor of potential energy for
% S4RRPR9
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
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:03
% EndTime: 2019-12-31 17:10:03
% DurationCPUTime: 0.06s
% Computational Cost: add. (44->31), mult. (76->48), div. (0->0), fcn. (81->8), ass. (0->19)
t128 = sin(qJ(2));
t139 = g(3) * t128;
t129 = sin(qJ(1));
t130 = cos(qJ(2));
t138 = t129 * t130;
t125 = pkin(7) + qJ(4);
t123 = sin(t125);
t131 = cos(qJ(1));
t137 = t131 * t123;
t124 = cos(t125);
t136 = t131 * t124;
t126 = sin(pkin(7));
t135 = t131 * t126;
t127 = cos(pkin(7));
t134 = t131 * t127;
t133 = g(1) * t131 + g(2) * t129;
t132 = pkin(2) * t130 + qJ(3) * t128 + pkin(1);
t122 = -g(3) * t130 + t133 * t128;
t1 = [0, -t133, g(1) * t129 - g(2) * t131, 0, 0, 0, 0, 0, -t133 * t130 - t139, t122, -g(1) * (t129 * t126 + t130 * t134) - g(2) * (t127 * t138 - t135) - t127 * t139, -g(1) * (t129 * t127 - t130 * t135) - g(2) * (-t126 * t138 - t134) + t126 * t139, -t122, -g(3) * (t128 * pkin(2) - t130 * qJ(3) + pkin(4)) + (g(2) * pkin(5) - g(1) * t132) * t131 + (-g(1) * pkin(5) - g(2) * t132) * t129, 0, 0, 0, 0, 0, -g(1) * (t129 * t123 + t130 * t136) - g(2) * (t124 * t138 - t137) - t124 * t139, -g(1) * (t129 * t124 - t130 * t137) - g(2) * (-t123 * t138 - t136) + t123 * t139;];
U_reg = t1;
