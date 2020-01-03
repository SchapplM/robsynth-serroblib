% Calculate minimal parameter regressor of potential energy for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:47
% EndTime: 2019-12-31 16:36:47
% DurationCPUTime: 0.06s
% Computational Cost: add. (47->29), mult. (119->57), div. (0->0), fcn. (152->10), ass. (0->24)
t128 = sin(pkin(4));
t132 = sin(qJ(3));
t142 = t128 * t132;
t133 = sin(qJ(2));
t141 = t128 * t133;
t135 = cos(qJ(3));
t140 = t128 * t135;
t136 = cos(qJ(2));
t139 = t128 * t136;
t130 = cos(pkin(4));
t138 = t130 * t133;
t137 = t130 * t136;
t134 = cos(qJ(4));
t131 = sin(qJ(4));
t129 = cos(pkin(8));
t127 = sin(pkin(8));
t126 = t130 * t132 + t133 * t140;
t125 = -t127 * t138 + t129 * t136;
t124 = t127 * t137 + t129 * t133;
t123 = t127 * t136 + t129 * t138;
t122 = t127 * t133 - t129 * t137;
t121 = t125 * t135 + t127 * t142;
t120 = t123 * t135 - t129 * t142;
t1 = [-g(3) * qJ(1), 0, -g(1) * t125 - g(2) * t123 - g(3) * t141, g(1) * t124 + g(2) * t122 - g(3) * t139, 0, 0, 0, 0, 0, -g(1) * t121 - g(2) * t120 - g(3) * t126, -g(1) * (-t125 * t132 + t127 * t140) - g(2) * (-t123 * t132 - t129 * t140) - g(3) * (t130 * t135 - t132 * t141), 0, 0, 0, 0, 0, -g(1) * (t121 * t134 + t124 * t131) - g(2) * (t120 * t134 + t122 * t131) - g(3) * (t126 * t134 - t131 * t139), -g(1) * (-t121 * t131 + t124 * t134) - g(2) * (-t120 * t131 + t122 * t134) - g(3) * (-t126 * t131 - t134 * t139);];
U_reg = t1;
