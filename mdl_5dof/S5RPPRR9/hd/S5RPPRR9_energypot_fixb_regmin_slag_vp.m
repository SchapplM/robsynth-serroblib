% Calculate minimal parameter regressor of potential energy for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR9_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:47
% EndTime: 2019-12-31 18:02:47
% DurationCPUTime: 0.06s
% Computational Cost: add. (45->26), mult. (90->44), div. (0->0), fcn. (104->8), ass. (0->19)
t136 = sin(qJ(1));
t124 = sin(qJ(4));
t135 = g(3) * t124;
t134 = cos(pkin(8));
t133 = sin(pkin(8));
t123 = sin(qJ(5));
t126 = cos(qJ(4));
t132 = t123 * t126;
t125 = cos(qJ(5));
t131 = t125 * t126;
t127 = cos(qJ(1));
t130 = t127 * pkin(1) + t136 * qJ(2);
t129 = t136 * pkin(1) - t127 * qJ(2);
t113 = -t127 * t134 - t136 * t133;
t114 = t127 * t133 - t136 * t134;
t128 = g(1) * t113 + g(2) * t114;
t116 = -g(1) * t127 - g(2) * t136;
t115 = g(1) * t136 - g(2) * t127;
t1 = [0, t116, t115, t116, -t115, -g(3) * pkin(5) - g(1) * t130 - g(2) * t129, t128, g(1) * t114 - g(2) * t113, -g(1) * (t127 * pkin(2) + t130) - g(2) * (t136 * pkin(2) + t129) - g(3) * (-qJ(3) + pkin(5)), 0, 0, 0, 0, 0, t128 * t126 + t135, g(3) * t126 - t128 * t124, 0, 0, 0, 0, 0, -g(1) * (-t113 * t131 + t114 * t123) - g(2) * (-t113 * t123 - t114 * t131) + t125 * t135, -g(1) * (t113 * t132 + t114 * t125) - g(2) * (-t113 * t125 + t114 * t132) - t123 * t135;];
U_reg = t1;
