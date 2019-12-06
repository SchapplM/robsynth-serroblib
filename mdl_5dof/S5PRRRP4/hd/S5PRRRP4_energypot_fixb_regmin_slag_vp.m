% Calculate minimal parameter regressor of potential energy for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:46:19
% EndTime: 2019-12-05 16:46:20
% DurationCPUTime: 0.07s
% Computational Cost: add. (85->34), mult. (106->45), div. (0->0), fcn. (114->8), ass. (0->25)
t123 = qJ(2) + qJ(3);
t122 = cos(t123);
t129 = cos(qJ(2));
t139 = t129 * pkin(2) + pkin(3) * t122 + pkin(1);
t121 = sin(t123);
t137 = g(3) * t121;
t124 = sin(pkin(8));
t126 = sin(qJ(4));
t136 = t124 * t126;
t128 = cos(qJ(4));
t135 = t124 * t128;
t125 = cos(pkin(8));
t134 = t125 * t126;
t133 = t125 * t128;
t132 = g(1) * t125 + g(2) * t124;
t115 = t122 * t136 + t133;
t117 = t122 * t134 - t135;
t131 = g(1) * t117 + g(2) * t115 + t126 * t137;
t130 = -pkin(6) - pkin(5);
t127 = sin(qJ(2));
t118 = t122 * t133 + t136;
t116 = t122 * t135 - t134;
t114 = -g(3) * t122 + t132 * t121;
t113 = -g(1) * t118 - g(2) * t116 - t128 * t137;
t1 = [-g(3) * qJ(1), 0, -g(3) * t127 - t132 * t129, -g(3) * t129 + t132 * t127, 0, -t132 * t122 - t137, t114, 0, 0, 0, 0, 0, t113, t131, t113, -t114, -t131, -g(1) * (t118 * pkin(4) + t117 * qJ(5) - t124 * t130 + t139 * t125) - g(2) * (t116 * pkin(4) + t115 * qJ(5) + t139 * t124 + t125 * t130) - g(3) * (t127 * pkin(2) - t122 * pkin(7) + qJ(1)) + (-g(3) * (pkin(4) * t128 + qJ(5) * t126 + pkin(3)) - t132 * pkin(7)) * t121;];
U_reg = t1;
