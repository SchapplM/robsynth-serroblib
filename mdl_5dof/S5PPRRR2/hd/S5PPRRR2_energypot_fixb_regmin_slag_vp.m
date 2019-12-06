% Calculate minimal parameter regressor of potential energy for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:48
% EndTime: 2019-12-05 15:14:48
% DurationCPUTime: 0.06s
% Computational Cost: add. (53->26), mult. (58->41), div. (0->0), fcn. (62->8), ass. (0->22)
t142 = g(3) * qJ(1);
t126 = pkin(9) + qJ(3);
t122 = sin(t126);
t141 = g(3) * t122;
t127 = qJ(4) + qJ(5);
t124 = sin(t127);
t128 = sin(pkin(8));
t140 = t128 * t124;
t125 = cos(t127);
t139 = t128 * t125;
t130 = sin(qJ(4));
t138 = t128 * t130;
t131 = cos(qJ(4));
t137 = t128 * t131;
t129 = cos(pkin(8));
t136 = t129 * t124;
t135 = t129 * t125;
t134 = t129 * t130;
t133 = t129 * t131;
t132 = g(1) * t129 + g(2) * t128;
t123 = cos(t126);
t1 = [-t142, -g(1) * (t129 * pkin(1) + t128 * qJ(2)) - g(2) * (t128 * pkin(1) - t129 * qJ(2)) - t142, 0, -t132 * t123 - t141, -g(3) * t123 + t132 * t122, 0, 0, 0, 0, 0, -g(1) * (t123 * t133 + t138) - g(2) * (t123 * t137 - t134) - t131 * t141, -g(1) * (-t123 * t134 + t137) - g(2) * (-t123 * t138 - t133) + t130 * t141, 0, 0, 0, 0, 0, -g(1) * (t123 * t135 + t140) - g(2) * (t123 * t139 - t136) - t125 * t141, -g(1) * (-t123 * t136 + t139) - g(2) * (-t123 * t140 - t135) + t124 * t141;];
U_reg = t1;
