% Calculate minimal parameter regressor of potential energy for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:12
% EndTime: 2019-12-05 17:10:12
% DurationCPUTime: 0.06s
% Computational Cost: add. (53->24), mult. (61->39), div. (0->0), fcn. (68->10), ass. (0->23)
t129 = qJ(2) + qJ(3);
t125 = sin(t129);
t145 = g(3) * t125;
t128 = qJ(4) + qJ(5);
t124 = sin(t128);
t130 = sin(pkin(9));
t144 = t130 * t124;
t126 = cos(t128);
t143 = t130 * t126;
t132 = sin(qJ(4));
t142 = t130 * t132;
t134 = cos(qJ(4));
t141 = t130 * t134;
t131 = cos(pkin(9));
t140 = t131 * t124;
t139 = t131 * t126;
t138 = t131 * t132;
t137 = t131 * t134;
t136 = g(1) * t131 + g(2) * t130;
t135 = cos(qJ(2));
t133 = sin(qJ(2));
t127 = cos(t129);
t1 = [-g(3) * qJ(1), 0, -g(3) * t133 - t136 * t135, -g(3) * t135 + t136 * t133, 0, -t136 * t127 - t145, -g(3) * t127 + t136 * t125, 0, 0, 0, 0, 0, -g(1) * (t127 * t137 + t142) - g(2) * (t127 * t141 - t138) - t134 * t145, -g(1) * (-t127 * t138 + t141) - g(2) * (-t127 * t142 - t137) + t132 * t145, 0, 0, 0, 0, 0, -g(1) * (t127 * t139 + t144) - g(2) * (t127 * t143 - t140) - t126 * t145, -g(1) * (-t127 * t140 + t143) - g(2) * (-t127 * t144 - t139) + t124 * t145;];
U_reg = t1;
