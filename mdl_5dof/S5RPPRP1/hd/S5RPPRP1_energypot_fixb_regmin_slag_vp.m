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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:36:24
% EndTime: 2019-12-05 17:36:24
% DurationCPUTime: 0.08s
% Computational Cost: add. (78->36), mult. (81->52), div. (0->0), fcn. (76->8), ass. (0->24)
t136 = cos(qJ(4));
t125 = t136 * pkin(4) + pkin(3);
t130 = sin(pkin(8));
t131 = cos(pkin(8));
t132 = -qJ(5) - pkin(6);
t149 = t125 * t131 - t130 * t132;
t148 = g(1) * t130;
t133 = qJ(2) + pkin(5);
t147 = g(1) * t133;
t129 = qJ(1) + pkin(7);
t127 = cos(t129);
t134 = sin(qJ(4));
t145 = t127 * t134;
t143 = t131 * t134;
t142 = t131 * t136;
t126 = sin(t129);
t137 = cos(qJ(1));
t141 = t137 * pkin(1) + t127 * pkin(2) + t126 * qJ(3);
t135 = sin(qJ(1));
t140 = -t135 * pkin(1) + t127 * qJ(3);
t139 = g(2) * t126 - g(3) * t127;
t138 = g(2) * t135 - g(3) * t137;
t121 = g(1) * t131 + t139 * t130;
t1 = [0, t138, g(2) * t137 + g(3) * t135, t138 * pkin(1) - t147, t139 * t131 - t148, -t121, -g(2) * t127 - g(3) * t126, -t147 - g(2) * (-t126 * pkin(2) + t140) - g(3) * t141, 0, 0, 0, 0, 0, -t136 * t148 - g(2) * (-t126 * t142 + t145) - g(3) * (t126 * t134 + t127 * t142), t134 * t148 - g(2) * (t126 * t143 + t127 * t136) - g(3) * (t126 * t136 - t127 * t143), t121, -g(1) * (t130 * t125 + t131 * t132 + t133) - g(2) * (pkin(4) * t145 + t140) - g(3) * (t149 * t127 + t141) + (-g(2) * (-pkin(2) - t149) - g(3) * pkin(4) * t134) * t126;];
U_reg = t1;
