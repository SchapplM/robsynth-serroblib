% Calculate minimal parameter regressor of potential energy for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:15
% EndTime: 2019-12-31 18:43:15
% DurationCPUTime: 0.06s
% Computational Cost: add. (64->29), mult. (70->41), div. (0->0), fcn. (68->8), ass. (0->22)
t136 = sin(qJ(3));
t149 = g(3) * t136;
t148 = qJ(2) + pkin(5);
t135 = sin(qJ(4));
t139 = cos(qJ(3));
t147 = t135 * t139;
t138 = cos(qJ(4));
t146 = t138 * t139;
t145 = pkin(4) * t135 + pkin(6);
t133 = qJ(1) + pkin(8);
t131 = sin(t133);
t132 = cos(t133);
t144 = g(1) * t132 + g(2) * t131;
t137 = sin(qJ(1));
t140 = cos(qJ(1));
t143 = -g(1) * t140 - g(2) * t137;
t142 = t143 * pkin(1);
t130 = t138 * pkin(4) + pkin(3);
t134 = -qJ(5) - pkin(7);
t141 = t130 * t139 - t134 * t136 + pkin(2);
t129 = -g(3) * t139 + t144 * t136;
t1 = [0, t143, g(1) * t137 - g(2) * t140, -g(3) * t148 + t142, 0, 0, 0, 0, 0, -t144 * t139 - t149, t129, 0, 0, 0, 0, 0, -g(1) * (t131 * t135 + t132 * t146) - g(2) * (t131 * t146 - t132 * t135) - t138 * t149, -g(1) * (t131 * t138 - t132 * t147) - g(2) * (-t131 * t147 - t132 * t138) + t135 * t149, -t129, -g(3) * (t136 * t130 + t139 * t134 + t148) + t142 + (-g(1) * t141 + g(2) * t145) * t132 + (-g(1) * t145 - g(2) * t141) * t131;];
U_reg = t1;
