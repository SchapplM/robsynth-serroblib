% Calculate minimal parameter regressor of potential energy for
% S5PPRRR3
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:58
% EndTime: 2019-12-05 15:16:58
% DurationCPUTime: 0.08s
% Computational Cost: add. (53->34), mult. (100->63), div. (0->0), fcn. (120->10), ass. (0->24)
t152 = g(3) * qJ(1);
t135 = sin(pkin(9));
t136 = sin(pkin(8));
t151 = t135 * t136;
t138 = cos(pkin(8));
t150 = t135 * t138;
t139 = sin(qJ(4));
t149 = t135 * t139;
t141 = cos(qJ(4));
t148 = t135 * t141;
t142 = cos(qJ(3));
t147 = t135 * t142;
t140 = sin(qJ(3));
t146 = t136 * t140;
t145 = t136 * t142;
t144 = t138 * t140;
t143 = t138 * t142;
t137 = cos(pkin(9));
t134 = qJ(4) + qJ(5);
t133 = cos(t134);
t132 = sin(t134);
t131 = t137 * t143 + t146;
t130 = t137 * t145 - t144;
t1 = [-t152, -g(1) * (t138 * pkin(1) + t136 * qJ(2)) - g(2) * (t136 * pkin(1) - t138 * qJ(2)) - t152, 0, -g(1) * t131 - g(2) * t130 - g(3) * t147, -g(1) * (-t137 * t144 + t145) - g(2) * (-t137 * t146 - t143) + g(3) * t135 * t140, 0, 0, 0, 0, 0, -g(1) * (t131 * t141 + t138 * t149) - g(2) * (t130 * t141 + t136 * t149) - g(3) * (-t137 * t139 + t141 * t147), -g(1) * (-t131 * t139 + t138 * t148) - g(2) * (-t130 * t139 + t136 * t148) - g(3) * (-t137 * t141 - t139 * t147), 0, 0, 0, 0, 0, -g(1) * (t131 * t133 + t132 * t150) - g(2) * (t130 * t133 + t132 * t151) - g(3) * (-t137 * t132 + t133 * t147), -g(1) * (-t131 * t132 + t133 * t150) - g(2) * (-t130 * t132 + t133 * t151) - g(3) * (-t132 * t147 - t137 * t133);];
U_reg = t1;
