% Calculate minimal parameter regressor of potential energy for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:47:31
% EndTime: 2019-12-05 15:47:31
% DurationCPUTime: 0.06s
% Computational Cost: add. (52->29), mult. (61->45), div. (0->0), fcn. (65->10), ass. (0->24)
t132 = qJ(2) + pkin(9);
t150 = g(3) * sin(t132);
t133 = qJ(4) + qJ(5);
t130 = sin(t133);
t134 = sin(pkin(8));
t149 = t134 * t130;
t131 = cos(t133);
t148 = t134 * t131;
t137 = sin(qJ(4));
t147 = t134 * t137;
t139 = cos(qJ(4));
t146 = t134 * t139;
t135 = cos(pkin(8));
t145 = t135 * t130;
t144 = t135 * t131;
t143 = t135 * t137;
t142 = t135 * t139;
t141 = g(1) * t135 + g(2) * t134;
t140 = cos(qJ(2));
t138 = sin(qJ(2));
t136 = -qJ(3) - pkin(5);
t129 = cos(t132);
t127 = t140 * pkin(2) + pkin(1);
t1 = [-g(3) * qJ(1), 0, -g(3) * t138 - t141 * t140, -g(3) * t140 + t141 * t138, -g(1) * (t135 * t127 - t134 * t136) - g(2) * (t134 * t127 + t135 * t136) - g(3) * (t138 * pkin(2) + qJ(1)), 0, 0, 0, 0, 0, -g(1) * (t129 * t142 + t147) - g(2) * (t129 * t146 - t143) - t139 * t150, -g(1) * (-t129 * t143 + t146) - g(2) * (-t129 * t147 - t142) + t137 * t150, 0, 0, 0, 0, 0, -g(1) * (t129 * t144 + t149) - g(2) * (t129 * t148 - t145) - t131 * t150, -g(1) * (-t129 * t145 + t148) - g(2) * (-t129 * t149 - t144) + t130 * t150;];
U_reg = t1;
