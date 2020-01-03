% Calculate minimal parameter regressor of potential energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:32
% EndTime: 2020-01-03 11:20:32
% DurationCPUTime: 0.10s
% Computational Cost: add. (94->43), mult. (94->66), div. (0->0), fcn. (93->10), ass. (0->26)
t136 = sin(pkin(8));
t152 = g(1) * t136;
t139 = qJ(2) + pkin(5);
t151 = g(1) * t139;
t150 = qJ(4) * t136;
t134 = qJ(1) + pkin(7);
t129 = sin(t134);
t138 = cos(pkin(8));
t149 = t129 * t138;
t131 = cos(t134);
t148 = t131 * t138;
t135 = sin(pkin(9));
t147 = t135 * t138;
t137 = cos(pkin(9));
t146 = t137 * t138;
t140 = sin(qJ(1));
t145 = t140 * pkin(1) + t129 * pkin(2);
t144 = -g(2) * t129 + g(3) * t131;
t141 = cos(qJ(1));
t143 = -g(2) * t140 + g(3) * t141;
t142 = -t141 * pkin(1) - t129 * qJ(3);
t133 = pkin(9) + qJ(5);
t130 = cos(t133);
t128 = sin(t133);
t126 = g(1) * t138 + t144 * t136;
t1 = [0, t143, -g(2) * t141 - g(3) * t140, t143 * pkin(1) - t151, t144 * t138 - t152, -t126, g(2) * t131 + g(3) * t129, -t151 - g(2) * (-t131 * qJ(3) + t145) - g(3) * (-t131 * pkin(2) + t142), -t137 * t152 - g(2) * (t129 * t146 - t131 * t135) - g(3) * (-t129 * t135 - t131 * t146), t135 * t152 - g(2) * (-t129 * t147 - t131 * t137) - g(3) * (-t129 * t137 + t131 * t147), t126, -g(1) * (t136 * pkin(3) - t138 * qJ(4) + t139) - g(2) * (pkin(3) * t149 + t129 * t150 + t145) - g(3) * t142 + (g(2) * qJ(3) - g(3) * (-pkin(3) * t138 - pkin(2) - t150)) * t131, 0, 0, 0, 0, 0, -t130 * t152 - g(2) * (-t131 * t128 + t130 * t149) - g(3) * (-t129 * t128 - t130 * t148), t128 * t152 - g(2) * (-t128 * t149 - t131 * t130) - g(3) * (t128 * t148 - t129 * t130);];
U_reg = t1;
