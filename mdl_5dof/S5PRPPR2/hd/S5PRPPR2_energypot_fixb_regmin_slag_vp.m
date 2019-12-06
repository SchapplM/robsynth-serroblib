% Calculate minimal parameter regressor of potential energy for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% U_reg [1x16]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:41
% EndTime: 2019-12-05 15:24:41
% DurationCPUTime: 0.08s
% Computational Cost: add. (78->37), mult. (86->56), div. (0->0), fcn. (87->10), ass. (0->29)
t131 = qJ(2) + pkin(8);
t126 = sin(t131);
t152 = g(3) * t126;
t130 = pkin(9) + qJ(5);
t125 = sin(t130);
t133 = sin(pkin(7));
t151 = t133 * t125;
t127 = cos(t130);
t150 = t133 * t127;
t132 = sin(pkin(9));
t149 = t133 * t132;
t134 = cos(pkin(9));
t148 = t133 * t134;
t135 = cos(pkin(7));
t147 = t135 * t125;
t146 = t135 * t127;
t145 = t135 * t132;
t144 = t135 * t134;
t138 = cos(qJ(2));
t124 = t138 * pkin(2) + pkin(1);
t136 = -qJ(3) - pkin(5);
t143 = t133 * t124 + t135 * t136;
t137 = sin(qJ(2));
t142 = t137 * pkin(2) + qJ(1);
t141 = t135 * t124 - t133 * t136;
t140 = g(1) * t135 + g(2) * t133;
t128 = cos(t131);
t139 = pkin(3) * t128 + qJ(4) * t126;
t1 = [-g(3) * qJ(1), 0, -g(3) * t137 - t140 * t138, -g(3) * t138 + t140 * t137, -g(1) * t141 - g(2) * t143 - g(3) * t142, -g(1) * (t128 * t144 + t149) - g(2) * (t128 * t148 - t145) - t134 * t152, -g(1) * (-t128 * t145 + t148) - g(2) * (-t128 * t149 - t144) + t132 * t152, g(3) * t128 - t140 * t126, -g(1) * (t139 * t135 + t141) - g(2) * (t139 * t133 + t143) - g(3) * (t126 * pkin(3) - t128 * qJ(4) + t142), 0, 0, 0, 0, 0, -g(1) * (t128 * t146 + t151) - g(2) * (t128 * t150 - t147) - t127 * t152, -g(1) * (-t128 * t147 + t150) - g(2) * (-t128 * t151 - t146) + t125 * t152;];
U_reg = t1;
