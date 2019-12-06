% Calculate minimal parameter regressor of potential energy for
% S5PRRRP6
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
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:16
% EndTime: 2019-12-05 16:52:17
% DurationCPUTime: 0.08s
% Computational Cost: add. (89->40), mult. (120->60), div. (0->0), fcn. (132->8), ass. (0->27)
t139 = sin(qJ(2));
t151 = g(3) * t139;
t136 = sin(pkin(8));
t141 = cos(qJ(2));
t150 = t136 * t141;
t137 = cos(pkin(8));
t149 = t137 * t141;
t138 = sin(qJ(3));
t148 = t138 * t141;
t140 = cos(qJ(3));
t147 = t140 * t141;
t146 = pkin(3) * t138 + pkin(5);
t145 = g(1) * t137 + g(2) * t136;
t132 = t140 * pkin(3) + pkin(2);
t142 = -pkin(7) - pkin(6);
t144 = t132 * t141 - t139 * t142 + pkin(1);
t135 = qJ(3) + qJ(4);
t133 = sin(t135);
t134 = cos(t135);
t126 = t133 * t150 + t137 * t134;
t128 = t133 * t149 - t136 * t134;
t143 = g(1) * t128 + g(2) * t126 + t133 * t151;
t130 = -g(3) * t141 + t145 * t139;
t129 = t136 * t133 + t134 * t149;
t127 = -t137 * t133 + t134 * t150;
t125 = -g(1) * t129 - g(2) * t127 - t134 * t151;
t1 = [-g(3) * qJ(1), 0, -t145 * t141 - t151, t130, 0, 0, 0, 0, 0, -g(1) * (t136 * t138 + t137 * t147) - g(2) * (t136 * t147 - t137 * t138) - t140 * t151, -g(1) * (t136 * t140 - t137 * t148) - g(2) * (-t136 * t148 - t137 * t140) + t138 * t151, 0, 0, 0, 0, 0, t125, t143, t125, -t130, -t143, -g(1) * (t129 * pkin(4) + t128 * qJ(5)) - g(2) * (t127 * pkin(4) + t126 * qJ(5)) - g(3) * (t141 * t142 + qJ(1)) - (pkin(4) * t134 + qJ(5) * t133 + t132) * t151 + (-g(1) * t144 + g(2) * t146) * t137 + (-g(1) * t146 - g(2) * t144) * t136;];
U_reg = t1;
