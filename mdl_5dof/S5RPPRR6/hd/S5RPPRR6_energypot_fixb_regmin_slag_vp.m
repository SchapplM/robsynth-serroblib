% Calculate minimal parameter regressor of potential energy for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:07
% EndTime: 2019-12-31 17:58:07
% DurationCPUTime: 0.05s
% Computational Cost: add. (63->26), mult. (60->40), div. (0->0), fcn. (58->10), ass. (0->21)
t137 = pkin(9) + qJ(4);
t133 = sin(t137);
t153 = g(3) * t133;
t152 = g(3) * (qJ(2) + pkin(5));
t138 = qJ(1) + pkin(8);
t134 = sin(t138);
t142 = sin(qJ(5));
t151 = t134 * t142;
t144 = cos(qJ(5));
t150 = t134 * t144;
t136 = cos(t138);
t149 = t136 * t142;
t148 = t136 * t144;
t147 = g(1) * t136 + g(2) * t134;
t143 = sin(qJ(1));
t145 = cos(qJ(1));
t146 = -g(1) * t145 - g(2) * t143;
t140 = cos(pkin(9));
t139 = sin(pkin(9));
t135 = cos(t137);
t1 = [0, t146, g(1) * t143 - g(2) * t145, t146 * pkin(1) - t152, -g(3) * t139 - t147 * t140, -g(3) * t140 + t147 * t139, -g(1) * t134 + g(2) * t136, -g(1) * (t145 * pkin(1) + t136 * pkin(2) + t134 * qJ(3)) - g(2) * (t143 * pkin(1) + t134 * pkin(2) - t136 * qJ(3)) - t152, 0, 0, 0, 0, 0, -t147 * t135 - t153, -g(3) * t135 + t147 * t133, 0, 0, 0, 0, 0, -g(1) * (t135 * t148 + t151) - g(2) * (t135 * t150 - t149) - t144 * t153, -g(1) * (-t135 * t149 + t150) - g(2) * (-t135 * t151 - t148) + t142 * t153;];
U_reg = t1;
