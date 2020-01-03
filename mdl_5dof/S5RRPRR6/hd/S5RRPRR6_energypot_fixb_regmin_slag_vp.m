% Calculate minimal parameter regressor of potential energy for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:54
% EndTime: 2020-01-03 12:05:54
% DurationCPUTime: 0.06s
% Computational Cost: add. (72->33), mult. (69->53), div. (0->0), fcn. (74->10), ass. (0->20)
t145 = sin(pkin(9));
t156 = g(1) * t145;
t144 = qJ(1) + qJ(2);
t140 = sin(t144);
t146 = cos(pkin(9));
t155 = t140 * t146;
t142 = cos(t144);
t154 = t142 * t146;
t147 = sin(qJ(4));
t153 = t146 * t147;
t149 = cos(qJ(4));
t152 = t146 * t149;
t151 = g(2) * t140 - g(3) * t142;
t150 = cos(qJ(1));
t148 = sin(qJ(1));
t143 = qJ(4) + qJ(5);
t141 = cos(t143);
t139 = sin(t143);
t138 = g(2) * t142 + g(3) * t140;
t1 = [0, -g(2) * t148 + g(3) * t150, -g(2) * t150 - g(3) * t148, 0, -t151, -t138, -t151 * t146 - t156, -g(1) * t146 + t151 * t145, t138, -g(1) * (pkin(6) + pkin(5)) - g(2) * (t148 * pkin(1) + t140 * pkin(2) - t142 * qJ(3)) - g(3) * (-t150 * pkin(1) - t142 * pkin(2) - t140 * qJ(3)), 0, 0, 0, 0, 0, -t149 * t156 - g(2) * (t140 * t152 - t142 * t147) - g(3) * (-t140 * t147 - t142 * t152), t147 * t156 - g(2) * (-t140 * t153 - t142 * t149) - g(3) * (-t140 * t149 + t142 * t153), 0, 0, 0, 0, 0, -t141 * t156 - g(2) * (-t142 * t139 + t141 * t155) - g(3) * (-t140 * t139 - t141 * t154), t139 * t156 - g(2) * (-t139 * t155 - t142 * t141) - g(3) * (t139 * t154 - t140 * t141);];
U_reg = t1;
