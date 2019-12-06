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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:36:10
% EndTime: 2019-12-05 18:36:10
% DurationCPUTime: 0.06s
% Computational Cost: add. (72->32), mult. (69->53), div. (0->0), fcn. (74->10), ass. (0->20)
t147 = sin(pkin(9));
t158 = g(1) * t147;
t146 = qJ(1) + qJ(2);
t142 = sin(t146);
t148 = cos(pkin(9));
t157 = t142 * t148;
t144 = cos(t146);
t156 = t144 * t148;
t149 = sin(qJ(4));
t155 = t148 * t149;
t151 = cos(qJ(4));
t154 = t148 * t151;
t153 = g(2) * t142 - g(3) * t144;
t152 = cos(qJ(1));
t150 = sin(qJ(1));
t145 = qJ(4) + qJ(5);
t143 = cos(t145);
t141 = sin(t145);
t140 = g(2) * t144 + g(3) * t142;
t1 = [0, g(2) * t150 - g(3) * t152, g(2) * t152 + g(3) * t150, 0, t153, t140, t153 * t148 - t158, -g(1) * t148 - t153 * t147, -t140, -g(1) * (pkin(6) + pkin(5)) - g(2) * (-t150 * pkin(1) - t142 * pkin(2) + t144 * qJ(3)) - g(3) * (t152 * pkin(1) + t144 * pkin(2) + t142 * qJ(3)), 0, 0, 0, 0, 0, -t151 * t158 - g(2) * (-t142 * t154 + t144 * t149) - g(3) * (t142 * t149 + t144 * t154), t149 * t158 - g(2) * (t142 * t155 + t144 * t151) - g(3) * (t142 * t151 - t144 * t155), 0, 0, 0, 0, 0, -t143 * t158 - g(2) * (t144 * t141 - t143 * t157) - g(3) * (t142 * t141 + t143 * t156), t141 * t158 - g(2) * (t141 * t157 + t144 * t143) - g(3) * (-t141 * t156 + t142 * t143);];
U_reg = t1;
