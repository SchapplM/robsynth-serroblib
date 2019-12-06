% Calculate minimal parameter regressor of potential energy for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:12:58
% EndTime: 2019-12-05 17:12:58
% DurationCPUTime: 0.07s
% Computational Cost: add. (59->30), mult. (71->53), div. (0->0), fcn. (82->10), ass. (0->19)
t146 = sin(qJ(2));
t154 = g(3) * t146;
t143 = sin(pkin(9));
t148 = cos(qJ(2));
t153 = t143 * t148;
t144 = cos(pkin(9));
t152 = t144 * t148;
t145 = sin(qJ(3));
t151 = t145 * t148;
t147 = cos(qJ(3));
t150 = t147 * t148;
t142 = qJ(3) + qJ(4);
t149 = g(1) * t144 + g(2) * t143;
t141 = qJ(5) + t142;
t140 = cos(t142);
t139 = sin(t142);
t138 = cos(t141);
t137 = sin(t141);
t1 = [-g(3) * qJ(1), 0, -t149 * t148 - t154, -g(3) * t148 + t149 * t146, 0, 0, 0, 0, 0, -g(1) * (t143 * t145 + t144 * t150) - g(2) * (t143 * t150 - t144 * t145) - t147 * t154, -g(1) * (t143 * t147 - t144 * t151) - g(2) * (-t143 * t151 - t144 * t147) + t145 * t154, 0, 0, 0, 0, 0, -g(1) * (t143 * t139 + t140 * t152) - g(2) * (-t144 * t139 + t140 * t153) - t140 * t154, -g(1) * (-t139 * t152 + t143 * t140) - g(2) * (-t139 * t153 - t144 * t140) + t139 * t154, 0, 0, 0, 0, 0, -g(1) * (t143 * t137 + t138 * t152) - g(2) * (-t144 * t137 + t138 * t153) - t138 * t154, -g(1) * (-t137 * t152 + t143 * t138) - g(2) * (-t137 * t153 - t144 * t138) + t137 * t154;];
U_reg = t1;
