% Calculate minimal parameter regressor of potential energy for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:26
% EndTime: 2019-12-05 16:04:26
% DurationCPUTime: 0.08s
% Computational Cost: add. (70->42), mult. (169->73), div. (0->0), fcn. (207->10), ass. (0->27)
t151 = sin(pkin(5));
t168 = pkin(6) * t151;
t155 = sin(qJ(4));
t167 = t151 * t155;
t156 = sin(qJ(2));
t166 = t151 * t156;
t158 = cos(qJ(4));
t165 = t151 * t158;
t159 = cos(qJ(2));
t164 = t151 * t159;
t153 = cos(pkin(5));
t163 = t153 * t156;
t162 = t153 * t159;
t150 = sin(pkin(9));
t152 = cos(pkin(9));
t144 = t150 * t156 - t152 * t162;
t146 = t150 * t162 + t152 * t156;
t161 = -g(1) * t146 - g(2) * t144 + g(3) * t164;
t145 = t150 * t159 + t152 * t163;
t147 = -t150 * t163 + t152 * t159;
t160 = g(1) * t147 + g(2) * t145 + g(3) * t166;
t157 = cos(qJ(5));
t154 = sin(qJ(5));
t148 = t153 * t158 - t155 * t164;
t143 = t144 * t155 - t152 * t165;
t142 = t146 * t155 + t150 * t165;
t1 = [-g(3) * qJ(1), 0, -t160, -t161, t160, t161, -g(1) * (t152 * pkin(1) + t147 * pkin(2) + t146 * qJ(3) + t150 * t168) - g(2) * (t150 * pkin(1) + t145 * pkin(2) + t144 * qJ(3) - t152 * t168) - g(3) * (t153 * pkin(6) + qJ(1) + (pkin(2) * t156 - qJ(3) * t159) * t151), 0, 0, 0, 0, 0, -g(1) * t142 - g(2) * t143 - g(3) * t148, -g(1) * (t146 * t158 - t150 * t167) - g(2) * (t144 * t158 + t152 * t167) - g(3) * (-t153 * t155 - t158 * t164), 0, 0, 0, 0, 0, -g(1) * (t142 * t157 + t147 * t154) - g(2) * (t143 * t157 + t145 * t154) - g(3) * (t148 * t157 + t154 * t166), -g(1) * (-t142 * t154 + t147 * t157) - g(2) * (-t143 * t154 + t145 * t157) - g(3) * (-t148 * t154 + t157 * t166);];
U_reg = t1;
