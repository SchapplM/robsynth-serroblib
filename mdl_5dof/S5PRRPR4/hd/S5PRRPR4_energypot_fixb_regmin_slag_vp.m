% Calculate minimal parameter regressor of potential energy for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:23:21
% EndTime: 2019-12-05 16:23:21
% DurationCPUTime: 0.07s
% Computational Cost: add. (61->33), mult. (80->52), div. (0->0), fcn. (84->8), ass. (0->21)
t147 = sin(qJ(2));
t157 = g(3) * t147;
t143 = sin(pkin(8));
t149 = cos(qJ(2));
t156 = t143 * t149;
t144 = cos(pkin(8));
t155 = t144 * t149;
t146 = sin(qJ(3));
t154 = t146 * t149;
t148 = cos(qJ(3));
t153 = t148 * t149;
t152 = pkin(3) * t146 + pkin(5);
t151 = g(1) * t144 + g(2) * t143;
t141 = t148 * pkin(3) + pkin(2);
t145 = -qJ(4) - pkin(6);
t150 = t141 * t149 - t145 * t147 + pkin(1);
t142 = qJ(3) + pkin(9) + qJ(5);
t140 = cos(t142);
t139 = sin(t142);
t138 = -g(3) * t149 + t151 * t147;
t1 = [-g(3) * qJ(1), 0, -t151 * t149 - t157, t138, 0, 0, 0, 0, 0, -g(1) * (t143 * t146 + t144 * t153) - g(2) * (t143 * t153 - t144 * t146) - t148 * t157, -g(1) * (t143 * t148 - t144 * t154) - g(2) * (-t143 * t154 - t144 * t148) + t146 * t157, -t138, -g(3) * (t147 * t141 + t149 * t145 + qJ(1)) + (-g(1) * t150 + g(2) * t152) * t144 + (-g(1) * t152 - g(2) * t150) * t143, 0, 0, 0, 0, 0, -g(1) * (t143 * t139 + t140 * t155) - g(2) * (-t144 * t139 + t140 * t156) - t140 * t157, -g(1) * (-t139 * t155 + t143 * t140) - g(2) * (-t139 * t156 - t144 * t140) + t139 * t157;];
U_reg = t1;
