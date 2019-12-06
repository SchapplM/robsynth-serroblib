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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:29:05
% EndTime: 2019-12-05 17:29:05
% DurationCPUTime: 0.08s
% Computational Cost: add. (94->41), mult. (94->63), div. (0->0), fcn. (93->10), ass. (0->27)
t140 = sin(pkin(8));
t157 = g(1) * t140;
t143 = qJ(2) + pkin(5);
t156 = g(1) * t143;
t138 = qJ(1) + pkin(7);
t133 = sin(t138);
t155 = g(2) * t133;
t142 = cos(pkin(8));
t154 = t133 * t142;
t135 = cos(t138);
t153 = t135 * t142;
t139 = sin(pkin(9));
t152 = t139 * t142;
t141 = cos(pkin(9));
t151 = t141 * t142;
t145 = cos(qJ(1));
t150 = t145 * pkin(1) + t135 * pkin(2) + t133 * qJ(3);
t144 = sin(qJ(1));
t149 = -t144 * pkin(1) + t135 * qJ(3);
t148 = -g(3) * t135 + t155;
t147 = g(2) * t144 - g(3) * t145;
t146 = pkin(3) * t142 + qJ(4) * t140;
t137 = pkin(9) + qJ(5);
t134 = cos(t137);
t132 = sin(t137);
t128 = g(1) * t142 + t148 * t140;
t1 = [0, t147, g(2) * t145 + g(3) * t144, t147 * pkin(1) - t156, t148 * t142 - t157, -t128, -g(2) * t135 - g(3) * t133, -t156 - g(2) * (-t133 * pkin(2) + t149) - g(3) * t150, -t141 * t157 - g(2) * (-t133 * t151 + t135 * t139) - g(3) * (t133 * t139 + t135 * t151), t139 * t157 - g(2) * (t133 * t152 + t135 * t141) - g(3) * (t133 * t141 - t135 * t152), t128, -g(1) * (t140 * pkin(3) - t142 * qJ(4) + t143) - g(2) * t149 - g(3) * (t146 * t135 + t150) - (-pkin(2) - t146) * t155, 0, 0, 0, 0, 0, -t134 * t157 - g(2) * (t135 * t132 - t134 * t154) - g(3) * (t133 * t132 + t134 * t153), t132 * t157 - g(2) * (t132 * t154 + t135 * t134) - g(3) * (-t132 * t153 + t133 * t134);];
U_reg = t1;
