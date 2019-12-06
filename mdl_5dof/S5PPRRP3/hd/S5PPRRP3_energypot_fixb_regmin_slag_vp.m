% Calculate minimal parameter regressor of potential energy for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% U_reg [1x16]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:15
% EndTime: 2019-12-05 15:11:15
% DurationCPUTime: 0.11s
% Computational Cost: add. (81->45), mult. (183->67), div. (0->0), fcn. (216->8), ass. (0->33)
t142 = cos(pkin(8));
t160 = pkin(2) * t142;
t159 = g(3) * qJ(1);
t140 = sin(pkin(8));
t144 = sin(qJ(4));
t158 = t140 * t144;
t146 = cos(qJ(4));
t157 = t140 * t146;
t147 = cos(qJ(3));
t156 = t140 * t147;
t141 = sin(pkin(7));
t145 = sin(qJ(3));
t155 = t141 * t145;
t154 = t141 * t147;
t143 = cos(pkin(7));
t153 = t143 * t145;
t152 = t143 * t147;
t151 = t143 * pkin(1) + t141 * qJ(2);
t150 = t141 * pkin(1) - qJ(2) * t143;
t129 = t142 * t154 - t153;
t124 = t129 * t144 - t141 * t157;
t131 = t142 * t152 + t155;
t126 = t131 * t144 - t143 * t157;
t132 = t142 * t146 + t144 * t156;
t149 = g(1) * t126 + g(2) * t124 + g(3) * t132;
t128 = t142 * t155 + t152;
t130 = t142 * t153 - t154;
t148 = g(3) * t140 * t145 + g(1) * t130 + g(2) * t128;
t133 = -t142 * t144 + t146 * t156;
t127 = t131 * t146 + t143 * t158;
t125 = t129 * t146 + t141 * t158;
t123 = -g(1) * t127 - g(2) * t125 - g(3) * t133;
t1 = [-t159, -g(1) * t151 - g(2) * t150 - t159, 0, -g(1) * t131 - g(2) * t129 - g(3) * t156, t148, 0, 0, 0, 0, 0, t123, t149, t123, -t148, -t149, -g(1) * (pkin(3) * t131 + pkin(4) * t127 + pkin(6) * t130 + qJ(5) * t126 + t143 * t160 + t151) - g(2) * (pkin(3) * t129 + pkin(4) * t125 + pkin(6) * t128 + qJ(5) * t124 + t141 * t160 + t150) - g(3) * (t133 * pkin(4) - t142 * pkin(5) + t132 * qJ(5) + qJ(1)) + (-g(3) * (pkin(3) * t147 + pkin(6) * t145 + pkin(2)) + (-g(1) * t143 - g(2) * t141) * pkin(5)) * t140;];
U_reg = t1;
