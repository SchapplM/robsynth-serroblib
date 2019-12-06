% Calculate minimal parameter regressor of potential energy for
% S5PRPPR1
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
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:08
% EndTime: 2019-12-05 15:22:08
% DurationCPUTime: 0.08s
% Computational Cost: add. (96->40), mult. (90->59), div. (0->0), fcn. (91->10), ass. (0->23)
t129 = sin(pkin(8));
t141 = g(3) * t129;
t140 = pkin(5) + qJ(1);
t127 = pkin(7) + qJ(2);
t121 = sin(t127);
t131 = cos(pkin(8));
t139 = t121 * t131;
t123 = cos(t127);
t138 = t123 * t131;
t128 = sin(pkin(9));
t137 = t128 * t131;
t130 = cos(pkin(9));
t136 = t130 * t131;
t135 = t123 * pkin(2) + t121 * qJ(3) + cos(pkin(7)) * pkin(1);
t134 = g(1) * t123 + g(2) * t121;
t133 = t121 * pkin(2) - t123 * qJ(3) + sin(pkin(7)) * pkin(1);
t132 = pkin(3) * t131 + qJ(4) * t129;
t126 = pkin(9) + qJ(5);
t122 = cos(t126);
t120 = sin(t126);
t116 = g(1) * t121 - g(2) * t123;
t115 = -g(3) * t131 + t134 * t129;
t1 = [-g(3) * qJ(1), 0, -t134, t116, -t134 * t131 - t141, t115, -t116, -g(1) * t135 - g(2) * t133 - g(3) * t140, -g(1) * (t121 * t128 + t123 * t136) - g(2) * (t121 * t136 - t123 * t128) - t130 * t141, -g(1) * (t121 * t130 - t123 * t137) - g(2) * (-t121 * t137 - t123 * t130) + t128 * t141, -t115, -g(1) * (t132 * t123 + t135) - g(2) * (t132 * t121 + t133) - g(3) * (t129 * pkin(3) - t131 * qJ(4) + t140), 0, 0, 0, 0, 0, -g(1) * (t121 * t120 + t122 * t138) - g(2) * (-t123 * t120 + t122 * t139) - t122 * t141, -g(1) * (-t120 * t138 + t121 * t122) - g(2) * (-t120 * t139 - t123 * t122) + t120 * t141;];
U_reg = t1;
