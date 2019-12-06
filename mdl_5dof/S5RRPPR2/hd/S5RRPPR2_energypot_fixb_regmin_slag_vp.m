% Calculate minimal parameter regressor of potential energy for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:21
% EndTime: 2019-12-05 18:20:21
% DurationCPUTime: 0.06s
% Computational Cost: add. (77->30), mult. (58->45), div. (0->0), fcn. (56->10), ass. (0->20)
t139 = g(1) * (qJ(3) + pkin(6) + pkin(5));
t127 = sin(pkin(9));
t138 = g(1) * t127;
t128 = cos(pkin(9));
t129 = sin(qJ(5));
t137 = t128 * t129;
t131 = cos(qJ(5));
t136 = t128 * t131;
t126 = qJ(1) + qJ(2);
t123 = cos(t126);
t132 = cos(qJ(1));
t135 = t132 * pkin(1) + pkin(2) * t123;
t122 = sin(t126);
t130 = sin(qJ(1));
t134 = -t130 * pkin(1) - pkin(2) * t122;
t121 = pkin(8) + t126;
t118 = sin(t121);
t119 = cos(t121);
t133 = g(2) * t118 - g(3) * t119;
t1 = [0, g(2) * t130 - g(3) * t132, g(2) * t132 + g(3) * t130, 0, g(2) * t122 - g(3) * t123, g(2) * t123 + g(3) * t122, -g(2) * t134 - g(3) * t135 - t139, t133 * t128 - t138, -g(1) * t128 - t133 * t127, -g(2) * t119 - g(3) * t118, -t139 - g(2) * (-t118 * pkin(3) + t119 * qJ(4) + t134) - g(3) * (t119 * pkin(3) + t118 * qJ(4) + t135), 0, 0, 0, 0, 0, -t131 * t138 - g(2) * (-t118 * t136 + t119 * t129) - g(3) * (t118 * t129 + t119 * t136), t129 * t138 - g(2) * (t118 * t137 + t119 * t131) - g(3) * (t118 * t131 - t119 * t137);];
U_reg = t1;
