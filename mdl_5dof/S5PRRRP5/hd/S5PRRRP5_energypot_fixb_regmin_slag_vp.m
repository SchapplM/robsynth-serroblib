% Calculate minimal parameter regressor of potential energy for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:17
% EndTime: 2019-12-05 16:49:17
% DurationCPUTime: 0.08s
% Computational Cost: add. (64->35), mult. (85->54), div. (0->0), fcn. (89->8), ass. (0->21)
t130 = sin(qJ(2));
t140 = g(3) * t130;
t126 = qJ(3) + qJ(4);
t123 = sin(t126);
t129 = sin(qJ(3));
t139 = t129 * pkin(3) + pkin(4) * t123 + pkin(5);
t127 = sin(pkin(8));
t132 = cos(qJ(2));
t138 = t127 * t132;
t128 = cos(pkin(8));
t137 = t128 * t132;
t136 = t129 * t132;
t131 = cos(qJ(3));
t135 = t131 * t132;
t134 = g(1) * t128 + g(2) * t127;
t124 = cos(t126);
t121 = t131 * pkin(3) + pkin(4) * t124 + pkin(2);
t125 = -qJ(5) - pkin(7) - pkin(6);
t133 = t121 * t132 - t125 * t130 + pkin(1);
t120 = -g(3) * t132 + t134 * t130;
t1 = [-g(3) * qJ(1), 0, -t134 * t132 - t140, t120, 0, 0, 0, 0, 0, -g(1) * (t127 * t129 + t128 * t135) - g(2) * (t127 * t135 - t128 * t129) - t131 * t140, -g(1) * (t127 * t131 - t128 * t136) - g(2) * (-t127 * t136 - t128 * t131) + t129 * t140, 0, 0, 0, 0, 0, -g(1) * (t127 * t123 + t124 * t137) - g(2) * (-t128 * t123 + t124 * t138) - t124 * t140, -g(1) * (-t123 * t137 + t127 * t124) - g(2) * (-t123 * t138 - t128 * t124) + t123 * t140, -t120, -g(3) * (t130 * t121 + t132 * t125 + qJ(1)) + (-g(1) * t133 + g(2) * t139) * t128 + (-g(1) * t139 - g(2) * t133) * t127;];
U_reg = t1;
