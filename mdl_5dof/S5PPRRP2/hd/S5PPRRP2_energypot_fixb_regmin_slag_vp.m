% Calculate minimal parameter regressor of potential energy for
% S5PPRRP2
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:09:09
% EndTime: 2019-12-05 15:09:09
% DurationCPUTime: 0.08s
% Computational Cost: add. (85->36), mult. (103->47), div. (0->0), fcn. (108->8), ass. (0->24)
t121 = pkin(8) + qJ(3);
t120 = cos(t121);
t136 = cos(pkin(8)) * pkin(2) + pkin(1) + pkin(3) * t120;
t134 = g(3) * qJ(1);
t119 = sin(t121);
t133 = g(3) * t119;
t122 = sin(pkin(7));
t125 = sin(qJ(4));
t132 = t122 * t125;
t126 = cos(qJ(4));
t131 = t122 * t126;
t123 = cos(pkin(7));
t130 = t123 * t125;
t129 = t123 * t126;
t128 = g(1) * t123 + g(2) * t122;
t113 = t120 * t132 + t129;
t115 = t120 * t130 - t131;
t127 = g(1) * t115 + g(2) * t113 + t125 * t133;
t124 = -pkin(5) - qJ(2);
t116 = t120 * t129 + t132;
t114 = t120 * t131 - t130;
t112 = -g(3) * t120 + t128 * t119;
t111 = -g(1) * t116 - g(2) * t114 - t126 * t133;
t1 = [-t134, -g(1) * (t123 * pkin(1) + t122 * qJ(2)) - g(2) * (t122 * pkin(1) - t123 * qJ(2)) - t134, 0, -t128 * t120 - t133, t112, 0, 0, 0, 0, 0, t111, t127, t111, -t112, -t127, -g(1) * (t116 * pkin(4) + t115 * qJ(5) - t122 * t124 + t136 * t123) - g(2) * (t114 * pkin(4) + t113 * qJ(5) + t136 * t122 + t123 * t124) - g(3) * (sin(pkin(8)) * pkin(2) - t120 * pkin(6) + qJ(1)) + (-g(3) * (pkin(4) * t126 + qJ(5) * t125 + pkin(3)) - t128 * pkin(6)) * t119;];
U_reg = t1;
