% Calculate minimal parameter regressor of potential energy for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP13_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP13_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:59:36
% EndTime: 2019-12-31 18:59:36
% DurationCPUTime: 0.06s
% Computational Cost: add. (53->36), mult. (109->47), div. (0->0), fcn. (115->6), ass. (0->24)
t126 = sin(qJ(3));
t139 = pkin(3) * t126;
t129 = cos(qJ(3));
t138 = g(3) * t129;
t125 = sin(qJ(4));
t127 = sin(qJ(1));
t137 = t127 * t125;
t128 = cos(qJ(4));
t136 = t127 * t128;
t130 = cos(qJ(1));
t135 = t130 * t125;
t134 = t130 * t128;
t133 = t130 * pkin(1) + t127 * qJ(2);
t132 = t127 * pkin(1) - t130 * qJ(2);
t119 = g(1) * t127 - g(2) * t130;
t115 = t126 * t137 - t134;
t117 = t126 * t135 + t136;
t131 = g(1) * t115 - g(2) * t117 + t125 * t138;
t120 = g(1) * t130 + g(2) * t127;
t118 = -t126 * t134 + t137;
t116 = t126 * t136 + t135;
t114 = -g(3) * t126 + t119 * t129;
t113 = -g(1) * t116 - g(2) * t118 - t128 * t138;
t1 = [0, -t120, t119, t120, -t119, -g(3) * pkin(5) - g(1) * t133 - g(2) * t132, 0, 0, 0, 0, 0, -t119 * t126 - t138, -t114, 0, 0, 0, 0, 0, t113, t131, t113, t114, -t131, -g(1) * (t116 * pkin(4) + t130 * pkin(6) + t115 * qJ(5) + t127 * t139 + t133) - g(2) * (t118 * pkin(4) + t127 * pkin(6) - t117 * qJ(5) - t130 * t139 + t132) - g(3) * (t126 * pkin(7) + pkin(2) + pkin(5)) + (-g(3) * (pkin(4) * t128 + qJ(5) * t125 + pkin(3)) + t119 * pkin(7)) * t129;];
U_reg = t1;
