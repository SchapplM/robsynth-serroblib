% Calculate minimal parameter regressor of potential energy for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:00
% EndTime: 2019-12-31 17:44:00
% DurationCPUTime: 0.07s
% Computational Cost: add. (77->27), mult. (87->39), div. (0->0), fcn. (84->8), ass. (0->22)
t123 = qJ(1) + pkin(7);
t119 = sin(t123);
t120 = cos(t123);
t136 = g(1) * t120 + g(2) * t119;
t126 = qJ(2) + pkin(5);
t138 = g(3) * t126;
t130 = cos(qJ(1));
t137 = t130 * pkin(1) + t120 * pkin(2) + t119 * qJ(3);
t128 = sin(qJ(1));
t135 = -g(1) * t130 - g(2) * t128;
t134 = t128 * pkin(1) + t119 * pkin(2) - t120 * qJ(3);
t124 = sin(pkin(8));
t125 = cos(pkin(8));
t133 = pkin(3) * t125 + qJ(4) * t124;
t127 = sin(qJ(5));
t129 = cos(qJ(5));
t132 = t124 * t129 - t125 * t127;
t131 = t124 * t127 + t125 * t129;
t115 = -g(1) * t119 + g(2) * t120;
t114 = -g(3) * t124 - t136 * t125;
t113 = -g(3) * t125 + t136 * t124;
t1 = [0, t135, g(1) * t128 - g(2) * t130, t135 * pkin(1) - t138, t114, t113, t115, -g(1) * t137 - g(2) * t134 - t138, t114, t115, -t113, -g(1) * (t133 * t120 + t137) - g(2) * (t133 * t119 + t134) - g(3) * (t124 * pkin(3) - t125 * qJ(4) + t126), 0, 0, 0, 0, 0, -g(3) * t132 - t136 * t131, g(3) * t131 - t136 * t132;];
U_reg = t1;
