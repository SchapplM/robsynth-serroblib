% Calculate minimal parameter regressor of potential energy for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% U_reg [1x14]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:33:40
% EndTime: 2019-12-05 15:33:40
% DurationCPUTime: 0.07s
% Computational Cost: add. (62->33), mult. (73->47), div. (0->0), fcn. (70->8), ass. (0->24)
t124 = cos(qJ(4));
t112 = t124 * pkin(4) + pkin(3);
t117 = qJ(2) + pkin(8);
t114 = sin(t117);
t115 = cos(t117);
t120 = -qJ(5) - pkin(6);
t136 = t112 * t115 - t114 * t120;
t135 = g(3) * t114;
t118 = sin(pkin(7));
t122 = sin(qJ(4));
t132 = t118 * t122;
t131 = t118 * t124;
t119 = cos(pkin(7));
t130 = t119 * t122;
t129 = t119 * t124;
t125 = cos(qJ(2));
t113 = t125 * pkin(2) + pkin(1);
t121 = -qJ(3) - pkin(5);
t128 = t118 * t113 + t119 * t121;
t123 = sin(qJ(2));
t127 = t123 * pkin(2) + qJ(1);
t126 = g(1) * t119 + g(2) * t118;
t110 = t119 * t113;
t1 = [-g(3) * qJ(1), 0, -g(3) * t123 - t126 * t125, -g(3) * t125 + t126 * t123, -g(1) * (-t118 * t121 + t110) - g(2) * t128 - g(3) * t127, 0, 0, 0, 0, 0, -g(1) * (t115 * t129 + t132) - g(2) * (t115 * t131 - t130) - t124 * t135, -g(1) * (-t115 * t130 + t131) - g(2) * (-t115 * t132 - t129) + t122 * t135, g(3) * t115 - t126 * t114, -g(1) * (t136 * t119 + t110) - g(2) * (-pkin(4) * t130 + t128) - g(3) * (t114 * t112 + t115 * t120 + t127) + (-g(1) * (pkin(4) * t122 - t121) - g(2) * t136) * t118;];
U_reg = t1;
