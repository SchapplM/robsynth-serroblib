% Calculate minimal parameter regressor of potential energy for
% S5PRPRP2
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
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:03
% EndTime: 2019-12-05 15:31:03
% DurationCPUTime: 0.08s
% Computational Cost: add. (80->36), mult. (77->48), div. (0->0), fcn. (74->8), ass. (0->21)
t125 = cos(qJ(4));
t115 = t125 * pkin(4) + pkin(3);
t121 = sin(pkin(8));
t122 = cos(pkin(8));
t123 = -qJ(5) - pkin(6);
t136 = t115 * t122 - t121 * t123;
t135 = g(3) * t121;
t134 = pkin(5) + qJ(1);
t120 = pkin(7) + qJ(2);
t116 = sin(t120);
t124 = sin(qJ(4));
t132 = t116 * t124;
t130 = t122 * t124;
t129 = t122 * t125;
t128 = t116 * pkin(2) + sin(pkin(7)) * pkin(1);
t117 = cos(t120);
t127 = t117 * pkin(2) + t116 * qJ(3) + cos(pkin(7)) * pkin(1);
t126 = g(1) * t117 + g(2) * t116;
t111 = g(1) * t116 - g(2) * t117;
t110 = -g(3) * t122 + t126 * t121;
t1 = [-g(3) * qJ(1), 0, -t126, t111, -t126 * t122 - t135, t110, -t111, -g(1) * t127 - g(2) * (-t117 * qJ(3) + t128) - g(3) * t134, 0, 0, 0, 0, 0, -g(1) * (t117 * t129 + t132) - g(2) * (t116 * t129 - t117 * t124) - t125 * t135, -g(1) * (t116 * t125 - t117 * t130) - g(2) * (-t116 * t130 - t117 * t125) + t124 * t135, -t110, -g(1) * (pkin(4) * t132 + t127) - g(2) * (t136 * t116 + t128) - g(3) * (t121 * t115 + t122 * t123 + t134) + (-g(1) * t136 - g(2) * (-pkin(4) * t124 - qJ(3))) * t117;];
U_reg = t1;
