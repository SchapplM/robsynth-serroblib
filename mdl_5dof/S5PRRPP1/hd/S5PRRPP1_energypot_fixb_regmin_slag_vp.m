% Calculate minimal parameter regressor of potential energy for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:51
% EndTime: 2019-12-05 16:06:51
% DurationCPUTime: 0.08s
% Computational Cost: add. (83->29), mult. (63->33), div. (0->0), fcn. (56->8), ass. (0->17)
t124 = sin(qJ(3));
t130 = t124 * pkin(3) + pkin(5) + qJ(1);
t125 = cos(qJ(3));
t113 = t125 * pkin(3) + pkin(2);
t121 = pkin(7) + qJ(2);
t114 = sin(t121);
t116 = cos(t121);
t123 = -pkin(6) - qJ(4);
t129 = t114 * t113 + t116 * t123 + sin(pkin(7)) * pkin(1);
t128 = g(1) * t116 + g(2) * t114;
t127 = t116 * t113 - t114 * t123 + cos(pkin(7)) * pkin(1);
t122 = qJ(3) + pkin(8);
t115 = sin(t122);
t117 = cos(t122);
t126 = pkin(4) * t117 + qJ(5) * t115;
t108 = g(1) * t114 - g(2) * t116;
t1 = [-g(3) * qJ(1), 0, -t128, t108, 0, 0, 0, 0, 0, -g(3) * t124 - t128 * t125, -g(3) * t125 + t128 * t124, -t108, -g(1) * t127 - g(2) * t129 - g(3) * t130, -g(3) * t115 - t128 * t117, -t108, g(3) * t117 - t128 * t115, -g(1) * (t126 * t116 + t127) - g(2) * (t126 * t114 + t129) - g(3) * (t115 * pkin(4) - t117 * qJ(5) + t130);];
U_reg = t1;
