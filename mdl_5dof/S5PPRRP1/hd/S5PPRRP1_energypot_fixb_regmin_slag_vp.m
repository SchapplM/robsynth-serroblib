% Calculate minimal parameter regressor of potential energy for
% S5PPRRP1
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
% U_reg [1x14]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:07:13
% EndTime: 2019-12-05 15:07:13
% DurationCPUTime: 0.06s
% Computational Cost: add. (63->32), mult. (70->42), div. (0->0), fcn. (67->8), ass. (0->20)
t123 = g(3) * qJ(1);
t108 = pkin(8) + qJ(3);
t106 = sin(t108);
t122 = g(3) * t106;
t109 = sin(pkin(7));
t113 = sin(qJ(4));
t121 = t109 * t113;
t114 = cos(qJ(4));
t120 = t109 * t114;
t110 = cos(pkin(7));
t119 = t110 * t113;
t118 = t110 * t114;
t117 = pkin(4) * t113 + pkin(5) + qJ(2);
t116 = g(1) * t110 + g(2) * t109;
t105 = t114 * pkin(4) + pkin(3);
t107 = cos(t108);
t111 = -qJ(5) - pkin(6);
t115 = t105 * t107 - t106 * t111 + cos(pkin(8)) * pkin(2) + pkin(1);
t103 = -g(3) * t107 + t116 * t106;
t1 = [-t123, -g(1) * (t110 * pkin(1) + t109 * qJ(2)) - g(2) * (t109 * pkin(1) - t110 * qJ(2)) - t123, 0, -t116 * t107 - t122, t103, 0, 0, 0, 0, 0, -g(1) * (t107 * t118 + t121) - g(2) * (t107 * t120 - t119) - t114 * t122, -g(1) * (-t107 * t119 + t120) - g(2) * (-t107 * t121 - t118) + t113 * t122, -t103, -g(3) * (t106 * t105 + t107 * t111 + sin(pkin(8)) * pkin(2) + qJ(1)) + (-g(1) * t115 + g(2) * t117) * t110 + (-g(1) * t117 - g(2) * t115) * t109;];
U_reg = t1;
