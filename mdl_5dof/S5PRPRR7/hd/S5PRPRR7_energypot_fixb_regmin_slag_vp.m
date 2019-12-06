% Calculate minimal parameter regressor of potential energy for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:00:44
% EndTime: 2019-12-05 16:00:44
% DurationCPUTime: 0.07s
% Computational Cost: add. (45->31), mult. (78->50), div. (0->0), fcn. (82->8), ass. (0->19)
t118 = cos(qJ(2));
t125 = g(3) * t118;
t113 = sin(pkin(8));
t116 = sin(qJ(2));
t124 = t113 * t116;
t114 = cos(pkin(8));
t123 = t114 * t116;
t115 = sin(qJ(4));
t122 = t115 * t116;
t117 = cos(qJ(4));
t121 = t116 * t117;
t120 = g(1) * t114 + g(2) * t113;
t119 = pkin(2) * t118 + qJ(3) * t116 + pkin(1);
t112 = qJ(4) + qJ(5);
t111 = cos(t112);
t110 = sin(t112);
t109 = g(3) * t116 + t120 * t118;
t108 = t120 * t116 - t125;
t1 = [-g(3) * qJ(1), 0, -t109, t108, t109, -t108, -g(3) * (t116 * pkin(2) - t118 * qJ(3) + qJ(1)) + (g(2) * pkin(5) - g(1) * t119) * t114 + (-g(1) * pkin(5) - g(2) * t119) * t113, 0, 0, 0, 0, 0, -g(1) * (t113 * t117 + t114 * t122) - g(2) * (t113 * t122 - t114 * t117) + t115 * t125, -g(1) * (-t113 * t115 + t114 * t121) - g(2) * (t113 * t121 + t114 * t115) + t117 * t125, 0, 0, 0, 0, 0, -g(1) * (t110 * t123 + t113 * t111) - g(2) * (t110 * t124 - t114 * t111) + t110 * t125, -g(1) * (-t113 * t110 + t111 * t123) - g(2) * (t114 * t110 + t111 * t124) + t111 * t125;];
U_reg = t1;
