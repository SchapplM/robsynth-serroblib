% Calculate minimal parameter regressor of potential energy for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% U_reg [1x13]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPPRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:58:06
% EndTime: 2019-12-05 14:58:06
% DurationCPUTime: 0.06s
% Computational Cost: add. (57->31), mult. (81->52), div. (0->0), fcn. (88->8), ass. (0->22)
t121 = g(3) * qJ(1);
t106 = sin(pkin(8));
t120 = g(3) * t106;
t110 = sin(qJ(5));
t119 = t106 * t110;
t111 = cos(qJ(5));
t118 = t106 * t111;
t107 = sin(pkin(7));
t108 = cos(pkin(8));
t117 = t107 * t108;
t105 = pkin(9) + qJ(4);
t100 = sin(t105);
t109 = cos(pkin(7));
t116 = t109 * t100;
t101 = cos(t105);
t115 = t109 * t101;
t114 = t109 * pkin(1) + t107 * qJ(2);
t113 = t107 * pkin(1) - t109 * qJ(2);
t112 = pkin(2) * t108 + qJ(3) * t106;
t99 = t107 * t100 + t108 * t115;
t98 = t101 * t117 - t116;
t1 = [-t121, -g(1) * t114 - g(2) * t113 - t121, -g(1) * (t112 * t109 + t114) - g(2) * (t112 * t107 + t113) - g(3) * (t106 * pkin(2) - t108 * qJ(3) + qJ(1)), 0, -g(1) * t99 - g(2) * t98 - t101 * t120, -g(1) * (t107 * t101 - t108 * t116) - g(2) * (-t100 * t117 - t115) + t100 * t120, 0, 0, 0, 0, 0, -g(1) * (t109 * t119 + t99 * t111) - g(2) * (t107 * t119 + t98 * t111) - g(3) * (t101 * t118 - t108 * t110), -g(1) * (t109 * t118 - t99 * t110) - g(2) * (t107 * t118 - t98 * t110) - g(3) * (-t101 * t119 - t108 * t111);];
U_reg = t1;
