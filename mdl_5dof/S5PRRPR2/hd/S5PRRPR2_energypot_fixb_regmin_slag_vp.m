% Calculate minimal parameter regressor of potential energy for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:34
% EndTime: 2019-12-05 16:17:34
% DurationCPUTime: 0.05s
% Computational Cost: add. (78->29), mult. (52->40), div. (0->0), fcn. (52->10), ass. (0->16)
t108 = sin(pkin(9));
t115 = g(3) * t108;
t109 = cos(pkin(9));
t110 = sin(qJ(5));
t114 = t109 * t110;
t111 = cos(qJ(5));
t113 = t109 * t111;
t107 = pkin(8) + qJ(2);
t106 = qJ(3) + t107;
t102 = sin(t106);
t103 = cos(t106);
t112 = g(1) * t103 + g(2) * t102;
t105 = cos(t107);
t104 = sin(t107);
t101 = g(1) * t102 - g(2) * t103;
t1 = [-g(3) * qJ(1), 0, -g(1) * t105 - g(2) * t104, g(1) * t104 - g(2) * t105, 0, -t112, t101, -t112 * t109 - t115, -g(3) * t109 + t112 * t108, -t101, -g(1) * (t103 * pkin(3) + t102 * qJ(4) + pkin(2) * t105 + cos(pkin(8)) * pkin(1)) - g(2) * (t102 * pkin(3) - t103 * qJ(4) + pkin(2) * t104 + sin(pkin(8)) * pkin(1)) - g(3) * (pkin(6) + pkin(5) + qJ(1)), 0, 0, 0, 0, 0, -g(1) * (t102 * t110 + t103 * t113) - g(2) * (t102 * t113 - t103 * t110) - t111 * t115, -g(1) * (t102 * t111 - t103 * t114) - g(2) * (-t102 * t114 - t103 * t111) + t110 * t115;];
U_reg = t1;
