% Calculate minimal parameter regressor of potential energy for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_energypot_fixb_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:03:22
% EndTime: 2019-12-05 17:03:22
% DurationCPUTime: 0.04s
% Computational Cost: add. (31->17), mult. (45->27), div. (0->0), fcn. (48->8), ass. (0->16)
t103 = qJ(3) + qJ(4);
t101 = sin(t103);
t115 = g(2) * t101;
t104 = sin(qJ(5));
t106 = sin(qJ(2));
t114 = t106 * t104;
t107 = cos(qJ(5));
t113 = t106 * t107;
t109 = cos(qJ(2));
t112 = t109 * t104;
t111 = t109 * t107;
t110 = g(1) * t109 + g(3) * t106;
t108 = cos(qJ(3));
t105 = sin(qJ(3));
t102 = cos(t103);
t1 = [-g(3) * qJ(1), 0, -t110, g(1) * t106 - g(3) * t109, 0, 0, 0, 0, 0, g(2) * t105 - t110 * t108, g(2) * t108 + t110 * t105, 0, 0, 0, 0, 0, -t110 * t102 + t115, g(2) * t102 + t110 * t101, 0, 0, 0, 0, 0, -g(1) * (t102 * t111 + t114) + t107 * t115 - g(3) * (t102 * t113 - t112), -g(1) * (-t102 * t112 + t113) - t104 * t115 - g(3) * (-t102 * t114 - t111);];
U_reg = t1;
