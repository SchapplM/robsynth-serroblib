% Calculate minimal parameter regressor of potential energy for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRP6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:18
% EndTime: 2019-12-05 15:41:18
% DurationCPUTime: 0.08s
% Computational Cost: add. (62->37), mult. (127->49), div. (0->0), fcn. (132->6), ass. (0->23)
t118 = sin(qJ(2));
t120 = cos(qJ(2));
t131 = pkin(2) * t120 + qJ(3) * t118 + pkin(1);
t129 = g(3) * t120;
t117 = sin(qJ(4));
t127 = t117 * t118;
t119 = cos(qJ(4));
t126 = t118 * t119;
t125 = t118 * pkin(2) + qJ(1);
t115 = sin(pkin(7));
t116 = cos(pkin(7));
t124 = t115 * pkin(5) + t131 * t116;
t123 = g(1) * t116 + g(2) * t115;
t122 = -t116 * pkin(5) + t131 * t115;
t102 = t115 * t117 - t116 * t126;
t104 = t115 * t126 + t116 * t117;
t121 = g(1) * t102 - g(2) * t104 + t119 * t129;
t105 = t115 * t127 - t116 * t119;
t103 = t115 * t119 + t116 * t127;
t101 = g(3) * t118 + t123 * t120;
t100 = t123 * t118 - t129;
t99 = -g(1) * t103 - g(2) * t105 + t117 * t129;
t1 = [-g(3) * qJ(1), 0, -t101, t100, t101, -t100, -g(1) * t124 - g(2) * t122 - g(3) * (-t120 * qJ(3) + t125), 0, 0, 0, 0, 0, t99, t121, t99, -t101, -t121, -g(1) * (t115 * pkin(3) + t103 * pkin(4) + t102 * qJ(5) + t124) - g(2) * (-t116 * pkin(3) + t105 * pkin(4) - t104 * qJ(5) + t122) - g(3) * (t118 * pkin(6) + t125) + (-g(3) * (-pkin(4) * t117 + qJ(5) * t119 - qJ(3)) - t123 * pkin(6)) * t120;];
U_reg = t1;
