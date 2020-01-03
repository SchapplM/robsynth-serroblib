% Calculate minimal parameter regressor of potential energy for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:04
% EndTime: 2019-12-31 16:54:04
% DurationCPUTime: 0.07s
% Computational Cost: add. (35->21), mult. (53->33), div. (0->0), fcn. (54->8), ass. (0->17)
t113 = pkin(7) + qJ(3);
t111 = sin(t113);
t125 = g(3) * t111;
t116 = sin(qJ(4));
t117 = sin(qJ(1));
t124 = t117 * t116;
t118 = cos(qJ(4));
t123 = t117 * t118;
t119 = cos(qJ(1));
t122 = t119 * t116;
t121 = t119 * t118;
t120 = g(1) * t119 + g(2) * t117;
t115 = cos(pkin(7));
t114 = sin(pkin(7));
t112 = cos(t113);
t110 = g(1) * t117 - g(2) * t119;
t1 = [0, -t120, t110, -g(3) * t114 - t115 * t120, -g(3) * t115 + t114 * t120, -t110, -g(1) * (t119 * pkin(1) + t117 * qJ(2)) - g(2) * (t117 * pkin(1) - t119 * qJ(2)) - g(3) * pkin(4), 0, 0, 0, 0, 0, -t112 * t120 - t125, -g(3) * t112 + t111 * t120, 0, 0, 0, 0, 0, -g(1) * (t112 * t121 + t124) - g(2) * (t112 * t123 - t122) - t118 * t125, -g(1) * (-t112 * t122 + t123) - g(2) * (-t112 * t124 - t121) + t116 * t125;];
U_reg = t1;
