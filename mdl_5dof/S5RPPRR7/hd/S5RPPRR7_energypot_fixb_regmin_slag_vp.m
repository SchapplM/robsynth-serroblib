% Calculate minimal parameter regressor of potential energy for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:52
% EndTime: 2019-12-31 17:59:53
% DurationCPUTime: 0.05s
% Computational Cost: add. (46->23), mult. (52->38), div. (0->0), fcn. (50->8), ass. (0->16)
t125 = g(3) * (qJ(2) + pkin(5));
t118 = cos(qJ(4));
t124 = g(3) * t118;
t114 = sin(qJ(5));
t115 = sin(qJ(4));
t123 = t114 * t115;
t117 = cos(qJ(5));
t122 = t115 * t117;
t112 = qJ(1) + pkin(8);
t110 = sin(t112);
t111 = cos(t112);
t121 = -g(1) * t110 + g(2) * t111;
t116 = sin(qJ(1));
t119 = cos(qJ(1));
t120 = -g(1) * t119 - g(2) * t116;
t1 = [0, t120, g(1) * t116 - g(2) * t119, t120 * pkin(1) - t125, g(1) * t111 + g(2) * t110, t121, -g(1) * (t119 * pkin(1) + t111 * pkin(2) + t110 * qJ(3)) - g(2) * (t116 * pkin(1) + t110 * pkin(2) - t111 * qJ(3)) - t125, 0, 0, 0, 0, 0, t121 * t115 - t124, g(3) * t115 + t121 * t118, 0, 0, 0, 0, 0, -g(1) * (t110 * t122 + t111 * t114) - g(2) * (t110 * t114 - t111 * t122) - t117 * t124, -g(1) * (-t110 * t123 + t111 * t117) - g(2) * (t110 * t117 + t111 * t123) + t114 * t124;];
U_reg = t1;
