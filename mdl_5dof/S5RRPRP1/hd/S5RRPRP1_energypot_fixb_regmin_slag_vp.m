% Calculate minimal parameter regressor of potential energy for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x16]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:22:22
% EndTime: 2019-12-05 18:22:22
% DurationCPUTime: 0.04s
% Computational Cost: add. (58->24), mult. (41->32), div. (0->0), fcn. (35->8), ass. (0->17)
t109 = qJ(1) + qJ(2);
t107 = cos(t109);
t114 = cos(qJ(1));
t118 = t114 * pkin(1) + pkin(2) * t107;
t117 = qJ(3) + pkin(6) + pkin(5);
t106 = sin(t109);
t112 = sin(qJ(1));
t116 = -t112 * pkin(1) - pkin(2) * t106;
t105 = pkin(8) + t109;
t101 = sin(t105);
t102 = cos(t105);
t115 = g(2) * t101 - g(3) * t102;
t113 = cos(qJ(4));
t111 = sin(qJ(4));
t110 = -qJ(5) - pkin(7);
t104 = t113 * pkin(4) + pkin(3);
t1 = [0, g(2) * t112 - g(3) * t114, g(2) * t114 + g(3) * t112, 0, g(2) * t106 - g(3) * t107, g(2) * t107 + g(3) * t106, -g(1) * t117 - g(2) * t116 - g(3) * t118, 0, 0, 0, 0, 0, -g(1) * t111 + t115 * t113, -g(1) * t113 - t115 * t111, -g(2) * t102 - g(3) * t101, -g(1) * (t111 * pkin(4) + t117) - g(2) * (-t101 * t104 - t102 * t110 + t116) - g(3) * (-t101 * t110 + t102 * t104 + t118);];
U_reg = t1;
