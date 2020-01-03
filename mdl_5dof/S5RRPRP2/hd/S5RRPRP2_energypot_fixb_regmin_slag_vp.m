% Calculate minimal parameter regressor of potential energy for
% S5RRPRP2
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
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:42
% EndTime: 2019-12-31 19:49:42
% DurationCPUTime: 0.05s
% Computational Cost: add. (79->25), mult. (58->34), div. (0->0), fcn. (52->8), ass. (0->17)
t117 = qJ(1) + qJ(2);
t125 = qJ(3) + pkin(6) + pkin(5);
t112 = pkin(8) + t117;
t108 = sin(t112);
t109 = cos(t112);
t124 = g(1) * t109 + g(2) * t108;
t118 = sin(qJ(4));
t120 = cos(qJ(4));
t123 = pkin(4) * t120 + qJ(5) * t118 + pkin(3);
t113 = sin(t117);
t114 = cos(t117);
t119 = sin(qJ(1));
t121 = cos(qJ(1));
t122 = -g(1) * (t121 * pkin(1) + pkin(2) * t114) - g(2) * (t119 * pkin(1) + pkin(2) * t113);
t107 = -g(3) * t118 - t124 * t120;
t106 = -g(3) * t120 + t124 * t118;
t1 = [0, -g(1) * t121 - g(2) * t119, g(1) * t119 - g(2) * t121, 0, -g(1) * t114 - g(2) * t113, g(1) * t113 - g(2) * t114, -g(3) * t125 + t122, 0, 0, 0, 0, 0, t107, t106, t107, -g(1) * t108 + g(2) * t109, -t106, -g(3) * (t118 * pkin(4) - t120 * qJ(5) + t125) + (g(2) * pkin(7) - g(1) * t123) * t109 + (-g(1) * pkin(7) - g(2) * t123) * t108 + t122;];
U_reg = t1;
