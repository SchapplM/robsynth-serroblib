% Calculate minimal parameter regressor of potential energy for
% S5RPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:23
% EndTime: 2019-12-31 17:49:23
% DurationCPUTime: 0.06s
% Computational Cost: add. (90->32), mult. (74->41), div. (0->0), fcn. (65->8), ass. (0->22)
t125 = qJ(2) + pkin(5);
t132 = g(3) * t125;
t122 = qJ(1) + pkin(7);
t116 = sin(t122);
t118 = cos(t122);
t131 = g(1) * t118 + g(2) * t116;
t127 = sin(qJ(1));
t128 = cos(qJ(1));
t130 = -g(1) * t128 - g(2) * t127;
t121 = pkin(8) + qJ(4);
t115 = sin(t121);
t117 = cos(t121);
t124 = cos(pkin(8));
t129 = t124 * pkin(3) + pkin(4) * t117 + qJ(5) * t115 + pkin(2);
t126 = -pkin(6) - qJ(3);
t123 = sin(pkin(8));
t120 = t128 * pkin(1);
t119 = t127 * pkin(1);
t113 = -g(1) * t116 + g(2) * t118;
t112 = -g(3) * t115 - t131 * t117;
t111 = -g(3) * t117 + t131 * t115;
t1 = [0, t130, g(1) * t127 - g(2) * t128, t130 * pkin(1) - t132, -g(3) * t123 - t131 * t124, -g(3) * t124 + t131 * t123, t113, -g(1) * (t118 * pkin(2) + t116 * qJ(3) + t120) - g(2) * (t116 * pkin(2) - t118 * qJ(3) + t119) - t132, 0, 0, 0, 0, 0, t112, t111, t112, t113, -t111, -g(1) * t120 - g(2) * t119 - g(3) * (t123 * pkin(3) + t115 * pkin(4) - t117 * qJ(5) + t125) + (-g(1) * t129 - g(2) * t126) * t118 + (g(1) * t126 - g(2) * t129) * t116;];
U_reg = t1;
