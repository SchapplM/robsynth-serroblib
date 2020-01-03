% Calculate minimal parameter regressor of potential energy for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP12_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:19
% EndTime: 2019-12-31 18:57:20
% DurationCPUTime: 0.06s
% Computational Cost: add. (43->33), mult. (76->43), div. (0->0), fcn. (74->6), ass. (0->21)
t120 = cos(qJ(3));
t129 = g(3) * t120;
t116 = sin(qJ(4));
t118 = sin(qJ(1));
t128 = t118 * t116;
t119 = cos(qJ(4));
t127 = t118 * t119;
t121 = cos(qJ(1));
t126 = t121 * t116;
t125 = t121 * t119;
t124 = pkin(4) * t116 + pkin(6);
t123 = g(1) * (t121 * pkin(1) + t118 * qJ(2));
t109 = g(1) * t118 - g(2) * t121;
t111 = t119 * pkin(4) + pkin(3);
t115 = -qJ(5) - pkin(7);
t117 = sin(qJ(3));
t122 = t111 * t117 + t115 * t120;
t113 = t118 * pkin(1);
t110 = g(1) * t121 + g(2) * t118;
t108 = -g(3) * t117 + t109 * t120;
t1 = [0, -t110, t109, t110, -t109, -t123 - g(2) * (-t121 * qJ(2) + t113) - g(3) * pkin(5), 0, 0, 0, 0, 0, -t109 * t117 - t129, -t108, 0, 0, 0, 0, 0, -g(1) * (t117 * t127 + t126) - g(2) * (-t117 * t125 + t128) - t119 * t129, -g(1) * (-t117 * t128 + t125) - g(2) * (t117 * t126 + t127) + t116 * t129, t108, -t123 - g(2) * t113 - g(3) * (t120 * t111 - t117 * t115 + pkin(2) + pkin(5)) + (-g(1) * t122 - g(2) * t124) * t118 + (-g(1) * t124 - g(2) * (-qJ(2) - t122)) * t121;];
U_reg = t1;
