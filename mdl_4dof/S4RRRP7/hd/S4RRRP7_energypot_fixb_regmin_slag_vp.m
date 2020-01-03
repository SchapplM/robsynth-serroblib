% Calculate minimal parameter regressor of potential energy for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRP7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:03
% EndTime: 2019-12-31 17:21:03
% DurationCPUTime: 0.06s
% Computational Cost: add. (44->30), mult. (96->42), div. (0->0), fcn. (105->6), ass. (0->19)
t118 = sin(qJ(2));
t128 = g(3) * t118;
t119 = sin(qJ(1));
t121 = cos(qJ(2));
t127 = t119 * t121;
t117 = sin(qJ(3));
t122 = cos(qJ(1));
t126 = t122 * t117;
t120 = cos(qJ(3));
t125 = t122 * t120;
t124 = g(1) * t122 + g(2) * t119;
t112 = t117 * t127 + t125;
t114 = -t119 * t120 + t121 * t126;
t123 = g(1) * t114 + g(2) * t112 + t117 * t128;
t115 = t119 * t117 + t121 * t125;
t113 = t120 * t127 - t126;
t111 = -g(3) * t121 + t124 * t118;
t110 = -g(1) * t115 - g(2) * t113 - t120 * t128;
t1 = [0, -t124, g(1) * t119 - g(2) * t122, 0, 0, 0, 0, 0, -t124 * t121 - t128, t111, 0, 0, 0, 0, 0, t110, t123, t110, -t111, -t123, -g(1) * (t115 * pkin(3) + t119 * pkin(5) + t114 * qJ(4) + (pkin(2) * t121 + pkin(1)) * t122) - g(2) * (t119 * pkin(1) + pkin(2) * t127 + t113 * pkin(3) - t122 * pkin(5) + t112 * qJ(4)) - g(3) * (-t121 * pkin(6) + pkin(4)) + (-g(3) * (pkin(3) * t120 + qJ(4) * t117 + pkin(2)) - t124 * pkin(6)) * t118;];
U_reg = t1;
