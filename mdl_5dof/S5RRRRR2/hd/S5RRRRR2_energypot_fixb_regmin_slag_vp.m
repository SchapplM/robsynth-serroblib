% Calculate minimal parameter regressor of potential energy for
% S5RRRRR2
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
% U_reg [1x27]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_energypot_fixb_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:53:42
% EndTime: 2019-12-05 18:53:42
% DurationCPUTime: 0.04s
% Computational Cost: add. (52->19), mult. (48->30), div. (0->0), fcn. (52->10), ass. (0->19)
t115 = qJ(3) + qJ(4);
t111 = sin(t115);
t128 = g(3) * t111;
t116 = qJ(1) + qJ(2);
t112 = sin(t116);
t117 = sin(qJ(5));
t127 = t112 * t117;
t120 = cos(qJ(5));
t126 = t112 * t120;
t114 = cos(t116);
t125 = t114 * t117;
t124 = t114 * t120;
t123 = g(1) * t114 + g(2) * t112;
t122 = cos(qJ(1));
t121 = cos(qJ(3));
t119 = sin(qJ(1));
t118 = sin(qJ(3));
t113 = cos(t115);
t1 = [0, -g(1) * t122 - g(2) * t119, g(1) * t119 - g(2) * t122, 0, -t123, g(1) * t112 - g(2) * t114, 0, 0, 0, 0, 0, -g(3) * t118 - t123 * t121, -g(3) * t121 + t123 * t118, 0, 0, 0, 0, 0, -t123 * t113 - t128, -g(3) * t113 + t123 * t111, 0, 0, 0, 0, 0, -g(1) * (t113 * t124 + t127) - g(2) * (t113 * t126 - t125) - t120 * t128, -g(1) * (-t113 * t125 + t126) - g(2) * (-t113 * t127 - t124) + t117 * t128;];
U_reg = t1;
