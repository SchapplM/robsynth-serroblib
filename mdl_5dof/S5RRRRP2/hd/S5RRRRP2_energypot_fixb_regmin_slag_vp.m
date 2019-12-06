% Calculate minimal parameter regressor of potential energy for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:48:12
% EndTime: 2019-12-05 18:48:12
% DurationCPUTime: 0.04s
% Computational Cost: add. (57->23), mult. (45->29), div. (0->0), fcn. (42->8), ass. (0->15)
t124 = qJ(1) + qJ(2);
t119 = sin(t124);
t121 = cos(t124);
t129 = g(2) * t119 - g(3) * t121;
t128 = cos(qJ(1));
t127 = cos(qJ(3));
t126 = sin(qJ(1));
t125 = sin(qJ(3));
t123 = qJ(3) + qJ(4);
t122 = -qJ(5) - pkin(8) - pkin(7);
t120 = cos(t123);
t118 = sin(t123);
t117 = t127 * pkin(3) + pkin(4) * t120 + pkin(2);
t116 = g(2) * t121 + g(3) * t119;
t1 = [0, g(2) * t126 - g(3) * t128, g(2) * t128 + g(3) * t126, 0, t129, t116, 0, 0, 0, 0, 0, -g(1) * t125 + t129 * t127, -g(1) * t127 - t129 * t125, 0, 0, 0, 0, 0, -g(1) * t118 + t129 * t120, -g(1) * t120 - t129 * t118, -t116, -g(1) * (t125 * pkin(3) + pkin(4) * t118 + pkin(5) + pkin(6)) - g(2) * (-t126 * pkin(1) - t119 * t117 - t121 * t122) - g(3) * (t128 * pkin(1) + t121 * t117 - t119 * t122);];
U_reg = t1;
