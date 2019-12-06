% Calculate minimal parameter regressor of potential energy for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:03
% EndTime: 2019-12-05 15:43:03
% DurationCPUTime: 0.04s
% Computational Cost: add. (63->21), mult. (46->26), div. (0->0), fcn. (42->10), ass. (0->14)
t124 = pkin(9) + qJ(4);
t125 = pkin(8) + qJ(2);
t120 = sin(t125);
t122 = cos(t125);
t128 = g(1) * t122 + g(2) * t120;
t127 = cos(pkin(9));
t126 = sin(pkin(9));
t123 = qJ(5) + t124;
t121 = cos(t124);
t119 = sin(t124);
t118 = cos(t123);
t117 = sin(t123);
t116 = g(1) * t120 - g(2) * t122;
t1 = [-g(3) * qJ(1), 0, -t128, t116, -g(3) * t126 - t128 * t127, -g(3) * t127 + t128 * t126, -t116, -g(1) * (t122 * pkin(2) + t120 * qJ(3) + cos(pkin(8)) * pkin(1)) - g(2) * (t120 * pkin(2) - t122 * qJ(3) + sin(pkin(8)) * pkin(1)) - g(3) * (pkin(5) + qJ(1)), 0, 0, 0, 0, 0, -g(3) * t119 - t128 * t121, -g(3) * t121 + t128 * t119, 0, 0, 0, 0, 0, -g(3) * t117 - t128 * t118, -g(3) * t118 + t128 * t117;];
U_reg = t1;
