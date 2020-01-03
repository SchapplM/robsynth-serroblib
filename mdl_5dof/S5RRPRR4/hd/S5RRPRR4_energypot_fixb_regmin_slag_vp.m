% Calculate minimal parameter regressor of potential energy for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:02:05
% EndTime: 2020-01-03 12:02:05
% DurationCPUTime: 0.03s
% Computational Cost: add. (46->18), mult. (35->25), div. (0->0), fcn. (32->10), ass. (0->13)
t123 = qJ(1) + qJ(2);
t117 = pkin(9) + t123;
t128 = g(2) * sin(t117) - g(3) * cos(t117);
t127 = cos(qJ(1));
t126 = cos(qJ(4));
t125 = sin(qJ(1));
t124 = sin(qJ(4));
t122 = qJ(4) + qJ(5);
t121 = cos(t123);
t120 = cos(t122);
t119 = sin(t123);
t118 = sin(t122);
t1 = [0, -g(2) * t125 + g(3) * t127, -g(2) * t127 - g(3) * t125, 0, -g(2) * t119 + g(3) * t121, -g(2) * t121 - g(3) * t119, -g(1) * (qJ(3) + pkin(6) + pkin(5)) - g(2) * (t125 * pkin(1) + pkin(2) * t119) - g(3) * (-t127 * pkin(1) - pkin(2) * t121), 0, 0, 0, 0, 0, -g(1) * t124 - t128 * t126, -g(1) * t126 + t128 * t124, 0, 0, 0, 0, 0, -g(1) * t118 - t128 * t120, -g(1) * t120 + t128 * t118;];
U_reg = t1;
