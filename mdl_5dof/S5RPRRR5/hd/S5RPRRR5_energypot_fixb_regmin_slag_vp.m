% Calculate minimal parameter regressor of potential energy for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:27
% EndTime: 2019-12-05 18:16:27
% DurationCPUTime: 0.03s
% Computational Cost: add. (45->12), mult. (33->18), div. (0->0), fcn. (30->8), ass. (0->13)
t118 = qJ(1) + pkin(9) + qJ(3);
t116 = sin(t118);
t117 = cos(t118);
t127 = g(2) * t116 - g(3) * t117;
t123 = sin(qJ(1));
t125 = cos(qJ(1));
t126 = g(2) * t123 - g(3) * t125;
t124 = cos(qJ(4));
t122 = sin(qJ(4));
t121 = qJ(4) + qJ(5);
t120 = cos(t121);
t119 = sin(t121);
t1 = [0, t126, g(2) * t125 + g(3) * t123, -g(1) * (qJ(2) + pkin(5)) + t126 * pkin(1), 0, t127, g(2) * t117 + g(3) * t116, 0, 0, 0, 0, 0, -g(1) * t122 + t127 * t124, -g(1) * t124 - t127 * t122, 0, 0, 0, 0, 0, -g(1) * t119 + t127 * t120, -g(1) * t120 - t127 * t119;];
U_reg = t1;
