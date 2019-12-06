% Calculate minimal parameter regressor of potential energy for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:29
% EndTime: 2019-12-05 18:58:29
% DurationCPUTime: 0.03s
% Computational Cost: add. (48->13), mult. (32->20), div. (0->0), fcn. (32->10), ass. (0->15)
t121 = qJ(1) + qJ(2);
t119 = qJ(3) + t121;
t113 = sin(t119);
t114 = cos(t119);
t126 = g(2) * t113 - g(3) * t114;
t125 = cos(qJ(1));
t124 = cos(qJ(4));
t123 = sin(qJ(1));
t122 = sin(qJ(4));
t120 = qJ(4) + qJ(5);
t118 = cos(t121);
t117 = cos(t120);
t116 = sin(t121);
t115 = sin(t120);
t1 = [0, g(2) * t123 - g(3) * t125, g(2) * t125 + g(3) * t123, 0, g(2) * t116 - g(3) * t118, g(2) * t118 + g(3) * t116, 0, t126, g(2) * t114 + g(3) * t113, 0, 0, 0, 0, 0, -g(1) * t122 + t126 * t124, -g(1) * t124 - t126 * t122, 0, 0, 0, 0, 0, -g(1) * t115 + t126 * t117, -g(1) * t117 - t126 * t115;];
U_reg = t1;
