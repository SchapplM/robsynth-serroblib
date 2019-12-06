% Calculate minimal parameter regressor of potential energy for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:42
% EndTime: 2019-12-05 17:47:42
% DurationCPUTime: 0.04s
% Computational Cost: add. (41->24), mult. (51->28), div. (0->0), fcn. (45->6), ass. (0->14)
t129 = sin(qJ(3));
t134 = pkin(3) * t129;
t130 = sin(qJ(1));
t132 = cos(qJ(1));
t133 = t132 * pkin(1) + t130 * qJ(2);
t120 = g(1) * t130 - g(2) * t132;
t131 = cos(qJ(3));
t128 = -qJ(4) - pkin(6);
t126 = t130 * pkin(1);
t124 = qJ(3) + pkin(8) + qJ(5);
t123 = cos(t124);
t122 = sin(t124);
t121 = g(1) * t132 + g(2) * t130;
t1 = [0, -t121, t120, t121, -t120, -g(1) * t133 - g(2) * (-t132 * qJ(2) + t126) - g(3) * pkin(5), 0, 0, 0, 0, 0, -g(3) * t131 - t120 * t129, g(3) * t129 - t120 * t131, -t121, -g(1) * (-t132 * t128 + t130 * t134 + t133) - g(2) * (-t130 * t128 + t126 + (-qJ(2) - t134) * t132) - g(3) * (t131 * pkin(3) + pkin(2) + pkin(5)), 0, 0, 0, 0, 0, -g(3) * t123 - t120 * t122, g(3) * t122 - t120 * t123;];
U_reg = t1;
