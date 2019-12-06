% Calculate minimal parameter regressor of potential energy for
% S5PRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:40
% EndTime: 2019-12-05 16:19:40
% DurationCPUTime: 0.04s
% Computational Cost: add. (54->21), mult. (39->24), div. (0->0), fcn. (35->8), ass. (0->13)
t128 = pkin(8) + qJ(2);
t125 = sin(t128);
t126 = cos(t128);
t132 = g(1) * t126 + g(2) * t125;
t131 = cos(qJ(3));
t130 = sin(qJ(3));
t129 = -pkin(6) - qJ(4);
t127 = qJ(3) + pkin(9) + qJ(5);
t124 = t131 * pkin(3) + pkin(2);
t123 = cos(t127);
t122 = sin(t127);
t121 = g(1) * t125 - g(2) * t126;
t1 = [-g(3) * qJ(1), 0, -t132, t121, 0, 0, 0, 0, 0, -g(3) * t130 - t132 * t131, -g(3) * t131 + t132 * t130, -t121, -g(1) * (t126 * t124 - t125 * t129 + cos(pkin(8)) * pkin(1)) - g(2) * (t125 * t124 + t126 * t129 + sin(pkin(8)) * pkin(1)) - g(3) * (t130 * pkin(3) + pkin(5) + qJ(1)), 0, 0, 0, 0, 0, -g(3) * t122 - t132 * t123, -g(3) * t123 + t132 * t122;];
U_reg = t1;
