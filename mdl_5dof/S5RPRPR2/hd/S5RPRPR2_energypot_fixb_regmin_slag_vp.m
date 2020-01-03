% Calculate minimal parameter regressor of potential energy for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:34:00
% EndTime: 2020-01-03 11:34:00
% DurationCPUTime: 0.04s
% Computational Cost: add. (70->24), mult. (46->29), div. (0->0), fcn. (40->10), ass. (0->16)
t125 = qJ(2) + pkin(5);
t118 = qJ(1) + pkin(8);
t116 = qJ(3) + t118;
t112 = sin(t116);
t113 = cos(t116);
t124 = g(2) * t112 - g(3) * t113;
t121 = sin(qJ(1));
t122 = cos(qJ(1));
t123 = -g(2) * t121 + g(3) * t122;
t120 = cos(pkin(9));
t119 = sin(pkin(9));
t117 = pkin(9) + qJ(5);
t115 = cos(t117);
t114 = sin(t117);
t111 = g(2) * t113 + g(3) * t112;
t1 = [0, t123, -g(2) * t122 - g(3) * t121, t123 * pkin(1) - g(1) * t125, 0, -t124, -t111, -g(1) * t119 - t124 * t120, -g(1) * t120 + t124 * t119, t111, -g(1) * (pkin(6) + t125) - g(2) * (t112 * pkin(3) - t113 * qJ(4) + pkin(2) * sin(t118) + t121 * pkin(1)) - g(3) * (-t113 * pkin(3) - t112 * qJ(4) - pkin(2) * cos(t118) - t122 * pkin(1)), 0, 0, 0, 0, 0, -g(1) * t114 - t124 * t115, -g(1) * t115 + t124 * t114;];
U_reg = t1;
