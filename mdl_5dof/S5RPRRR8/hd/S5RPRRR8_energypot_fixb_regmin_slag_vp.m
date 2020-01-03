% Calculate minimal parameter regressor of potential energy for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:06:05
% EndTime: 2019-12-31 19:06:05
% DurationCPUTime: 0.04s
% Computational Cost: add. (36->16), mult. (63->27), div. (0->0), fcn. (72->8), ass. (0->15)
t129 = cos(qJ(3));
t128 = sin(qJ(3));
t124 = sin(qJ(1));
t126 = cos(qJ(1));
t114 = -t124 * t128 - t126 * t129;
t115 = t124 * t129 - t126 * t128;
t127 = g(1) * t114 - g(2) * t115;
t125 = cos(qJ(4));
t123 = sin(qJ(4));
t122 = qJ(4) + qJ(5);
t121 = cos(t122);
t120 = sin(t122);
t117 = -g(1) * t126 - g(2) * t124;
t116 = g(1) * t124 - g(2) * t126;
t1 = [0, t117, t116, t117, -t116, -g(1) * (t126 * pkin(1) + t124 * qJ(2)) - g(2) * (t124 * pkin(1) - t126 * qJ(2)) - g(3) * pkin(5), 0, t127, -g(1) * t115 - g(2) * t114, 0, 0, 0, 0, 0, g(3) * t123 + t127 * t125, g(3) * t125 - t127 * t123, 0, 0, 0, 0, 0, g(3) * t120 + t127 * t121, g(3) * t121 - t127 * t120;];
U_reg = t1;
