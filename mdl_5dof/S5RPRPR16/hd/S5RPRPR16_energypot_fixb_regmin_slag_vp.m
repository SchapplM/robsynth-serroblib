% Calculate minimal parameter regressor of potential energy for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR16_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:32
% EndTime: 2019-12-31 18:39:32
% DurationCPUTime: 0.05s
% Computational Cost: add. (38->33), mult. (76->43), div. (0->0), fcn. (74->6), ass. (0->18)
t122 = sin(qJ(3));
t132 = pkin(3) * t122;
t131 = g(3) * t122;
t123 = sin(qJ(1));
t125 = cos(qJ(3));
t130 = t123 * t125;
t121 = sin(qJ(5));
t126 = cos(qJ(1));
t129 = t126 * t121;
t124 = cos(qJ(5));
t128 = t126 * t124;
t127 = t126 * pkin(1) + t123 * qJ(2);
t116 = g(1) * t123 - g(2) * t126;
t119 = t123 * pkin(1);
t117 = g(1) * t126 + g(2) * t123;
t115 = t116 * t125 - t131;
t114 = g(3) * t125 + t116 * t122;
t1 = [0, -t117, t116, t117, -t116, -g(1) * t127 - g(2) * (-t126 * qJ(2) + t119) - g(3) * pkin(5), 0, 0, 0, 0, 0, -t114, -t115, -t117, t114, t115, -g(1) * (-qJ(4) * t130 + t123 * t132 + t127) - g(2) * (t123 * pkin(6) + t119) - g(3) * (t125 * pkin(3) + t122 * qJ(4) + pkin(2) + pkin(5)) + (-g(1) * pkin(6) - g(2) * (qJ(4) * t125 - qJ(2) - t132)) * t126, 0, 0, 0, 0, 0, -g(1) * (-t121 * t130 + t128) - g(2) * (t123 * t124 + t125 * t129) - t121 * t131, -g(1) * (-t124 * t130 - t129) - g(2) * (-t123 * t121 + t125 * t128) - t124 * t131;];
U_reg = t1;
