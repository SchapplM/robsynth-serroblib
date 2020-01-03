% Calculate minimal parameter regressor of potential energy for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x24]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:26:14
% EndTime: 2019-12-31 17:26:14
% DurationCPUTime: 0.04s
% Computational Cost: add. (30->16), mult. (44->26), div. (0->0), fcn. (48->8), ass. (0->16)
t117 = qJ(2) + qJ(3);
t115 = sin(t117);
t129 = g(3) * t115;
t118 = sin(qJ(4));
t120 = sin(qJ(1));
t128 = t120 * t118;
t121 = cos(qJ(4));
t127 = t120 * t121;
t123 = cos(qJ(1));
t126 = t123 * t118;
t125 = t123 * t121;
t124 = g(1) * t123 + g(2) * t120;
t122 = cos(qJ(2));
t119 = sin(qJ(2));
t116 = cos(t117);
t1 = [0, -t124, g(1) * t120 - g(2) * t123, 0, 0, 0, 0, 0, -g(3) * t119 - t124 * t122, -g(3) * t122 + t124 * t119, 0, 0, 0, 0, 0, -t124 * t116 - t129, -g(3) * t116 + t124 * t115, 0, 0, 0, 0, 0, -g(1) * (t116 * t125 + t128) - g(2) * (t116 * t127 - t126) - t121 * t129, -g(1) * (-t116 * t126 + t127) - g(2) * (-t116 * t128 - t125) + t118 * t129;];
U_reg = t1;
