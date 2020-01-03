% Calculate minimal parameter regressor of potential energy for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPPR5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:29
% EndTime: 2019-12-31 17:46:30
% DurationCPUTime: 0.05s
% Computational Cost: add. (61->28), mult. (99->38), div. (0->0), fcn. (108->8), ass. (0->21)
t130 = sin(qJ(1));
t129 = g(3) * (-qJ(3) + pkin(5));
t128 = cos(pkin(7));
t127 = sin(pkin(7));
t120 = cos(qJ(1));
t126 = t120 * pkin(1) + t130 * qJ(2);
t125 = t120 * pkin(2) + t126;
t124 = t130 * pkin(1) - t120 * qJ(2);
t102 = -t120 * t128 - t130 * t127;
t103 = t120 * t127 - t130 * t128;
t123 = g(1) * t103 - g(2) * t102;
t122 = g(1) * t102 + g(2) * t103;
t121 = t130 * pkin(2) + t124;
t118 = cos(pkin(8));
t117 = sin(pkin(8));
t116 = pkin(8) + qJ(5);
t110 = cos(t116);
t109 = sin(t116);
t105 = -g(1) * t120 - g(2) * t130;
t104 = g(1) * t130 - g(2) * t120;
t1 = [0, t105, t104, t105, -t104, -g(3) * pkin(5) - g(1) * t126 - g(2) * t124, t122, t123, -g(1) * t125 - g(2) * t121 - t129, g(3) * t117 + t122 * t118, g(3) * t118 - t122 * t117, -t123, -g(1) * (-t102 * pkin(3) + t103 * qJ(4) + t125) - g(2) * (-t103 * pkin(3) - t102 * qJ(4) + t121) - t129, 0, 0, 0, 0, 0, g(3) * t109 + t122 * t110, g(3) * t110 - t122 * t109;];
U_reg = t1;
