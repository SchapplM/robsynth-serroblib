% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPP7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_energykin_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:05:53
% EndTime: 2019-12-31 21:05:53
% DurationCPUTime: 0.08s
% Computational Cost: add. (183->36), mult. (394->76), div. (0->0), fcn. (219->4), ass. (0->30)
t136 = pkin(3) + pkin(4);
t123 = qJD(1) ^ 2;
t135 = t123 / 0.2e1;
t134 = cos(qJ(3));
t122 = cos(qJ(2));
t133 = t122 * t123;
t121 = sin(qJ(2));
t108 = (-pkin(2) * t122 - pkin(7) * t121 - pkin(1)) * qJD(1);
t130 = t122 * qJD(1);
t114 = pkin(6) * t130 + qJD(2) * pkin(7);
t120 = sin(qJ(3));
t132 = t120 * t108 + t134 * t114;
t131 = qJD(1) * t121;
t129 = qJD(1) * qJD(2);
t116 = -qJD(3) + t130;
t105 = -t116 * qJ(4) + t132;
t128 = t121 * t129;
t127 = t122 * t129;
t113 = -qJD(2) * pkin(2) + pkin(6) * t131;
t126 = t134 * t108 - t120 * t114;
t125 = qJD(4) - t126;
t110 = t120 * qJD(2) + t134 * t131;
t124 = t110 * qJ(4) - t113;
t109 = -t134 * qJD(2) + t120 * t131;
t106 = t109 * pkin(3) - t124;
t104 = t116 * pkin(3) + t125;
t103 = -t136 * t109 + qJD(5) + t124;
t102 = t109 * qJ(5) + t105;
t101 = -t110 * qJ(5) + t136 * t116 + t125;
t1 = [t135, 0, 0, t121 ^ 2 * t135, t121 * t133, t128, t127, qJD(2) ^ 2 / 0.2e1, pkin(1) * t133 - pkin(6) * t128, -t123 * pkin(1) * t121 - pkin(6) * t127, t110 ^ 2 / 0.2e1, -t110 * t109, -t110 * t116, t109 * t116, t116 ^ 2 / 0.2e1, t113 * t109 - t126 * t116, t113 * t110 + t132 * t116, t104 * t116 + t106 * t109, t104 * t110 - t105 * t109, -t105 * t116 - t106 * t110, t105 ^ 2 / 0.2e1 + t106 ^ 2 / 0.2e1 + t104 ^ 2 / 0.2e1, t101 * t116 - t103 * t109, -t102 * t116 + t103 * t110, -t101 * t110 + t102 * t109, t102 ^ 2 / 0.2e1 + t101 ^ 2 / 0.2e1 + t103 ^ 2 / 0.2e1;];
T_reg = t1;
