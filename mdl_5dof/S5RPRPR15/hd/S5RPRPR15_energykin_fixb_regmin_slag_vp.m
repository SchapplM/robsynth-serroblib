% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPR15_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:23
% EndTime: 2019-12-31 18:37:23
% DurationCPUTime: 0.06s
% Computational Cost: add. (144->35), mult. (300->76), div. (0->0), fcn. (163->6), ass. (0->29)
t123 = qJD(1) ^ 2;
t130 = t123 / 0.2e1;
t129 = cos(pkin(8));
t128 = t123 * qJ(2);
t120 = sin(qJ(3));
t122 = cos(qJ(3));
t110 = (pkin(3) * t120 - qJ(4) * t122 + qJ(2)) * qJD(1);
t114 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t111 = qJD(3) * qJ(4) + t120 * t114;
t118 = sin(pkin(8));
t101 = t118 * t110 + t129 * t111;
t127 = qJD(1) * t122;
t126 = qJD(3) * t114;
t125 = t120 * qJD(1);
t124 = qJD(1) * qJD(3);
t100 = t129 * t110 - t118 * t111;
t109 = -qJD(3) * pkin(3) - t122 * t114 + qJD(4);
t121 = cos(qJ(5));
t119 = sin(qJ(5));
t116 = -qJD(1) * pkin(1) + qJD(2);
t115 = qJD(5) + t125;
t113 = t118 * qJD(3) + t129 * t127;
t112 = -t129 * qJD(3) + t118 * t127;
t104 = t112 * pkin(4) + t109;
t103 = -t119 * t112 + t121 * t113;
t102 = t121 * t112 + t119 * t113;
t99 = -t112 * pkin(7) + t101;
t98 = pkin(4) * t125 - t113 * pkin(7) + t100;
t1 = [t130, 0, 0, t116 * qJD(1), t128, qJ(2) ^ 2 * t130 + t116 ^ 2 / 0.2e1, t122 ^ 2 * t130, -t122 * t123 * t120, t122 * t124, -t120 * t124, qJD(3) ^ 2 / 0.2e1, t120 * t128 + t122 * t126, -t120 * t126 + t122 * t128, t100 * t125 + t109 * t112, -t101 * t125 + t109 * t113, -t100 * t113 - t101 * t112, t101 ^ 2 / 0.2e1 + t100 ^ 2 / 0.2e1 + t109 ^ 2 / 0.2e1, t103 ^ 2 / 0.2e1, -t103 * t102, t103 * t115, -t102 * t115, t115 ^ 2 / 0.2e1, (-t119 * t99 + t121 * t98) * t115 + t104 * t102, -(t119 * t98 + t121 * t99) * t115 + t104 * t103;];
T_reg = t1;
