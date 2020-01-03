% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRRPR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:38
% EndTime: 2019-12-31 21:11:39
% DurationCPUTime: 0.11s
% Computational Cost: add. (138->32), mult. (201->73), div. (0->0), fcn. (92->6), ass. (0->29)
t129 = pkin(3) + pkin(4);
t111 = qJD(1) + qJD(2);
t108 = t111 ^ 2;
t128 = t108 / 0.2e1;
t127 = pkin(1) * qJD(1);
t114 = sin(qJ(3));
t126 = t111 * t114;
t117 = cos(qJ(3));
t125 = t111 * t117;
t121 = sin(qJ(2)) * t127;
t106 = t111 * pkin(7) + t121;
t103 = qJD(3) * qJ(4) + t117 * t106;
t124 = qJD(3) * t114;
t123 = qJD(3) * t117;
t122 = t114 * t106 + qJD(4);
t120 = cos(qJ(2)) * t127;
t119 = qJ(4) * t114 + pkin(2);
t116 = cos(qJ(5));
t113 = sin(qJ(5));
t109 = qJD(3) - qJD(5);
t107 = -t111 * pkin(2) - t120;
t102 = (-t113 * t117 + t114 * t116) * t111;
t101 = (t113 * t114 + t116 * t117) * t111;
t100 = -qJD(3) * pkin(3) + t122;
t99 = -pkin(8) * t125 + t103;
t98 = -t120 + (-pkin(3) * t117 - t119) * t111;
t97 = -pkin(8) * t126 - t129 * qJD(3) + t122;
t96 = t120 + (t129 * t117 + t119) * t111;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t128, t111 * t120, -t111 * t121, t114 ^ 2 * t128, t114 * t108 * t117, t111 * t124, t111 * t123, qJD(3) ^ 2 / 0.2e1, -t106 * t124 - t107 * t125, -t106 * t123 + t107 * t126, -qJD(3) * t100 - t98 * t125, (t100 * t114 + t103 * t117) * t111, qJD(3) * t103 - t98 * t126, t103 ^ 2 / 0.2e1 + t98 ^ 2 / 0.2e1 + t100 ^ 2 / 0.2e1, t102 ^ 2 / 0.2e1, -t102 * t101, -t102 * t109, t101 * t109, t109 ^ 2 / 0.2e1, t96 * t101 - (-t113 * t99 + t116 * t97) * t109, t96 * t102 + (t113 * t97 + t116 * t99) * t109;];
T_reg = t1;
