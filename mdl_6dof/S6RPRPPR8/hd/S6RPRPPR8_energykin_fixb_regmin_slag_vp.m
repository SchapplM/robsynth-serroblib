% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% T_reg [1x28]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPRPPR8_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:59:59
% EndTime: 2019-03-09 02:59:59
% DurationCPUTime: 0.11s
% Computational Cost: add. (142->46), mult. (264->89), div. (0->0), fcn. (98->4), ass. (0->31)
t127 = -pkin(3) - pkin(4);
t117 = qJD(1) ^ 2;
t126 = t117 / 0.2e1;
t125 = t117 * qJ(2);
t108 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t114 = sin(qJ(3));
t104 = qJD(3) * qJ(4) + t114 * t108;
t124 = qJD(1) * t114;
t123 = qJD(3) * t108;
t116 = cos(qJ(3));
t122 = t116 * qJD(1);
t110 = qJ(4) * t122;
t121 = qJD(5) + t110;
t120 = -pkin(8) + t127;
t119 = qJD(1) * qJD(3);
t101 = -qJ(5) * t124 - t104;
t118 = qJD(4) + (-qJ(5) * qJD(1) - t108) * t116;
t115 = cos(qJ(6));
t113 = sin(qJ(6));
t111 = -qJD(1) * pkin(1) + qJD(2);
t109 = qJD(6) + t122;
t106 = -t113 * qJD(3) + t115 * t124;
t105 = t115 * qJD(3) + t113 * t124;
t103 = -t110 + (pkin(3) * t114 + qJ(2)) * qJD(1);
t102 = -qJD(3) * pkin(3) - t116 * t108 + qJD(4);
t100 = (t127 * t114 - qJ(2)) * qJD(1) + t121;
t99 = qJD(3) * pkin(5) - t101;
t98 = t127 * qJD(3) + t118;
t97 = t120 * qJD(3) + t118;
t96 = (pkin(5) * t116 + t120 * t114 - qJ(2)) * qJD(1) + t121;
t1 = [t126, 0, 0, t111 * qJD(1), t125, qJ(2) ^ 2 * t126 + t111 ^ 2 / 0.2e1, t116 ^ 2 * t126, -t116 * t117 * t114, t116 * t119, -t114 * t119, qJD(3) ^ 2 / 0.2e1, t114 * t125 + t116 * t123, -t114 * t123 + t116 * t125, -t102 * qJD(3) + t103 * t124 (t102 * t116 - t104 * t114) * qJD(1), t104 * qJD(3) - t103 * t122, t104 ^ 2 / 0.2e1 + t103 ^ 2 / 0.2e1 + t102 ^ 2 / 0.2e1, -t101 * qJD(3) + t100 * t122, t98 * qJD(3) + t100 * t124 (-t101 * t114 - t116 * t98) * qJD(1), t98 ^ 2 / 0.2e1 + t101 ^ 2 / 0.2e1 + t100 ^ 2 / 0.2e1, t106 ^ 2 / 0.2e1, -t106 * t105, t106 * t109, -t105 * t109, t109 ^ 2 / 0.2e1 (-t113 * t97 + t115 * t96) * t109 + t99 * t105 -(t113 * t96 + t115 * t97) * t109 + t99 * t106;];
T_reg  = t1;
