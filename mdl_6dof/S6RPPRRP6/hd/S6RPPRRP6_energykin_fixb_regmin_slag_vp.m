% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRP6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:11:17
% EndTime: 2019-03-09 02:11:17
% DurationCPUTime: 0.10s
% Computational Cost: add. (154->35), mult. (269->76), div. (0->0), fcn. (118->4), ass. (0->27)
t121 = qJD(1) ^ 2;
t129 = t121 / 0.2e1;
t128 = -pkin(1) - qJ(3);
t118 = sin(qJ(4));
t120 = cos(qJ(4));
t103 = -qJD(2) + (pkin(4) * t118 - pkin(8) * t120 - t128) * qJD(1);
t113 = qJ(2) * qJD(1) + qJD(3);
t109 = -pkin(7) * qJD(1) + t113;
t105 = qJD(4) * pkin(8) + t109 * t118;
t117 = sin(qJ(5));
t119 = cos(qJ(5));
t127 = t117 * t103 + t119 * t105;
t126 = qJD(1) * t120;
t125 = qJD(4) * t109;
t110 = -t128 * qJD(1) - qJD(2);
t124 = t110 * qJD(1);
t123 = qJD(1) * qJD(4);
t106 = -qJD(4) * pkin(4) - t109 * t120;
t122 = t103 * t119 - t117 * t105;
t114 = -qJD(1) * pkin(1) + qJD(2);
t112 = qJD(1) * t118 + qJD(5);
t108 = qJD(4) * t117 + t119 * t126;
t107 = -t119 * qJD(4) + t117 * t126;
t101 = pkin(5) * t107 - qJ(6) * t108 + t106;
t100 = qJ(6) * t112 + t127;
t99 = -pkin(5) * t112 + qJD(6) - t122;
t1 = [t129, 0, 0, t114 * qJD(1), t121 * qJ(2), qJ(2) ^ 2 * t129 + t114 ^ 2 / 0.2e1, t113 * qJD(1), t124, t110 ^ 2 / 0.2e1 + t113 ^ 2 / 0.2e1, t120 ^ 2 * t129, -t120 * t121 * t118, t120 * t123, -t118 * t123, qJD(4) ^ 2 / 0.2e1, t118 * t124 + t120 * t125, -t118 * t125 + t120 * t124, t108 ^ 2 / 0.2e1, -t108 * t107, t108 * t112, -t107 * t112, t112 ^ 2 / 0.2e1, t106 * t107 + t122 * t112, t106 * t108 - t127 * t112, t101 * t107 - t112 * t99, -t100 * t107 + t108 * t99, t100 * t112 - t101 * t108, t100 ^ 2 / 0.2e1 + t101 ^ 2 / 0.2e1 + t99 ^ 2 / 0.2e1;];
T_reg  = t1;
