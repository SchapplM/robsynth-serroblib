% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRRP5
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
% T_reg [1x25]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRP5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:08:46
% EndTime: 2019-03-09 02:08:46
% DurationCPUTime: 0.10s
% Computational Cost: add. (114->33), mult. (213->72), div. (0->0), fcn. (91->4), ass. (0->27)
t113 = qJD(1) ^ 2;
t122 = t113 / 0.2e1;
t121 = cos(qJ(5));
t120 = -pkin(1) - qJ(3);
t110 = sin(qJ(5));
t111 = sin(qJ(4));
t112 = cos(qJ(4));
t96 = -qJD(2) + (pkin(4) * t111 - pkin(8) * t112 - t120) * qJD(1);
t106 = qJ(2) * qJD(1) + qJD(3);
t102 = -pkin(7) * qJD(1) + t106;
t98 = qJD(4) * pkin(8) + t102 * t111;
t119 = t110 * t96 + t121 * t98;
t118 = qJD(1) * t112;
t117 = qJD(4) * t102;
t103 = -t120 * qJD(1) - qJD(2);
t116 = t103 * qJD(1);
t115 = qJD(1) * qJD(4);
t114 = -t110 * t98 + t121 * t96;
t99 = -qJD(4) * pkin(4) - t102 * t112;
t107 = -qJD(1) * pkin(1) + qJD(2);
t105 = qJD(1) * t111 + qJD(5);
t101 = t110 * qJD(4) + t121 * t118;
t100 = -t121 * qJD(4) + t110 * t118;
t93 = pkin(5) * t100 + qJD(6) + t99;
t92 = -qJ(6) * t100 + t119;
t91 = pkin(5) * t105 - qJ(6) * t101 + t114;
t1 = [t122, 0, 0, t107 * qJD(1), t113 * qJ(2), qJ(2) ^ 2 * t122 + t107 ^ 2 / 0.2e1, t106 * qJD(1), t116, t103 ^ 2 / 0.2e1 + t106 ^ 2 / 0.2e1, t112 ^ 2 * t122, -t112 * t113 * t111, t112 * t115, -t111 * t115, qJD(4) ^ 2 / 0.2e1, t111 * t116 + t112 * t117, -t111 * t117 + t112 * t116, t101 ^ 2 / 0.2e1, -t101 * t100, t101 * t105, -t100 * t105, t105 ^ 2 / 0.2e1, t99 * t100 + t114 * t105, t99 * t101 - t119 * t105, -t100 * t92 - t101 * t91, t92 ^ 2 / 0.2e1 + t91 ^ 2 / 0.2e1 + t93 ^ 2 / 0.2e1;];
T_reg  = t1;
