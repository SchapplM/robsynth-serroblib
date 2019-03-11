% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRPR6_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:51:32
% EndTime: 2019-03-09 01:51:32
% DurationCPUTime: 0.07s
% Computational Cost: add. (107->38), mult. (195->79), div. (0->0), fcn. (69->4), ass. (0->30)
t114 = qJD(1) ^ 2;
t125 = t114 / 0.2e1;
t124 = -pkin(1) - qJ(3);
t111 = sin(qJ(4));
t123 = qJD(1) * t111;
t122 = qJD(4) * qJ(5);
t107 = qJD(1) * qJ(2) + qJD(3);
t101 = -qJD(1) * pkin(7) + t107;
t121 = qJD(4) * t101;
t102 = -t124 * qJD(1) - qJD(2);
t120 = t102 * qJD(1);
t113 = cos(qJ(4));
t119 = t113 * qJD(1);
t118 = pkin(4) * t123 - qJD(2);
t117 = qJD(1) * qJD(4);
t116 = pkin(5) * qJD(1) - t101;
t115 = -qJ(5) * t113 - t124;
t112 = cos(qJ(6));
t110 = sin(qJ(6));
t108 = -qJD(1) * pkin(1) + qJD(2);
t105 = qJD(6) + t119;
t100 = t112 * qJD(4) + t110 * t123;
t99 = t110 * qJD(4) - t112 * t123;
t98 = -t111 * t101 - t122;
t97 = -qJD(4) * pkin(4) - t113 * t101 + qJD(5);
t96 = t115 * qJD(1) + t118;
t95 = -t116 * t111 + t122;
t94 = qJD(5) + t116 * t113 + (-pkin(4) - pkin(8)) * qJD(4);
t93 = (pkin(8) * t111 + t115) * qJD(1) + t118;
t1 = [t125, 0, 0, t108 * qJD(1), t114 * qJ(2), qJ(2) ^ 2 * t125 + t108 ^ 2 / 0.2e1, t107 * qJD(1), t120, t102 ^ 2 / 0.2e1 + t107 ^ 2 / 0.2e1, t113 ^ 2 * t125, -t113 * t114 * t111, t113 * t117, -t111 * t117, qJD(4) ^ 2 / 0.2e1, t111 * t120 + t113 * t121, t102 * t119 - t111 * t121 (t111 * t98 + t113 * t97) * qJD(1), t97 * qJD(4) - t96 * t123, -t98 * qJD(4) - t96 * t119, t96 ^ 2 / 0.2e1 + t98 ^ 2 / 0.2e1 + t97 ^ 2 / 0.2e1, t100 ^ 2 / 0.2e1, -t100 * t99, t100 * t105, -t99 * t105, t105 ^ 2 / 0.2e1 (-t110 * t93 + t112 * t94) * t105 + t95 * t99 -(t110 * t94 + t112 * t93) * t105 + t95 * t100;];
T_reg  = t1;
