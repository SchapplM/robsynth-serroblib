% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% T_reg [1x30]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPRRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:29:00
% EndTime: 2019-03-09 02:29:00
% DurationCPUTime: 0.08s
% Computational Cost: add. (148->36), mult. (288->82), div. (0->0), fcn. (159->6), ass. (0->31)
t117 = (pkin(1) + qJ(3)) * qJD(1) - qJD(2);
t131 = qJD(1) ^ 2;
t139 = t131 / 0.2e1;
t130 = cos(qJ(4));
t121 = qJD(1) * qJ(2) + qJD(3);
t116 = -qJD(1) * pkin(7) + t121;
t133 = -pkin(8) * qJD(1) + t116;
t110 = qJD(4) * pkin(4) + t133 * t130;
t127 = sin(qJ(4));
t111 = t133 * t127;
t126 = sin(qJ(5));
t129 = cos(qJ(5));
t137 = t126 * t110 + t129 * t111;
t136 = qJD(4) * t116;
t135 = t117 * qJD(1);
t134 = qJD(1) * qJD(4);
t132 = t129 * t110 - t126 * t111;
t115 = t127 * qJD(1) * pkin(4) + t117;
t113 = (t126 * t130 + t127 * t129) * qJD(1);
t128 = cos(qJ(6));
t125 = sin(qJ(6));
t123 = qJD(4) + qJD(5);
t122 = -qJD(1) * pkin(1) + qJD(2);
t114 = (-t126 * t127 + t129 * t130) * qJD(1);
t112 = qJD(6) + t113;
t107 = t128 * t114 + t125 * t123;
t106 = t125 * t114 - t128 * t123;
t105 = t113 * pkin(5) - t114 * pkin(9) + t115;
t104 = t123 * pkin(9) + t137;
t103 = -t123 * pkin(5) - t132;
t1 = [t139, 0, 0, t122 * qJD(1), t131 * qJ(2), qJ(2) ^ 2 * t139 + t122 ^ 2 / 0.2e1, t121 * qJD(1), t135, t117 ^ 2 / 0.2e1 + t121 ^ 2 / 0.2e1, t130 ^ 2 * t139, -t130 * t131 * t127, t130 * t134, -t127 * t134, qJD(4) ^ 2 / 0.2e1, t127 * t135 + t130 * t136, -t127 * t136 + t130 * t135, t114 ^ 2 / 0.2e1, -t114 * t113, t114 * t123, -t113 * t123, t123 ^ 2 / 0.2e1, t115 * t113 + t132 * t123, t115 * t114 - t137 * t123, t107 ^ 2 / 0.2e1, -t107 * t106, t107 * t112, -t106 * t112, t112 ^ 2 / 0.2e1 (-t125 * t104 + t128 * t105) * t112 + t103 * t106 -(t128 * t104 + t125 * t105) * t112 + t103 * t107;];
T_reg  = t1;
