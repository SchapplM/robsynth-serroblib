% Calculate minimal parameter regressor of fixed base kinetic energy for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S6RPPPRR4_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_energykin_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_energykin_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:51
% EndTime: 2019-03-09 01:35:51
% DurationCPUTime: 0.06s
% Computational Cost: add. (126->39), mult. (223->75), div. (0->0), fcn. (89->6), ass. (0->31)
t118 = qJD(1) ^ 2;
t126 = t118 / 0.2e1;
t115 = sin(qJ(5));
t117 = cos(qJ(5));
t103 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t111 = sin(pkin(9));
t112 = cos(pkin(9));
t123 = qJ(2) * qJD(1);
t99 = t112 * t103 - t111 * t123;
t96 = qJD(1) * pkin(3) + qJD(4) - t99;
t95 = qJD(1) * pkin(7) + t96;
t125 = t117 * qJD(3) + t115 * t95;
t124 = t111 * t103;
t122 = qJD(1) * t117;
t121 = t115 * qJD(1);
t120 = qJD(1) * qJD(5);
t119 = -t115 * qJD(3) + t117 * t95;
t100 = t112 * t123 + t124;
t116 = cos(qJ(6));
t114 = sin(qJ(6));
t110 = qJD(1) * qJ(4);
t109 = qJD(3) ^ 2 / 0.2e1;
t107 = -qJD(1) * pkin(1) + qJD(2);
t104 = -qJD(6) + t121;
t102 = t114 * qJD(5) - t116 * t122;
t101 = t116 * qJD(5) + t114 * t122;
t97 = t100 - t110;
t93 = t124 - t110 + (-pkin(5) * t115 + pkin(8) * t117 + qJ(2) * t112) * qJD(1);
t92 = qJD(5) * pkin(8) + t125;
t91 = -qJD(5) * pkin(5) - t119;
t1 = [t126, 0, 0, -t107 * qJD(1), t118 * qJ(2), qJ(2) ^ 2 * t126 + t107 ^ 2 / 0.2e1, -t99 * qJD(1), t100 * qJD(1), t100 ^ 2 / 0.2e1 + t99 ^ 2 / 0.2e1 + t109, -t96 * qJD(1), -t97 * qJD(1), t109 + t97 ^ 2 / 0.2e1 + t96 ^ 2 / 0.2e1, t117 ^ 2 * t126, -t117 * t118 * t115, -t117 * t120, t115 * t120, qJD(5) ^ 2 / 0.2e1, t119 * qJD(5) - t97 * t121, -t125 * qJD(5) - t97 * t122, t102 ^ 2 / 0.2e1, t102 * t101, -t102 * t104, -t101 * t104, t104 ^ 2 / 0.2e1 -(-t114 * t92 + t116 * t93) * t104 - t91 * t101 (t114 * t93 + t116 * t92) * t104 + t91 * t102;];
T_reg  = t1;
