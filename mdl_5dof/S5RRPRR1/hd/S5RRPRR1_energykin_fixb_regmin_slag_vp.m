% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% T_reg [1x26]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_energykin_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:22:25
% EndTime: 2019-07-18 17:22:25
% DurationCPUTime: 0.11s
% Computational Cost: add. (101->30), mult. (255->71), div. (0->0), fcn. (157->6), ass. (0->30)
t118 = qJD(1) ^ 2;
t125 = t118 / 0.2e1;
t124 = pkin(3) + qJ(3);
t117 = cos(qJ(2));
t120 = qJD(1) * t117;
t104 = t124 * t120;
t113 = sin(qJ(4));
t116 = cos(qJ(4));
t111 = qJD(2) * pkin(1);
t114 = sin(qJ(2));
t121 = qJD(1) * t114;
t98 = qJD(2) * pkin(2) - t124 * t121 + t111;
t123 = t116 * t104 + t113 * t98;
t122 = t117 ^ 2 * t118;
t119 = qJD(1) * qJD(2);
t101 = t113 * t121 - t116 * t120;
t103 = qJD(3) + (-pkin(1) - pkin(2)) * t120;
t115 = cos(qJ(5));
t112 = sin(qJ(5));
t109 = qJD(2) + qJD(4);
t106 = -pkin(1) * t120 + qJD(3);
t105 = -qJ(3) * t121 + t111;
t102 = (t113 * t117 + t114 * t116) * qJD(1);
t99 = qJD(5) + t101;
t96 = t102 * t115 + t109 * t112;
t95 = t102 * t112 - t115 * t109;
t94 = -pkin(4) * t102 + t103;
t93 = t104 * t113 - t116 * t98;
t92 = pkin(4) * t109 + t123;
t1 = [t125, 0, 0, t114 ^ 2 * t125, t114 * t118 * t117, t114 * t119, t117 * t119, qJD(2) ^ 2 / 0.2e1, 0, 0, qJ(3) * t122 - t105 * t121, qJ(3) ^ 2 * t122 / 0.2e1 + t105 ^ 2 / 0.2e1 + t106 ^ 2 / 0.2e1, t102 ^ 2 / 0.2e1, -t102 * t101, t102 * t109, -t101 * t109, t109 ^ 2 / 0.2e1, t101 * t103 - t109 * t93, t103 * t102 - t123 * t109, t96 ^ 2 / 0.2e1, -t96 * t95, t96 * t99, -t95 * t99, t99 ^ 2 / 0.2e1, (-t112 * t92 + t115 * t94) * t99 + t93 * t95, -(t112 * t94 + t115 * t92) * t99 + t93 * t96;];
T_reg  = t1;
