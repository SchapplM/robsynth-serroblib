% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5PRRPP2_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:10:08
% EndTime: 2019-12-05 16:10:08
% DurationCPUTime: 0.11s
% Computational Cost: add. (114->27), mult. (263->66), div. (0->0), fcn. (151->6), ass. (0->28)
t121 = qJD(2) ^ 2;
t129 = t121 / 0.2e1;
t117 = sin(qJ(3));
t118 = sin(qJ(2));
t112 = qJD(2) * pkin(6) + t118 * qJD(1);
t122 = qJ(4) * qJD(2) + t112;
t107 = qJD(3) * pkin(3) - t122 * t117;
t119 = cos(qJ(3));
t108 = t122 * t119;
t115 = sin(pkin(8));
t116 = cos(pkin(8));
t104 = t115 * t107 + t116 * t108;
t128 = qJD(2) * t117;
t127 = qJD(2) * t119;
t126 = qJD(3) * t112;
t120 = cos(qJ(2));
t125 = t120 * qJD(1);
t124 = qJD(1) * qJD(2);
t123 = qJD(2) * qJD(3);
t103 = t116 * t107 - t115 * t108;
t111 = -t125 + qJD(4) + (-pkin(3) * t119 - pkin(2)) * qJD(2);
t113 = -qJD(2) * pkin(2) - t125;
t110 = (t115 * t119 + t116 * t117) * qJD(2);
t109 = t115 * t128 - t116 * t127;
t102 = qJD(3) * qJ(5) + t104;
t101 = t109 * pkin(4) - t110 * qJ(5) + t111;
t100 = -qJD(3) * pkin(4) + qJD(5) - t103;
t1 = [qJD(1) ^ 2 / 0.2e1, t129, t120 * t124, -t118 * t124, t117 ^ 2 * t129, t117 * t121 * t119, t117 * t123, t119 * t123, qJD(3) ^ 2 / 0.2e1, -t113 * t127 - t117 * t126, t113 * t128 - t119 * t126, -t103 * t110 - t104 * t109, t104 ^ 2 / 0.2e1 + t103 ^ 2 / 0.2e1 + t111 ^ 2 / 0.2e1, -t100 * qJD(3) + t101 * t109, t100 * t110 - t102 * t109, t102 * qJD(3) - t101 * t110, t102 ^ 2 / 0.2e1 + t101 ^ 2 / 0.2e1 + t100 ^ 2 / 0.2e1;];
T_reg = t1;
