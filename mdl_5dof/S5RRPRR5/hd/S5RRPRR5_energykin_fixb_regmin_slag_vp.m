% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% T_reg [1x24]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RRPRR5_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_energykin_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:57
% EndTime: 2020-01-03 12:03:57
% DurationCPUTime: 0.11s
% Computational Cost: add. (210->33), mult. (311->73), div. (0->0), fcn. (197->8), ass. (0->32)
t129 = cos(qJ(4));
t114 = cos(pkin(9));
t112 = qJD(1) + qJD(2);
t127 = pkin(1) * qJD(1);
t124 = sin(qJ(2)) * t127;
t107 = t112 * qJ(3) + t124;
t122 = pkin(7) * t112 + t107;
t100 = t122 * t114;
t116 = sin(qJ(4));
t113 = sin(pkin(9));
t99 = t122 * t113;
t128 = t129 * t100 - t116 * t99;
t126 = t112 * t113;
t125 = t112 * t114;
t123 = cos(qJ(2)) * t127;
t121 = -t116 * t100 - t129 * t99;
t120 = qJD(3) - t123;
t102 = (-pkin(3) * t114 - pkin(2)) * t112 + t120;
t118 = cos(qJ(5));
t115 = sin(qJ(5));
t111 = qJD(4) + qJD(5);
t110 = t114 ^ 2;
t109 = t113 ^ 2;
t105 = -t112 * pkin(2) + t120;
t104 = (t129 * t113 + t114 * t116) * t112;
t103 = t116 * t126 - t129 * t125;
t95 = t103 * pkin(4) + t102;
t94 = -t115 * t103 + t118 * t104;
t93 = t118 * t103 + t115 * t104;
t92 = -t103 * pkin(8) + t128;
t91 = qJD(4) * pkin(4) - t104 * pkin(8) + t121;
t1 = [qJD(1) ^ 2 / 0.2e1, 0, 0, t112 ^ 2 / 0.2e1, t112 * t123, -t112 * t124, -t105 * t125, t105 * t126, (t109 + t110) * t112 * t107, t105 ^ 2 / 0.2e1 + (t110 / 0.2e1 + t109 / 0.2e1) * t107 ^ 2, t104 ^ 2 / 0.2e1, -t104 * t103, t104 * qJD(4), -t103 * qJD(4), qJD(4) ^ 2 / 0.2e1, qJD(4) * t121 + t102 * t103, -t128 * qJD(4) + t102 * t104, t94 ^ 2 / 0.2e1, -t94 * t93, t94 * t111, -t93 * t111, t111 ^ 2 / 0.2e1, t95 * t93 + (-t115 * t92 + t118 * t91) * t111, t95 * t94 - (t115 * t91 + t118 * t92) * t111;];
T_reg = t1;
