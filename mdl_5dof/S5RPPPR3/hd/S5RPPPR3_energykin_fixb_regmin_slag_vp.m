% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% 
% Output:
% T_reg [1x19]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPPPR3_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:44:00
% EndTime: 2019-12-31 17:44:00
% DurationCPUTime: 0.06s
% Computational Cost: add. (82->33), mult. (203->65), div. (0->0), fcn. (102->6), ass. (0->25)
t101 = sin(pkin(8));
t103 = cos(pkin(8));
t102 = sin(pkin(7));
t99 = (pkin(1) * t102 + qJ(3)) * qJD(1);
t93 = t101 * qJD(2) + t103 * t99;
t112 = qJD(1) * t101;
t111 = qJD(1) * t103;
t104 = cos(pkin(7));
t110 = -pkin(1) * t104 - pkin(2);
t92 = t103 * qJD(2) - t101 * t99;
t91 = qJD(4) - t92;
t109 = qJ(4) * t101 - t110;
t107 = qJD(1) ^ 2;
t106 = cos(qJ(5));
t105 = sin(qJ(5));
t98 = t110 * qJD(1) + qJD(3);
t95 = (t101 * t106 - t103 * t105) * qJD(1);
t94 = (t101 * t105 + t103 * t106) * qJD(1);
t90 = t93 ^ 2 / 0.2e1;
t89 = qJD(3) + (-pkin(3) * t103 - t109) * qJD(1);
t88 = t93 * t111;
t87 = -pkin(6) * t111 + t93;
t86 = -pkin(6) * t112 + t91;
t85 = -qJD(3) + ((pkin(3) + pkin(4)) * t103 + t109) * qJD(1);
t1 = [t107 / 0.2e1, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t102 ^ 2 / 0.2e1 + t104 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t107, -t98 * t111, t98 * t112, -t92 * t112 + t88, t90 + t92 ^ 2 / 0.2e1 + t98 ^ 2 / 0.2e1, -t89 * t111, t91 * t112 + t88, -t89 * t112, t90 + t89 ^ 2 / 0.2e1 + t91 ^ 2 / 0.2e1, t95 ^ 2 / 0.2e1, -t95 * t94, t95 * qJD(5), -t94 * qJD(5), qJD(5) ^ 2 / 0.2e1, t85 * t94 + (-t105 * t87 + t106 * t86) * qJD(5), t85 * t95 - (t105 * t86 + t106 * t87) * qJD(5);];
T_reg = t1;
