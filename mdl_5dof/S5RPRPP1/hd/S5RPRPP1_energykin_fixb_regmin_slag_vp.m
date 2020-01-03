% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% T_reg [1x17]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRPP1_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:13
% EndTime: 2019-12-31 18:09:13
% DurationCPUTime: 0.06s
% Computational Cost: add. (137->31), mult. (321->73), div. (0->0), fcn. (173->6), ass. (0->27)
t116 = qJD(1) ^ 2;
t123 = t116 / 0.2e1;
t115 = cos(qJ(3));
t120 = qJD(1) * t115;
t111 = sin(pkin(7));
t105 = (pkin(1) * t111 + pkin(6)) * qJD(1);
t114 = sin(qJ(3));
t122 = t114 * qJD(2) + t115 * t105;
t100 = qJ(4) * t120 + t122;
t110 = sin(pkin(8));
t112 = cos(pkin(8));
t109 = t115 * qJD(2);
t99 = qJD(3) * pkin(3) + t109 + (-qJ(4) * qJD(1) - t105) * t114;
t95 = t112 * t100 + t110 * t99;
t121 = qJD(1) * t114;
t119 = qJD(1) * qJD(3);
t113 = cos(pkin(7));
t118 = -pkin(1) * t113 - pkin(2);
t94 = -t110 * t100 + t112 * t99;
t101 = qJD(4) + (-pkin(3) * t115 + t118) * qJD(1);
t106 = t118 * qJD(1);
t103 = (t110 * t115 + t112 * t114) * qJD(1);
t102 = t110 * t121 - t112 * t120;
t96 = t102 * pkin(4) - t103 * qJ(5) + t101;
t93 = qJD(3) * qJ(5) + t95;
t92 = -qJD(3) * pkin(4) + qJD(5) - t94;
t1 = [t123, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t111 ^ 2 / 0.2e1 + t113 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t116, t114 ^ 2 * t123, t114 * t116 * t115, t114 * t119, t115 * t119, qJD(3) ^ 2 / 0.2e1, -t106 * t120 + (-t114 * t105 + t109) * qJD(3), -t122 * qJD(3) + t106 * t121, -t95 * t102 - t94 * t103, t95 ^ 2 / 0.2e1 + t94 ^ 2 / 0.2e1 + t101 ^ 2 / 0.2e1, -t92 * qJD(3) + t96 * t102, -t93 * t102 + t92 * t103, t93 * qJD(3) - t96 * t103, t93 ^ 2 / 0.2e1 + t96 ^ 2 / 0.2e1 + t92 ^ 2 / 0.2e1;];
T_reg = t1;
