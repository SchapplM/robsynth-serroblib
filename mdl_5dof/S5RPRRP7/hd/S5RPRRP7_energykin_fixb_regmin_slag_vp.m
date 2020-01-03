% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% T_reg [1x22]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRP7_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:45:36
% EndTime: 2019-12-31 18:45:36
% DurationCPUTime: 0.07s
% Computational Cost: add. (133->31), mult. (300->74), div. (0->0), fcn. (160->6), ass. (0->28)
t124 = qJD(1) ^ 2;
t134 = t124 / 0.2e1;
t118 = sin(pkin(8));
t112 = (pkin(1) * t118 + pkin(6)) * qJD(1);
t121 = sin(qJ(3));
t123 = cos(qJ(3));
t132 = t121 * qJD(2) + t123 * t112;
t106 = qJD(3) * pkin(7) + t132;
t119 = cos(pkin(8));
t128 = -pkin(1) * t119 - pkin(2);
t107 = (-pkin(3) * t123 - pkin(7) * t121 + t128) * qJD(1);
t120 = sin(qJ(4));
t122 = cos(qJ(4));
t133 = t122 * t106 + t120 * t107;
t131 = qJD(1) * t121;
t130 = t123 * qJD(1);
t129 = qJD(1) * qJD(3);
t127 = t123 * qJD(2) - t121 * t112;
t126 = -t120 * t106 + t122 * t107;
t105 = -qJD(3) * pkin(3) - t127;
t114 = -qJD(4) + t130;
t113 = t128 * qJD(1);
t111 = t120 * qJD(3) + t122 * t131;
t110 = -t122 * qJD(3) + t120 * t131;
t102 = t110 * pkin(4) - t111 * qJ(5) + t105;
t101 = -t114 * qJ(5) + t133;
t100 = t114 * pkin(4) + qJD(5) - t126;
t1 = [t134, 0, 0, qJD(2) ^ 2 / 0.2e1 + (t118 ^ 2 / 0.2e1 + t119 ^ 2 / 0.2e1) * pkin(1) ^ 2 * t124, t121 ^ 2 * t134, t121 * t124 * t123, t121 * t129, t123 * t129, qJD(3) ^ 2 / 0.2e1, t127 * qJD(3) - t113 * t130, -t132 * qJD(3) + t113 * t131, t111 ^ 2 / 0.2e1, -t111 * t110, -t111 * t114, t110 * t114, t114 ^ 2 / 0.2e1, t105 * t110 - t126 * t114, t105 * t111 + t133 * t114, t100 * t114 + t102 * t110, t100 * t111 - t101 * t110, -t101 * t114 - t102 * t111, t101 ^ 2 / 0.2e1 + t102 ^ 2 / 0.2e1 + t100 ^ 2 / 0.2e1;];
T_reg = t1;
