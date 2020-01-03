% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% T_reg [1x27]
%   minimal parameter regressor of kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR13_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:15:27
% EndTime: 2019-12-31 19:15:27
% DurationCPUTime: 0.06s
% Computational Cost: add. (124->34), mult. (263->75), div. (0->0), fcn. (154->6), ass. (0->29)
t121 = qJD(1) ^ 2;
t129 = t121 / 0.2e1;
t128 = cos(qJ(4));
t127 = t121 * qJ(2);
t118 = sin(qJ(3));
t120 = cos(qJ(3));
t106 = (pkin(3) * t118 - pkin(7) * t120 + qJ(2)) * qJD(1);
t112 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t107 = qJD(3) * pkin(7) + t118 * t112;
t117 = sin(qJ(4));
t126 = t117 * t106 + t128 * t107;
t125 = qJD(1) * t120;
t124 = qJD(3) * t112;
t123 = qJD(1) * qJD(3);
t122 = t128 * t106 - t117 * t107;
t113 = t118 * qJD(1) + qJD(4);
t108 = -qJD(3) * pkin(3) - t120 * t112;
t119 = cos(qJ(5));
t116 = sin(qJ(5));
t114 = -pkin(1) * qJD(1) + qJD(2);
t111 = qJD(5) + t113;
t110 = t117 * qJD(3) + t128 * t125;
t109 = -t128 * qJD(3) + t117 * t125;
t101 = t109 * pkin(4) + t108;
t100 = -t116 * t109 + t119 * t110;
t99 = t119 * t109 + t116 * t110;
t98 = -t109 * pkin(8) + t126;
t97 = t113 * pkin(4) - t110 * pkin(8) + t122;
t1 = [t129, 0, 0, t114 * qJD(1), t127, qJ(2) ^ 2 * t129 + t114 ^ 2 / 0.2e1, t120 ^ 2 * t129, -t120 * t121 * t118, t120 * t123, -t118 * t123, qJD(3) ^ 2 / 0.2e1, t118 * t127 + t120 * t124, -t118 * t124 + t120 * t127, t110 ^ 2 / 0.2e1, -t110 * t109, t110 * t113, -t109 * t113, t113 ^ 2 / 0.2e1, t108 * t109 + t122 * t113, t108 * t110 - t126 * t113, t100 ^ 2 / 0.2e1, -t100 * t99, t100 * t111, -t99 * t111, t111 ^ 2 / 0.2e1, (-t116 * t98 + t119 * t97) * t111 + t101 * t99, -(t116 * t97 + t119 * t98) * t111 + t101 * t100;];
T_reg = t1;
