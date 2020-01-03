% Calculate minimal parameter regressor of fixed base kinetic energy for
% S5RPRRR12
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
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_reg = S5RPRRR12_energykin_fixb_regmin_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_energykin_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_energykin_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_energykin_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_kinetic_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:13:08
% EndTime: 2019-12-31 19:13:08
% DurationCPUTime: 0.07s
% Computational Cost: add. (128->33), mult. (270->74), div. (0->0), fcn. (159->6), ass. (0->29)
t116 = qJD(1) ^ 2;
t123 = t116 / 0.2e1;
t111 = sin(qJ(4));
t114 = cos(qJ(4));
t115 = cos(qJ(3));
t104 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t118 = -pkin(7) * qJD(1) + t104;
t98 = qJD(3) * pkin(3) + t118 * t115;
t112 = sin(qJ(3));
t99 = t118 * t112;
t122 = t111 * t98 + t114 * t99;
t121 = t116 * qJ(2);
t103 = (pkin(3) * t112 + qJ(2)) * qJD(1);
t120 = qJD(3) * t104;
t119 = qJD(1) * qJD(3);
t117 = -t111 * t99 + t114 * t98;
t101 = (t111 * t115 + t112 * t114) * qJD(1);
t113 = cos(qJ(5));
t110 = sin(qJ(5));
t108 = qJD(3) + qJD(4);
t107 = -qJD(1) * pkin(1) + qJD(2);
t102 = (-t111 * t112 + t114 * t115) * qJD(1);
t100 = qJD(5) + t101;
t95 = t113 * t102 + t110 * t108;
t94 = t110 * t102 - t113 * t108;
t93 = t101 * pkin(4) - t102 * pkin(8) + t103;
t92 = t108 * pkin(8) + t122;
t91 = -t108 * pkin(4) - t117;
t1 = [t123, 0, 0, t107 * qJD(1), t121, qJ(2) ^ 2 * t123 + t107 ^ 2 / 0.2e1, t115 ^ 2 * t123, -t115 * t116 * t112, t115 * t119, -t112 * t119, qJD(3) ^ 2 / 0.2e1, t112 * t121 + t115 * t120, -t112 * t120 + t115 * t121, t102 ^ 2 / 0.2e1, -t102 * t101, t102 * t108, -t101 * t108, t108 ^ 2 / 0.2e1, t103 * t101 + t117 * t108, t103 * t102 - t122 * t108, t95 ^ 2 / 0.2e1, -t95 * t94, t95 * t100, -t94 * t100, t100 ^ 2 / 0.2e1, (-t110 * t92 + t113 * t93) * t100 + t91 * t94, -(t110 * t93 + t113 * t92) * t100 + t91 * t95;];
T_reg = t1;
